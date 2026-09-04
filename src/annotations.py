# vim: set expandtab shiftwidth=4 softtabstop=4:

# === UCSF ChimeraX Copyright ===
# Copyright 2016 Regents of the University of California.
# All rights reserved.  This software provided pursuant to a
# license agreement containing restrictions on its disclosure,
# duplication and use.  For details see:
# http://www.rbvi.ucsf.edu/chimerax/docs/licensing.html
# This notice must be embedded in or attached to all copies,
# including partial copies, of the software or any revisions
# or derivations thereof.
# === UCSF ChimeraX Copyright ===

"""
Free-text residue annotations for the Investigate tool - see claude_idee.md for
the design. A JSON file (<structure-name>.chopchop.json, in ChopChopMF's shared
download folder by default - get_settings(session).download_dir, same folder
every other ChopChopMF tool already remembers - overridable per model via
set_store_path()) is the source of truth, not the live ChimeraX session: it
survives a plain close/reopen of the structure, which ChimeraX custom residue
attributes do not (those only survive an explicit ChimeraX session
save/restore). One AnnotationStore per file path is cached per ChimeraX
session so every tool/dialog touching the same structure shares one in-memory
copy instead of racing separate writers.
"""

from datetime import datetime
from pathlib import Path
import json

_SESSION_STORES = {}       # session -> {path: AnnotationStore}
_SESSION_MODEL_PATH = {}   # session -> {model: path}  (explicit user overrides only)
_CLEANUP_REGISTERED = set()  # sessions whose REMOVE_MODELS cleanup handler is already registered


def _ensure_cleanup_handler(session):
    """A closed model kept as a key in _SESSION_MODEL_PATH would otherwise
    stay referenced (and alive) for the rest of the ChimeraX session, even
    though nothing else needs it once the model is gone - drop its entry as
    soon as ChimeraX reports it removed."""
    if session in _CLEANUP_REGISTERED:
        return
    _CLEANUP_REGISTERED.add(session)
    from chimerax.core.models import REMOVE_MODELS

    def _on_models_removed(_trigger_name, removed_models, session=session):
        model_paths = _SESSION_MODEL_PATH.get(session)
        if model_paths:
            for model in removed_models:
                model_paths.pop(model, None)

    session.triggers.add_handler(REMOVE_MODELS, _on_models_removed)


def _file_stem(model):
    model_path = getattr(model, "filename", None)
    if model_path:
        return Path(model_path).stem
    safe_name = "".join(c if c.isalnum() or c in "-_" else "_" for c in model.name)
    return safe_name or f"model_{model.id_string}"


def _default_path(session, model):
    from .utils import get_settings
    download_dir = Path(get_settings(session).download_dir)
    return download_dir / f"{_file_stem(model)}.chopchop.json"


class AnnotationStore:
    """One current note per residue, keyed by "chain_id:residue_number" - a
    single editable value (not an append-only log), so it can be edited
    directly in the Chart tab's table as well as from the Residue tab.

    Also holds a durable copy of whatever per-residue *values* any ChopChopMF
    tool has computed (PDBePISA's delta_g_score, PAE Analysis's pdockq, etc.) -
    a live ChimeraX residue attribute only exists while that model is open in
    this session (it does not survive a plain close/reopen of the structure
    file), so Investigate's Chart snapshots whatever is currently live into
    this file every time it builds its rows. That happens generically, once,
    in Investigate - not something every individual tool needs to remember to
    do itself."""

    def __init__(self, path):
        self.path = path
        self._data = {}
        self._values = {}
        if path.exists():
            try:
                with open(path) as f:
                    raw = json.load(f)
                self._data = raw.get("residues", {})
                self._values = raw.get("values", {})
            except (OSError, ValueError):
                self._data = {}
                self._values = {}

    @staticmethod
    def _key(chain_id, number):
        return f"{chain_id}:{number}"

    def note_for(self, chain_id, number):
        """The current note dict ({"tool", "text", "ts"}) for a residue, or None."""
        return self._data.get(self._key(chain_id, number))

    def set_note(self, chain_id, number, tool_name, text):
        """Replace the note for a residue and save to disk immediately
        (annotations represent the user's own analysis, not a regeneratable
        score, so autosave beats risking loss to a forgotten Save click).
        An empty/whitespace-only text removes the note."""
        key = self._key(chain_id, number)
        text = text.strip()
        if not text:
            self._data.pop(key, None)
        else:
            self._data[key] = {"tool": tool_name, "text": text, "ts": datetime.now().isoformat(timespec="seconds")}
        self._save()

    def all_notes(self):
        return dict(self._data)

    def values_for(self, chain_id, number):
        """Dict of {attr_name: value} last recorded for this residue (from
        any tool), or {} if none - the durable fallback used when a residue's
        live ChimeraX attribute is gone (e.g. after closing and reopening the
        structure) but was recorded here at some earlier point."""
        return dict(self._values.get(self._key(chain_id, number), {}))

    def record_values(self, residue_values):
        """Bulk-update many residues' values in a single save.
        `residue_values` is {(chain_id, number): {attr_name: value, ...}, ...}.
        Only actually writes to disk if something changed, so calling this on
        every Chart refresh doesn't thrash the file when nothing's new."""
        changed = False
        for (chain_id, number), values in residue_values.items():
            if not values:
                continue
            key = self._key(chain_id, number)
            if self._values.get(key) != values:
                self._values[key] = dict(values)
                changed = True
        if changed:
            self._save()

    def clear_stale_values(self, stale):
        """Bulk-remove specific attribute names from many residues' stored
        values in a single save. `stale` is {(chain_id, number): {attr_name,
        ...}, ...}. Used when a tool's current live state proves a
        previously-recorded value is no longer applicable for a residue (e.g.
        it fell out of the current chain pair's interface, or fell out of
        scope after a scoping bugfix) - without this, a value recorded once
        (even mistakenly, by an earlier bug) stays in the file forever, since
        record_values() only ever adds/replaces, never removes."""
        changed = False
        for (chain_id, number), names in stale.items():
            if not names:
                continue
            key = self._key(chain_id, number)
            values = self._values.get(key)
            if not values:
                continue
            for name in names:
                if name in values:
                    del values[name]
                    changed = True
            if not values:
                self._values.pop(key, None)
        if changed:
            self._save()

    def known_value_names(self):
        """Every attribute name ever recorded for any residue in this store -
        lets a Chart column stay visible even when no residue currently has
        the live ChimeraX attribute (e.g. right after reopening a structure)."""
        names = set()
        for values in self._values.values():
            names.update(values.keys())
        return names

    def _save(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with open(self.path, "w") as f:
            json.dump({"residues": self._data, "values": self._values}, f, indent=2, default=str)


def get_store(session, model):
    """The shared AnnotationStore for `model`, for this session - at
    set_store_path()'s path if one was chosen for this model, else the default
    (get_settings(session).download_dir)/<structure-name>.chopchop.json."""
    _ensure_cleanup_handler(session)
    model_paths = _SESSION_MODEL_PATH.setdefault(session, {})
    path = model_paths.get(model) or _default_path(session, model)
    return _store_at(session, path)


def set_store_path(session, model, path):
    """Point `model`'s annotations at `path` (e.g. chosen via a file dialog)
    for the rest of this session, loading any existing notes already there."""
    _ensure_cleanup_handler(session)
    _SESSION_MODEL_PATH.setdefault(session, {})[model] = path
    return _store_at(session, path)


def _store_at(session, path):
    stores = _SESSION_STORES.setdefault(session, {})
    if path not in stores:
        stores[path] = AnnotationStore(path)
    return stores[path]
