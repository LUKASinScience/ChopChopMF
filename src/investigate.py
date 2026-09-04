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
Investigate Tool
Cross-tool residue "dossier": pick a model, then Ctrl+click a residue in the 3D
view (or use any other ChopChopMF tool's Select action) to see every custom
ChimeraX residue attribute any tool has recorded for it (e.g. PDBePISA's
residue_score, ChopMissense's MissenseScores), plus your own free-text note.
The Chart tab is popup-only (a big table gets cramped in a docked panel): it
lists every residue of the model at once, with pLDDT color-coded, an editable
Notes column, CSV export, and a copy-out of the notes file.
See claude_idee.md for the design this implements (v1 scope).
"""

import csv
import io
from pathlib import Path

import numpy as np
from matplotlib.colors import LinearSegmentedColormap

from chimerax.core.tools import ToolInstance
from chimerax.core.commands import run
from chimerax.atomic import selected_residues
from chimerax.ui import MainToolWindow
from Qt.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QWidget, QComboBox,
    QTabWidget, QLineEdit, QTableWidget, QTableWidgetItem, QFileDialog, QDialog,
    QGroupBox, QInputDialog,
)
from Qt.QtGui import QFont, QColor
from Qt.QtCore import Qt

from .utils import (
    make_scrollable, make_guide_button, show_error, open_figure_dialog, GUIDE_BASE_URL, track_dialog,
    PDOCKQ_THRESHOLDS, PDOCKQ_LOW, ILIS_THRESHOLD, ILIS_HIGH, threshold_badge, pae_value_colors,
    get_settings,
)
from . import annotations

# Same pLDDT confidence bands as ChimeraX's own AlphaFold palette / the PAE
# Analysis tool's Plots tab (0-50 red, 50-70 orange, 70-90 yellow, 90-100 blue).
_PLDDT_BANDS = [
    (0, 50, "#FF0000", "#FFFFFF"),
    (50, 70, "#FFA500", "#000000"),
    (70, 90, "#FFFF00", "#000000"),
    (90, 100, "#6495ED", "#000000"),
]

# Single-hue, light->dark sequential ramp for the coverage heatmap (magnitude,
# not identity/category, so one hue is correct here) - light gray for "no
# value yet" (NaN) reads clearly against the teal brand color for "present".
_COVERAGE_CMAP = LinearSegmentedColormap.from_list("chopchop_coverage", ["#F2F2F2", "#00695C"])
_COVERAGE_CMAP.set_bad(color="#DDDDDD")


# Residue attributes the *other* ChopChopMF analysis tools are known to write,
# with a clear display label - always shown as Chart columns, blank if that
# tool hasn't been run on this model yet, so the table itself shows the user
# which analyses are still worth running (not just whatever already happens
# to be set). Anything else found on a residue still gets its own column too
# (raw attribute name as the label), so nothing is ever silently hidden.
_KNOWN_ATTRS = {
    "residue_score": "PDBePISA class",
    "delta_g_score": "PDBePISA ΔG",
    "MissenseScores": "AlphaMissense",
    # "Interface", not "PAE", for the three that are purely distance/pLDDT-based
    # and never actually read the PAE matrix (pDockQ, buried area, H-bonds all
    # come from contact geometry - only LIS/cLIS/iLIS/ipSAE below use PAE
    # values themselves) - a "PAE" prefix on these implied a PAE dependency
    # that isn't there and was reported as confusing.
    "chopchop_pae_contact": "Interface contact",
    "chopchop_pae_clis_confident": "PAE cLIS",
    "chopchop_pae_pdockq": "Interface pDockQ",
    "chopchop_pae_buried_area": "Interface Buried area",
    "chopchop_pae_hbonds": "Interface H-bonds",
    "chopchop_pae_lis": "PAE LIS",
    "chopchop_pae_clis": "PAE cLIS score",
    "chopchop_pae_ilis": "PAE iLIS",
    "chopchop_pae_ipsae": "PAE ipSAE (d0chn)",
    "chopchop_pae_contact_pae": "PAE (Å)",
    "chopchop_alphasync_asa": "AlphaSync SASA",
    "chopchop_alphasync_rsa": "AlphaSync RSA",
    "chopchop_alphasync_surface": "AlphaSync Surface",
    "chopchop_alphasync_disorder": "AlphaSync Disorder",
    "chopchop_alphasync_secstruct": "AlphaSync Sec. Str.",
    "chopchop_phospho_kinase": "Phospho kinase(s)",
    "chopchop_signal_organelle": "Signal/targeting organelle",
    "chopchop_tm_helix": "TM helix",
}

# ChimeraX core itself (chimerax.atomic.structure.py, "Local model scoring")
# registers a residue attribute named "<metric>_score" for any per-residue QA
# metric it finds embedded in an opened mmCIF file - for AlphaFold-style files
# whose local metric is literally named "pLDDT", that comes out as the exact
# attribute "pLDDT_score". It's not from ChopChopMF and not fake test data -
# it's a real, ChimeraX-registered duplicate of the pLDDT already shown from
# the structure's B-factor column, so it's excluded from Chart/Plots discovery
# here rather than shown as a second, confusingly-similar column. Any other
# "<something>_score" metric a file might embed is a different metric and
# still shown normally - only this exact, guaranteed-duplicate name is hidden.
_HIDDEN_ATTRS = {"pLDDT_score"}

# One-line, plain-language gloss per known column, for the "Export for AI
# analysis" file - lets a pasted-in LLM interpret the table without the user
# having to explain every column by hand first.
_COLUMN_DESCRIPTIONS = {
    "pLDDT": "AlphaFold's per-residue confidence (0-100). >90 very high, "
             "70-90 confident, 50-70 low, <50 very low/likely disordered.",
    "PDBePISA class": "PISA interface classification for this residue (PDBePISA tool).",
    "PDBePISA ΔG": "Solvation free energy gain (ΔG, kcal/mol) on interface formation "
                   "for this residue, from PDBePISA - more negative usually means a "
                   "more stabilizing contribution to the interface.",
    "AlphaMissense": "AlphaMissense pathogenicity score (0-1). Higher = more likely "
                      "pathogenic if this residue were mutated.",
    "Interface contact": "True if this residue was part of a contact pair with the other "
                    "chain in PAE Analysis (pDockQ's contact definition: C-beta to C-beta, "
                    "or C-alpha for glycine, distance <= 8 Å - distance-only, no PAE filter; "
                    "see PAE cLIS for a PAE-filtered contact).",
    "PAE cLIS": "True if this residue is part of a cLIS-qualifying contact pair "
                "(contributes to the confident-interface score iLIS).",
    "Interface pDockQ": "pDockQ score (Bryant, Pozzati & Elofsson 2022) for the chain-pair "
                  "interface this residue was last a contact-residue of in PAE Analysis - "
                  ">0.5 high confidence, 0.23-0.5 weak/medium, <0.23 poor. One shared value "
                  "for the whole interface (not residue-specific), only set on the actual "
                  "contact residues (see Interface contact) - not on every residue of both "
                  "chains. Derived from contact count and pLDDT only, not from PAE values.",
    "Interface Buried area": "Buried interface area (Å²) for the chain-pair interface this "
                        "residue was last a contact-residue of in PAE Analysis (ChimeraX's "
                        "measure buriedarea) - one shared value for the whole interface, "
                        "only set on the actual contact residues.",
    "Interface H-bonds": "Hydrogen bond count for the chain-pair interface this residue was "
                   "last a contact-residue of in PAE Analysis (ChimeraX's hbonds command) - "
                   "one shared value for the whole interface, only set on the actual contact "
                   "residues.",
    "PAE LIS": "LIS score (Kim et al. 2024) for the chain-pair interface this residue was "
               "last a contact-residue of in PAE Analysis - one shared value for the whole "
               "interface, only set on the actual contact residues.",
    "PAE cLIS score": "cLIS score (Kim et al. 2024, distance-filtered LIS) for the chain-pair "
                       "interface this residue was last a contact-residue of in PAE Analysis - "
                       "one shared value for the whole interface, only set on the actual "
                       "contact residues.",
    "PAE iLIS": "iLIS score (Kim et al. 2024) for the chain-pair interface this residue was "
                "last a contact-residue of in PAE Analysis - >=0.223 is a high-confidence "
                "interaction. One shared value for the whole interface, only set on the "
                "actual contact residues.",
    "PAE ipSAE (d0chn)": "ipSAE, fixed chain-length variant (Dunbrack, DunbrackLab/IPSAE) for "
                          "the chain-pair interface this residue was last a contact-residue "
                          "of in PAE Analysis - one shared value for the whole interface, "
                          "only set on the actual contact residues.",
    "PAE (Å)": "The actual predicted aligned error (Å), not a derived score - averaged over "
               "this residue's contacts (<=8 Å) with the other chain, same direction "
               "ChimeraX itself uses to color its PAE contact pseudobonds. Roughly <=12 Å is "
               "confident (PAE Analysis Tab 1's own default contact cutoff), ~20 Å or higher "
               "is uncertain. Only set for residues that are actually in contact.",
    "AlphaSync SASA": "Solvent-accessible surface area in Å² (AlphaSync, DSSP algorithm).",
    "AlphaSync RSA": "Relative solvent-accessible surface area (0-1 fraction, AlphaSync).",
    "AlphaSync Surface": "Surface/core classification from AlphaSync (<=25% RSA = buried).",
    "AlphaSync Disorder": "AlphaSync disorder call for this residue ('*' disordered, '.' structured).",
    "AlphaSync Sec. Str.": "AlphaSync secondary structure assignment for this residue.",
    "Notes": "The user's own free-text note for this residue, written in Investigate.",
    "Phospho kinase(s)": "Kinase consensus motif(s) matched at this Ser/Thr/Tyr residue "
                          "(Phospho Sites tool) - a heuristic motif match filtered by pLDDT/SASA, "
                          "not a calibrated prediction; see docs/acknowledgements.md for citations.",
    "Signal/targeting organelle": "Subcellular targeting-signal motif(s) matched at this residue "
                                   "(Signal Peptide tool) - a heuristic motif match filtered by SASA, "
                                   "not a calibrated prediction; see docs/acknowledgements.md for citations.",
    "TM helix": "True if this residue is part of a candidate transmembrane helix (TM Helix tool: "
                "DSSP alpha-helix + hydrophobicity threshold) - a heuristic call, not a calibrated "
                "TMHMM-grade prediction.",
}


def _attr_columns(model, store=None):
    """(attribute names, display labels), same order - known ChopChopMF
    attributes first (always present, even if no residue has them set yet),
    then any other custom attribute actually found live on the model, or
    (if `store` is given) previously recorded in its annotation store - so a
    column stays visible even right after reopening a structure, before any
    tool has recomputed anything in the new session (except _HIDDEN_ATTRS -
    see its comment)."""
    discovered = {name for r in model.residues for name, _ in r.custom_attrs}
    if store is not None:
        discovered |= store.known_value_names()
    extra = sorted(n for n in discovered if n not in _KNOWN_ATTRS and n not in _HIDDEN_ATTRS)
    names = list(_KNOWN_ATTRS) + extra
    labels = [_KNOWN_ATTRS.get(n, n) for n in names]
    return names, labels


def _plddt_colors(value):
    for lo, hi, bg, fg in _PLDDT_BANDS:
        if lo <= value <= hi:
            return bg, fg
    return None, None


# Attributes with an established, citable confidence threshold (same bands as
# PAE Analysis's own Scores tab, via the shared utils.py helpers) - colored the
# same way here so a score means the same thing wherever it's shown. Every
# other attribute (buried area, H-bonds, LIS/cLIS, ipSAE d0chn, ...) has no
# literature-backed good/bad cutoff and intentionally stays uncolored.
_ATTR_BADGE_THRESHOLDS = {
    "chopchop_pae_pdockq": (PDOCKQ_THRESHOLDS, PDOCKQ_LOW),
    "chopchop_pae_ilis": ([(ILIS_THRESHOLD, *ILIS_HIGH)], None),
}

# Attributes that are a real PAE value (Angstrom) rather than a derived score -
# colored on the same continuous blue-to-white scale as the PAE Matrix plot and
# ChimeraX's own PAE contact pseudobonds, instead of a two/three-band badge.
_ATTR_PAE_VALUE = {"chopchop_pae_contact_pae"}


def _attr_badge_colors(attr_name, value):
    if attr_name in _ATTR_PAE_VALUE:
        try:
            return pae_value_colors(float(value))
        except (TypeError, ValueError):
            return None, None
    spec = _ATTR_BADGE_THRESHOLDS.get(attr_name)
    if spec is None:
        return None, None
    try:
        numeric_value = float(value)
    except (TypeError, ValueError):
        return None, None
    thresholds, low = spec
    badge = threshold_badge(numeric_value, thresholds, low)
    if badge is None:
        return None, None
    bg, fg, _label = badge
    return bg, fg


class Investigate(ToolInstance):

    SESSION_ENDURING = False
    SESSION_SAVE = False

    def __init__(self, session, tool_name):
        super().__init__(session, tool_name)
        self.display_name = "Investigate"
        self._current_residue = None
        self._open_dialogs = []
        self._scatter_residues = []
        self.tool_window = MainToolWindow(self)
        self._build_ui()
        self._selection_handler = session.triggers.add_handler(
            "selection changed", self._on_selection_changed
        )
        self.tool_window.manage("side")

    def delete(self):
        if self._selection_handler is not None:
            self.session.triggers.remove_handler(self._selection_handler)
            self._selection_handler = None
        super().delete()

    def _build_ui(self):
        layout = QVBoxLayout()
        layout.setSpacing(4)
        layout.setContentsMargins(6, 6, 6, 6)

        model_row = QHBoxLayout()
        model_row.addWidget(QLabel("Model:"))
        self.model_selector = QComboBox()
        self.model_selector.currentIndexChanged.connect(self._refresh_annotations_path_label)
        self.model_selector.currentIndexChanged.connect(self._refresh_plot_axis_choices)
        model_row.addWidget(self.model_selector, stretch=1)
        refresh_button = QPushButton("↻")
        refresh_button.setToolTip("Refresh model list")
        refresh_button.clicked.connect(self._refresh_models)
        model_row.addWidget(refresh_button)
        layout.addLayout(model_row)

        file_row = QHBoxLayout()
        self.annotations_path_label = QLabel("")
        self.annotations_path_label.setStyleSheet("color: gray; font-size: 10px;")
        self.annotations_path_label.setWordWrap(True)
        self.annotations_path_label.setToolTip(
            "To save a timestamped snapshot of this file, or load an earlier one back, "
            "use ChopChopMF's Settings toolbar tool."
        )
        file_row.addWidget(self.annotations_path_label, stretch=1)
        refresh_annotations_button = QPushButton("↻")
        refresh_annotations_button.setToolTip(
            "Refresh this label - useful after changing this model's annotations file "
            "from ChopChopMF's Settings toolbar tool."
        )
        refresh_annotations_button.clicked.connect(self._refresh_annotations_path_label)
        file_row.addWidget(refresh_annotations_button)
        layout.addLayout(file_row)

        tabs = QTabWidget()
        tabs.addTab(make_scrollable(self._create_residue_tab()), "Residue")
        tabs.addTab(make_scrollable(self._create_chart_tab()), "Chart")
        tabs.addTab(make_scrollable(self._create_plots_tab()), "Plots")
        layout.addWidget(tabs)

        self.tool_window.ui_area.setLayout(layout)
        self._refresh_models()

    def _refresh_models(self):
        self.model_selector.clear()
        for m in self.session.models.list():
            if hasattr(m, "residues"):
                self.model_selector.addItem(f"#{m.id_string} {m.name}", m.id_string)
        self._refresh_annotations_path_label()

    def _refresh_annotations_path_label(self):
        model = self._current_model()
        if model is None:
            self.annotations_path_label.setText("No model selected.")
            return
        store = annotations.get_store(self.session, model)
        self.annotations_path_label.setText(f"Annotations file: {store.path}")

    def _current_model(self):
        id_string = self.model_selector.currentData()
        if not id_string:
            return None
        return next((m for m in self.session.models.list() if m.id_string == id_string), None)

    def _clear_layout(self, layout):
        while layout.count():
            item = layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()

    # ---------- Residue tab ----------

    def _create_residue_tab(self):
        widget = QWidget()
        layout = QVBoxLayout()
        layout.setAlignment(Qt.AlignTop)
        layout.setSpacing(3)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.addWidget(make_guide_button("4-analyze-structure"))
        click_hint = QLabel("Ctrl+click an atom to select it (plain click just rotates the view).")
        click_hint.setStyleSheet("color: gray; font-size: 10px;")
        layout.addWidget(click_hint)

        self.residue_heading = QLabel("No residue selected.")
        bold_font = QFont()
        bold_font.setBold(True)
        self.residue_heading.setFont(bold_font)
        self.residue_heading.setStyleSheet(
            "background-color: palette(mid); color: palette(text); border-radius: 4px; padding: 5px;"
        )
        layout.addWidget(self.residue_heading)

        self.attrs_container = QWidget()
        self.attrs_layout = QVBoxLayout()
        self.attrs_layout.setAlignment(Qt.AlignTop)
        self.attrs_layout.setSpacing(1)
        self.attrs_layout.setContentsMargins(0, 0, 0, 0)
        self.attrs_container.setLayout(self.attrs_layout)
        layout.addWidget(self.attrs_container)

        note_row = QHBoxLayout()
        note_row.addWidget(QLabel("Note:"))
        self.note_input = QLineEdit()
        self.note_input.setPlaceholderText("Add a note for this residue…")
        note_row.addWidget(self.note_input, stretch=1)
        save_note_button = QPushButton("Save")
        save_note_button.clicked.connect(self._save_note)
        note_row.addWidget(save_note_button)
        layout.addLayout(note_row)

        self.note_meta_label = QLabel("")
        self.note_meta_label.setStyleSheet("color: gray; font-size: 10px;")
        layout.addWidget(self.note_meta_label)

        widget.setLayout(layout)
        return widget

    def _on_selection_changed(self, trigger_name, _data):
        residues = selected_residues(self.session)
        if len(residues) != 1:
            self._current_residue = None
            self.residue_heading.setStyleSheet(
                "background-color: palette(mid); color: palette(text); border-radius: 4px; padding: 5px;"
            )
            if len(residues) == 0:
                self.residue_heading.setText("No residue selected.")
            else:
                self.residue_heading.setText(f"{len(residues)} residues selected - Ctrl+click just one to inspect.")
            self._clear_layout(self.attrs_layout)
            self._refresh_note()
            return
        self._show_residue(residues[0])

    def _add_attr_row(self, name, value_text, bg=None, fg=None):
        row_widget = QWidget()
        row = QHBoxLayout()
        row.setContentsMargins(4, 1, 4, 1)
        name_label = QLabel(name)
        bold_font = QFont()
        bold_font.setBold(True)
        name_label.setFont(bold_font)
        row.addWidget(name_label)
        value_label = QLabel(value_text)
        if bg:
            value_label.setStyleSheet(f"background-color: {bg}; color: {fg}; padding: 1px 6px; border-radius: 3px;")
        row.addWidget(value_label)
        row.addStretch()
        row_widget.setLayout(row)
        self.attrs_layout.addWidget(row_widget)

    def _show_residue(self, residue):
        self._current_residue = residue
        self.residue_heading.setStyleSheet(
            "background-color: #2a78d6; color: white; border-radius: 4px; padding: 5px;"
        )
        self.residue_heading.setText(f"Chain {residue.chain_id}, residue {residue.number} ({residue.name})")

        self._clear_layout(self.attrs_layout)
        atom = residue.principal_atom
        if atom is not None:
            bg, fg = _plddt_colors(atom.bfactor)
            self._add_attr_row("pLDDT", f"{atom.bfactor:.1f}", bg, fg)

        attrs = [(n, v) for n, v in residue.custom_attrs if n not in _HIDDEN_ATTRS]
        for name, value in sorted(attrs):
            bg, fg = _attr_badge_colors(name, value)
            self._add_attr_row(name, str(value), bg, fg)
        if not attrs and atom is None:
            self.attrs_layout.addWidget(QLabel("No tool has recorded a value for this residue yet."))

        self._refresh_note()

    def _refresh_note(self):
        model = self._current_model()
        if self._current_residue is None or model is None:
            self.note_input.setText("")
            self.note_meta_label.setText("")
            return
        store = annotations.get_store(self.session, model)
        note = store.note_for(self._current_residue.chain_id, self._current_residue.number)
        if note:
            self.note_input.setText(note["text"])
            self.note_meta_label.setText(f"Last edited {note['ts']} via {note['tool']}")
        else:
            self.note_input.setText("")
            self.note_meta_label.setText("No note yet.")

    def _save_note(self):
        if self._current_residue is None:
            show_error(self.tool_window.ui_area, "ChopChopMF", "Ctrl+click a residue in the 3D view first.")
            return
        model = self._current_model()
        if model is None:
            show_error(self.tool_window.ui_area, "ChopChopMF", "Select a model first.")
            return
        store = annotations.get_store(self.session, model)
        try:
            store.set_note(self._current_residue.chain_id, self._current_residue.number, "Investigate", self.note_input.text())
        except OSError as e:
            show_error(self.tool_window.ui_area, "ChopChopMF", f"Failed to save note:\n{e}")
            return
        self._refresh_note()

    # ---------- Chart tab (popup-only - a full residue table doesn't fit a docked panel) ----------

    def _create_chart_tab(self):
        widget = QWidget()
        layout = QVBoxLayout()
        layout.setAlignment(Qt.AlignTop)
        layout.setSpacing(4)
        layout.setContentsMargins(6, 6, 6, 6)
        info_label = QLabel(
            "Every residue of the selected model, with pLDDT, every value any tool has "
            "recorded, and an editable Notes column."
        )
        info_label.setWordWrap(True)
        layout.addWidget(info_label)

        open_button = QPushButton("Open Chart")
        open_button.clicked.connect(self._open_chart_window)
        layout.addWidget(open_button)

        widget.setLayout(layout)
        return widget

    def _merge_stored_and_live(self, model, store):
        """{residue: merged_values_dict} for every residue of `model` - live
        ChimeraX attributes merged with whatever's durably stored, live wins.

        Also reconciles the store, and this is the part that matters: an
        attribute name currently live on AT LEAST ONE residue of `model`
        means the tool that owns it has (re)computed in this ChimeraX
        session, so its current live state is authoritative for every
        residue - a residue that no longer carries that attribute (fell out
        of the current chain pair's interface, or fell out of scope after a
        scoping fix) has the stale leftover actively removed from the store,
        not just silently omitted from what's displayed. Without this, a
        value written once - even by an earlier bug - stays in the durable
        file forever and keeps reappearing here, merged back in as soon as
        the live attribute disappears, no matter how many times the live
        computation is fixed and rerun.

        An attribute name that is NOT live anywhere in `model` right now was
        only ever recorded in an earlier session (this one just reopened the
        structure and nothing has recomputed yet) - left untouched, so
        durable persistence across a plain close/reopen still works."""
        live_attr_names = {name for r in model.residues for name, _ in r.custom_attrs}
        snapshot, stale, merged_by_residue = {}, {}, {}
        for r in model.residues:
            key = (r.chain_id, r.number)
            live_attrs = dict(r.custom_attrs)  # list of (name, value) pairs, not a dict
            stored = store.values_for(r.chain_id, r.number)
            stale_names = (live_attr_names & stored.keys()) - live_attrs.keys()
            if stale_names:
                stale[key] = stale_names
                for name in stale_names:
                    stored.pop(name, None)
            if live_attrs:
                snapshot[key] = live_attrs
            merged_by_residue[r] = {**stored, **live_attrs}
        if stale:
            store.clear_stale_values(stale)
        if snapshot:
            store.record_values(snapshot)
        return merged_by_residue

    def _build_chart_rows(self, model=None):
        """Every residue of the model, as (headers, rows, row_keys) - rows are
        plain strings, ready for a QTableWidget/CSV; row_keys are the matching
        (chain_id, number) pairs, same order, used to write an edited Notes
        cell back to the store. Always includes pLDDT (the structure's own
        B-factor column - no other tool needs to have run first) plus one
        column per _KNOWN_ATTRS entry (blank if that tool hasn't run on this
        model yet - the point is to show which analyses are still available,
        not just what already happened to be set) plus Notes.

        `model`, if given, overrides the main panel's currently-selected model -
        the Chart popup keeps its own model selector so it can show a
        different model than what's selected in the background Residue tab.

        See _merge_stored_and_live() for how a value is resolved (live wins
        over stored) and how stale stored values get cleaned up."""
        if model is None:
            model = self._current_model()
        if model is None:
            return [], [], []
        store = annotations.get_store(self.session, model)
        attr_names, attr_labels = _attr_columns(model, store)
        merged_by_residue = self._merge_stored_and_live(model, store)

        # Populated columns (this model actually has a value for them) grouped
        # together on the left, ahead of the still-blank ones - lets you scan
        # "what's actually been run on this model" at a glance instead of
        # hunting across a wide, mostly-empty fixed column order.
        has_value = {
            name: any(merged.get(name) not in (None, "") for merged in merged_by_residue.values())
            for name in attr_names
        }
        ordered = sorted(range(len(attr_names)), key=lambda i: not has_value[attr_names[i]])
        attr_names = [attr_names[i] for i in ordered]
        attr_labels = [attr_labels[i] for i in ordered]

        headers = ["Chain", "Residue #", "Name", "pLDDT"] + attr_labels + ["Notes"]
        rows, row_keys = [], []
        for r in sorted(model.residues, key=lambda r: (r.chain_id, r.number)):
            note = store.note_for(r.chain_id, r.number)
            merged = merged_by_residue[r]
            atom = r.principal_atom
            plddt = f"{atom.bfactor:.1f}" if atom is not None else ""
            row = [r.chain_id, str(r.number), r.name, plddt]
            row.extend(str(merged.get(name, "")) for name in attr_names)
            row.append(note["text"] if note else "")
            rows.append(row)
            row_keys.append((r.chain_id, r.number))
        return headers, rows, row_keys

    def _open_chart_window(self):
        model = self._current_model()
        if model is None:
            show_error(self.tool_window.ui_area, "ChopChopMF", "Select a model first.")
            return

        dialog = QDialog(self.tool_window.ui_area)
        dialog.setAttribute(Qt.WA_DeleteOnClose)
        dialog.resize(1100, 700)
        layout = QVBoxLayout()

        # The Chart popup keeps its own model selector, independent of the
        # Residue tab's - so you can flip between models to compare their
        # Charts side by side without closing/reopening this window or
        # disturbing whatever's selected in the background panel.
        model_row = QHBoxLayout()
        model_row.addWidget(QLabel("Model:"))
        popup_model_selector = QComboBox()
        for m in self.session.models.list():
            if hasattr(m, "residues"):
                popup_model_selector.addItem(f"#{m.id_string} {m.name}", m.id_string)
        idx = popup_model_selector.findData(model.id_string)
        if idx >= 0:
            popup_model_selector.setCurrentIndex(idx)
        model_row.addWidget(popup_model_selector, stretch=1)
        layout.addLayout(model_row)

        def current_popup_model():
            id_string = popup_model_selector.currentData()
            if not id_string:
                return None
            return next((m for m in self.session.models.list() if m.id_string == id_string), None)

        search_row = QHBoxLayout()
        search_row.addWidget(QLabel("Search residue #:"))
        search_input = QLineEdit()
        search_input.setPlaceholderText("e.g. 117 - filters as you type")
        search_row.addWidget(search_input)
        layout.addLayout(search_row)

        table = QTableWidget(0, 0)
        table.verticalHeader().setVisible(False)
        layout.addWidget(table)
        layout.addWidget(QLabel("Only the Notes column (last) is editable - edit a cell and press Enter to save."))

        state = {"row_keys": [], "notes_col": None, "plddt_col": 3}

        def populate():
            m = current_popup_model()
            dialog.setWindowTitle(f"ChopChopMF - Investigate Chart - {m.name if m else 'no model'}")
            if m is None:
                table.setRowCount(0)
                table.setColumnCount(0)
                return
            headers, rows, row_keys = self._build_chart_rows(m)
            state["row_keys"] = row_keys
            state["notes_col"] = len(headers) - 1
            badge_cols = {
                headers.index(_KNOWN_ATTRS[attr_name]): attr_name
                for attr_name in (*_ATTR_BADGE_THRESHOLDS, *_ATTR_PAE_VALUE)
                if _KNOWN_ATTRS[attr_name] in headers
            }
            table.blockSignals(True)
            table.clear()
            table.setColumnCount(len(headers))
            table.setRowCount(len(rows))
            table.setHorizontalHeaderLabels(headers)
            for i, row in enumerate(rows):
                for col, value in enumerate(row):
                    item = QTableWidgetItem(value)
                    if col != state["notes_col"]:
                        item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                    bg = fg = None
                    if col == state["plddt_col"] and value:
                        bg, fg = _plddt_colors(float(value))
                    elif col in badge_cols and value:
                        bg, fg = _attr_badge_colors(badge_cols[col], value)
                    if bg:
                        item.setBackground(QColor(bg))
                        item.setForeground(QColor(fg))
                    table.setItem(i, col, item)
            table.resizeColumnsToContents()
            table.blockSignals(False)
            apply_search_filter()

        def apply_search_filter():
            text = search_input.text().strip()
            resnum_col = 1  # Chain, Residue #, Name, pLDDT, ...
            for row in range(table.rowCount()):
                item = table.item(row, resnum_col)
                match = not text or (item is not None and text in item.text())
                table.setRowHidden(row, not match)

        search_input.textChanged.connect(apply_search_filter)
        popup_model_selector.currentIndexChanged.connect(populate)

        def on_item_changed(item):
            notes_col = state["notes_col"]
            row_keys = state["row_keys"]
            if notes_col is None or item.column() != notes_col or item.row() >= len(row_keys):
                return
            m = current_popup_model()
            if m is None:
                return
            chain_id, number = row_keys[item.row()]
            store = annotations.get_store(self.session, m)
            try:
                store.set_note(chain_id, number, "Investigate Chart", item.text())
            except OSError as e:
                show_error(dialog, "ChopChopMF", f"Failed to save note:\n{e}")
            if self._current_residue is not None and (self._current_residue.chain_id, self._current_residue.number) == (chain_id, number):
                self._refresh_note()

        table.itemChanged.connect(on_item_changed)
        populate()

        # Widen the window (once, on open) so every column - including any
        # ChopChopMF adds in the future - is visible without horizontal
        # scrolling, capped to the screen so it never overflows the display.
        from Qt.QtWidgets import QApplication
        content_width = table.horizontalHeader().length() + table.verticalHeader().width() + 50
        screen = QApplication.primaryScreen()
        max_width = screen.availableGeometry().width() - 100 if screen else 1800
        dialog.resize(min(content_width, max_width), dialog.height())

        button_row = QHBoxLayout()
        refresh_button = QPushButton("↻ Refresh")
        refresh_button.clicked.connect(populate)
        button_row.addWidget(refresh_button)

        export_button = QPushButton("Export as CSV…")
        export_button.clicked.connect(lambda: self._export_chart_csv(dialog, current_popup_model()))
        button_row.addWidget(export_button)

        export_all_button = QPushButton("Export as CSV… (all models)")
        export_all_button.setToolTip(
            "Every residue of every open model in one CSV, using the fixed known-tool "
            "columns only (not a model's own extra discovered attributes, which can "
            "differ model to model) - so every row lines up under the same columns."
        )
        export_all_button.clicked.connect(lambda: self._export_combined_csv(dialog))
        button_row.addWidget(export_all_button)

        ai_export_button = QPushButton("Export for AI analysis…")
        ai_export_button.setToolTip(
            "Save a Markdown file with column explanations and your notes up front - "
            "ready to paste into Claude/ChatGPT/Gemini and ask it to analyze."
        )
        ai_export_button.clicked.connect(lambda: self._export_for_ai_analysis(dialog, current_popup_model()))
        button_row.addWidget(ai_export_button)
        layout.addLayout(button_row)

        dialog.setLayout(layout)
        dialog.show()
        track_dialog(self._open_dialogs, dialog)

    def _build_combined_chart_rows(self):
        """Every residue of every open model, as (headers, rows) - one combined
        table across models, prefixed with a Model column. Uses only the fixed
        _KNOWN_ATTRS columns, not each model's own extra discovered attributes
        (those can differ from model to model, which would misalign columns
        across rows from different models) - for full per-model detail
        including any extra attributes, export a single model's Chart instead.
        Uses _merge_stored_and_live() for each model, same as
        _build_chart_rows - including its stale-value cleanup."""
        headers = ["Model", "Chain", "Residue #", "Name", "pLDDT"] + list(_KNOWN_ATTRS.values()) + ["Notes"]
        rows = []
        for model in self.session.models.list():
            if not hasattr(model, "residues"):
                continue
            store = annotations.get_store(self.session, model)
            model_label = f"#{model.id_string} {model.name}"
            merged_by_residue = self._merge_stored_and_live(model, store)
            for r in sorted(model.residues, key=lambda r: (r.chain_id, r.number)):
                note = store.note_for(r.chain_id, r.number)
                merged = merged_by_residue[r]
                atom = r.principal_atom
                plddt = f"{atom.bfactor:.1f}" if atom is not None else ""
                row = [model_label, r.chain_id, str(r.number), r.name, plddt]
                row.extend(str(merged.get(name, "")) for name in _KNOWN_ATTRS)
                row.append(note["text"] if note else "")
                rows.append(row)
        return headers, rows

    def _export_combined_csv(self, dialog):
        headers, rows = self._build_combined_chart_rows()
        if not rows:
            show_error(dialog, "ChopChopMF", "No open models with residues found.")
            return
        export_dir = get_settings(self.session).export_dir
        file_path, _ = QFileDialog.getSaveFileName(
            dialog, "Save Combined Table As", str(Path(export_dir) / "investigate_chart_all_models.csv"), "CSV Files (*.csv)"
        )
        if not file_path:
            return
        try:
            with open(file_path, "w", newline="", encoding="utf-8-sig") as f:
                writer = csv.writer(f)
                writer.writerow(headers)
                writer.writerows(rows)
        except OSError as e:
            show_error(dialog, "ChopChopMF", f"Failed to save table:\n{e}")

    def _export_chart_csv(self, dialog, model=None):
        headers, rows, _row_keys = self._build_chart_rows(model)
        export_dir = get_settings(self.session).export_dir
        file_path, _ = QFileDialog.getSaveFileName(
            dialog, "Save Table As", str(Path(export_dir) / "investigate_chart.csv"), "CSV Files (*.csv)"
        )
        if not file_path:
            return
        try:
            with open(file_path, "w", newline="", encoding="utf-8-sig") as f:
                writer = csv.writer(f)
                writer.writerow(headers)
                writer.writerows(rows)
        except OSError as e:
            show_error(dialog, "ChopChopMF", f"Failed to save table:\n{e}")

    def _build_ai_export_text(self, question="", model=None):
        """Markdown text for the 'Export for AI analysis' button: the user's
        analysis question up front (asked for at export time - a generic
        review prompt if left blank), a short preamble explaining what
        ChopChopMF is and what each present column means, the user's own
        per-residue notes (usually extra context for that question), then the
        full residue table - as a Markdown table if it's short enough to stay
        readable, otherwise as a fenced CSV block."""
        if model is None:
            model = self._current_model()
        if model is None:
            return None
        headers, rows, row_keys = self._build_chart_rows(model)

        ran_tools = [
            h for idx, h in enumerate(headers)
            if h in _KNOWN_ATTRS.values() and any(row[idx] for row in rows)
        ]

        lines = [f"# ChopChopMF residue analysis export - {model.name}", ""]
        lines.append("## What to analyze")
        lines.append("")
        if question.strip():
            lines.append(question.strip())
        else:
            lines.append(
                "No specific question was given - do a general review instead: point "
                "out residues that look unusual or noteworthy given the data below "
                "(e.g. low confidence, high pathogenicity, or several tools agreeing "
                "on the same residue)."
            )
        lines.append("")

        lines.append("## References - please read these before analyzing, if you can browse")
        lines.append("")
        lines.append(
            "These are the canonical background sources for the data and terms used in "
            "this file. If you have web access, actually open and read them (not just "
            "list them back) before answering:"
        )
        lines.append("")
        lines.append(
            f"- ChopChopMF usage guide - {GUIDE_BASE_URL} - how each ChopChopMF tool/column "
            "mentioned below actually works."
        )
        lines.append(
            "- AlphaFold Guide - https://lukasinscience.github.io/AlphaFold-Guide/ - background on "
            "what pLDDT, PAE, ipTM, and AlphaMissense scores mean and how to interpret them."
        )
        lines.append(
            "- ChimeraX FigureStyle guide - https://lukasinscience.github.io/ChimeraX-FigureStyle/ - "
            "conventions for rendering/coloring structures in ChimeraX, useful if the analysis "
            "touches how a residue or region is visualized."
        )
        lines.append("")

        lines.append(
            "ChopChopMF is a ChimeraX plug-in for protein structure analysis. "
            "This file lists every residue of the structure below, together with "
            "whatever ChopChopMF analysis tools have been run on it so far."
        )
        lines.append("")
        if ran_tools:
            lines.append(f"Analyses already run on this structure: {', '.join(ran_tools)}.")
        else:
            lines.append("No ChopChopMF analysis tool has been run on this structure yet - only pLDDT is available.")
        lines.append("")

        lines.append("## Column meanings")
        lines.append("")
        for header in headers:
            desc = _COLUMN_DESCRIPTIONS.get(header)
            if desc:
                lines.append(f"- **{header}**: {desc}")
        lines.append("")

        notes = [(row_keys[i], row[-1]) for i, row in enumerate(rows) if row[-1]]
        if notes:
            lines.append("## Notes (likely the most important context)")
            lines.append("")
            for (chain_id, number), text in notes:
                lines.append(f"- Chain {chain_id}, residue {number}: {text}")
            lines.append("")

        lines.append(f"## Full residue table ({len(rows)} residues)")
        lines.append("")
        if len(rows) <= 50:
            lines.append("| " + " | ".join(headers) + " |")
            lines.append("|" + "|".join(["---"] * len(headers)) + "|")
            for row in rows:
                cells = [(c.replace("|", "/") or "-") for c in row]
                lines.append("| " + " | ".join(cells) + " |")
        else:
            lines.append(f"({len(rows)} rows is too many for a readable Markdown table, given as CSV instead)")
            lines.append("")
            lines.append("```csv")
            buf = io.StringIO()
            writer = csv.writer(buf)
            writer.writerow(headers)
            writer.writerows(rows)
            lines.append(buf.getvalue().rstrip("\n"))
            lines.append("```")
        lines.append("")
        return "\n".join(lines)

    def _export_for_ai_analysis(self, dialog, model=None):
        if model is None:
            model = self._current_model()
        if model is None:
            show_error(dialog, "ChopChopMF", "Select a model first.")
            return
        question, ok = QInputDialog.getMultiLineText(
            dialog, "Export for AI analysis",
            "What do you want the AI to look at or answer? (optional - written into "
            "the file itself so you don't have to explain it again when you paste it in)",
            ""
        )
        if not ok:
            return
        text = self._build_ai_export_text(question, model)
        export_dir = get_settings(self.session).export_dir
        file_path, _ = QFileDialog.getSaveFileName(
            dialog, "Save AI Analysis Export", str(Path(export_dir) / f"{model.name}_ai_export.md"), "Markdown Files (*.md)"
        )
        if not file_path:
            return
        try:
            Path(file_path).write_text(text, encoding="utf-8")
        except OSError as e:
            show_error(dialog, "ChopChopMF", f"Failed to save export:\n{e}")

    # ---------- Plots tab ----------

    def _create_plots_tab(self):
        widget = QWidget()
        layout = QVBoxLayout()
        layout.setAlignment(Qt.AlignTop)
        layout.setSpacing(4)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.addWidget(make_guide_button("4-analyze-structure"))

        self._add_plot_section(
            layout, "Attribute Coverage", self._draw_coverage_heatmap,
            "Every residue vs. every value any ChopChopMF tool can produce - gray means "
            "that tool hasn't recorded anything for this residue yet. Good for spotting "
            "residues with evidence from several tools at once."
        )

        scatter_group = QGroupBox("Compare Two Values")
        scatter_layout = QVBoxLayout()
        scatter_desc = QLabel(
            "Pick any two values with numeric data to plot one dot per residue - "
            "e.g. pLDDT vs. AlphaMissense score. Click a point in the opened figure "
            "to select that residue (jumps straight to its Residue-tab dossier)."
        )
        scatter_desc.setWordWrap(True)
        scatter_layout.addWidget(scatter_desc)
        axis_row = QHBoxLayout()
        axis_row.addWidget(QLabel("X:"))
        self.scatter_x = QComboBox()
        axis_row.addWidget(self.scatter_x, stretch=1)
        axis_row.addWidget(QLabel("Y:"))
        self.scatter_y = QComboBox()
        axis_row.addWidget(self.scatter_y, stretch=1)
        scatter_layout.addLayout(axis_row)
        scatter_button = QPushButton("Open Figure")
        scatter_button.clicked.connect(self._open_scatter_figure)
        scatter_layout.addWidget(scatter_button)
        scatter_group.setLayout(scatter_layout)
        layout.addWidget(scatter_group)

        widget.setLayout(layout)
        self._refresh_plot_axis_choices()
        return widget

    def _add_plot_section(self, layout, heading, draw_func, description):
        group = QGroupBox(heading)
        group_layout = QVBoxLayout()
        desc_label = QLabel(description)
        desc_label.setWordWrap(True)
        group_layout.addWidget(desc_label)
        open_button = QPushButton("Open Figure")
        open_button.clicked.connect(lambda: self._open_figure(draw_func, heading))
        group_layout.addWidget(open_button)
        group.setLayout(group_layout)
        layout.addWidget(group)

    def _open_figure(self, draw_func, title):
        dialog = open_figure_dialog(self.tool_window.ui_area, draw_func, title)
        track_dialog(self._open_dialogs, dialog)

    def _open_scatter_figure(self):
        dialog = open_figure_dialog(
            self.tool_window.ui_area, self._draw_scatter_figure, "Compare Values",
            on_pick=self._on_scatter_pick
        )
        track_dialog(self._open_dialogs, dialog)

    def _on_scatter_pick(self, event):
        residues = getattr(self, "_scatter_residues", None)
        if not residues or not event.ind:
            return
        residue = residues[event.ind[0]]
        # Select by direct atom flags (not an atomspec string) - unambiguous even
        # with several open models that happen to share a chain letter/number.
        run(self.session, "select clear")
        for atom in residue.atoms:
            atom.selected = True

    def _no_data_message(self, figure, message):
        ax = figure.add_subplot(111)
        ax.text(0.5, 0.5, message, ha="center", va="center", wrap=True)
        ax.axis("off")

    def _numeric_columns(self, model):
        """(label, extractor) pairs for every value with at least one numeric
        reading on this model - pLDDT plus any known/discovered residue
        attribute (live or durably recorded) that parses as a float for at
        least one residue. Feeds the scatter plot's X/Y axis pickers."""
        residues = model.residues
        store = annotations.get_store(self.session, model)

        def plddt(r):
            atom = r.principal_atom
            return atom.bfactor if atom is not None else None

        columns = [("pLDDT", plddt)]
        attr_names, attr_labels = _attr_columns(model, store)
        for name, label in zip(attr_names, attr_labels):
            def extractor(r, name=name):
                merged = {**store.values_for(r.chain_id, r.number), **dict(r.custom_attrs)}
                try:
                    return float(merged.get(name))
                except (TypeError, ValueError):
                    return None
            if any(extractor(r) is not None for r in residues):
                columns.append((label, extractor))
        return columns

    def _refresh_plot_axis_choices(self):
        if not hasattr(self, "scatter_x"):
            return  # Plots tab not built yet
        model = self._current_model()
        self.scatter_x.blockSignals(True)
        self.scatter_y.blockSignals(True)
        self.scatter_x.clear()
        self.scatter_y.clear()
        if model is not None:
            labels = [label for label, _ in self._numeric_columns(model)]
            self.scatter_x.addItems(labels)
            self.scatter_y.addItems(labels)
            if len(labels) > 1:
                self.scatter_y.setCurrentIndex(1)
        self.scatter_x.blockSignals(False)
        self.scatter_y.blockSignals(False)

    def _draw_coverage_heatmap(self, figure):
        model = self._current_model()
        if model is None:
            self._no_data_message(figure, "No model selected.")
            return
        residues = sorted(model.residues, key=lambda r: (r.chain_id, r.number))
        if not residues:
            self._no_data_message(figure, "No residues on this model.")
            return
        store = annotations.get_store(self.session, model)
        attr_names, attr_labels = _attr_columns(model, store)
        columns = ["pLDDT"] + attr_labels
        data = np.full((len(residues), len(columns)), np.nan)
        for i, r in enumerate(residues):
            atom = r.principal_atom
            if atom is not None:
                data[i, 0] = atom.bfactor / 100.0
            attrs = {**store.values_for(r.chain_id, r.number), **dict(r.custom_attrs)}
            for j, name in enumerate(attr_names, start=1):
                if name in attrs:
                    try:
                        data[i, j] = float(attrs[name])
                    except (TypeError, ValueError):
                        data[i, j] = 1.0  # non-numeric but present -> "has a value"

        for j in range(1, data.shape[1]):
            col = data[:, j]
            finite = col[~np.isnan(col)]
            if finite.size == 0:
                continue
            lo, hi = finite.min(), finite.max()
            if hi > lo:
                data[:, j] = (col - lo) / (hi - lo)
            else:
                data[:, j] = np.where(np.isnan(col), np.nan, 0.5)

        ax = figure.add_subplot(111)
        im = ax.imshow(data, aspect="auto", cmap=_COVERAGE_CMAP, vmin=0, vmax=1, interpolation="nearest")
        ax.set_xticks(range(len(columns)))
        ax.set_xticklabels(columns, rotation=45, ha="right", fontsize=8)
        step = max(1, len(residues) // 25)
        yticks = list(range(0, len(residues), step))
        ax.set_yticks(yticks)
        ax.set_yticklabels([f"{residues[i].chain_id}{residues[i].number}" for i in yticks], fontsize=7)
        ax.set_title("Attribute coverage per residue (gray = no value yet)")
        figure.colorbar(im, ax=ax, fraction=0.03, pad=0.04, label="normalized value / present")
        figure.tight_layout()

    def _draw_scatter_figure(self, figure):
        model = self._current_model()
        if model is None:
            self._no_data_message(figure, "No model selected.")
            return
        columns = dict(self._numeric_columns(model))
        x_label, y_label = self.scatter_x.currentText(), self.scatter_y.currentText()
        if x_label not in columns or y_label not in columns:
            self._no_data_message(figure, "Pick two values with numeric data first.")
            return
        x_extract, y_extract = columns[x_label], columns[y_label]
        xs, ys, residues = [], [], []
        for r in model.residues:
            x, y = x_extract(r), y_extract(r)
            if x is not None and y is not None:
                xs.append(x)
                ys.append(y)
                residues.append(r)
        if not xs:
            self._no_data_message(figure, "No residue has both values yet.")
            return
        self._scatter_residues = residues
        ax = figure.add_subplot(111)
        ax.scatter(xs, ys, s=18, color="#00695C", alpha=0.7, edgecolors="none", picker=5)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(f"{y_label} vs. {x_label} ({len(xs)} residues) - click a point to select it")
