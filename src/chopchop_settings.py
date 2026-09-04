#!/usr/bin/env python3

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
Settings Tool
Central place for every persistent location ChopChopMF tools save output to:
the shared download folder (AlphaMissense fetch/Sequence/ChopMissense/PDBePISA
XML lookups), the shared export folder (suggested starting folder for CSV/
Markdown "Save As" dialogs across PAE Analysis/Batch Analysis/Investigate),
PDBePISA's .defattr output folder, and - per open model - Investigate's
durable `.chopchop.json` annotations file, including saving a timestamped
snapshot ("Save Session As...") and loading an earlier one back ("Change...").
Nothing computed here; this tool only points other tools at where to read/write.
"""

from datetime import datetime
from pathlib import Path
import shutil

from chimerax.core.tools import ToolInstance
from chimerax.ui import MainToolWindow
from Qt.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QComboBox, QLineEdit,
    QFileDialog, QGroupBox, QWidget,
)

from .utils import make_guide_button, make_scrollable, show_error, get_settings
from . import annotations


class ChopChopSettings(ToolInstance):

    SESSION_ENDURING = False
    SESSION_SAVE = False

    def __init__(self, session, tool_name):
        super().__init__(session, tool_name)
        self.display_name = "Setup"
        self.tool_window = MainToolWindow(self)
        self._build_ui()
        self.tool_window.manage("side")

    def _build_ui(self):
        layout = QVBoxLayout()
        layout.setSpacing(6)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.addWidget(make_guide_button("7-setup"))

        locations_group = QGroupBox("Shared Locations")
        locations_layout = QVBoxLayout()
        self._download_field = self._add_folder_setting_row(
            locations_layout, "Download folder", "download_dir",
            "Used by AlphaMissense fetch, Sequence, ChopMissense, and PDBePISA lookups - "
            "changing it here is the same as changing it in any of those tools' own fields."
        )
        self._export_field = self._add_folder_setting_row(
            locations_layout, "Export folder (CSV/Markdown suggestions)", "export_dir",
            "Suggested starting folder for PAE Analysis/Batch Analysis/Investigate's "
            "CSV and Markdown \"Save As\" dialogs - they still ask every time, this only "
            "changes where the dialog starts."
        )
        self._add_defattr_row(locations_layout)
        locations_group.setLayout(locations_layout)
        layout.addWidget(locations_group)

        annotations_group = QGroupBox("Model Annotations (Investigate Sessions)")
        annotations_layout = QVBoxLayout()
        model_row = QHBoxLayout()
        model_row.addWidget(QLabel("Model:"))
        self.model_selector = QComboBox()
        self.model_selector.currentIndexChanged.connect(self._refresh_annotations_path_label)
        model_row.addWidget(self.model_selector, stretch=1)
        refresh_models_button = QPushButton("↻")
        refresh_models_button.setToolTip("Refresh model list")
        refresh_models_button.clicked.connect(self._refresh_models)
        model_row.addWidget(refresh_models_button)
        annotations_layout.addLayout(model_row)

        self.annotations_path_label = QLabel("")
        self.annotations_path_label.setStyleSheet("color: gray; font-size: 10px;")
        self.annotations_path_label.setWordWrap(True)
        annotations_layout.addWidget(self.annotations_path_label)

        button_row = QHBoxLayout()
        save_session_button = QPushButton("Save Session As…")
        save_session_button.setToolTip(
            "Snapshot this model's current .chopchop.json to a new, timestamped file - the "
            "live file keeps being used/updated as normal, so this is purely a copy to come "
            "back to later without it ever getting silently overwritten by further work."
        )
        save_session_button.clicked.connect(self._save_annotations_copy)
        button_row.addWidget(save_session_button)
        change_path_button = QPushButton("Change…")
        change_path_button.setToolTip(
            "Point this model's annotations at a different file - pick an existing "
            "*.chopchop.json (e.g. one saved earlier with 'Save Session As…') to load that "
            "session back and keep working in it, or a new filename to start fresh without "
            "touching the current file."
        )
        change_path_button.clicked.connect(self._choose_annotations_path)
        button_row.addWidget(change_path_button)
        annotations_layout.addLayout(button_row)
        annotations_group.setLayout(annotations_layout)
        layout.addWidget(annotations_group)

        layout.addStretch(1)
        container = QWidget()
        container.setLayout(layout)
        outer_layout = QVBoxLayout()
        outer_layout.setContentsMargins(0, 0, 0, 0)
        outer_layout.addWidget(make_scrollable(container))
        self.tool_window.ui_area.setLayout(outer_layout)
        self._refresh_models()

    # ---------- Shared Locations ----------

    def _add_folder_setting_row(self, layout, label_text, attr_name, tooltip):
        """One row: <label> [path field] [Browse...] - reads/writes
        get_settings(session).<attr_name> directly. Typed edits commit on
        losing focus/Enter (QLineEdit.editingFinished), not per keystroke;
        picking via Browse commits immediately, since that's an unambiguous
        choice. Used for the two plain shared folders (download_dir,
        export_dir) - defattr_output_dir has its own row (see
        _add_defattr_row) since an empty value has a different, non-folder
        meaning there."""
        layout.addWidget(QLabel(f"{label_text}:"))
        row = QHBoxLayout()
        field = QLineEdit(getattr(get_settings(self.session), attr_name))
        field.setToolTip(tooltip)
        row.addWidget(field, stretch=1)

        def commit():
            value = field.text().strip()
            if value:
                setattr(get_settings(self.session), attr_name, value)
            else:
                field.setText(getattr(get_settings(self.session), attr_name))

        field.editingFinished.connect(commit)

        def browse():
            folder = QFileDialog.getExistingDirectory(
                self.tool_window.ui_area, f"Choose {label_text}", field.text()
            )
            if folder:
                field.setText(folder)
                setattr(get_settings(self.session), attr_name, folder)

        browse_button = QPushButton("Browse…")
        browse_button.clicked.connect(browse)
        row.addWidget(browse_button)
        layout.addLayout(row)
        return field

    def _add_defattr_row(self, layout):
        layout.addWidget(QLabel("PDBePISA .defattr output folder:"))
        row = QHBoxLayout()
        self._defattr_field = QLineEdit(get_settings(self.session).defattr_output_dir)
        self._defattr_field.setPlaceholderText("(default: next to the loaded PISA XML file)")
        self._defattr_field.setToolTip(
            "Where PDBePISA writes both its interface-class and ΔG-coloring .defattr files. "
            "Empty means the default: next to whichever PISA XML file was loaded."
        )
        row.addWidget(self._defattr_field, stretch=1)

        def commit():
            get_settings(self.session).defattr_output_dir = self._defattr_field.text().strip()

        self._defattr_field.editingFinished.connect(commit)

        def browse():
            folder = QFileDialog.getExistingDirectory(
                self.tool_window.ui_area, "Choose .defattr Output Folder", self._defattr_field.text()
            )
            if folder:
                self._defattr_field.setText(folder)
                get_settings(self.session).defattr_output_dir = folder

        browse_button = QPushButton("Browse…")
        browse_button.clicked.connect(browse)
        row.addWidget(browse_button)
        layout.addLayout(row)

        reset_button = QPushButton("Reset to default (next to XML)")
        reset_button.clicked.connect(self._reset_defattr)
        layout.addWidget(reset_button)

    def _reset_defattr(self):
        self._defattr_field.setText("")
        get_settings(self.session).defattr_output_dir = ""

    # ---------- Model Annotations (moved here from Investigate) ----------

    def _refresh_models(self):
        self.model_selector.clear()
        for m in self.session.models.list():
            if hasattr(m, "residues"):
                self.model_selector.addItem(f"#{m.id_string} {m.name}", m.id_string)
        self._refresh_annotations_path_label()

    def _current_model(self):
        id_string = self.model_selector.currentData()
        if not id_string:
            return None
        return next((m for m in self.session.models.list() if m.id_string == id_string), None)

    def _refresh_annotations_path_label(self):
        model = self._current_model()
        if model is None:
            self.annotations_path_label.setText("No model selected.")
            return
        store = annotations.get_store(self.session, model)
        self.annotations_path_label.setText(f"Annotations file: {store.path}")

    def _choose_annotations_path(self):
        model = self._current_model()
        if model is None:
            show_error(self.tool_window.ui_area, "ChopChopMF", "Select a model first.")
            return
        current_path = annotations.get_store(self.session, model).path
        file_path, _ = QFileDialog.getSaveFileName(
            self.tool_window.ui_area, "Choose Annotations File", str(current_path), "JSON Files (*.json)"
        )
        if not file_path:
            return
        annotations.set_store_path(self.session, model, Path(file_path))
        self._refresh_annotations_path_label()

    def _save_annotations_copy(self):
        """Snapshot the selected model's current annotations file to a new
        path the user picks - a "Save Session As" that never touches the
        live file (which keeps accumulating changes as normal), so a
        snapshot from today can't be silently overwritten by tomorrow's
        session. Suggests a timestamped name right next to the live file,
        so old snapshots are easy to spot and don't collide with each other
        or the live file."""
        model = self._current_model()
        if model is None:
            show_error(self.tool_window.ui_area, "ChopChopMF", "Select a model first.")
            return
        store = annotations.get_store(self.session, model)
        if not store.path.exists():
            show_error(
                self.tool_window.ui_area, "ChopChopMF",
                "No annotations file exists for this model yet (it's created on the first note "
                "or the first tool-computed value)."
            )
            return
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        suggested_name = f"{store.path.stem}_session_{timestamp}{store.path.suffix}"
        suggested_path = store.path.parent / suggested_name
        file_path, _ = QFileDialog.getSaveFileName(
            self.tool_window.ui_area, "Save Session As", str(suggested_path), "JSON Files (*.json)"
        )
        if not file_path:
            return
        try:
            shutil.copy(store.path, file_path)
        except OSError as e:
            show_error(self.tool_window.ui_area, "ChopChopMF", f"Failed to save session:\n{e}")
