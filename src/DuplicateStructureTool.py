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
DuplicateStructureTool.py

ChimeraX GUI tool providing three tabs:
1. Duplicate Structure
2. Measure Center
3. Symmetry Copies
"""

from chimerax.core.tools import ToolInstance
from chimerax.core.commands import run
from chimerax.ui import MainToolWindow
from Qt.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QLabel, QComboBox, QPushButton, QWidget,
    QCheckBox, QDoubleSpinBox, QGroupBox, QTabWidget, QLineEdit
)
from Qt.QtGui import QFont
from chimerax.map import Volume

ALL_CHAINS_LABEL = "*(All chains)"
BLANK_CHAIN_LABEL = "(blank)"


class DuplicateStructureTool(ToolInstance):
    SESSION_ENDURING = False
    SESSION_SAVE = True
    help = "help:user/tools/duplicatestructure.html"

    def __init__(self, session, tool_name):
        super().__init__(session, tool_name)
        self.display_name = "Duplicate Structure"
        self.tool_window = MainToolWindow(self)
        self._build_ui()
        self.tool_window.manage('side')

    def _build_ui(self):
        self.tabs = QTabWidget()
        self.tabs.addTab(self._build_duplicate_tab(), "Duplicate Structure")
        self.tabs.addTab(self._build_measure_center_tab(), "Measure Center")
        self.tabs.addTab(self._build_symmetry_tab(), "Symmetry Copies")

        layout = QVBoxLayout()
        layout.addWidget(self.tabs)
        self.tool_window.ui_area.setLayout(layout)

    # -------------------------
    # Duplicate Structure Tab
    # -------------------------
    def _build_duplicate_tab(self):
        layout = QVBoxLayout()
        title = QLabel("Duplicate Structure Tool")
        tf = QFont(); tf.setPointSize(14); tf.setBold(True)
        title.setFont(tf)
        layout.addWidget(title)

        layout.addWidget(QLabel("Select model:"))
        self.model_selector = QComboBox()
        layout.addWidget(self.model_selector)

        layout.addWidget(QLabel("Select chain:"))
        self.chain_selector = QComboBox()
        layout.addWidget(self.chain_selector)

        refresh_btn = QPushButton("↻ Refresh model list & chains")
        refresh_btn.clicked.connect(self._refresh_models_and_chains)
        layout.addWidget(refresh_btn)

        offset_group = QGroupBox("Offset duplicate")
        offset_layout = QVBoxLayout()
        self.offset_enable = QCheckBox("Apply offset to duplicate")
        self.offset_enable.stateChanged.connect(self._toggle_offset_enabled)
        offset_layout.addWidget(self.offset_enable)

        vec_row = QHBoxLayout()
        vec_row.addWidget(QLabel("ΔX (Å):"))
        self.dx = QDoubleSpinBox(); self._setup_spin(self.dx, 10.0); vec_row.addWidget(self.dx)
        vec_row.addWidget(QLabel("ΔY (Å):"))
        self.dy = QDoubleSpinBox(); self._setup_spin(self.dy, 0.0); vec_row.addWidget(self.dy)
        vec_row.addWidget(QLabel("ΔZ (Å):"))
        self.dz = QDoubleSpinBox(); self._setup_spin(self.dz, 0.0); vec_row.addWidget(self.dz)

        offset_layout.addLayout(vec_row)
        offset_group.setLayout(offset_layout)
        layout.addWidget(offset_group)

        btn_row = QHBoxLayout()
        duplicate_btn = QPushButton("ChopChop Double")
        duplicate_btn.clicked.connect(self.duplicate_structure)
        btn_row.addWidget(duplicate_btn)
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.delete)
        btn_row.addWidget(close_btn)
        layout.addLayout(btn_row)

        widget = QWidget()
        widget.setLayout(layout)
        self._refresh_models_and_chains()
        self.model_selector.currentIndexChanged.connect(self._refresh_chain_list)
        self._toggle_offset_enabled()
        return widget

    # -------------------------
    # Measure Center Tab
    # -------------------------
    def _build_measure_center_tab(self):
        layout = QVBoxLayout()
        title = QLabel("Measure Center of Map")
        tf = QFont(); tf.setPointSize(14); tf.setBold(True)
        title.setFont(tf)
        layout.addWidget(title)

        layout.addWidget(QLabel("Select map (volume model):"))
        self.map_selector = QComboBox()
        layout.addWidget(self.map_selector)

        refresh_btn = QPushButton("↻ Refresh maps")
        refresh_btn.clicked.connect(self._refresh_maps)
        layout.addWidget(refresh_btn)

        measure_btn = QPushButton("ChopChop Measure Center")
        measure_btn.clicked.connect(self.measure_center)
        layout.addWidget(measure_btn)

        layout.addStretch(1)
        widget = QWidget()
        widget.setLayout(layout)
        self._refresh_maps()
        return widget

    # -------------------------
    # Symmetry Copies Tab
    # -------------------------
    def _build_symmetry_tab(self):
        layout = QVBoxLayout()
        title = QLabel("Symmetry Copies")
        tf = QFont(); tf.setPointSize(14); tf.setBold(True)
        title.setFont(tf)
        layout.addWidget(title)

        layout.addWidget(QLabel("Select structure model:"))
        self.sym_model_selector = QComboBox()
        layout.addWidget(self.sym_model_selector)

        layout.addWidget(QLabel("Symmetry group (e.g. C2, C3, D2):"))
        self.sym_group_box = QComboBox()
        for sym in ["C2","C3","C4","C5","C6","D2","D3","D4","T","O","I"]:
            self.sym_group_box.addItem(sym)
        layout.addWidget(self.sym_group_box)

        paste_layout = QHBoxLayout()
        paste_layout.addWidget(QLabel("Paste XYZ for center:"))
        self.xyz_paste = QLineEdit()
        self.xyz_paste.setPlaceholderText("353.79, 353.79, 333.95")
        paste_layout.addWidget(self.xyz_paste)
        paste_btn = QPushButton("Apply")
        paste_btn.clicked.connect(self._apply_xyz_paste)
        paste_layout.addWidget(paste_btn)
        layout.addLayout(paste_layout)

        coord_layout = QHBoxLayout()
        coord_layout.addWidget(QLabel("Center X:"))
        self.cx = QDoubleSpinBox(); self._setup_spin(self.cx, 0); coord_layout.addWidget(self.cx)
        coord_layout.addWidget(QLabel("Y:"))
        self.cy = QDoubleSpinBox(); self._setup_spin(self.cy, 0); coord_layout.addWidget(self.cy)
        coord_layout.addWidget(QLabel("Z:"))
        self.cz = QDoubleSpinBox(); self._setup_spin(self.cz, 0); coord_layout.addWidget(self.cz)
        layout.addLayout(coord_layout)

        refresh_btn = QPushButton("↻ Refresh models")
        refresh_btn.clicked.connect(self._refresh_sym_models)
        layout.addWidget(refresh_btn)

        sym_btn = QPushButton("ChopChop Symmetry Copies")
        sym_btn.clicked.connect(self.make_symmetry_copies)
        layout.addWidget(sym_btn)

        self.sym_status = QLabel("")
        layout.addWidget(self.sym_status)
        layout.addStretch(1)

        widget = QWidget()
        widget.setLayout(layout)
        self._refresh_sym_models()
        return widget

    # -------------------------
    # Utility Functions
    # -------------------------
    def _setup_spin(self, spin, val):
        spin.setRange(-1e6, 1e6)
        spin.setDecimals(3)
        spin.setSingleStep(1.0)
        spin.setValue(float(val))

    def _toggle_offset_enabled(self):
        enabled = self.offset_enable.isChecked()
        self.dx.setEnabled(enabled)
        self.dy.setEnabled(enabled)
        self.dz.setEnabled(enabled)

    def _refresh_models_and_chains(self):
        self._populate_models()
        self._refresh_chain_list()

    def _populate_models(self):
        self.model_selector.clear()
        for m in self.session.models.list():
            if hasattr(m, "residues"):
                self.model_selector.addItem(m.id_string)

    def _refresh_chain_list(self):
        self.chain_selector.clear()
        model = self._model_from_id(self.model_selector.currentText().strip())
        if model is None:
            return
        chains = sorted(set(r.chain_id for r in model.residues))
        self.chain_selector.addItem(ALL_CHAINS_LABEL)
        for c in chains:
            self.chain_selector.addItem(c if c else BLANK_CHAIN_LABEL)

    def _model_from_id(self, id_str):
        for m in self.session.models.list():
            if getattr(m, "id_string", None) == id_str:
                return m
        return None

    # -------------------------
    # Duplicate Structure Logic
    # -------------------------
    def duplicate_structure(self):
        model_id = self.model_selector.currentText().strip()
        chain_choice = self.chain_selector.currentText().strip()
        if not model_id:
            print("Please select a model.")
            return

        specific_chain = chain_choice not in ("", ALL_CHAINS_LABEL)
        if specific_chain:
            chain_id = "" if chain_choice == BLANK_CHAIN_LABEL else chain_choice
            what = f"Model #{model_id}, Chain '{chain_id}'"
        else:
            what = f"Model #{model_id} (all chains)"

        before = {m.id_string for m in self.session.models.list() if hasattr(m, "residues")}
        run(self.session, f"combine #{model_id}")
        after = {m.id_string for m in self.session.models.list() if hasattr(m, "residues")}
        new_ids = sorted(after - before, key=lambda x: [int(y) for y in x.split(".") if y.isdigit()])
        if not new_ids:
            return

        created = ", ".join(f"#{nid}" for nid in new_ids)
        print(f"Duplicated {what} -> {created}")

        if specific_chain:
            for nid in new_ids:
                if chain_choice == BLANK_CHAIN_LABEL:
                    self._delete_non_blank_chains_in_model(nid)
                else:
                    run(self.session, f"delete #{nid} & ~ #{nid}/{chain_id}")

        if self.offset_enable.isChecked():
            dx, dy, dz = self.dx.value(), self.dy.value(), self.dz.value()
            id_list = ",".join("#" + i for i in new_ids)
            run(self.session,
                f"move x {dx} models {id_list} ; move y {dy} models {id_list} ; move z {dz} models {id_list}")

    def _delete_non_blank_chains_in_model(self, model_id):
        model = self._model_from_id(model_id)
        if not model:
            return
        chains = [c for c in set(r.chain_id for r in model.residues) if c]
        if chains:
            run(self.session, "delete " + " ".join(f"#{model_id}/{c}" for c in chains))
        run(self.session, "select clear")

    # -------------------------
    # Measure Center Logic
    # -------------------------
    def measure_center(self):
        map_id = self.map_selector.currentData()
        if not map_id or map_id.startswith("("):
            msg = "Please select a valid density map."
            print(msg)
            self.session.logger.warning(msg)
            return

        cmd = f"measure center #{map_id}"
        print(f"Executing: {cmd}")
        run(self.session, cmd)

    def _refresh_maps(self):
        self.map_selector.clear()
        for m in self.session.models.list():
            if isinstance(m, Volume):
                self.map_selector.addItem(f"{m.id_string} ({m.name})", m.id_string)
        if self.map_selector.count() == 0:
            self.map_selector.addItem("(No maps loaded)")

    # -------------------------
    # Symmetry Copies Logic
    # -------------------------
    def make_symmetry_copies(self):
        model_id = self.sym_model_selector.currentText().strip()
        if not model_id:
            return
        sym = self.sym_group_box.currentText().strip()
        x, y, z = self.cx.value(), self.cy.value(), self.cz.value()
        cmd = f"sym #{model_id} {sym} center {x:.2f},{y:.2f},{z:.2f}"
        print(f"Executing: {cmd}")
        run(self.session, cmd)
        self.sym_status.setText(f"Created symmetry copies for #{model_id} ({sym})")

    def _apply_xyz_paste(self):
        text = self.xyz_paste.text().strip().strip("()")
        try:
            vals = [float(x) for x in text.replace(",", " ").split() if x]
            if len(vals) != 3:
                raise ValueError
        except Exception:
            self.sym_status.setText(" Could not parse XYZ input.")
            return
        x, y, z = vals
        self.cx.setValue(x); self.cy.setValue(y); self.cz.setValue(z)
        self.sym_status.setText(f"XYZ center set to ({x:.2f}, {y:.2f}, {z:.2f})")

    def _refresh_sym_models(self):
        self.sym_model_selector.clear()
        for m in self.session.models.list():
            if hasattr(m, "residues"):
                self.sym_model_selector.addItem(m.id_string)
        if self.sym_model_selector.count() == 0:
            self.sym_model_selector.addItem("(No models loaded)")


def show_tool(session):
    return DuplicateStructureTool(session, "duplicate structure")
