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
Transmembrane Helix Tool
Runs DSSP, groups consecutive alpha-helices >= a minimum length, and flags
those whose mean Kyte & Doolittle (1982) hydrophobicity is above a threshold
as candidate transmembrane helices. When ChimeraX's own `mlp` lipophilicity
table is available, also reports a SASA-weighted Fauchere & Pliska (1983)
lipophilicity score per helix (the rbvi/chimerax-recipes "helixmlp" method) -
membrane-facing surface should itself be lipophilic, a stronger signal than
residue identity alone. Results are a heuristic motif match, not a
calibrated prediction - see docs/acknowledgements.md for citations.
"""

from chimerax.core.tools import ToolInstance
from chimerax.core.commands import run
from chimerax.atomic import Residue
from chimerax.ui import MainToolWindow
from Qt.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QComboBox,
    QDoubleSpinBox, QSpinBox, QDialog, QTableWidget, QTableWidgetItem, QHeaderView,
    QAbstractItemView,
)
from Qt.QtCore import Qt

from .utils import make_guide_button, busy_cursor, show_error
from . import cellbio


class TransmembraneHelix(ToolInstance):

    SESSION_ENDURING = False
    SESSION_SAVE = True

    def __init__(self, session, tool_name):
        super().__init__(session, tool_name)
        self.display_name = "TM Helix"
        self.tool_window = MainToolWindow(self)
        self._hits = []
        self._build_ui()
        self.tool_window.manage('side')

    def _build_ui(self):
        layout = QVBoxLayout()
        layout.setAlignment(Qt.AlignTop)
        layout.addWidget(make_guide_button("5-cell-biology"))
        description = QLabel("Heuristic prediction only - see guide for methodology/citations.")
        description.setWordWrap(True)
        layout.addWidget(description)

        self.chain_selector = QComboBox()
        layout.addWidget(QLabel("Model : Chain"))
        layout.addWidget(self.chain_selector)
        refresh_button = QPushButton("↻ Refresh model list")
        refresh_button.clicked.connect(self._refresh_chain_list)
        layout.addWidget(refresh_button)

        row = QHBoxLayout()
        row.addWidget(QLabel("Min. helix length"))
        self.length_spin = QSpinBox()
        self.length_spin.setRange(5, 60)
        self.length_spin.setValue(15)
        row.addWidget(self.length_spin)
        row.addWidget(QLabel("Hydrophobicity ≥"))
        self.hydro_spin = QDoubleSpinBox()
        self.hydro_spin.setRange(-5.0, 5.0)
        self.hydro_spin.setValue(1.0)
        row.addWidget(self.hydro_spin)
        layout.addLayout(row)

        self.run_button = QPushButton("ChopChop predict TM Helices")
        self.run_button.clicked.connect(self._run)
        layout.addWidget(self.run_button)

        self.results_button = QPushButton("Show Results…")
        self.results_button.clicked.connect(self._show_results)
        self.results_button.setEnabled(False)
        layout.addWidget(self.results_button)

        self.tool_window.ui_area.setLayout(layout)
        self._refresh_chain_list()

    def _refresh_chain_list(self):
        self.chain_selector.clear()
        for model in self.session.models.list():
            if not hasattr(model, "chains"):
                continue
            for chain in model.chains:
                self.chain_selector.addItem(f"{model.id_string}:{chain.chain_id}")

    def _selected_chain(self):
        text = self.chain_selector.currentText()
        if not text:
            return None
        model_id, chain_id = text.split(":", 1)
        for model in self.session.models.list():
            if model.id_string == model_id and hasattr(model, "chains"):
                for chain in model.chains:
                    if chain.chain_id == chain_id:
                        return chain
        return None

    def _run(self):
        with busy_cursor(self.run_button):
            chain = self._selected_chain()
            if chain is None:
                show_error(self.tool_window.ui_area, "ChopChopMF", "Select a model:chain first.")
                return
            structure_id = chain.structure.id_string
            try:
                self.session.logger.status("TM Helix: running DSSP...")
                run(self.session, f"dssp #{structure_id}")
                self.session.logger.status("TM Helix: measuring SASA...")
                run(self.session, f"measure sasa #{structure_id}")
            except Exception as e:
                show_error(self.tool_window.ui_area, "ChopChopMF", f"DSSP/SASA failed: {e}")
                return
            _sequence, residues = cellbio.sequence_and_residues(chain)
            self._hits = cellbio.scan_tm_helices(
                residues, self.length_spin.value(), self.hydro_spin.value())
            self._write_attributes(residues)

            if self._hits:
                specs = [f"/{h['start_residue'].chain_id}:{h['start_residue'].number}-{h['end_residue'].number}"
                         for h in self._hits]
                run(self.session, f"select {' '.join(specs)}")
                run(self.session, f"color {' '.join(specs)} gold")
                run(self.session, f"style {' '.join(specs)} cartoon")
            self.results_button.setEnabled(True)
            self.session.logger.info(f"TM Helix: {len(self._hits)} candidate helix(es) found.")

    def _write_attributes(self, candidates):
        Residue.register_attr(self.session, "chopchop_tm_helix", "TM Helix", attr_type=bool)
        hit_residues = set()
        for h in self._hits:
            hit_residues.update(cellbio.residue_range(h["start_residue"], h["end_residue"], candidates))
        for res in candidates:
            res.chopchop_tm_helix = res in hit_residues

    def _show_results(self):
        dialog = QDialog(self.tool_window.ui_area)
        dialog.setWindowTitle("TM Helix - Results")
        dialog.setAttribute(Qt.WA_DeleteOnClose)
        dialog.resize(600, 500)
        layout = QVBoxLayout()

        hint = QLabel("Click a row to select that helix in the 3D view.")
        hint.setWordWrap(True)
        layout.addWidget(hint)

        table = QTableWidget(len(self._hits), 5)
        table.setHorizontalHeaderLabels(["Chain", "Residues", "Length", "Hydrophobicity (K-D)", "MLP score (Å²-weighted)"])
        for i, h in enumerate(self._hits):
            r1, r2 = h["start_residue"], h["end_residue"]
            mlp = f"{h['mlp_score']:.3f}" if h["mlp_score"] is not None else "n/a"
            values = [r1.chain_id, f"{r1.number}-{r2.number}", str(h["length"]), f"{h['hydrophobicity']:.2f}", mlp]
            for col, value in enumerate(values):
                item = QTableWidgetItem(value)
                item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                table.setItem(i, col, item)
        table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        table.verticalHeader().setVisible(False)
        table.setSortingEnabled(True)
        table.setSelectionBehavior(QAbstractItemView.SelectRows)
        table.cellClicked.connect(lambda row, _col: self._select_hit(row))
        layout.addWidget(table)

        close_button = QPushButton("Close")
        close_button.clicked.connect(dialog.close)
        layout.addWidget(close_button)
        dialog.setLayout(layout)
        dialog.finished.connect(lambda _result: run(self.session, "view"))
        dialog.show()

    def _select_hit(self, row):
        if row >= len(self._hits):
            return
        h = self._hits[row]
        r1, r2 = h["start_residue"], h["end_residue"]
        spec = f"/{r1.chain_id}:{r1.number}-{r2.number}"
        run(self.session, f"select {spec}")
        run(self.session, f"view {spec}")
