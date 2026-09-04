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
Phospho Sites Tool
Scans a chosen chain's sequence for kinase phosphorylation-site consensus
motifs (see cellbio.py::PHOSPHO_MOTIFS), keeping only S/T/Y residues that are
both disordered (pLDDT below a threshold) and solvent-exposed (SASA above a
threshold) - the standard "flexible and accessible" filter for a real
phospho-acceptor site. Results are a heuristic motif match, not a calibrated
prediction - see docs/acknowledgements.md for per-motif citations.
"""

from chimerax.core.tools import ToolInstance
from chimerax.core.commands import run
from chimerax.atomic import Residue
from chimerax.ui import MainToolWindow
from Qt.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QComboBox, QDoubleSpinBox,
    QDialog, QTableWidget, QTableWidgetItem, QHeaderView, QAbstractItemView,
)
from Qt.QtCore import Qt

from .utils import make_guide_button, busy_cursor, show_error, CELLBIO_PALETTE
from . import cellbio

_KINASE_COLORS = {
    name: CELLBIO_PALETTE[i % len(CELLBIO_PALETTE)]
    for i, (name, *_rest) in enumerate(cellbio.PHOSPHO_MOTIFS)
}


class PhosphoSites(ToolInstance):

    SESSION_ENDURING = False
    SESSION_SAVE = True

    def __init__(self, session, tool_name):
        super().__init__(session, tool_name)
        self.display_name = "Phospho Sites"
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
        row.addWidget(QLabel("pLDDT <"))
        self.plddt_spin = QDoubleSpinBox()
        self.plddt_spin.setRange(0.0, 100.0)
        self.plddt_spin.setValue(70.0)
        row.addWidget(self.plddt_spin)
        row.addWidget(QLabel("SASA (Å²) >"))
        self.sasa_spin = QDoubleSpinBox()
        self.sasa_spin.setRange(0.0, 1000.0)
        self.sasa_spin.setValue(5.0)
        row.addWidget(self.sasa_spin)
        layout.addLayout(row)

        self.run_button = QPushButton("ChopChop predict Phospho Sites")
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
            try:
                self.session.logger.status("Phospho Sites: measuring SASA...")
                run(self.session, f"measure sasa #{chain.structure.id_string}")
            except Exception as e:
                show_error(self.tool_window.ui_area, "ChopChopMF", f"SASA calculation failed: {e}")
                return
            sequence, residues = cellbio.sequence_and_residues(chain)
            self._hits = cellbio.scan_phospho_sites(
                sequence, residues, self.plddt_spin.value(), self.sasa_spin.value())
            self._write_attributes(chain.structure, residues)
            self.results_button.setEnabled(True)

            if self._hits:
                sel = " ".join(f"/{h['residue'].chain_id}:{h['residue'].number}" for h in self._hits)
                run(self.session, f"select {sel}")
                run(self.session, f"style {sel} sphere")
                for h in self._hits:
                    kinase = h["kinases"][0][0]
                    r = h["residue"]
                    run(self.session, f"color /{r.chain_id}:{r.number} {_KINASE_COLORS[kinase]}")
            self.session.logger.info(f"Phospho Sites: {len(self._hits)} residue(s) flagged.")

    def _write_attributes(self, structure, candidates):
        Residue.register_attr(self.session, "chopchop_phospho_kinase", "Phospho Sites", attr_type=str)
        hit_map = {h["residue"]: h for h in self._hits}
        for res in candidates:
            hit = hit_map.get(res)
            if hit is not None:
                res.chopchop_phospho_kinase = ", ".join(k for k, _c, _a in hit["kinases"])
            else:
                try:
                    delattr(res, "chopchop_phospho_kinase")
                except AttributeError:
                    pass

    def _show_results(self):
        dialog = QDialog(self.tool_window.ui_area)
        dialog.setWindowTitle("Phospho Sites - Results")
        dialog.setAttribute(Qt.WA_DeleteOnClose)
        dialog.resize(650, 550)
        layout = QVBoxLayout()

        legend = QLabel(" &nbsp; ".join(
            f'<span style="background-color:{color}; color:white; padding:1px 5px; '
            f'border-radius:3px;">{kinase}</span>'
            for kinase, color in _KINASE_COLORS.items()))
        legend.setWordWrap(True)
        layout.addWidget(legend)
        hint = QLabel("Click a row to select that residue in the 3D view.")
        hint.setWordWrap(True)
        layout.addWidget(hint)

        table = QTableWidget(len(self._hits), 5)
        table.setHorizontalHeaderLabels(["Chain", "Residue", "Kinase(s)", "pLDDT", "SASA (Å²)"])
        for i, h in enumerate(self._hits):
            r = h["residue"]
            kinases = ", ".join(f"{k}{' (approx.)' if a else ''}" for k, _c, a in h["kinases"])
            values = [r.chain_id, f"{r.name}{r.number}", kinases, f"{h['plddt']:.1f}", f"{h['sasa']:.1f}"]
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
        r = self._hits[row]["residue"]
        run(self.session, f"select /{r.chain_id}:{r.number}")
        run(self.session, f"view /{r.chain_id}:{r.number}")
