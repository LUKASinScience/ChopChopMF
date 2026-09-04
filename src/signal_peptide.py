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
Signal Peptide Tool
Scans a chosen chain's sequence for subcellular targeting-signal motifs (NLS,
NES, ER retention, peroxisomal PTS1/PTS2 - see cellbio.py::SIGNAL_MOTIFS)
plus an N-terminal mitochondrial-presequence charge heuristic, keeping only
motifs that are surface-exposed (SASA above a threshold). Results are a
heuristic motif match, not a calibrated prediction - see
docs/acknowledgements.md for per-motif citations.
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

_ORGANELLES = sorted({organelle for _n, _p, organelle, _c, _a in cellbio.SIGNAL_MOTIFS} | {"Mitochondrion (import)"})
_ORGANELLE_COLORS = {name: CELLBIO_PALETTE[i % len(CELLBIO_PALETTE)] for i, name in enumerate(_ORGANELLES)}


class SignalPeptide(ToolInstance):

    SESSION_ENDURING = False
    SESSION_SAVE = True

    def __init__(self, session, tool_name):
        super().__init__(session, tool_name)
        self.display_name = "Signal Peptide"
        self.tool_window = MainToolWindow(self)
        self._hits = []
        self._mito_hit = None
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
        row.addWidget(QLabel("SASA (Å²) >"))
        self.sasa_spin = QDoubleSpinBox()
        self.sasa_spin.setRange(0.0, 1000.0)
        self.sasa_spin.setValue(5.0)
        row.addWidget(self.sasa_spin)
        layout.addLayout(row)

        self.run_button = QPushButton("ChopChop predict localization sequences")
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
                self.session.logger.status("Signal Peptide: measuring SASA...")
                run(self.session, f"measure sasa #{chain.structure.id_string}")
            except Exception as e:
                show_error(self.tool_window.ui_area, "ChopChopMF", f"SASA calculation failed: {e}")
                return
            sequence, residues = cellbio.sequence_and_residues(chain)
            self._hits = cellbio.scan_signal_motifs(sequence, residues, self.sasa_spin.value())
            is_mito, pos, neg = cellbio.scan_mitochondrial_presequence(sequence)
            self._mito_hit = {"positive": pos, "negative": neg} if is_mito else None
            self._write_attributes(residues)

            select_specs = []
            for h in self._hits:
                r1, r2 = h["start_residue"], h["end_residue"]
                spec = f"/{r1.chain_id}:{r1.number}-{r2.number}"
                select_specs.append(spec)
                run(self.session, f"color {spec} {_ORGANELLE_COLORS[h['organelle']]}")
            if select_specs:
                run(self.session, f"select {' '.join(select_specs)}")
                run(self.session, f"style {' '.join(select_specs)} sphere")
            self.results_button.setEnabled(True)
            self.session.logger.info(
                f"Signal Peptide: {len(self._hits)} motif hit(s)"
                f"{', N-terminal mitochondrial-import signature present' if self._mito_hit else ''}.")

    def _write_attributes(self, candidates):
        Residue.register_attr(self.session, "chopchop_signal_organelle", "Signal Peptide", attr_type=str)
        hit_map = {}
        for h in self._hits:
            for res in cellbio.residue_range(h["start_residue"], h["end_residue"], candidates):
                hit_map.setdefault(res, []).append(h["organelle"])
        for res in candidates:
            organelles = hit_map.get(res)
            if organelles:
                res.chopchop_signal_organelle = ", ".join(sorted(set(organelles)))
            else:
                try:
                    delattr(res, "chopchop_signal_organelle")
                except AttributeError:
                    pass

    def _show_results(self):
        dialog = QDialog(self.tool_window.ui_area)
        dialog.setWindowTitle("Signal Peptide - Results")
        dialog.setAttribute(Qt.WA_DeleteOnClose)
        dialog.resize(650, 500)
        layout = QVBoxLayout()

        if self._mito_hit:
            mito_label = QLabel(
                f"N-terminal mitochondrial-import signature: {self._mito_hit['positive']} positive / "
                f"{self._mito_hit['negative']} negative charges in first {cellbio.MTS_WINDOW} residues "
                f"({cellbio.MTS_CITATION}).")
            mito_label.setWordWrap(True)
            layout.addWidget(mito_label)

        organelles_present = {h["organelle"] for h in self._hits}
        legend = QLabel(" &nbsp; ".join(
            f'<span style="background-color:{_ORGANELLE_COLORS[org]}; color:white; padding:1px 5px; '
            f'border-radius:3px;">{org}</span>'
            for org in _ORGANELLES if org in organelles_present))
        legend.setWordWrap(True)
        layout.addWidget(legend)
        hint = QLabel("Click a row to select that region in the 3D view.")
        hint.setWordWrap(True)
        layout.addWidget(hint)

        table = QTableWidget(len(self._hits), 6)
        table.setHorizontalHeaderLabels(["Signal", "Organelle", "Chain", "Residues", "Sequence", "SASA (Å²)"])
        for i, h in enumerate(self._hits):
            r1, r2 = h["start_residue"], h["end_residue"]
            name = h["name"] + (" (approx.)" if h["approximate"] else "")
            values = [name, h["organelle"], r1.chain_id, f"{r1.number}-{r2.number}", h["seq"], f"{h['sasa']:.1f}"]
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
