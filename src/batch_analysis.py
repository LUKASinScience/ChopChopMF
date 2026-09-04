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
Batch Analysis Tool
Runs PAE Analysis's interface scores (pDockQ, buried area, H-bonds, LIS/cLIS/
iLIS, ipSAE d0chn) across every currently open model at once, instead of one
chain pair at a time in PAE Analysis itself. Pick a chain-role rule once -
either the same chain letters on every model, or "first chain(s) vs. last
chain" (borrowed from PPIScreenML's own batch convention, for batches whose
models don't all share one chain-naming scheme) - and get one aggregated,
exportable table.

A model missing PAE data is never silently skipped without a look: for the
AlphaFold3-server file naming convention (`..._model_<N>.cif` next to
`..._full_data_<N>.json`, the exact pattern AlphaFold3-server output uses),
the matching PAE file is found and loaded automatically - verified by the
file actually existing at that exact name, not guessed at. This deliberately
does NOT try to guess other tools' naming conventions (ColabFold, local
AlphaFold-Multimer, etc.) - a model in one of those formats still needs its
PAE `.json` loaded once via PAE Analysis's own "Load .json file" first, and
is reported as skipped (not silently left out) if it isn't.

Two input modes: score whatever's already open in ChimeraX, or point at a
folder of many `.cif`/`.pdb` files - each one is opened, scored, and closed
again before the next, so a batch of many predictions never needs them all
open in ChimeraX at once (PAE Analysis itself only ever works with exactly
one open model - this second mode exists specifically so a folder of many
predictions doesn't require fighting that restriction by hand).
"""

import csv
import re
from pathlib import Path

from chimerax.core.tools import ToolInstance
from chimerax.core.commands import run
from chimerax.atomic import Residue
from chimerax.ui import MainToolWindow
from Qt.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QWidget, QComboBox,
    QLineEdit, QTableWidget, QTableWidgetItem, QFileDialog, QTabWidget,
)
from Qt.QtCore import Qt

from .utils import make_scrollable, make_guide_button, show_error, busy_cursor, get_settings
from . import scoring

_HEADERS = [
    "Model", "Name", "Chain(s) 1", "Chain(s) 2", "Status", "pDockQ", "Contacts",
    "Buried area (Å²)", "H-bonds", "LIS", "cLIS", "iLIS", "ipSAE (d0chn)",
]


class BatchAnalysis(ToolInstance):

    SESSION_ENDURING = False
    SESSION_SAVE = False

    def __init__(self, session, tool_name):
        super().__init__(session, tool_name)
        self.display_name = "Batch Analysis"
        self._rows = []
        self.tool_window = MainToolWindow(self)
        self._build_ui()
        self.tool_window.manage("side")

    def _build_ui(self):
        outer_layout = QVBoxLayout()
        tabs = QTabWidget()
        tabs.addTab(make_scrollable(self._create_batch_tab()), "Batch")
        outer_layout.addWidget(tabs)
        self.tool_window.ui_area.setLayout(outer_layout)
        self._update_rule_visibility()
        self._update_mode_visibility()

    def _create_batch_tab(self):
        widget = QWidget()
        layout = QVBoxLayout()
        layout.setAlignment(Qt.AlignTop)
        layout.setSpacing(4)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.addWidget(make_guide_button("4-analyze-structure"))

        info = QLabel(
            "Runs PAE Analysis's interface scores across many models at once. "
            "For AlphaFold3-server files (“…_model_N.cif” next to “…_full_data_N.json”), "
            "the matching PAE file is found and loaded automatically. Other formats "
            "(ColabFold, local AlphaFold-Multimer) still need PAE Analysis's own "
            "“Load .json file” first, once per model - a model with no PAE data "
            "either way is reported as skipped, never guessed at."
        )
        info.setWordWrap(True)
        layout.addWidget(info)

        mode_row = QHBoxLayout()
        mode_row.addWidget(QLabel("Input:"))
        self.mode_combo = QComboBox()
        self.mode_combo.addItems([
            "Already-open models",
            "Folder of files (opens, scores, and closes each one automatically)",
        ])
        self.mode_combo.currentIndexChanged.connect(self._update_mode_visibility)
        mode_row.addWidget(self.mode_combo, stretch=1)
        layout.addLayout(mode_row)

        self.folder_row = QWidget()
        folder_layout = QHBoxLayout()
        folder_layout.setContentsMargins(0, 0, 0, 0)
        folder_layout.addWidget(QLabel("Folder:"))
        self.folder_edit = QLineEdit()
        self.folder_edit.setPlaceholderText("Folder containing the .cif/.pdb files")
        folder_layout.addWidget(self.folder_edit, stretch=1)
        browse_button = QPushButton("Browse…")
        browse_button.clicked.connect(self._choose_folder)
        folder_layout.addWidget(browse_button)
        self.folder_row.setLayout(folder_layout)
        layout.addWidget(self.folder_row)

        rule_row = QHBoxLayout()
        rule_row.addWidget(QLabel("Chain rule:"))
        self.rule_combo = QComboBox()
        self.rule_combo.addItems([
            "Same chain letters on every model",
            "First chain(s) = side 1, last chain = side 2 (PPIScreenML-style)",
        ])
        self.rule_combo.currentIndexChanged.connect(self._update_rule_visibility)
        rule_row.addWidget(self.rule_combo, stretch=1)
        layout.addLayout(rule_row)

        self.letters_row = QWidget()
        letters_layout = QHBoxLayout()
        letters_layout.setContentsMargins(0, 0, 0, 0)
        letters_layout.addWidget(QLabel("Chain 1:"))
        self.chain1_edit = QLineEdit("A")
        self.chain1_edit.setMaximumWidth(50)
        letters_layout.addWidget(self.chain1_edit)
        letters_layout.addWidget(QLabel("Chain 2:"))
        self.chain2_edit = QLineEdit("B")
        self.chain2_edit.setMaximumWidth(50)
        letters_layout.addWidget(self.chain2_edit)
        letters_layout.addStretch()
        self.letters_row.setLayout(letters_layout)
        layout.addWidget(self.letters_row)

        self.run_button = QPushButton("Run Batch")
        self.run_button.clicked.connect(self._run_batch)
        layout.addWidget(self.run_button)

        self.table = QTableWidget(0, len(_HEADERS))
        self.table.setHorizontalHeaderLabels(_HEADERS)
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.setMinimumHeight(220)
        layout.addWidget(self.table)

        export_button = QPushButton("Export as CSV…")
        export_button.clicked.connect(self._export_csv)
        layout.addWidget(export_button)

        widget.setLayout(layout)
        return widget

    def _update_rule_visibility(self):
        self.letters_row.setVisible(self.rule_combo.currentIndex() == 0)

    def _update_mode_visibility(self):
        self.folder_row.setVisible(self.mode_combo.currentIndex() == 1)

    def _choose_folder(self):
        folder = QFileDialog.getExistingDirectory(self.tool_window.ui_area, "Choose Folder")
        if folder:
            self.folder_edit.setText(folder)

    def _resolve_chain_groups(self, model):
        """(group1, group2, residues1, residues2) for this model under the
        current chain-role rule, or None if the rule doesn't apply to it
        (e.g. the requested letters aren't present). Only protein (amino
        acid) chains are considered - an AlphaFold3 prediction with a bound
        ligand/ion commonly gets those as their own separate "chain" IDs
        (e.g. an ATP or Mg2+ as chain C/D/E), which would otherwise silently
        get swept into the "first chain(s) vs. last chain" rule's groups and
        produce meaningless (all-zero) scores."""
        chain_ids = sorted({
            r.chain_id for r in model.residues if r.polymer_type == Residue.PT_AMINO
        })
        if len(chain_ids) < 2:
            return None
        if self.rule_combo.currentIndex() == 0:
            c1 = self.chain1_edit.text().strip()
            c2 = self.chain2_edit.text().strip()
            if not c1 or not c2 or c1 == c2 or c1 not in chain_ids or c2 not in chain_ids:
                return None
            group1, group2 = {c1}, {c2}
        else:
            group2 = {chain_ids[-1]}
            group1 = set(chain_ids[:-1])
        residues1 = [r for r in model.residues if r.chain_id in group1 and r.polymer_type == Residue.PT_AMINO]
        residues2 = [r for r in model.residues if r.chain_id in group2 and r.polymer_type == Residue.PT_AMINO]
        if not residues1 or not residues2:
            return None
        return group1, group2, residues1, residues2

    def _find_af3_pae_json(self, model):
        """For an AlphaFold3-server-style model file '..._model_<N>.cif' (or
        .pdb), look for the matching '..._full_data_<N>.json' the server
        writes next to it - the exact naming convention, not a fuzzy guess.
        Returns the Path if that exact file exists, else None (other tools'
        PAE file naming, e.g. ColabFold or local AlphaFold-Multimer, isn't
        covered here - those still need PAE Analysis's own file picker)."""
        path = getattr(model, "filename", None)
        if not path:
            return None
        p = Path(path)
        m = re.match(r"^(.*)_model_(\d+)$", p.stem)
        if not m:
            return None
        candidate = p.parent / f"{m.group(1)}_full_data_{m.group(2)}.json"
        return candidate if candidate.exists() else None

    def _ensure_pae_loaded(self, model):
        """The model's alphafold_pae if already loaded, or if not, whatever
        loading its auto-detected AlphaFold3-server PAE file (see
        _find_af3_pae_json) achieves - None if neither works."""
        pae_obj = getattr(model, "alphafold_pae", None)
        if pae_obj is not None:
            return pae_obj
        json_path = self._find_af3_pae_json(model)
        if json_path is None:
            return None
        try:
            run(self.session, f'alphafold pae #{model.id_string} file "{json_path}" plot false')
        except Exception as e:
            self.session.logger.warning(f"Batch Analysis: found {json_path.name} for model #{model.id_string} but loading it failed: {e}")
            return None
        pae_obj = getattr(model, "alphafold_pae", None)
        if pae_obj is not None:
            self.session.logger.info(f"Batch Analysis: auto-loaded {json_path.name} for model #{model.id_string}.")
        return pae_obj

    def _score_model(self, model):
        """One results row for this already-open model: chain groups, then
        PAE-load, then the actual scores - or a skip/error row with the same
        column count if any step along the way doesn't work out."""
        groups = self._resolve_chain_groups(model)
        if groups is None:
            return [model.id_string, model.name, "", "", "skipped: chain rule didn't match"] + [""] * 8
        group1, group2, residues1, residues2 = groups
        c1_label, c2_label = ",".join(sorted(group1)), ",".join(sorted(group2))
        pae_obj = self._ensure_pae_loaded(model)
        if pae_obj is None:
            return [model.id_string, model.name, c1_label, c2_label, "skipped: no PAE data loaded"] + [""] * 8
        try:
            pdockq, n_contacts = scoring.compute_pdockq(residues1, residues2)
            buried_area = scoring.compute_buried_area(residues1, residues2)
            hbonds = scoring.compute_hbond_count(self.session, model, group1, group2)
            lis, clis, ilis = scoring.compute_lis_clis_ilis(residues1, residues2, pae_obj)
            ipsae = scoring.compute_ipsae_d0chn(residues1, residues2, pae_obj)
        except Exception as e:
            return [model.id_string, model.name, c1_label, c2_label, f"error: {e}"] + [""] * 8
        return [
            model.id_string, model.name, c1_label, c2_label, "ok",
            f"{pdockq:.3f}", str(n_contacts), f"{buried_area:.0f}", str(hbonds),
            f"{lis:.3f}", f"{clis:.3f}", f"{ilis:.3f}", f"{ipsae:.3f}",
        ]

    def _run_batch(self):
        if self.mode_combo.currentIndex() == 0:
            self._run_batch_open_models()
        else:
            self._run_batch_folder()

    def _run_batch_open_models(self):
        rows = [
            self._score_model(model) for model in self.session.models.list()
            if hasattr(model, "residues")
        ]
        self._rows = rows
        self._populate_table()
        if not rows:
            show_error(self.tool_window.ui_area, "ChopChopMF", "No open models with residues found.")

    def _run_batch_folder(self):
        folder = Path(self.folder_edit.text().strip())
        if not folder.is_dir():
            show_error(self.tool_window.ui_area, "ChopChopMF", "Pick a valid folder first.")
            return
        structure_files = sorted(p for p in folder.iterdir() if p.suffix.lower() in (".cif", ".pdb"))
        if not structure_files:
            show_error(self.tool_window.ui_area, "ChopChopMF", "No .cif/.pdb files found in that folder.")
            return

        rows = []
        with busy_cursor(self.run_button):
            for i, path in enumerate(structure_files, 1):
                self.session.logger.info(f"Batch Analysis: {i}/{len(structure_files)} - {path.name}")
                try:
                    opened = run(self.session, f'open "{path}"')
                except Exception as e:
                    rows.append(["", path.name, "", "", f"error: failed to open ({e})"] + [""] * 8)
                    continue
                model = opened[0] if opened else None
                if model is None or not hasattr(model, "residues"):
                    rows.append(["", path.name, "", "", "error: did not open as an atomic structure"] + [""] * 8)
                    continue
                try:
                    rows.append(self._score_model(model))
                finally:
                    # Always close before the next file, success or not - the
                    # whole point of folder mode is never needing every
                    # prediction open in ChimeraX at once.
                    try:
                        run(self.session, f"close #{model.id_string}")
                    except Exception as e:
                        self.session.logger.warning(f"Batch Analysis: failed to close #{model.id_string} ({path.name}): {e}")
        self._rows = rows
        self._populate_table()

    def _populate_table(self):
        self.table.setRowCount(len(self._rows))
        for i, row in enumerate(self._rows):
            for col, value in enumerate(row):
                self.table.setItem(i, col, QTableWidgetItem(value))
        self.table.resizeColumnsToContents()

    def _export_csv(self):
        if not self._rows:
            show_error(self.tool_window.ui_area, "ChopChopMF", "Run the batch first.")
            return
        export_dir = get_settings(self.session).export_dir
        file_path, _ = QFileDialog.getSaveFileName(
            self.tool_window.ui_area, "Save Batch Results As",
            str(Path(export_dir) / "batch_analysis.csv"), "CSV Files (*.csv)"
        )
        if not file_path:
            return
        try:
            with open(file_path, "w", newline="", encoding="utf-8-sig") as f:
                writer = csv.writer(f)
                writer.writerow(_HEADERS)
                writer.writerows(self._rows)
        except OSError as e:
            show_error(self.tool_window.ui_area, "ChopChopMF", f"Failed to save results:\n{e}")
