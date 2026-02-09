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
CropStructureTool.py
A ChimeraX tool for:
- Cropping a structure by removing residues not in a specified range
- Deleting an entire chain from the model
It uses a tabbed interface for easy selection and consistent style with other tools.
"""

from chimerax.core.tools import ToolInstance
from chimerax.ui import MainToolWindow
from chimerax.core.commands import run
from Qt.QtWidgets import (QVBoxLayout, QLineEdit, QLabel, QPushButton, QComboBox,
                           QWidget, QTabWidget)
from Qt.QtGui import QFont


class CropStructureTool(ToolInstance):
    SESSION_ENDURING = False
    SESSION_SAVE = True
    help = "help:user/tools/cropstructure.html"

    def __init__(self, session, tool_name):
        super().__init__(session, tool_name)
        self.display_name = "Crop Structure"
        self.tool_window = MainToolWindow(self)
        self._build_ui()

    def _build_ui(self):
        layout = QVBoxLayout()

        title_label = QLabel("Crop Structure Tool")
        title_font = QFont()
        title_font.setPointSize(14)
        title_font.setBold(True)
        title_label.setFont(title_font)
        layout.addWidget(title_label)

        tabs = QTabWidget()
        tabs.addTab(self._create_crop_tab(), "Crop Residues")
        tabs.addTab(self._create_delete_tab(), "Delete Chain")
        layout.addWidget(tabs)

        container = QWidget()
        container.setLayout(layout)
        self.tool_window.ui_area.setLayout(layout)
        self.tool_window.manage('side')

    def _create_crop_tab(self):
        widget = QWidget()
        layout = QVBoxLayout()

        layout.addWidget(QLabel("Select model:"))
        self.crop_model_selector = QComboBox()
        layout.addWidget(self.crop_model_selector)

        layout.addWidget(QLabel("Select chain:"))
        self.crop_chain_selector = QComboBox()
        layout.addWidget(self.crop_chain_selector)

        refresh_button = QPushButton("↻ Refresh model list and chains")
        refresh_button.clicked.connect(self._refresh_crop_models)
        layout.addWidget(refresh_button)

        layout.addWidget(QLabel("Residue range to keep (e.g., 1-50,60):"))
        self.residue_input = QLineEdit()
        layout.addWidget(self.residue_input)

        crop_button = QPushButton("ChopChop Crop")
        crop_button.clicked.connect(self.crop_structure)
        layout.addWidget(crop_button)

        self._refresh_crop_models()
        widget.setLayout(layout)
        return widget

    def _create_delete_tab(self):
        widget = QWidget()
        layout = QVBoxLayout()

        layout.addWidget(QLabel("Select model:"))
        self.delete_model_selector = QComboBox()
        layout.addWidget(self.delete_model_selector)

        layout.addWidget(QLabel("Select chain:"))
        self.delete_chain_selector = QComboBox()
        layout.addWidget(self.delete_chain_selector)

        refresh_button = QPushButton("Refresh models and chains")
        refresh_button.clicked.connect(self._refresh_delete_models)
        layout.addWidget(refresh_button)

        delete_button = QPushButton("Delete Chain")
        delete_button.clicked.connect(self.delete_chain)
        layout.addWidget(delete_button)

        self._refresh_delete_models()
        widget.setLayout(layout)
        return widget

    def _refresh_crop_models(self):
        self._populate_models(self.crop_model_selector)
        self._refresh_chain_list(self.crop_model_selector, self.crop_chain_selector)

        self.crop_model_selector.currentIndexChanged.connect(
            lambda: self._refresh_chain_list(self.crop_model_selector, self.crop_chain_selector))

    def _refresh_delete_models(self):
        self._populate_models(self.delete_model_selector)
        self._refresh_chain_list(self.delete_model_selector, self.delete_chain_selector)

        self.delete_model_selector.currentIndexChanged.connect(
            lambda: self._refresh_chain_list(self.delete_model_selector, self.delete_chain_selector))

    def _populate_models(self, selector):
        selector.clear()
        for mdl in self.session.models.list():
            if hasattr(mdl, "residues"):
                selector.addItem(mdl.id_string)

    def _refresh_chain_list(self, model_selector, chain_selector):
        chain_selector.clear()
        model_id = model_selector.currentText().strip()
        model = next((m for m in self.session.models if m.id_string == model_id), None)
        if model:
            chains = sorted(set(r.chain_id for r in model.residues))
            for cid in chains:
                chain_selector.addItem(cid)

    def crop_structure(self):
        model_id = self.crop_model_selector.currentText().strip()
        chain_id = self.crop_chain_selector.currentText().strip()
        residue_range = self.residue_input.text().strip()
        if not (model_id and chain_id and residue_range):
            self.session.logger.warning("Fill in Model, Chain, and Residue Range fields.")
            return
        try:
            self.session.logger.info(f"Parsing residue range: '{residue_range}'")
            residues_to_keep = self.parse_residue_range(residue_range)
            total = self.get_total_residues(model_id, chain_id)
            to_remove = [r for r in range(1, total+1) if r not in residues_to_keep]
            if to_remove:
                sel_cmd = f"select #{model_id}/{chain_id}:{','.join(map(str, to_remove))}"
                self.session.logger.info(f"Selecting residues to delete: {sel_cmd}")
                run(self.session, sel_cmd)
                run(self.session, "delete sel")
                run(self.session, "select clear")
            self.session.logger.info(f"Cropped Model {model_id} Chain {chain_id}.")
        except Exception as e:
            self.session.logger.error(f"Error cropping structure: {e}")

    def delete_chain(self):
        model_id = self.delete_model_selector.currentText().strip()
        chain_id = self.delete_chain_selector.currentText().strip()
        if not (model_id and chain_id):
            self.session.logger.warning("Please select a model and chain.")
            return
        try:
            self.session.logger.info(f"Deleting Chain {chain_id} from Model {model_id}...")
            run(self.session, f"select #{model_id}/{chain_id}")
            run(self.session, "delete sel")
            run(self.session, "select clear")
            self.session.logger.info(f"Deleted Chain {chain_id} from Model {model_id}.")
        except Exception as e:
            self.session.logger.error(f"Error deleting chain: {e}")

    def get_total_residues(self, model_id, chain_id):
        model = next((m for m in self.session.models if m.id_string == model_id), None)
        if not model:
            raise ValueError(f"Model {model_id} not found.")
        residues = [r for r in model.residues if r.chain_id == chain_id]
        return len(residues)

    @staticmethod
    def parse_residue_range(residue_range):
        residues = []
        for part in residue_range.split(','):
            part = part.strip()
            if not part:
                continue
            if '-' in part:
                start, end = map(int, part.split('-'))
                residues.extend(range(start, end+1))
            else:
                residues.append(int(part))
        return residues
