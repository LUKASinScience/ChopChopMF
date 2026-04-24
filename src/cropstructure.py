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
- Previewing the crop via "cartoon hide sel" 
"""

from chimerax.core.tools import ToolInstance
from chimerax.ui import MainToolWindow
from chimerax.core.commands import run
from Qt.QtWidgets import (QVBoxLayout, QLineEdit, QLabel, QPushButton, QComboBox,
                           QWidget, QTabWidget, QHBoxLayout)
from Qt.QtGui import QFont


class CropStructureTool(ToolInstance):
    SESSION_ENDURING = False
    SESSION_SAVE = True
    help = "help:user/tools/cropstructure.html"

    def __init__(self, session, tool_name):
        super().__init__(session, tool_name)
        self.display_name = "Crop Structure"
        self.tool_window = MainToolWindow(self)
        
        # Status for the preview logic
        self.preview_active = False
        
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

        # Button Area
        btn_layout = QHBoxLayout()
        
        self.preview_btn = QPushButton("Hide Deletion Preview")
        self.preview_btn.clicked.connect(self.toggle_preview)
        btn_layout.addWidget(self.preview_btn)

        crop_button = QPushButton("ChopChop Crop")
        crop_button.setStyleSheet("font-weight: bold; background-color: #d32f2f; color: white;")
        crop_button.clicked.connect(self.crop_structure)
        btn_layout.addWidget(crop_button)
        
        layout.addLayout(btn_layout)

        self._refresh_crop_models()
        widget.setLayout(layout)
        return widget

    def toggle_preview(self):
        """Executes preview via 'select' and 'hide' or performs an 'undo'."""
        if not self.preview_active:
            model_id = self.crop_model_selector.currentText().strip()
            chain_id = self.crop_chain_selector.currentText().strip()
            residue_range = self.residue_input.text().strip()
            
            if not (model_id and chain_id and residue_range):
                self.session.logger.warning("Select Model, Chain and Range first.")
                return

            try:
                # 1. Calculate residues that should NOT be kept
                residues_to_keep = self.parse_residue_range(residue_range)
                total_residues = self.get_all_residue_numbers(model_id, chain_id)
                to_remove = [n for n in total_residues if n not in residues_to_keep]

                if to_remove:
                    # 2. Command: Select and Hide
                    spec = f"#{model_id}/{chain_id}:{','.join(map(str, to_remove))}"
                    run(self.session, f"select {spec}")
                    run(self.session, "hide sel; cartoon hide sel")
                    
                    self.preview_active = True
                    self.preview_btn.setText("Reset Preview (Undo)")
                    self.session.logger.info(f"Preview active: Hidden residues {spec}")
            except Exception as e:
                self.session.logger.error(f"Preview Error: {e}")
        else:
            # 3. Undo function
            run(self.session, "undo")
            self.preview_active = False
            self.preview_btn.setText("Hide Deletion Preview")

    def crop_structure(self):
        """Deletes the currently selected (hidden) residues or recalculates them."""
        model_id = self.crop_model_selector.currentText().strip()
        chain_id = self.crop_chain_selector.currentText().strip()
        residue_range = self.residue_input.text().strip()

        if not (model_id and chain_id and residue_range):
            return

        try:
            # If preview is active, the residues are already selected
            if not self.preview_active:
                residues_to_keep = self.parse_residue_range(residue_range)
                total_residues = self.get_all_residue_numbers(model_id, chain_id)
                to_remove = [n for n in total_residues if n not in residues_to_keep]
                if to_remove:
                    spec = f"#{model_id}/{chain_id}:{','.join(map(str, to_remove))}"
                    run(self.session, f"select {spec}")
            
            # Execute the deletion process
            run(self.session, "delete sel")
            run(self.session, "select clear")
            self.preview_active = False
            self.preview_btn.setText("Hide Deletion Preview")
            self.session.logger.info(f"Cropped Model {model_id} Chain {chain_id}.")
        except Exception as e:
            self.session.logger.error(f"Error cropping structure: {e}")

    # --- Helper Functions ---

    def get_all_residue_numbers(self, model_id, chain_id):
        """Returns a list of all residue numbers for the chain."""
        model = next((m for m in self.session.models if m.id_string == model_id), None)
        if model:
            return [r.number for r in model.residues if r.chain_id == chain_id]
        return []

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

    def delete_chain(self):
        model_id = self.delete_model_selector.currentText().strip()
        chain_id = self.delete_chain_selector.currentText().strip()
        if model_id and chain_id:
            run(self.session, f"delete #{model_id}/{chain_id}")

    @staticmethod
    def parse_residue_range(residue_range):
        residues = []
        for part in residue_range.split(','):
            part = part.strip()
            if not part: continue
            if '-' in part:
                start, end = map(int, part.split('-'))
                residues.extend(range(start, end+1))
            else:
                residues.append(int(part))
        return residues
