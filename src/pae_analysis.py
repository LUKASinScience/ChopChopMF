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
PAE Analysis Tool
A ChimeraX tool for performing AlphaFold error analysis using PAE.
The tool has two tabs:
1. Standard PAE Contact Analysis
2. Residue Selection from Pseudobond Model
Standard output messages are printed to the terminal.
"""

import webbrowser
from chimerax.core.tools import ToolInstance
from chimerax.core.commands import run
from chimerax.ui import MainToolWindow
from Qt.QtWidgets import (QVBoxLayout, QLabel, QPushButton, QWidget, QComboBox, QTabWidget)
from Qt.QtGui import QFont

class PAEAnalysis(ToolInstance):

    SESSION_ENDURING = False
    SESSION_SAVE = True
    help = "https://www.cgl.ucsf.edu/chimerax/docs/user/tools/pae_analysis.html"

    def __init__(self, session, tool_name):
        super().__init__(session, tool_name)
        self.json_loaded = False
        self.display_name = "PAE Analysis Tool"
        self.tool_window = MainToolWindow(self)
        self._build_ui()
        self.tool_window.manage('side')

    def _build_ui(self):
        layout = QVBoxLayout()

        title_label = QLabel("PAE Analysis Tool")
        bold_font = QFont()
        bold_font.setBold(True)
        title_label.setFont(bold_font)
        layout.addWidget(title_label)

        tabs = QTabWidget()
        tabs.addTab(self._create_contacts_tab(), "1. PAE Contacts")
        tabs.addTab(self._create_residue_tab(), "2. PAE Contact Residues")
        layout.addWidget(tabs)

        container = QWidget()
        container.setLayout(layout)
        self.tool_window.ui_area.setLayout(layout)

        self._refresh_chain_list()

    def _create_contacts_tab(self):
        widget = QWidget()
        layout = QVBoxLayout()

        layout.addWidget(QLabel("PAE (Predicted Aligned Error):"))

        instructions = QLabel(
            "PAE Analysis requires exactly one open model and its corresponding .json file.\n"
            "Click the 'Load .json file' button to show the AlphaFold Error Plot tool."
        )
        instructions.setWordWrap(True)
        layout.addWidget(instructions)

        refresh_button = QPushButton("↻ Refresh model list")
        refresh_button.clicked.connect(self._refresh_chain_list)
        layout.addWidget(refresh_button)

        layout.addWidget(QLabel("Select first chain:"))
        self.chain1_selector = QComboBox()
        layout.addWidget(self.chain1_selector)

        layout.addWidget(QLabel("Select second chain:"))
        self.chain2_selector = QComboBox()
        layout.addWidget(self.chain2_selector)

        layout.addWidget(QLabel("Select distance (Å):"))
        self.distance_selector = QComboBox()
        for d in range(3, 11):
            self.distance_selector.addItem(str(d))
        layout.addWidget(self.distance_selector)

        load_json_button = QPushButton("Load .json file")
        load_json_button.clicked.connect(self._load_json_file)
        layout.addWidget(load_json_button)

        run_button = QPushButton("ChopChop PAE")
        run_button.clicked.connect(self._run_pae_analysis)
        layout.addWidget(run_button)

        widget.setLayout(layout)
        return widget

    def _create_residue_tab(self):
        widget = QWidget()
        layout = QVBoxLayout()

        instructions = QLabel(
            "This tab selects residues that are connected by pseudobonds in the PAE Contacts model.\n"
            "It allows you to investigate the interactions between the two chains more thoroughly\n"
            "based on the settings from the ChopChop PAE analysis."
        )
        instructions.setWordWrap(True)
        layout.addWidget(instructions)

        run_button = QPushButton("ChopChop PAE interaction Residues")
        run_button.clicked.connect(self._run_pae_selection)
        layout.addWidget(run_button)

        explanation = QLabel(
            "A more precise analysis can be done by using the Webbrowser 'PAE Viewer'\n"
            "developed by Christoph Elfmann and Jörg Stülke."
        )
        layout.addWidget(explanation)

        pae_viewer_button = QPushButton("Open PAE Viewer")
        pae_viewer_button.clicked.connect(lambda: webbrowser.open("https://subtiwiki.uni-goettingen.de/v4/paeViewerDemo"))
        layout.addWidget(pae_viewer_button)

        widget.setLayout(layout)
        return widget

    def _refresh_chain_list(self):
        self.chain1_selector.clear()
        self.chain2_selector.clear()
        models = [m for m in self.session.models.list() if hasattr(m, "residues")]
        if len(models) != 1:
            self.session.logger.error("PAE Analysis requires exactly one open model with residue information.")
            return
        model = models[0]
        chains = sorted(set(res.chain_id for res in model.residues))
        for chain in chains:
            self.chain1_selector.addItem(chain)
            self.chain2_selector.addItem(chain)
        if "A" in chains:
            self.chain1_selector.setCurrentText("A")
        if "B" in chains:
            self.chain2_selector.setCurrentText("B")
        self.distance_selector.setCurrentText("5")
        self.session.logger.info("Chain lists refreshed.")

    def _load_json_file(self):
        run(self.session, 'ui tool show "AlphaFold Error Plot"')
        self.json_loaded = True

    def _run_pae_analysis(self):
        models = [m for m in self.session.models.list() if hasattr(m, "residues")]
        if len(models) != 1:
            self.session.logger.error("PAE Analysis requires exactly one open model with residue information.")
            return
        if not self.json_loaded:
            self.session.logger.error("Please load the corresponding .json file before running the analysis.")
            return
        chain1 = self.chain1_selector.currentText()
        chain2 = self.chain2_selector.currentText()
        if not chain1 or not chain2 or chain1 == chain2:
            self.session.logger.error("Both chains must be selected and different.")
            return
        distance = self.distance_selector.currentText()
        command = f"alphafold contacts /{chain1} to /{chain2} distance {distance}"
        self.session.logger.info(f"Running command: {command}")
        run(self.session, command)

    def _run_pae_selection(self):
        pb_model = next((m for m in self.session.models.list() if m.name == "PAE Contacts"), None)
        if pb_model is None:
            self.session.logger.error("No pseudobond model named 'PAE Contacts' found.")
            return
        residues = set()
        for pb in pb_model.pseudobonds:
            for atom in pb.atoms:
                residues.add(atom.residue)
        for res in residues:
            self.session.logger.info(f"Residue: chain {res.chain_id}, number {res.number}, name {res.name}")
        run(self.session, "select clear")
        for res in residues:
            for atom in res.atoms:
                atom.selected = True
        run(self.session, "show sel")
        run(self.session, "style sel stick")
        run(self.session, "color sel bychain")
        run(self.session, "color sel byhetero")
        pb_model.delete()
        self.session.logger.info("Pseudobond model deleted. Residues remain selected and styled by heteroatom.")

from chimerax.core.toolshed import BundleAPI

class _PAEAnalysisBundleAPI(BundleAPI):
    @staticmethod
    def start_tool(session, tool_name):
        return PAEAnalysis(session, tool_name)

bundle_api = _PAEAnalysisBundleAPI()

