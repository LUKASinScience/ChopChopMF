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
Foldseek Analysis Tool
A ChimeraX tool to perform structural alignment and analysis using Foldseek.
This updated version includes a selector for the target database (PDB or AlphaFold DB).
Standard output messages are printed to the terminal.
"""

from chimerax.core.tools import ToolInstance
from chimerax.core.commands import run
from chimerax.ui import MainToolWindow
from Qt.QtWidgets import QVBoxLayout, QLabel, QPushButton, QWidget, QComboBox, QTabWidget
from Qt.QtGui import QFont, QDesktopServices
from Qt.QtCore import QUrl, Qt

class FoldseekAnalysis(ToolInstance):

    SESSION_ENDURING = False
    SESSION_SAVE = True
    help = "https://www.cgl.ucsf.edu/chimerax/docs/user/tools/foldseek.html"

    def __init__(self, session, tool_name):
        super().__init__(session, tool_name)

        # Set name displayed on title bar
        self.display_name = "Foldseek Analysis Tool"
        self.tool_window = MainToolWindow(self)
        self._build_ui()
        self.tool_window.manage('side')

    def _build_ui(self):
        # Create tab widget
        tabs = QTabWidget()

        # First tab: Main functionality
        main_tab = QWidget()
        main_layout = QVBoxLayout()

        # Title
        title = QLabel("Foldseek: Structural Alignment and Analysis")
        title_font = QFont()
        title_font.setBold(True)
        title_font.setPointSize(14)
        title.setFont(title_font)
        main_layout.addWidget(title)

        # Description
        description = QLabel("Foldseek provides tools for structural alignment and analysis. Select the query structure and the target database.")
        description.setWordWrap(True)
        main_layout.addWidget(description)
        
        # --- Database Selection (NEW) ---
        database_label = QLabel("Select Target Database:")
        main_layout.addWidget(database_label)

        self.database_selector = QComboBox()
        self.database_selector.addItem("PDB (default)")
        self.database_selector.addItem("AlphaFold DB (afdb50)")
        main_layout.addWidget(self.database_selector)
        # ---------------------------------

        # Model Selection
        model_label = QLabel("Select Model:")
        main_layout.addWidget(model_label)

        self.model_selector = QComboBox()
        self.chain_selector = QComboBox()
        self._populate_models_and_chains()
        main_layout.addWidget(self.model_selector)

        # Chain Selection
        chain_label = QLabel("Select Chain:")
        main_layout.addWidget(chain_label)
        main_layout.addWidget(self.chain_selector)
        
        # Refresh button for model list
        refresh_button = QPushButton("↻ Refresh model list")
        refresh_button.clicked.connect(self._refresh_model_list)
        main_layout.addWidget(refresh_button)

        # Run Button
        run_button = QPushButton("ChopChop Foldseek")
        run_button.clicked.connect(self._run_foldseek)
        main_layout.addWidget(run_button)

        main_tab.setLayout(main_layout)

        # Second tab: References (top-aligned, consistent 1-line spacing)
        ref_tab = QWidget()
        ref_layout = QVBoxLayout()
        ref_layout.setContentsMargins(0, 0, 0, 0)
        ref_layout.setSpacing(5)  # consistent single-line spacing
        ref_layout.setAlignment(Qt.AlignTop)

        # Reference section title
        ref_label = QLabel("References:")
        ref_label.setFont(QFont("", weight=QFont.Bold))
        ref_layout.addWidget(ref_label, alignment=Qt.AlignTop)

        # Reference description
        ref_description = QLabel("Foldseek is a rapid structural alignment tool, developed by the Steinegger Lab, that searches and compares protein models to identify similarities in their folds and architectures.")
        ref_description.setWordWrap(True)
        ref_layout.addWidget(ref_description, alignment=Qt.AlignTop)

        # Foldseek Paper Button
        paper_button = QPushButton("Foldseek Paper")
        paper_button.clicked.connect(lambda: QDesktopServices.openUrl(QUrl("https://www.nature.com/articles/s41587-023-01773-0")))
        ref_layout.addWidget(paper_button, alignment=Qt.AlignTop)

        # Steinegger Lab Button
        lab_button = QPushButton("Steinegger Lab")
        lab_button.clicked.connect(lambda: QDesktopServices.openUrl(QUrl("https://steineggerlab.com/en/")))
        ref_layout.addWidget(lab_button, alignment=Qt.AlignTop)

        ref_tab.setLayout(ref_layout)

        # Add tabs to tab widget
        tabs.addTab(main_tab, "Analysis")
        tabs.addTab(ref_tab, "References")

        # Set layout
        container = QWidget()
        container_layout = QVBoxLayout()
        container_layout.setContentsMargins(0, 0, 0, 0)
        container_layout.addWidget(tabs)
        container.setLayout(container_layout)
        self.tool_window.ui_area.setLayout(container_layout)

    def _populate_models_and_chains(self):
        self.model_selector.clear()
        self.chain_selector.clear()
        models, chains = self.select_models_and_chains(self.session)

        for model in models:
            self.model_selector.addItem(model)
        for chain in chains:
            self.chain_selector.addItem(chain)

    def select_models_and_chains(self, session):
        """Fetch all available models and chains."""
        all_models = []
        all_chains = []
        self.session.logger.info("Fetching available models and chains...")
        for model in session.models.list():
            if hasattr(model, "residues"):
                all_models.append(model.id_string)
                for chain_id in set(res.chain_id for res in model.residues):
                    all_chains.append(f"{model.id_string}:{chain_id}")
        self.session.logger.info(f"Found models: {all_models}")
        self.session.logger.info(f"Found chains: {all_chains}")
        return all_models, all_chains
        
    def _refresh_model_list(self):
        """Refresh the dropdown to list currently open models in ChimeraX."""
        print("Refreshing model list...")
        self._populate_models_and_chains()
        print("Model list refreshed.")
        
    # Removed redundant get_models_and_chains, using select_models_and_chains instead

    def _run_foldseek(self):
        model = self.model_selector.currentText()
        chain = self.chain_selector.currentText()
        # NEW: Get the selected database text
        selected_db = self.database_selector.currentText()

        if not model or not chain:
            self.session.logger.error("Model and chain must be selected.")
            return

        # Extract chain ID (e.g., '1:A' -> 'A')
        chain_id = chain.split(':')[1] 

        # NEW: Determine the database parameter based on selection
        if "AlphaFold DB" in selected_db:
            db_param = "database afdb50"
        else:
            # No explicit parameter needed for the default PDB search
            db_param = "" 

        # Construct the final command string
        # .strip() handles case where db_param is empty (PDB default)
        command = f"foldseek #{model} /{chain_id} {db_param}".strip() 
        
        self.session.logger.info(f"Running command: {command}")
        run(self.session, command)

# Register the tool in ChimeraX
from chimerax.core.toolshed import BundleAPI

class _FoldseekAnalysisBundleAPI(BundleAPI):

    @staticmethod
    def start_tool(session, tool_name):
        return FoldseekAnalysis(session, tool_name)

bundle_api = _FoldseekAnalysisBundleAPI()
