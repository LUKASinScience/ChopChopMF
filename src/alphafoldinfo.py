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
AlphaFoldInfoTool.py / Integrated Tool
ChimeraX tool to display AlphaFold resources, color structure by pLDDT score, 
provide UniProt linking and annotation commands, and integrate AlphaSync data analysis.
"""

import re
import webbrowser
import requests
import json

from chimerax.core.tools import ToolInstance
from chimerax.core.commands import run
from chimerax.ui import MainToolWindow
from chimerax.core.models import Model
from Qt.QtWidgets import (
    QVBoxLayout, QLabel, QPushButton, QWidget, QHBoxLayout, QFrame, QTabWidget, 
    QComboBox, QLineEdit, QListWidget, QListWidgetItem,
    QTableWidget, QTableWidgetItem, QHeaderView, QScrollArea
)
from Qt.QtGui import QFont
from Qt.QtCore import Qt


# ---------- AlphaSync Constants ----------
BASE_URL = "https://alphasync.stjude.org/api/v1/"
HEADERS = {"Accept": "application/json"}
# -----------------------------------------


class AlphaFold2(ToolInstance):

    SESSION_ENDURING = False
    SESSION_SAVE = True
    help = "https://www.ebi.ac.uk/training/online/courses/alphafold/"
    
    # --- UniProt Regex for Extraction ---
    UNIPROT_REGEX = re.compile(
        r"""
        (?:
            ^AF\-([A-Z0-9]{6,10})     |      # AF-P04637...
            AlphaFold[\s\-\_]+([A-Z0-9]{6,10})  # AlphaFold Q7LBC6
        )
        """,
        re.VERBOSE | re.IGNORECASE
    )
    # ------------------------------------

    def __init__(self, session, tool_name):
        super().__init__(session, tool_name)

        # Automatically show the tool pane
        run(session, "ui tool show AlphaFold")

        self.display_name = "AlphaFold Info"
        self.tool_window = MainToolWindow(self)
        
        # --- UI elements for existing tabs ---
        self.model_combo = QComboBox()      # Used for pLDDT Coloring tab
        self.chain_combo = QComboBox()      # Used for UniProt tab (chain selector)
        self.uniprot_input = QLineEdit()    # Used for UniProt tab (ID input)
        self.tab_structure_list = QListWidget() # Used for UniProt tab (model list)
        self.warning_label = QLabel()       # Warning label for UniProt tab
        self.tabs = QTabWidget() 
        
        # --- UI elements for AlphaSync tab ---
        self.protein_acc = None             # extracted UniProt ID for AlphaSync
        self.alphasync_input_edit = QLineEdit() # Input for AlphaSync tab
        self.alphasync_status_label = QLabel("") # Status for AlphaSync tab
        self.alphasync_structure_list = QListWidget() # Model list for AlphaSync
        self.tab_residue_table = QTableWidget()     # Table for AlphaSync data


        # Populate initial lists and guess UniProt ID
        self._get_models_and_chains() 
        self._try_guess_uniprot_id()

        self._build_ui()
        self.tool_window.manage('side')
        
    def _extract_uniprot_from_model_name(self, name: str):
        """Extracts UniProt ID using the defined regex."""
        m = self.UNIPROT_REGEX.search(name)
        if not m:
            return None
        return m.group(1) or m.group(2)

    # --- Helper: refresh list of open models and chains ---
    def _get_models_and_chains(self):
        """Updates the internal lists of models and model-chains."""
        self._available_models = []
        self._available_chains = []
        
        for m in self.session.models.list():
            if isinstance(m, Model) and hasattr(m, "residues"):
                model_id_str = m.id_string
                # Label format: #ID Model Name (e.g., #1 AF-P04637-F1)
                label = f"#{model_id_str} {m.name if m.name else '(unnamed)'}"
                self._available_models.append((model_id_str, label))
                
                # Collect chains for the UniProt tab (e.g., #1/A)
                chain_ids = set(res.chain_id for res in m.residues)
                for chain_id in chain_ids:
                    chain_label = f"#{model_id_str}/{chain_id}"
                    self._available_chains.append(chain_label)


    def _refresh_model_list(self):
        """Refresh the list of open models for the pLDDT Coloring tab (QComboBox)."""
        self._get_models_and_chains()
        
        self.model_combo.clear()
        for _, label in self._available_models:
            self.model_combo.addItem(label)
        
        self._refresh_uniprot_selectors()
        self._refresh_alphasync_model_list() 
            
    def _refresh_uniprot_selectors(self):
        """
        Refreshes the structure list (QListWidget) and chain selector (QComboBox) 
        for the UniProt tab.
        """
        self._get_models_and_chains()
        
        # 1. Populate AlphaFold Structure List (QListWidget)
        self.tab_structure_list.clear()
        
        for model_id_str, label in self._available_models:
            model_name = label.split(" ", 1)[-1] 
            
            if self._extract_uniprot_from_model_name(model_name): 
                self.tab_structure_list.addItem(QListWidgetItem(label)) 
        
        # 2. Populate Chain Selector (QComboBox)
        self.chain_combo.clear()
        for chain_label in self._available_chains:
            self.chain_combo.addItem(chain_label)
            
        # 3. Try to guess ID and chain based on the newly filtered list (for initial load)
        self._try_guess_uniprot_id()


    def _try_guess_uniprot_id(self):
        """
        Attempt to pre-fill the UniProt ID and chain from an open AlphaFold model 
        upon tool start or refresh.
        """
        self._get_models_and_chains()

        for model_id_str, label in self._available_models:
            model_name = label.split(" ", 1)[-1] 
            uniprot = self._extract_uniprot_from_model_name(model_name)
            
            if uniprot:
                # Set UniProt ID for UniProt tab
                self.uniprot_input.setText(uniprot.upper())
                # Set UniProt ID for AlphaSync tab
                self.alphasync_input_edit.setText(uniprot.upper())
                self.protein_acc = uniprot.upper()
                
                # Try to auto-select the first chain for the guessed model
                first_chain_spec = next((c for c in self._available_chains if c.startswith(model_id_str + '/')), None)
                if first_chain_spec:
                    index = self.chain_combo.findText(first_chain_spec)
                    if index >= 0:
                        self.chain_combo.setCurrentIndex(index)
                
                self.session.logger.info(f"Guessed UniProt ID: {uniprot.upper()} from model {label}")
                return
        
        self.uniprot_input.setText("")
        self.alphasync_input_edit.setText("")

    def _use_selected_model_for_uniprot(self):
        """
        Action for the 'Use Selected Model' button (UniProt tab): 
        Sets UniProt ID and auto-selects the first chain based on the selected list item.
        """
        item = self.tab_structure_list.currentItem()
        if not item:
            self.session.logger.info("No model selected in the UniProt list.")
            return

        selected_label = item.text() # e.g., "#1 AF-P04637-F1"
        parts = selected_label.split(" ", 1)
        # Check for malformed label, should always have at least two parts
        if len(parts) < 2:
            self.session.logger.error(f"Could not parse model label: {selected_label}")
            return
            
        model_spec = parts[0] # e.g., "#1"
        model_name = parts[1] # e.g., "AF-P04637-F1"
        
        # 1. Set UniProt ID
        uniprot = self._extract_uniprot_from_model_name(model_name)
        if not uniprot:
            self.session.logger.error(f"Could not extract UniProt ID from: {model_name}")
            return
            
        uniprot_upper = uniprot.upper()
        self.uniprot_input.setText(uniprot_upper)
        self.alphasync_input_edit.setText(uniprot_upper) # Sync with AlphaSync tab

        # 2. Auto-select the first chain associated with this model
        first_chain_spec = next((c for c in self._available_chains if c.startswith(model_spec + '/')), None)
        
        if first_chain_spec:
            index = self.chain_combo.findText(first_chain_spec)
            if index >= 0:
                self.chain_combo.setCurrentIndex(index)
                self.session.logger.info(f"Chain selector set to {first_chain_spec}")
            else:
                self.session.logger.warning(f"Chain {first_chain_spec} not found in selector list after refresh.")

        self.session.logger.info(f"UniProt ID and chain set automatically from model selection: {model_spec}")

    def _use_selected_model_for_alphasync(self):
        """
        Action for the 'Use Selected Model' button (AlphaSync tab): 
        Sets UniProt ID for the AlphaSync input field.
        
        FIX: Added parsing for the selected label to remove the leading model ID.
        """
        item = self.alphasync_structure_list.currentItem()
        if not item:
            self._alphasync_log("No model selected in the AlphaSync list.")
            return

        # The label is in the format "#ID Model Name (e.g., #1 AF-P04637-F1)"
        selected_label = item.text() 
        
        # Extract the model name string only (e.g., "AF-P04637-F1")
        parts = selected_label.split(" ", 1)
        if len(parts) < 2:
            self._alphasync_log(f"Could not parse model label: {selected_label}")
            return
            
        model_name = parts[1]

        uniprot = self._extract_uniprot_from_model_name(model_name)

        if not uniprot:
            self._alphasync_log(f"Could not extract UniProt ID from: {model_name}")
            return

        uniprot_upper = uniprot.upper()
        self.protein_acc = uniprot_upper
        self.uniprot_input.setText(uniprot_upper) # Sync with UniProt tab
        self.alphasync_input_edit.setText(uniprot_upper)
        self._alphasync_log(f"UniProt ID set to {self.protein_acc}")
        
    def _build_ui(self):
        self.tabs = QTabWidget() 

        # --- Tab 1: Confidence score info and coloring ---
        main_tab = self._build_main_tab()
        self.tabs.addTab(main_tab, "pLDDT Coloring")
        
        # --- Tab 2: UniProt Tab ---
        uniprot_tab = self._build_uniprot_tab()
        self.tabs.addTab(uniprot_tab, "UniProt")

        # --- Tab 3: Databases (Quick Links) ---
        db_tab = self._build_database_tab()
        self.tabs.addTab(db_tab, "Databases")
        
        # --- Tab 4: AlphaSync (NEW) ---
        alphasync_tab = self._build_alphasync_tab()
        self.tabs.addTab(alphasync_tab, "AlphaSync")


        # Set layout
        container = QWidget()
        layout = QVBoxLayout()
        layout.addWidget(self.tabs)
        container.setLayout(layout)
        self.tool_window.ui_area.setLayout(layout)

    def _build_main_tab(self):
        main_tab = QWidget()
        main_layout = QVBoxLayout()

        title = QLabel("AlphaFold Confidence Score Information")
        title_font = QFont()
        title_font.setPointSize(12)
        title_font.setBold(True)
        title.setFont(title_font)
        main_layout.addWidget(title)

        bold_font = QFont()
        bold_font.setBold(True)
        plddt_label = QLabel("pLDDT (per-residue confidence):")
        plddt_label.setFont(bold_font)
        main_layout.addWidget(plddt_label)
        main_layout.addWidget(QLabel("Confidence level of individual residues in a predicted structure."))

        main_layout.addWidget(QLabel("Color Coding:"))
        legend_layout = QVBoxLayout()
        colors = [
            ("0–50: Very Low", "darkorange"),
            ("50–70: Low", "yellow"),
            ("70–90: Confident", "deepskyblue"),
            ("90–100: Very High", "blue")
        ]
        for label_text, color in colors:
            row = QHBoxLayout()
            box = QFrame()
            box.setFixedSize(20, 20)
            box.setStyleSheet(f"background-color: {color}; border: 1px solid black;")
            row.addWidget(box)
            row.addWidget(QLabel(label_text))
            legend_layout.addLayout(row)
        main_layout.addLayout(legend_layout)

        main_layout.addWidget(QLabel("Select model to color:"))
        main_layout.addWidget(self.model_combo) 
        self._refresh_model_list() 

        refresh_btn = QPushButton("↻ Refresh model list")
        refresh_btn.clicked.connect(self._refresh_model_list)
        main_layout.addWidget(refresh_btn)

        color_btn = QPushButton("Color selected model by AlphaFold2 pLDDT score")
        color_btn.clicked.connect(self._color_selected_model)
        main_layout.addWidget(color_btn)
        
        main_layout.addStretch()
        main_tab.setLayout(main_layout)
        return main_tab

    def _build_uniprot_tab(self):
        uniprot_tab = QWidget()
        uniprot_layout = QVBoxLayout()

        title = QLabel("UniProt Annotation & Association")
        title_font = QFont()
        title_font.setPointSize(12)
        title_font.setBold(True)
        title.setFont(title_font)
        uniprot_layout.addWidget(title)
        uniprot_layout.addWidget(QLabel("1. Select an open AlphaFold model:"))
        
        # --- Structure Selection List (QListWidget) ---
        uniprot_layout.addWidget(self.tab_structure_list)
        
        btn_layout = QHBoxLayout()
        
        refresh_btn = QPushButton("Refresh Model List")
        refresh_btn.clicked.connect(self._refresh_uniprot_selectors)
        btn_layout.addWidget(refresh_btn)

        # Explicit button to set ID and Chain
        use_btn = QPushButton("Use Selected Model")
        use_btn.clicked.connect(self._use_selected_model_for_uniprot)
        btn_layout.addWidget(use_btn)
        
        uniprot_layout.addLayout(btn_layout)

        uniprot_layout.addWidget(QLabel("2. Check/Edit UniProt ID (Pre-set after using model):"))
        uniprot_layout.addWidget(self.uniprot_input)
        
        uniprot_layout.addWidget(QLabel("3. Select Model Chain to Associate (Pre-set after using model):"))
        uniprot_layout.addWidget(self.chain_combo)
        
        # --- WARNING LABEL ---
        self.warning_label = QLabel("⚠️ Please ensure matching UniProt ID with selected model/chain.")
        uniprot_layout.addWidget(self.warning_label)
        # --- END WARNING LABEL ---

        associate_btn = QPushButton("4. Fetch UniProt Annotation & Associate")
        associate_btn.clicked.connect(self._associate_uniprot_annotation)
        uniprot_layout.addWidget(associate_btn)
        
       
        self._refresh_uniprot_selectors() 
        
        uniprot_layout.addStretch()
        uniprot_tab.setLayout(uniprot_layout)
        return uniprot_tab

    def _build_database_tab(self):
        db_tab = QWidget()
        db_layout = QVBoxLayout()

        db_layout.addWidget(QLabel("Quick Links to Databases:"))

        uniprot_btn = QPushButton("UniProt")
        uniprot_btn.clicked.connect(lambda: webbrowser.open("https://www.uniprot.org/"))
        db_layout.addWidget(uniprot_btn)
        
        db_layout.addWidget(QLabel("Developed by Google DeepMind and EMBL-EBI"))

        afdb_btn = QPushButton("AlphaFold Protein Structure Database")
        afdb_btn.clicked.connect(lambda: webbrowser.open("https://alphafold.com/"))
        db_layout.addWidget(afdb_btn)

        db_layout.addStretch()
        db_tab.setLayout(db_layout)
        return db_tab
    
    #-------------------------AlphaSync Tab Logic------------------
  
    
    def _build_alphasync_tab(self):
        alphasync_tool_tabs = QTabWidget()
        
        # -------------------- Main Layout --------------------
        main_widget = QWidget()
        main_layout = QVBoxLayout()
        
        title = QLabel("AlphaSync Data Integration")
        font = QFont()
        font.setPointSize(12)
        font.setBold(True)
        title.setFont(font)
        main_layout.addWidget(title)

        # -------------------- Reference/Browse Buttons --------------------
        button_row_layout = QHBoxLayout()
        
        reference_btn = QPushButton("AlphaSync reference")
        reference_btn.clicked.connect(lambda: webbrowser.open("https://www.nature.com/articles/s41594-025-01719-x"))
        button_row_layout.addWidget(reference_btn)
        
        browse_btn = QPushButton("Browse AlphaSync webside")
        browse_btn.clicked.connect(lambda: webbrowser.open("https://alphasync.stjude.org"))
        button_row_layout.addWidget(browse_btn)

        main_layout.addLayout(button_row_layout)
        
        # -------------------- Input + Fetch Button --------------------
        input_layout = QHBoxLayout()

        input_layout.addWidget(QLabel("UniProt ID:"))
        input_layout.addWidget(self.alphasync_input_edit)

        go_btn = QPushButton("ChopChop AlphaSync Residue Data")
        go_btn.clicked.connect(self._fetch_all_data)
        input_layout.addWidget(go_btn)

        main_layout.addLayout(input_layout)

        # Status label
        main_layout.addWidget(self.alphasync_status_label)

        #--------STRUCTURE SELECTION TAB (AlphaSync Sub-tab)--------------
       
        tab_structure_widget = QWidget()
        tab_structure_layout = QVBoxLayout()
        tab_structure_widget.setLayout(tab_structure_layout)

        label = QLabel("AlphaFold Structures Loaded in ChimeraX")
        font_b = QFont()
        font_b.setBold(True)
        label.setFont(font_b)
        tab_structure_layout.addWidget(label)

        self.alphasync_structure_list.setSelectionMode(QListWidget.SingleSelection)
        tab_structure_layout.addWidget(self.alphasync_structure_list)

        btn_layout = QHBoxLayout()

        # -------------------- REFRESH BUTTON --------------------
        refresh_btn = QPushButton("Refresh Model List")
        refresh_btn.clicked.connect(self._refresh_alphasync_model_list)
        btn_layout.addWidget(refresh_btn)

        # -------------------- USE SELECTED MODEL --------------------
        use_btn = QPushButton("Use Selected Model")
        use_btn.clicked.connect(self._use_selected_model_for_alphasync)
        btn_layout.addWidget(use_btn)

        tab_structure_layout.addLayout(btn_layout)
        alphasync_tool_tabs.addTab(tab_structure_widget, "Structure Selection")
        
        self._refresh_alphasync_model_list()

        #------------------RESIDUE DATA TAB (AlphaSync Sub-tab)-------------
  
        self.tab_residue_table.setColumnCount(9)
        self.tab_residue_table.setHorizontalHeaderLabels([
            "Site", "Residue", "pLDDT", "SASA", "RSA",
            "Surface", "Disorder", "Sec. str.", "Contacts"
        ])
        self.tab_residue_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        alphasync_tool_tabs.addTab(self.tab_residue_table, "Residue Data")

        # -----------------EXPLANATION TAB (AlphaSync Sub-tab)--------------
     
        tab_explanation = QWidget()
        exp_layout = QVBoxLayout()
        tab_explanation.setLayout(exp_layout)

        explanations = {
            "SASA": "The solvent-accessible surface area in Å² (DSSP algorithm).",
            "RSA": "Relative solvent-accessible surface area scaled to percentage (Tien et al, 2013).",
            "Surface": "Surface/core classification. ≤25% RSA = buried (Levy et al., 2010).",
            "Disorder": "An asterisk (*) indicates a disordered residue, while a period (.) indicates a structured residue according to the AlphaFold 2 structure prediction. This value was averaged over a ±10 aa window.)."
        }

        for key, text in explanations.items():
            btn = QPushButton(key)
            btn.setCheckable(True)

            scroll = QScrollArea()
            scroll.setWidgetResizable(True)

            lbl = QLabel(text)
            lbl.setWordWrap(True)
            lbl.setContentsMargins(5, 5, 5, 5)
            scroll.setWidget(lbl)
            scroll.setVisible(False)

            btn.toggled.connect(lambda checked, s=scroll: s.setVisible(checked))

            exp_layout.addWidget(btn)
            exp_layout.addWidget(scroll)

        exp_layout.addStretch()
        alphasync_tool_tabs.addTab(tab_explanation, "Explanation")
        
        main_layout.addWidget(alphasync_tool_tabs)
        main_widget.setLayout(main_layout)
        return main_widget
        
    def _refresh_alphasync_model_list(self):
        """Refreshes the model list specifically for the AlphaSync Structure Selection tab."""
        self.alphasync_structure_list.clear()
        
        for m in self.session.models.list():
            n = m.name

            # Check if the name matches an AlphaFold pattern
            if n.startswith("AF-") or n.lower().startswith("alphafold"):
                if self._extract_uniprot_from_model_name(n):
                    # Use model label format to be consistent with UniProt tab
                    model_id_str = m.id_string
                    label = f"#{model_id_str} {m.name if m.name else '(unnamed)'}"
                    self.alphasync_structure_list.addItem(QListWidgetItem(label))

  
    #------------------------------LOGGING------------------------------
  

    def _alphasync_log(self, msg):
        self.alphasync_status_label.setText(msg)

    
    #-----------------------------DATA FETCHING-------------------------
  

    def _fetch_all_data(self):
        user_input = self.alphasync_input_edit.text().strip()

        if not user_input and not self.protein_acc:
            self._alphasync_log("No UniProt ID provided.")
            return

        if user_input:
            self.protein_acc = user_input

        self._alphasync_log(f"Using UniProt ID: {self.protein_acc}")

        self._fetch_residue_data()

    def _fetch_residue_data(self):
        try:
            url = f"{BASE_URL}protein/{self.protein_acc}"
            self._alphasync_log(f"Fetching Residue Data: {url}")
            r = requests.get(url, headers=HEADERS)
            r.raise_for_status()
            data = r.json()
            self._populate_residue_table(data)
            self._alphasync_log("Residue data successfully loaded.")
        except Exception as e:
            self._alphasync_log(f"Error fetching Residue Data: {e}")
            self.tab_residue_table.clearContents()
            self.tab_residue_table.setRowCount(0)

    #---------------------------TABLE POPULATION-----------------------

    def _populate_residue_table(self, data):
        self.tab_residue_table.setRowCount(len(data))

        for i, res in enumerate(data):
            self.tab_residue_table.setItem(i, 0, QTableWidgetItem(str(res.get("site", ""))))
            self.tab_residue_table.setItem(i, 1, QTableWidgetItem(res.get("aa", "")))

            plddt = res.get("plddt")
            self.tab_residue_table.setItem(i, 2, QTableWidgetItem(f"{plddt:.1f}" if plddt is not None else ""))

            asa = res.get("asa")
            self.tab_residue_table.setItem(i, 3, QTableWidgetItem(f"{asa:.1f} Å²" if asa is not None else ""))

            rsa = res.get("relasa10")
            rsa_str = f"{rsa*100:.0f}%" if isinstance(rsa, (int, float)) else ""
            self.tab_residue_table.setItem(i, 4, QTableWidgetItem(rsa_str))

            self.tab_residue_table.setItem(i, 5, QTableWidgetItem(res.get("surf", "")))
            self.tab_residue_table.setItem(i, 6, QTableWidgetItem(res.get("dis", "")))
            self.tab_residue_table.setItem(i, 7, QTableWidgetItem(res.get("sec", "")))

            contacts = res.get("contacts")
            if isinstance(contacts, list):
                contacts_str = ", ".join(str(c) for c in contacts)
            else:
                contacts_str = str(contacts) if contacts else ""
            self.tab_residue_table.setItem(i, 8, QTableWidgetItem(contacts_str))

        self.tab_residue_table.resizeColumnsToContents()

    #-------------------------Existing Tab Logic------------------------

    def _color_selected_model(self):
        """Select the chosen model and color by AlphaFold pLDDT score."""
        if self.model_combo.count() == 0:
            self.session.logger.info("No models available to color.")
            return

        selected_text = self.model_combo.currentText()
        if not selected_text.startswith("#"):
            self.session.logger.info("No model selected.")
            return

        model_id = selected_text.split()[0]

        try:
            run(self.session, f"color byattribute bfactor {model_id} palette alphafold")
            self.session.logger.info(f"Colored {model_id} by AlphaFold pLDDT score.")
        except Exception as e:
            self.session.logger.warning(f"Failed to color {model_id}: {e}")
            
    def _associate_uniprot_annotation(self):
        """Run the ChimeraX command to fetch and associate UniProt data."""
        uniprot_id = self.uniprot_input.text().strip().upper()
        chain_spec = self.chain_combo.currentText()

        # 1. Input checks
        if not uniprot_id:
            self.session.logger.error("UniProt ID is required.")
            return

        if not chain_spec:
            self.session.logger.error("A model chain must be selected.")
            return

        # 2. Execute Command
        command = f"open {uniprot_id} fromDatabase uniprot associate {chain_spec}"
        
        self.session.logger.info(f"Running UniProt command: {command}")
        
        try:
            run(self.session, command)
            
            self.session.logger.info(f"Successfully requested UniProt annotation for {uniprot_id} and association with {chain_spec}.")
        except Exception as e:
            self.session.logger.error(f"Failed to run UniProt command: {e}")


# --- Register the tool ---
from chimerax.core.toolshed import BundleAPI


class _AlphaFoldInfoBundleAPI(BundleAPI):
    @staticmethod
    def start_tool(session, tool_name):
        return AlphaFold2(session, tool_name)


bundle_api = _AlphaFoldInfoBundleAPI()