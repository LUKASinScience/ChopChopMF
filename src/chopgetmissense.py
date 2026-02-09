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
chopget_missense.py
A ChimeraX tool to fetch AlphaMissense data, process TSV files, and generate a defattr file for coloring.

Features:
- Automatic download and extraction of hotspot data for a given UniProt ID.
- Automatic chain length detection from the downloaded PDB file.
- Option to upload a custom TSV file; if chosen, a user can pick which currently open ChimeraX model (PDB) to use for chain-length determination via a dropdown.
- A "Refresh model list" button allows updating the model list if new models are opened in ChimeraX.
"""

from chimerax.core.tools import ToolInstance
from chimerax.core.commands import run
from pathlib import Path
from chimerax.ui import MainToolWindow
from Qt.QtWidgets import (
    QVBoxLayout, QLabel, QLineEdit, QPushButton, QWidget, QHBoxLayout, QRadioButton,
    QFileDialog, QFrame, QComboBox, QTabWidget
)
from Qt.QtCore import Qt
import webbrowser
import os
import requests
import zipfile
import pandas as pd


class ChopChopGetMissense(ToolInstance):

    SESSION_ENDURING = False
    SESSION_SAVE = True
    help = "help:user/tools/chopchop_get_missense.html"

    def __init__(self, session, tool_name):
        super().__init__(session, tool_name)
        self.display_name = "AlphaMissense"
        self.tool_window = MainToolWindow(self)
        self._build_ui()
        self.tool_window.manage('side')

    def _build_ui(self):
        tabs = QTabWidget()

        # First tab: Main functionality
        main_tab = QWidget()
        main_layout = QVBoxLayout()

        # UniProt ID
        self.uniprot_input = QLineEdit()
        main_layout.addWidget(QLabel("Enter UniProt ID of Human Protein:"))
        main_layout.addWidget(self.uniprot_input)

        # Download folder path
        self.download_path_input = QLineEdit()
        self.download_path_input.setText(str(Path.home() / "Downloads"))
        main_layout.addWidget(QLabel("Enter download folder path for fetched Files and Attribute file for Scoring Scheme:"))
        main_layout.addWidget(self.download_path_input)

        # Option to upload custom TSV
        self.use_uploaded_tsv_checkbox = QRadioButton("Use uploaded AlphaMissense TSV file")
        self.use_uploaded_tsv_checkbox.setChecked(False)
        main_layout.addWidget(self.use_uploaded_tsv_checkbox)

        # TSV file path
        self.file_path_input = QLineEdit()
        self.file_path_input.setPlaceholderText("Choose a TSV file")
        self.file_path_input.setDisabled(True)
        main_layout.addWidget(QLabel("Upload AlphaMissense TSV file:"))
        main_layout.addWidget(self.file_path_input)

        browse_button = QPushButton("Browse")
        browse_button.clicked.connect(self.browse_file)
        browse_button.setDisabled(True)
        main_layout.addWidget(browse_button)

        # Model selector for chain length
        self.model_selector = QComboBox()
        self.model_selector.setDisabled(True)
        main_layout.addWidget(QLabel("For uploaded TSV: Select open PDB model for chain length:"))
        main_layout.addWidget(self.model_selector)

        # Refresh button for model list
        refresh_button = QPushButton("↻ Refresh model list")
        refresh_button.clicked.connect(self._refresh_model_list)
        refresh_button.setDisabled(True)
        main_layout.addWidget(refresh_button)

        # Toggle function for enabling/disabling relevant widgets
        def toggle_inputs():
            state = self.use_uploaded_tsv_checkbox.isChecked()
            self.file_path_input.setDisabled(not state)
            browse_button.setDisabled(not state)
            self.model_selector.setDisabled(not state)
            refresh_button.setDisabled(not state)

        self.use_uploaded_tsv_checkbox.toggled.connect(toggle_inputs)

        # Execute button
        execute_button = QPushButton("ChopChop Missense PDB")
        execute_button.clicked.connect(self.run_script)
        main_layout.addWidget(execute_button)

        # Color bar legend
        main_layout.addWidget(QLabel("Scoring Scheme:"))
        color_bar = QVBoxLayout()
        colors = [
            ("1: Benign                  AlphaMissense Score: >0.34", "navy"),
            ("2: Ambiguous           AlphaMissense Score: 0.34-0.564", "cornflowerblue"),
            ("3: Pathogenic low     AlphaMissense Score: 0.564-0.78", "red"),
            ("4: Pathogenic high    AlphaMissense Score: 0.78-1.0", "darkred"),
            ("5: No score", "white")
        ]
        for label, color in colors:
            row = QHBoxLayout()
            box = QFrame()
            box.setFixedSize(20, 20)
            box.setStyleSheet(f"background-color: {color}; border: 1px solid black;")
            row.addWidget(box)
            row.addWidget(QLabel(label))
            color_bar.addLayout(row)
        main_layout.addLayout(color_bar)

        main_tab.setLayout(main_layout)

        # Second tab: References
        ref_tab = QWidget()
        ref_layout = QVBoxLayout()
        ref_layout.setContentsMargins(0, 0, 0, 0)
        ref_layout.setSpacing(5)
        ref_layout.setAlignment(Qt.AlignTop)

        ref_layout.addWidget(QLabel("AlphaMissense PDB and Scores are fetched from"))
        hegelab_button = QPushButton("Hegelab AlphaMissense")
        hegelab_button.clicked.connect(lambda: webbrowser.open("https://alphamissense.hegelab.org/hotspot"))
        ref_layout.addWidget(hegelab_button)

        ref_layout.addWidget(QLabel("Reference Google DeepMind's AlphaMissense"))
        alphamissense_button = QPushButton("AlphaMissense Paper")
        alphamissense_button.clicked.connect(lambda: webbrowser.open("https://www.science.org/doi/10.1126/science.adg7492"))
        ref_layout.addWidget(alphamissense_button)

        ref_tab.setLayout(ref_layout)

        # Add tabs
        tabs.addTab(main_tab, "Analysis")
        tabs.addTab(ref_tab, "References")

        container = QWidget()
        container_layout = QVBoxLayout()
        container_layout.addWidget(tabs)
        container.setLayout(container_layout)
        self.tool_window.ui_area.setLayout(container_layout)
    

    def browse_file(self):
        file_path, _ = QFileDialog.getOpenFileName(None, "Select AlphaMissense TSV File", str(Path.home()), "TSV Files (*.tsv)")
        if file_path:
            self.file_path_input.setText(file_path)

    def _refresh_model_list(self):
        """Refresh the dropdown to list currently open models in ChimeraX."""
        self.model_selector.clear()
        # List all open models that contain atomic data
        for model in self.session.models.list():
            if getattr(model, "atoms", None):
                label = f"Model {model.id} - {model.name}"
                # Store the actual model object in the QComboBox item
                self.model_selector.addItem(label, model)

    def run_script(self):
        download_dir = Path(self.download_path_input.text().strip())

        if not download_dir.exists():
            try:
                download_dir.mkdir(parents=True, exist_ok=True)
            except Exception as e:
                self.session.logger.warning(f"Failed to create directory {download_dir}: {e}")
                return

        uniprot_id = self.uniprot_input.text().strip()
        if not uniprot_id:
            self.session.logger.warning("Please provide a UniProt ID.")
            return

        self.session.logger.info(f"Starting download for UniProt ID: {uniprot_id}")
        zip_file_path, extracted_folder = self.download_hotspot_file(uniprot_id, download_dir)

        if zip_file_path and extracted_folder:
            # Try to get the TSV file from the extracted folder
            downloaded_tsv = next(extracted_folder.glob("*.tsv"), None)

            if self.use_uploaded_tsv_checkbox.isChecked():
                # Use the user-uploaded TSV
                tsv_file = Path(self.file_path_input.text().strip())
                if not tsv_file.exists():
                    self.session.logger.warning("Provided TSV file does not exist. Aborting.")
                    return

                # Determine chain length from the selected open model
                model_idx = self.model_selector.currentIndex()
                if model_idx < 0:
                    self.session.logger.warning("No open PDB model selected. Aborting.")
                    return

                open_model = self.model_selector.itemData(model_idx)
                if not open_model:
                    self.session.logger.warning("Selected model is invalid. Aborting.")
                    return

                chain_length = self.get_chain_length_from_model(open_model)
                self.session.logger.info(f"Determined chain length from open model: {chain_length}")

                # Generate the defattr file
                self.process_uploaded_tsv(tsv_file, extracted_folder / "AM_score.defattr", chain_length)
            else:
                # Automatic mode: use downloaded TSV, parse PDB from extracted folder
                if not downloaded_tsv:
                    self.session.logger.warning("No TSV file found in the extracted folder.")
                    return
                tsv_file = downloaded_tsv
                pdb_files = list(extracted_folder.glob("*.pdb"))
                if not pdb_files:
                    self.session.logger.warning("No PDB files found in the extracted folder.")
                    return
                pdb_file = pdb_files[0]
                chain_length = self.get_chain_length_from_pdb_file(pdb_file)
                self.session.logger.info(f"Determined chain length from downloaded PDB: {chain_length}")

                # Generate the defattr file
                self.generate_defattr_file(tsv_file, extracted_folder / "AM_score.defattr", chain_length)

            # Step 4: Load PDB files into ChimeraX
            pdb_files = list(extracted_folder.glob("*.pdb"))
            if not pdb_files:
                self.session.logger.warning("No PDB files found in the extracted folder.")
                return
            for pdb_file in pdb_files:
                run(self.session, f"open {pdb_file}")

            # Step 5: Apply coloring using defattr file
            run(self.session, f"open {extracted_folder / 'AM_score.defattr'}")
            self.session.logger.info("Applied coloring using AM_score.defattr")

            # Step 6: Apply coloring with palette
            run(self.session, "color byattribute AM_score palette 1,navy:2,cornflowerblue:3,red:4,darkred:5,white")
            self.session.logger.info("Applied color byattribute with custom palette")
            self.session.logger.info(f"Process completed successfully for UniProt ID {uniprot_id}")

    def download_hotspot_file(self, uniprot_id, download_dir):
        base_url = "https://alphamissense.hegelab.org/download_hot"
        download_url = f"{base_url}?uniprot_id={uniprot_id}"
        zip_file_path = download_dir / f"{uniprot_id}_hotspots.zip"
        extracted_folder = download_dir / f"{uniprot_id}_hotspots"

        try:
            self.session.logger.info(f"Downloading hotspot dataset from {download_url}...")
            response = requests.get(download_url)
            response.raise_for_status()
            with open(zip_file_path, "wb") as f:
                f.write(response.content)
            self.session.logger.info(f"Downloaded to {zip_file_path}")

            with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
                zip_ref.extractall(extracted_folder)
            self.session.logger.info(f"Extracted to {extracted_folder}")
            return zip_file_path, extracted_folder
        except Exception as e:
            self.session.logger.warning(f"Failed during processing: {e}")
            return None, None

    def process_uploaded_tsv(self, tsv_file, defattr_file, chain_length):
        """Process the uploaded TSV file to generate a defattr file."""
        try:
            self.session.logger.info(f"Processing uploaded TSV file: {tsv_file}")
            df = pd.read_csv(tsv_file, sep="\t", engine='python')
            position_col = next((col for col in df.columns if "position" in col.lower()), None)
            score_col = next((col for col in df.columns if "mean" in col.lower() and "am score" in col.lower()), None)
            if not position_col or not score_col:
                raise ValueError("Required columns not found in the TSV file.")
            df = df[[position_col, score_col]]
            df.columns = ["position", "mean AM score from all"]
            self.session.logger.info(f"TSV file loaded with {len(df)} rows and {len(df.columns)} columns.")
            self.generate_defattr_file_from_dataframe(df, defattr_file, chain_length)
        except Exception as e:
            self.session.logger.warning(f"Failed to process uploaded TSV file: {e}")

    def generate_defattr_file(self, tsv_file, defattr_file, chain_length):
        """Generate a defattr file for coloring based on AM_score from a TSV file."""
        try:
            df = pd.read_csv(tsv_file, sep="\t", usecols=["position", "mean AM score from all"], engine='python')
            self.generate_defattr_file_from_dataframe(df, defattr_file, chain_length)
        except Exception as e:
            self.session.logger.warning(f"Error generating defattr file from TSV: {e}")

    def generate_defattr_file_from_dataframe(self, df, defattr_file, max_residue_number):
        """Generate a defattr file from a DataFrame using residues 1..max_residue_number.
        If a residue is missing in the TSV, a default score of 5 is assigned.
        """
        ranges = {
            "benign": (0, 0.34, 1),
            "ambiguous": (0.34, 0.564, 2),
            "patho0": (0.564, 0.78, 3),
            "patho1": (0.78, 1.0, 4)
        }
        header = [
            "# Coloring scheme:",
            "# 1: Benign (navy)",
            "# 2: Ambiguous (cornflowerblue)",
            "# 3: Pathogenic low (red)",
            "# 4: Pathogenic high (darkred)",
            "# 5: No score (white)"
        ]
        lines = header + ["attribute: AM_score", "match mode: 1-to-1", "recipient: residues"]

        residue_scores = {}
        for index, row in df.iterrows():
            residue_id = int(row["position"])
            score = float(row["mean AM score from all"])
            # Assign a default of 5, then override if it falls in a known range
            assigned_val = 5
            for category, (low, high, assigned_score) in ranges.items():
                if low <= score < high:
                    assigned_val = assigned_score
                    break
            residue_scores[residue_id] = assigned_val

        # Fill up from 1..max_residue_number
        for residue_id in range(1, max_residue_number + 1):
            assigned_score = residue_scores.get(residue_id, 5)
            lines.append(f"\t:{residue_id}\t{assigned_score}")
        lines.append("")

        with open(defattr_file, "w") as f:
            f.write("\n".join(lines))
        self.session.logger.info(f"Generated defattr file: {defattr_file}")

    def get_chain_length_from_pdb_file(self, pdb_file):
        """Parse the PDB file to determine the maximum residue number (chain length)."""
        max_res = 0
        try:
            with open(pdb_file, "r") as f:
                for line in f:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        try:
                            resnum = int(line[22:26].strip())
                            if resnum > max_res:
                                max_res = resnum
                        except:
                            continue
        except Exception as e:
            self.session.logger.warning(f"Error parsing PDB file {pdb_file}: {e}")
        return max_res

    def get_chain_length_from_model(self, model):
        """Determine the maximum residue number from an open ChimeraX model."""
        max_res = 0
        try:
            # Assuming model.atoms is iterable and each atom has 'resid'
            for atom in model.atoms:
                if atom.resid > max_res:
                    max_res = atom.resid
        except Exception as e:
            self.session.logger.warning(f"Error reading model {model.id}: {e}")
        return max_res



