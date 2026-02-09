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
MissenseAlignment.py
A ChimeraX tool for fetching AlphaMissense data (via a UniProt ID), aligning sequences,
and coloring structures based on the alignment. The tool uses the AlphaMissense scores
from the human protein (from the TSV file) and maps them onto the selected model chain.
Only residues that perfectly match the human sequence (i.e. the same residue letter)
receive their corresponding score; all other residues (mismatches or gaps) are assigned
the default "no score" color (yellow).
User Interface:
 - A dropdown to select "modelID:chainID" for alignment.
 - A QLineEdit for entering the human protein UniProt ID.
 - A button labeled "ChopChop Missense" to run the process.
"""

from chimerax.core.tools import ToolInstance
from chimerax.core.commands import run
from pathlib import Path
from chimerax.ui import MainToolWindow
from Qt.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, QComboBox, QWidget, QFrame
)
import os
import requests
import pandas as pd
import zipfile

class ChopChopMissense(ToolInstance):
    """
    A ChimeraX tool for fetching AlphaMissense data, aligning sequences, and coloring structures.
    The alignment remaps scores from the human protein (from the TSV) onto the chosen model chain.
    Only residues that perfectly match the human sequence get their score; all other residues receive
    the default "no score" color (yellow).
    """
    SESSION_ENDURING = False
    SESSION_SAVE = True

    def __init__(self, session, tool_name):
        super().__init__(session, tool_name)
        self.display_name = "Missense"
        self.tool_window = MainToolWindow(self)
        self._build_ui()
        self.tool_window.manage('side')

    def _build_ui(self):
        # Collect open models and chains using the initial method.
        all_models, all_chains = self._select_models_and_chains(self.session)

        layout = QVBoxLayout()

        # Dropdown for selecting "modelID:chainID"
        self.combobox = QComboBox()
        for chain in all_chains:
            self.combobox.addItem(chain)
        layout.addWidget(QLabel("Choose model ID and chain ID for alignment of the non Human Protein"))
        layout.addWidget(self.combobox)
        
        # Refresh button for model list
        refresh_button = QPushButton("↻ Refresh model list")
        refresh_button.clicked.connect(self._refresh_model_list)
        layout.addWidget(refresh_button)

        # Input for the human protein UniProt ID
        self.uniprot_input = QLineEdit()
        layout.addWidget(QLabel("Enter Human Protein UniProt ID to fetch AlphaMissense data"))
        layout.addWidget(self.uniprot_input)

        # Button to execute the process
        execute_button = QPushButton("ChopChop Missense Alignment")
        layout.addWidget(execute_button)
        execute_button.clicked.connect(self.fetch_align_and_color)
        
        # Color bar legend
        layout.addWidget(QLabel("Scoring Scheme:"))
        color_bar = QVBoxLayout()

        colors = [
            ("1: Benign                     AlphaMissense Score: >0.34", "navy"),
            ("2: Ambiguous             AlphaMissense Score: 0.34-0.564", "cornflowerblue"),
            ("3: Pathogenic low      AlphaMissense Score: 0.564-0.78", "red"),
            ("4: Pathogenic high    AlphaMissense Score: 0.78-1.0", "darkred"),
            ("5: No score               no Alignment to Human Protein", "yellow")
        ]

        for label_text, color in colors:
            color_item = QHBoxLayout()
            color_box = QFrame()
            color_box.setFixedSize(20, 20)
            color_box.setStyleSheet(f"background-color: {color}; border: 1px solid black;")
            color_label = QLabel(label_text)
            color_item.addWidget(color_box)
            color_item.addWidget(color_label)
            color_bar.addLayout(color_item)

        layout.addLayout(color_bar)

        container = QWidget()
        container.setLayout(layout)
        self.tool_window.ui_area.setLayout(layout)

    def _select_models_and_chains(self, session):
        """
        Collects all open models and their chains.
        Returns two lists:
          - List of model id strings.
          - List of "modelID:chainID" strings.
        """
        self.session.logger.info("Collecting open models and chains.")
        all_models = []
        all_chains = []

        for model in session.models:
            model_id_str = model.id_string
            all_models.append(model_id_str)
            if hasattr(model, 'chains'):
                for chain in model.chains:
                    chain_id = chain.chain_id
                    all_chains.append(f"{model_id_str}:{chain_id}")
        self.session.logger.info(f"Found models: {all_models}, chains: {all_chains}")
        return all_models, all_chains

    def _refresh_model_list(self):
        """Refresh the dropdown to list currently open models in ChimeraX."""
        self.populate_combobox()

    def populate_combobox(self):
        """Fetch available models and chains and populate the dropdown."""
        all_models, all_chains = self.get_models_and_chains(self.session)
        self.combobox.clear()
        for chain in all_chains:
            self.combobox.addItem(chain)

    def get_models_and_chains(self, session):
        """Retrieve all available models and chains using the same simple format."""
        all_models = []
        all_chains = []
        for model in session.models.list():
            if hasattr(model, "residues"):
                all_models.append(model.id_string)
                for chain_id in set(res.chain_id for res in model.residues):
                    # Use the simple format: "modelID:chainID"
                    all_chains.append(f"{model.id_string}:{chain_id}")
        return all_models, all_chains

    def fetch_align_and_color(self):
        """
        Main workflow:
         1. Retrieve the selected model/chain and human UniProt ID.
         2. Extract the chain sequence from ChimeraX.
         3. Download the AlphaMissense hotspot data and extract the human protein sequence from the TSV.
         4. Perform sequence alignment using Muscle.
         5. Parse the alignment file to obtain mappings of residues for both sequences.
         6. Create a defattr file that maps the human AlphaMissense scores to the model residues.
            Only if the aligned residues are identical is the corresponding score assigned.
            Otherwise, the default "no score" color (yellow, code 5) is used.
         7. Apply coloring in ChimeraX.
        """
        self.session.logger.info("Starting fetch, align, and color process.")
        selected_chain = self.combobox.currentText()
        uniprot_id = self.uniprot_input.text().strip()

        if not selected_chain or not uniprot_id:
            self.session.logger.warning("Please provide both a model/chain and a UniProt ID.")
            return

        try:
            model_id, chain_id = selected_chain.split(":")
        except ValueError:
            self.session.logger.warning("Invalid chain selection format. Please re-check.")
            return

        self.session.logger.info(f"Selected model: {model_id}, chain: {chain_id}")

        # 1. Extract model chain sequence from ChimeraX
        model_sequence = self.get_chain_sequence_as_string(model_id, chain_id)
        if not model_sequence:
            self.session.logger.warning("Failed to get chain sequence.")
            return

        # 2. Download hotspot data and extract human protein sequence from TSV
        zip_file_path, extracted_folder = self.download_hotspot_file(uniprot_id, Path.home() / "Downloads")
        if not extracted_folder:
            self.session.logger.warning("Failed to download and extract AlphaMissense data.")
            return

        tsv_file = extracted_folder / f"AlphaMissense-Hotspot-{uniprot_id}.tsv"
        human_sequence = self.extract_tsv_sequence_as_string(tsv_file)
        if not human_sequence:
            self.session.logger.warning("Failed to extract sequence from AlphaMissense TSV.")
            return

        # 3. Perform alignment between the model sequence and the human protein sequence
        alignment_file = self.perform_alignment(model_sequence, human_sequence, extracted_folder)
        if not alignment_file:
            self.session.logger.warning("Failed to perform sequence alignment.")
            return

        # 4. Parse the alignment file to obtain mappings for both sequences
        model_aln, human_aln, model_res_numbers, human_res_numbers = self.parse_alignment_with_residues(alignment_file)
        if model_aln is None or human_aln is None:
            self.session.logger.warning("Failed to parse alignment file.")
            return

        # 5. Map human scores to model residues only if the residues match perfectly.
        self.color_model(model_aln, human_aln, model_res_numbers, human_res_numbers, tsv_file, extracted_folder)

    def get_chain_sequence_as_string(self, model_id, chain_id):
        """
        Uses ChimeraX commands to display the sequence for the specified chain,
        saves it as FASTA, and returns the sequence string (without header).
        """
        self.session.logger.info(f"Fetching sequence for model {model_id}, chain {chain_id}.")
        output_dir = Path("/tmp")
        output_dir.mkdir(parents=True, exist_ok=True)
        sequence_path = output_dir / "missense_sequence1.fasta"

        try:
            run(self.session, f"sequence chain #{model_id}/{chain_id}")
            run(self.session, f"save {sequence_path} format fasta")
            self.session.logger.info(f"Sequence saved to {sequence_path}.")

            with open(sequence_path, "r", encoding="utf-8-sig") as file:
                lines = file.readlines()

            seq_lines = [line.strip() for line in lines if not line.startswith(">")]
            sequence = "".join(seq_lines).replace("\n", "").replace("\r", "")
            return sequence

        except Exception as e:
            self.session.logger.warning(f"Error fetching chain sequence: {e}")
            return None

    def download_hotspot_file(self, uniprot_id, download_dir):
        """
        Downloads the AlphaMissense hotspot dataset (ZIP file) for the given UniProt ID,
        extracts it, and returns the ZIP file path and the extraction folder.
        """
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

    def extract_tsv_sequence_as_string(self, tsv_file):
        """
        Extracts the complete amino acid sequence from the TSV file.
        Checks for the column 'a.a' or 'a.a.'.
        """
        self.session.logger.info(f"Extracting sequence from AlphaMissense TSV file: {tsv_file}.")
        try:
            df = pd.read_csv(tsv_file, sep="\t", engine="python")
            possible_cols = ["a.a", "a.a."]
            sequence_column = None
            for col in possible_cols:
                if col in df.columns:
                    sequence_column = col
                    break
            if not sequence_column:
                raise ValueError("Could not find 'a.a' or 'a.a.' column in the TSV file.")
            sequence = "".join(df[sequence_column].dropna().astype(str))
            sequence = sequence.replace("-", "").strip()
            return sequence

        except Exception as e:
            self.session.logger.warning(f"Error extracting sequence from TSV: {e}")
            return None

    def perform_alignment(self, sequence1, sequence2, extracted_folder):
        """
        Uses ChimeraX to align the two sequences using MUSCLE and saves the alignment to a FASTA file.
        Returns the alignment file path.
        """
        self.session.logger.info("Performing alignment for sequences.")
        alignment_file_path = extracted_folder / "missense_alignment.fasta"

        try:
            run(self.session, f"sequence align {sequence1},{sequence2} program muscle")
            run(self.session, f"save {alignment_file_path} format fasta alignment 1")
            self.session.logger.info(f"Alignment completed and saved to {alignment_file_path}.")
            return alignment_file_path

        except Exception as e:
            self.session.logger.warning(f"Error during alignment: {e}")
            return None

    def parse_alignment_with_residues(self, alignment_file):
        """
        Parses the FASTA alignment file (assumed to contain exactly two sequences).
        Returns:
          - model_aln: aligned model sequence (string)
          - human_aln: aligned human sequence (string)
          - model_res_numbers: list mapping each alignment column to the model residue number (or None if gap)
          - human_res_numbers: list mapping each alignment column to the human residue number (or None if gap)
        """
        self.session.logger.info(f"Parsing alignment file: {alignment_file}.")
        try:
            with open(alignment_file, "r") as f:
                lines = f.readlines()

            sequences = []
            current_seq = ""
            for line in lines:
                if line.startswith(">"):
                    if current_seq:
                        sequences.append(current_seq)
                        current_seq = ""
                else:
                    current_seq += line.strip()
            if current_seq:
                sequences.append(current_seq)

            if len(sequences) < 2:
                self.session.logger.warning("Alignment file does not contain two sequences.")
                return None, None, None, None

            # Assume first sequence is from the model and second from the human protein
            model_aln = sequences[0]
            human_aln = sequences[1]

            model_res_numbers = []
            human_res_numbers = []
            model_counter = 1
            human_counter = 1
            aln_length = len(model_aln)
            for i in range(aln_length):
                if model_aln[i] == "-":
                    model_res_numbers.append(None)
                else:
                    model_res_numbers.append(model_counter)
                    model_counter += 1
                if human_aln[i] == "-":
                    human_res_numbers.append(None)
                else:
                    human_res_numbers.append(human_counter)
                    human_counter += 1

            self.session.logger.info("Alignment parsed successfully with residue mappings.")
            return model_aln, human_aln, model_res_numbers, human_res_numbers

        except Exception as e:
            self.session.logger.warning(f"Error parsing alignment file: {e}")
            return None, None, None, None

    def color_model(self, model_aln, human_aln, model_res_numbers, human_res_numbers, tsv_file, extracted_folder):
        """
        Creates a .defattr file that maps the AlphaMissense scores (from the human protein)
        to the corresponding model chain residues using the alignment.
        For each alignment column where both sequences have a residue and the residue letters match exactly,
        the human residue number is used to look up its score (from the TSV), and that score is assigned to
        the corresponding model residue number. If the residue letters differ or a gap is present,
        the default "no score" color (yellow, code 5) is assigned.
        """
        self.session.logger.info("Coloring model based on alignment and AlphaMissense data.")
        try:
            df = pd.read_csv(tsv_file, sep="\t", engine="python")
            score_col = "mean AM score from all"
            if score_col not in df.columns:
                raise ValueError(f"Column '{score_col}' not found in TSV file.")
            scores = df[score_col].tolist()  # human residue 1 corresponds to scores[0]

            # Build a mapping from model residue number to the human score,
            # but only if the aligned residue letters match exactly.
            model_score_mapping = {}
            for i in range(len(model_aln)):
                m_res = model_res_numbers[i]
                h_res = human_res_numbers[i]
                if m_res is not None and h_res is not None and model_aln[i] == human_aln[i]:
                    if 1 <= h_res <= len(scores):
                        model_score_mapping[m_res] = scores[h_res - 1]

            model_length = max([num for num in model_res_numbers if num is not None])

            attribute_lines = [
                "attribute: MissenseScores",
                "match mode: 1-to-1",
                "recipient: residues"
            ]
            for r in range(1, model_length + 1):
                score = model_score_mapping.get(r, None)
                if score is None:
                    color_code = 5  # Default "no score" color (yellow)
                else:
                    if score > 0.78:
                        color_code = 4
                    elif score > 0.564:
                        color_code = 3
                    elif score > 0.34:
                        color_code = 2
                    else:
                        color_code = 1
                attribute_lines.append(f"\t:{r}\t{color_code}")

            attribute_file = extracted_folder / "MissenseScores.defattr"
            with open(attribute_file, "w") as f:
                f.write("\n".join(attribute_lines))
                f.write("\n")
            self.session.logger.info(f"Attributes saved to {attribute_file}.")

            run(self.session, f"open {attribute_file}")
            run(self.session, "color byattribute MissenseScores palette 1,navy:2,cornflowerblue:3,firebrick:4,darkred:5,yellow")
            self.session.logger.info("Coloring completed successfully.")

        except Exception as e:
            self.session.logger.warning(f"Error during coloring: {e}")
