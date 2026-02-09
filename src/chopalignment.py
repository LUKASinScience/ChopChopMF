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

from chimerax.core.tools import ToolInstance
from chimerax.core.commands import run
from pathlib import Path
from chimerax.ui import MainToolWindow
from Qt.QtWidgets import (
    QLineEdit, QComboBox, QWidget, QVBoxLayout, QLabel,
    QPushButton, QHBoxLayout, QFileDialog, QFrame
)
from Qt.QtGui import QAction
import csv
import os
import requests
import distutils.dir_util
import shutil


class ChopChopMFalignment(ToolInstance):
    SESSION_ENDURING = False
    SESSION_SAVE = True
    help = "help:user/tools/tutorial.html"

    def __init__(self, session, tool_name):
        super().__init__(session, tool_name)
        self.display_name = "Sequence"
        self.tool_window = MainToolWindow(self)
        self.tool_window.fill_context_menu = self.fill_context_menu
        self._build_ui()

    def _build_ui(self):
        all_models, all_chains, all_modl_index, all_chains_with_index = self.select_models_and_chains(self.session)

        layout = QVBoxLayout()
        self.combobox = QComboBox()
        for x in all_chains:
            model_id, chain_id = x.split(":")
            self.combobox.addItem(f"Model {model_id}, Chain {chain_id}")

        layout.addWidget(QLabel("Choose model ID and chain ID for alignment"))
        layout.addWidget(self.combobox)

        refresh_button = QPushButton("↻ Refresh model list")
        refresh_button.clicked.connect(self._refresh_model_list)
        layout.addWidget(refresh_button)

        self.text1 = QLineEdit()
        layout.addWidget(QLabel("Enter an amino acid sequence or UniProt ID to be aligned to your model"))
        layout.addWidget(self.text1)

        button1 = QPushButton("ChopChop SequenceAlignment")
        layout.addWidget(button1)

        self.download_path_input = QLineEdit()
        self.download_path_input.setText(str(Path.home() / "Downloads"))
        layout.addWidget(QLabel("Enter download folder path:"))
        layout.addWidget(self.download_path_input)

        browse_button = QPushButton("Browse")
        browse_button.clicked.connect(self.browse_folder)
        layout.addWidget(browse_button)

        button2 = QPushButton("Download")
        layout.addWidget(button2)

        button1.clicked.connect(self.align)
        button2.clicked.connect(self.save)

        layout.addWidget(QLabel("Scoring Scheme:"))
        color_bar = QVBoxLayout()
        colors = [
            ("conserved", "teal"),
            ("semi-conserved", "palegoldenrod"),
            ("not conserved", "lightsalmon"),
            ("gap", "dimgrey"),
        ]
        for label, color in colors:
            color_item = QHBoxLayout()
            color_box = QFrame()
            color_box.setFixedSize(20, 20)
            color_box.setStyleSheet(f"background-color: {color}; border: 1px solid black;")
            color_label = QLabel(label)
            color_item.addWidget(color_box)
            color_item.addWidget(color_label)
            color_bar.addLayout(color_item)
        layout.addLayout(color_bar)

        layout.addWidget(QLabel("Customize Colors after running ChopChop SequenceAlignment if needed:"))
        self.color_inputs = {}
        color_entries = [
            ("1 (Conserved):", "teal"),
            ("2 (Semi-Conserved):", "palegoldenrod"),
            ("3 (Not Conserved):", "lightsalmon"),
            ("4 (Gap):", "dimgrey"),
        ]
        for label, default_color in color_entries:
            entry_layout = QHBoxLayout()
            color_label = QLabel(label)
            color_input = QLineEdit(default_color)
            entry_layout.addWidget(color_label)
            entry_layout.addWidget(color_input)
            layout.addLayout(entry_layout)
            self.color_inputs[label.split(" ")[0]] = color_input

        apply_color_button = QPushButton("Apply New Color Scheme")
        apply_color_button.clicked.connect(self.apply_custom_colors)
        layout.addWidget(apply_color_button)

        container = QWidget()
        container.setLayout(layout)
        self.tool_window.ui_area.setLayout(layout)
        self.tool_window.manage('side')

    def apply_custom_colors(self):
        palette = []
        for key, input_field in self.color_inputs.items():
            color = input_field.text().strip()
            palette.append(f"{key},{color}")
        color_command = "color byattribute MUSCLEalignment palette " + ":".join(palette)
        try:
            run(self.session, color_command)
            print("Custom color scheme applied:", color_command)
        except Exception as e:
            print("Error applying custom color scheme:", str(e))

    def browse_folder(self):
        folder_path = QFileDialog.getExistingDirectory(self.tool_window.ui_area, "Select Download Folder")
        if folder_path:
            self.download_path_input.setText(folder_path)

    def _refresh_model_list(self):
        print("Refreshing model list...")
        self.populate_combobox()
        print("Model list refreshed.")

    def populate_combobox(self):
        all_models, all_chains = self.get_models_and_chains(self.session)
        self.combobox.clear()
        for chain in all_chains:
            self.combobox.addItem(chain)

    def get_models_and_chains(self, session):
        all_models, all_chains = [], []
        for model in session.models.list():
            if hasattr(model, "residues"):
                all_models.append(model.id_string)
                for chain_id in set(res.chain_id for res in model.residues):
                    all_chains.append(f"Model {model.id_string}, Chain {chain_id}")
        return all_models, all_chains

    def fill_context_menu(self, menu, x, y):
        action = QAction("Sample Action", menu)
        action.triggered.connect(lambda: print("Sample Action triggered"))
        menu.addAction(action)

    def is_semi_conserved(self, residue1, residue2):
        # Define the criteria for semi-conserved residues based on similarity of amino acid properties
        property_similarities = {
            'G': ['G'],
            'P': ['P'],
            'A': ['A'],
            'S': ['S'],
            'T': ['T'],
            'V': ['V', 'I'],
            'L': ['L', 'I', 'V'],
            'I': ['I', 'L', 'V'],
            'M': ['M'],
            'F': ['F', 'Y', 'W'],
            'Y': ['Y', 'F', 'W'],
            'W': ['W', 'F', 'Y'],
            'H': ['H', 'N', 'Q'],
            'N': ['N', 'H', 'Q'],
            'Q': ['Q', 'H', 'N'],
            'R': ['R', 'K'],
            'K': ['K', 'R'],
            'D': ['D', 'E'],
            'E': ['E', 'D'],
            'C': ['C'],
        }
        return residue1.upper() in property_similarities.get(residue2.upper(), [])

    def add_chain_id(self, file_location, file_name, chain_id, model_id):
        try:
            # Construct the file path
            file_path = file_location + file_name
            # Load the file contents
            with open(file_path) as file:
                lines = file.readlines()
            # Add "/" and chain ID before ":" in the first column starting from line 10
            modified_lines = []
            for i, line in enumerate(lines):
                if i >= 9 and line.startswith("\t:"):
                    modified_lines.append(f"\t#{model_id}/{chain_id}{line[1:]}")
                else:
                    modified_lines.append(line)
            # Construct the new file path
            new_file_path = file_location + "/Chain{}_{}".format(chain_id, file_name[1:])
            # Save the modified contents to the new file
            with open(new_file_path, "w") as file:
                file.writelines(modified_lines)
            print(f"Chain ID added successfully! Modified file saved as: {new_file_path}")
        except FileNotFoundError:
            print("File not found. Please enter a valid file location and name.")
        except Exception as e:
            print("An error occurred:", str(e))
    
    
    def align(self):
	    # --- Create output directory ---
	    mypath = Path().home()
	    output_dir = f"{mypath}/ChopChop_output"
	    if not os.path.exists(output_dir):
	        os.mkdir(output_dir)

	    # --- Get input model and chain IDs from combobox ---
	    model_chain = str(self.combobox.currentText())
	    if len(model_chain) == 16:
	        model = model_chain[6]
	        chain = model_chain[15]
	    elif len(model_chain) == 17:
	        model = model_chain[6] + model_chain[7]
	        chain = model_chain[16]
	    elif len(model_chain) == 18:
	        model = model_chain[6] + model_chain[7] + model_chain[8]
	        chain = model_chain[17]
	    else:
	        print(" Could not parse model/chain from combobox selection.")
	        return

	    # --- NEW: Detect actual starting residue number ---
	    start_resnum = 1
	    model_obj = None
	    for m in self.session.models:
	        if hasattr(m, "chains") and m.id_string == model:
	            model_obj = m
	            break
	    if model_obj:
	        for ch in model_obj.chains:
	            if ch.chain_id == chain:
	                start_resnum = ch.residues[0].number
	                break
	    print(f"Detected chain {chain} in model #{model} starting at residue {start_resnum}")

	    # --- Get FASTA sequence of chosen chain ---
	    run(self.session, f"sequence chain #{model}/{chain}")
	    run(self.session, f"save {output_dir}/chopchop_sequence1.fasta format fasta alignment {model}/{chain}")

	    with open(f"{output_dir}/chopchop_sequence1.fasta", "r", encoding="utf-8-sig") as file:
	        lines = file.readlines()

	    header = True
	    sequence1 = ""
	    for line in lines:
	        if header:
	            header = False
	            continue
	        sequence1 += line.strip()
	    sequence1 = sequence1.replace(" ", "").replace("\n", "")

	    # --- Format FASTA sequence from input text box ---
	    sequence2 = str(self.text1.text()).strip().replace(" ", "").replace("\n", "")

	    # --- Get FASTA from UniProt ID if input contains digits ---
	    if any(char.isdigit() for char in sequence2):
	        url = f"https://www.uniprot.org/uniprot/{sequence2}.fasta"
	        response = requests.get(url)
	        if response.status_code == 200:
	            fasta_file_path = os.path.join(output_dir, f"{sequence2}.fasta")
	            with open(fasta_file_path, "w") as fasta_file:
	                fasta_file.write(response.text)
	            print(f"FASTA file for {sequence2} downloaded successfully.")
	        else:
	            raise ValueError(f"Failed to download FASTA file for UniProt ID: {sequence2}")

	        with open(fasta_file_path, "r", encoding="utf-8-sig") as file:
	            lines = file.readlines()
	        header = True
	        sequence2 = ""
	        for line in lines:
	            if header:
	                header = False
	                continue
	            sequence2 += line.strip()
	        sequence2 = sequence2.replace(" ", "").replace("\n", "")

	    # --- Perform sequence alignment ---
	    protein_name1 = "proteinA"
	    protein_name2 = "proteinB"
	    run(self.session, f"sequence align {sequence1},{sequence2} program muscle")

	    # --- Save last alignment in FASTA format ---
	    alignment_Id = 1
	    while True:
	        try:
	            run(self.session, f"save {output_dir}/chopchop_alignment.fasta format fasta alignment {alignment_Id}")
	            alignment_Id += 1
	        except Exception:
	            break

	    # --- Parse alignment results ---
	    alignment = []
	    with open(f"{output_dir}/chopchop_alignment.fasta") as alignment_file:
	        lines = alignment_file.readlines()

	        start_index = 1
	        end_index = lines.index("\n")
	        sequence1_alignment = "".join(lines[start_index:end_index]).strip().replace("\n", "")

	        start_index = lines.index("\n") + 2
	        sequence2_alignment = "".join(lines[start_index:]).strip().replace("\n", "")

	        print("Sequence 1 Alignment:", sequence1_alignment)
	        print("Sequence 2 Alignment:", sequence2_alignment)

	        for residue1, residue2 in zip(sequence1_alignment, sequence2_alignment):
	            alignment.append([residue1, residue2])

	    # --- Define scoring ---
	    conserved_score = 1
	    semi_conserved_score = 2
	    non_conserved_score = 3
	    gap_score = 4

	    # --- Calculate alignment scores ---
	    scores = []
	    for residue1, residue2 in alignment:
	        if residue1 == residue2:
	            score = conserved_score
	        elif residue1 != '-' and residue2 != '-':
	            if self.is_semi_conserved(residue1, residue2):
	                score = semi_conserved_score
	            else:
	                score = non_conserved_score
	        elif residue1 == '-' or residue2 == '-':
	            score = gap_score
	        else:
	            score = non_conserved_score
	        scores.append(score)

	    # --- Write alignment results to TSV ---
	    output_file = f"{output_dir}/{protein_name1}_{protein_name2}_alignment.csv"
	    with open(output_file, "w", newline="") as file:
	        writer = csv.writer(file, delimiter="\t")
	        writer.writerow([protein_name1, protein_name2, "Alignment Score"])
	        writer.writerow(["", "", "1: Conserved, 2: Semi-Conserved, 3: Not Conserved, 4: Gap"])
	        for (residue1, residue2), score in zip(alignment, scores):
	            writer.writerow([residue1, residue2, score])

	    print("Alignment results saved to:", output_file)

	    # --- Remove double quotes and save final CSV ---
	    with open(output_file, "r") as inp:
	        lines = inp.readlines()
	    output_file_final = f"{output_dir}/{protein_name1}_{protein_name2}_alignment_final.csv"
	    with open(output_file_final, "w") as out:
	        for line in lines:
	            if '"' not in line:
	                out.write(line)

	    # --- Parse alignment into column CSV ---
	    input_file_path = output_file_final
	    with open(input_file_path, "r") as tsv_file:
	        tsv_reader = csv.reader(tsv_file, delimiter="\t")
	        next(tsv_reader)  # header
	        next(tsv_reader)  # scoring row

	        protein_names_column1 = []
	        protein_names_column2 = []
	        alignment_scores_column3 = []
	        for row in tsv_reader:
	            protein_names_column1.append(row[0])
	            protein_names_column2.append(row[1])
	            alignment_scores_column3.append(row[2])

	    output_file_path_column1 = f"{output_dir}/column1_proteins.csv"
	    with open(output_file_path_column1, "w", newline="") as output_file_column1:
	        csv_writer = csv.writer(output_file_column1, delimiter="\t")
	        csv_writer.writerow(["Residue", "Alignment Score"])
	        for i in range(len(protein_names_column1)):
	            csv_writer.writerow([protein_names_column1[i], alignment_scores_column3[i]])
	        output_file_column1.write("# Scoring: 1 = Conserved, 2 = Semi-Conserved, 3 = Not Conserved, 4 = Gap\n")

	    # --- Load the CSV and clean it ---
	    data = []
	    with open(output_file_path_column1, "r") as file:
	        reader = csv.reader(file, delimiter="\t")
	        for row in reader:
	            if not row[0].startswith('-'):
	                data.append(row)

	    # --- Build defattr header ---
	    header = (
	        "#  \n"
	        "#  Use this file to assign the attribute in Chimera with the Define Attribute tool.\n"
	        "#  open .defattr then use the command:\n"
	        "#  color byattribute MUSCLEalignment palette 1,teal:2,palegoldenrod:3,lightsalmon:4,dimgrey\n"
	        "#  \n"
	        "#  Scoring: 1 = Conserved, 2 = Semi-Conserved, 3 = Not Conserved, 4 = Gap\n"
	        "#  \n"
	        "attribute: MUSCLEalignment\n"
	        "match mode: 1-to-1\n"
	        "recipient: residues\n"
	    )

	    # --- Number residues correctly using detected start_resnum ---
	    for i in range(1, len(data) - 1):
	        residue = data[i][0]
	        residue_parts = residue.split(":")
	        residue_parts[0] = residue_parts[0].lstrip("ARNDCEQGHILKMFPSTWYV")
	        residue_number = start_resnum + (i - 1)
	        data[i][0] = f"\t{residue_parts[0]}:{residue_number}"

	    output_filename = f"{os.path.splitext('/column1_proteins.csv')[0]}_MUSCLEalignment4ChimeraX.defattr"
	    output_filepath = f"{output_dir}{output_filename}"

	    with open(output_filepath, "w", newline="") as file:
	        file.write(header)
	        for row in data[1:-1]:
	            file.write("\t".join(row).strip('"') + "\n")

	    print("The updated data has been saved as", output_filename)

	    # --- Color the model in ChimeraX ---
	    self.add_chain_id(output_dir, output_filename, chain, model)
	    run(self.session, f"open {output_dir}/Chain{chain}_{output_filename[1:]}")
	    run(self.session, "color byattribute MUSCLEalignment palette 1,teal:2,palegoldenrod:3,lightsalmon:4,dimgrey")

	    # --- Apply lighting ---
	    run(self.session, "lighting flat")
	    run(self.session, "lighting soft")
    
        

        
    def save(self):
        #copy output folder to download directory
        download_dir = Path(self.download_path_input.text().strip())
        mypath = Path().home()
        if not download_dir.exists():
            try:
                download_dir.mkdir(parents=True, exist_ok=True)
            except Exception as e:
                self.session.logger.warning(f"Failed to create directory {download_dir}: {e}")
                return
        output_dir = "{}/ChopChop_output".format(mypath)
        distutils.dir_util.copy_tree(output_dir, str(download_dir))
        #delete original output folder
        shutil.rmtree(output_dir)
        print("Alignment downloaded in {}".format(download_dir))
        

    def select_models_and_chains(self, session):
        all_chains = []
        all_chains_with_index = {}
        all_models = []
        all_modl_index = {}

        for index, model in enumerate(session.models):
            # Skip models without atomic data (e.g., volumes, maps)
            if not hasattr(model, "chains"):
                continue

            # Skip submodels if needed (optional, but keeps it tidy)
            if '.' in model.id_string:
                continue

            # Store atomic models only
            all_models.append(model.id_string)
            all_modl_index[model.id_string] = index

            # Add each chain in the atomic model
            for chain in model.chains:
                info = f"{model.id_string}:{chain.chain_id}"
                all_chains.append(info)
                all_chains_with_index[info] = index

        return all_models, all_chains, all_modl_index, all_chains_with_index
    
      
    def fill_context_menu(self, menu, x, y):
        # Add any tool-specific items to the given context menu (a QMenu instance).
        # The menu will then be automatically fijlled out with generic tool-related actions
        # (e.g. Hide Tool, Help, Dockable Tool, etc.) 
        #
        # The x,y args are the x() and y() values of QContextMenuEvent, in the rare case
        # where the items put in the menu depends on where in the tool interface the menu
        # was raised.
        #from Qt.QtGui import QAction
        clear_action = QAction("Clear", menu)
        clear_action.triggered.connect(lambda *args: self.line_edit.clear())
        menu.addAction(clear_action)


    
    

