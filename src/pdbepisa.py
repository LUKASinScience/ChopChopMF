

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
PDBePISA ChimeraX Tool (ChopChopMF)

This tool preserves the original two-tab interface:
1) Interface Scoring (parses PDBePISA XML to tag residues by interface class and writes a defattr)
2) ΔG Filter (reads SOLVATIONENERGY from PDBePISA XML and colors by per-residue ΔG)

Added functionality (without removing any existing function):
- Load multiple ΔG XML files.
- Dropdown to select the "active" XML for ΔG coloring.
- "Append mode" checkbox to apply ΔG coloring cumulatively across all loaded files.
- Color swatches next to ΔG palette fields (visual cue for the entered color).
- Neutral band & baseline color controls (GUI-only). When enabled, paints a baseline (default lightgrey),
  and keeps ΔG==0, BSA==0, and |ΔG| ≤ ε residues neutral (i.e., they stay grey).
- "Plot ΔG Values" button opens a window showing a bar chart, scatter plot, and a text list
  for all residues that currently have a delta_g_score attribute. Residue labels include ChainID:ResNum.
  The plot uses the same color scheme as the structure (current GUI palette).

Launches from the current folder by default and logs actions to the ChimeraX logger.
"""

from chimerax.core.tools import ToolInstance
from chimerax.core.commands import run
from chimerax.ui import MainToolWindow
from Qt.QtWidgets import (
    QVBoxLayout, QLabel, QPushButton, QWidget, QFileDialog,
    QComboBox, QFrame, QHBoxLayout, QLineEdit, QTabWidget,
    QSlider, QCheckBox, QTextEdit, QSizePolicy, QSpacerItem
)
from Qt.QtCore import Qt
from Qt.QtGui import QFont
import xml.etree.ElementTree as ET
import webbrowser
import os

# ---- plotting ----
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


class DGPlotWindow(QWidget):
    """Popup window that plots all residues with delta_g_score (bar + scatter) and lists exact values."""
    def __init__(self, session, triplets, color_func):
        super().__init__()
        self.setWindowTitle("ΔG Values")
        self.session = session
        # triplets = [(chain_id, resnum_str, dg_value_float), ...]
        self.triplets = triplets
        self.color_func = color_func  # callable: float -> color string

        layout = QVBoxLayout(self)

        # Build colors array based on current palette
        values = [dg for (_, _, dg) in triplets]
        colors = [self.color_func(v) for v in values]

        # Matplotlib figure with two axes
        fig = Figure(figsize=(9, 7))
        canvas = FigureCanvas(fig)
        layout.addWidget(canvas)

        ax_bar = fig.add_subplot(2, 1, 1)
        ax_scat = fig.add_subplot(2, 1, 2)

        labels = [f"{c}:{r}" for (c, r, _) in triplets]

        # Bar chart (top) — use palette colors
        ax_bar.bar(range(len(values)), values, color=colors)
        ax_bar.set_title("ΔG per Residue (Bar)")
        ax_bar.set_ylabel("ΔG (kcal/mol)")
        # Show fewer x tick labels to avoid clutter; list is in the text area anyway
        if len(labels) <= 50:
            ax_bar.set_xticks(range(len(labels)))
            ax_bar.set_xticklabels(labels, rotation=90, fontsize=7)
        else:
            ax_bar.set_xticks([])

        # Scatter plot (bottom) — use palette colors
        ax_scat.scatter(range(len(values)), values, s=14, c=colors)
        ax_scat.set_title("ΔG Distribution (Scatter)")
        ax_scat.set_xlabel("Residue Index")
        ax_scat.set_ylabel("ΔG (kcal/mol)")

        fig.tight_layout()
        canvas.draw()

        # Text list of values
        txt = QTextEdit()
        txt.setReadOnly(True)
        lines = [f"{c}:{r}\t{dg:.3f}" for (c, r, dg) in triplets]
        txt.setText("\n".join(lines))
        layout.addWidget(txt)

        self.resize(900, 800)


class PDBePISA(ToolInstance):

    SESSION_ENDURING = False
    SESSION_SAVE = True
    help = "https://www.ebi.ac.uk/pdbe/pisa/pistart.html"

    def __init__(self, session, tool_name):
        super().__init__(session, tool_name)
        self.display_name = "PDBePISA Tool"
        self.tool_window = MainToolWindow(self)

        # ΔG state for multi-file handling
        self.dg_xml_data = {}        # {basename.xml: parsed XML root}
        self.active_dg_file = None   # name of the active XML in the dropdown
        self.append_mode = False     # cumulative vs. single-file coloring

        # GUI-driven neutral controls (user can change in-tab)
        self._neutral_on = True
        self._epsilon = 0.01
        self._neutral_color = "lightgrey"

        self._build_tabs()
        self.tool_window.manage('side')

    # ---- window / tabs ----
    def _build_tabs(self):
        self.tabs = QTabWidget()

        tab1 = self._build_ui_tab1()
        tab2 = self._build_ui_tab2()
        self.tabs.addTab(tab1, "Interface Scoring")
        self.tabs.addTab(tab2, "ΔG Filter")

        self.tabs.setCurrentIndex(0)

        layout = QVBoxLayout()
        layout.addWidget(self.tabs)
        container = QWidget()
        container.setLayout(layout)
        self.tool_window.ui_area.setLayout(layout)

    # ---- TAB 1: Interface Scoring (unchanged) ----
    def _build_ui_tab1(self):
        layout = QVBoxLayout()

        pisa_button = QPushButton("Open PDBePISA Website")
        pisa_button.clicked.connect(lambda: webbrowser.open(self.help))
        layout.addWidget(pisa_button)

        title = QLabel("PDBePISA: Macromolecular Interface Analysis")
        title_font = QFont()
        title_font.setBold(True)
        title_font.setPointSize(14)
        title.setFont(title_font)
        layout.addWidget(title)

        description = QLabel(
            "How to get the XML file\n"
            "Open PDBePISA and enter your PDB ID or upload your .pdb file, analyse the structure for interfaces. "
            "At the PISA Interface List, press the 'Details' button under Interfaces, "
            "scroll down to Interfacing residues at the PISA Interface page and press the 'XML' button with Display level Residues. "
            "Save the .xml file and open it in ChopChopMF."
        )
        description.setWordWrap(True)
        layout.addWidget(description)

        file_button = QPushButton("Select PDBePISA XML File")
        file_button.clicked.connect(self.select_file)
        layout.addWidget(file_button)

        layout.addWidget(QLabel("Select defattr file for coloring interfaces:"))
        chop_button = QPushButton("ChopChop PISA Interfaces")
        chop_button.clicked.connect(self.chopchop_interfaces)
        layout.addWidget(chop_button)

        layout.addWidget(QLabel("Scoring Scheme:"))
        color_bar = QVBoxLayout()
        colors = [
            ("Buried", "darkorange"),
            ("Hydrogen Bond", "cornflowerblue"),
            ("Salt Bridge", "purple"),
        ]
        for label, color in colors:
            row = QHBoxLayout()
            swatch = QFrame()
            swatch.setFixedSize(20, 20)
            swatch.setStyleSheet(f"background-color: {color}; border: 1px solid black;")
            row.addWidget(swatch)
            row.addWidget(QLabel(label))
            color_bar.addLayout(row)
        layout.addLayout(color_bar)

        layout.addWidget(QLabel("Customize Colors for Scoring Scheme:"))
        self.color_inputs = {}
        color_entries = [
            ("1 (Buried):", "darkorange"),
            ("2 (Hydrogen Bond):", "cornflowerblue"),
            ("3 (Salt Bridge):", "purple"),
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

        tab1 = QWidget()
        tab1.setLayout(layout)
        return tab1

    def select_file(self):
        xml_file, _ = QFileDialog.getOpenFileName(
            caption="Select PDBePISA XML File",
            filter="XML Files (*.xml)"
        )
        if xml_file:
            residues_by_chain, defattr_lines = self.parse_pdbepisa_xml(xml_file)
            commands = self.generate_commands(residues_by_chain)
            self.run_chimerax_commands(commands)
            self.write_defattr_file(defattr_lines, xml_file)

    def chopchop_interfaces(self):
        defattr_file, _ = QFileDialog.getOpenFileName(
            caption="Select Generated Defattr File",
            filter="Defattr Files (*.defattr)"
        )
        if defattr_file:
            run(self.session, f"open {defattr_file}")
            run(self.session, "color byattribute residue_score palette 1,darkorange:2,cornflowerblue:3,purple")
            self.session.logger.info("ChopChop Interfaces coloring applied.")

    def apply_custom_colors(self):
        palette = []
        for key, input_field in self.color_inputs.items():
            color = input_field.text().strip()
            palette.append(f"{key},{color}")
        color_command = "color byattribute residue_score palette " + ":".join(palette)
        try:
            run(self.session, color_command)
            self.session.logger.info("Custom color scheme applied.")
        except Exception as e:
            self.session.logger.error(f"Error applying custom color scheme: {e}")

    def parse_pdbepisa_xml(self, xml_file):
        try:
            tree = ET.parse(xml_file)
            root = tree.getroot()
        except ET.ParseError as e:
            self.session.logger.error(f"XML parsing error: {e}")
            return {}, []

        residues_by_chain = {}
        defattr_lines = [
            "#\n",
            "# Use this file to assign the attribute in ChimeraX with the Define Attribute tool.\n",
            "# Use the command: color byattribute residue_score palette 1,darkorange:2,cornflowerblue:3,purple\n",
            "#\n",
            "# Scoring: 1 = Buried, 2 = Buried with H-bond, 3 = Buried with Salt Bridge\n",
            "#\n",
            "attribute: residue_score\n",
            "match mode: 1-to-1\n",
            "recipient: residues\n"
        ]

        for residue in root.findall('.//RESIDUE'):
            structure_element = residue.find('STRUCTURE')
            bsa_element = residue.find('BURIEDSURFACEAREA')
            hsdc_element = residue.find('HSDC')

            if structure_element is not None and bsa_element is not None and bsa_element.text:
                structure_text = structure_element.text.strip()
                bsa_value = float(bsa_element.text.strip())
                hsdc_value = hsdc_element.text.strip() if hsdc_element is not None else ""

                chain_id, _, residue_info = structure_text.partition(':')
                chain_id = chain_id.strip()
                res_num = ''.join(filter(str.isdigit, residue_info))

                if bsa_value > 0:
                    if chain_id not in residues_by_chain:
                        residues_by_chain[chain_id] = []
                    residues_by_chain[chain_id].append(res_num)

                    score = 1
                    if hsdc_value == "H":
                        score = 2
                    elif hsdc_value == "HS":
                        score = 3

                    defattr_lines.append(f"/{chain_id}:{res_num}\t{score}\n")

        return residues_by_chain, defattr_lines

    def generate_commands(self, residues_by_chain):
        if not residues_by_chain:
            self.session.logger.info("No interacting residues found or XML file is empty/incorrect.")
            return None

        chimerax_command = []
        for chain, residues in residues_by_chain.items():
            residues_str = ','.join(residues)
            chimerax_command.append(f"#1/{chain}:{residues_str}")
        return 'select ' + ' '.join(chimerax_command)

    def run_chimerax_commands(self, commands):
        if commands:
            run(self.session, commands)
            self.session.logger.info(f"Executed ChimeraX command: {commands}")
            run(self.session, "color sel darkorange")
            self.session.logger.info("Colored selected residues in darkorange.")

    def write_defattr_file(self, defattr_lines, xml_file):
        output_file = os.path.splitext(xml_file)[0] + "_output.defattr"
        try:
            with open(output_file, 'w') as file:
                for line in defattr_lines:
                    if line.startswith("/"):
                        file.write(f"\t{line}")
                    else:
                        file.write(line)
            self.session.logger.info(f"Defattr file written to {output_file}")
        except Exception as e:
            self.session.logger.error(f"Error writing defattr file: {e}")

    # ---- TAB 2: ΔG Filter (extended with multi-XML support + neutral + plotting) ----
    def _build_ui_tab2(self):
        layout = QVBoxLayout()

        file_btn = QPushButton("Load PDBePISA XML File")
        file_btn.clicked.connect(self.load_xml_for_dg)
        layout.addWidget(file_btn)

        # Active XML selector
        layout.addWidget(QLabel("Active XML file:"))
        self.dg_selector = QComboBox()
        self.dg_selector.currentIndexChanged.connect(self.set_active_dg_file)
        layout.addWidget(self.dg_selector)

        # Append mode checkbox
        self.append_checkbox = QCheckBox("Append mode (cumulative across all loaded XMLs)")
        self.append_checkbox.stateChanged.connect(self.toggle_append_mode)
        layout.addWidget(self.append_checkbox)

        # Neutral band controls (GUI-owned state)
        row_neutral = QHBoxLayout()
        self.neutral_checkbox = QCheckBox("Neutral band (±ε) & baseline color")
        self.neutral_checkbox.setChecked(True)
        self.neutral_checkbox.stateChanged.connect(self._sync_neutral_state)
        row_neutral.addWidget(self.neutral_checkbox)

        row_neutral.addWidget(QLabel("ε (kcal/mol):"))
        self.eps_edit = QLineEdit("0.01")
        self.eps_edit.setFixedWidth(60)
        self.eps_edit.textChanged.connect(lambda txt, w=self.eps_edit: self._validate_float_lineedit(w, txt, 0.0, 0.5))
        self.eps_edit.textChanged.connect(lambda _: self._sync_neutral_state())
        row_neutral.addWidget(self.eps_edit)

        row_neutral.addWidget(QLabel("Neutral color:"))
        self.neutral_color = QLineEdit("lightgrey")
        self.neutral_color.setFixedWidth(120)
        self.neutral_color.setToolTip("lightgrey")
        self.neutral_color.textChanged.connect(lambda txt, w=self.neutral_color: self._validate_color_lineedit(w, txt))
        self.neutral_swatch = QFrame(); self.neutral_swatch.setFixedSize(20, 20)
        self.neutral_swatch.setStyleSheet("background-color: lightgrey; border: 1px solid black;")
        self.neutral_color.textChanged.connect(lambda c: self.neutral_swatch.setStyleSheet(
            f"background-color: {c.strip()}; border: 1px solid black;"))
        row_neutral.addWidget(self.neutral_color)
        row_neutral.addWidget(self.neutral_swatch)
        row_neutral.addItem(QSpacerItem(10, 10, QSizePolicy.Expanding, QSizePolicy.Minimum))
        layout.addLayout(row_neutral)

        # Cutoff controls (unchanged)
        self.slider_label = QLabel("ΔG cutoff: 0.50 kcal/mol")
        layout.addWidget(self.slider_label)

        self.dg_slider = QSlider(Qt.Horizontal)
        self.dg_slider.setMinimum(-100)
        self.dg_slider.setMaximum(250)
        self.dg_slider.setValue(50)
        self.dg_slider.valueChanged.connect(self.update_slider_label)
        layout.addWidget(self.dg_slider)

        self.only_above_cutoff = QCheckBox("Only show residues ≥ cutoff")
        self.only_above_cutoff.setChecked(True)
        layout.addWidget(self.only_above_cutoff)

        info_label = QLabel(
            "ΔG values are read from the 'SOLVATIONENERGY' field in the PDBePISA XML.\n"
            "Residues with ΔG = 0.0 kcal/mol and residues with BURIEDSURFACEAREA = 0 are excluded."
        )
        info_label.setWordWrap(True)
        layout.addWidget(info_label)

        self.apply_btn = QPushButton("ChopChop ΔG Coloring")
        self.apply_btn.clicked.connect(self.apply_dg_coloring)
        layout.addWidget(self.apply_btn)

        # Plot ΔG Values button directly under ChopChop button
        plot_btn = QPushButton("Plot ΔG Values")
        plot_btn.clicked.connect(self.plot_dg_values)
        layout.addWidget(plot_btn)

        layout.addWidget(QLabel("ΔG Palette (Ranges excluding ~0 kcal/mol):"))
        self.dg_palette_inputs = {}

        color_ranges = [
            ("-1.0 to -0.51", "navy"),
            ("-0.5 to -0.01", "deepskyblue"),
            ("0.01 to 0.5", "moccasin"),
            ("0.51 to 0.99", "goldenrod"),
            ("1.0 to 1.49", "orange"),
            ("1.5 to 1.99", "darkorange"),
            ("2.0 to 2.49", "red"),
            ("2.5 and higher", "firebrick"),
        ]

        for label, default_color in color_ranges:
            row = QHBoxLayout()
            row.addWidget(QLabel(f"{label} kcal/mol:"))
            entry = QLineEdit(default_color)
            self.dg_palette_inputs[label] = entry
            row.addWidget(entry)

            # live color swatch
            swatch = QFrame()
            swatch.setFixedSize(20, 20)
            swatch.setStyleSheet(f"background-color: {default_color}; border: 1px solid black;")
            row.addWidget(swatch)

            def _sync_swatch(line_edit=entry, frame=swatch):
                c = line_edit.text().strip()
                frame.setStyleSheet(f"background-color: {c}; border: 1px solid black;")
            entry.textChanged.connect(_sync_swatch)

            layout.addLayout(row)

        apply_color_btn = QPushButton("Apply New Color Scheme")
        apply_color_btn.clicked.connect(self.apply_custom_dg_colors)
        layout.addWidget(apply_color_btn)

        dg_tab = QWidget()
        dg_tab.setLayout(layout)
        return dg_tab

    # ---- helpers for UI validation / state ----
    def _validate_color_lineedit(self, le: QLineEdit, text: str):
        ok = bool(text.strip())
        le.setStyleSheet("" if ok else "border: 1px solid red;")

    def _validate_float_lineedit(self, le: QLineEdit, text: str, lo: float, hi: float):
        try:
            v = float(text.strip())
            ok = (lo <= v <= hi)
        except Exception:
            ok = False
        le.setStyleSheet("" if ok else "border: 1px solid red;")

    def _sync_neutral_state(self):
        self._neutral_on = self.neutral_checkbox.isChecked()
        try:
            v = float(self.eps_edit.text().strip())
        except Exception:
            v = 0.01
        # clamp
        if v < 0.0: v = 0.0
        if v > 0.5: v = 0.5
        self._epsilon = v
        c = self.neutral_color.text().strip()
        self._neutral_color = c if c else "lightgrey"

    # ---- palette helper for plotting ----
    def _color_for_value(self, v: float) -> str:
        """
        Map a ΔG value to the current GUI palette color (same bins as structure coloring).
        Bins:
          [-inf,-1.0) use first bin color (navy by default)
          [-1.0,-0.51], [-0.5,-0.01], [0.01,0.5], [0.51,0.99], [1.0,1.49], [1.5,1.99], [2.0,2.49], [2.5, +inf)
        Values in the neutral band (|ΔG| ≤ ε) are treated as neutral baseline color.
        """
        # Neutral check
        try:
            eps_val = float(self.eps_edit.text().strip())
        except Exception:
            eps_val = self._epsilon
        eps_val = max(0.0, min(0.5, eps_val))
        if self.neutral_checkbox.isChecked() and abs(v) <= eps_val:
            return self.neutral_color.text().strip() or "lightgrey"

        # Read colors from GUI fields
        # Order must match label order added in _build_ui_tab2
        labels = [
            "-1.0 to -0.51",
            "-0.5 to -0.01",
            "0.01 to 0.5",
            "0.51 to 0.99",
            "1.0 to 1.49",
            "1.5 to 1.99",
            "2.0 to 2.49",
            "2.5 and higher",
        ]
        colors = [(self.dg_palette_inputs[l].text().strip() or "black") for l in labels]

        # Map value to bin
        if v < -1.0:
            return colors[0]
        if -1.0 <= v <= -0.51:
            return colors[0]
        if -0.5 <= v <= -0.01:
            return colors[1]
        if 0.01 <= v <= 0.5:
            return colors[2]
        if 0.51 <= v <= 0.99:
            return colors[3]
        if 1.0 <= v <= 1.49:
            return colors[4]
        if 1.5 <= v <= 1.99:
            return colors[5]
        if 2.0 <= v <= 2.49:
            return colors[6]
        # v >= 2.5 or in ( -0.01, 0.01 ) falls through
        if v >= 2.5:
            return colors[7]
        # Anything extremely close to zero but outside neutral band gets nearest side:
        return colors[2] if v > 0 else colors[1]

    # ---- plotting callback ----
    def plot_dg_values(self):
        """Collect all residues with delta_g_score and open a plotting window (uses current GUI palette)."""
        triplets = []  # (chain_id, resnum_str, dg_val)
        try:
            from chimerax.atomic import AtomicStructure
            models = [m for m in self.session.models.list() if isinstance(m, AtomicStructure)]
            if not models:
                self.session.logger.warning("No atomic models open; cannot plot ΔG.")
                return
            for m in models:
                for r in m.residues:
                    if hasattr(r, 'delta_g_score'):
                        try:
                            val = float(getattr(r, 'delta_g_score'))
                        except Exception:
                            continue
                        # Respect neutral baseline choice (omit pure zeros from plot if you prefer; here we include them)
                        chain_id = r.chain_id if r.chain_id is not None else ""
                        resnum = str(r.number)
                        if getattr(r, 'insertion_code', None) and r.insertion_code.strip():
                            resnum += r.insertion_code.strip()
                        triplets.append((chain_id, resnum, val))
        except Exception as e:
            self.session.logger.error(f"Failed to collect ΔG values for plotting: {e}")
            return

        if not triplets:
            self.session.logger.info("No residues with delta_g_score found. Apply ΔG coloring first.")
            return

        try:
            win = DGPlotWindow(self.session, triplets, self._color_for_value)
            if not hasattr(self, "_plot_windows"):
                self._plot_windows = []
            self._plot_windows.append(win)  # keep ref
            win.show()
            self.session.logger.info(f"Plotted {len(triplets)} residues with ΔG values (palette-matched).")
        except Exception as e:
            self.session.logger.error(f"Failed to create ΔG plot window: {e}")

    # ---- multi-XML handling (ΔG tab) ----
    def load_xml_for_dg(self):
        fname, _ = QFileDialog.getOpenFileName(
            caption="Select PDBePISA XML File",
            filter="XML Files (*.xml)"
        )
        if not fname:
            return
        try:
            tree = ET.parse(fname)
            root = tree.getroot()
        except ET.ParseError as e:
            self.session.logger.error(f"XML parsing error: {e}")
            return

        base = os.path.basename(fname)
        self.dg_xml_data[base] = root

        # Update dropdown if new
        existing = [self.dg_selector.itemText(i) for i in range(self.dg_selector.count())]
        if base not in existing:
            self.dg_selector.addItem(base)

        if self.active_dg_file is None:
            self.active_dg_file = base
            self.dg_selector.setCurrentText(base)

        self.session.logger.info(f"Loaded PDBePISA XML for ΔG: {base}")

    def set_active_dg_file(self):
        current = self.dg_selector.currentText()
        if current in self.dg_xml_data:
            self.active_dg_file = current
            self.session.logger.info(f"Active ΔG XML set to: {current}")

    def toggle_append_mode(self, state):
        self.append_mode = (state == Qt.Checked)
        self.session.logger.info("ΔG coloring mode: Append" if self.append_mode else "ΔG coloring mode: Reset")

    # ---- ΔG coloring (kept original flow; now with neutral baseline & band) ----
    def apply_dg_coloring(self):
        if not self.dg_xml_data:
            self.session.logger.error("No XML file loaded for ΔG filtering.")
            return

        cutoff = self.dg_slider.value() / 100
        filter_active = self.only_above_cutoff.isChecked()

        # Apply neutral baseline first (overwrites previous schemes)
        if self._neutral_on:
            try:
                run(self.session, f"color all {self._neutral_color}")
                self.session.logger.info(f"Applied neutral baseline: {self._neutral_color}")
            except Exception as e:
                self.session.logger.warning(f"Baseline color failed ('{self._neutral_color}'): {e}")

        attr_lines = [
            "#\n",
            "# ΔG/solvation energy per residue\n",
            "# Source: SOLVATIONENERGY field from PDBePISA XML\n",
            "attribute: delta_g_score\n",
            "match mode: 1-to-1\n",
            "recipient: residues\n"
        ]

        to_select = []

        # Choose which XML roots to process
        if self.append_mode:
            xml_roots = list(self.dg_xml_data.values())
        else:
            if not self.active_dg_file:
                self.session.logger.error("No active XML file selected for ΔG coloring.")
                return
            xml_roots = [self.dg_xml_data[self.active_dg_file]]

        # Build defattr and selection from chosen roots
        seen = set()  # avoid duplicate entries when appending
        skipped_neutral = 0
        for root in xml_roots:
            for residue in root.findall('.//RESIDUE'):
                struct = residue.find('STRUCTURE')
                dg = residue.find('SOLVATIONENERGY')
                bsa = residue.find('BURIEDSURFACEAREA')

                if struct is None or dg is None or dg.text is None:
                    continue

                try:
                    dg_val = float(dg.text.strip())
                    if bsa is None or float(bsa.text.strip()) == 0.0:
                        # no buried area -> keep neutral
                        continue
                except Exception:
                    continue

                # Keep ΔG==0 or |ΔG| ≤ ε neutral (omit from defattr)
                if dg_val == 0.0:
                    continue
                try:
                    eps_val = float(self.eps_edit.text().strip())
                except Exception:
                    eps_val = self._epsilon
                if self.neutral_checkbox.isChecked() and abs(dg_val) <= max(0.0, min(0.5, eps_val)):
                    skipped_neutral += 1
                    continue

                if filter_active and dg_val < cutoff:
                    continue

                chain_id, _, residue_info = struct.text.strip().partition(':')
                res_num = ''.join(filter(str.isdigit, residue_info))
                key = (chain_id, res_num)
                if key in seen:
                    continue
                seen.add(key)

                attr_lines.append(f"\t/{chain_id}:{res_num}\t{dg_val:.3f}\n")
                to_select.append(f"#1/{chain_id}:{res_num}")

        # Write defattr and apply
        output_file = "dg_coloring_output.defattr"
        try:
            with open(output_file, 'w') as f:
                f.writelines(attr_lines)
            self.session.logger.info(
                f"ΔG defattr file written to {output_file} (neutral-kept residues: {skipped_neutral})"
            )
        except Exception as e:
            self.session.logger.error(f"Error writing ΔG defattr file: {e}")
            return

        if to_select:
            run(self.session, 'select ' + ' '.join(to_select))
        else:
            run(self.session, 'select clear')

        run(self.session, f"open {output_file}")
        # Use default palette; user can override with Apply New Color Scheme (plot uses GUI palette too)
        run(self.session, "color byattribute delta_g_score palette navy:deepskyblue:white:moccasin:goldenrod:orange:darkorange:red:firebrick range -1.0,2.5")
        self.session.logger.info("ΔG-based coloring applied.")

    def apply_custom_dg_colors(self):
        palette_mapping = {
            "-1.0 to -0.51": "-1.0",
            "-0.5 to -0.01": "-0.5",
            "0.01 to 0.5": "0.01",
            "0.51 to 0.99": "0.51",
            "1.0 to 1.49": "1.0",
            "1.5 to 1.99": "1.5",
            "2.0 to 2.49": "2.0",
            "2.5 and higher": "2.5"
        }

        palette = []
        for label in palette_mapping.keys():
            color = self.dg_palette_inputs[label].text().strip()
            palette.append(color)

        try:
            palette_str = ":".join(palette)
            run(self.session, f"color byattribute delta_g_score palette {palette_str} range -1.0,2.5")
            self.session.logger.info("Custom ΔG color scheme applied.")
        except Exception as e:
            self.session.logger.error(f"Error applying ΔG color scheme: {e}")

    def update_slider_label(self):
        val = self.dg_slider.value() / 100
        self.slider_label.setText(f"ΔG cutoff: {val:.2f} kcal/mol")


# ---- ChimeraX bundle API glue (unchanged) ----
from chimerax.core.toolshed import BundleAPI
class _PDBePISABundleAPI(BundleAPI):
    @staticmethod
    def start_tool(session, tool_name):
        return PDBePISA(session, tool_name)

bundle_api = _PDBePISABundleAPI()
