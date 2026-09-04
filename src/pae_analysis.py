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
The tool has four tabs:
1. Contacts - load PAE data and generate distance/PAE-filtered contacts
   between two chains
2. Scores - pDockQ, LIS/cLIS/iLIS, buried area and H-bond count for the
   selected chain pair, with citable confidence badges where established
3. Plots - pLDDT, ipTM/pTM, and PAE plots, each drawn fresh in its own
   floating, savable window
4. Residues - select/style the residues behind the generated contacts
Standard output messages are printed to the terminal.
"""

import csv
import json
import webbrowser
from pathlib import Path

from matplotlib.colors import LinearSegmentedColormap

from chimerax.core.tools import ToolInstance
from chimerax.core.commands import run
from chimerax.core.errors import UserError
from chimerax.atomic import Residue, selected_residues
from chimerax.ui import MainToolWindow
from Qt.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QWidget, QComboBox,
    QTabWidget, QFileDialog, QCheckBox, QDoubleSpinBox, QGroupBox,
)
from Qt.QtGui import QFont
from Qt.QtCore import Qt
from .utils import (
    make_scrollable, make_guide_button, show_error, open_figure_dialog, open_residue_table_dialog,
    track_dialog,
    STATUS_GOOD, STATUS_WARNING, STATUS_CRITICAL, PDOCKQ_THRESHOLDS, PDOCKQ_LOW,
    ILIS_THRESHOLD, ILIS_HIGH, threshold_badge, PAE_CMAP as _PAE_CMAP, pae_value_colors,
    get_settings, CATEGORICAL_PALETTE,
)
from . import scoring

_PLDDT_BANDS = [(0, 50, "#FF0000"), (50, 70, "#FFA500"), (70, 90, "#FFFF00"), (90, 100, "#6495ED")]

# Per-bond coloring within a found set of H-bonds: shorter/stronger -> good, longer/weaker -> critical.
_HBOND_CMAP = LinearSegmentedColormap.from_list("hbond", [STATUS_GOOD[0], STATUS_CRITICAL[0]])

# ipTM/pTM bar coloring - same >0.8 / 0.6-0.8 / <0.6 traffic-light convention used
# throughout the AlphaFold-Multimer literature and in the ChopChopMF AlphaFold guide's
# own Confidence Metrics page (lukasinscience.github.io/AlphaFold-Guide).
_IPTM_THRESHOLDS = [(0.8, STATUS_GOOD[0]), (0.6, STATUS_WARNING[0])]
_IPTM_LOW = STATUS_CRITICAL[0]

# Colors selected residues per chain, replacing ChimeraX's default "bychain"
# palette of muted named colors - see utils.CATEGORICAL_PALETTE.
_CHAIN_PALETTE = CATEGORICAL_PALETTE

# Uniform color a residue is reset to once it's no longer part of the current
# highlight, so results from a previous Select action never linger on screen
# (dataviz skill's "muted ink" tone - distinct from every chain/status color above).
_NEUTRAL_RESIDUE_COLOR = "#898781"


def _iptm_bar_color(value):
    for min_value, color in _IPTM_THRESHOLDS:
        if value >= min_value:
            return color
    return _IPTM_LOW


class PAEAnalysis(ToolInstance):

    SESSION_ENDURING = False
    SESSION_SAVE = True
    help = "https://www.cgl.ucsf.edu/chimerax/docs/user/tools/pae_analysis.html"

    def __init__(self, session, tool_name):
        super().__init__(session, tool_name)
        self.json_loaded = False
        self.confidence_scores = {}
        self._open_dialogs = []
        self._highlighted_residues = set()
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
        tabs.addTab(make_scrollable(self._create_contacts_tab()), "1. Contacts")
        tabs.addTab(make_scrollable(self._create_scores_tab()), "2. Scores")
        tabs.addTab(make_scrollable(self._create_confidence_tab()), "3. Plots")
        tabs.addTab(make_scrollable(self._create_residue_tab()), "4. Residues")
        layout.addWidget(tabs)

        container = QWidget()
        container.setLayout(layout)
        self.tool_window.ui_area.setLayout(layout)

        self._refresh_chain_list()

    def _create_contacts_tab(self):
        widget = QWidget()
        layout = QVBoxLayout()
        layout.setAlignment(Qt.AlignTop)
        layout.addWidget(make_guide_button("4-analyze-structure"))

        layout.addWidget(QLabel("PAE (Predicted Aligned Error):"))

        refresh_button = QPushButton("↻ Refresh model list")
        refresh_button.clicked.connect(self._refresh_chain_list)
        layout.addWidget(refresh_button)

        layout.addWidget(QLabel("Select first chain:"))
        self.chain1_selector = QComboBox()
        layout.addWidget(self.chain1_selector)

        layout.addWidget(QLabel("Select second chain:"))
        self.chain2_selector = QComboBox()
        layout.addWidget(self.chain2_selector)

        distance_row = QHBoxLayout()
        distance_row.addWidget(QLabel("Contact distance (Å):"))
        self.distance_spinbox = QDoubleSpinBox()
        self.distance_spinbox.setRange(1.0, 15.0)
        self.distance_spinbox.setSingleStep(0.5)
        self.distance_spinbox.setValue(5.0)
        distance_row.addWidget(self.distance_spinbox)
        layout.addLayout(distance_row)

        pae_row = QHBoxLayout()
        self.max_pae_checkbox = QCheckBox("Limit to residue pairs with PAE ≤")
        self.max_pae_checkbox.setChecked(True)
        self.max_pae_spinbox = QDoubleSpinBox()
        self.max_pae_spinbox.setRange(0.0, 30.0)
        self.max_pae_spinbox.setSingleStep(0.5)
        self.max_pae_spinbox.setValue(12.0)
        self.max_pae_checkbox.toggled.connect(self.max_pae_spinbox.setEnabled)
        pae_row.addWidget(self.max_pae_checkbox)
        pae_row.addWidget(self.max_pae_spinbox)
        pae_row.addWidget(QLabel("Å"))
        layout.addLayout(pae_row)
        pae_hint = QLabel("PAE = how confident AlphaFold is about this pair's relative position; 12 Å is a common cutoff for a confident contact.")
        pae_hint.setWordWrap(True)
        pae_hint.setStyleSheet("color: gray; font-size: 11px;")
        layout.addWidget(pae_hint)

        load_json_button = QPushButton("Load .json file")
        load_json_button.clicked.connect(self._load_json_file)
        layout.addWidget(load_json_button)

        run_button = QPushButton("ChopChop PAE")
        run_button.clicked.connect(self._run_pae_analysis)
        layout.addWidget(run_button)

        widget.setLayout(layout)
        return widget

    def _create_confidence_tab(self):
        widget = QWidget()
        layout = QVBoxLayout()
        layout.setAlignment(Qt.AlignTop)
        layout.addWidget(make_guide_button("4-analyze-structure"))

        self._add_plot_section(
            layout, "pLDDT per Residue", self._draw_plddt_figure,
            "Per-residue confidence from the structure itself (works without a loaded PAE file)."
        )
        self._add_plot_section(
            layout, "ipTM / pTM", self._draw_iptm_figure,
            "Global confidence scores, if present in the loaded .json."
        )
        self._add_plot_section(
            layout, "PAE Matrix", self._draw_pae_heatmap_figure,
            "The full predicted aligned error heatmap."
        )

        widget.setLayout(layout)
        return widget

    def _create_scores_tab(self):
        widget = QWidget()
        layout = QVBoxLayout()
        layout.setAlignment(Qt.AlignTop)
        layout.addWidget(make_guide_button("4-analyze-structure"))

        refresh_button = QPushButton("↻ Refresh Scores")
        refresh_button.clicked.connect(self._refresh_scores)
        layout.addWidget(refresh_button)

        self.scores_container = QWidget()
        self.scores_container_layout = QVBoxLayout()
        self.scores_container_layout.setAlignment(Qt.AlignTop)
        self.scores_container.setLayout(self.scores_container_layout)
        layout.addWidget(self.scores_container)

        button_row = QHBoxLayout()
        deselect_button = QPushButton("Deselect")
        deselect_button.setToolTip("Clear the current selection and reset highlighted residues to the neutral color.")
        deselect_button.clicked.connect(self._deselect_and_reset_style)
        button_row.addWidget(deselect_button)
        table_button = QPushButton("Open residue table")
        table_button.setToolTip("Show the currently selected residues in a sortable, exportable table.")
        table_button.clicked.connect(self._show_residue_table)
        button_row.addWidget(table_button)
        export_button = QPushButton("Export scores as CSV…")
        export_button.clicked.connect(self._export_scores_csv)
        button_row.addWidget(export_button)
        layout.addLayout(button_row)

        widget.setLayout(layout)
        self._refresh_scores()
        return widget

    def _make_score_row(self, name, value_text, badge=None, button_text=None, button_callback=None):
        row = QWidget()
        row_layout = QHBoxLayout()
        row_layout.setContentsMargins(4, 4, 4, 4)
        row_layout.setSpacing(10)
        name_label = QLabel(name)
        bold_font = QFont()
        bold_font.setBold(True)
        name_label.setFont(bold_font)
        name_label.setMinimumWidth(90)
        row_layout.addWidget(name_label)
        row_layout.addWidget(QLabel(value_text))
        if badge is not None:
            bg_color, text_color, text = badge
            badge_label = QLabel(text)
            badge_label.setStyleSheet(
                f"background-color: {bg_color}; color: {text_color}; "
                "border-radius: 8px; padding: 3px 10px; font-weight: bold;"
            )
            row_layout.addWidget(badge_label)
        row_layout.addStretch(1)
        if button_text is not None:
            button = QPushButton(button_text)
            button.clicked.connect(button_callback)
            row_layout.addWidget(button)
        row.setLayout(row_layout)
        return row

    def _clear_layout(self, layout):
        while layout.count():
            item = layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()

    def _current_interface(self):
        """Return (model, chain1, chain2, residues1, residues2), or None if not ready."""
        model = self._get_single_model(warn=False)
        if model is None:
            return None
        chain1 = self.chain1_selector.currentText()
        chain2 = self.chain2_selector.currentText()
        if not chain1 or not chain2 or chain1 == chain2:
            return None
        residues1 = [r for r in model.residues if r.chain_id == chain1]
        residues2 = [r for r in model.residues if r.chain_id == chain2]
        if not residues1 or not residues2:
            return None
        return model, chain1, chain2, residues1, residues2

    def _compute_scores(self):
        """Compute the verified interface scores for the current model/chain pair.

        Returns a dict of results, or None if there is no single model with two
        distinct chains selected. Missing individual values (e.g. LIS/cLIS/iLIS
        without a loaded PAE file) are omitted from the dict rather than guessed.
        """
        interface = self._current_interface()
        if interface is None:
            return None
        model, chain1, chain2, residues1, residues2 = interface

        # Rounded here, once, at the source - so the Scores tab, the residue
        # attributes it writes below (-> Investigate's Chart/Residue panel),
        # and the CSV export all show the same clean numbers instead of raw
        # floating-point precision nobody needs.
        results = {"model": model.id_string, "chain1": chain1, "chain2": chain2}
        pdockq, n_contacts = scoring.compute_pdockq(residues1, residues2)
        results["pdockq"] = round(pdockq, 3)
        results["n_contacts"] = n_contacts
        results["buried_area"] = round(scoring.compute_buried_area(residues1, residues2), 1)
        results["hbonds"] = scoring.compute_hbond_count(self.session, model, chain1, chain2)
        # The residues these pair-level scores actually get written onto below -
        # the real interface (CB-CB <= 8 A contacts), not every residue of both
        # entire chains. See _write_score_attributes for why this matters.
        contact_residues, _ = scoring.compute_contact_residues(residues1, residues2)

        contact_pae_map = {}
        pae_obj = getattr(model, "alphafold_pae", None)
        if pae_obj is not None:
            lis, clis, ilis = scoring.compute_lis_clis_ilis(residues1, residues2, pae_obj)
            results["lis"] = round(lis, 3)
            results["clis"] = round(clis, 3)
            results["ilis"] = round(ilis, 3)
            results["ipsae_d0chn"] = round(scoring.compute_ipsae_d0chn(residues1, residues2, pae_obj), 3)
            contact_pae_map = scoring.compute_contact_pae(residues1, residues2, pae_obj)
            if contact_pae_map:
                results["contact_pae"] = round(sum(contact_pae_map.values()) / len(contact_pae_map), 2)

        self._write_score_attributes(model, results, contact_residues, contact_pae_map)
        return results

    def _set_scalar_attribute(self, attr_name, candidates, value_map):
        """Register (if needed) and set `attr_name` for every residue in
        `value_map` to its given value; every other residue in `candidates`
        has any leftover value from a previous chain pair removed (a scalar
        has no natural "unset" value the way a bool has False, so a stale
        score from a chain that's no longer part of the current pair is
        actually deleted rather than left behind) - mirrors
        _set_bool_attribute's reset behavior for the two boolean attributes."""
        Residue.register_attr(self.session, attr_name, "PAE Analysis", attr_type=float)
        for r in candidates:
            if r in value_map:
                setattr(r, attr_name, value_map[r])
            else:
                try:
                    delattr(r, attr_name)
                except AttributeError:
                    pass

    def _write_score_attributes(self, model, results, contact_residues, contact_pae_map):
        """Write the just-computed interface scores as real ChimeraX residue
        attributes - runs every time scores are (re)computed (Scores tab
        refresh, a loaded .json, CSV export), not just when the user clicks a
        Select button. This is what makes a score show up in Investigate's
        Chart/Plots without an extra manual step, and (via Investigate's own
        snapshotting) survive a close/reopen of the structure durably in its
        .chopchop.json file.

        Written only onto `contact_residues` (the actual CB-CB <= 8 A
        interface, same set as the "chopchop_pae_contact" flag) - NOT every
        residue of both selected chains. These are pair-level aggregate
        scores (one shared number for the whole interface), so writing them
        onto e.g. a floppy, uninvolved N-terminus 300 residues away from the
        actual contact would silently claim that residue is part of the
        scored interface when it isn't - the same value would appear on
        every single residue of both chains, which is exactly what made this
        look broken/meaningless when inspected per-residue in Investigate."""
        pair = contact_residues
        self._set_scalar_attribute("chopchop_pae_pdockq", model.residues, {r: results["pdockq"] for r in pair})
        self._set_scalar_attribute("chopchop_pae_buried_area", model.residues, {r: results["buried_area"] for r in pair})
        self._set_scalar_attribute("chopchop_pae_hbonds", model.residues, {r: float(results["hbonds"]) for r in pair})
        if "lis" in results:
            self._set_scalar_attribute("chopchop_pae_lis", model.residues, {r: results["lis"] for r in pair})
            self._set_scalar_attribute("chopchop_pae_clis", model.residues, {r: results["clis"] for r in pair})
            self._set_scalar_attribute("chopchop_pae_ilis", model.residues, {r: results["ilis"] for r in pair})
            self._set_scalar_attribute("chopchop_pae_ipsae", model.residues, {r: results["ipsae_d0chn"] for r in pair})
            # Unlike the four calls just above (one shared value for the
            # whole pair), this is a genuine per-residue value - each contact
            # residue gets its own mean real-PAE reading to the other side,
            # not the pair's aggregate - so a residue without a qualifying
            # contact correctly ends up with no value at all. Called
            # unconditionally here (even with an empty `contact_pae_map`) so
            # switching to a chain pair with zero contacts still clears any
            # stale per-residue value left over from a previous pair.
            self._set_scalar_attribute(
                "chopchop_pae_contact_pae", model.residues,
                {r: round(v, 2) for r, v in contact_pae_map.items()}
            )

    def _refresh_scores(self):
        self._clear_layout(self.scores_container_layout)
        results = self._compute_scores()
        if results is None:
            self.scores_container_layout.addWidget(
                QLabel("Select two different chains in Tab 1 to see interface scores.")
            )
            return
        has_pae = "lis" in results

        pdockq_badge = threshold_badge(results["pdockq"], PDOCKQ_THRESHOLDS, PDOCKQ_LOW)
        self.scores_container_layout.addWidget(self._make_score_row(
            "pDockQ", f"{results['pdockq']:.3f}  ({results['n_contacts']} contacts ≤ 8 Å)", pdockq_badge,
            "Select contacts", self._select_score_residues
        ))

        self.scores_container_layout.addWidget(self._make_score_row(
            "Buried area", f"{results['buried_area']:.0f} Å²"
        ))
        self.scores_container_layout.addWidget(self._make_score_row(
            "H-bonds", str(results["hbonds"]), None, "Select H-bonds", self._select_hbonds
        ))

        if has_pae:
            ilis_badge = ILIS_HIGH if results["ilis"] >= ILIS_THRESHOLD else None
            self.scores_container_layout.addWidget(self._make_score_row("LIS", f"{results['lis']:.3f}"))
            self.scores_container_layout.addWidget(self._make_score_row("cLIS", f"{results['clis']:.3f}"))
            self.scores_container_layout.addWidget(self._make_score_row(
                "iLIS", f"{results['ilis']:.3f}", ilis_badge,
                "Select confident pairs", self._select_clis_residues
            ))
            self.scores_container_layout.addWidget(self._make_score_row(
                "ipSAE (d0chn)", f"{results['ipsae_d0chn']:.3f}"
            ))
            if "contact_pae" in results:
                bg, fg = pae_value_colors(results["contact_pae"])
                self.scores_container_layout.addWidget(self._make_score_row(
                    "Contact PAE", f"{results['contact_pae']:.2f} Å (avg. over contact residues)",
                    (bg, fg, f"{results['contact_pae']:.2f} Å")
                ))
                pae_note = QLabel(
                    "The actual PAE value (not a derived score) averaged per contact residue, "
                    "same blue-to-white scale as the pseudobonds/PAE Matrix - roughly ≤12 Å is "
                    "confident (Tab 1's own default cutoff), values near/above 20 Å are uncertain."
                )
                pae_note.setWordWrap(True)
                pae_note.setStyleSheet("color: gray; font-size: 11px;")
                self.scores_container_layout.addWidget(pae_note)
        else:
            self.scores_container_layout.addWidget(
                QLabel("LIS / cLIS / iLIS: load a PAE .json file in Tab 1 first.")
            )

    def _delete_pae_contacts_model(self):
        """Interface Scores works on chain coordinates directly, not on pseudobonds -
        remove any leftover 'PAE Contacts' model from Tab 1 so it doesn't linger here."""
        for m in list(self.session.models.list()):
            if m.name == "PAE Contacts":
                m.delete()

    def _set_bool_attribute(self, attr_name, candidates, true_residues):
        """Register (if needed) and set a boolean residue attribute - True for
        residues in `true_residues`, False for the rest of `candidates` (so a
        stale True from a previous chain pair/run doesn't linger). Makes the
        result show up automatically in the Investigate tool's Chart, which
        discovers any residue.custom_attrs generically - no Investigate-side
        code needed per tool."""
        Residue.register_attr(self.session, attr_name, "PAE Analysis", attr_type=bool)
        true_set = set(true_residues)
        for r in candidates:
            setattr(r, attr_name, r in true_set)

    def _select_score_residues(self):
        interface = self._current_interface()
        if interface is None:
            show_error(self.tool_window.ui_area, "ChopChopMF", "Select two different chains in Tab 1 first.")
            return
        model, _, _, residues1, residues2 = interface
        residues, _ = scoring.compute_contact_residues(residues1, residues2)
        self._set_bool_attribute("chopchop_pae_contact", model.residues, residues)
        if not residues:
            message = "No contact residues found within 8 Å for this chain pair."
            self.session.logger.error(message)
            show_error(self.tool_window.ui_area, "ChopChopMF", message)
            return
        self._delete_pae_contacts_model()
        self._select_and_style_residues(residues)

    def _select_hbonds(self):
        interface = self._current_interface()
        if interface is None:
            show_error(self.tool_window.ui_area, "ChopChopMF", "Select two different chains in Tab 1 first.")
            return
        _, chain1, chain2, _, _ = interface
        for m in list(self.session.models.list()):
            if m.name == "ChopChop H-Bonds":
                m.delete()
        command = f'hbonds /{chain1} restrict /{chain2} reveal true select true name "ChopChop H-Bonds"'
        self.session.logger.info(f"Running command: {command}")
        try:
            run(self.session, command)
        except UserError as e:
            self.session.logger.error(str(e))
            show_error(self.tool_window.ui_area, "ChopChopMF", str(e))
            return
        hb_model = next((m for m in self.session.models.list() if m.name == "ChopChop H-Bonds"), None)
        if hb_model is None or len(hb_model.pseudobonds) == 0:
            message = "No hydrogen bonds found between these chains."
            self.session.logger.error(message)
            show_error(self.tool_window.ui_area, "ChopChopMF", message)
            return
        lengths = [pb.length for pb in hb_model.pseudobonds]
        min_len, max_len = min(lengths), max(lengths)
        span = max_len - min_len
        for pb in hb_model.pseudobonds:
            t = (pb.length - min_len) / span if span > 0 else 0.0
            r, g, b, _a = _HBOND_CMAP(t)
            pb.color = (int(r * 255), int(g * 255), int(b * 255), 255)
        self.session.logger.info(
            f"{len(hb_model.pseudobonds)} H-bond(s) selected, colored by length "
            "(green = shorter/stronger, red = longer/weaker)."
        )

    def _select_clis_residues(self):
        interface = self._current_interface()
        if interface is None:
            show_error(self.tool_window.ui_area, "ChopChopMF", "Select two different chains in Tab 1 first.")
            return
        model, _, _, residues1, residues2 = interface
        pae_obj = getattr(model, "alphafold_pae", None)
        if pae_obj is None:
            message = "Please load the corresponding .json file in Tab 1 first."
            self.session.logger.error(message)
            show_error(self.tool_window.ui_area, "ChopChopMF", message)
            return
        residues = scoring.compute_clis_confident_residues(residues1, residues2, pae_obj)
        self._set_bool_attribute("chopchop_pae_clis_confident", model.residues, residues)
        if not residues:
            message = "No cLIS-confident residue pairs found for this chain pair."
            self.session.logger.error(message)
            show_error(self.tool_window.ui_area, "ChopChopMF", message)
            return
        self._delete_pae_contacts_model()
        self._select_and_style_residues(residues)

    def _deselect_and_reset_style(self):
        self._reset_highlighted_residues(self._highlighted_residues)
        self._highlighted_residues = set()
        run(self.session, "select clear")

    def _show_residue_table(self):
        residues = selected_residues(self.session)
        if len(residues) == 0:
            show_error(
                self.tool_window.ui_area, "ChopChopMF",
                "Nothing is selected. Click one of the Select buttons above first."
            )
            return
        dialog = open_residue_table_dialog(self.tool_window.ui_area, residues, "ChopChopMF - Selected Residues")
        track_dialog(self._open_dialogs, dialog)

    def _export_scores_csv(self):
        results = self._compute_scores()
        if results is None:
            show_error(
                self.tool_window.ui_area, "ChopChopMF",
                "Select two different chains in Tab 1 first."
            )
            return
        export_dir = get_settings(self.session).export_dir
        file_path, _ = QFileDialog.getSaveFileName(
            self.tool_window.ui_area, "Save Interface Scores As",
            str(Path(export_dir) / "interface_scores.csv"), "CSV Files (*.csv)"
        )
        if not file_path:
            return
        fieldnames = [
            "model", "chain1", "chain2", "pdockq", "n_contacts", "buried_area",
            "hbonds", "lis", "clis", "ilis", "ipsae_d0chn", "contact_pae",
        ]
        try:
            with open(file_path, "w", newline="", encoding="utf-8-sig") as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerow(results)
        except OSError as e:
            show_error(self.tool_window.ui_area, "ChopChopMF", f"Failed to save scores:\n{e}")

    def _add_plot_section(self, layout, heading, draw_func, description):
        group = QGroupBox(heading)
        group_layout = QVBoxLayout()
        desc_label = QLabel(description)
        desc_label.setWordWrap(True)
        group_layout.addWidget(desc_label)
        open_button = QPushButton("Open Figure")
        open_button.clicked.connect(lambda: self._open_figure(draw_func, heading))
        group_layout.addWidget(open_button)
        group.setLayout(group_layout)
        layout.addWidget(group)

    def _open_figure(self, draw_func, title):
        dialog = open_figure_dialog(self.tool_window.ui_area, draw_func, title)
        track_dialog(self._open_dialogs, dialog)

    def _no_data_message(self, figure, message):
        ax = figure.add_subplot(111)
        ax.text(0.5, 0.5, message, ha="center", va="center", wrap=True)
        ax.axis("off")

    def _pae_row_chain_id(self, residue_or_atom):
        if isinstance(residue_or_atom, Residue):
            return residue_or_atom.chain_id
        return residue_or_atom.residue.chain_id

    def _pae_chain_boundaries(self, pae_obj):
        rows = pae_obj.row_residues_or_atoms()
        chain_ids = [self._pae_row_chain_id(ra) for ra in rows]
        boundaries = [i for i in range(1, len(chain_ids)) if chain_ids[i] != chain_ids[i - 1]]
        return chain_ids, boundaries

    def _draw_plddt_figure(self, figure):
        model = self._get_single_model(warn=False)
        if model is None:
            self._no_data_message(figure, "No single open model.")
            return
        residues = [r for r in model.residues if r.principal_atom is not None]
        if not residues:
            self._no_data_message(figure, "No residues with a representative atom found.")
            return
        numbers = [r.number for r in residues]
        plddt = [r.principal_atom.bfactor for r in residues]

        ax = figure.add_subplot(111)
        for low, high, color in _PLDDT_BANDS:
            ax.axhspan(low, high, color=color, alpha=0.12, linewidth=0)
        chain_ids = [r.chain_id for r in residues]
        for i in range(1, len(chain_ids)):
            if chain_ids[i] != chain_ids[i - 1]:
                ax.axvline(numbers[i], color="black", linewidth=0.5, linestyle="--")
        ax.plot(numbers, plddt, color="black", linewidth=1)
        ax.set_xlabel("Residue number")
        ax.set_ylabel("pLDDT")
        ax.set_ylim(0, 100)
        ax.set_title("Per-residue pLDDT (from B-factor column)")

    def _draw_iptm_figure(self, figure):
        ax = figure.add_subplot(111)
        iptm = self.confidence_scores.get("iptm")
        ptm = self.confidence_scores.get("ptm")
        if iptm is None and ptm is None:
            self._no_data_message(
                figure,
                "No ipTM/pTM data found in the loaded .json file.\n"
                "(Common for AlphaFold3 - check its separate summary_confidences file.)"
            )
            return
        labels, values = [], []
        if iptm is not None:
            labels.append("ipTM")
            values.append(iptm)
        if ptm is not None:
            labels.append("pTM")
            values.append(ptm)
        ax.barh(labels, values, color=[_iptm_bar_color(v) for v in values])
        ax.axvline(0.6, color=STATUS_WARNING[0], linewidth=0.8, linestyle="--")
        ax.axvline(0.8, color=STATUS_GOOD[0], linewidth=0.8, linestyle="--")
        ax.set_xlim(0, 1)
        ax.set_xlabel("Score")
        ax.set_title("Confidence scores (green >0.8, yellow 0.6-0.8, red <0.6)")
        for i, v in enumerate(values):
            ax.text(v + 0.02, i, f"{v:.2f}", va="center")

    def _draw_pae_heatmap_figure(self, figure):
        model = self._get_single_model(warn=False)
        pae_obj = getattr(model, "alphafold_pae", None) if model is not None else None
        if pae_obj is None:
            self._no_data_message(figure, "Load a PAE .json file in Tab 1 first.")
            return
        matrix = pae_obj.pae_matrix
        ax = figure.add_subplot(111)
        im = ax.imshow(matrix, cmap=_PAE_CMAP, vmin=0, vmax=30, origin="upper")
        figure.colorbar(im, ax=ax, label="PAE (Å)")
        _, boundaries = self._pae_chain_boundaries(pae_obj)
        for b in boundaries:
            ax.axhline(b - 0.5, color="black", linewidth=0.5)
            ax.axvline(b - 0.5, color="black", linewidth=0.5)
        ax.set_xlabel("Scored residue")
        ax.set_ylabel("Aligned residue")
        ax.set_title("Predicted Aligned Error")

    def _create_residue_tab(self):
        widget = QWidget()
        layout = QVBoxLayout()
        layout.setAlignment(Qt.AlignTop)
        layout.addWidget(make_guide_button("4-analyze-structure"))

        self.remove_pb_model_checkbox = QCheckBox("Delete pseudobonds after selecting residues")
        self.remove_pb_model_checkbox.setToolTip(
            "The blue/red 'PAE Contacts' pseudobond lines from Tab 1 are normally deleted once their "
            "residues are selected below, to keep the scene tidy. Uncheck this to keep the pseudobonds "
            "visible for further inspection."
        )
        self.remove_pb_model_checkbox.setChecked(True)
        layout.addWidget(self.remove_pb_model_checkbox)

        run_button = QPushButton("ChopChop PAE interaction Residues")
        run_button.clicked.connect(self._run_pae_selection)
        layout.addWidget(run_button)

        pae_viewer_button = QPushButton("Open PAE Viewer")
        pae_viewer_button.setToolTip("PAE Viewer by Christoph Elfmann and Jörg Stülke (subtiwiki.uni-goettingen.de)")
        pae_viewer_button.clicked.connect(lambda: webbrowser.open("https://subtiwiki.uni-goettingen.de/v4/paeViewerDemo"))
        layout.addWidget(pae_viewer_button)

        widget.setLayout(layout)
        return widget

    def _get_single_model(self, warn=True):
        models = [m for m in self.session.models.list() if hasattr(m, "residues")]
        if len(models) != 1:
            message = (
                "PAE Analysis requires exactly one open atomic model with residue information. "
                "If your complex was opened as separate chains, combine them into one model first "
                "(e.g. the ChimeraX 'combine' command), or open the multimer prediction as a single "
                "structure."
            )
            self.session.logger.error(message)
            if warn:
                show_error(self.tool_window.ui_area, "ChopChopMF", message)
            return None
        return models[0]

    def _refresh_chain_list(self):
        self.chain1_selector.clear()
        self.chain2_selector.clear()
        model = self._get_single_model(warn=False)
        if model is None:
            return
        chains = sorted(set(res.chain_id for res in model.residues))
        for chain in chains:
            self.chain1_selector.addItem(chain)
            self.chain2_selector.addItem(chain)
        if "A" in chains:
            self.chain1_selector.setCurrentText("A")
        if "B" in chains:
            self.chain2_selector.setCurrentText("B")
        self.session.logger.info("Chain lists refreshed.")

    def _load_json_file(self):
        model = self._get_single_model()
        if model is None:
            return
        start_dir = str(Path.home())
        model_path = getattr(model, "filename", None)
        if model_path:
            model_dir = Path(model_path).parent
            if model_dir.is_dir():
                start_dir = str(model_dir)
        file_path, _ = QFileDialog.getOpenFileName(
            self.tool_window.ui_area, "Select AlphaFold PAE .json file",
            start_dir, "JSON Files (*.json)"
        )
        if not file_path:
            return
        # plot false: ChimeraX's own "alphafold pae" command opens its own
        # separate "AlphaFold Predicted Aligned Error" plot window by default
        # (plot=True unless colorDomains is set) - redundant clutter here since
        # ChopChopMF's own Plots tab already draws the PAE heatmap on demand.
        command = f'alphafold pae #{model.id_string} file "{file_path}" plot false'
        self.session.logger.info(f"Running command: {command}")
        try:
            run(self.session, command)
        except UserError as e:
            self.json_loaded = False
            self.session.logger.error(str(e))
            show_error(self.tool_window.ui_area, "ChopChopMF", f"Failed to load PAE data:\n{e}")
            return
        if getattr(model, "alphafold_pae", None) is None:
            self.json_loaded = False
            message = (
                "The PAE data could not be linked to the open model. Make sure the .json file "
                "matches the currently open structure."
            )
            self.session.logger.error(message)
            show_error(self.tool_window.ui_area, "ChopChopMF", message)
            return
        self.json_loaded = True
        self.confidence_scores = self._read_confidence_scores(file_path)
        self.session.logger.info(f"PAE data loaded for model #{model.id_string}.")
        self._refresh_scores()

    def _read_confidence_scores(self, file_path):
        """Best-effort extraction of ipTM/pTM from the loaded PAE .json, if present.

        ChimeraX itself has no concept of ipTM/pTM (only PAE), but ColabFold's scores.json -
        the very file used for PAE above - commonly also carries these as top-level keys.
        AlphaFold3 keeps them in a separate summary_confidences file, so nothing is found
        there; this is deliberately tolerant and never raises.
        """
        try:
            with open(file_path) as f:
                data = json.load(f)
        except (OSError, ValueError):
            return {}
        if isinstance(data, list) and data and isinstance(data[0], dict):
            data = data[0]
        if not isinstance(data, dict):
            return {}
        return {key: data[key] for key in ("iptm", "ptm") if key in data}

    def _run_pae_analysis(self):
        model = self._get_single_model()
        if model is None:
            return
        if not self.json_loaded or getattr(model, "alphafold_pae", None) is None:
            message = "Please load the corresponding .json file before running the analysis."
            self.session.logger.error(message)
            show_error(self.tool_window.ui_area, "ChopChopMF", message)
            return
        chain1 = self.chain1_selector.currentText()
        chain2 = self.chain2_selector.currentText()
        if not chain1 or not chain2 or chain1 == chain2:
            message = "Both chains must be selected and different."
            self.session.logger.error(message)
            show_error(self.tool_window.ui_area, "ChopChopMF", message)
            return

        distance = self.distance_spinbox.value()
        command = f"alphafold contacts /{chain1} to /{chain2} distance {distance}"
        if self.max_pae_checkbox.isChecked():
            command += f" maxPae {self.max_pae_spinbox.value()}"
        self.session.logger.info(f"Running command: {command}")
        try:
            run(self.session, command)
        except UserError as e:
            self.session.logger.error(str(e))
            show_error(self.tool_window.ui_area, "ChopChopMF", str(e))

    def _reset_highlighted_residues(self, residues):
        """Revert previously highlighted residues to the neutral color/style, so a
        new Select action never leaves the old one's residues colored on screen."""
        if not residues:
            return
        run(self.session, "select clear")
        for res in residues:
            for atom in res.atoms:
                atom.selected = True
        run(self.session, f"color sel {_NEUTRAL_RESIDUE_COLOR}")
        run(self.session, "hide sel atoms")
        run(self.session, "select clear")

    def _select_and_style_residues(self, residues):
        residues = set(residues)
        for res in residues:
            self.session.logger.info(f"Residue: chain {res.chain_id}, number {res.number}, name {res.name}")
        self._reset_highlighted_residues(self._highlighted_residues - residues)
        run(self.session, "select clear")
        for res in residues:
            for atom in res.atoms:
                atom.selected = True
        run(self.session, "show sel")
        run(self.session, "style sel stick")
        chain_ids = sorted({res.chain_id for res in residues})
        for i, chain_id in enumerate(chain_ids):
            color = _CHAIN_PALETTE[i % len(_CHAIN_PALETTE)]
            run(self.session, f'color /{chain_id} & sel {color}')
        run(self.session, "color sel byhetero")
        self._highlighted_residues = residues

    def _run_pae_selection(self):
        pb_model = next((m for m in self.session.models.list() if m.name == "PAE Contacts"), None)
        if pb_model is None:
            message = "No pseudobond model named 'PAE Contacts' found. Run the PAE Contacts analysis in the first tab first."
            self.session.logger.error(message)
            show_error(self.tool_window.ui_area, "ChopChopMF", message)
            return
        residues = set()
        for pb in pb_model.pseudobonds:
            for atom in pb.atoms:
                residues.add(atom.residue)
        if not residues:
            message = "The 'PAE Contacts' model has no pseudobonds to select residues from."
            self.session.logger.error(message)
            show_error(self.tool_window.ui_area, "ChopChopMF", message)
            return
        self._select_and_style_residues(residues)
        if self.remove_pb_model_checkbox.isChecked():
            pb_model.delete()
            self.session.logger.info("Pseudobond model deleted. Residues remain selected and styled by heteroatom.")
        else:
            self.session.logger.info("Residues remain selected and styled by heteroatom. 'PAE Contacts' model kept.")
