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

"""
Shared helpers used across the ChopChopMF tools:
- make_scrollable: wrap a tab's content widget so it can be scrolled.
- make_guide_button: button that opens the online usage guide at a specific section.
- busy_cursor: wait cursor + button disable during long-running actions.
- show_error: modal error dialog to accompany logger warnings.
- safe_extractall: zip-slip-safe extraction of downloaded archives.
- get_settings: persistent, session-cached ChopChopMF settings (download_dir).
- open_residue_table_dialog: floating, sortable table of residues with a CSV export button.
- open_table_dialog: same, but for an arbitrary table (any headers/rows), used by Investigate.
"""

from contextlib import contextmanager
from pathlib import Path
import csv
import os
import webbrowser

from Qt.QtWidgets import (
    QScrollArea, QApplication, QMessageBox, QPushButton, QDialog, QVBoxLayout,
    QFileDialog, QTableWidget, QTableWidgetItem, QHeaderView,
)
from Qt.QtCore import Qt, QSize
from matplotlib.colors import LinearSegmentedColormap

from chimerax.core.settings import Settings

MIN_TAB_HEIGHT = 700
GUIDE_BASE_URL = "https://lukasinscience.github.io/ChopChopMF/usage/"


class _TallSizeHintScrollArea(QScrollArea):
    """QScrollArea whose sizeHint() has a taller minimum, without setting an
    actual minimumSize.

    ChimeraX sizes a newly-docked tool panel from `ui_area.sizeHint()` (see
    chimerax.ui.gui.ToolWindow.manage()), but does NOT clamp the panel's
    minimum size to it afterwards - the dock can still be shrunk by the user
    or by ChimeraX's own resize-fallback logic. `setMinimumHeight()` instead
    installs a hard floor that Qt enforces everywhere, which fights the
    scroll area's whole purpose (being smaller than its content) and was
    found to actually break the scrollbar in practice. Overriding sizeHint()
    only influences that one-time initial-size request.
    """
    def __init__(self, min_height):
        super().__init__()
        self._min_height = min_height

    def sizeHint(self):
        hint = super().sizeHint()
        return QSize(hint.width(), max(hint.height(), self._min_height))


def make_scrollable(widget, min_height=MIN_TAB_HEIGHT):
    """Wrap a tab's content widget in a scroll area so it can be scrolled when it doesn't fit.

    Tabs open at a comfortable default height (via sizeHint, not a hard
    minimum size - see `_TallSizeHintScrollArea`); combined with
    `layout.setAlignment(Qt.AlignTop)` on each tab's own layout, shorter tabs
    stay packed at the top instead of being stretched apart, while taller
    tabs keep scrolling normally.
    """
    scroll = _TallSizeHintScrollArea(min_height)
    scroll.setWidgetResizable(True)
    scroll.setWidget(widget)
    return scroll


def make_guide_button(anchor, label="📖 Open Guide / Tutorial"):
    """Create a button that opens the ChopChopMF usage guide at a specific section in the browser."""
    button = QPushButton(label)
    button.setToolTip("Opens the ChopChopMF usage guide for this tool in your browser.")
    button.clicked.connect(lambda: webbrowser.open(f"{GUIDE_BASE_URL}#{anchor}"))
    return button


@contextmanager
def busy_cursor(button=None):
    """Show a wait cursor and disable a trigger button for the duration of a long-running action."""
    if button is not None:
        button.setEnabled(False)
    QApplication.setOverrideCursor(Qt.WaitCursor)
    try:
        yield
    finally:
        QApplication.restoreOverrideCursor()
        if button is not None:
            button.setEnabled(True)


def show_error(parent, title, message):
    """Show a modal error dialog (in addition to logging)."""
    QMessageBox.warning(parent, title, message)


# Fixed status palette (good/warning/serious/critical) - shared by every tool
# that shows a confidence badge (PAE Analysis's Scores tab, Investigate's
# Chart/Residue panel), so a color always means the same verdict wherever it
# appears. Never themed, never reused for anything but a confidence verdict.
# Each entry is (background, text) - text color picked per background for
# readability (dark backgrounds get white text, bright ones get dark text).
STATUS_GOOD = ("#0ca30c", "#FFFFFF")
STATUS_WARNING = ("#fab219", "#0b0b0b")
STATUS_SERIOUS = ("#ec835a", "#0b0b0b")
STATUS_CRITICAL = ("#d03b3b", "#FFFFFF")

# Interpretation bands with a citable source - see docs/acknowledgements.md.
# Metrics without an established "good/bad" threshold in the literature (buried
# area, H-bonds, LIS, cLIS, ipSAE d0chn) intentionally get no badge anywhere.
PDOCKQ_THRESHOLDS = [
    (0.5, *STATUS_GOOD, "High confidence"),
    (0.23, *STATUS_WARNING, "Weak / medium"),
]
PDOCKQ_LOW = (*STATUS_CRITICAL, "Poor")
ILIS_THRESHOLD = 0.223
ILIS_HIGH = (*STATUS_GOOD, "High-confidence interaction")


# Categorical palette (validated CVD-safe order, dataviz skill) - fixed
# assignment order, shared by every tool that colors a fixed set of named
# categories (PAE Analysis's per-chain residue coloring, Cell Biology's
# per-kinase/per-organelle badges). Never reassigned based on a filter/
# selection - a category always gets the same color wherever it appears.
CATEGORICAL_PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#4a3aa7"]

# Second, distinct categorical palette - shares no hue with CATEGORICAL_PALETTE
# above (or with STATUS_*/PAE_CMAP), so a Cell Biology tool's per-kinase/per-
# organelle coloring never gets confused with a PAE Analysis chain or a
# confidence badge elsewhere in the bundle. Fixed assignment order, same rule
# as CATEGORICAL_PALETTE - never reassigned based on a filter/selection.
CELLBIO_PALETTE = [
    "#0f9b8e",  # teal
    "#c0399f",  # magenta
    "#9a8b1f",  # olive
    "#7856d1",  # violet
    "#8a5a2b",  # brown
    "#17a3c4",  # cyan
    "#a33b52",  # maroon
    "#6fae1f",  # lime
    "#4c6a8a",  # slate
]


def threshold_badge(value, thresholds, low):
    """thresholds: [(min_value, bg_color, text_color, label), ...] sorted
    highest-first. Returns the first (bg_color, text_color, label) whose
    min_value `value` meets, or `low` (itself a badge tuple, or None) if it's
    below every threshold."""
    for min_value, bg_color, text_color, label in thresholds:
        if value >= min_value:
            return bg_color, text_color, label
    return low


# Mirrors ChimeraX's own built-in 'pae'/'paecontacts' colormap (chimerax.core.colors)
# - blue = confident/low error, red-ish/gray/white = uncertain/high error, 0-30 Å -
# so a real PAE value looks the same everywhere it's shown: the pseudobonds ChimeraX
# itself colors in PAE Analysis's Contacts tab, the PAE Matrix plot, and (via
# pae_value_colors below) any residue attribute derived from actual PAE values.
PAE_CMAP = LinearSegmentedColormap.from_list(
    "pae", list(zip(
        [0, 5 / 30, 10 / 30, 15 / 30, 20 / 30, 25 / 30, 1.0],
        ["#0000FF", "#6495ED", "#FFFF00", "#FFA500", "#808080", "#D3D3D3", "#FFFFFF"],
    ))
)


def pae_value_colors(value, vmax=30.0):
    """(bg_hex, fg_hex) for a raw PAE value in Angstrom, using the same
    blue-to-white scale as ChimeraX's own PAE coloring (see PAE_CMAP above).
    Text color is picked from the background's perceived luminance."""
    t = max(0.0, min(1.0, value / vmax))
    r, g, b, _a = PAE_CMAP(t)
    bg = f"#{int(r * 255):02x}{int(g * 255):02x}{int(b * 255):02x}"
    luminance = 0.299 * r + 0.587 * g + 0.114 * b
    fg = "#FFFFFF" if luminance < 0.5 else "#0b0b0b"
    return bg, fg


def track_dialog(open_dialogs, dialog):
    """Append `dialog` to a tool's `_open_dialogs` list, and remove it again
    once the dialog closes (its `finished` signal) - without this, every
    floating plot/table dialog opened over a session's lifetime stays
    referenced forever, even long after the user closed it (WA_DeleteOnClose
    only frees the underlying Qt widget, not this Python-side list entry)."""
    open_dialogs.append(dialog)
    dialog.finished.connect(lambda _result: open_dialogs.remove(dialog) if dialog in open_dialogs else None)


def open_figure_dialog(parent, draw_func, title="Plot", on_pick=None):
    """Show a fresh copy of a plot in its own floating, resizable window with a Save button.

    A matplotlib Figure can only be rendered by one canvas at a time, so this creates a new
    Figure/canvas here rather than reusing one already embedded in a tab. `draw_func(figure)`
    must (re)draw the desired content onto the given Figure - the same function used to draw
    the compact, embedded preview is meant to be reused here.

    `on_pick`, if given, is connected to the canvas's matplotlib `pick_event` (the drawn
    artist(s) must set their own `picker` for anything to actually fire) - lets a caller
    make a plot clickable (e.g. jump to the clicked data point's underlying object) without
    every caller needing its own canvas/dialog plumbing. Optional and backward-compatible -
    existing callers that don't pass it behave exactly as before.
    """
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.figure import Figure

    dialog = QDialog(parent)
    dialog.setWindowTitle(title)
    dialog.setAttribute(Qt.WA_DeleteOnClose)
    dialog.resize(800, 600)
    layout = QVBoxLayout()

    figure = Figure(figsize=(8, 6))
    canvas = FigureCanvas(figure)
    draw_func(figure)
    canvas.draw()
    if on_pick is not None:
        canvas.mpl_connect("pick_event", on_pick)
    layout.addWidget(canvas)

    def _save():
        file_path, _ = QFileDialog.getSaveFileName(
            dialog, "Save Plot As", str(Path.home() / "plot.png"),
            "PNG Image (*.png);;PDF Document (*.pdf);;SVG Image (*.svg)"
        )
        if not file_path:
            return
        try:
            figure.savefig(file_path, dpi=200, bbox_inches="tight")
        except OSError as e:
            show_error(dialog, "ChopChopMF", f"Failed to save plot:\n{e}")

    save_button = QPushButton("Save…")
    save_button.clicked.connect(_save)
    layout.addWidget(save_button)

    dialog.setLayout(layout)
    dialog.show()
    return dialog


def open_residue_table_dialog(parent, residues, title="Selected Residues"):
    """Show a sortable table (chain, number, name, pLDDT) of `residues` in its own
    floating window, with a "Export as CSV..." button. `residues` is any iterable
    of chimerax.atomic.Residue; pLDDT is read from each residue's principal atom
    B-factor (the standard AlphaFold convention) and left blank if unavailable.
    """
    rows = []
    for r in sorted(residues, key=lambda r: (r.chain_id, r.number)):
        atom = r.principal_atom
        plddt = f"{atom.bfactor:.1f}" if atom is not None else ""
        rows.append((r.chain_id, str(r.number), r.name, plddt))

    dialog = QDialog(parent)
    dialog.setWindowTitle(title)
    dialog.setAttribute(Qt.WA_DeleteOnClose)
    dialog.resize(420, 500)
    layout = QVBoxLayout()

    table = QTableWidget(len(rows), 4)
    table.setHorizontalHeaderLabels(["Chain", "Residue #", "Name", "pLDDT"])
    table.setSortingEnabled(False)
    for i, (chain_id, number, name, plddt) in enumerate(rows):
        for col, value in enumerate((chain_id, number, name, plddt)):
            item = QTableWidgetItem(value)
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            table.setItem(i, col, item)
    table.setSortingEnabled(True)
    table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
    table.verticalHeader().setVisible(False)
    layout.addWidget(table)

    def _export():
        file_path, _ = QFileDialog.getSaveFileName(
            dialog, "Save Residue Table As", str(Path.home() / "residues.csv"), "CSV Files (*.csv)"
        )
        if not file_path:
            return
        try:
            with open(file_path, "w", newline="", encoding="utf-8-sig") as f:
                writer = csv.writer(f)
                writer.writerow(["chain", "residue_number", "name", "plddt"])
                writer.writerows(rows)
        except OSError as e:
            show_error(dialog, "ChopChopMF", f"Failed to save residue table:\n{e}")

    export_button = QPushButton("Export as CSV…")
    export_button.clicked.connect(_export)
    layout.addWidget(export_button)

    dialog.setLayout(layout)
    dialog.show()
    return dialog


def open_table_dialog(parent, headers, rows, title="Table", csv_filename="table.csv", extra_buttons=None):
    """Show an arbitrary read-only table (any headers/string rows) in a large
    floating window, with an "Export as CSV..." button. `extra_buttons` is an
    optional list of (label, callback) pairs added below it - each callback is
    called with the dialog as its only argument (e.g. to save a related file).
    """
    dialog = QDialog(parent)
    dialog.setWindowTitle(title)
    dialog.setAttribute(Qt.WA_DeleteOnClose)
    dialog.resize(1000, 700)
    layout = QVBoxLayout()

    table = QTableWidget(len(rows), len(headers))
    table.setHorizontalHeaderLabels(headers)
    table.setSortingEnabled(False)
    for i, row in enumerate(rows):
        for col, value in enumerate(row):
            item = QTableWidgetItem(value)
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            table.setItem(i, col, item)
    table.setSortingEnabled(True)
    table.verticalHeader().setVisible(False)
    layout.addWidget(table)

    def _export():
        file_path, _ = QFileDialog.getSaveFileName(
            dialog, "Save Table As", str(Path.home() / csv_filename), "CSV Files (*.csv)"
        )
        if not file_path:
            return
        try:
            with open(file_path, "w", newline="", encoding="utf-8-sig") as f:
                writer = csv.writer(f)
                writer.writerow(headers)
                writer.writerows(rows)
        except OSError as e:
            show_error(dialog, "ChopChopMF", f"Failed to save table:\n{e}")

    export_button = QPushButton("Export as CSV…")
    export_button.clicked.connect(_export)
    layout.addWidget(export_button)

    for label, callback in (extra_buttons or []):
        button = QPushButton(label)
        button.clicked.connect(lambda checked=False, cb=callback: cb(dialog))
        layout.addWidget(button)

    dialog.setLayout(layout)
    dialog.show()
    return dialog


def safe_extractall(zip_ref, dest_dir):
    """Extract a ZipFile, rejecting any member that would land outside dest_dir (zip-slip guard)."""
    dest_dir = Path(dest_dir).resolve()
    for member in zip_ref.namelist():
        target = (dest_dir / member).resolve()
        if target != dest_dir and not str(target).startswith(str(dest_dir) + os.sep):
            raise ValueError(f"Unsafe path in zip file: {member}")
    zip_ref.extractall(dest_dir)


class ChopChopMFSettings(Settings):
    AUTO_SAVE = {
        "download_dir": str(Path.home() / "Downloads"),
        # Suggested starting folder for one-off CSV/Markdown "Save As" dialogs
        # (PAE Analysis/Batch Analysis/Investigate exports) - the dialog still
        # asks every time, this only replaces the previous plain Path.home().
        "export_dir": str(Path.home() / "Downloads"),
        # PDBePISA's .defattr output folder - "" (falsy) means the existing
        # default behavior (next to the loaded PISA XML file), same meaning
        # as the None it replaces. Persistent (not a per-tool-instance
        # attribute) so it can be shown/changed centrally in the Settings
        # tool even while PDBePISA itself isn't open.
        "defattr_output_dir": "",
    }


_settings_cache = {}


def get_settings(session):
    """Return a session-cached, persistent ChopChopMF settings object."""
    if session not in _settings_cache:
        _settings_cache[session] = ChopChopMFSettings(session, "ChopChopMF")
    return _settings_cache[session]
