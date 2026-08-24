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
"""

from contextlib import contextmanager
from pathlib import Path
import os
import webbrowser

from Qt.QtWidgets import QScrollArea, QApplication, QMessageBox, QPushButton
from Qt.QtCore import Qt, QSize

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
    }


_settings_cache = {}


def get_settings(session):
    """Return a session-cached, persistent ChopChopMF settings object."""
    if session not in _settings_cache:
        _settings_cache[session] = ChopChopMFSettings(session, "ChopChopMF")
    return _settings_cache[session]
