# Changelog

All notable changes to ChopChopMF are documented here. This file was introduced with version 1.3, so entries for 1.2 and earlier are limited to what could be reliably reconstructed — see the [ChimeraX Toolshed page](https://cxtoolshed.rbvi.ucsf.edu/apps/chimeraxchopchopmf) for the historical release listing.

## 1.3

### Fixed
- `chopalignment.py`: removed a duplicate `fill_context_menu` method whose "Clear" action referenced a non-existent widget and would crash on click.
- `chopalignment.py`: replaced `distutils.dir_util.copy_tree` (removed in Python 3.12) with `shutil.copytree`.
- `chopalignment.py`: model/chain selection in the alignment tool no longer relies on fragile string-length parsing; it now reads the value directly from the dropdown.
- `chopmissense.py`: replaced a hardcoded `/tmp` path with `tempfile.gettempdir()`, fixing the tool on Windows.
- `pdbepisa.py`: interface selection and ΔG coloring no longer hardcode model `#1`; both tabs now have a model selector.

### Security / Robustness
- Added request timeouts to all network calls (UniProt, AlphaMissense, AlphaSync).
- Quoted file paths passed into ChimeraX `run()` commands.
- Added input validation for the sequence/UniProt ID field.
- Added zip-slip protection when extracting downloaded AlphaMissense archives (`utils.safe_extractall`).

### Added
- Tab content in all 9 tools is now scrollable and opens at a comfortable default height; each tab's layout is top-aligned so extra space collects as blank space at the bottom instead of stretching widgets apart.
- Every tab of every tool now has a "📖 Open Guide / Tutorial" button at the top, linking directly to the relevant section of the [online usage guide](https://lukasinscience.github.io/ChopChopMF/usage/) in the user's browser.
- Busy cursor and button-disable feedback during network-bound actions (sequence alignment, AlphaMissense fetch, AlphaSync fetch).
- Error dialogs (in addition to log messages) for the main failure points of each tool's primary action.
- Download folder path is now remembered between sessions (`chopalignment.py`, `chopgetmissense.py`).
- New shared `src/utils.py` module (`make_scrollable`, `make_guide_button`, `busy_cursor`, `show_error`, `safe_extractall`, `get_settings`) to remove duplication introduced by the scrolling/feedback work.

### Changed
- Unified the version string across `bundle_info.xml`, `setup.py`, `src/__init__.py`, and the docs (previously inconsistently at 1.1/1.2).
- Removed dead per-file `BundleAPI` classes in `pae_analysis.py`, `foldseekanalysis.py`, `alphafoldinfo.py`, `pdbepisa.py` (never loaded by ChimeraX; only `src/__init__.py`'s `bundle_api` is used).
- Declared `requests` and `matplotlib` as explicit bundle dependencies.

## 1.2 and earlier

Pre-dates this changelog. See the Toolshed release notes for details.
