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
ChopChopMF bundle API wiring for ChimeraX.
Maps toolbar/tool names to their concrete tool classes. 

"""

from chimerax.core.toolshed import BundleAPI
from chimerax.core.commands import run

# Subclass from chimerax.core.toolshed.BundleAPI
class ChopChopMFAPI(BundleAPI):
    api_version = 1  # Register commands with BundleInfo and ToolInfo
    _version_logged = False  # internal flag to avoid duplicate log messages

    @staticmethod
    def start_tool(session, bi, ti):
        # Log version only once per ChimeraX session
        if not ChopChopMFAPI._version_logged:
            session.logger.info("ChopChopMF-1.1")
            ChopChopMFAPI._version_logged = True

        # Start the corresponding tool based on the name in bundle_info.xml.
        tool_name = ti.name  # Get the tool name from the ToolInfo instance

        # Check the tool name and instantiate the appropriate tool class
        if tool_name == "Sequence":
            from .chopalignment import ChopChopMFalignment
            return ChopChopMFalignment(session, tool_name)
        
        elif tool_name == "AlphaMissense":
            from .chopgetmissense import ChopChopGetMissense
            return ChopChopGetMissense(session, tool_name)
        
        elif tool_name == "Missense":
            from .chopmissense import ChopChopMissense
            return ChopChopMissense(session, tool_name)
        
        elif tool_name == "Crop Structure":
            from .cropstructure import CropStructureTool
            return CropStructureTool(session, tool_name)
            
        elif tool_name == "Duplicate Structure":
            from .DuplicateStructureTool import DuplicateStructureTool
            return DuplicateStructureTool(session, tool_name)
        
        elif tool_name == "AlphaFold2":
            from .alphafoldinfo import AlphaFold2
            return AlphaFold2(session, tool_name)
        
        elif tool_name == "PAE Analysis":
            from .pae_analysis import PAEAnalysis
            return PAEAnalysis(session, tool_name)
                
        elif tool_name == "PDBePISA":
            from .pdbepisa import PDBePISA
            return PDBePISA(session, tool_name)
            
        elif tool_name == "Foldseek Analysis":
            from .foldseekanalysis import FoldseekAnalysis
            return FoldseekAnalysis(session, tool_name)
            
        elif tool_name == "Undo":
            session.logger.info("Undoing last action...")
            run(session, "undo")
            return None  
        
        else:
            raise ValueError(f"Unknown tool name: {tool_name}")

    @staticmethod
    def get_class(class_name):
        # Given a class name, return the corresponding tool class.
        if class_name == "ChopChopMFalignment":
            from .chopalignment import ChopChopMFalignment
            return ChopChopMFalignment
        
        elif class_name == "ChopChopGetMissense":
            from .chopgetmissense import ChopChopGetMissense
            return ChopChopGetMissense
        
        elif class_name == "ChopChopMissense":
            from .chopmissense import ChopChopMissense
            return ChopChopMissense
            
        elif class_name == "CropStructureTool":
            from .cropstructure import CropStructureTool
            return CropStructureTool
            
        elif class_name == "DuplicateStructureTool":
            from .DuplicateStructureTool import DuplicateStructureTool
            return DuplicateStructureTool
        
        elif class_name == "AlphaFold2":
            from .alphafoldinfo import AlphaFold2
            return AlphaFold2

        elif class_name == "PAEAnalysis":
            from .pae_analysis import PAEAnalysis
            return PAEAnalysis
            
        elif class_name == "PDBePISA":
            from .pdbepisa import PDBePISA
            return PDBePISA
            
        elif class_name == "FoldseekAnalysis":
            from .foldseekanalysis import FoldseekAnalysis
            return FoldseekAnalysis
        
        else:
            raise ValueError(f"Unknown class name '{class_name}'")

    @staticmethod
    def run_provider(session, name, mgr, **kw):
        # This method is called to invoke the functionality of the provider.
        if mgr == session.toolbar:
            if name == "Sequence":
                from .chopalignment import ChopChopMFalignment
                return ChopChopMFalignment(session, name)
            
            elif name == "AlphaMissense":
                from .chopgetmissense import ChopChopGetMissense
                return ChopChopGetMissense(session, name)
            
            elif name == "Missense":
                from .chopmissense import ChopChopMissense
                return ChopChopMissense(session, name)
                
            elif name == "Crop Structure":
                from .cropstructure import CropStructureTool
                return CropStructureTool(session, name)
            
            elif name == "Duplicate Structure":
                from .DuplicateStructureTool import DuplicateStructureTool
                return DuplicateStructureTool(session, name)
                
            elif name == "AlphaFold2":
                from .alphafoldinfo import AlphaFold2
                return AlphaFold2(session, name)

            elif name == "PAE Analysis":
                from .pae_analysis import PAEAnalysis
                return PAEAnalysis(session, name)
            
            elif name == "PDBePISA":
                from .pdbepisa import PDBePISA
                return PDBePISA(session, name)
                
            elif name == "Foldseek Analysis":
                from .foldseekanalysis import FoldseekAnalysis
                return FoldseekAnalysis(session, name)
            
            elif name == "Undo":
                run(session, "undo")
                return None

            else:
                raise ValueError(f"Unknown toolbar provider: {name}")

        else:
            raise ValueError(f"Unknown manager: {mgr}")

# Create the ``bundle_api`` object that ChimeraX expects.
bundle_api = ChopChopMFAPI()


