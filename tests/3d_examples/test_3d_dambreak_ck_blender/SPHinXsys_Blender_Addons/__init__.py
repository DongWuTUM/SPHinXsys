bl_info = {
    "name"        : "SPHinXsys",
    "author"      : "Dong Wu & Xiangyu Hu",
    "description" : "SPHinXsys Blender plugin",
    "blender"     : (4, 4, 3),
    "version"     : (0, 0, 2),
    "category"    : "Generic",
}

import bpy
from .props     import SPHinXsysGlobalProps, SPHinXsysProps, registerSPHinXsysProperties, unregisterSPHinXsysProperties
from .saveMesh  import Mesh_OT_ClearSTLFiles, Mesh_OT_SaveSelectdMeshes2File, Mesh_OT_SaveAllMeshes2File
from .panel     import SPHinXsys_PT_Panel
from .SPHinXsys import SPHinXsys_OT_runCase
from .sphinxsys_color_panel import SPH_OT_ApplySPHColorMaterial, SPH_OT_UpdateColorField
from .particle_color_operator import SPH_OT_SetupParticleVisualization

# List of classes to register
classes = [
    SPHinXsysGlobalProps,
    SPHinXsysProps,
    Mesh_OT_ClearSTLFiles,
    Mesh_OT_SaveSelectdMeshes2File,
    Mesh_OT_SaveAllMeshes2File,
    SPHinXsys_PT_Panel,
    SPHinXsys_OT_runCase,
    SPH_OT_ApplySPHColorMaterial,
    SPH_OT_UpdateColorField,
    SPH_OT_SetupParticleVisualization,
]

def register():
    # Unregister in case of reload
    try:
        unregister()
    except Exception:
        pass
    # Register all
    for cls in classes:
        try:
            bpy.utils.register_class(cls)
        except ValueError:
            pass
    # Register properties
    try:
        registerSPHinXsysProperties()
    except Exception:
        pass

def unregister():
    # Unregister properties
    try:
        unregisterSPHinXsysProperties()
    except Exception:
        pass
    # Unregister classes in reverse
    for cls in reversed(classes):
        try:
            bpy.utils.unregister_class(cls)
        except ValueError:
            pass

if __name__ == "__main__":
    register()
