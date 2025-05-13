import bpy, os
from . import props

class Mesh_OT_ClearSTLFiles(bpy.types.Operator):
    bl_idname = "mesh.clearstlfiles"
    bl_label = "Clear STL Files"
    bl_options = {'REGISTER'}

    def execute(self, context):
        # Ensure the .blend file is saved
        blend_path = bpy.data.filepath
        if not blend_path:
            self.report({'ERROR'}, "Please save the .blend file first")
            return {'CANCELLED'}

        # Determine the directory of the .blend file
        blend_dir = os.path.dirname(blend_path)
        mesh_dir  = os.path.join(blend_dir, "mesh")

        # Clear out each of the subfolders
        for sub in ("Fluids", "Rigids", "Elastics"):
            folder = os.path.join(mesh_dir, sub)
            if os.path.isdir(folder):
                for filename in os.listdir(folder):
                    file_path = os.path.join(folder, filename)
                    if os.path.isfile(file_path):
                        os.remove(file_path)

        # Reset the lists of block names
        props.clearAllBlocksName()
        return {'FINISHED'}

class Mesh_OT_SaveSelectdMeshes2File(bpy.types.Operator):
    bl_idname = "mesh.saveselectedmeshes2file"
    bl_label = "Save selected meshes to file"
    bl_options = {'REGISTER'}

    def execute(self, context):
        # Ensure the .blend file is saved
        blend_path = bpy.data.filepath
        if not blend_path:
            self.report({'ERROR'}, "Please save the .blend file first")
            return {'CANCELLED'}

        # Compute output directories relative to the .blend
        blend_dir = os.path.dirname(blend_path)
        mesh_dir  = os.path.join(blend_dir, "mesh")
        fluids_dir   = os.path.join(mesh_dir, "Fluids")
        rigids_dir   = os.path.join(mesh_dir, "Rigids")
        elastics_dir = os.path.join(mesh_dir, "Elastics")

        # Create directories if they don’t exist
        for d in (fluids_dir, rigids_dir, elastics_dir):
            os.makedirs(d, exist_ok=True)
            
        # Clear and then build up the lists of block names
        props.clearAllBlocksName()
        view_layer = context.view_layer
        obj_active = view_layer.objects.active

        # Iterate over all objects in the scene
        for obj in bpy.data.objects:
            if obj.type != 'MESH':
                continue
            if not obj.data.m_SPHinXsysSettings.m_Export:
                continue

            # Select object and set it active
            obj.select_set(True)
            view_layer.objects.active = obj

            media = obj.data.m_SPHinXsysSettings.m_TypeofMedia
            if media == 'Fluid':
                out_dir = fluids_dir
                props.appendBlocksName('Fluid', obj.name)
            elif media == 'RigidSolid':
                out_dir = rigids_dir
                props.appendBlocksName('RigidSolid', obj.name)
            else:
                out_dir = elastics_dir
                props.appendBlocksName('ElasticSolid', obj.name)

            # Construct filename and export to STL
            filename = bpy.path.clean_name(obj.name) + ".stl"
            filepath = os.path.join(out_dir, filename)
            bpy.ops.export_mesh.stl(filepath=filepath, ascii=True, use_selection=True)

            # Deselect after export
            obj.select_set(False)

        # Restore the original active object
        view_layer.objects.active = obj_active
        return {'FINISHED'}
    
class Mesh_OT_SaveAllMeshes2File(bpy.types.Operator):
    bl_idname = "mesh.saveallmeshes2file"
    bl_label = "Save All Meshes"
    bl_options = {'REGISTER'}

    def execute(self, context):
        # Ensure the .blend file is saved
        blend_path = bpy.data.filepath
        if not blend_path:
            self.report({'ERROR'}, "Please save the .blend file first")
            return {'CANCELLED'}

        # Compute output directories relative to the .blend
        blend_dir = os.path.dirname(blend_path)
        mesh_dir  = os.path.join(blend_dir, "mesh")
        fluids_dir   = os.path.join(mesh_dir, "Fluids")
        rigids_dir   = os.path.join(mesh_dir, "Rigids")
        elastics_dir = os.path.join(mesh_dir, "Elastics")

        # Create directories if they don’t exist
        for d in (fluids_dir, rigids_dir, elastics_dir):
            os.makedirs(d, exist_ok=True)

        # Reset and build lists of block names
        props.clearAllBlocksName()
        view_layer = context.view_layer
        obj_active = view_layer.objects.active

        # Iterate over all objects in the scene
        for obj in bpy.data.objects:
            if obj.type != 'MESH':
                continue

            # Select object and set it active
            obj.select_set(True)
            view_layer.objects.active = obj

            media = obj.data.m_SPHinXsysSettings.m_TypeofMedia
            if media == 'Fluid':
                out_dir = fluids_dir
                props.appendBlocksName('Fluid', obj.name)
            elif media == 'RigidSolid':
                out_dir = rigids_dir
                props.appendBlocksName('RigidSolid', obj.name)
            else:
                out_dir = elastics_dir
                props.appendBlocksName('ElasticSolid', obj.name)

            # Construct filename and export all meshes
            filename = bpy.path.clean_name(obj.name) + ".stl"
            filepath = os.path.join(out_dir, filename)
            bpy.ops.export_mesh.stl(filepath=filepath, ascii=True, use_selection=True)

            # Deselect after export
            obj.select_set(False)

        # Restore the original active object
        view_layer.objects.active = obj_active
        return {'FINISHED'}
