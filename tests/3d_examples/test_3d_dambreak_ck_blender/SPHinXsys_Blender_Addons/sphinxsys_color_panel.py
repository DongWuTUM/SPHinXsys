import bpy
from .SPHinXsys import water_obj_name, frame_update_handler

# Operator to apply vertex color material
class SPH_OT_ApplySPHColorMaterial(bpy.types.Operator):
    bl_idname = "sphinxsys.apply_color_material"
    bl_label = "Show Color"
    bl_description = "Apply material that reads vertex color layer 'SPH_color' for particle coloring"

    @classmethod
    def poll(cls, context):
        return bpy.data.objects.get(water_obj_name) is not None

    def execute(self, context):
        obj = bpy.data.objects.get(water_obj_name)
        if not obj or obj.type != 'MESH':
            self.report({'ERROR'}, f"Mesh object '{water_obj_name}' not found")
            return {'CANCELLED'}

        # Create or get material
        mat_name = f"{water_obj_name}_ColorMat"
        mat = bpy.data.materials.get(mat_name) or bpy.data.materials.new(mat_name)
        mat.use_nodes = True
        nodes = mat.node_tree.nodes
        links = mat.node_tree.links
        # Clear existing nodes
        for node in list(nodes):
            nodes.remove(node)
        # Output node
        out_node = nodes.new(type="ShaderNodeOutputMaterial")
        out_node.location = (300, 0)
        # Principled BSDF
        bsdf_node = nodes.new(type="ShaderNodeBsdfPrincipled")
        bsdf_node.location = (0, 0)
        # Attribute node
        vc_node = nodes.new(type="ShaderNodeAttribute")
        vc_node.location = (-300, 0)
        vc_node.attribute_name = "SPH_color"
        # Link nodes
        links.new(vc_node.outputs["Color"], bsdf_node.inputs["Base Color"])
        links.new(bsdf_node.outputs["BSDF"], out_node.inputs["Surface"])

        # Assign material to object
        obj.data.materials.clear()
        obj.data.materials.append(mat)

        # Refresh colors for current frame
        frame_update_handler(context.scene)
        context.scene.frame_set(context.scene.frame_current)

        self.report({'INFO'}, f"Applied vertex color material '{mat_name}' and refreshed colors")
        return {'FINISHED'}

# Operator to refresh color field without rerunning simulation
class SPH_OT_UpdateColorField(bpy.types.Operator):
    bl_idname = "sphinxsys.update_color"
    bl_label = "Refresh Color"
    bl_description = "Refresh vertex colors based on selected field (Density/Velocity)"

    @classmethod
    def poll(cls, context):
        return bpy.data.objects.get(water_obj_name) is not None

    def execute(self, context):
        frame_update_handler(context.scene)
        context.scene.frame_set(context.scene.frame_current)
        self.report({'INFO'}, "Refreshed SPH particle colors")
        return {'FINISHED'}

def register():
    bpy.utils.register_class(SPH_OT_ApplySPHColorMaterial)
    bpy.utils.register_class(SPH_OT_UpdateColorField)

def unregister():
    bpy.utils.unregister_class(SPH_OT_UpdateColorField)
    bpy.utils.unregister_class(SPH_OT_ApplySPHColorMaterial)
