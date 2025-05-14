import bpy
from .SPHinXsys import water_obj_name

class SPH_OT_SetupParticleVisualization(bpy.types.Operator):
    """Instance point cloud as colored spheres via Geometry Nodes"""
    bl_idname = "sphinxsys.setup_particle_viz"
    bl_label = "Setup Particle Viz"
    bl_description = "Use Geometry Nodes to instance spheres on particles with vertex colors"

    @classmethod
    def poll(cls, context):
        return bpy.data.objects.get(water_obj_name) is not None

    def execute(self, context):
        # 1. Create or get sphere prototype
        if 'ParticleProto' not in bpy.data.objects:
            bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=1, radius=1)
            proto = bpy.context.active_object
            proto.name = 'ParticleProto'
            for p in proto.data.polygons:
                p.use_smooth = True
        else:
            proto = bpy.data.objects['ParticleProto']

        # 2. Create or get material
        mat_name = 'M_Particle'
        if mat_name not in bpy.data.materials:
            mat = bpy.data.materials.new(mat_name)
            mat.use_nodes = True
            nodes = mat.node_tree.nodes
            links = mat.node_tree.links
            nodes.clear()
            out = nodes.new('ShaderNodeOutputMaterial'); out.location = (300, 0)
            bsdf = nodes.new('ShaderNodeBsdfPrincipled'); bsdf.location = (0, 0)
            attr = nodes.new('ShaderNodeAttribute'); attr.location = (-200, 0)
            attr.attribute_name = 'SPH_color'
            links.new(attr.outputs['Color'], bsdf.inputs['Base Color'])
            links.new(bsdf.outputs['BSDF'], out.inputs['Surface'])
        else:
            mat = bpy.data.materials[mat_name]

        proto.data.materials.clear()
        proto.data.materials.append(mat)
        proto.hide_viewport = True
        proto.hide_render = True

        # 3. Get point cloud object
        obj = bpy.data.objects.get(water_obj_name)
        if not obj:
            self.report({'ERROR'}, f"Object '{water_obj_name}' not found")
            return {'CANCELLED'}

        # 4. Build Geometry Nodes node tree
        gn = bpy.data.node_groups.new(name='GN_Particles', type='GeometryNodeTree')
        # Node group interface
        gn.interface.new_socket(name='Geometry', in_out='INPUT', socket_type='NodeSocketGeometry', description='')
        gn.interface.new_socket(name='Scale',    in_out='INPUT', socket_type='NodeSocketFloat',    description='')
        gn.interface.new_socket(name='Geometry', in_out='OUTPUT',socket_type='NodeSocketGeometry', description='')
        gn_in  = gn.nodes.new('NodeGroupInput');  gn_in.location  = (-600, 0)
        gn_out = gn.nodes.new('NodeGroupOutput'); gn_out.location = (600,  0)

        # Nodes
        cap     = gn.nodes.new('GeometryNodeCaptureAttribute');  cap.location = (-300,  0)
        cap.domain = 'POINT'; cap.capture_items.new('RGBA', 'SPH_color')
        inst    = gn.nodes.new('GeometryNodeInstanceOnPoints');  inst.location = (0,    0)
        obj_info= gn.nodes.new('GeometryNodeObjectInfo');        obj_info.location = (-300, -200)
        obj_info.inputs['Object'].default_value     = proto
        obj_info.inputs['As Instance'].default_value= True
        real    = gn.nodes.new('GeometryNodeRealizeInstances'); inst.location = (300,  0)
        setm    = gn.nodes.new('GeometryNodeSetMaterial');        setm.location = (300, -200)
        setm.inputs['Material'].default_value = mat

        # Links
        links = gn.links
        links.new(gn_in.outputs['Geometry'], cap.inputs['Geometry'])
        links.new(cap.outputs['Geometry'], inst.inputs['Points'])
        links.new(obj_info.outputs['Geometry'], inst.inputs['Instance'])
        links.new(gn_in.outputs['Scale'],   inst.inputs['Scale'])
        links.new(inst.outputs[0],          real.inputs[0])
        links.new(real.outputs[0],          setm.inputs['Geometry'])
        links.new(setm.outputs[0],          gn_out.inputs['Geometry'])

        # 5. Apply modifier to the point cloud object
        mod = obj.modifiers.new(name='GN_Particles', type='NODES')
        mod.node_group = gn

        self.report({'INFO'}, 'Particle visualization setup complete')
        return {'FINISHED'}

def register():
    bpy.utils.register_class(SPH_OT_SetupParticleVisualization)

def unregister():
    bpy.utils.unregister_class(SPH_OT_SetupParticleVisualization)
