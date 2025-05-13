import bpy

class SPHinXsys_PT_Panel(bpy.types.Panel):
    bl_idname = "SPHinXsys_PT_Panel"
    bl_label = "SPHinXsys"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "SPHinXsys"

    def draw(self, context):
        layout = self.layout

        col = layout.column()
        col.operator('mesh.clearstlfiles', text="Clear STL Files")
        col.operator('mesh.saveselectedmeshes2file', text="Save Selected Meshes")
        col.operator('mesh.saveallmeshes2file', text="Save All Meshes")
        # Modified Run Case button to match new Operator ID
        col.operator('sphinxsys.run_case', text='Run Case', icon='PLAY')
        # Add Show Color, Refresh Color and Setup Particle Viz buttons
        col.operator('sphinxsys.apply_color_material', text='Show Color', icon='COLOR')
        col.operator('sphinxsys.update_color', text='Refresh Color', icon='FILE_REFRESH')
        col.operator('sphinxsys.setup_particle_viz', text='Setup Particle Viz', icon='PARTICLES')

        col = layout.column()
        gs = context.scene.m_SPHinXsysGlobalSettings
        col.prop(gs, 'm_WaterBlockFilePath')
        col.prop(gs, 'm_RigidBlockFilePath')
        col.prop(gs, 'm_ElasticBlockFilePath')

        col.prop(gs, 'm_ColorField')

        split = layout.split(factor=0.5)
        split.prop(gs, 'm_SimDomainLower')
        split.prop(gs, 'm_SimDomainUpper')

        col.prop(gs, 'm_Gravity')
        col.prop(gs, 'm_EndTime')

        col.prop(gs, 'm_ParticleReferenceSpacing')
        col.prop(gs, 'm_Rho0f')
        col.prop(gs, 'm_URef')

        # Display SPHinXsys media settings for each mesh object
        for mesh in bpy.data.meshes:
            box = layout.box()
            box.use_property_split = True
            box.use_property_decorate = False
            box.label(text=mesh.name)

            box.prop(mesh.m_SPHinXsysSettings, 'm_TypeofMedia', text='Media Type')
            box.prop(mesh.m_SPHinXsysSettings, 'm_Export', text='Exported?')

            if mesh.m_SPHinXsysSettings.m_TypeofMedia == 'Fluid':
                for prop in dir(mesh.m_SPHinXsysSettings):
                    if prop.startswith('m_Fluid'):
                        box.prop(mesh.m_SPHinXsysSettings, prop)
            else:
                for prop in dir(mesh.m_SPHinXsysSettings):
                    if prop.startswith('m_Solid'):
                        box.prop(mesh.m_SPHinXsysSettings, prop)

def register():
    bpy.utils.register_class(SPHinXsys_PT_Panel)

def unregister():
    bpy.utils.unregister_class(SPHinXsys_PT_Panel)
