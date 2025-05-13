import bpy

class SPHinXsysGlobalProps(bpy.types.PropertyGroup):
    m_SimDomainLower: bpy.props.FloatVectorProperty(
        name='Simulation Domain Lower',
        description='Lower bounds of the simulation domain',
        default=(-1.0, -1.0, -1.0)
    )
    m_SimDomainUpper: bpy.props.FloatVectorProperty(
        name='Simulation Domain Upper',
        description='Upper bounds of the simulation domain',
        default=(1.0, 1.0, 1.0)
    )
    m_ParticleReferenceSpacing: bpy.props.FloatProperty(
        name='Fluid Reference Spacing',
        description='Reference spacing between fluid particles',
        default=0.05
    )
    m_Gravity: bpy.props.FloatProperty(
        name='Gravity',
        description='Gravity acceleration',
        default=1.0
    )
    m_EndTime: bpy.props.FloatProperty(
        name='End Time',
        description='End time of the simulation',
        default=20.0
    )
    m_Rho0f: bpy.props.FloatProperty(
        name='Reference Density',
        description='Reference density for fluid',
        default=1.0
    )
    m_URef: bpy.props.FloatProperty(
        name='Reference Velocity',
        description='Reference velocity for fluid',
        default=20.0
    )
    m_Cf: bpy.props.FloatProperty(
        name='Reference Sound Speed',
        description='Reference speed of sound for fluid',
        default=200.0
    )
    m_WaterBlockFilePath: bpy.props.StringProperty(
        name='Water Block File Path',
        description='Path to water block STL files',
        default="./Mesh/Fluids/"
    )
    m_RigidBlockFilePath: bpy.props.StringProperty(
        name='Rigid Block File Path',
        description='Path to rigid block STL files',
        default="./Mesh/Rigids/"
    )
    m_ElasticBlockFilePath: bpy.props.StringProperty(
        name='Elastic Block File Path',
        description='Path to elastic block STL files',
        default="./Mesh/Elastics/"
    )

    # Color field selection (Density or Velocity)
    m_ColorField: bpy.props.EnumProperty(
        name='Color Field',
        description='Field used for particle coloring',
        items=[
            ('Density',  'Density',  'Color by density'),
            ('Velocity', 'Velocity', 'Color by velocity magnitude'),
        ],
        default='Density'
    )


class SPHinXsysProps(bpy.types.PropertyGroup):
    m_TypeofMedia: bpy.props.EnumProperty(
        name='Media Type',
        description='Type of media for this mesh',
        items=[
            ('RigidSolid',   'Rigid Solid',   'Rigid solid'),
            ('ElasticSolid', 'Elastic Solid', 'Elastic solid'),
            ('Fluid',        'Fluid',         'Fluid'),
        ]
    )
    m_Export: bpy.props.BoolProperty(
        name='Export',
        description='Whether to export this mesh',
        default=True
    )
    # Fluid-specific properties
    m_FluidDensity: bpy.props.FloatProperty(name='Fluid Density', default=1.0)
    m_FluidCharacteristicVelocity: bpy.props.FloatProperty(name='Characteristic Velocity', default=10.0)
    m_FluidSpeedofSound: bpy.props.FloatProperty(name='Speed of Sound', default=10.0)
    m_FluidReynoldsNumber: bpy.props.FloatProperty(name='Reynolds Number', default=10.0)
    m_FluidDynamicsViscosity: bpy.props.FloatProperty(name='Dynamic Viscosity', default=10.0)
    # Solid-specific properties
    m_SolidReferenceDensity: bpy.props.FloatProperty(name='Reference Density', default=10.0)
    m_SolidYoungsModulus: bpy.props.FloatProperty(name='Young’s Modulus', default=1.4e3)
    m_SolidPoissonRatio: bpy.props.FloatProperty(name='Poisson Ratio', default=0.4)


def registerSPHinXsysProperties():
    bpy.types.Scene.m_SPHinXsysGlobalSettings = bpy.props.PointerProperty(
        name="SPHinXsys Settings",
        type=SPHinXsysGlobalProps
    )
    bpy.types.Mesh.m_SPHinXsysSettings = bpy.props.PointerProperty(
        name="SPHinXsys Mesh Settings",
        type=SPHinXsysProps
    )


def unregisterSPHinXsysProperties():
    del bpy.types.Scene.m_SPHinXsysGlobalSettings
    del bpy.types.Mesh.m_SPHinXsysSettings

# The following functions maintain lists of block names for Save Mesh operations
m_FluidBlocksName   = []
m_RigidBlocksName   = []
m_ElasticBlocksName = []

def appendBlocksName(blockType: str, blockName: str):
    """Add a mesh name to the corresponding block list."""
    if blockType == 'Fluid':
        if blockName not in m_FluidBlocksName:
            m_FluidBlocksName.append(blockName)
    elif blockType == 'RigidSolid':
        if blockName not in m_RigidBlocksName:
            m_RigidBlocksName.append(blockName)
    elif blockType == 'ElasticSolid':
        if blockName not in m_ElasticBlocksName:
            m_ElasticBlocksName.append(blockName)

def deleteBlocksName(blockType: str, blockName: str):
    """Remove a mesh name from the corresponding block list."""
    if blockType == 'Fluid' and blockName in m_FluidBlocksName:
        m_FluidBlocksName.remove(blockName)
    elif blockType == 'RigidSolid' and blockName in m_RigidBlocksName:
        m_RigidBlocksName.remove(blockName)
    elif blockType == 'ElasticSolid' and blockName in m_ElasticBlocksName:
        m_ElasticBlocksName.remove(blockName)

def clearAllBlocksName():
    """Clear all block name lists."""
    m_FluidBlocksName.clear()
    m_RigidBlocksName.clear()
    m_ElasticBlocksName.clear()
