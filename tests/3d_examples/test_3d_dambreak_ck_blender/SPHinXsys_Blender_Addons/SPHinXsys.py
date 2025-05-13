import bpy
import os
import sys
import glob
import traceback

from . import props
import pyvista as pv       # requires `pip install vtk pyvista`
import numpy as np

# --- Add the add-on folder to the DLL/Python search path (Windows) ---
plugin_dir = os.path.dirname(os.path.abspath(__file__))
try:
    os.add_dll_directory(plugin_dir)
except Exception:
    pass
if plugin_dir not in sys.path:
    sys.path.insert(0, plugin_dir)

# --- Try to import the SPHinXsys core module ---
try:
    import test_3d_dambreak_ck_blender as SPHinXsys
    isSPHinXsysInstalled = True
except Exception:
    traceback.print_exc()
    isSPHinXsysInstalled = False

# --- Global state ---
water_vtp_files = []
wall_obj_created = False
water_obj_name = "WaterBodyParticles"
wall_obj_name  = "WallParticles"

@bpy.app.handlers.persistent
def frame_update_handler(scene):
    """On each frame change: reload VTP, update point positions, SPH_value and vertex color."""
    if not water_vtp_files:
        return

    idx = scene.frame_current - 1
    if idx < 0 or idx >= len(water_vtp_files):
        return

    fp = water_vtp_files[idx]
    try:
        grid = pv.read(fp)
        pts  = np.array(grid.points)
        field = scene.m_SPHinXsysGlobalSettings.m_ColorField
        arr   = grid.point_data[field]
        # If vector field, take magnitude
        if arr.ndim == 2:
            arr = np.linalg.norm(arr, axis=1)
        else:
            arr = np.array(arr)
    except Exception as e:
        print(f"[SPHinXsys handler] VTP read error: {e}")
        return

    obj = bpy.data.objects.get(water_obj_name)
    if not obj or obj.type != 'MESH':
        return
    mesh = obj.data

    # 1) Update positions
    for i, v in enumerate(mesh.vertices):
        v.co = pts[i]

    # 2) Write into SPH_value attribute
    val_attr = mesh.attributes.get("SPH_value")
    if val_attr:
        for i, d in enumerate(val_attr.data):
            d.value = float(arr[i])

    # 3) Write into vertex-color layer "SPH_color"
    vmin = float(arr.min())
    vmax = float(arr.max())
    span = vmax - vmin
    if span == 0.0:
        span = 1.0  # Prevent division by zero

    col_attr = mesh.attributes.get("SPH_color")
    if col_attr:
        for i, d in enumerate(col_attr.data):
            t = (arr[i] - vmin) / span
            d.color = (t, 0.0, 1.0 - t, 1.0)
    mesh.update()

class SPHinXsys_OT_runCase(bpy.types.Operator):
    """Read a folder of VTP files and animate them as a particle cloud."""
    bl_idname  = "sphinxsys.run_case"
    bl_label   = "Run Case"
    bl_options = {'REGISTER'}

    def execute(self, context):
        global water_vtp_files, wall_obj_created

        # 1) Ensure the .blend is saved
        blend_fp = bpy.data.filepath
        if not blend_fp:
            self.report({'ERROR'}, "Please save the .blend file first")
            return {'CANCELLED'}
        blend_dir = os.path.dirname(blend_fp)
        os.chdir(blend_dir)

        # 2) Run SPHinXsys simulation
        if not isSPHinXsysInstalled:
            self.report({'ERROR'}, "Failed to load SPHinXsys core module")
            return {'CANCELLED'}
        try:
            P  = context.scene.m_SPHinXsysGlobalSettings
            bp = SPHinXsys.BasicParameters()
            bp.sim_domain_lower        = P.m_SimDomainLower
            bp.sim_domain_upper        = P.m_SimDomainUpper
            bp.particle_spacing_ref    = P.m_ParticleReferenceSpacing
            bp.gravity_g               = P.m_Gravity
            bp.end_time                = P.m_EndTime
            bp.rho0_f                  = P.m_Rho0f
            bp.U_ref                   = P.m_URef
            bp.water_block_file_path   = P.m_WaterBlockFilePath
            bp.rigid_block_file_path   = P.m_RigidBlockFilePath
            bp.elastic_block_file_path = P.m_ElasticBlockFilePath
            bp.water_block_file_name   = props.m_FluidBlocksName   or ["water.stl"]
            bp.rigid_block_file_name   = props.m_RigidBlocksName   or ["rigid.stl"]
            bp.elastic_block_file_name = props.m_ElasticBlocksName or ["elastic.stl"]
            sim = SPHinXsys.Simulator(bp)
            sim.run()
        except Exception as e:
            traceback.print_exc()
            self.report({'ERROR'}, f"Simulation failed: {e}")
            return {'CANCELLED'}

        # 3) Collect output VTPs and build the particle mesh
        outdir  = os.path.join(blend_dir, "output")
        pattern = os.path.join(outdir, "WaterBody_*.vtp")
        water_vtp_files = sorted(glob.glob(pattern))
        if not water_vtp_files:
            self.report({'WARNING'}, "No WaterBody_*.vtp files found")
        else:
            g0  = pv.read(water_vtp_files[0])
            pts = np.array(g0.points)
            mesh = bpy.data.meshes.new(water_obj_name)
            mesh.from_pydata(pts.tolist(), [], [])

            # Create SPH_value (scalar) and SPH_color (vertex color) layers
            mesh.attributes.new(name="SPH_value",   type='FLOAT',       domain='POINT')
            mesh.attributes.new(name="SPH_color",   type='FLOAT_COLOR', domain='POINT')
            mesh.update()

            # Remove old object if present
            if water_obj_name in bpy.data.objects:
                old = bpy.data.objects[water_obj_name]
                bpy.data.objects.remove(old, do_unlink=True)

            obj = bpy.data.objects.new(water_obj_name, mesh)
            context.collection.objects.link(obj)

            # Set frame range to match number of files
            scene = context.scene
            scene.frame_start = 1
            scene.frame_end   = len(water_vtp_files)

        # 4) Load wall particles on first run
        if not wall_obj_created:
            wf = os.path.join(outdir, "Wall_0000000000.vtp")
            if os.path.isfile(wf):
                try:
                    gw = pv.read(wf)
                    pw = np.array(gw.points)
                    mesh2 = bpy.data.meshes.new(wall_obj_name)
                    mesh2.from_pydata(pw.tolist(), [], [])
                    mesh2.update()

                    if wall_obj_name in bpy.data.objects:
                        old2 = bpy.data.objects[wall_obj_name]
                        bpy.data.objects.remove(old2, do_unlink=True)

                    obj2 = bpy.data.objects.new(wall_obj_name, mesh2)
                    context.collection.objects.link(obj2)
                except Exception as e:
                    print(f"[SPHinXsys] Failed to load wall particles: {e}")
                    self.report({'WARNING'}, f"Failed to load wall particles: {e}")
            else:
                self.report({'WARNING'}, f"Wall particle file not found: {wf}")
            wall_obj_created = True

        # 5) Register frame handler and start playback
        handlers = bpy.app.handlers.frame_change_post
        if frame_update_handler not in handlers:
            handlers.append(frame_update_handler)
        frame_update_handler(context.scene)
        try:
            bpy.ops.screen.animation_play()
        except Exception:
            pass

        self.report({'INFO'}, "Setup complete, playing animation…")
        return {'FINISHED'}
