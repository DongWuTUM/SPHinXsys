import bpy

# —— 1. 创建或获取小球原型 ParticleProto ——
if 'ParticleProto' not in bpy.data.objects:
    bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=3, radius=1)
    particle = bpy.context.active_object
    particle.name = 'ParticleProto'
    # 平滑着色
    for poly in particle.data.polygons:
        poly.use_smooth = True
else:
    particle = bpy.data.objects['ParticleProto']

# —— 2. 创建或获取材质 M_Particle ——
if 'M_Particle' not in bpy.data.materials:
    mat = bpy.data.materials.new(name='M_Particle')
    mat.use_nodes = True
    nodes = mat.node_tree.nodes
    links = mat.node_tree.links
    nodes.clear()
    # 输出节点
    out = nodes.new('ShaderNodeOutputMaterial'); out.location = (300, 0)
    # BSDF
    bsdf = nodes.new('ShaderNodeBsdfPrincipled'); bsdf.location = (0, 0)
    # Attribute 节点读取 SPH_color
    attr = nodes.new('ShaderNodeAttribute'); attr.location = (-200, 0)
    attr.attribute_name = 'SPH_color'
    links.new(attr.outputs['Color'], bsdf.inputs['Base Color'])
    links.new(bsdf.outputs['BSDF'], out.inputs['Surface'])
else:
    mat = bpy.data.materials['M_Particle']

# 隐藏并赋材质给原型
particle.data.materials.clear()
particle.data.materials.append(mat)
particle.hide_viewport = True
particle.hide_render = True

# —— 3. 获取点云对象 ——
obj = bpy.data.objects.get('WaterBodyParticles')
if obj is None:
    raise Exception("找不到对象 WaterBodyParticles，请确认名称一致")

# —— 4. 新建 Geometry Nodes 节点树 ——

# 创建节点树
gn = bpy.data.node_groups.new('GN_Particles', 'GeometryNodeTree')

# —— 5. 创建组接口插槽 ——
gn.interface.new_socket(name='Geometry', in_out='INPUT',  socket_type='NodeSocketGeometry')
gn.interface.new_socket(name='Scale',    in_out='INPUT',  socket_type='NodeSocketFloat')
gn.interface.new_socket(name='Geometry', in_out='OUTPUT', socket_type='NodeSocketGeometry')

# 添加组输入/输出节点
gn_in  = gn.nodes.new('NodeGroupInput');  gn_in.location  = (-600,  0)
gn_out = gn.nodes.new('NodeGroupOutput'); gn_out.location = ( 600,  0)

# —— 6. 捕获颜色属性 ——
cap_color = gn.nodes.new('GeometryNodeCaptureAttribute'); cap_color.location = (-300,  0)
cap_color.domain = 'POINT'
cap_color.capture_items.new('RGBA', 'SPH_color')

# —— 7. 实例化 & 原型信息 ——
inst = gn.nodes.new('GeometryNodeInstanceOnPoints'); inst.location = (0,   0)
obj_info = gn.nodes.new('GeometryNodeObjectInfo'); obj_info.location = (-300, -200)
obj_info.inputs['Object'].default_value = particle
obj_info.inputs['As Instance'].default_value = True

# —— 8. 实现化 & 赋材质 ——
real = gn.nodes.new('GeometryNodeRealizeInstances'); real.location = (300,  0)
setm = gn.nodes.new('GeometryNodeSetMaterial');        setm.location = (300, -200)
setm.inputs['Material'].default_value = mat

# —— 9. 连接节点 ——
links = gn.links
links.new(gn_in.outputs['Geometry'],    cap_color.inputs['Geometry'])
links.new(cap_color.outputs['Geometry'], inst.inputs['Points'])
links.new(obj_info.outputs['Geometry'], inst.inputs['Instance'])
links.new(gn_in.outputs['Scale'],       inst.inputs['Scale'])
links.new(inst.outputs[0],              real.inputs[0])
links.new(real.outputs[0],              setm.inputs['Geometry'])
links.new(setm.outputs[0],              gn_out.inputs['Geometry'])

# —— 10. 应用 Modifier ——
mod = obj.modifiers.new('GN_Particles', 'NODES')
mod.node_group = gn

print("Geometry Nodes setup complete (uniform scale) ✅")
