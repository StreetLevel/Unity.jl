import Unity: UnityMesh, UnityText
using ColorTypes, GeometryTypes

#define vertices
verts = [Point3f0(0,0,0),Point3f0(1,0,0),Point3f0(1,1,0),Point3f0(0,1,0)]
#define indices for surface mesh
surf_inds = Int32[0,2,1,0,3,2]
#define indices for line mesh
line_inds = Int32[0,1,1,2,2,3,3,0]
#define indices for vertex mesh
vert_inds = Int32[0,1,2,3]
#define colors
colors = rand(RGBA{Float32},length(verts))
#colors = [RGBA{Float32}(1,1,1,.001) for vert in verts]

#define optons
options = ["surface_shader = flat"]
transparent = ["surface_shader = transparent"]
point_size = ["point_size = 0.5"]

#define text
unity_text = UnityText("Hello World!", Point3f0(1.5,2,0), Point3f0(0.1,0.1,.15), Point3f0(1,1,1), rand(RGBA{Float32}))

#create meshes
surface_mesh = UnityMesh("Mesh 1", verts, Int32[], Int32[], surf_inds, colors, options, UnityText[unity_text])
line_mesh = UnityMesh("Mesh 2", map(x->Point3f0(x[1]+3,x[2],x[3]),verts), Int32[], line_inds, Int32[], colors, String[])
point_mesh = UnityMesh("Mesh 3", map(x->Point3f0(x[1],x[2]+3,x[3]),verts), vert_inds, Int32[], Int32[], colors, point_size)
mesh = UnityMesh("Mesh 4", map(x->Point3f0(x[1]+3,x[2]+3,x[3]),verts), vert_inds, line_inds, surf_inds, colors, vcat(point_size,transparent))

#send meshes to unity
socket = connect(8052)
write(socket, surface_mesh)
write(socket, line_mesh)
write(socket, point_mesh)
write(socket, mesh)
close(socket)