include("/Users/maximilianbittens/Documents/Unity/Unity.jl/src/MeshReaders.jl")
using GeometryTypes
using ColorTypes
import Unity: UnityMesh, PyramidMesh, UnityVector3, UnityText
import MeshReaders: GmshReader, load

file1 = "/Users/maximilianbittens/Documents/IBNM/Dissertation/Modelle/Knochen\ heile.msh"
file2 = "/Users/maximilianbittens/Documents/IBNM/Dissertation/Modelle/Mayo.msh"

bone_mesh = MeshReaders.load(GmshReader, file1)
mayo_mesh = MeshReaders.load(GmshReader, file2)
bone_triangles = vcat(map(x->x[4:end],filter(x->x[1:3]==[2,3,1],bone_mesh.elements))-1...)
mayo_triangles = vcat(map(x->x[4:end],filter(x->x[1:3]==[2,3,1],mayo_mesh.elements))-1...)
bone_lines = vcat(map(x->x[4:end],filter(x->x[1:3]==[1,3,1],bone_mesh.elements))-1...)
mayo_lines = vcat(map(x->x[4:end],filter(x->x[1:3]==[1,3,1],mayo_mesh.elements))-1...)
bone_vertices = UnityVector3[UnityVector3(bone_mesh.nodes[1,i],bone_mesh.nodes[2,i],bone_mesh.nodes[3,i]) for i = 1:size(bone_mesh.nodes,2)]
#mayo_vertices = UnityVector3[UnityVector3(mayo_mesh.nodes[1,i],mayo_mesh.nodes[2,i],mayo_mesh.nodes[3,i]) for i = 1:size(mayo_mesh.nodes,2)]
mayo_vertices = [Point3f0(mayo_mesh.nodes[1,i],mayo_mesh.nodes[2,i],mayo_mesh.nodes[3,i]) for i = 1:size(mayo_mesh.nodes,2)]
#mayo_vertices = map(x->pcamat*x,mayo_vertices)
bone_color = [RGBA{Float32}(.8,.8,.8,.3) for vert in bone_vertices]
mayo_color = [RGBA{Float32}(.8,0,.3,.3) for vert in mayo_vertices]
#mayo_color = [rand(RGBA{Float32}) for vert in mayo_vertices]
#mayo_color = [rand(RGBA{Float32}) for i = 1:3:length(mayo_triangles)]
bone_options = ["surface_shader = wireframe"]
mayo_options = ["surface_shader = transparent"]


unity_text = UnityText("",UnityVector3(-30.,400.,-1418.),UnityVector3(5,5,1.),UnityVector3(90.,180.,1.))

unity_bone_mesh = UnityMesh("Bone",bone_vertices,Int32[],bone_lines,bone_triangles,bone_color,bone_options,UnityText[unity_text])
unity_mayo_mesh = UnityMesh("Mayo",mayo_vertices,Int32[],mayo_lines,mayo_triangles,mayo_color,mayo_options)
unity_mayo_mesh = UnityMesh("Mayo",mayo_vertices,Int32[],Int32[],mayo_triangles,mayo_color,mayo_options)



function plot(um::UnityMesh)
	socket = connect(8052)
	write(socket,um)
	close(socket)
end

function rand!(um::UnityMesh)
	
	for vert in um.vertices
		vert.x += randn()/2
		vert.y += randn()/2
		vert.z += randn()/2
	end
	um.colors = rand(RGBA{Float32},length(um.vertices))
end


#sleep(1)
#socket = connect(8052)
#for i = 1:10
#	unity_bone_mesh.text[1].text = "$i"
#	rand!(unity_bone_mesh)
#	#plot(unity_bone_mesh)
#	write(socket,unity_bone_mesh)
#	sleep(.75)
#end
#close(socket)


