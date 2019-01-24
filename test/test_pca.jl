using GeometryTypes

const AInt = Int32  
const AFloat = Float32

include("/Users/maximilianbittens/.julia/v0.6/Unity/src/MeshReaders.jl")
include("/Users/maximilianbittens/.julia/v0.6/Unity/src/AnsysReader.jl")

using GeometryTypes
using ColorTypes
import Unity: UnityMesh, PyramidMesh, UnityVector3, UnityText
using Colors

function triangle_area(vertices::Vector{Point3f0}, faces::Vector{Face{3,AInt}})
	areas = Vector{AFloat}();
	for face in faces
		push!(areas,triangle_area(vertices[face]...))
	end
	return areas
end

function triangle_area(P1::Point, P2::Point, P3::Point)
	L = ( 	
		sqrt(sum((P1-P2).^2)), 
		sqrt(sum((P2-P3).^2)), 
		sqrt(sum((P3-P1).^2)) 
		)
	s = sum(L)/2.0
	return sqrt( s*(s-L[1])*(s-L[2])*(s-L[3]) )
end

function centroids(vertices::Vector{Point3f0}, faces::Vector{Face{3,AInt}})
	#c = Vector{Point{3,Float64}}()
	c = Vector{Point{3,Float32}}()
	for face in faces
		vs = vertices[face]
		#push!(c,sum(convert(Tuple{Point{3,Float64},Point{3,Float64},Point{3,Float64}},vs))/3.)
		push!(c,sum(vs)/3.)
	end
	return c
end

#
function vertices_and_faces(mesh)
	vertices = Vector{Point3f0}()
	for i = 1:size(mesh.nodes,2)
		push!(vertices,Point3f0(mesh.nodes[:,i]))
	end
	faces = map(x->Face{3,AInt}(x[4:end]),filter(x->x[1:3]==[2,3,1],mesh.elements))
	return vertices,faces
end
#


function covariance_matrix_and_weighted_center(mesh)
	vertices, faces = vertices_and_faces(mesh)
	return covariance_matrix_and_weighted_center(vertices, faces)
end

function covariance_matrix_and_weighted_center(vertices::Vector{Point3f0}, faces::Vector{Face{3,AInt}})
	face_centroids = centroids(vertices,faces)
	info("face_centroids")
	display(face_centroids)
	tri_areas = triangle_area(vertices,faces)
	weighted_center = sum(face_centroids.*tri_areas)/sum(tri_areas)
	#info("weighted_center")
	#display(weighted_center)
	zero_centered_centroids = map((x,y)->(x-weighted_center)*sqrt(y),face_centroids,tri_areas)
	#info("zero_centeredweighted_center")
	#display(zero_centered_centroids)
	wcm = hcat(map(x->[x...],zero_centered_centroids)...)
	#info("wcm")
	#display(wcm)
	return (wcm*wcm')./(length(zero_centered_centroids)-1), weighted_center
end

#
function minarc(v1,coords=[Point3f0(1,0,0),Point3f0(0,1,0),Point3f0(0,0,1)])
	arc(v1,v2) = acos( dot(v1,v2)/(norm(v1)*norm(v2)))
	a = [ arc(v1,coords[1]),arc(v1,coords[2]),arc(v1,coords[3])]
	ind = a[1]<min(a[2],a[3]) ? 1 : (a[2]<min(a[1],a[3]) ? 2 : ( (a[3]<min(a[1],a[2]) ? 3 : error() )))
	return a[ind], ind
end

function pca_weighted(mesh)
	vertices, faces = vertices_and_faces(mesh)
	return pca_weighted(vertices,faces)
end
#

function pca_weighted(vertices::Vector{Point3f0}, faces::Vector{Face{3,AInt}})
	try
	cov_mat, wc = covariance_matrix_and_weighted_center(vertices,faces)
	#info("cov")
	#display(cov_mat)
	λ, ϕ = eigs(cov_mat,nev=2)
	#display(ϕ)
	#pca = align(Point3f0(ϕ[:,1]/norm(ϕ[:,1])), Point3f0(ϕ[:,2]/norm(ϕ[:,2])))
	pca = [Point3f0(ϕ[:,1]/norm(ϕ[:,1])),Point3f0(ϕ[:,2]/norm(ϕ[:,2])),Point3f0(cross(ϕ[:,1],ϕ[:,2])/norm(cross(ϕ[:,1],ϕ[:,2])))]
	return pca,wc
	catch
		display(faces)
		error()
	end
end

function mul!(res::Vector{Float32},mat::StaticArrays.SArray{Tuple{3,3},Float32,2,9},vec::Point3f0)
	fill!(res,0.0f0)
	for i = 1:3
		for j = 1:3
			res[i] += mat[i,j]*vec[j]
		end
	end
end

import StaticArrays

type OBB #Oriented Bounding Box
	pca::StaticArrays.SArray{Tuple{3,3},Float32,2,9}
	minb::Vector{Float32}
	maxb::Vector{Float32}
	OBB(pca,min,max) = new(pca,min,max)
	OBB(pca) = new(hcat(pca...),[Inf32,Inf32,Inf32],[-Inf32,-Inf32,-Inf32])
end

function Base.copy(obb::OBB)
	return OBB(copy(obb.pca),copy(obb.minb),copy(obb.maxb))
end

function oriented_bounding_box(vertices::Vector{Point3f0},pca::Vector{Point3f0})
	obb = OBB(pca)
	verts = map(x->inv(obb.pca)*x,vertices)
	obb.minb,obb.maxb = boundingbox(verts)
	return obb
end

function boundingbox(vertices::Vector{Point3f0})
	return Point3f0(minimum(map(x->x[1],vertices)),minimum(map(x->x[2],vertices)),minimum(map(x->x[3],vertices))),Point3f0(maximum(map(x->x[1],vertices)),maximum(map(x->x[2],vertices)),maximum(map(x->x[3],vertices)))
end

function boundingbox_vertices(minb,maxb)
	Δ = maxb-minb
	vertices = Point3f0[minb, minb+[Δ[1],0,0],minb+[0,Δ[2],0],minb+[0,0,Δ[3]],minb+[Δ[1],Δ[2],0],minb+[Δ[1],0,Δ[3]],minb+[0,Δ[2],Δ[3]],maxb]
	return vertices
end

function GeometryTypes.volume(obb::OBB)
	Δ = obb.maxb-obb.minb
	return reduce(*,Δ)
end

#
function pos(P::Union{Point3f0,Vector{Float32}},A,n)
	return dot(P-A,n)>0
end

function pos(P::Vector{Point3f0},A,n)
	return all(map(x->pos(x,A,n),P))
end
#

function inside(obb::OBB, x::Union{Point3f0,Vector{Float32}})
	#return all(obb.minb.<(x+Float32(eps()))) && all(obb.maxb.>(x-Float32(eps())))
	Δ = obb.maxb-obb.minb
	return all(obb.minb.<(x+Δ/100)) && all(obb.maxb.>(x-Δ/100))
end

function inside(obb::OBB, x::Vector{Point3f0})
	return sum(map(xx->inside(obb,xx),x))>1
end

function Base.isvalid(obb::OBB)
	return all(obb.minb .< obb.maxb)
end

function inside_vertices(obb::OBB, vertices)
	inds = [inside(obb,inv(obb.pca)*vert) for vert in vertices]
	return inds,vertices[inds]
end

function inside_faces(obb::OBB, vertices, faces)
	inds = [inside(obb,[map(x->inv(obb.pca)*x,vertices[face])...]) for face in faces]
	return inds,faces[inds]
end

function oriented_bounding_box(obb::OBB,vertices,faces,vinds,finds)
	isd_fcs_inds,isd_fcs = inside_faces(obb,vertices,faces[finds])
	faceverts = unique(vcat(convert(Vector{Vector{Int32}},faces[finds[isd_fcs_inds]])...))
	isd_verts_inds,isd_verts = inside_vertices(obb,vertices[faceverts])
	if length(faceverts) > 3 && length(isd_fcs) > 1
		pca, wc = pca_weighted(vertices,isd_fcs)
		#newobb = oriented_bounding_box(isd_verts,pca)
		newobb = oriented_bounding_box(vertices[faceverts],pca)
	else
		info("length(faceverts) ")
		display(length(faceverts))
		info("length(isd_fcs) ")
		display(length(isd_fcs))
		newobb = OBB(obb.pca)
		wc = Point3f0(0)
	end
	return newobb,wc,vinds,finds[isd_fcs_inds]
end

function subdivide(obb::OBB, wc, direction)
	obb1,obb2 = copy(obb),copy(obb)
	len = (obb.maxb[direction]-obb.minb[direction])
	Δ = (inv(obb.pca)*wc)[direction] - obb.minb[direction]
	#Δ = len/2
	obb1.maxb[direction] -= (len-Δ)
	obb2.minb[direction] += Δ
	return obb1,obb2
end

function GeometryTypes.volume(elem::Tuple{Tuple{OBB,Point3f0,Vector{Int32},Vector{Int32}},Tuple{OBB,Point3f0,Vector{Int32},Vector{Int32}}})
	return volume(elem[1][1])+volume(elem[2][1])
end

function get_children(obbs::Vector{Tuple{Tuple{OBB,Point3f0,Vector{Int32},Vector{Int32}},Tuple{OBB,Point3f0,Vector{Int32},Vector{Int32}}}})
	#return obbs[1]
	#return obbs[findmin(map(x->volume(x),obbs))[2]]
	val = findmin(map(x->volume(x),obbs[1:2]))
	if abs(val[1]) < Inf32
		return obbs[val[2]]
	else
		return obbs[3]
	end
end

function refine(obb::OBB,vertices,faces,wc,vinds,finds,obbtree)
	directions = sortperm(obb.maxb-obb.minb,rev=true)
	#directions = directions[[1]]
	obbs = Vector{Tuple{Tuple{OBB,Point3f0,Vector{Int32},Vector{Int32}},Tuple{OBB,Point3f0,Vector{Int32},Vector{Int32}}}}()
	for direction in directions
		obb1,obb2 = subdivide(obb,wc,direction)

		#isd_fcs_inds1,isd_fcs1 = inside_faces(obb1,vertices,faces[finds])
		#isd_verts_inds1,isd_verts1 = inside_vertices(obb1,vertices[vinds])
		#isd_fcs_inds2,isd_fcs2 = inside_faces(obb2,vertices,faces[finds])
		#isd_verts_inds2,isd_verts2 = inside_vertices(obb2,vertices[vinds])

		obb1tmp,wc1tmp,vinds1tmp,finds1tmp = oriented_bounding_box(obb1,vertices,faces,vinds,finds)
		obb2tmp,wc2tmp,vinds2tmp,finds2tmp = oriented_bounding_box(obb2,vertices,faces,vinds,finds)

		#obbtreeelement1 = OBBTreeElement(1,obb1,Nullable{OBBTreeElement}(),wc1tmp,vinds[isd_verts_inds1],finds[isd_fcs_inds1])
		#obbtreeelement2 = OBBTreeElement(1,obb2,Nullable{OBBTreeElement}(),wc2tmp,vinds[isd_verts_inds2],finds[isd_fcs_inds2])
		#println("OBB Mesh 1")
		#plot("OBB mesh 1", obbtree, obbtreeelement1, rand(RGBA{Float32}))
		#println("OBB Mesh 2")
		#plot("OBB mesh 2", obbtree, obbtreeelement2, rand(RGBA{Float32}))
		#println("OBB 1")
		#plot("Bounding Box 1", obbtreeelement1.obb, rand(RGBA{Float32}))
		#println("OBB 2")		
		#plot("Bounding Box 2", obbtreeelement2.obb, rand(RGBA{Float32}))
		#println("done")
		#sleep(.1)

		obb1,wc1,vinds1,finds1 = obb1tmp,wc1tmp,vinds1tmp,finds1tmp 
		obb2,wc2,vinds2,finds2 = obb2tmp,wc2tmp,vinds2tmp,finds2tmp 

		push!(obbs,((obb1,wc1,vinds1,finds1),(obb2,wc2,vinds2,finds2)))
	end
	return get_children(obbs)
end

abstract type AbstractTreeElement end
type OBBTreeElement <: AbstractTreeElement
	level::Int32
	obb::OBB
	children::Nullable{Tuple{AbstractTreeElement,AbstractTreeElement}}
	parent::Nullable{AbstractTreeElement}
	wc::Point3f0
	vinds::Vector{Int32}
	finds::Vector{Int32}
	OBBTreeElement(level,obb,wc) = new(level,obb,Nullable{Tuple{AbstractTreeElement,AbstractTreeElement}}(),Nullable{AbstractTreeElement}(),wc,Vector{Int32}(),Vector{Int32}())
	OBBTreeElement(level,obb,parent,wc) = new(level,obb,Nullable{Tuple{AbstractTreeElement,AbstractTreeElement}}(),Nullable{AbstractTreeElement}(parent),wc,Vector{Int32}(),Vector{Int32}())
	OBBTreeElement(level,obb,parent,wc,vinds,finds) = new(level,obb,Nullable{Tuple{AbstractTreeElement,AbstractTreeElement}}(),Nullable{AbstractTreeElement}(parent),wc,vinds,finds)
end

type OBBTree
	maxdepth::Int32
	tree::Vector{Set{OBBTreeElement}}
	vertices::Vector{Point3f0}
	faces::Vector{Face{3,AInt}}
	OBBTree(maxdepth,vertices,faces) = new(maxdepth,Vector{Set{OBBTreeElement}}(maxdepth),vertices,faces)
end

function Base.show(io::IO, obt::OBBTree)
    print(io,"OBBTree")
end

function Base.show(io::IO, obt::OBBTreeElement)
    print(io,"OBBTreeElement Level $(obt.level)")
end

function GeometryTypes.volume(elem::OBBTreeElement)
	Δ = elem.obb.maxb-elem.obb.minb
	return reduce(*,Δ)
end

function Base.push!(obbt::OBBTree, elem::OBBTreeElement)
	push!(obbt.tree[elem.level],elem)
end

function intialize!(obbtree::OBBTree)
	for i = 1:obbtree.maxdepth; obbtree.tree[i] = Set{OBBTreeElement}(); end
	pca, wc = pca_weighted(obbtree.vertices,obbtree.faces)
	obb = oriented_bounding_box(obbtree.vertices,pca)
	elem1 = OBBTreeElement(1,obb,wc)
	elem1.vinds = collect(Int32(1):Int32(length(obbtree.vertices)))
	elem1.finds = collect(Int32(1):Int32(length(obbtree.faces)))
	push!(obbtree,elem1)
	#refine!(obbtree,elem1)
	return nothing
end

function refine!(obbtree::OBBTree, elem::OBBTreeElement)
	#((obb1,wc1,vinds1,finds1),(obb2,wc2,vinds2,finds2)) = refine(elem.obb,obbtree.vertices,obbtree.faces,elem.wc,elem.vinds,elem.finds)
	((obb1,wc1,vinds1,finds1),(obb2,wc2,vinds2,finds2)) = refine(elem.obb,obbtree.vertices,obbtree.faces,elem.wc,elem.vinds,elem.finds,obbtree)
	if isvalid(obb1) && isvalid(obb2)
		child1 = OBBTreeElement(elem.level+1,obb1,elem,wc1,vinds1,finds1)
		child2 = OBBTreeElement(elem.level+1,obb2,elem,wc2,vinds2,finds2)
		elem.children = Nullable{Tuple{OBBTreeElement,OBBTreeElement}}((child1,child2))
		push!(obbtree,child1)
		push!(obbtree,child2)
	else
		child1 = OBBTreeElement(elem.level+1,obb1,elem,wc1,vinds1,finds1)
		child2 = OBBTreeElement(elem.level+1,obb2,elem,wc2,vinds2,finds2)
	end
	return child1,child2
end

function get_highest_level_descendants(obbtr::OBBTree)
	ret = Set{OBBTreeElement}()
	tmp = Set{OBBTreeElement}()
	push!(tmp,first(obbtr.tree[1]))
	while !isempty(tmp)
		elem = pop!(tmp)
		if elem.children.hasvalue
			union!(tmp,elem.children.value)
		else
			push!(ret,elem)
		end
	end
	return ret
end

function Base.isvalid(elem::OBBTreeElement)
	return isvalid(elem.obb)
end

function Base.delete!(obbtree::OBBTree, child1::OBBTreeElement, child2::OBBTreeElement)
	if !child1.children.hasvalue && !child2.children.hasvalue
		delete!(obbtree.tree[child1.level],child1)
		delete!(obbtree.tree[child2.level],child2)
		child1.parent.value.children = Nullable{Tuple{OBBTreeElement,OBBTreeElement}}()
		child1.parent = Nullable{OBBTreeElement}()
		child2.parent = Nullable{OBBTreeElement}()
		return true
	else
		return false
	end
end

function adaptive_refine!(obbtr::OBBTree, tol::Float32)
	ret = get_highest_level_descendants(obbtr)
	while !isempty(ret)
		elem = pop!(ret)
		if !elem.children.hasvalue && elem.level < obbtr.maxdepth
			child1,child2 = refine!(obbtr,elem)
			println((volume(child1)+volume(child2))/volume(elem))
			if isvalid(child1) && isvalid(child2) && (volume(child1)+volume(child2))/volume(elem) < tol
				info("Refinement on Level: $(child1.parent.value.level)")
				push!(ret,child1)
				push!(ret,child2)
			else
				delete!(obbtr,child1,child2)
			end
		end
	end
end


function collect_parents!(ret::Set{OBBTreeElement},elem::OBBTreeElement)
	push!(ret,elem)
	if elem.parent.hasvalue
		collect_parents!(ret,elem.parent.value)
	end
	return nothing
end

function collect_parents(elem::OBBTreeElement)
	ret = Set{OBBTreeElement}()
	collect_parents!(ret,elem)
	return ret
end


function plot_parents(obbtree::OBBTree,leaf::OBBTreeElement)
	ret = collect_parents(leaf)
	for (i,elem) in enumerate(ret)
		plot("OBB fitted $i",elem.obb,RGBA{Float32}(rand(),rand(),rand(),.4))
	end
end

function plot(obbtr::OBBTree)
	ret = get_highest_level_descendants(obbtr)
	#colors = linspace(RGBA{Float32}(1,0,0,1),RGBA{Float32}(0,0,1,1),length(ret))
	colors = rand(RGB{Float32},length(ret))
	for (i,elem) in enumerate(ret)
		plot("OBB fitted $i",elem.obb,RGBA{Float32}(colors[i].r,colors[i].g,colors[i].b,.4))

		#sleep(.1)
		#plot("OBB wc $i",elem.obb.pca, elem.wc)
		
		plot("OBB mesh fitted $i", obbtr, elem, RGBA{Float32}(colors[i].r,colors[i].g,colors[i].b,.4))
	end
end

import Colors


function vertices_and_faces(vertices,faces,finds)	
	ins_fcs = faces[finds]
	ins_fcs_new = similar(faces[finds])
	vertinds = convert(Vector{Vector{Int32}},ins_fcs)
	unq_vert_inds = unique(vcat(vertinds...))
	verts = vertices[unq_vert_inds]
	for i = 1:length(ins_fcs_new)
		ins_fcs_new[i] = map(y->findfirststop(x->x==ins_fcs[i][y],unq_vert_inds),1:3)
	end
	return verts,ins_fcs_new
end


function plot(tag::String,obbtree::OBBTree,elem::OBBTreeElement,color::RGBA{Float32})
	#test = [inside(elem.obb,[map(x->inv(elem.obb.pca)*x,obbtree.vertices[face])...]) for face in obbtree.faces]
	verts,fcs = vertices_and_faces(obbtree.vertices,obbtree.faces,elem.finds)
	#colors = [color for i = 1:length(obbtree.vertices)]
	colors = [color for i = 1:length(verts)]
	options = ["surface_shader = transparent"]
	unity_faces = vcat(map(xx->convert(Vector{UInt32},xx-1),fcs)...)
	utm = UnityMesh(tag,verts,Int32[],Int32[],unity_faces,colors,options)
	socket = connect(8052)
	write(socket,utm)
	close(socket)
	#sleep(.1)
end

function plot(tag::String,obb::OBB, color::RGBA{Float32}=RGBA{Float32}(1,0,0,1))
	vertices = boundingbox_vertices(obb.minb,obb.maxb)
	info("boundingbox size: ",volume(obb))
	lines =  [0,1,0,2,0,3,3,6,3,5,2,6,1,4,5,7,1,5,4,2,4,7,7,6]
	colors = Colors.linspace(RGBA{Float32}(1,1,1,1),color,8)
	verts = map(x->obb.pca*x,vertices)
	#println(obb.minb)
	#println(obb.maxb)
	#utm = UnityMesh(tag, verts, Int32[0,7], lines, Int32[], colors, [" = "])
	utm = UnityMesh(tag, verts, Int32[], lines, Int32[], colors, [" = "])
	socket = connect(8052)
	#println(socket)
	a = write(socket,utm)
	#println(a)
	close(socket)
end



function plot(points::Vector{Point3f0})
	red = RGBA{Float32}(1,0,0)
	unity_point_mesh = UnityMesh("Point Mesh",points,collect(0:length(points)-1), Int32[], Int32[], [red for i=1:length(points)], [" = "])
	socket = connect(8052)
	write(socket,unity_point_mesh)
	close(socket)
end

function plot(tag::String, pca::StaticArrays.SArray{Tuple{3,3},Float32,2,9}, weighted_center::GeometryTypes.Point{3,Float32})
	_pca = Point3f0[pca[:,i] for i = 1:3]
	plot(tag,_pca,weighted_center)
end

function plot(tag::String,pca::Vector{Point3f0}, weighted_center::Point3f0)
	black = RGBA{Float32}(0,0,0)
	red = RGBA{Float32}(1,0,0)
	green = RGBA{Float32}(0,1,0)
	blue = RGBA{Float32}(0,0,1)
	vertices = vcat([weighted_center],map(x->weighted_center+1*x,pca))
	unity_mesh = UnityMesh(tag,vertices, Int32[0,1,2,3], Int32[0,1,0,2,0,3], Int32[], [black, red, green, blue ], ["point_size = 0.02"])
	socket = connect(8052)
	write(socket,unity_mesh)
	close(socket)
end

function test()
	#vertices,faces = vertices_and_faces(mayo_mesh)
	vertices,_faces = vertices_and_faces(bone_mesh)
	pca, wc = pca_weighted(vertices,_faces)
	
	centered_vertices = map(x->x-wc,vertices)
	faces = Vector{GeometryTypes.Face{3,Int32}}()
	for face in _faces
		if all(map(x->x[3],centered_vertices[face]).>30)
			push!(faces,face)
		end
	end
	#centered_mesh = UnityMesh("Mayo Center",centered_vertices,Int32[],Int32[],mayo_triangles,mvayo_color,mayo_options)
	centered_mesh = UnityMesh("Mayo Center",centered_vertices,Int32[],Int32[],bone_triangles,bone_color,["surface_shader = transparent"])
	plot(centered_mesh)
	
	coords = [Point3f0(1,0,0),Point3f0(0,1,0),Point3f0(0,0,1)]
	plot("Coordinate System",coords,Point3f0(0,0,0))
	
	
	obbtree = OBBTree(5,centered_vertices,faces)
	intialize!(obbtree)
	adaptive_refine!(obbtree,1f0)
	plot(obbtree)
	return obbtree
end


function plot(um::UnityMesh)
	socket = connect(8052)
	write(socket,um)
	close(socket)
end

function findfirststop{T<:Function}(f::T,vec::Vector)
	for (i,v) in enumerate(vec)
	(f(v) == true) ? (return i) : (nothing)
	end
	return -1
end


function d(mesh)
	vertices,boundary_faces = vertices_and_faces(mesh)
	vertinds = convert(Vector{Vector{Int32}},boundary_faces)
	unq_vert_inds = unique(vcat(vertinds...))
	boundary_vertices = vertices[unq_vert_inds]
	for i = 1:length(boundary_faces)
		boundary_faces[i] = map(y->findfirststop(x->x==boundary_faces[i][y],unq_vert_inds),1:3)
	end
	return boundary_vertices,boundary_faces
end


function simpletest2()
	file1 = "/Users/maximilianbittens/.julia/v0.6/Unity/models/gmsh/box.msh"
	vertices,eldat = AnsysReader.read(file1)
	triangles = vcat(convert(Vector{Vector{Int32}},eldat)...)-Int32(1)
	options = ["surface_shader = transparent"]
	color = [RGBA{Float32}(1,0,0,1) for i = 1:length(vertices)]
	unity_mesh = UnityMesh("Mesh",vertices,Int32[],Int32[],triangles,color,options)
	plot(unity_mesh)
end

function simpletest()
	file1 = "/Users/maximilianbittens/.julia/v0.6/Unity/models/gmsh/box.msh"
	box_mesh = MeshReaders.load(MeshReaders.GmshReader, file1)
	unity_triangles = vcat(map(x->x[4:end],filter(x->x[1:3]==[2,3,1],box_mesh.elements))-1...)
	vertices = Point3f0[Point3f0(box_mesh.nodes[1,i],box_mesh.nodes[2,i],box_mesh.nodes[3,i]) for i = 1:size(box_mesh.nodes,2)]
	options = ["surface_shader = transparent","point_size = 0.02"]
	color = [RGBA{Float32}(1,1,1,1) for i = 1:length(vertices)]
	unity_mesh = UnityMesh("Mesh",vertices,collect(UInt32(0):UInt32(length(vertices)-1)),Int32[],unity_triangles,color,options)
	plot(unity_mesh)
	vertices,faces = vertices_and_faces(box_mesh)
	obbtree = OBBTree(5,vertices,faces)
	intialize!(obbtree)
	elem = first(obbtree.tree[1])
	plot("elem coords",elem.obb.pca, elem.wc)
	pca,wc = pca_weighted(vertices,faces)
	plot("elem coords2",pca, wc)
	#child1,child2 = refine!(obbtree,elem)
	plot(obbtree)
end

# load mesh
file1 = "/Users/maximilianbittens/.julia/v0.6/Unity/models/gmsh/Box.msh"
box_mesh = MeshReaders.load(MeshReaders.GmshReader, file1)
# get boundary
boundary_vertices,boundary_faces = d(box_mesh)
# get all vertices and faces
vertices,faces = vertices_and_faces(box_mesh)
# get weighted center
pca,wc = pca_weighted(vertices,faces)
#shift vertices and boundary_vertices
#vertices = map(x->x-wc,vertices)
#boundary_vertices = map(x->x-wc,boundary_vertices)

# plot mesh
unity_triangles = vcat(convert(Vector{Vector{Int32}},faces)...)-Int32(1)
options = ["surface_shader = transparent","point_size = 0.02"]
color = [RGBA{Float32}(.5,.5,.5,.3) for i = 1:length(vertices)]
unity_mesh = UnityMesh("Mesh",vertices,collect(UInt32(0):UInt32(length(vertices)-1)),Int32[],unity_triangles,color,options)
plot(unity_mesh)

# plot coordinate system
coords = [Point3f0(1,0,0),Point3f0(0,1,0),Point3f0(0,0,1)]
plot("Coordinate System",coords,Point3f0(0,0,0))

#pcac,wcc = pca_weighted(boundary_vertices,boundary_faces)
#plot("elem coordsc",pcac, wcc+Point3f0(.1,0,0))

obbtree = OBBTree(3,vertices,faces)
intialize!(obbtree)
elem = first(obbtree.tree[1])
#plot("elem coords",elem.obb.pca, elem.wc)
adaptive_refine!(obbtree,2.5f0)
plot(obbtree) 
elem = first(filter(!isempty,obbtree.tree)[end])
#plot_parents(obbtree,elem)

#=
pca,wc = pca_weighted(boundary_vertices, boundary_faces[child1.finds])
obbnew = oriented_bounding_box(boundary_vertices[child1.vinds],pca)
obbtreeelem = OBBTreeElement(1,obbnew,Nullable{OBBTreeElement}(),wc,child1.vinds,child1.finds)

println("OBB Mesh new")
plot("OBB meshnew", obbtree, obbtreeelem, rand(RGBA{Float32}))
plot("Bounding Box new", obbtreeelem.obb, rand(RGBA{Float32}))
plot("elem coordsc",pca, wc)
#simpletest()
=#



