module Unity
import JSON
import ColorTypes
import ColorTypes: RGB, RGBA
import Colors: weighted_color_mean
import GeometryTypes
import GeometryTypes: SimpleRectangle, Point3f0, Face
import Sockets.TCPSocket

include(joinpath("palette.jl"))
include(joinpath("colormaps.jl"))

mutable struct UnityVector3
    x::Float32
    y::Float32
    z::Float32
end

mutable struct UnityColor
    r::Float32
    g::Float32
    b::Float32
    a::Float32
end

mutable struct UnityText
    text::String
    pos::UnityVector3
    scale::UnityVector3
    rot::UnityVector3
    color::UnityColor
end

function Base.convert(::Type{UnityColor},c::ColorTypes.RGBA{Float32})
    return UnityColor(c.r,c.g,c.b,c.alpha)
end
function Base.convert(::Type{UnityVector3},x::GeometryTypes.Point3f0)
    return UnityVector3(x[1],x[2],x[3])
end
#function Base.convert(::Vector{UInt32},x::Array{GeometryTypes.Face{3,Int32},1})
#    vcat(map(xx->convert(Vector{UInt32},xx-1),x)...)
#end

#Unity mesh with c-like indexing
mutable struct UnityMesh
    id::String
    vertices::Vector{UnityVector3}
    points::Vector{UInt32}
    lines::Vector{UInt32}
    triangles::Vector{UInt32}
    colors::Vector{UnityColor}
    options::Vector{String}
    text::Vector{UnityText}
    visible::Vector{Bool}
end

function UnityMesh(id::String, vertices::Vector, points::Vector,  lines::Vector, triangles::Vector, colors::Vector, options::Vector, visible=Bool[false])
    return UnityMesh(id,vertices,points,lines,triangles,colors,options,UnityText[],visible)
end


function Base.write(socket::TCPSocket, um::UnityMesh)
    jum = JSON.json(um)
    retval = write(socket, jum*"UNITY_MESH_JSON_FORMATTED")
    return retval
end

function clear(socket::TCPSocket)
    retval = write(socket, "UNITY_RESET_ALL")
    return retval
end

function clear(socket::TCPSocket, id::Int)
    @assert 0 <= id <= 9
    retval = write(socket, "UNITY_RESET_$id")
    return retval
end

function screenshot(socket::TCPSocket,filename::String)
    retval = write(socket, filename*"UNITY_SCREENSHOT")
    return retval
end

mutable struct UnityCameraSettings
    id::String
    main_camera_position::Vector{UnityVector3} #should be all vectors of length 0 or 1
    main_camera_rotation::Vector{UnityVector3}
    main_camera_scale::Vector{UnityVector3}
    background_color::Vector{UnityColor}
    perspective::Vector{Bool}
end

function Base.write(socket::TCPSocket, ucs::UnityCameraSettings)
    jum = JSON.json(ucs)
    retval = write(socket, jum*"UNITY_CAMERA_SETTINGS")
    return retval
end

#Unity Pyramid mesh with c-like indices
mutable struct PyramidMesh
    id::String
    vertices::Vector{UnityVector3}
    pyramids::Vector{UInt32}
    colors::Vector{ColorTypes.RGBA{Float32}}
end

begin
local pattern = [ 0,2,1,0,3,2,2,3,1,0,1,3 ]
function Base.convert(::Type{UnityMesh},msh::PyramidMesh,dublic_vert::Bool=false)
    if dublic_vert
        return convert_and_duplicate(UnityMesh,msh,pattern)
    else
        return convert(UnityMesh,msh,pattern)
    end
end
end

function Base.convert(::Type{UnityMesh},msh::PyramidMesh,pattern::Vector{Int})
    triangles = Vector{UInt32}()
    inds = msh.pyramids
    @assert mod(length(inds),4) == 0
    for i = 1:4:length(inds)
        append!(triangles,
            [
            inds[i+pattern[1]],inds[i+pattern[2]],inds[i+pattern[3]],
            inds[i+pattern[4]],inds[i+pattern[5]],inds[i+pattern[6]],
            inds[i+pattern[7]],inds[i+pattern[8]],inds[i+pattern[9]],
            inds[i+pattern[10]],inds[i+pattern[11]],inds[i+pattern[12]]
            ]
            )
    end
    return UnityMesh(msh.id,msh.vertices,UInt32[],UInt32[],triangles,msh.colors)
end

function convert_and_duplicate(::Type{UnityMesh},msh::PyramidMesh,pattern::Vector{Int})
    id = msh.id
    verts = msh.vertices
    lines = UInt32[]
    points = UInt32[]
    inds = msh.pyramids
    clrs = msh.colors
    @assert mod(length(inds),4) == 0

    #new mesh
    #triangles = Vector{Int32}()
    vertices = Vector{UnityVector3}()
    colors = Vector{ColorTypes.RGBA{Float32}}()

    for i = 1:4:length(inds)

        tmpinds = inds[i:i+3]+1
        tmpverts = verts[tmpinds]
        tmpclrs = clrs[tmpinds]

        append!(vertices, tmpverts[pattern+1])
        append!(colors, tmpclrs[pattern+1])

    end

    triangles = UInt32[i-1 for i = 1:length(vertices)]

    return UnityMesh(id,vertices,points,lines,triangles,colors)

end

function Base.write(tcpstream::TCPSocket, um::PyramidMesh)
    jum = JSON.json(convert(UnityMesh,um))
    return write(tcpstream, jum)
end

import Combinatorics

function boundary(upm::PyramidMesh)
    verts = Vector{GeometryTypes.Point3f0}()
    pyramids = Vector{UInt32}()
    clrs = Vector{ColorTypes.RGBA{Float32}}()
    for i = 1:4:length(upm.pyramids)
        # todo
    end
end


function ColorBar(clrmap::BoundedColorMap)

    add_verts_faces_color!(vertices,faces,lines,colors,a,b,c,d,ca,cb,flip) = begin
        push!(vertices,a) #0
        push!(vertices,b) #1
        push!(vertices,c) #2
        push!(vertices,d) #3
        push!(colors,ca)
        push!(colors,cb)
        push!(colors,ca)
        push!(colors,cb)
        #flip ? push!(faces,Face{3,UInt32}(length(vertices)-4,length(vertices)-1,length(vertices)-2)) : push!(faces,Face{3,UInt32}(length(vertices)-1,length(vertices)-4,length(vertices)-2))
        #flip ? push!(faces,Face{3,UInt32}(length(vertices)-4,length(vertices)-3,length(vertices)-1)) : push!(faces,Face{3,UInt32}(length(vertices)-3,length(vertices)-4,length(vertices)-1))
        #push!(lines,length(vertices)-4)
        #push!(lines,length(vertices)-2)
        #push!(lines,length(vertices)-3)
        #push!(lines,length(vertices)-1)
        flip ? push!(faces,Face{3,UInt32}(length(vertices)-3,length(vertices),length(vertices)-1)) : push!(faces,Face{3,UInt32}(length(vertices),length(vertices)-3,length(vertices)-1))
        flip ? push!(faces,Face{3,UInt32}(length(vertices)-3,length(vertices)-2,length(vertices))) : push!(faces,Face{3,UInt32}(length(vertices)-2,length(vertices)-3,length(vertices)))
        push!(lines,length(vertices)-3)
        push!(lines,length(vertices)-1)
        push!(lines,length(vertices)-2)
        push!(lines,length(vertices))
    end

    offset = 0.
    n = 10

    rel = 10.
    sr = SimpleRectangle(0.,0.,rel,.5)

    a = map(Point3f0,convert(Vector{Float64},linspace(0,sr.w,n)),repmat([offset+sr.h],n),repmat([0],n))
    b = map(Point3f0,convert(Vector{Float64},linspace(0,sr.w,n)),repmat([offset],n),repmat([0],n))

    vals = linspace(clrmap.mindta.value, clrmap.maxdta.value, n )
    clrmap = apply(clrmap.clrmap, convert( Vector{Float32}, vals ) )


    faces = Face{3,UInt32}[]
    lines = UInt32[]
    vertices = Point3f0[]
    colors = Vector{RGBA{Float32}}()

    for i = 1:n-1
      add_verts_faces_color!(vertices,faces,lines,colors,a[i],a[i+1],b[i],b[i+1],clrmap[i],clrmap[i+1],false)
    end

    return vertices,faces,lines,colors,vals,b
end
import Printf.@sprintf
function ColorBar(Id::String, clrmap::BoundedColorMap, alpha=0.4)

    vertices,faces,lines,colors,vals,b = ColorBar(clrmap)
    unity_vertices = map(x->UnityVector3(x[1]-8000,x[2]-8000,x[3]-105), vertices )
    unity_faces = vcat([[UInt32(x[1]),UInt32(x[2]),UInt32(x[3])] for x in faces]...)
    unity_colors = map(x->RGBA{Float32}(x.r,x.g,x.b,alpha),colors)
    text = map((x,y)->UnityText(@sprintf("%.3e",x),y-GeometryTypes.Point3f0(7999.6,8000.1,105),UnityVector3(.1,.1,.05),UnityVector3(1,180,1),RGBA{Float32}(0.0,0.0,0.0,1.0)), vals,b)
    #plot(ub,Id,unity_vertices,unity_faces,lines,unity_colors,ub.config[:wireframe_options],text)
    return UnityMesh(Id,unity_vertices,Int32[],Int32[],unity_faces,unity_colors,["surface_shader = transparent"],text)

end



end #module Unity
