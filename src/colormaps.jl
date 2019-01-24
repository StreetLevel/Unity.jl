import FixedPointNumbers: N0f8

#defined in Plotting/palettes.jl
colormapdict = Dict{String,Vector{RGB{Float32}}}(

	 "Spectral"  =>  Spectral,
	 "RdYlGn"    =>  RdYlGn,
	 "GnYlRd"    =>  reverse(RdYlGn),
	 "RdBu"      =>  RdBu,
	 "PiYG"      =>  PiYG,
	 "PRGn"      =>  PRGn,
	 "RdYlBu"    =>  RdYlBu,
	 "BuYlRd"    =>  BuYlRd,#reverse(RdYlBu), 
	 "BrBG"      =>  BrBG,
	 "RdGy"      =>  RdGy,
	 "GyRd"      =>  reverse(RdGy),
	 "PuOr"      =>  PuOr,
	 "Arnold"    =>  reverse(PiYG),
	 "BuRd"      =>  reverse(RdBu),
	 "GnYlRd"    =>  reverse(RdYlGn),
	 "Paired"    =>  Paired,
	 "Purples"   =>  Purples,
	 "BlPu"		 =>  BlPu,
	 "BkOr"   =>  BkOr,
	 "GrBl" => GrBl,
	 "distgclrs" => distgclrs,
	 "YG" => YG)


abstract AbstractColorMap

type ColorMap <: AbstractColorMap
	 map::Vector{RGB{Float32}}
	 function ColorMap(str::String)
		  clrmap = new()
		  clrmap.map = generate(colormapdict[str],1:256)
		  return clrmap
	 end
end

type BoundedColorMap
	 clrmap::ColorMap
	 mindta::Nullable{Float32}
	 maxdta::Nullable{Float32}
	 function BoundedColorMap(clrmap::ColorMap)
		  return new(clrmap,Nullable{Float32}(),Nullable{Float32}())
	 end
end

type test
	 x::Float32
	 y::Nullable{Float32}
end

function generate(colormap::Vector{RGB{Float32}},data::AbstractVector)
	 inds = map(x->Int(floor(x)),linspace(1,length(data),length(colormap)))
	 clrs = map(mean,map(x->colormap[find(y->y==x,inds)],unique(inds)))
	 inds[1] += sum(map(i->inds[i+1]-inds[i],1:(length(inds)-1)))-length(data)
	 return length(data) == length(clrs) ? clrs : vcat(map(i->(inds[i+1]-inds[i]) > 1 ? linspace(clrs[i],clrs[i+1],inds[i+1]-inds[i]) : clrs[i] ,1:(length(inds)-1))...)
end

import Base:mean
function mean(clrs::Vector{RGB{Float32}})
	 r,g,b,a = zeros(Float32,4)
	 n = length(clrs)
	 for clr in clrs
		  r+=clr.r
		  g+=clr.g
		  b+=clr.b
	 end
	 r/=n;g/=n;b/=n
	 return RGB{Float32}(r,g,b)
end

type JetColorMap <: AbstractColorMap
	 map::Vector{RGBA{Float32}}
	 function JetColorMap()
		  cmap = Array(RGBA{Float32}, 256)
		  for i = 0:255
				n=4*i/256
				cmap[i+1] = RGBA(min(max(min(n-1.5,-n+4.5),0.),1.),
									 min(max(min(n-0.5,-n+3.5),0.),1.),
									 min(max(min(n+0.5,-n+2.5),0.),1.),1.)
		  end
		  return new(cmap)
	 end
end

type HSVColorMap <: AbstractColorMap
	 map::Vector{RGBA{Float32}}
	 function HSVColorMap()
		  return new(linspace(Colors.HSV(0,1,1),Colors.HSV(330,1,1),256))
	 end
end

function apply!(cm::BoundedColorMap, data::Union{Array{Float32,1},Array{Float64,1},AbstractVector{Int}})
	 cm.mindta = minimum(data)
	 cm.maxdta = maximum(data)
	 return apply(cm.clrmap,data)
end

function apply(cm::AbstractColorMap,data::Union{Array{Float64,1},Array{Float32,1}})
	 color = Array{RGBA{Float32}}(length(data))
	 mindta = minimum(data)
	 maxdta = maximum(data)
	 interv = maxdta-mindta
	 if isapprox(abs(interv),0.,atol=1e-10)
		  fill!(color,cm.map[floor(Int,length(cm.map)/2)])
	 else
		  pos = 0.0
		  for i = 1:length(data)
				pos = (data[i]-mindta)/interv * 255.0 + 1.0
				clr_1 = cm.map[floor(Int64,pos)]
				clr_2 = cm.map[ceil(Int64,pos)]
				weight = pos%1
				clr = weighted_color_mean(weight,clr_2,clr_1)
				#clr = (1-weight)*clr_1+weight*clr_2
				color[i] = clr
		  end
	 end
	 return color
end

function apply(cm::AbstractColorMap,data::AbstractVector{Int})
	 color = Array{RGBA{Float32}}(length(data))
	 mindta = minimum(data)
	 maxdta = maximum(data)
	 interv = maxdta-mindta
	 if interv < 1
		  fill!(color,cm.map[floor(Int,length(cm.map)/2)])
	 else
		  pos = 0.0
		  for i = 1:length(data)
				pos = (data[i]-mindta)/interv * 255.0 + 1.0
				clr_1 = cm.map[floor(Int64,pos)]
				clr_2 = cm.map[ceil(Int64,pos)]
				weight = pos%1
				clr = weighted_color_mean(weight,clr_2,clr_1)
				#clr = (1-weight)*clr_1+weight*clr_2
				color[i] = clr
		  end
	 end
	 return color
end
