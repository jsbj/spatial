using DataStructures
using PyPlot
using PyCall
using StatsBase
using Distributions
using Glob
using JLD
using Shell
using LocallyWeightedRegression
using LinearAlgebra
xr = pyimport("xarray")

include("helpers.jl")
include("config.jl")
include("transforms.jl")
include("creators.jl")
include("figures.jl")

# cd("/project2/abbot/jsbj/spatial")
# cd("/Users/jonah/repos/papers/active/spatial_from_control")
const models = create_model_dict()
const grid = "12x24_equaldist"

# if length(ARGS) > 0
#   args_list = split(ARGS[1],":")
#   str = args_list[1] * "("
#
#   str_list = []
#
#
#   for i in 2:length(args_list)
#     if match(r"\b\d+\b",args_list[i]) != nothing
#       push!(str_list,args_list[i])
#     else
#       if args_list[i] in ["false","true"]
#         push!(str_list,args_list[i])
#       else
#         if i in [2,3,5,length(args_list)]
#           push!(str_list,"\"" * args_list[i] * "\"")
#         else
#           push!(str_list,":$(Symbol(args_list[i]))")
#         end
#       end
#     end
#   end
#
#   str = str * join(str_list,",") * ")"
#   println(str)
#   eval(parse(str))
# end
