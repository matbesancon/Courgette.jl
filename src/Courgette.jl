module Courgette

using SparseArrays
using JuMP
using LinearAlgebra: dot, I, Diagonal
import Random

include("cut_generating_cone.jl")
include("normalizing.jl")

end # module
