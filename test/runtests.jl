using Courgette
using Test

import Random
import GLPK
import LinearAlgebra
using JuMP
using SparseArrays: sprand

@testset "CGLP cuts fractional point" begin
    Random.seed!(33)
    A = round.(Int, 10 .* sprand(100, 50, 0.5))
    b = round.(Int, rand(100) .* 20)
    Ah = vcat(A, -Int.(LinearAlgebra.Diagonal(LinearAlgebra.I, 50)), Int.(LinearAlgebra.Diagonal(LinearAlgebra.I, 50)))
    bh = vcat(b, -ones(Int, 50), zeros(Int, 50))
    primal_model = Model(
        with_optimizer(GLPK.Optimizer, msg_lev = GLPK.OFF)
    )
    @variable(primal_model, 0 <= x[1:50] <= 1)
    @constraint(primal_model, A * x .>= b)
    @objective(primal_model, Min, sum(x))
    optimize!(primal_model)
    @test termination_status(primal_model) == MOI.OPTIMAL
    xvalue = JuMP.value.(x)
    jfrac = 1
    while jfrac <= length(xvalue)
        if xvalue[jfrac] > 5 * eps(Float64) && xvalue[jfrac] < 1 - 5 * eps(Float64)
            break
        end
        jfrac += 1
    end
    π = eachindex(xvalue) .== jfrac
    π0 = 0

    (cglp, variables) = Courgette.cut_generating_cone(with_optimizer(GLPK.Optimizer, msg_lev = GLPK.OFF), Ah, xvalue, bh, π, π0)
    Courgette.normalize_cglp!(Courgette.StandardNormalization(), cglp::JuMP.Model, variables)
    optimize!(cglp)
    @test termination_status(cglp) == MOI.OPTIMAL
    @test LinearAlgebra.dot(JuMP.value.(variables[:γ]), xvalue) - JuMP.value(variables[:γ0]) < 0
end
