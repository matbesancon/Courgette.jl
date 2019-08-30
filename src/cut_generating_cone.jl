
"""
    cut_generating_cone(optimizer, A, xfrac, b, π, π0, lower_bounds = false)

Given the primal problem:
```
max c ⋅ x
A x ≥ b
```
and some integrality constraints, a branching rule can be deduced,
of the form `π ⋅ x ≤ π0 ⋁ π ⋅ x ≥ π0 + 1`.
Build the Cut Generating Linear Problem (CGLP) and returns the problem,
dual multipliers and cutting coefficients in a named tuple.
If `lower_bounds`, dual multipliers are added for `x ≥ 0`.
"""
function cut_generating_cone(optimizer, A, xfrac, b, π, π0, lower_bounds = false)
    cglp = Model(optimizer)
    n = length(xfrac)
    size(A, 2) == n == length(π) || throw(DimensionError("A and x must be of fitting dimensions"))
    m = size(A, 1)
    @variables(cglp, begin
        u[1:m] >= 0
        u0 >= 0
        v[1:m] >= 0
        v0 >= 0
        γ[1:n]
        γ0
    end)

    @constraints(cglp, begin
        γ .== A' * u .- u0 .* π
        γ .== A' * v .+ v0 .* π
        γ0 == dot(u, b) - u0 * π0
        γ0 == dot(v, b) + v0 * (π0 + 1)
    end)
    @objective(cglp, Min, dot(γ, xfrac) - γ0)
    return (cglp, (γ = γ, γ0 = γ0, u = u, u0 = u0, v = v, v0 = v0))
end
