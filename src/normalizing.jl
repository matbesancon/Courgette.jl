
abstract type ConeNormalization end

struct TrivialNormalization <: ConeNormalization end
struct StandardNormalization <: ConeNormalization end

"""
    normalize_cglp!(::ConeNormalization, cglp, variables)

Add a constraint or set of constraints to the CGLP using the
`ConeNormalization` technique. Return the `ConstraintRef`
"""
function normalize_cglp! end

function normalize_cglp!(::TrivialNormalization, cglp::JuMP.Model, variables::NamedTuple)
    @constraint(cglp, variables.u0 + variables.v0 == 1)
end

function normalize_cglp!(::StandardNormalization, cglp::JuMP.Model, variables::NamedTuple)
    @constraint(cglp, sum(variables.u + variables.v) + variables.u0 + variables.v0 == 1)
end
