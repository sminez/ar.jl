#=
A selection of different implementations for numerically computing
the 4-vector 4-differentail Dμ.

In Cartesian coordinates this is:
    Dμ = ∂μ/αμ = ∂ / (αμ∂xμ) = α0∂0 - αi∂i = α0∂0 - ∇i
=#
const ξ_GROUPS = [
    ("Ξi", ["1", "2", "3"]),
    ("Ξi0", [a for a in ALLOWED if length(a) == 2 && contains(a, "0")]),
    ("Ξjk", [a for a in ALLOWED if length(a) == 2 && ~contains(a, "0")]),
    ("Ξ0jk", [a for a in ALLOWED if length(a) == 3 && contains(a, "0")])
]
const α_TO_GROUP = Dict(
    vcat([[(μ, group) for μ in αs] for (group, αs) in ξ_GROUPS]...)
)
const ∂_RIGHT = Dict(1=>"2", 2=>"3", 3=>"1")
const ∂_LEFT = Dict(1=>"3", 2=>"1", 3=>"2")


############################################
# .: Derivative operations on Ξ vectors :. #
############################################
"""Differnetiation of a symbolic_ξα just records the operation"""
function diff_ξ(wrt::α, j::symbolic_ξ)
    new_ξ = symbolic_ξ(j.val, j.unit, [p for p in j.partials])
    push!(new_ξ.partials, wrt)
    return new_ξ
end

function ∂(component::ξα, ∂_wrt::α)
    # Using div so that we can chose division by/into in algebra.jl
    new_α = div(component.alpha, ∂_wrt)
    new_ξ = diff_ξ(∂_wrt, component.xi)
    return symbolic_ξα(new_α, new_ξ)
end

∂(vec::Vector, wrt::α) = [∂(comp, wrt) for comp in vec]
∇i(vec::Vector) = vcat([∂(vec, α(i)) for i in ["1" "2" "3"]]...)
Dμ(vec::Vector) = vcat(∂(vec, α("0")), ∇i(vec))
DG(vec::Vector) = vcat([∂(vec, α(i)) for i in ALLOWED]...)

# Non-unicode versions
partial(comp::ξα, wrt::α) = ∂(comp, wrt)
partial(vec::Vector, wrt::α) = ∂(vec, wrt)
del_i(vec::Vector) = ∇i(vec)
D_mu(vec::Vector) = Dμ(vec)


##############################################
# .: Operations on differentiated vectors :. #
##############################################
"""__by_α__ :: Vector{symbolic_ξα} -> Iterators.GroupBy

Group derivative terms according to their unit element.
If vector_groups=true then this will group the elements into
"""
function by_α(vec::Vector{symbolic_ξα}, vector_groups=false)
    sort!(vec, lt=(x,y) -> ix(x.alpha) < ix(y.alpha))
    if vector_groups
        #XXX :: default to grouping by index if not part of a group
        return groupby(x -> get(α_TO_GROUP, x.alpha.index, x.alpha.index), vec)
    else
        return groupby(x -> x.alpha.index, vec)
    end
end

function by_ξ(vec::Vector{symbolic_ξα})
    sort!(vec, lt=(x,y) -> ix(α(x.xi.unit)) < ix(α(y.xi.unit)))
    return groupby(x -> x.xi.unit, vec)
end


"""__by_∇__ :: Vector{symbolic_ξα} -> Vector{symbolic_ξα}

Attempt to identify vector derivatives in each group
TODO:: Allow for grouping of second derivatives using the level flag"""
function by_∇(vec::Vector{symbolic_ξα}, level=1)
    output = Vector{symbolic_ξα}()

    for components in by_α(vec, true)
        for (group_name, group_ξs) in ξ_GROUPS
            # Define helpers
            function div_check(j::symbolic_ξα)
                ∂_index = parse(j.xi.partials[level].index)
                return ∂_index == 0 ? false : group_ξs[∂_index] == j.xi.unit
            end

            function curl_check(j::symbolic_ξα)
                ∂_index = parse(j.xi.partials[level].index)

            end

            # Check for ∇•group_name
            matches = filter(div_check, components)
            if matches != []
                first_α = matches[1].alpha
                if all([m.alpha == first_α for m in matches])
                    # XXX:: This is either an α index or the name of a Ξ group
                    lead_unit = matches[1].xi.unit
                    div_component = symbolic_ξ("∇•$group_name", lead_unit)
                    push!(output, symbolic_ξα(first_α, div_component))
                    filter!(x -> !(x in matches), components)
                end
            end

            # Check for ∇xgroup_name
            matches = filter(curl_check, components)
        end
        # XXX:: add remaining elements back in
        append!(output, components)
    end
    return output
end


"""__show_by_α__
Pretty print an α_grouped derivative.
"""
function show_by_α(vec::Vector{symbolic_ξα}, vector_notation=false)
    grouped  = vector_notation ? by_α(vec, true) : by_α(vec)
    for group in grouped
        print("α$(group[1].alpha.index)(")
        formatted = [
            "$(e.alpha.sign == 1 ? "+" : "-")$(e.xi)"
            for e in group
        ]
        print(join(formatted, " "))
        println(")")
    end
end
