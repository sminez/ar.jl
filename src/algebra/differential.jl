#=
A selection of different implementations for numerically computing
the 4-vector 4-differentail Dμ.

In Cartesian coordinates this is:
    Dμ = ∂μ/αμ = ∂ / (αμ∂xμ) = α0∂0 - αi∂i = α0∂0 - ∇i
=#
const ξ_GROUPS = [
    ("i", ["1", "2", "3"]),
    ("i0", [a for a in ALLOWED if length(a) == 2 && contains(a, "0")]),
    ("jk", [a for a in ALLOWED if length(a) == 2 && ~contains(a, "0")]),
    ("0jk", [a for a in ALLOWED if length(a) == 3 && contains(a, "0")])
]
const α_TO_GROUP = Dict(
    vcat([[(μ, group) for μ in αs] for (group, αs) in ξ_GROUPS]...)
)
# These two dictionaries allow for iterating over the definition of curl by
# associating each element in a Ξ group with its clockwise and counter-clockwise
# elements of the 3-cycle formed by the elements: 1 → 2 → 3
#                                                 ↑ ← ← ← ↓
const CW = Dict(1=>2, 2=>3, 3=>1)
const ACW = Dict(1=>3, 2=>1, 3=>2)


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


###################################################
# .: Grouping terms for display of derivatives :. #
###################################################
"""__by_α__ :: Vector{symbolic_ξα} -> Iterators.GroupBy

Group derivative terms according to their unit element.
If vector_groups=true then this will group the elements into"""
function by_α(vec::Vector{symbolic_ξα}, vector_groups=false)
    sort!(vec, lt=(x,y) -> ix(x.alpha) < ix(y.alpha))
    if vector_groups
        return groupby(x -> get(α_TO_GROUP, x.alpha.index, x.alpha.index), vec)
    else
        return groupby(x -> x.alpha.index, vec)
    end
end


function replace_grad(components, group_name, groups_ξs)
    output = Vector{symbolic_ξα}()

    function grad_check(j::symbolic_ξα, comp::symbolic_ξα, wrt::String)
        comp_match = j.xi.unit == comp.xi.unit
        ∂_match = j.xi.partials[level].index == wrt
        return comp_match && ∂_match
    end

    for comp in components
        grad_elements = []
        for wrt in ["1","2","3"]
            elem = filter(j->grad_check(j,comp,wrt), components)
            elem != [] && push!(grad_elements, elem[1])
        end

        if length(grad_elements) == 3
            first_α = grad_elements[1].alpha
            α_group = get(α_TO_GROUP, first_α.index, first_α.index)
            grad_component = symbolic_ξ("∇Ξ$group_name")
            grad_α = α(Symbol(α_group), first_α.sign)
            push!(output, symbolic_ξα(grad_α, grad_component))
            filter!(x -> !(x in grad_elements), components)
        end
        return (components, output)
    end
end


function replace_div(components, group_name, groups_ξs)
    output = Vector{symbolic_ξα}()

    function div_check(j::symbolic_ξα)
        ∂_index = parse(j.xi.partials[level].index)
        return ∂_index == 0 ? false : group_ξs[∂_index] == j.xi.unit
    end

    div_elements = filter(div_check, components)
    if div_elements != []
        first_α = div_elements[1].alpha

        if all([d.alpha == first_α for d in div_elements])
            # XXX:: This is either an α index or the name of a Ξ group
            α_group = get(α_TO_GROUP, first_α.index, first_α.index)
            div_component = symbolic_ξ("∇•Ξ$group_name")
            div_α = α(Symbol(α_group), first_α.sign)
            push!(output, symbolic_ξα(div_α, div_component))
            filter!(x -> !(x in div_elements), components)
        end
    end
    return (components, output)
end


function replace_curl(components, group_name, group_ξs)
    output = Vector{symbolic_ξα}()

    function curl_check(j::symbolic_ξα, term::Symbol, wrt::α)
        diff = j.xi.partials[level]
        return (j.xi.val == term) && (diff == wrt)
    end

    curl_elements = Vector{Tuple{symbolic_ξα,Symbol}}()
    curl_element_missing = false

    for (i, index) in enumerate(group_ξs)
        # curl is [∂cw(ξacw) - ∂acw(ξcw)] ∀ indices in group_ξs
        # NOTE:: This is [pos, neg] below until I work out some better
        #        names for everything!
        ∂_cw = α(string(CW[i]))
        ∂_acw = α(string(ACW[i]))
        ξ_cw = Symbol("ξ" * group_ξs[CW[i]])
        ξ_acw = Symbol("ξ" * group_ξs[ACW[i]])
        pos = filter(j->curl_check(j, ξ_acw, ∂_cw), components)
        neg = filter(j->curl_check(j, ξ_cw, ∂_acw), components)

        if pos != [] && neg != []
            # extract the filtered elements from the arrays
            pos, neg = pos[1], neg[1]

            if pos.alpha.sign == 1 && neg.alpha.sign == -1
                append!(curl_elements, [(pos, :+), (neg, :+)])
            elseif pos.alpha.sign == -1 && neg.alpha.sign == 1
                append!(curl_elements, [(pos, :-), (neg, :-)])
            else
                curl_element_missing = true
            end
        else
            curl_element_missing = true
        end
        curl_element_missing && break
    end

    if length(curl_elements) == 6
        sign = 0
        if all(comp[2] == :+ for comp in curl_elements)
            sign = 1
        elseif all(comp[2] == :- for comp in curl_elements)
            sign = -1
        end

        curl_component = symbolic_ξ("∇xΞ$group_name")
        first_α = curl_elements[1][1].alpha
        α_group = get(α_TO_GROUP, first_α.index, first_α.index)
        curl_α = α(Symbol(α_group), sign)
        push!(output, symbolic_ξα(curl_α, curl_component))
        filter!(x -> !(x in [e[1] for e in curl_elements]), components)
    end
    return (components, output)
end


function replace_group_partials(components, group_name, group_ξs)
    output = Vector{symbolic_ξα}()

    function dmu_check(j::symbolic_ξα, comp::String, wrt::String)
        j.xi.unit == comp && j.xi.partials[level].index == wrt
    end

    for wrt in ["0","1","2","3"]
        xi_elements = []
        for comp in group_ξs
            elem = filter(j->dmu_check(j,comp,wrt), components)
            if elem != []
                push!(xi_elements, elem[1])
            end
        end

        if length(xi_elements) == 3
            first_α = xi_elements[1].alpha
            α_group = get(α_TO_GROUP, first_α.index, first_α.index)
            xi_component = symbolic_ξ("∂$(wrt)Ξ$group_name")
            xi_α = α(Symbol(α_group), first_α.sign)
            push!(output, symbolic_ξα(xi_α, xi_component))
            filter!(x -> !(x in xi_elements), components)
        end
    end
    return (components, output)
end


"""__by_∇__ :: Vector{symbolic_ξα} -> Vector{symbolic_ξα}

Attempt to identify vector derivatives in each group
TODO:: Allow for grouping of second derivatives using the level flag"""
function by_∇(vec::Vector{symbolic_ξα}, level=1)
    output = Vector{symbolic_ξα}()

    for components in by_α(vec, true)
        for (group, ξs) in ξ_GROUPS
            components, grads = replace_grad(components, group, ξs)
            append!(output, grads)

            components, divs = replace_div(components, group, ξs)
            append!(output, divs)

            components, curls = replace_curl(components, group, ξs)
            append!(output, curls)

            components, group_∂s = replace_group_partials(components, group, ξs)
            append!(output, group_∂s)
        end

        # Add remaining elements under the correct grouped α value
        for component in components
            i = component.alpha.index
            component.alpha.index = get(α_TO_GROUP, i, i)
        end
        append!(output, components)
    end
    return output
end


"""__show_by_α__
Pretty print an α_grouped derivative.
"""
function show_by_α(vec::Vector{symbolic_ξα}, vector_notation=false)
    if vector_notation
        function vec_sort(x,y)
            ix_x = findin(ALLOWED_GROUPS, [Symbol(x.alpha.index)])[1]
            ix_y = findin(ALLOWED_GROUPS, [Symbol(y.alpha.index)])[1]
            return ix_x < ix_y
        end
        sort!(vec, lt=vec_sort)
        grouped = groupby(x -> get(α_TO_GROUP, x.alpha.index, x.alpha.index), vec)
    else
        grouped = by_α(vec)
    end

    for group in grouped
        i = group[1].alpha.index
        if vector_notation
            print("α$(get(α_TO_GROUP, i, i))(")
        else
            print("α$i(")
        end
        formatted = [
            "$(e.alpha.sign == 1 ? "+" : "-")$(e.xi)"
            for e in group
        ]
        print(join(formatted, " "))
        println(")")
    end
end
