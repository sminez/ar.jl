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
"""__by_α__

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


"""__to_del__

Replace the elements found by one of the replacement functions with an
alternate form using ∇ notation.
The replacement functions all use this as final step to do the replacement."""
function to_del(elements::Vector{symbolic_ξα}, replacement::String, sign=0)
    first_α = elements[1].alpha
    α_group = get(α_TO_GROUP, first_α.index, first_α.index)
    sign = (sign != 0) ? sign : first_α.sign
    element_α = α(Symbol(α_group), sign)
    new_component = symbolic_ξ(replacement)
    return symbolic_ξα(element_α, new_component)
end


function replace_grad(terms, group_name, groups_ξs, level)
    output = Vector{symbolic_ξα}()

    function grad_check(j::symbolic_ξα, term::symbolic_ξα, wrt::String)
        term_match = j.xi.unit == term.xi.unit
        ∂_match = j.xi.partials[level].index == wrt
        return term_match && ∂_match
    end

    for term in terms
        grad_elements = Vector{symbolic_ξα}()
        for wrt in ["1","2","3"]
            elem = filter(j->grad_check(j, term, wrt), terms)
            elem != [] && push!(grad_elements, elem[1])
        end

        if length(grad_elements) == 3
            push!(output, to_del(grad_elements, "∇Ξ$group_name"))
            filter!(x -> !(x in grad_elements), terms)
        end
    end
    return (terms, output)
end


function replace_div(terms, group_name, group_ξs, level)
    output = Vector{symbolic_ξα}()

    function div_check(j::symbolic_ξα)
        ∂_index = parse(j.xi.partials[level].index)
        return ∂_index == 0 ? false : group_ξs[∂_index] == j.xi.unit
    end

    div_elements = filter(div_check, terms)

    if length(div_elements) == 3
        if all([d.alpha == div_elements[1].alpha for d in div_elements])
            push!(output, to_del(div_elements, "∇Ξ•$group_name"))
            filter!(x -> !(x in div_elements), terms)
        end
    end
    return (terms, output)
end


function replace_curl(terms, group_name, group_ξs, level)
    output = nothing
    output = Vector{symbolic_ξα}()

    function curl_check(j::symbolic_ξα, term::Symbol, wrt::α)
        diff = j.xi.partials[level]
        return (j.xi.val == term) && (diff == wrt)
    end

    curl_elements = Vector{symbolic_ξα}()
    signs = Vector{Integer}()
    curl_element_missing = false

    for (i, index) in enumerate(group_ξs)
        ∂_cw = α(string(CW[i]))
        ∂_acw = α(string(ACW[i]))
        ξ_cw = Symbol("ξ" * group_ξs[CW[i]])
        ξ_acw = Symbol("ξ" * group_ξs[ACW[i]])

        positive_term = filter(j->curl_check(j, ξ_acw, ∂_cw), terms)
        negative_term = filter(j->curl_check(j, ξ_cw, ∂_acw), terms)

        if positive_term != [] && negative_term != []
            # NOTE:: filter returns an array so we need to extract the value
            #        only when we know that the filter returned a value
            pos, neg = positive_term[1].alpha.sign, negative_term[1].alpha.sign
            # If both terms have the same sign this isn't a Curl
            curl_element_missing = (pos == neg)

            if !curl_element_missing
                append!(curl_elements, [positive_term[1], negative_term[1]])
                push!(signs, pos)
            end
        else
            # We are missing one or both terms
            curl_element_missing = true
        end
        # Stop looking for terms if anything is missing
        curl_element_missing && break
    end

    if length(curl_elements) == 6
        sign = all(signs .== 1) ? 1 : 0
        sign = all(signs .== -1) ? -1 : sign
        push!(output, to_del(curl_elements, "∇Ξx$group_name", sign))
        filter!(x -> !(x in curl_elements), terms)
    end
    return (terms, output)
end


function replace_group_partials(terms, group_name, group_ξs, level)
    output = Vector{symbolic_ξα}()

    function partial_check(j::symbolic_ξα, term::String, wrt::String)
        j.xi.unit == term && j.xi.partials[level].index == wrt
    end

    for wrt in ["0","1","2","3"]
        group_elements = Vector{symbolic_ξα}()
        for term in group_ξs
            elem = filter(j->partial_check(j, term, wrt), terms)
            if elem != []
                push!(group_elements, elem[1])
            end
        end

        if length(group_elements) == 3
            push!(output, to_del(group_elements, "∂$(wrt)Ξ$group_name"))
            filter!(x -> !(x in group_elements), terms)
        end
    end
    return (terms, output)
end


"""__by_∇__

Attempt to identify vector derivatives in each group
TODO:: Allow for grouping of second derivatives using the level flag"""
function by_∇(vec::Vector{symbolic_ξα}, level=1)
    output = Vector{symbolic_ξα}()

    for terms in by_α(vec, true)
        for (group, ξs) in ξ_GROUPS
            (terms, replacements) = replace_grad(terms, group, ξs, level)
            replacements != [] && append!(output, replacements)

            (terms, replacements) = replace_div(terms, group, ξs, level)
            replacements != [] && append!(output, replacements)

            (terms, replacements) = replace_curl(terms, group, ξs, level)
            replacements != [] && append!(output, replacements)

            (terms, replacements) = replace_group_partials(terms, group, ξs, level)
            replacements != [] && append!(output, replacements)
        end

        # Add remaining elements under the correct grouped α value
        for term in terms
            i = term.alpha.index
            term.alpha.index = get(α_TO_GROUP, i, i)
        end
        append!(output, terms)
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
        grouped = groupby(x->get(α_TO_GROUP, x.alpha.index, x.alpha.index), vec)
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
