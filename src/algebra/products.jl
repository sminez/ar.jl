#=
All operations on AR types.

Multiplying αs
==============
This is based on a set of simplification rules based on allowed
manipulations of elements in the algebra.
(NOTE:: In all notation, αμ.αν is simplified to αμν)

(1)   αpμ == αμp == αμ
      "Multiplication by αp (r-point) is idempotent. (αp is the identity)"
(2i)  α0^2 == αp
      "Repeated α0 indices can just be removed."
(2ii) αi^2 == -αp
      "Repeated αi indices can be removed by negating"
(2iii) α^2 == +-αp
      "All elements square to either +αp or -αp"
(3)   αμν == -ανμ
      "Adjacent indices can be popped by negating."


I am converting the current product into an array of integers in order
to allow for the different orderings of each final product in a flexible
way. Ordering is a mapping of index (0,1,2,3) to position in the final
product. This should be stable regardless of how we define the 16
elements of the algebra.

The algorithm makes use of the fact that for any ordering we can
dermine the whether the total number of pops is odd or even by looking
at the first element alone and then recursing on the rest of the
ordering as a sub-problem.
If the final position of the first element is even then it will take an
odd number of pops to correctly position it. We can then look only at
the remaining elements and re-label them with indices 1->(n-1) and
repeat the process until we are done.
NOTE:: I proved this by brute force. (i.e. listing all options and
showing that the proposition holds...!)
=#
const DIVISION_TYPE = "into"  # One of "by" or "into"


##########################
# .: Operations on αs :. #
##########################
"""__find_prod__
Pre-compute the product of two αs in the algebra by first checking for special
case such as multiplication by αp and squaring and then using the set of rules
defined above to pop and eliminate indices in order to determine the final
product.
_NOTE_:: The implementation of this is based on the paramaters at the top of this
       file (algebra.jl). These can be modified in order to change the algebra
       and see how the resulting equations are affected.
"""
function find_prod(i::α, j::α; metric=METRIC, allowed=ALLOWED)
    # set the paramaters being used
    # TODO:: Once the paramaters of the algebra have been finalised this should
    #        be moved back to the top of the file.
    metric = Dict(zip(prod([string(m) for m in 0:length(metric)-1]), metric))
    targets = Dict([(Set(a), a) for a in allowed])

    # Rule (1) :: Multiplication by αp is idempotent
    i.index == "p" && return α(j.index, (i.sign * j.sign))
    j.index == "p" && return α(i.index, (i.sign * j.sign))

    # Rule (2) :: Squaring and popping
    sign = i.sign * j.sign
    components = i.index * j.index
    intersection = intersect(Set(i.index), Set(j.index))

    for repeated in intersection
        first, second = find([c for c in components] .== repeated)
        # Distance - 1 as we only need to get adjacent to the first occurance
        n_pops = second - first - 1
        # Only a total odd number of pops will negate
        if n_pops % 2 == 1
            sign *= -1
        end
        # Cancelling unit elements negates based on the metric being used
        sign *= metric[repeated]
        components = String(filter(μ -> μ != repeated, [c for c in components]))
    end

    # If everything cancelled then i == j and we are left with αp (r-point)
    length(components) == 0 && return α("p", sign)

    # Rule (3) :: Popping to the correct order
    target = targets[Set(components)]

    # Allow for immediate return if the product is already in the correct order
    target == components && return α(target, sign)

    ordering = Dict([(c,i) for (i,c) in enumerate(target)])
    current = [ordering[c] for c in components]

    while length(current) > 1
        sign *= iseven(current[1]) ? -1 : 1
        shift!(current)
        new_order = Dict([(j,i) for (i,j) in enumerate(sort(current))])
        current = [new_order[k] for k in current]
    end

    #= Alternate version:: not working yet
    original = [ordering[c] for c in components]

    while length(original) > 0
        k = original[1]
        if iseven(k); sign *= -1 end
        shift!(original)
        for comp in original
            if comp > k; comp -= 1 end
        end
    end
    =#
    return α(target, sign)
end


############################################################
# .: The pre-computed Cayley table and helper functions :. #
############################################################
const CAYLEY = permutedims(
    [find_prod(α(i), α(j)) for i in ALLOWED, j in ALLOWED],
    [1,2]
)

"""Helper to quickly find an α in the Cayley table"""
ix(a::α) = findin(ALLOWED, [a.index])[1]

"""Find the multiplicative inverse of an α element"""
inv(a::α) = α(a.index, (a^2).sign * a.sign)

################################
# .: Operations on αs alone :. #
################################
*(i::α, j::α) = find_prod(i, j)
\(i::α, j::α) = inv(i) * j
/(i::α, j::α) = i * inv(j)

#= NOTE:: This version should be faster as it is a memoisation of find_prod
          but I am getting some very odd behaviour when i = αp that I can't
          work out so I am simply running everything through find_prod for now.
function *(i::α, j::α)
    prod = CAYLEY[ix(i),ix(j)]
    prod.sign *= i.sign * j.sign
    return prod
end
=#


##########################
# .: Operations on ξs :. #
##########################
function *(i::ξ, j::ξ)
    val = Symbol(string(i.val) * string(j.val))
    return ξ(val, j.unit, j.partials)
end


#################################
# .: Operations on ξαs and αs:. #
#################################
*(i::ξα, a::α) = ξα(i.alpha * a, i.xi)
*(a::α, i::ξα) = ξα(a * i.alpha, i.xi)
/(i::ξα, a::α) = ξα(i.alpha / a, i.xi)
\(a::α, i::ξα) = ξα(a \ i.alpha, i.xi)


###########################
# .: Operations on ξαs :. #
###########################
*(i::ξα, j::ξα) = ξα(i.alpha * j.alpha, i.xi * j.xi)


###################################
# .: Operations on Vector{ξα}s :. #
###################################
"""The wedge product of two multi-vectors (symbol is `wedge`)"""
function ∧(v1::Vector{ξα}, v2::Vector{ξα})
    vec = [p[1] * p[2] for p in Iterators.product(v1, v2)]
    return vcat([g for g in groupby(x -> x.alpha.index, vec)]...)
end

"""The dot product of two multi-vectors (symbol is `circ`)"""
function ∘(v1::Vector{ξα}, v2::Vector{ξα})
    product = 0.0
    sort!(v1, lt=(x,y) -> ix(x.alpha) < ix(y.alpha))
    sort!(v2, lt=(x,y) -> ix(x.alpha) < ix(y.alpha))
    shortest = length(v1) < length(v2) ? v1 : v2

    for component in shortest
        a = component.alpha
    end
    # iterate over the shorter and calculate the scalar product of ξs for
    # matching αs
    return ξα(α("p"), product)
end

wedge(v1::Vector{ξα}, v2::Vector{ξα}) = ∧(v1, v2)
dot(v1::Vector{ξα}, v2::Vector{ξα}) = ∘(v1, v2)
full(v1::Vector{ξα}, v2::Vector{ξα}) = vcat([∘(v1, v2)], ∧(v1, v2))

# NOTE:: Also defined as the bare product of two multivectors
*(v1::Vector{ξα}, v2::Vector{ξα}) = vcat([∘(v1, v2)], ∧(v1, v2))


################################################################################
# Elsewhere in the code, you should use div() so that changing between the two
# is managed by changing the value of DIVISION_TYPE in this file alone.
if DIVISION_TYPE == "by"
    div(comp::ξα, a::α) = comp / a
    div(i::ξα, j::ξα) = i / j
    div(i::α, j::α) = i / j
elseif DIVISION_TYPE == "into"
    div(comp::ξα, a::α) = a \ comp
    div(i::ξα, j::ξα) = j \ i
    div(i::α, j::α) = j \ i
else
    throw(ArgumentError("Invalid division type: $DIVISION_TYPE"))
end
