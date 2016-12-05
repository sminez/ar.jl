#=
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
=#

####################
# .: Parameters :. #
####################
const ALLOWED = ["p","0","1","2","3","10","20","30","23",
                 "31","12","023","031","012","123","0123"]
const TARGETS = Dict([(Set(a), a) for a in ALLOWED])
const METRIC = Dict(zip("0123", [1 -1 -1 -1]))
const DIVISION_TYPE = "into"  # One of "by" or "into"


#######################
# .: Unit Elements :. #
#######################
"""One of the 16 unit element of the algebra listed in ALLOWED"""
type α
    index::String
    sign::Int8

    function α(index::String, sign::Integer)
        sign in [1, -1]  || error("invalid α: $index, $sign")
        index in ALLOWED || error("invalid α: $index, $sign")
        new(index, Int8(sign))
    end

    function α(index::String)
        sign = '-' in index ? -1 : 1
        val = sign > 0 ? index : index[2:end]
        val in ALLOWED || error("invalid α: $index")
        α(index, sign)
    end
end

# Non-unicode version
typealias alpha α

"""αs are equal if their indices and signs are equal"""
==(i::α, j::α) = i.index == j.index && i.sign == j.sign
"""negation swaps the sign"""
-(i::α) = α(i.index, (-1*i.sign))


################################
# .: Vectors and components :. #
################################
#=
    NOTE:: While αs are allowed to have either sign, when they are bound as
    part of an AR_pair their sign information is transferred to the paramater
    they are paired with in order to simplify implementations of calculations
    and improve performance.

    As these operations require manipulation of α values, they are defined at
    the bottom of the file after the generation of the Cayley table.
=#
typealias AR_pair Tuple{Any,α}

typealias ξα Tuple{Float64,α}
typealias xi_component ξα

typealias Ξ Vector{ξα}
typealias xi Ξ


###############################
# .: Formatting & printing :. #
###############################
show(io::IO, a::α) = print(io, "$(a.sign > 0 ? "" : "-")α$(a.index)")
show(io::IO, j::ξα) = print(io, "α$(j[2].index) $(j[1])")


###################################################
# .: Initialisors for the general multivectors :. #
###################################################
Ξμ(init::Float64)   = [(init, α(a)) for a in "0123"]
Ξμν(init::Float64)  = [(init, α(a)) for a in ALLOWED if length(a) == 2]
Ξμνρ(init::Float64) = [(init, α(a)) for a in ALLOWED if length(a) == 3]
ΞG(init::Float64)   = [(init, α(a)) for a in ALLOWED]
# Non-unicode versions
xi_1(init::Float64) = Ξμ(init::Float64)
xi_2(init::Float64) = Ξμν(init::Float64)
xi_3(init::Float64) = Ξμνρ(init::Float64)
xi_G(init::Float64) = ΞG(init::Float64)


##########################
# .: Operations on αs :. #
##########################
"""
__find_prod__

Pre-compute the product of two αs in the algebra by first checking for special
case such as multiplication by αp and squaring and then using the set of rules
defined above to pop and eliminate indices in order to determine the final
product.
_NOTE_:: The implementation of this is based on the paramaters at the top of this
       file (algebra.jl). These can be modified in order to change the algebra
       and see how the resulting equations are affected.
"""
function find_prod(i::α, j::α)
    # Rule (1) :: Multiplication by αp is idempotent
    i.index == "p" && return α(j.index, (i.sign * j.sign))
    j.index == "p" && return α(i.index, (i.sign * j.sign))

    # Rule (2) :: Squaring and popping
    components = i.index * j.index
    intersection = intersect(Set(i.index), Set(j.index))

    # Determine if the sign of the initial product is +ve or -ve
    sign = i.sign * j.sign

    for repeated in intersection
        first, second = find([c for c in components] .== repeated)
        # Distance - 1 as we only need to get adjacent to the first occurance
        n_pops = second - first - 1
        # Only a total odd number of pops will negate
        sign *= (n_pops % 2 == 1 ? -1 : 1)
        # Cancelling unit elements negates based on the metric being used
        sign *= METRIC[repeated]
        # Remove the duplicate elements
        components = String(filter(μ -> μ != repeated, [c for c in components]))
    end

    # If everything cancelled then i == j and we are left with αp (r-point)
    length(components) == 0 && return α("p", sign)

    # Rule (3) :: Sorting of elements via bubble sort to track pops
    target = TARGETS[Set(components)]

    # Allow for immediate return if the product is already in the correct order
    target == components && return α(target, sign)

    #=
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
    ordering = Dict([(c,i) for (i,c) in enumerate(target)])
    current = [ordering[c] for c in components]
    while length(current) > 0
        sign *= iseven(current[1]) ? -1 : 1
        shift!(current)
        new_order = Dict([(j,i) for (i,j) in enumerate(sort(current))])
        current = [new_order[k] for k in current]
    end

    return α(target, sign)
end


############################################################
# .: The pre-computed Cayley table and helper functions :. #
############################################################
const CAYLEY = permutedims(
    [find_prod(α(i), α(j)) for j in ALLOWED, i in ALLOWED],
    [2,1]
)

"""Helper to quickly find an α in the Cayley table"""
ix(a::α) = find(μ -> μ == a.index, ALLOWED)[1]

"""Helper to aid in transferring sign information from αs to ξs"""
extract_sign(a::α) = (a.sign == -1) ? (α(a.index, 1), -1) : (a, 1)


################################
# .: Operations on αs alone :. #
################################
"""Lookup in the Cayley table and determine the new sign"""
function *(i::α, j::α)
    prod = CAYLEY[ix(i),ix(j)]
    prod.sign *= i.sign * j.sign
    return prod
end

"""Find the inverse of the first argument and then multiply"""
function \(i::α, j::α)
    i_inverse = α(i.index, (i * i).sign * i.sign)
    return i_inverse * j
end

"""Find the inverse of the second argument and then multiply"""
function /(i::α, j::α)
    j_inverse = α(j.index, (j * j).sign * j.sign)
    return i * j_inverse
end


################################
# .: Operations on ξα pairs :. #
################################
#=
    NOTE:: Addition/subtraction of ξα pairs is only defined for matching αs.
    It is also deliberate that there is no definition given for multiplication
    of an ξα pair by a scalar: in accordance with the principle of Absolute
    Relativity, this is not allowed and instead we must multiply by a
    (scalar, αp) pair.
=#
function +(i::ξα, j::ξα)
    (iξ, iα), (jξ, jα) = i, j
    iα == jα || error("can only add components with the same α: $iα != $jα")
    return (iξ+jξ, iα)
end

function -(i::ξα, j::ξα)
    jξ, jα = j
    return i + (-jξ, jα)
end

function *(i::ξα, j::ξα)
    (iξ, iα), (jξ, jα) = i, j
    new_α, sign = extract_sign(iα*jα)
    ξ = sign * iξ * jξ
    return (ξ, new_α)
end

function \(i::ξα, j::ξα)
    (iξ, iα), (jξ, jα) = i, j
    new_α, sign = extract_sign(iα \ jα)
    ξ = sign * (iξ \ jξ)
    return (ξ, new_α)
end

function /(i::ξα, j::ξα)
    (iξ, iα), (jξ, jα) = i, j
    new_α, sign = extract_sign(iα / jα)
    ξ = sign * (iξ / jξ)
    return (ξ, new_α)
end


######################################
# .: Operations on ξα pairs and αs:. #
######################################
*(i::ξα, a::α) = (i[1], a * i[2])
*(a::α, i::ξα) = i * a
\(a::α, i::ξα) = (i[1], a \ i[2])
/(i::ξα, a::α) = (i[1], i[2] / a)

# Elsewhere in the code, you should use div() so that changing between the two
# is managed by changing the value of DIVISION_TYPE in this file alone.
if DIVISION_TYPE == "by"
    div(comp::ξα, a::α) = comp / a
    div(i::ξα, j::ξα) = i / j
elseif DIVISION_TYPE == "into"
    div(comp::ξα, a::α) = a \ comp
    div(i::ξα, j::ξα) = j \ i
else
    error("Invalid division type: $DIVISION_TYPE")
end
