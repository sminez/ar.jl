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
const ALLOWED_GROUPS = [Symbol(g) for g in ["p","0","i","i0","jk","0jk","123","0123"]]
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
        sign in [1, -1]  || throw(TypeError("invalid α: $index, $sign"))
        index in ALLOWED || throw(TypeError("invalid α: $index, $sign"))
        new(index, Int8(sign))
    end

    function α(index::String)
        sign = '-' in index ? -1 : 1
        val = sign > 0 ? index : index[2:end]
        val in ALLOWED || throw(TypeError("invalid α: $index"))
        new(val, sign)
    end

    function α(group::Symbol, sign::Integer)
        sign in [1, -1]  || error("invalid α: $index, $sign")
        group in ALLOWED_GROUPS || error("invalid α: $group, $sign")
        new(string(group), sign)
    end

    function α(a::α)
        throw(TypeError("αs can not be initialised with another α"))
    end
end

# Non-unicode version
typealias alpha α

==(i::α, j::α) = (i.index == j.index) && (i.sign == j.sign)
-(a::α) = α(a.index, (-1 * a.sign))


"""A simple representation of a ξ component to manipulate symbolically"""
type symbolic_ξ
    val::Symbol
    unit::String
    partials::Vector{α}

    # TODO:: provide some validation on s
    symbolic_ξ(s::String, g::String) = new(Symbol(s), g, Vector{α}())
    symbolic_ξ(s::String) = new(Symbol(s), s, Vector{α}())
    symbolic_ξ(a::α) = new(Symbol("ξ" * a.index), a.index, Vector{α}())
    symbolic_ξ(a::α, p::Vector{α}) = new(Symbol("ξ" * a.index), a.index, p)
    symbolic_ξ(s::Symbol, u::String, p::Vector{α}) = new(s, u, p)
end

==(i::symbolic_ξ, j::symbolic_ξ) = (i.val == j.val) && (i.partials == j.partials)


################################
# .: Vectors and components :. #
################################
"""__ξα vector components__
(component, alpha) pairs are stored as ξα data structures that also keep track
of which partial derivaties we need to caclulate.
This _must_ be specified directly in the case of a floating point array representing
the component values of Ξ at each point in spacetime but it will be inferred for
functions and symbolic values.

The intended way of creating function_ξα values is via the f-string macro defined
in `parse.jl` in the utils directory. Helper functions for initialising the
symbolic representations of the general multi-vectors are provided below.
"""
abstract ξα

immutable symbolic_ξα <: ξα
    alpha::α
    xi::symbolic_ξ
    wrt::Set{α}

    """Symbols track the initial component and differentiate wrt everything"""
    function symbolic_ξα(alpha::α)
        xi = symbolic_ξ(alpha)
        wrt = Set([α(i) for i in ["0","1","2","3"]])
        new(alpha, xi, wrt)
    end

    function symbolic_ξα(alpha::α, xi::symbolic_ξ)
        new(alpha, xi, Set([α(i) for i in ["0","1","2","3"]]))
    end
end

==(i::ξα, j::ξα) = (i.alpha == j.alpha)&&(i.xi == j.xi)&&(i.wrt == j.wrt)


type function_ξα <: ξα
    alpha::α
    xi::Expr
    wrt::Set{α}
end

type array_ξα <: ξα
    alpha::α
    xi::Array{Float64}
    wrt::Set{α}
end

# Non-unicode
typealias AR_pair ξα
typealias symbolic_AR_pair symbolic_ξα
typealias function_AR_pair function_ξα
typealias array_AR_pair array_ξα


###############################
# .: Formatting & printing :. #
###############################
function show(io::IO, a::α)
    print(io, "$(a.sign > 0 ? "" : "-")α$(a.index)")
end

function show(io::IO, j::symbolic_ξ)
    print(io, "$(join(["∂"*p.index for p in reverse(j.partials)]))$(j.val)")
    #print(io, "$(join(["∂"*p.index for p in j.partials]))ξ$(j.initial.index)")
end

function show(io::IO, j::ξα)
    print(io, "$(j.alpha) $(j.xi)")
end


###################################################
# .: Initialisors for the general multivectors :. #
###################################################
Ξp   = [symbolic_ξα(α("p"))]
Ξμ   = [symbolic_ξα(α(a)) for a in ["0", "1", "2", "3"]]
Ξμν  = [symbolic_ξα(α(a)) for a in ALLOWED if length(a) == 2]
Ξμνρ = [symbolic_ξα(α(a)) for a in ALLOWED if length(a) == 3]
ΞG   = [symbolic_ξα(α(a)) for a in ALLOWED]

# Non-unicode versions
xi_p = Ξp
xi_1 = Ξμ
xi_2 = Ξμν
xi_3 = Ξμνρ
xi_G = ΞG

"""Validator for custom Ξ definitions"""
function check_Ξ(vec::Vector)
    length(Set(v.alpha for v in vec)) < length(vec) && error("Repeated α in Ξ")
    return vec
end

check_xi(vec::Vector) = check_Ξ(vec)


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
function find_prod(i::α, j::α)
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
        sign *= (n_pops % 2 == 1 ? -1 : 1)
        # Cancelling unit elements negates based on the metric being used
        sign *= METRIC[repeated]
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
    [find_prod(α(i), α(j)) for i in ALLOWED, j in ALLOWED],
    [1,2]
)

"""Helper to quickly find an α in the Cayley table"""
ix(a::α) = findin(ALLOWED, [a.index])[1]

"""Helper to aid in transferring sign information from αs to ξs"""
extract_sign(a::α) = (a.sign == -1) ? (α(a.index, 1), -1) : (a, 1)

"""Find the multiplicative inverse of an α element"""
inv(a::α) = α(a.index, (a^2).sign * a.sign)

################################
# .: Operations on αs alone :. #
################################
#= NOTE:: This version should be faster as it is a memoisation of find_prod
          but I am getting some very odd behaviour when i = αp that I can't
          work out so I am simply running everything through find_prod for now.

function *(i::α, j::α)
    prod = CAYLEY[ix(i),ix(j)]
    prod.sign *= i.sign * j.sign
    return prod
end
=#
"""Find the product of two α values"""
*(i::α, j::α) = find_prod(i, j)

"""Find the inverse of the first argument and then multiply"""
\(i::α, j::α) = inv(i) * j

"""Find the inverse of the second argument and then multiply"""
/(i::α, j::α) = i * inv(j)


###########################
# .: Operations on ξαs :. #
###########################
#=

TODO:: Get this working with the new ξα definition

NOTE:: Addition/subtraction of ξα pairs is only defined for matching αs.
It is also deliberate that there is no definition given for multiplication
of an ξα pair by a scalar: in accordance with the principle of Absolute
Relativity, this is not allowed and instead we must multiply by a
(scalar, αp) pair.

function +(i::ξα, j::ξα)
    (iξ, iα), (jξ, jα) = (i.xi, i.alpha), (j.xi, j.alpha)
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
=#


#################################
# .: Operations on ξαs and αs:. #
#################################
function *(i::ξα, a::α)
    new_ξα = copy(i)
    new_ξα.alpha = new_ξα.alpha * a
    return new_ξα
end

function *(a::α, i::ξα)
    new_ξα = copy(i)
    new_ξα.alpha = a * new_ξα.alpha
    return new_ξα
end

function \(a::α, i::ξα)
    new_ξα = copy(i)
    new_ξα.alpha = a \ new_ξα.alpha
    return new_ξα
end

function /(i::ξα, a::α)
    new_ξα = copy(i)
    new_ξα.alpha = new_ξα.alpha / a
    return new_ξα
end

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
