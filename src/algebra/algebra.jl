#=
This package defines the types and operations required for the
modified Dirac-Clifford Algebra to be computed numerically
according to the principle of Absolute Relativity.

As in the convention used in the theory, greek indices run 0,1,2,3
and latin indices run 1,2,3.

The definition of the 4-differntial Dμ is given in differentail.jl
as this must be computed numerically and the way in which that is
performed is allowed to vary.


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
import Base.*, Base./, Base.\, Base.==, Base.show


####################
# .: Parameters :. #
####################
const ALLOWED = ["p","0","1","2","3","10","20","30","23",
                 "31","12","023","031","012","123","0123"]
const TARGETS = Dict([(Set(a), a) for a in ALLOWED])
const METRIC = Dict(zip("0123", [1 -1 -1 -1]))

#######################
# .: Unit Elements :. #
#######################
"""One of the 16 unit element of the algebra listed in `allowed`"""
type α
    index::String
    sign::Int8

    function α(index::String, sign::Integer)
        sign in [1, -1]  || error("invalid α: $index, $sign")
        index in ALLOWED || error("invalid α: $index, $sign")
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

"""Pretty print αs when working in the REPL"""
show(io::IO, a::alpha) = print(io, (a.sign > 0 ? "" : "-") * "α" * a.index)

"""αs are equal if their indices and signs are equal"""
==(i::alpha, j::alpha) = i.index == j.index && i.sign == j.sign

################################
# .: Vectors and components :. #
################################
typealias AR_pair Tuple{Any,alpha}

typealias ξα Tuple{Real,alpha}
typealias xi_component ξα

typealias Ξ Vector{AR_pair}
typealias xi Ξ



##########################
# .: Operations on αs :. #
##########################
"""Pre-compute the product of two αs in the algebra by first checking for special
case such as multiplication by αp and squaring and then using the set of rules
defined above to pop and eliminate indices in order to determine the final
product.
NOTE:: The implementation of this is based on the paramaters at the top of this
       file (algebra.jl). These can be modified in order to change the algebra
       and see how the resulting equations are affected."""
function find_prod(i::alpha, j::alpha)
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


#################################################################
# .: The pre-computed Cayley table and associated operations :. #
#################################################################
# .' is the transposition operator as multi-arrays in Julia are column major...
const CAYLEY = permutedims(
    [find_prod(α(i), α(j)) for j in allowed, i in allowed],
    [2,1]
)

"""Helper to quickly find an α in the Cayley table"""
ix(a::α) = find(μ -> μ == a.index, allowed)

"""Lookup in the Cayley table and determine the new sign"""
function *(i::α, j::α)
    prod = CAYLEY[ix(i),ix(j)][1]
    prod.sign *= i.sign * j.sign
    return prod
end

"""Find the inverse of the first argument and then multiply"""
function /(i::α, j::α)
    i_inverse = α(i.index, (i * i).sign * i.sign)
    return i_inverse * j
end

"""Find the inverse of the second argument and then multiply"""
function \(i::α, j::α)
    j_inverse = α(j.index, (j * j).sign * j.sign)
    return i * j_inverse
end

################################
# .: Viewing a Cayley Table :. #
################################
"""This is a quick and dirty way to print out the Cayley table in
    a couple of different ways so that it can be pasted into excel
    for visualising and looking at properties such as sign and symmetry."""
function view_cayley(output="indices")
    s = join(CAYLEY[1,:], ",")
    println(",$s")
    for i in 1:16
        c = CAYLEY[i,:]
        print("$(c[1]),")
        if output == "indices"
            println(join([a for a in c], ","))
        elseif output == "strindices"
            println(join([a.index for a in c], ","))
        elseif output == "colmap"
            println(join([find(x -> x == a.index, allowed)[1] for a in c], ","))
        elseif output == "sign"
            println(join([a.sign for a in c], ","))
        else
            error("Invalid output specification")
        end
    end
end
