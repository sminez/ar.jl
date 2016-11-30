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
import Base.*, Base.==, Base.show


####################
# .: Paramaters :. #
####################
const allowed = ["p","0","1","2","3","10","20","30","23",
                 "31","12","023","031","012","123","0123"]
const targets = Dict([(Set(a), a) for a in allowed])
const αp_squares =  ["0", "10", "20", "30", "123"]

#######################
# .: Unit Elements :. #
#######################
abstract alpha

"One of the 16 unit element of the algebra"
type α <: alpha
    index::String
    sign::Int8

    function α(index::String, sign::Integer)
        sign in [1, -1]  || error("invalid α: $index, $sign")
        index in allowed || error("invalid α: $index, $sign")
        index in allowed || error("invalid α: $index, $sign")
        new(index, Int8(sign))
    end

    function α(index::String)
        sign = '-' in index ? -1 : 1
        val = sign > 0 ? index : index[2:end]
        val in allowed || error("invalid α: $index")
        α(index, sign)
    end
end

"Pretty print αs when working in the REPL"
function show(io::IO, a::alpha)
    str = (a.sign > 0 ? "" : "-") * "α" * a.index
    print(io, str)
end

################################
# .: Vectors and components :. #
################################
typealias AR_pair Tuple{Any,α}
typealias ξμ Tuple{Real,α}

typealias Ξ Vector{AR_pair}
typealias xi Ξ


#########################
# .: Functions on αs :. #
#########################
function ==(i::alpha, j::alpha)
    i.index == j.index && i.sign == j.sign
end

"""Compute the product of two αs in the algebra by first checking for special
case such as multiplication by αp and squaring and then using the set of rules
defined above to pop and eliminate indices in order to determine the final
product.
NOTE:: The implementation of this is based on the paramaters at the top of this
       file (algebra.jl). These can be modified in order to change the algebra
       and see how the resulting equations are affected.
       The choice of division by or division into is made in differntial.jl
       as this is implemented as a flag on the division implementations there"""
function *(i::alpha, j::alpha)
    # Rule (1) :: Multiplication by αp
    if contains(i.index, "p")
        return α(j.index, (i.sign * j.sign))
    elseif contains(j.index, "p")
        return α(i.index, (i.sign * j.sign))
    end

    # Rule (2) :: Squaring and popping
    if i.index == j.index
        sign = (i.index in αp_squares) ? 1 : -1
        return α("p", sign * i.sign * j.sign)
    end

    components = i.index * j.index
    intersection = intersect(Set(i.index), Set(j.index))

    # Determine if the initial sign is +ve or -ve
    sign = i.sign * j.sign

    # Remove duplicate indices and handle negation from associated pops
    for repeated in intersection
        first, second = find([c for c in components] .== repeated)
        # Distance - 1 as we only need to get adjacent to the first occurance
        n_pops = second - first - 1
        if n_pops % 2 == 1
            # Only a total odd number of pops will negate
            sign *= -1
        end
        if repeated != '0'
            # Cancelling α0 is idempotend, cancelling αi negates
            sign *= -1
        end
        filtered = filter(μ -> μ != repeated, [c for c in components])
        components = String(filtered)
    end

    # Rule (3) :: Sorting of elements via bubble sort to track pops
    target = targets[Set(components)]

    # Allow for immediate return if the product is already in the correct order
    target == components && return α(target, sign)

    #=  NOTE:: I am converting the current product into an array of integers
        in order to allow for the different orderings of each final product in
        a flexible way. Ordering is a mapping of index (0,1,2,3) to position in
        the final product. This should be stable regardless of how we define the
        16 elements of the algebra.
    =#
    ordering = Dict([(c,i) for (i,c) in enumerate(target)])
    current = [ordering[c] for c in components]
    _, n_pops = count_pops(current)
    if n_pops % 2 == 1
        sign *= -1
    end

    return α(target, sign)
end


"""Count the number of pops required to sort a source array into ascending
order. This is really just a modified version of the standard recursive
merge sort with a count of the number of pops are used during the merge step.
In the main implementation of * we don't care about the final sorted array but
we need to return it in order to allow the recursive implementation."""
function count_pops(source::Array{Int, 1})
    # Base case to terminate recursion
    if length(source) == 1
        return (source, 0)
    end

    # Split the input into two halves
    mid = div(length(source), 2)
    left, right = source[1:mid], source[mid+1:end]

    # Recurse of each half to get sub counts
    sorted_left, left_pops = count_pops(left)
    sorted_right, right_pops = count_pops(right)

    # Merge the halves together and count pops
    total_pops = left_pops + right_pops
    i, j = 0, 0
    merged = []

    while (i < length(sorted_left)) && (j < length(sorted_right))
        # Taking from the left if both are equal
        if sorted_left[i+1] <= sorted_right[j+1]
            append!(merged, sorted_left[i+1])
            i += 1
        else
            append!(merged, sorted_right[j+1])
            j += 1
            total_pops += (length(sorted_left) - i)
        end
    end

    # Append any remaing element from the array that wasn't drained
    i == length(sorted_left) && append!(merged, sorted_right[j+1:end])
    j == length(sorted_right) && append!(merged, sorted_left[i+1:end])

    return merged, total_pops
end

#################################
# .: Generate a Cayley Table :. #
#################################
#=  This is a quick and dirty way to print out the Cayley table in
    a couple of different ways so that it can be pasted into excel
    for visualising and looking at properties such as sign and symmetry.
=#
cayley = [[α(i) * α(j) for j in allowed] for i in allowed]
s = join(cayley[1], ",")
println(",$s")
for c in cayley
    print("$(c[1]),")
    # Indices only for quick colour scale view
    println(join([find(x -> x == a.index, allowed)[1] for a in c], ","))
    # String indices
    #println(join([a.index for a in c], ","))
    # Sign
    #println(join([a.sign for a in c], ","))
end
