#=
This package defines the types and operations required for the
modified Dirac-Clifford Algebra to be computed numerically
according to the principle of Absolute Relativity.

As in the convention used in the theory, greek indices run 0,1,2,3
and latin indices run 1,2,3.

The definition of the 4-differntial Dμ is given in differentail.jl
as this must be computed numerically and the way in which that is
performed is allowed to vary.
=#
import Base.*, Base.<=, Base.==, Base.show

####################
# .: Paramaters :. #
####################
const allowed = ["p","0","1","2","3","10","20","30", "23",
                 "31","12","023","031","012", "123","0123"]
const targets = Dict([(Set(a), a) for a in allowed])
const unit_ordering = Dict("1" => "2", "2" => "3", "3" => "1")
const αp_squares =  ["0", "10", "20", "30", "123"]

#######################
# .: Unit Elements :. #
#######################
abstract alpha

"Single element Indices for computing and ordering products"
type unit_element <: alpha
    index::String

    function unit_element(i::Char)
        if i ∈ ['0','1','2','3']
            new([i])
        else
            error("α must be ∈ [0,1,2,3]")
        end
    end

    function unit_element(i::String)
        if contains("0123", i) && length(i) == 1
            new(i)
        else
            error("α must be ∈ [0,1,2,3]")
        end
    end
end

"One of the 16 unit element of the algebra"
type α <: alpha
    index::String

    "Confirm that index is actually valid"
    function α(index::String)
        if index[1] == '-'
            val = index[2:end]
        else
            val = index
        end

        if val in allowed
            new(index)
        else
            error("invalid α: $index")
        end
    end
end

"Pretty print αs when working in the REPL"
function show(io::IO, a::alpha)
    if a.index[1] == '-'
        print(io, "-α" * a.index[2:end])
    else
        print(io, "α" * a.index)
    end
end

#########################
# .: Functions on αs :. #
#########################
function <=(i::unit_element, j::unit_element)
    if i.index == j.index
        return true
    elseif i.index == "0"
        return true
    elseif j.index == "0"
        return false
    # See the Paramaters section at the top of the file
    elseif unit_ordering[i.index] == j.index
        return true
    else
        return false
    end
end

function ==(i::alpha, j::alpha)
    i.index == j.index
end

"Compute the product of two αs in the algebra"
function *(i::alpha, j::alpha)
    #=
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
    # Rule (1) :: Multiplication by αp
    if contains(i.index, "p")
        if contains(i.index, "-")
            return α("-" * j.index)
        else
            return j
        end
    elseif contains(j.index, "p")
        if contains(j.index, "-")
            return α("-" * i.index)
        else
            return i
        end
    end

    # Rule (2) :: Squaring and popping
    if i.index == j.index
        if i.index in αp_squares
            return α("p")
        else
            return α("-p")
        end
    end

    #components = append!([c for c in i.index], [c for c in j.index])
    components = i.index * j.index
    intersection = intersect(Set(i.index), Set(j.index))

    # Determine if the initial product is +ve or -ve
    if ('-' in components) && ~('-' in intersection)
        # One of i or j is negative
        negate = true
    else
        # Both i and j are negative or neither is
        negate = false
    end

    # Remove duplicate indices and handle negation from associated pops
    for repeated in intersection
        # NOTE:: .<operator> is array broadcast syntax
        first, second = find([c for c in components] .== repeated)
        # Distance - 1 as we only need to get adjacent to the first occurance
        n_pops = second - first - 1
        if n_pops % 2 == 1
            # Only a total odd number of pops will negate
            negate = ~negate
        end
        filtered = filter(μ -> μ != repeated, [c for c in components])
        components = String(filtered)
    end

    # Rule (3) :: Sorting of elements via bubble sort to track pops
    target = targets[Set(components)]
    # Convert the components into an array of ints so we can sort easier
    ordering = Dict([(c,i) for (i,c) in enumerate(target)])
    current = [ordering[c] for c in components]
    _, n_pops = count_pops(current)
    if n_pops % 2 == 1
        negate = ~negate
    end

    # Deal with negation
    if negate
        return α("-" * target)
    else
        return α(target)
    end
end


"Count the number of pops required to sort a source array into ascending order"
function count_pops(source::Array{Int, 1})
    # NOTE:: This is just recursive merge sort but keeping track of the pops
    # Base case
    if length(source) == 1
        return (source, 0)
    end
    # Split
    mid = div(length(source), 2)
    left = source[1:mid]
    right = source[mid+1:end]
    # Recurse
    sorted_left, left_pops = count_pops(left)
    sorted_right, right_pops = count_pops(right)
    # Merge
    i, j = 1, 1
    total_pops = left_pops + right_pops
    merged = []
    while (i < length(sorted_left)) && (j < length(sorted_right))
        # Taking from the left if both are equal
        if sorted_left[i] <= sorted_right[j]
            append!(merged, sorted_left[i])
            i += 1
        else
            append!(merged, sorted_right[j])
            j += 1
            total_pops += (length(sorted_left) - i)
        end
    end
    if i == length(sorted_left)
        append!(merged, sorted_right[j:end])
    elseif j == length(sorted_right)
        append!(merged, sorted_left[i:end])
    end
    return merged, total_pops
end


##############################
# .: Generate Caley Table :. #
##############################
caley = [[α(i) * α(j) for j in allowed] for i in allowed]
for c in caley
    println(c)
end
