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
        allowed = [
            "p","0","1","2","3","10","20","30",
            "23","31","12","023","031","012",
            "123","0123"
        ]
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
"Check that unit_elements are correctly ordered"
function <=(i::unit_element, j::unit_element)
    d = Dict("1" => "2", "2" => "3", "3" => "1")
    if i.index == j.index
        return true
    elseif i.index == "0"
        return true
    elseif j.index == "0"
        return false
    elseif d[i.index] == j.index
        return true
    else
        return false
    end
end

function ==(i::unit_element, j::unit_element)
    i.index == j.index
end

"Multiplication of single index αs"
function *(i::unit_element, j::unit_element)
    if i.index == "0"
        return α("-" * j.index * "0")
    elseif j.index == "0"
        return α(i.index * "0")
    else
        # NOTE:: > is overloaded to give the correct ordering here
        if i < j
            return α(i.index * j.index)
        else
            return α("-" * j.index * i.index)
        end
    end
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
        "Adjacent indices can be poppped by negating."

    In addition to this, I am using a simple lookup table to simplify the
    squaring of elements as all square to either αp or -αp.
    =#
    print(i)
    println(j)

    # Rule (1) :: Multiplication by r-point
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

    # Rule (2) :: Squaring
    if i.index == j.index
        if i.index in ["0", "10", "20", "30", "123"]
            return α("p")
        else
            return α("-p")
        end
    end

    components = append!([c for c in i.index], [c for c in j.index])
    intersection = intersect(Set(i.index), Set(j.index))

    # Determine if the initial product is +ve or -ve
    if ('-' in components) && ~('-' in intersection)
        # Either i or j is negative
        negate = true
    else
        # Both i and j are negative or neither is
        negate = false
    end

    # Rule (3) :: Sorting of elements via bubble sort to track pops
    units = [unit_element(c) for c in components]
    needed_pops = true

    while needed_pops
        needed_pops = false
        for (i, unit) in enumerate(units[1:length(units)-1])
            if unit <= units[i+1]
                # Cancel repeated indices
                if unit == units[i+1]
                    filter!(μ -> μ != unit, units)
                    if unit.index > 0
                        negate = ~negate
                    end
                    needed_pops = true
                end
                break
            else
                # Bubble sort into the correct place and track pops
                while ~(units[i] <= units[i+1])
                    units[i:i+1] = units[i+1:-1:i]
                    needed_pops = true
                    negate = ~negate
                    i = i + 1
                    if i == length(units)
                        if units[end] == unit_element("1") && length(units) > 2
                            # NOTE:: This breaks the cycle of 123123...
                            units[i-1:i] = units[i:-1:i-1]
                        end
                        break
                    end
                end
            end
        end
    end

    product = join([u.index for u in units])

    # Deal with αi0 elements
    if length(product) == 2 && product[1] == '0'
        product = reverse(product)
        negate = ~negate
    end
    # Deal with negation
    if negate
        return α("-" * product)
    else
        return α(product)
    end
end


# .: Generate Caley Table :. #
allowed = [
    "p","0","1","2","3","10","20","30",
    "23","31","12","023","031","012",
    "123","0123"
]

caley = [[α(i) * α(j) for j in allowed] for i in allowed]
print(caley)
