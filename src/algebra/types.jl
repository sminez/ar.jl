#=
The following types are used for all calculations within the package:
    α
    ξ
    Ξ
=#

####################
# .: Parameters :. #
####################
# NOTE:: Original ordering from the paper
#const ALLOWED = ["p","0","1","2","3","10","20","30","23",
#                 "31","12","023","031","012","123","0123"]
const ALLOWED = [
    "p","23","31","12",
    "0","023","031","012",
    "123","1","2","3",
    "0123","10","20","30",
]
# NOTE:: This is the "natural" order for the table using my rule of:
#        Gn+1 = |Gn   Gn^n |
#               |Gn^n Gn   |
#        To build up the table...the sign distribution is...interesting!
#=
const ALLOWED = [
    "p","0","1","01",
    "2","02","12","012",
    "3","03","31","031",
    "23","023","123","0123"
]
=#
const ALLOWED_GROUPS = [Symbol(g) for g in ["p","0","i","i0","jk","0jk","123","0123"]]
const METRIC = [1 -1 -1 -1]


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
        new(val, sign)
    end

    function α(group::Symbol, sign::Integer)
        sign in [1, -1]  || error("invalid α: $index, $sign")
        #group in ALLOWED_GROUPS || error("invalid α: $group, $sign")
        new(string(group), sign)
    end

    function α(a::α)
        error("αs can not be initialised with another α")
    end
end

# Non-unicode version
typealias alpha α

==(i::α, j::α) = (i.index == j.index) && (i.sign == j.sign)
-(a::α) = α(a.index, (-1 * a.sign))


"""A simple representation of a ξ component to manipulate symbolically"""
type ξ
    val::Symbol
    unit::String
    partials::Vector{α}

    # TODO:: provide some validation on s
    ξ(s::String, g::String) = new(Symbol(s), g, Vector{α}())
    ξ(s::String) = new(Symbol(s), s, Vector{α}())
    ξ(a::α) = new(Symbol("ξ" * a.index), a.index, Vector{α}())
    ξ(a::α, p::Vector{α}) = new(Symbol("ξ" * a.index), a.index, p)
    ξ(s::Symbol, u::String, p::Vector{α}) = new(s, u, p)
end

==(i::ξ, j::ξ) = (i.val == j.val) && (i.partials == j.partials)


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
abstract Abstract_ξα

immutable ξα <: Abstract_ξα
    alpha::α
    xi::ξ
    wrt::Set{α}

    """Symbols track the initial component and differentiate wrt everything"""
    function ξα(alpha::α)
        xi = ξ(alpha)
        wrt = Set([α(i) for i in ["0","1","2","3"]])
        new(alpha, xi, wrt)
    end

    function ξα(alpha::α, xi::ξ)
        new(alpha, xi, Set([α(i) for i in ["0","1","2","3"]]))
    end
end


==(i::ξα, j::ξα) = (i.alpha == j.alpha) && (i.xi == j.xi) && (i.wrt == j.wrt)


type function_ξα <: Abstract_ξα
    alpha::α
    xi::Expr
    wrt::Set{α}
end

# Non-unicode
typealias AR_pair ξα


###############################
# .: Formatting & printing :. #
###############################
function show(io::IO, a::α)
    print(io, "$(a.sign > 0 ? "" : "-")α$(a.index)")
end

function show(io::IO, j::ξ)
    print(io, "$(join(["∂"*p.index for p in reverse(j.partials)]))$(j.val)")
end

function show(io::IO, j::ξα)
    print(io, "$(j.alpha) $(j.xi)")
end


###################################################
# .: Initialisors for the general multivectors :. #
###################################################
Ξp   = [ξα(α("p"))]
Ξμ   = [ξα(α(a)) for a in ["0", "1", "2", "3"]]
Ξμν  = [ξα(α(a)) for a in ALLOWED if length(a) == 2]
Ξμνρ = [ξα(α(a)) for a in ALLOWED if length(a) == 3]
ΞM = vcat([ξα(α("p"))], [ξα(α(a)) for a in ALLOWED if length(a) == 2 && !('0' in a)])
ΞT = vcat([ξα(α("0"))], [ξα(α(a)) for a in ALLOWED if length(a) == 3 && '0' in a])
ΞA = vcat([ξα(α("123"))], [ξα(α(a)) for a in ALLOWED if length(a) == 1 && !('0' in a)])
ΞE = vcat([ξα(α("0123"))], [ξα(α(a)) for a in ALLOWED if length(a) == 2 && '0' in a])
ΞG = [ξα(α(a)) for a in ALLOWED]

# Non-unicode versions
Xip = Ξp
Xi1 = Ξμ
Xi2 = Ξμν
Xi3 = Ξμνρ

XiM = ΞM
XiT = ΞT
XiA = ΞA
XiE = ΞE
XiG = ΞG
