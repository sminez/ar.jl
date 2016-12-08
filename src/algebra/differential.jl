#=
A selection of different implementations for numerically computing
the 4-vector 4-differentail Dμ.

In Cartesian coordinates this is:
    Dμ = ∂μ/αμ = ∂ / (αμ∂xμ) = α0∂0 - αi∂i = α0∂0 - ∇i
=#


############################################
# .: Derivative operations on Ξ vectors :. #
############################################
function ∂(component::symbolic_ξα, ∂_wrt::α)
    # Correct the sign for components that square to -αp
    # Using div so that we can chose division by/into in algebra.jl
    a = div(component.alpha, ∂_wrt)
    # Correct sign for components that square to -αp
    a.sign *= (∂_wrt^2).sign
    s = "∂" * ∂_wrt.index * string(component.xi)
    return symbolic_ξα(a, s, component.wrt)
end

∂(vec::Vector, wrt::α) = [∂(comp, wrt) for comp in vec]
∇i(vec::Vector) = vcat([∂(vec, α(i)) for i in ["1" "2" "3"]]...)
Dμ(vec::Vector) = vcat(∂(vec, α("0")), ∇i(vec))

# Non-unicode versions
partial(comp::symbolic_ξα, wrt::α) = ∂(comp, wrt)
partial(vec::Vector, wrt::α) = ∂(vec, wrt)
del_i(vec::Vector) = ∇i(vec)
D_mu(vec::Vector) = Dμ(vec)


###########################################
# .: Operations on derivatived vectors :. #
###########################################
"""
__by_α__

Group derivative terms according to their unit element.
This converts the result of a derivative - (∂,α-orig,ξα-new) - into tuple of
three α values - (α,∂,ξ) - representing the new unit element, the unit element
the derivative was taken with respect to and the origin of the ξ component in
the original Ξ multivector.

If the vector_notation flag is true then this will attempt to identify divs,
grads and curls in the resulting grouped derivative.
"""
function by_α(vec::Vector{symbolic_ξα})
    sort!(vec, lt=(x,y) -> ix(x.alpha) < ix(y.alpha))
    groupby(x -> x.alpha.index, vec)
end

"""
__show_by_α__

Pretty print an α_grouped derivative.
"""
function show_by_α(vec::Vector{symbolic_ξα}, vector_notation=false)
    for group in by_α(vec)
        print("α$(group[1].alpha.index)(")
        formatted = [
            "$(e.alpha.sign == 1 ? "+" : "-")$(e.xi)"
            for e in group
        ]
        print(join(formatted, " "))
        println(")")
    end
end
