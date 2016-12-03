#=
A selection of different implementations for numerically computing
the 4-vector 4-differentail Dμ.

In Cartesian coordinates this is:
    Dμ = ∂μ/αμ = ∂ / (αμ∂xμ) = α0∂0 - αi∂i = α0∂0 - ∇i

Partial Derivatives of multivector components are being defined as a vector
of (α, ξα) pairs so that we can keep track of which unit element each
component has been differnetiated with respect to in order to group at the
end of calculations. =#
typealias ∂μ Vector{Tuple{α,α,ξα}}
typealias d_mu ∂μ


############################################
# .: Derivative operations on Ξ vectors :. #
############################################
"""
**∂ | partial**

Find the partial derivative of component with respect to wrt.
This method allows you to specify a sign (+-1) in which represents calculating
∂μ vs-∂μ.
"""
function ∂(component::ξα, wrt::String)
    ∂_wrt = α(wrt)
    # Correct the sign for components that square to -αp
    ∂_wrt.sign = (∂_wrt^2).sign
    original_α = component[2]
    # Using div so that we can chose division by/into in algebra.jl
    derivative = div(component, ∂_wrt)
    return(∂_wrt::α, original_α::α, derivative::ξα)
end

∂(vec::Ξ, var::String) = [∂(comp, var) for comp in vec]
∇i(vec::Ξ) = vcat([∂(vec, i) for i in ["1" "2" "3"]]...)
Dμ(vec::Ξ) = vcat(∂(vec, "0"), ∇i(vec))

# Non-unicode versions
partial(comp::ξα, var::String) = ∂(comp::ξα, var::String)
partial(vec::Ξ, var::String) = ∂(vec::Ξ, var::String)
del_i(vec::Ξ) = ∇i(vec::Ξ)
D_mu(vec::Ξ) = Dμ(vec::Ξ)


###########################################
# .: Operations on derivatived vectors :. #
###########################################
"""
**by_α**

Group derivative terms according to their unit element.
This converts the result of a derivative - (∂,α-orig,ξα-new) - into tuple of
three α values - (α,∂,ξ) - representing the new unit element, the unit element
the derivative was taken with respect to and the origin of the ξ component in
the original Ξ multivector.

If the vector_notation flag is true then this will attempt to identify divs,
grads and curls in the resulting grouped derivative.
"""
function by_α(vec::∂μ)
    components = [(c[3][2], c[1], c[2]) for c in vec]
    sort!(components, lt=(x,y) -> ix(x[1]) < ix(y[1]))
    groupby(x -> x[1].index, components)
end

"""
**show_by_α**

Pretty print an α_grouped derivative.
"""
function show_by_α(vec::∂μ, vector_notation=false)
    for a in by_α(vec, vector_notation)
        print("α$(a[1][1].index)(")
        formatted = [
            # Sign is determine by both the sign of α and that of ∂
            # which is given by the metric in the algebra.
            "$(e[1].sign*e[2].sign > 0 ? "+" : "-")∂$(e[2].index)ξ$(e[3].index)"
            for e in a
        ]
        print(join(formatted, " "))
        println(")")
    end
end
