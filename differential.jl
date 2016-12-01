#=
A selection of different implementations for numerically computing
the 4-vector 4-differentail Dμ.

In Cartesian coordinates this is:
    Dμ = ∂μ/αμ = ∂ / (αμ∂xμ) = α0∂0 - αi∂i = α0∂0 - ∇i
=#
include("algebra.jl")


"Numeric Partial derivative of vec with respect to component var"
function ∂(vec::Ξ, var::Int)
    throw("Not yet implemented!")
end

"3-vector ∇"
function ∇(vec::Ξ)
    throw("Not yet implemented!")
end


"DμΞ4 = α0∂0Ξ4 - α1∂1Ξ4 - α2∂2Ξ4 - α3∂3Ξ4"
function Dμ(vec::Ξ)
    t = ∂(vec, 0)
    x, y, z = ∇i(vec)
    return Ξ(t, -x, -y, -z)
end
