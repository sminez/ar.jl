#=
A selection of different implementations for numerically computing
the 4-vector 4-differentail Dμ.

In Cartesian coordinates this is:
    Dμ = ∂μ/αμ = ∂ / (αμ∂xμ) = α0∂0 - αi∂i = α0∂0 - ∇i
=#
using algebra


"Numeric Partial derivative of vec with respect to component var"
function ∂(vec::Ξμ, var::Int)
    throw("Not yet implemented!")
end

"3-vector ∇"
function ∇(vec::Ξμ)
    throw("Not yet implemented!")
end


"D4Ξ4 = α0∂0 - α1∂1 - α2∂2 - α3∂3"
function Dμ(vec::Ξμ)
    t = ∂(vec, 0)
    x, y, z = ∇i(vec)
    return Ξμ(t, -x, -y, -z)
end

