VERSION >= v"0.4.0-dev+6521" && __precompile__(true)

##########################################
#= Calculating with Absolute Relativity =#
##########################################
module AR

######################
# .: Dependencies :. #
######################
import Base.*, Base./, Base.\, Base.==, Base.show


####################################
# .: Exported types and methods :. #
####################################
export  α, alpha,
        AR_pair,
        ξα, xi_component,
        Ξ, xi,
        view_cayley,

        ∂, partial,
        Dμ, D_mu,
        ∇, del


###########################
# .: Load module files :. #
###########################
for (dir, filename) in [
    ("algebra", "algebra.jl"),
    ("algebra", "differential.jl")
    ]
    include(joinpath(dir, filename))
end

# End of module AR
end
