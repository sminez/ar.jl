#=
AR :: Absolute Relativity
-------------------------
This package defines the types and operations required for the
principle of Absolute Relativity to be applied to calculations concerning
the dynamics of the relativistic fluid.

Some elements are purely analytic and allow for symbolic manipulation of
unit elements and multi-vectors while there is also support for numerical
calculations.

As in the convention used in the theory, greek indices run 0,1,2,3
and latin indices run 1,2,3.
=#
VERSION >= v"0.4" && __precompile__(true)

##########################################
#= Calculating with Absolute Relativity =#
##########################################
module AR

######################
# .: Dependencies :. #
######################
import Base.+, Base.-, Base.*, Base./, Base.\, Base.==, Base.show, Base.div
using Iterators
using Gadfly
using Colors


####################################
# .: Exported types and methods :. #
####################################
export  α, alpha,
        symbolic_ξ,
        ξα, AR_pair,
        symbolic_ξα, symbolic_AR_pair,
        function_ξα, function_AR_pair,
        array_ξα, array_AR_pair,
        check_Ξ, check_xi,

        Ξ, xi,
        Ξp, xi_p,
        Ξ, xi,
        ΞG, xi_G,
        Ξμ, xi_1,
        Ξμν, xi_2,
        Ξμνρ, xi_3,

        ∂μ, d_mu,
        ∂, partial,
        Dμ, D_mu,
        ∇i, del_i,
        by_α,
        by_∇,
        show_by_α,

        @f_str,

        convert_cayley,
        print_cayley,
        visualise_cayley



###########################
# .: Load module files :. #
###########################
for file_path in [
    "algebra/algebra.jl",
    "algebra/differential.jl",
    "utils/visualisation.jl",
    "utils/parse.jl"
    ]
    include(file_path)
end

# End of module AR
end
