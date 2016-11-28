#=
  This is an attempt to rewrite the fdtdlib.py file in Julia...
  Let's go!

  Useful Doc pages::
    The Manual: http://docs.julialang.org/en/latest/manual/
    Types in Julia: http://docs.julialang.org/en/release-0.5/manual/types/
    built-in functions: https://goo.gl/M1WEZh

  Curated list of packages and resources::
    https://github.com/svaksha/Julia.jl#index

  Plotting options::
    https://github.com/tbreloff/Plots.jl
    https://github.com/zyedidia/AnimatedPlots.jl
    https://github.com/jheinen/GR.jl
    http://gadflyjl.org/stable/
    https://github.com/JuliaPy/PyPlot.jl
    https://github.com/nolta/Winston.jl
=#

# Constants
const c0 = 299792458.
const μ0 = 4 * π * 1e-7
const ε0 = 1 / (μ0 * c0^2)

# Source pulse functions
"A simple Gaussian pulse with offset t0 and std-dev σ"
function gaussian(t::Float64, t0::Float64, σ::Float64)
  return exp(-((t - t0) ^ 2) / (2σ ^ 2))
end

abstract FDTD_state

type State <: FDTD_state
  n_dimensions::Int
  n_steps::Int
  boundary_values::Tuple{Number,Number}
  source_Fmax::Int
  λ_rmax::Int
  axis_lengths::Tuple{Vararg{Int}}
  # Array dimensions will be the same as the axes
  εr::Array{Float64}
  μr::Array{Float64}
  # Spacetime δs
  dt::Float64
  dx::Float64
  dy::Float64
  dz::Float64
end

#= These should be paramaters to the methods
medium::Tuple{Float64, Float64}
sources::Dict{String,Vector{Tuple{String,Any}}}
device_regions::Vector{Tuple{Float64,Float64,Any}}
=#


#= Primary Simulation Methods =#

