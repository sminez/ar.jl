using Plots
#=
Julia has both `import` and `using` for namespacing packages.
import is qualified (like Python) and using drops all exported
elements of a package into global namespace.
See here for more details:
    http://docs.julialang.org/en/release-0.5/manual/modules/

https://juliaplots.github.io/animations/
https://juliaplots.github.io/examples/gr/
https://juliaplots.github.io/output/
=#

gr()  # specify the GR backend for Plots

"A simple Gaussian pulse with offset t0 and std-dev σ"
function gauss(t::Float64, t0::Float64, σ::Float64)
  return exp(-((t - t0) ^ 2) / (2σ ^ 2))
end


# Collect will take a range and make it into a vector
x = collect(0:0.1:2π)

# The @animate macro builds an animation object that you can
# later decide what to do with.
anim = @animate for i = 0:60
    plot([gauss(k-(i/5.), 2., .5) for k in x])
    plot!((sin(x+i/5.)+1.)/2.)
    plot!((cos(x+i/5.)+1.)/2.)
end

#gif(anim, "example.gif", fps=15)
