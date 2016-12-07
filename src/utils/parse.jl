#=
Parsing Ξ expressions
=====================

This is a mix of utility functions that are intended to enable conversion
from traditional mathematical syntax to an internal representation that can
be manipulated.
=#
UNITS = Dict(zip("txyz", ["0", "1", "2", "3"]))

type eqn_component
    expr::Expr
    wrt::Set{α}
    alpha::α
end


"""
__equation"" strings__
Parse a string of the form `equation"αμν[f(t,x,y,z) _+_ g(t,x,y,z)]` into an
internal data structure that can be used to perform calculations
"""
macro equation_str(s)
    components = split(s, r"\s?_._\s?")
    binops = matchall(r"_(.)_", s)
    equation = []

    for component in components
        _wrt = Set{α}()
        match_α = match(r"α(?<alpha>\d*)\[(?<func>.*)\] ?", component)
        _α = match_α["alpha"]
        _func = match_α["func"]

        match_f_args = match(r"(?<mfunc>.*)\((?<strargs>.*)\)", _func)
        _mfunc = match_f_args["mfunc"]
        _strargs = match_f_args["strargs"]

        for arg in _strargs
            if arg in "txyz"
                # TODO:: negative alphas
                push!(_wrt, α(UNITS[arg]))
            end
        end
        push!(equation, eqn_component(parse(_func), _wrt, α(String(_α))))
        length(binops) > 0 && push!(equation, parse(shift!(binops)))
    end
    return equation
end
