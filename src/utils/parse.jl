#=
Parsing Ξ expressions
=====================

This is a mix of utility functions that are intended to enable conversion
from traditional mathematical syntax to an internal representation that can
be manipulated efficiently.

The main work being done so far is focusing on symbolic manipulation of the
algebra only whereas this would allow for analytic / numberic solutions to be
computed.
=#
const UNITS = Dict(zip("txyz", [α("0"), α("1"), α("2"), α("3")]))


""" __f_strings__
Parse a string of the form f"αμν[f(t,x,y,z) _+_ g(t,x,y,z)] into an
internal data structure that can be used to perform calculations.

Present limitations
- The function must be expressed as a linear sum of bracketed terms, each
  preceeded by an α value.
- The terms must be separated by a comma and no other commas may appear in the
  input string. (Open to suggestions for another separator!)
- The functions within the braces are currently expected to be a single
  function with arguments within brackets.

Future work
- Grouping of alpha terms shouldn't be too difficult: (α10 - α23)[sin(t)].
- Allow for polynomial functions.
- Allow for arbitrary complexity in functions
"""
function parse_function(s::String)
    components = split(s, r"\s?,\s?")
    func = []

    for component in components
        _wrt = Set{α}()
        match_α = match(
            r"(?<sign>[+-]?)α(?<alpha>\d+)\[(?<func>.*)\] ?",
            component
        )

        match_f_args = match(
            # TODO:: pull out mathematical functions and actually differentiate
            #r"(?<mfunc>.*)\((?<args>.*)\)",
            r".*\((?<args>.*)\)",
            match_α["func"]
        )

        _sign = match_α["sign"] == "-" ? -1 : 1
        _α = α(String(match_α["alpha"]), _sign)
        _comp = parse(match_α["func"])
        args = match_f_args["args"]

        for arg in args
            # Store alpha values that we need to calculate derivatives wrt
            arg in "txyz" && push!(_wrt, UNITS[arg])
        end
        push!(func, function_ξα(_α, _comp, _wrt))
    end
    return func
end

"""String macro version of parse_function :: f"<equation>"."""
macro f_str(s)
    parse_function(s)
end
