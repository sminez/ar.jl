#=
    A simple test runner script to collect tests and display results.

    Testing is being done using FactCheck.
        >>> https://github.com/JuliaArchive/FactCheck.jl

    Test contexts have been defined for logical concepts that are then
    iterated over for the possible input domain. Where possible I am
    trying to check every case, if not then at least every allowed
    input type.
=#
using FactCheck
include("src/AR.jl")
using AR

# #######################################################
# .: Bring in non-exported elements for use in tests :. #
# #######################################################
ALLOWED = AR.ALLOWED
METRIC = AR.METRIC
DIVISION_TYPE = AR.DIVISION_TYPE
units = ["0" "1" "2" "3"]

tests = []

for (root, dirs, files) in walkdir("tests/")
    for file in files
        startswith(file, "test_") && push!(tests, joinpath(root, file))
    end
end

println("Running test suite:...\n")

for test in tests
    println(".: $test :.")
    include(test)
end

# NOTE:: FactCheck suppresses the exit status so to re-raise it
#        at the end of the test run uncomment the following line.
# FactCheck.exitstatus()
