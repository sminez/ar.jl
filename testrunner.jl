#=
    A simple test runner script to collect tests and display results.

    NOTE:: Testsets have been defined in the tests directory and each additional
    test should print a short description of what the test is about before the
    call to @test so that it is easier to eyeball the output and see what is
    broken.
=#
using Base.Test
include("src/AR.jl")
using AR

# #######################################################
# .: Bring in non-exported elements for use in tests :. #
# #######################################################
ALLOWED = AR.ALLOWED
TARGETS = AR.TARGETS
METRIC = AR.METRIC
DIVISION_TYPE = AR.DIVISION_TYPE
# For printing the allowed elements in test descriptions
allowed = "{$(join(ALLOWED, ","))}"

tests = []

for (root, dirs, files) in walkdir("tests/")
    for file in files
        startswith(file, "test_") && push!(tests, joinpath(root, file))
    end
end

println("Running tests:...")

for test in tests
    include(test)
end
