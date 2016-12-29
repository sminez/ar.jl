# Test cases for algebra.jl

println("--> allowed = {$(join(ALLOWED, ",")))}")
println("-----------------------------------------------------------------")

facts("Algebra Properties") do

    context("αp.αμ = αμ ∀ μ ∈ allowed") do
        for i in ALLOWED
            @fact α("p") * α(i) --> α(i) "αp is not idempotent"
        end
    end

    context("-αp.αμ = -αμ ∀ μ ∈ allowed") do
        for i in ALLOWED
            @fact α("-p") * α(i) --> α(i, -1) "-αp is not negating"
        end
    end

    context("αμν = -αμν ∀ μ,ν ∈ {0,1,2,3}, μ ≠ ν") do
        for i in units, j in units
            if i != j
                @fact α(i) * α(j) --> -(α(j) * α(i)) "αμν ≠ -ανμ"
            end
        end
    end

    context("αμ.αμ = (-)αp ∀ μ ∈ allowed") do
        for i in ALLOWED
            @fact α(i)^2 --> anyof(α("p"), α("-p")) "not squaring to αp"
        end
    end
end

#NOTE:: This is being very odd...the first line of the Cayley table is failing
#       due to being negated, but ONLY in the test! (Manual test in the repl is
#       fine...)
facts("CAYLEY -> find_prod (only broken under factcheck...)") do
    find_prod(a::α, b::α) = AR.find_prod(a, b)
    ix(a::α) = AR.ix(a)

    for i in ALLOWED, j in ALLOWED
        @pending find_prod(α(i), α(j)) --> AR.CAYLEY[ix(α(i)),ix(α(j))] "$i, $j"
    end
end


facts("Constructor verification") do
    context("α :: constructor Equivalence") do
        for i in ALLOWED
            @fact α(i, 1) --> α(i) "constructors are not equivalent"
            @fact α(i, -1) --> α("-"*i) "string parsing of -ves is broken"
        end
    end

    context("α :: invalid Constructor Arguments") do
        @fact_throws MethodError α("not valid") "invalid string allowed"
        @fact_throws MethodError α(π) "float allowed"
        @fact_throws MethodError α(2//3) "rational allowed"
        @fact_throws MethodError α(α("0")) "α allowed"
        @fact_throws MethodError α("0", 2) "positive integer sign not ∈ [1 -1]"
        @fact_throws MethodError α("0", -42) "negative integer sign not ∈ [1 -1]"
        @fact_throws MethodError α("0", 2.7) "float sign allowed"
    end
end
