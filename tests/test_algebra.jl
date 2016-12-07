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


facts("Constructor verification :: α") do
    context("Constructor Equivalence") do
        for i in ALLOWED
            @fact α(i, 1) --> α(i) "constructors are not equivalent"
            @fact α(i, -1) --> α("-"*i) "string parsing of -ves is broken"
        end
    end

    context("Invalid Constructor Arguments") do
        @fact_throws MethodError α("not valid") "invalid string allowed"
        @fact_throws MethodError α(π) "float allowed"
        @fact_throws MethodError α(2//3) "rational allowed"
        @fact_throws MethodError α(α("0")) "α allowed"
        @fact_throws MethodError α("0", 2) "positive integer sign not ∈ [1 -1]"
        @fact_throws MethodError α("0", -42) "negative integer sign not ∈ [1 -1]"
        @fact_throws MethodError α("0", 2.7) "float sign allowed"
    end
end
