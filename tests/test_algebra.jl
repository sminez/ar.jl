# Test cases for algebra.jl

facts("AlgebraProperties") do
    println("--> allowed = {$(join(ALLOWED, ",")))}")
    println("-----------------------------------------------------------------")

    context("αp.αμ = αμ ∀ μ ∈ allowed") do
        for i in ALLOWED
            @fact α("p")*α(i) --> α(i) "αp is not idempotent"
        end
    end

    context("-αp.αμ = -αμ ∀ μ ∈ allowed") do
        for i in ALLOWED
            @fact α("-p")*α(i) --> α(i, -1) "-αp is not negating"
        end
    end

    context("αμν = -αμν ∀ μ,ν ∈ {0,1,2,3}, μ ≠ ν") do
        units = ["0" "1" "2" "3"]
        for i in units, j in units
            if i != j
                @fact α(i)*α(j) --> -(α(j)*α(i)) "αμν ≠ -ανμ"
            end
        end
    end

    context("aμ.αμ = (-)αp ∀ μ ∈ allowed") do
        for i in ALLOWED
            @fact α(i)^2 --> anyof(α("p"), α("-p")) "not squaring to αp"
        end
    end
end
