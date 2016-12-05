@testset "AlgebraProperties" begin
    println("Algebra Properties")
    println("-----------------------------------------------------------------")
    println("--> allowed = $allowed\n")

    #####################################

    println("αp.αμ = αμ ∀ μ ∈ allowed")
    @test all([α("p")*α(j) == α(j) for j in ALLOWED])

    #####################################

    println("-αp.αμ = -αμ ∀ μ ∈ allowed")
    @test all([α("p", -1)*α(j) == α(j, -1) for j in ALLOWED])

    #####################################

    println("αμν = -αμν ∀ μ,ν ∈ {0,1,2,3}, μ != ν")
    units = ["0" "1" "2" "3"]
    ij_ji = [(α(i)*α(j) == -(α(j)*α(i))) for i in units, j in units if i != j]
    @test all(ij_ji)

    #####################################

    println("aμ.αμ = (-)αp ∀ μ ∈ allowed")
    @test all([(α(i)^2).index == "p" for i in ALLOWED])

    println("\n-----------------------------------------------------------------")
end
