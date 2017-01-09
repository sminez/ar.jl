#=
Assorted operators that are neither products or differentials
    Commutators
=#

"""Group theoretic commutator: [a,b] = a.b.a-1.b-1"""
group_commutator(a::α, b::α) = a * b * inv(a) * inv(b)

"""Ring theoretic commutator: [a,b] = a.b - b.a"""
ring_commutator(a::ξα, b::ξα) = a*b - b*a

"""Ring theoretic anti-commutator: {a,b} = a.b + b.a"""
ring_anticommutator(a::ξα, b::ξα) = a*b + b*a
