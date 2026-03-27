import Pkg
Pkg.activate(".")

"""
TheoryOfGod

I = [ZERO < ○ < ONE] denotes a unit 1-dim space of information with origin ○ (no information) in its center including the corners ZERO and ONE.
Ω = I^I an ∞-dim metric and smooth vector space.
We have a Pretopology 𝕋 on Ω such that ϵᵢ ∈ 𝕋:
* ϵᵢ ⊆ Ω
* ϵ₂ ∈ ϵ₁.ϵ̃ => ϵ₂|ϵ₁ ⊆ ϵ₁ <=> ϵ₂ ⫉ ϵ₁ ⩓ ϵ₂ ∈ ϵ₃.ϵ̃ => ϵ₁ = ϵ₃
* x ∈ ϵᵢ ⊊ Ω: x.ρ = 0 => ϵᵢ.Φ(x) ∈ I is arbitrary, computable and smooth fuzzy existence potential towards ONE=true xor ZERO=false.

ϵ ⊊ Ω defines its existence inside a subset of Ω using an origin (μ), a radius (ρ) and a closed vs. open in each direction (∂) vector. These vectors are finite and all other dimensional coordinates of ϵ follow from linear interpolation.
If we use a horizontal axis for dimension and a vertical axis for coordinate in the dimension, for any ϵ, the chart looks like a stepwise linear function with finite non-zero radius intervals (active dimensions) and zero interval points within the interpolated regions.
Each child ϵ is a subset of its parent in the active dimensions declared by the parent.

god ⊊ God = Ω = I^I = I^(.) = [ZERO < ○ < ONE]^(.)
god observes or creates, God iterates.
"""
# module tog

# export ∃, ∃̇, ∃!

const T = Float32

using KernelAbstractions, StaticArrays, LinearAlgebra, Adapt
using HTTP, URIs, Sockets
using PNGFiles
using MiniFB
using Metal
const GPU_BACKEND = MetalBackend()
const GPU_BACKEND_WORKGROUPSIZE = 2^2^3

include("∃.jl")
const Ω = Ref(𝕋())
t() = t(Ω[])
const invϕ = one(T) / MathConstants.golden
# const name = Dict{∃, String}()
include("god.jl")
include("Octahedron.jl")
include("MiniFB.jl")
include("browser.jl")
include("godBrowser.jl")
const BROWSER_TASK = Threads.@spawn start(godbrowserstart, godbrowserkeypress)

g = godBROWSER[].g
browser=godBROWSER[].browser
g = god(
        t=zero(T),
        d=sort(SA[invϕ, invϕ^2, one(T)]), # t, x, y, z
        ẑeroμ=SA[○-T(0.0), ○-T(0.0), ○],
        ôneμ=SA[○+T(0.0), ○+T(0.0), ○+T(0.1)],
        ρ=(T(0.1), T(0.1), zero(T)),
        ♯=(10, 10))
ω = Ω[]
Ω[].Ο[Ω[]]
t()
g.ẑero.μ
g.ône.μ
g.ρ
g.θ
Ω[] === ω
Ω[].ϵ̃[Ω[]]
Ω[].ϵ̃[Ω[]][1]
Ω[].ϵ̃[Ω[]][2]
Ω[].ϵ̃[Ω[].ϵ̃[Ω[]][1]][1]

# dx, dy, d, μ, ρ, N=dxdy(g)
∃!(g, x -> T(0.1), Ω[])
∃!(g, x -> T(0.2), Ω[])
∃!(g, x -> T(0.3), Ω[])
ϕ̇ = Base.invokelatest() do
                        ∃̇(g, Ω[])
                    end
unique(ϕ̇)

focus!(g, 2, T(0.2))
move!(g, 2, T(0.2))
scale!(g, (T(0.05),T(0.05),zero(T)))

# ∃!(g, x -> begin
#         abs(x[2] - x[3]) < T(0.01) && return one(T)
#         # abs([2] - x[3]) < min(g.ρ[1],g.ρ[2])*T(0.01) && return one(T)
#         return ○
#     end, Ω[])

# @time Base.invokelatest() do
#                         ∃̇(g, Ω[]);
#                     end;

# # ∃!(g, (x...) -> prod(x), Ω[])
# ∃!(g, x -> prod(x), Ω[])
step!(g)
unique(ϕ̇)
# g.∂t₀=false
# focusup!(g, 2)
# focus!(g, SA[g.ône.μ[1],T(0.3),T(0.3),g.ône.μ[3]])
# focusup!(g, 3)
# focusdown!(g, 2)
# focusdown!(g, 3)
# moveup!(g, 2)
# moveup!(g, 1)
# movedown!(g, 2)
# movedown!(g, 1)
# scaleup!(g, 2)
# scaledown!(g, 1)
# scaleup!(g, 2)
# scaleup!(g, 3)
# g.ρ=(g.ρ[1],g.ρ[2],zero(T))
# scaledown!(g, 2)
# norm(g.ône.μ .- g.ẑero.μ)
# move!(g, SA[t(), ○, ○, ○])
# focus!(g, SA[t(), ○*exp(T(0.1)), ○*exp(T(0.1)), ○*exp(T(0.1))])
# scale!(g, (T(0.1), T(0.1), one(T)))
# jerkup!(g)
# speed!(g, T(0.001))
# Ω[].ϵ̃[Ω[]][1].Φ(SA[0.0,0.0,0.0,0.0])

include("Typst.jl")
# ∃!(g, typst("o"), Ω[])
∃!(g, typst("abcd"), Ω[])

# ∃!(g, typst("ii"), Ω[])

# g[], δ = step(g[], zero(T))
ω = Ω[]
# g=G[]
# ∃̇(G[], Ω[])

# const Ω = 𝕋()

# Ω.Ο[Ω.ϵ̃[Ω][1]]
# Ω[].Ο[Ω[]]
# t(Ω.ϵ̃[Ω][1])
# t(Ω.ϵ̃[Ω.ϵ̃[Ω][1]][1])

# G[].ẑero.μ
# G[].ẑero.ρ
# G[].ône.μ
# G[].ône.ρ
# G[].ρ
# G[], δ = step(G[])
# g = G[]
# ∃!(g, (x...)->T(0.1))
# g = move(g, SA[t(Ω.ϵ̃[Ω][1]), T(0.75), T(0.75), T(0.75)])
# ∃!(g, (x...)->T(0.2))
# g = move(g, SA[t(Ω.ϵ̃[Ω][1]), T(0.5), T(0.5), T(0.5)])

# println("idx=2 (tex): ", Array(out)[1])

# # μ(::𝕋) = SVector(ntuple(_ -> ○, length(d)))
# # μ(ϵ::∃) = ϵ.μ
# # ρ(Ω::𝕋) = μ(Ω)
# # ρ(ϵ::∃) = ϵ.ρ


