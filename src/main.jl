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
# using MiniFB
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
# include("MiniFB.jl")
include("Color.jl")
include("Typst.jl")
include("browser.jl")
include("godBrowser.jl")
const BROWSER_TASK = Threads.@spawn start(godbrowserstart, godbrowserkeypress)

g = godBROWSER[].g
browser = godBROWSER[].browser
# g = god(
#         t=zero(T),
#         d=sort(SA[invϕ, invϕ^2, one(T)]), # t, x, y, z
#         ẑeroμ=SA[○-T(0.0), ○-T(0.0), ○],
#         ôneμ=SA[○+T(0.0), ○+T(0.0), ○+T(0.1)],
#         ρ=(T(0.1), T(0.1), zero(T)),
#         ♯=(10, 10))
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
∃!(g, typst("abcd"), Ω[])
∃!(g, typst("imi"), Ω[])
∃!(g, x -> T(0.1), Ω[])
∃!(g, x -> T(0.2), Ω[])
# Φ=x -> T(0.2)
∃!(g, x -> T(0.3), Ω[])
ϕ̇ = Base.invokelatest() do
        ∃̇(g, Ω[])
end
unique(ϕ̇)

focus!(g, 2, T(0.2))
moveup!(g, 3)
movedown!(g, 3)
scale!(g, (T(0.025), T(0.025), zero(T)))
scale!(g, (T(0.05), T(0.05), zero(T)))
scale!(g, (T(0.1), T(0.1), zero(T)))
scale!(g, (T(0.2), T(0.2), zero(T)))
scale!(g, (T(0.2), T(0.2), T(0.01)))
scale!(g, (T(0.2), T(0.2), one(T)))
scale!(g, 3, one(T))
scale!(g, 3, T(0.05))
scale!(g, 3, zero(T))
g.∂t₀=false
∃!(g, x -> begin
                # c = (x[1], T(0.5), T(0.5), T(0.6))
                T(0.01)^2 < (x[2] .- T(0.5))^2 + (x[3] .- T(0.5))^2 + (x[4] .- T(0.6))^2 < T(0.02)^2 && return rgba2scalar(T(1.0),zero(T),zero(T),T(0.3))
                # T(0.01) < (x .- c).^2 < T(0.02) && return rgba2scalar(one(T),zero(T),zero(T),○)
                return ○
        end, Ω[])
