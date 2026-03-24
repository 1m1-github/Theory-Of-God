import Pkg
Pkg.activate(".")

"""
TheoryOfGod

I = [ZERO < ○ < ONE] denotes a unit 1-dim space of information with origin ○ (no information) in its center including the corners ZERO and ONE.
Ω = I^I an ∞-dim metric and smooth vector space.
We have a Pretopology 𝕋 on Ω such that ϵᵢ ∈ 𝕋:
* ϵᵢ ⊆ Ω
* ϵ₂ ∈ ϵ₁.ϵ̃ => ϵ₂|ϵ₁ ⊆ ϵ₁ <=> ϵ₂ ⫉ ϵ₁ ⩓ ϵ₂ ∈ ϵ₃.ϵ̃ => ϵ₁ = ϵ₃
* ϵ₁ ≠ ϵ₂ => ϵ₁ ∩ ϵ₂ = ∅
* x ∈ ϵᵢ ⊊ Ω: x.ρ = 0 => ϵᵢ.Φ(x) ∈ I is arbitrary, computable and smooth fuzzy existence potential towards ONE=true xor ZERO=false.

ϵ ⊊ Ω defines its existence inside a subset of Ω using an origin (μ), a radius (ρ) and a closed vs. open in each direction (∂) vector. These vectors are finite and all other dimensional coordinates of ϵ follow from linear interpolation.
If we use a horizontal axis for dimension and a vertical axis for coordinate in the dimension, for any ϵ, the chart looks like a stepwise linear function with finite non-zero radius intervals (active dimensions) and zero interval points within the interpolated regions.
Each child ϵ is a subset of its parent in the active dimensions declared by the parent.

god ⊊ God = Ω = I^I = I^(.) = [ZERO < ○ < ONE]^(.)
god observes or creates, God iterates.
"""
# module tog

# export ∃, ∃̇, ∃!

using KernelAbstractions, StaticArrays, LinearAlgebra, MiniFB
using Metal
const GPU_BACKEND = MetalBackend()
const GPU_BACKEND_WORKGROUPSIZE = 2^2^3

const T = Float32

include("00101_TheoryOfGod∃.jl")
const Ω = 𝕋()
# const name = Dict{∃, String}()
include("00103_TheoryOfGodgod.jl")
include("Octahedron.jl")

const invϕ = one(T) / MathConstants.golden
♯space = 10^1
const G = Ref(god(
    d=sort(SA[zero(T), invϕ, invϕ^2, one(T)]), # t, x, y, z
    # ẑeroμ=SA[t(), zero(T), zero(T), zero(T)],
    # f̂ocusμ=SA[t(), T(1.0), T(1.0), T(1.0)],
    ẑeroμ=SA[t(), ○, ○, ○],
    f̂ocusμ=SA[t(), T(0.6), T(0.6), T(0.6)],
    # ρ=(T(0.5),T(0.5),zero(T)),
    # ρ=(T(0.1),T(0.1),zero(T)),
    ρ=(T(0.1),T(0.1),T(0.1)),
    ♯=(♯space, ♯space)))
g=G[]
# G[] = flatten(g, g.f̂ocus.d[end])

# include("00102_TheoryOfGodMiniFB.jl")

∃!(g, x -> prod(x), Ω)
# g[], δ = step(g[], zero(T))
ω = Ω
∃̇(g, Ω)

# const Ω = 𝕋()
Ω.Ο[Ω]
# Ω.Ο[Ω.ϵ̃[Ω][1]]
# Ω.Ο[ϵ]
t()
# t(Ω.ϵ̃[Ω][1])
# t(Ω.ϵ̃[Ω.ϵ̃[Ω][1]][1])
# Ω.ϵ̃
g.ẑero.μ
g.ẑero.ρ
g.f̂ocus.μ
g.f̂ocus.ρ
g.ρ
G[] = step(G[]);
g=G[]
# ∃!(g, (x...)->T(0.1))
# g = move(g, SA[t(Ω.ϵ̃[Ω][1]), T(0.75), T(0.75), T(0.75)])
# ∃!(g, (x...)->T(0.2))
# g = move(g, SA[t(Ω.ϵ̃[Ω][1]), T(0.5), T(0.5), T(0.5)])

function start()
    # t = time()
    while true
        sleep(1)
        yield()
        
        @show "before pending_actions"
        while isready(pending_actions)
            action = take!(pending_actions)
            global g = action(g)
        end
        @show "after pending_actions"
        
        # t̂ = time()
        # dt = t̂ - t
        # t = t̂
        # g̃, _ = step(g, dt)
        # global g = g̃
        # println("δ=$δ")
        # δ || continue

        println("start() ẑero.μ=$(g.ẑero.μ)")
        println("start() f̂ocus.μ=$(g.f̂ocus.μ)")

        # p̂ixel = ∃̇(g)
        # @show p̂ixel
        # δ = Δ(pixel, p̂ixel)
        # @show δ
        # @show isempty(δ)
        # isempty(δ) && continue
        # todo
        # global buffer = p̂ixel
        @show "before updatebuffer"
        global buffer = updatebuffer(g)
        @show "after updatebuffer"
    end
end
# const godTASK = @async start()

# include("00104_TheoryOfGodTypst.jl")
# Φ_hi = Φ_typst(typst_to_matrix("hi"))
# Φ_hi = Φ_typst("hi")
# gpu_safe(○̂, 2)
# gpu_safe((x,y) -> ○, 2)
# φ_hi, mat_hi = Φ_typst("hi")
# gpu_safe(φ_hi, 2)
# ∃!(g, φ_hi)

# println("idx=2 (tex): ", Array(out)[1])

# # μ(::𝕋) = SVector(ntuple(_ -> ○, length(d)))
# # μ(ϵ::∃) = ϵ.μ
# # ρ(Ω::𝕋) = μ(Ω)
# # ρ(ϵ::∃) = ϵ.ρ


# # include("00090_BroadcastBrowser2Module.jl")
# # import Main.BroadcastBrowserModule: BroadcastBrowser, start
# # include("00105_TheoryOfGodgodBrowser.jl")
# # const BROWSERTASK = Threads.@spawn start(b -> godBrowser(b))
# # g=collect(values(godBROWSER[]))[1].g
