using LinearAlgebra

struct god
    ẑero::∃
    ône::∃
    ∂t₀::Bool
    v::T
    ρ::T
    Ω::𝕋
    ⚷::UInt
    ♯::NTuple
    ∇::UInt
    d::Function
end
function god(; d, μ, ρ, ⚷=zero(UInt), Φ=○̂, ♯=SA[1], ∇=typemax(UInt))
    N = length(d)
    ∂₀ = SVector(ntuple(_ -> (true, false), N))
    ẑero = ∃(Ω, d, μ, ρ, ∂₀, Φ)
    # μ₁ = SA[μ[1], ones(T, N - 1)...]
    ones = @SVector ones(T, N)
    ∂₁ = SVector(ntuple(_ -> (false, true), N))
    zeros = @SVector zeros(T, N)
    ône = ∃(Ω, d, ones, zeros, ∂₁, ○̂)
    god(ẑero, ône, true, zero(T), zero(T), 𝕋(), ⚷, ♯, ∇, (x, y) -> sqrt.(x .^ 2 .+ y .^ 2))
end
# isreal(ϵ::∃) = √(ϵ) === Ω
dh(⚷, g, n) = powermod(g, ⚷, n)
⚷⚷(⚷, dh) = powermod(dh, ⚷, DH_N)
⚷i(i, ⚷, ♯) =
    ntuple(length(♯)) do d
        mod1(i[d] + ⚷, ♯[d])
    end
i⚷(i, ⚷, ♯) =
    ntuple(length(♯)) do d
        î = mod(⚷, ♯[d])
        mod1(i[d] + ♯[d] - î, ♯[d])
    end
function ∃!(g::god, Φ, ω=Ω)
    ϵ = -(g.ône, g.ẑero, ω)
    ṫ = t(ω.Ο[ω] + 1)
    ρt = (one(T) - ṫ) * ○
    μt = ṫ + ρt
    μ = SA[μt, ϵ.μ[2:end]...]
    ρ = SA[ρt, ϵ.ρ[2:end]...]
    ϵ = ∃(ϵ, ϵ.d, μ, ρ, ϵ.∂, Φ)
    # ϵ = ∃(ϵ, ϵ.d, ϵ.μ, ϵ.ρ, ϵ.∂, ℼ(Φ))
    ∃!(ϵ, ω)
end

function ∃̇(g::god, ∇̄=1, ω=Ω)
    ϵ = β(-(g.ône, g.ẑero,ω), ω)
    # todo dim==2
    ẑeroμ = ϵ isa 𝕋 ? zeros(T, 3) : g.ẑero.μ[end-2:end]
    ôneμ = ϵ isa 𝕋 ? ones(T, 3) : g.ône.μ[end-2:end]
    ex, ey, wx, wy = calc_ew(ẑeroμ, ôneμ)
    i = fill(ϵ isa 𝕋 ? 0 : 1, g.♯..., GL_N)
    ΦΦ = []
    !(ϵ isa 𝕋) && push!(ΦΦ, ϵ.Φ)
    owners!(g, ϵ, i, ΦΦ, 0, ∇̄, ex, ey, wx, wy, ω)
    project(g, ΦΦ, x, i)
end
# function owners!(g, ϵ::𝕋, i, ΦΦ, ∇, ex, ey, wx, wy)
# end
# ϵ = Ω.ϵ̃[ϵ][1]
function owners!(g, ϵ, i, ΦΦ, ∇, ∇̄, ex, ey, wx, wy, ω)
    if 0 < ∇ && ϵ isa ∃
        intersects = pyramid_box_intersection(
            i, length(ΦΦ) + 1,
            g.ẑero.μ, g.ône.μ,
            ex, ey, wx, wy,
            ϵ.μ .- ϵ.ρ, ϵ.μ .+ ϵ.ρ,
            g.♯..., GL_N)
        intersects || return
        push!(ΦΦ, ϵ.Φ)
    end
    ∇ == ∇̄ && return
    for ϵ̃ = ω.ϵ̃[ϵ]
        owners!(g, ϵ̃, i, ΦΦ, ∇ + 1, ∇̄, ex, ey, wx, wy, ω)
    end
end
# z, o=ẑeroμ, ôneμ
function calc_ew(z, o)
    d = o - z
    u = SVector(d[1] == d[2] == 0 ? (0, d[3], -d[2]) : (d[2], -d[1], 0))
    u = u / norm(u)
    v = cross(d, u)
    v = v / norm(v)
    # dot(d, u)
    # dot(d, v)
    # dot(u, v)
    a = abs.(u)
    b = abs.(v)
    L = abs.(d)
    A = []
    # i = 1
    for i = 1:3
        0 < a[i] || continue
        0 < b[i] || continue
        wx = L[i] / a[i] / 4
        wy = L[i] / b[i] / 4
        a[i] * wx + b[i] * wy ≤ L[i] / 2 || continue
        j = (i + 1) % 3 + 1
        a[j] * wx + b[j] * wy ≤ L[j] / 2 || continue
        j = (j + 1) % 3 + 1
        a[j] * wx + b[j] * wy ≤ L[j] / 2 || continue
        push!(A, (wx, wy, wx * wy))
    end
    # i = 1
    # j = 2
    for i = 1:2, j = i+1:3
        D = a[i] * b[j] - a[j] * b[i]
        wx = (L[i] * b[j] - L[j] * b[i]) / D / 2
        wy = (L[j] * a[i] - L[i] * a[j]) / D / 2
        0 < wx || continue
        0 < wy || continue
        k = 6 - i - j
        a[k] * wx + b[k] * wy ≤ L[k] / 2 || continue
        push!(A, (wx, wy, wx * wy))
    end
    _, maxi = findmax(a -> a[3], A)
    wx = A[maxi][1]
    wy = A[maxi][2]
    ex = u * wx
    ey = v * wy
    ex, ey, wx, wy
end

struct ΦTuple{ΦΦ}
    ϕ::ΦΦ
end
@generated function Φ(φ::ΦTuple{ΦΦ}, i, x) where {ΦΦ}
    N = length(ΦΦ.parameters)
    branches = []
    for ĩ in 1:N
        push!(branches, quote
            if i == $ĩ
                return φ.φ[$ĩ](x)
            end
        end)
    end
    quote
        $(branches...)
        return zero(T)
    end
end
function project!(out, Φ, i, e1, e2)
    o = zero(T)
    for zi = 1:GL_N
        ĩ = @index(Global, NTuple)
        ϕi = i[ĩ..., zi]
        ϕ = Φ[ϕi]
        x =
            o += ϕ(x)
    end
    out[ĩ] = one(T) - exp(-o)
end
function project(g::god, Φ::ΦTuple{ΦΦ}, i) where {ΦΦ}
    out = KernelAbstractions.zeros(GPU_BACKEND, T, g.♯[2], g.♯[3])
    i̇ = KernelAbstractions.allocate(GPU_BACKEND, UInt16, size(i))
    copyto!(i̇, i)
    project!(GPU_BACKEND, GPU_BACKEND_WORKGROUPSIZE)(
        out,
        Φ, i̇,
        ndrange=(g.♯[2], g.♯[3])
    )
    KernelAbstractions.synchronize(GPU_BACKEND)
end

# const WHITE = (one(T), one(T), one(T), one(T))
# const BLACK = (zero(T), zero(T), zero(T), one(T))
# ∃̇(g::god) = ∃̇(g.ône - g.ẑero, g.♯, g.∇)
# # ϵ=g.ône - g.ẑero
# function ∃̇(ϵ::∃, ♯, ∇)
#     ϵ̂ = X(ϵ, ♯, ∇)
#     ϵ∃ = filter(ϵ -> ϵ !== God, vec(ϵ̂))
#     isempty(ϵ∃) && return fill(WHITE, ♯[2], ♯[3])
#     unique!(t, ϵ∃)
#     sort!(ϵ∃, by=t)
#     ϵt = map(t, ϵ∃)
#     t_to_tag = Dict(ϵt[i] => UInt32(i) for i in eachindex(ϵt))
#     Φi = zeros(UInt32, ♯[2], ♯[3], ♯[4])
#     for i in CartesianIndices(ϵ̂)
#         ϵ̂ᵢ = ϵ̂[i]
#         ϵ̂ᵢ === God && continue
#         Φi[i[2], i[3], i[4]] = t_to_tag[t(ϵ̂ᵢ)]
#     end
#     φ = ΦSet(ntuple(i -> ϵ∃[i].Φ, length(ϵ∃)))
#     rgba = render(φ, Φi, ♯)
#     ℼ̂(rgba)
# end
# ℼ̂(ϕ) = begin
#     pixel = fill(WHITE, size(ϕ, 2), size(ϕ, 3))
#     # i = collect(CartesianIndices(pixel))[1]
#     for i = CartesianIndices(pixel)
#         r = ϕ[1, Tuple(i)...]
#         g = ϕ[2, Tuple(i)...]
#         b = ϕ[3, Tuple(i)...]
#         a = ϕ[4, Tuple(i)...]
#         pixel[i] = r == g == b == a == ○ ? WHITE : (r, g, b, a)
#     end
#     pixel
# end
# # ρ(μ) = min(μ, 1 - μ)
function step(g::god, dt̂=one(T))
    if g.∂t₀
        ṫ = t()
        μ = SVector(ntuple(i -> i == 1 ? ṫ : g.ẑero.μ[i], length(g.ẑero.μ)))
    else
        δμ = g.ône.μ .- g.ẑero.μ
        all(d -> iszero(d), δμ) && return
        α = clamp(g.v * dt̂, zero(T), one(T))
        μ = g.ẑero.μ .+ α .* δμ
    end
    g = move(g, μ)
    # g.ẑero.Φ !== ○̂ && ∃!(g.ẑero)
    g
end
jerk(g::god, δ) = accelerate(g, g.v * exp(δ))
accelerate(g::god, δ) = speed(g, iszero(g.v) ? δ : g.v * exp(δ))
speed(g::god, v) = god(g.ẑero, g.ône, g.∂t₀, clamp(T(v), zero(T), one(T)), g.ρ, g.Ω, g.⚷, g.♯, g.∇, g.d)
# stop(g::god) = speed(g, zero(T))
# stoptime(g::god) = god(g.ẑero, g.ône, g.∂t₀, SA[zero(T), g.v[2:end]...], g.ρ, g.Ω, g.⚷, g.♯, g.∇, g.d)
# δ=T(0.01)
scale(g::god, δ) = begin
    ϵ = g.ône - g.ẑero
    ρ = min.(ϵ.ρ * exp(δ), ○)
    N = length(ϵ.μ)
    ône = min.(ϵ.μ .+ ρ, ones(T, N))
    ẑero = max.(zeros(T, N), ϵ.μ .- ρ)
    move(g, ône) # could be parallel
    move(g, ẑero) # could be parallel
end
# ẑeroμ=ône
# ẑeroμ=ẑero
move(g::god, ẑeroμ) =
    god(
        ∃(g.ẑero.ϵ̂, g.ẑero.d, ẑeroμ, g.ẑero.ρ, g.ẑero.∂, g.ẑero.Φ),
        # ∃(g.ône.ϵ̂, g.ône.d, SA[ẑeroμ[1], g.ône.μ[2:end]...], g.ône.ρ, g.ône.∂, g.ône.Φ),
        g.ône,
        g.∂t₀, g.v, g.ρ, g.Ω, g.⚷, g.♯, g.∇, g.d
    )
focus(g::god, ôneμ) =
    god(
        # ∃(g.ẑero.ϵ̂, g.ẑero.d, SA[ôneμ[1], g.ẑero.μ[2:end]...], g.ẑero.ρ, g.ẑero.∂, g.ẑero.Φ),
        g.ẑero,
        ∃(g.ône.ϵ̂, g.ône.d, ôneμ, g.ône.ρ, g.ône.∂, g.ône.Φ),
        g.∂t₀, g.v, g.ρ, g.Ω, g.⚷, g.♯, g.∇, g.d
    )
move(g::god, d, δ) = move(g, SVector(ntuple(i -> begin
        g.ẑero.d[i] == d && return δ(μ[i], T(0.01))
        μ[i]
    end, length(μ))))
focus(g::god, d, δ) = focus(g, SVector(ntuple(i -> begin
        g.ône.d[i] == d && return δ(μ[i], T(0.01))
        μ[i]
    end, length(μ))))
focusup(g, d) = focus(g, d, +)
focusdown(g, d) = focus(g, d, -)
moveup(g, d) = move(g, d, +)
movedown(g, d) = move(g, d, -)
jerkdown(g) = jerk(g, T(-0.01))
jerkup(g) = jerk(g, T(0.01))
scaledown(g) = scale(g, T(-0.01))
scaleup(g) = scale(g, T(0.01))

# home(g::god) = god(g.ẑero, g.ône, true, zero(T), g.ρ, g.Ω, g.⚷, g.♯, g.∇)

# struct ΦSet{Fs}
#     fs::Fs  # Tuple of Φ functions
# end
# # eval_Φ(φ,1,0.5,0.5,0.5,0.4)
# # @generated function eval_Φ(φ::ΦSet{Fs}, idx, t, x, y, z) where Fs
# @generated function eval_Φ(φ::ΦSet{Fs}, idx, x) where Fs
#     N = length(Fs.parameters)
#     branches = []
#     for i in 1:N
#         push!(branches, quote
#             if idx == $i
#                 # return φ.fs[$i](t, x, y, z)
#                 return φ.fs[$i](x)
#             end
#         end)
#     end
#     quote
#         $(branches...)
#         return (zero(T), zero(T), zero(T), zero(T))
#     end
# end
# @kernel function κ!(rgba, φ::ΦSet, Φi, ♯)
#     xi, yi = @index(Global, NTuple)
#     # xi, yi = 2,2
#     _, W, H, D = ♯
#     x = isone(W) ? ○ : (T(xi) - 1) / T(W - 1)
#     y = isone(H) ? ○ : (T(yi) - 1) / T(H - 1)
#     r, g, b, a = zero(T), zero(T), zero(T), zero(T)
#     # zi = collect(1:D)[2]
#     for zi = 1:D
#         one(T) ≤ a && break
#         z = isone(D) ? ○ : T(zi - 1) / T(D - 1)
#         idx = Φi[xi, yi, zi]
#         iszero(idx) && continue
#         # ṙ, ġ, ḃ, ȧ = eval_Φ(φ, idx, ○, x, y, z)
#         ṙ, ġ, ḃ, ȧ = eval_Φ(φ, idx, (○, x, y, z))
#         iszero(ȧ) && continue
#         rem = one(T) - a
#         r += ṙ * ȧ * rem
#         g += ġ * ȧ * rem
#         b += ḃ * ȧ * rem
#         a += ȧ * rem
#     end
#     rgba[1, xi, yi] = r
#     rgba[2, xi, yi] = g
#     rgba[3, xi, yi] = b
#     rgba[4, xi, yi] = a
# end
# # ♯=g.♯
# # φ, Φi, ♯
# function render(φ::ΦSet, Φi, ♯)
#     rgba = KernelAbstractions.zeros(GPU_BACKEND, T, 4, ♯[2], ♯[3])
#     i̇ = KernelAbstractions.allocate(GPU_BACKEND, UInt32, size(Φi))
#     copyto!(i̇, Φi)
#     Base.invokelatest() do
#         κ!(GPU_BACKEND, GPU_BACKEND_WORKGROUPSIZE)(
#             rgba, φ, i̇, ♯,
#             ndrange=(♯[2], ♯[3])
#         )
#     end
#     KernelAbstractions.synchronize(GPU_BACKEND)
#     Array(rgba)
# end

# # i = collect(CartesianIndices(Ξ))[2678]
# # Ξ[i].Φ(1)
# # ∇=typemax(UInt32)
# function X(ϵ::∃, ♯, ∇)
#     Ξ = Array{∀}(undef, ♯...)
#     ρ₀ = zero(ϵ.ρ)
#     # for i in CartesianIndices(Ξ)
#     # î = collect(1:length(Ξ))[1]
#     @time Threads.@threads for î in 1:length(Ξ)
#         # @time begin 
#         i = CartesianIndices(Ξ)[î]
#         x = X(i, ♯)
#         # xϵ = ∃(God, ϵ.d, x, ρ₀, ϵ.∂, ϵ.Φ)
#         xϵ = ∃(God, ϵ.d, SVector(x), ρ₀, ϵ.∂, ϵ.Φ)
#         Ξ[i], _ = X(xϵ, ∇)
#         # end;
#     end
#     Ξ
# end


# sum([0.707,-0.707,0] .* [1,1,1])
# sum([0.408,0.408,-0.816] .* [1,1,1])

const GL_N = 8
const GL_NODES = SVector{GL_N,T}(
    -0.9602898564975363,
    -0.7966664774136267,
    -0.5255324099163290,
    -0.1834346424956498,
    0.1834346424956498,
    0.5255324099163290,
    0.7966664774136267,
    0.9602898564975363
)
const GL_WEIGHTS = SVector{GL_N,T}(
    0.1012285362903763,
    0.2223810344533745,
    0.3137066458778873,
    0.3626837833783620,
    0.3626837833783620,
    0.3137066458778873,
    0.2223810344533745,
    0.1012285362903763,
)