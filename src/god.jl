"Octahedron from ẑero to focus."
mutable struct god
    ẑero::∃
    ône::∃
    ∂t₀::Bool
    v::T
    ρ::NTuple
    θ::T
    Ω::𝕋
    ⚷::UInt
    ♯::NTuple
    ∇̄::UInt
    norm::Function
end
function god(; d, t, ẑeroμ, ôneμ, ρ, θ=zero(T), ⚷=zero(UInt), ♯=(2, 2), ∇̄=typemax(UInt), n̂orm=x -> sqrt(sum(x̃ -> x̃^2, x)))
    @assert all(zero(T) .< d)
    N = length(d) + 1
    ∂ = SVector(ntuple(_ -> (true, true), N))
    z = @SVector zeros(T, N)
    ẑero = ∃(Ω[], SA[zero(T), d...], SA[t, ẑeroμ...], z, ∂, ○̂)
    ône = ∃(Ω[], SA[zero(T), d...], SA[t, ôneμ...], z, ∂, ○̂)
    god(ẑero, ône, true, zero(T), ρ, θ, 𝕋(), ⚷, ♯, ∇̄, n̂orm)
end

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

octahedron(g::god) = octahedron(g.ẑero.μ, g.ône.μ, g.ρ, g.θ, g.norm)
function octahedron(ẑeroμ, ôneμ, gρ, θ, n̂orm)
    N = length(ẑeroμ)
    d = ôneμ .- ẑeroμ
    D = n̂orm(d)^2
    perm = sortperm(abs.(d[end-2:end]))
    i1, i2 = N - 3 + perm[1], N - 3 + perm[2]
    e1 = zeros(SVector{N,T})
    e1 = setindex(e1, one(T), i1)
    d1 = e1 - dot(e1, d) / D * d
    d1 /= n̂orm(d1)
    e2 = zeros(SVector{N,T})
    e2 = setindex(e2, one(T), i2)
    d2 = e2 - dot(e2, d) / D * d
    d2 = d2 - dot(d2, d1) * d1
    d2 /= n̂orm(d2)
    dx = cos(θ) * d1 + sin(θ) * d2
    dy = -sin(θ) * d1 + cos(θ) * d2
    dx *= 2 * gρ[1]
    dy *= 2 * gρ[2]
    μ = (ôneμ .+ ẑeroμ) * ○
    ρ = (abs.(dx) + abs.(dy) + abs.(d)) * ○
    dx, dy, d, μ, ρ, N
end

function ∃!(g::god, Φ, ω=g.Ω, t0=t(ω), t1=one(T))
    # function ∃!(g::god, Φ, ω=g.Ω, t0=t(ω), t1=t(ω.Ο[ω]+1))
    (t0 < zero(T) || t0 < t(ω) || t1 < t0 || one(T) < t1) && return
    _, _, _, μ̃, ρ̃, _ = octahedron(g)
    try
        ρ̂ = (t1 - t0) * ○
        μ = SA[t0+ρ̂, μ̃[2:end]...]
        ρ = SA[ρ̂, ρ̃[2:end]...]
        ϵ = ∃(ω, g.ẑero.d, μ, ρ, g.ẑero.∂, Φ)
        ∃!(ϵ, ω)
    catch e
        bt = catch_backtrace()
        showerror(stderr, e, bt)
        nothing
    end
end

trivial(ϵ) = ϵ isa 𝕋 || ϵ.Φ === ○̂
function ∃̇(g::god, ω=g.Ω)
    try
        dx, dy, d, μ, ρ, N = octahedron(g)
        ϵ = ∃(ω, g.ẑero.d, μ, ρ, g.ẑero.∂, g.ẑero.Φ)
        ϵ = β(ϵ, ω, ω)
        istrivial = trivial(ϵ)
        ϵϵ = []
        !istrivial && push!(ϵϵ, ϵ)
        ôneρ = SA[zero(T), fill(g.ρ[end], N - 1)...]
        ône = g.ẑero.μ .+ d * ○ .* (one(T) .+ ôneρ)
        hasdepth = !iszero(g.ρ[end])
        nz = hasdepth ? GL_N - 1 : 1
        i = fill(istrivial ? 0 : 1, g.♯..., nz)
        owners!(g, ône, ϵ, i, ϵϵ, 0, dx, dy, nz, ω, istrivial)
        # unique(i)
        # count(x->x==1,i)/prod(size(i))
        # count(x->x==2,i)/prod(size(i))
        isempty(ϵϵ) && return fill(○, g.♯...)
        ΦΦ = ΦTuple(ntuple(i -> ϵϵ[i].Φ, length(ϵϵ)))
        μρϵϵ = map(ϵ -> μρΩ(ϵ), ϵϵ)
        ẑeros = SVector(ntuple(i -> μρϵϵ[i][1] .- μρϵϵ[i][2], length(μρϵϵ)))
        ônes = SVector(ntuple(i -> μρϵϵ[i][1] .+ μρϵϵ[i][2], length(μρϵϵ)))
        ∂z = SVector(ntuple(i -> ϵϵ[ceil(Int, i / N)].∂[((i-1)%N)+1][1], N * length(ϵϵ)))
        ∂o = SVector(ntuple(i -> ϵϵ[ceil(Int, i / N)].∂[((i-1)%N)+1][2], N * length(ϵϵ)))
        ôneϕ = if hasdepth
            z = @SVector zeros(T, N)
            ϵ = ∃(ω, g.ẑero.d, ône, z, g.ẑero.∂, ○̂)
            ϵ, found = X(ϵ, g.∇̄, ω)
            !found || trivial(ϵ) ? ○ : Metal.@allowscalar ϵ.Φ(ône)
        else
            zero(T)
        end
        godẑero = μ .- (dx .+ dy) * ○
        godẑeroône = ône .- godẑero
        Π(i, ΦΦ, ôneϕ, godẑero, godẑeroône, ẑeros, ônes, ∂z, ∂o, dx, dy, g.♯..., nz)
    catch e
        bt = catch_backtrace()
        showerror(stderr, e, bt)
        fill(○, g.♯...)
    end
end
function owners!(g, ône, ϵ, i, ϵϵ, ∇, dx, dy, nz, ω, istrivial)
    if 0 < ∇ && ϵ isa ∃ && !istrivial
        μ, ρ = μρΩ(ϵ)
        intersects = pyramid_box_intersection!(
            i, length(ϵϵ) + 1,
            g.ẑero.μ, ône,
            dx, dy,
            μ .- ρ, μ .+ ρ,
            g.♯..., nz)
        intersects || return
        push!(ϵϵ, ϵ)
    end
    ∇ == g.∇̄ && return
    # ϵ̃ = ω.ϵ̃[ϵ][1]
    # ϵ=ϵ̃
    for ϵ̃ = ω.ϵ̃[ϵ]
        owners!(g, ône, ϵ̃, i, ϵϵ, ∇ + 1, dx, dy, nz, ω, trivial(ϵ̃))
    end
end

struct ΦTuple{ΦT}
    ϕ::ΦT
end
@generated function Φ̇(Φ::ΦTuple{ΦT}, i, x) where ΦT
    N = length(ΦT.parameters)
    branches = []
    for ĩ = 1:N
        push!(branches, quote
            if i == $ĩ
                return Φ.ϕ[$ĩ](x)
            end
        end)
    end
    quote
        $(branches...)
        return zero(T)
    end
end
function Π(i, ΦΦ, ôneϕ, godẑero, godẑeroône, ẑeros, ônes, ∂z, ∂o, dx, dy, nx, ny, nz)
    out = KernelAbstractions.zeros(GPU_BACKEND, T, nx, ny)
    i̇ = KernelAbstractions.allocate(GPU_BACKEND, UInt16, size(i))
    copyto!(i̇, i)
    Π!(GPU_BACKEND, GPU_BACKEND_WORKGROUPSIZE)(
        out, i̇,
        ΦΦ, ôneϕ, godẑero, godẑeroône, ẑeros, ônes, ∂z, ∂o, dx, dy,
        nx, ny, nz, GL_NODES, GL_WEIGHTS,
        ndrange=(nx, ny)
    )
    KernelAbstractions.synchronize(GPU_BACKEND)
    Array(out)
end
# nx, ny=g.♯[1],g.♯[2]
@kernel function Π!(out, i, ΦΦ, ôneϕ, godẑero, godẑeroône, ẑeros, ônes, ∂z, ∂o, dx, dy, nx, ny, nz, gl_nodes, gl_weights)
    (ix, iy) = @index(Global, NTuple)
    ĩx = T(ix - 1) / T(nx - 1)
    ĩy = T(iy - 1) / T(ny - 1)
    ϕ = ôneϕ
    hasdepth = 1 < nz
    n = length(godẑero)
    for iz = 1:nz
        iϕ = i[ix, iy, iz]
        w = hasdepth ? gl_weights[iz] : one(T)
        if iszero(iϕ)
            ϕ += ○ * w
            continue
        end
        t = gl_nodes[iz]
        t̃ = one(T) - t
        z = godẑero .+ t * godẑeroône
        d̃x = t̃ * dx
        d̃y = t̃ * dy
        x = z .+ ĩx * d̃x .+ ĩy * d̃y
        ẋ = Base.setindex(x, godẑero[1], 1)
        zlocal = ẑeros[iϕ]
        olocal = ônes[iϕ]
        inner = true
        for k = 2:n
            ∂i = (iϕ - 1) * n + k
            if x[k] < zlocal[k] ||
               olocal[k] < x[k] ||
               olocal[k] == zlocal[k] ||
               (∂z[∂i] && x[k] == zlocal[k]) || (∂o[∂i] && x[k] == olocal[k])
                ϕ += ○ * w
                inner = false
                break
            end
            ẋ = Base.setindex(ẋ, (x[k] - zlocal[k]) / (olocal[k] - zlocal[k]), k) # todo case olocal < zlocal
        end
        if inner
            ϕ += Φ̇(ΦΦ, iϕ, ẋ) * w # todo clamp to [0,1]
        end
    end
    out[ix, iy] = hasdepth ? one(T) - exp(-ϕ) : ϕ
end

CHANGEΔ = T(0.01)
function step!(g::god, dt̂=one(T))
    if g.∂t₀
        ṫ = t()
        g.ẑero.μ[1] == ṫ && return false
        μ = SVector(ntuple(i -> i == 1 ? ṫ : g.ẑero.μ[i], length(g.ẑero.μ)))
    else
        δμ = g.ône.μ .- g.ẑero.μ
        all(d -> iszero(d), δμ) && return false
        α = T(clamp(g.v * dt̂, zero(T), one(T)))
        μ = g.ẑero.μ .+ α .* δμ
    end
    move!(g, μ)
    # g.ẑero.Φ !== ○̂ && ∃!(g.ẑero)
    true
end
jerk!(g::god, δ) = accelerate!(g, g.v * exp(δ))
accelerate!(g::god, δ) = speed!(g, iszero(g.v) ? δ : g.v * exp(δ))
speed!(g::god, v) = g.v = clamp(T(v), zero(T), one(T))
rotate!(g::god, θ) = begin
    try
        _, _, _, μ, ρ, _ = octahedron(g.ẑero.μ, g.ône.μ, g.ρ, θ, g.norm)
        ∃(g.Ω, g.ẑero.d, μ, ρ, g.ẑero.∂, g.ẑero.Φ)
        g.θ = θ
    catch
    end
end
scale!(g::god, ρ) = begin
    try
        _, _, _, μ, ρ̃, _ = octahedron(g.ẑero.μ, g.ône.μ, ρ, g.θ, g.norm)
        ∃(g.Ω, g.ẑero.d, μ, ρ̃, g.ẑero.∂, g.ẑero.Φ)
        g.ρ = ρ
    catch
    end
end
scale!(g::god, i, δ) = scale!(g, ntuple(ĩ -> begin
        ĩ == i && return g.ρ[ĩ] + δ
        g.ρ[ĩ]
    end, length(g.ρ)))
move!(g::god, ẑeroμ) = begin
    try
        ôneμ = SA[ẑeroμ[1], g.ône.μ[2:end]...]
        _, _, _, μ, ρ, _ = octahedron(ẑeroμ, ôneμ, g.ρ, g.θ, g.norm)
        ∃(g.Ω, g.ẑero.d, μ, ρ, g.ẑero.∂, g.ẑero.Φ)
        # any(g.ône.μ[2:end] .< ẑeroμ[2:end]) && return # todo handle ône<ẑero
        g.ẑero = ∃(g.ẑero.ϵ̂, g.ẑero.d, ẑeroμ, g.ẑero.ρ, g.ẑero.∂, g.ẑero.Φ)
        g.ône = ∃(g.ône.ϵ̂, g.ône.d, ôneμ, g.ône.ρ, g.ône.∂, g.ône.Φ)
    catch
    end
end
focus!(g::god, ôneμ) = begin
    try
        ẑeroμ = SA[ôneμ[1], g.ẑero.μ[2:end]...]
        _, _, _, μ, ρ, _ = octahedron(ẑeroμ, ôneμ, g.ρ, g.θ, g.norm)
        ∃(g.Ω, g.ẑero.d, μ, ρ, g.ẑero.∂, g.ẑero.Φ)
        # any(ôneμ[2:end] .< g.ẑero.μ[2:end]) && return # todo handle ône<ẑero
        g.ẑero = ∃(g.ẑero.ϵ̂, g.ẑero.d, ẑeroμ, g.ẑero.ρ, g.ẑero.∂, g.ẑero.Φ)
        g.ône = ∃(g.ône.ϵ̂, g.ône.d, ôneμ, g.ône.ρ, g.ône.∂, g.ône.Φ)
    catch
    end
end
move!(g::god, i, δ) =
    move!(g, SVector(ntuple(ĩ -> begin
            ĩ == i && return g.ẑero.μ[ĩ] + δ
            g.ẑero.μ[ĩ]
        end, length(g.ẑero.μ))))
focus!(g::god, i, δ) = focus!(g, SVector(ntuple(ĩ -> begin
        ĩ == i && return g.ône.μ[ĩ] + δ
        g.ône.μ[ĩ]
    end, length(g.ône.μ))))
focusup!(g, i) = focus!(g, i, CHANGEΔ)
focusdown!(g, i) = focus!(g, i, -CHANGEΔ)
moveup!(g, i) = move!(g, i, CHANGEΔ)
movedown!(g, i) = move!(g, i, -CHANGEΔ)
jerkup!(g) = jerk!(g, CHANGEΔ)
jerkdown!(g) = jerk!(g, -CHANGEΔ)
scaleup!(g, i) = scale!(g, i, CHANGEΔ)
scaledown!(g, i) = scale!(g, i, -CHANGEΔ)

const GL_NODES_RAW_8 = SVector{8,T}(
    -0.9602898564975363,
    -0.7966664774136267,
    -0.5255324099163290,
    -0.1834346424956498,
    0.1834346424956498,
    0.5255324099163290,
    0.7966664774136267,
    0.9602898564975363
)
const GL_WEIGHTS_8 = SVector{8,T}(
    0.1012285362903763,
    0.2223810344533745,
    0.3137066458778873,
    0.3626837833783620,
    0.3626837833783620,
    0.3137066458778873,
    0.2223810344533745,
    0.1012285362903763,
)
const GL_NODES_RAW_16 = SVector{16,T}(
    -0.9894009349916499,
    -0.9445750230732326,
    -0.8656312023878318,
    -0.7554044083550030,
    -0.6178762444026438,
    -0.4580167776572274,
    -0.2816035507792589,
    -0.0950125098376374,
    0.0950125098376374,
    0.2816035507792589,
    0.4580167776572274,
    0.6178762444026438,
    0.7554044083550030,
    0.8656312023878318,
    0.9445750230732326,
    0.9894009349916499
)
const GL_WEIGHTS_16 = SVector{16,T}(
    0.0271524594117541,
    0.0622535239386479,
    0.0951585116824928,
    0.1246289712555339,
    0.1495959888165767,
    0.1691565193950025,
    0.1826034150449236,
    0.1894506104550685,
    0.1894506104550685,
    0.1826034150449236,
    0.1691565193950025,
    0.1495959888165767,
    0.1246289712555339,
    0.0951585116824928,
    0.0622535239386479,
    0.0271524594117541
)
const GL_N = 8
const GL_WEIGHTS = GL_WEIGHTS_8
const GL_NODES_RAW = GL_NODES_RAW_8
const GL_NODES = ○ .+ GL_NODES_RAW ./ (GL_NODES_RAW[end] - GL_NODES_RAW[1]) # [0,1]

# using GLMakie
# fig = Figure()
# ax = Axis3(fig[1, 1])
# for t = 0.0:0.1:1.0
#     scatter!(ax, (g.ẑero.μ[2:end] .+ t * d[2:end])...; markersize=6, color=:red)
#     scatter!(ax, (ẑero[2:end] .+ t * dẑero[2:end])...; markersize=6, color=:red)
#     scatter!(ax, (ẑero[2:end] .+ t * dx[2:end])...; markersize=6, color=:red)
#     scatter!(ax, (ẑero[2:end] .+ t * dy[2:end])...; markersize=6, color=:pink)
# end
# scatter!(ax, g.ẑero.μ[2:end]..., ; markersize=6, color=:black)
# scatter!(ax, g.ône.μ[2:end]..., ; markersize=6, color=:black)
# scatter!(ax, (g.ẑero.μ[2:end] .+ ○ * d[2:end])...; markersize=6, color=:black)
# scatter!(ax, ône[2:end]..., ; markersize=6, color=:blue)
# scatter!(ax, ẑero[2:end]..., ; markersize=6, color=:blue)
# for i = 1:size(xout, 1)
#     for j = 1:size(xout, 2)
#         for k = 1:size(xout, 3)
#             scatter!(ax, xout[i, j, k, 2], xout[i, j, k, 3], xout[i, j, k, 4])
#         end
#     end
# end
# xout[1,1,1,:]
# xout[1,1,2,:]
# xout[1,1,3,:]
# xout[1,1,4,:]
# xout[1,1,5,:]
# xout[1,1,6,:]
# xout[1,1,7,:]
# dot(dx, dy)
# dot(dx, d)
# dot(dy, d)
# g.norm(dx)
# g.norm(dy)
