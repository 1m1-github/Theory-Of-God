"Octahedron from ẑero to focus."
mutable struct god
    ẑero::∃
    f̂ocus::∃
    ∂t₀::Bool
    v::T
    ρ::NTuple
    Ω::𝕋
    ⚷::UInt
    ♯::NTuple
    ∇̄::UInt
    norm::Function
end
function god(; d, ẑeroμ, f̂ocusμ, ρ, ⚷=zero(UInt), ♯=(2, 2), ∇̄=typemax(UInt), n̂orm=x -> sqrt(sum(x̃ -> x̃^2, x)))
    N = length(d)
    ∂ = SVector(ntuple(_ -> (true, true), N))
    z = @SVector zeros(T, N)
    ẑero = ∃(Ω[], d, ẑeroμ, z, ∂, ○̂)
    f̂ocus = ∃(Ω[], d, f̂ocusμ, z, ∂, ○̂)
    god(ẑero, f̂ocus, true, zero(T), ρ, 𝕋(), ⚷, ♯, ∇̄, n̂orm)
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

function dxdy(g::god)
    N = length(g.ẑero.μ)
    d = g.f̂ocus.μ .- g.ẑero.μ
    D = g.norm(d)^2
    ey = SA[zeros(T, N - 2)..., one(T), zero(T)]
    dy = ey .- d[end-1] / D * d
    dy /= g.norm(dy)
    ex = SA[zeros(T, N - 3)..., one(T), zeros(T, 2)...]
    dx = ex .- d[end-2] / D * d
    dx = dx .- dot(dx, dy) * dy
    dx /= g.norm(dx)
    dx *= 2 * g.ρ[1]
    dy *= 2 * g.ρ[2]
    μ = (g.f̂ocus.μ .+ g.ẑero.μ) * ○
    ρ = (abs.(dx) + abs.(dy) + abs.(d)) * ○
    dx, dy, d, μ, ρ, N
end

function ∃!(g::god, Φ, ω=g.Ω)
    _, _, _, μ̇, ρ̇, _ = dxdy(g)
    ṫ = t(ω.Ο[ω] + 1)
    ρt = (one(T) - ṫ) * ○ # todo create non-eternally
    μt = ṫ + ρt
    μ = SA[μt, μ̇[2:end]...]
    ρ = SA[ρt, ρ̇[2:end]...]
    try
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
        dx, dy, d, μ, ρ, N = dxdy(g)
        ϵ = ∃(ω, g.ẑero.d, μ, ρ, g.ẑero.∂, g.ẑero.Φ)
        ϵ = β(ϵ, ω, ω)
        istrivial = trivial(ϵ)
        i = fill(istrivial ? 0 : 1, g.♯..., GL_N - 1)
        Φ̃Φ̃ = []
        !istrivial && push!(Φ̃Φ̃, ϵ.Φ)
        f̂ocusρ = SA[zero(T), fill(g.ρ[end], N - 1)...]
        f̂ocus = g.ẑero.μ .+ d * ○ .* (one(T) .+ f̂ocusρ)
        hasdepth = !iszero(g.ρ[end])
        nz = hasdepth ? GL_N - 1 : 1
        owners!(g, f̂ocus, ϵ, i, Φ̃Φ̃, 0, dx, dy, nz, ω, istrivial)
        ΦΦ = ΦTuple(ntuple(i -> Φ̃Φ̃[i], length(Φ̃Φ̃)))
        Π̂, Π, f̂ocusϕ = if hasdepth
            z = @SVector zeros(T, N)
            ϵ = ∃(ω, g.ẑero.d, f̂ocus, z, g.ẑero.∂, ○̂)
            ϵ, found = X(ϵ, g.∇̄, ω)
            ϕ = !found || trivial(ϵ) ? ○ : ϵ.Φ(f̂ocus)
            project3d, project3d!, ϕ
        else
            project2d, project2d!, zero(T)
        end
        ẑero = μ .- (dx .+ dy) * ○
        dẑero = f̂ocus .- ẑero
        Π̂(ΦΦ, Π, i, ẑero, dẑero, dx, dy, g.♯..., f̂ocusϕ)
    catch e
        bt = catch_backtrace()
        showerror(stderr, e, bt)
        fill(○, g.♯...)
    end
end
function owners!(g, f̂ocus, ϵ, i, ΦΦ, ∇, dx, dy, nz, ω, istrivial)
    if 0 < ∇ && ϵ isa ∃ && !istrivial
        intersects = pyramid_box_intersection!(
            i, length(ΦΦ) + 1,
            g.ẑero.μ, f̂ocus,
            dx, dy,
            ϵ.μ .- ϵ.ρ, ϵ.μ .+ ϵ.ρ,
            g.♯..., nz)
        intersects || return
        push!(ΦΦ, ϵ.Φ)
    end
    ∇ == g.∇̄ && return
    for ϵ̃ = ω.ϵ̃[ϵ]
        owners!(g, f̂ocus, ϵ̃, i, ΦΦ, ∇ + 1, dx, dy, nz, ω, trivial(ϵ̃))
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
function project2d(ΦΦ, Π, i, ẑero, d, dx, dy, nx, ny, f̂ocusϕ)
    out = KernelAbstractions.zeros(GPU_BACKEND, T, nx, ny)
    i̇ = KernelAbstractions.allocate(GPU_BACKEND, UInt16, size(i))
    copyto!(i̇, i)
    Π(GPU_BACKEND, GPU_BACKEND_WORKGROUPSIZE)(
        out,
        ΦΦ, i̇, ẑero, dx, dy,
        nx, ny,
        ndrange=(nx, ny)
    )
    KernelAbstractions.synchronize(GPU_BACKEND)
    Array(out)
end
function project3d(ΦΦ, Π, i, ẑero, d, dx, dy, nx, ny, f̂ocusϕ)
    out = KernelAbstractions.zeros(GPU_BACKEND, T, nx, ny)
    i̇ = KernelAbstractions.allocate(GPU_BACKEND, UInt16, size(i))
    copyto!(i̇, i)
    Π(GPU_BACKEND, GPU_BACKEND_WORKGROUPSIZE)(
        out,
        ΦΦ, i̇, ẑero, d, dx, dy, f̂ocusϕ,
        nx, ny, GL_N - 1, GL_NODES, GL_WEIGHTS,
        ndrange=(nx, ny)
    )
    KernelAbstractions.synchronize(GPU_BACKEND)
    Array(out)
end
# todo does x actually belong to ϵ
@kernel function project2d!(out, ΦΦ, i, ẑero, dx, dy, nx, ny)
    (ix, iy) = @index(Global, NTuple)
    iϕ = i[ix, iy]
    if iszero(iϕ)
        out[ix, iy] = ○
    else
        ĩx = T(ix - 1) / T(nx - 1)
        ĩy = T(iy - 1) / T(ny - 1)
        x = ẑero .+ ĩx * dx .+ ĩy * dy
        out[ix, iy] = Φ̇(ΦΦ, iϕ, x) # todo clamp to [0,1]
    end
end
# todo does x actually belong to ϵ
@kernel function project3d!(out, ΦΦ, i, ẑero, d, dx, dy, f̂ocusϕ, nx, ny, nz, gl_nodes, gl_weights)
    (ix, iy) = @index(Global, NTuple)
    ĩx = T(ix - 1) / T(nx - 1)
    ĩy = T(iy - 1) / T(ny - 1)
    ϕ = f̂ocusϕ
    for iz = 1:nz
        iϕ = i[ix, iy, iz]
        if iszero(iϕ)
            ϕ += ○ * gl_weights[iz]
            continue
        end
        t = gl_nodes[iz]
        t̃ = one(T) - t
        z = ẑero .+ t * d
        d̃x = t̃ * dx
        d̃y = t̃ * dy
        x = z .+ ĩx * d̃x .+ ĩy * d̃y
        ϕ += Φ̇(ΦΦ, iϕ, x) * gl_weights[iz] # todo clamp to [0,1]
    end
    out[ix, iy] = one(T) - exp(-ϕ)
end

function step!(g::god, dt̂=one(T))
    if g.∂t₀
        ṫ = t()
        g.ẑero.μ[1] == ṫ && return false
        μ = SVector(ntuple(i -> i == 1 ? ṫ : g.ẑero.μ[i], length(g.ẑero.μ)))
    else
        δμ = g.f̂ocus.μ .- g.ẑero.μ
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
scale!(g::god, ρ) = g.ρ = ρ
scale!(g::god, i, δ) = scale!(g, ntuple(ĩ -> begin
        ĩ == i && return g.ρ[ĩ] + δ
        g.ρ[ĩ]
    end, length(g.ρ)))
move!(g::god, ẑeroμ) = begin
    try
        g.ẑero = ∃(g.ẑero.ϵ̂, g.ẑero.d, ẑeroμ, g.ẑero.ρ, g.ẑero.∂, g.ẑero.Φ)
        g.f̂ocus = ∃(g.f̂ocus.ϵ̂, g.f̂ocus.d, SA[ẑeroμ[1], g.f̂ocus.μ[2:end]...], g.f̂ocus.ρ, g.f̂ocus.∂, g.f̂ocus.Φ)
    catch
    end
end
focus!(g::god, ôneμ) = begin
    try
        g.ẑero = ∃(g.ẑero.ϵ̂, g.ẑero.d, SA[ôneμ[1], g.ẑero.μ[2:end]...], g.ẑero.ρ, g.ẑero.∂, g.ẑero.Φ)
        g.f̂ocus = ∃(g.f̂ocus.ϵ̂, g.f̂ocus.d, ôneμ, g.f̂ocus.ρ, g.f̂ocus.∂, g.f̂ocus.Φ)
    catch
    end
end
move!(g::god, i, δ) =
    move!(g, SVector(ntuple(ĩ -> begin
            ĩ == i && return g.ẑero.μ[ĩ] + δ
            g.ẑero.μ[ĩ]
        end, length(g.ẑero.μ))))
focus!(g::god, i, δ) = focus!(g, SVector(ntuple(ĩ -> begin
        ĩ == i && return g.f̂ocus.μ[ĩ] + δ
        g.f̂ocus.μ[ĩ]
    end, length(g.f̂ocus.μ))))
focusup!(g, i) = focus!(g, i, T(0.01))
focusdown!(g, i) = focus!(g, i, -T(0.01))
moveup!(g, i) = move!(g, i, T(0.01))
movedown!(g, i) = move!(g, i, -T(0.01))
jerkup!(g) = jerk!(g, T(0.01))
jerkdown!(g) = jerk!(g, T(-0.01))
scaleup!(g, i) = scale!(g, i, T(0.01))
scaledown!(g, i) = scale!(g, i, -T(0.01))

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
# scatter!(ax, g.f̂ocus.μ[2:end]..., ; markersize=6, color=:black)
# scatter!(ax, (g.ẑero.μ[2:end] .+ ○ * d[2:end])...; markersize=6, color=:black)
# scatter!(ax, f̂ocus[2:end]..., ; markersize=6, color=:blue)
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