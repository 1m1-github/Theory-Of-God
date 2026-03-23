"Octahedron from ẑero to focus."
struct god
    ẑero::∃
    f̂ocus::∃
    ∂t₀::Bool
    v::T
    # ρ::NTuple{2,T}
    ρ::NTuple
    Ω::𝕋
    ⚷::UInt
    ♯::NTuple
    # ♯::NTuple{2,UInt}
    ∇̄::UInt
    norm::Function
end
function god(; d, ẑeroμ, f̂ocusμ, ρ, ⚷=zero(UInt), ♯=(2, 2), ∇̄=typemax(UInt), n̂orm=x -> sqrt(sum(x̃ -> x̃^2, x)))
    N = length(d)
    ∂ = SVector(ntuple(_ -> (true, true), N))
    z = @SVector zeros(T, N)
    ẑero = ∃(Ω, d, ẑeroμ, z, ∂, ○̂)
    f̂ocus = ∃(Ω, d, f̂ocusμ, z, ∂, ○̂)
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
    dnorm = g.norm(d)
    dx = SA[zeros(T, N - 3)..., one(T), zeros(T, 2)...]
    dxort = dx .- d[end-2] / dnorm * d
    dxnorm = g.norm(dxort)
    dy = SA[zeros(T, N - 2)..., one(T), zero(T)]
    dyort = dy .- d[end-1] / dnorm * d
    dynorm = g.norm(dyort)
    dz = SA[zeros(T, N - 1)..., one(T)]
    dzort = dz .- d[end] / dnorm * d
    dznorm = g.norm(dzort)
    dx = dxort
    dy = dyort
    if dxnorm < dynorm && dxnorm < dznorm
        dx = dzort
    elseif dynorm < dxnorm && dynorm < dznorm
        dy = dzort
    end
    dx = dx / g.norm(dx)
    dy = dy .- dot(dy, dx)*dx
    dy = dy * 2 * g.ρ[2] / g.norm(dy)
    dx *= 2 * g.ρ[1]
    μ = (g.f̂ocus.μ .+ g.ẑero.μ) * ○
    ρ = (abs.(dx) + abs.(dy) + abs.(d)) * ○
    dx, dy, μ, d, ρ, N
end
function ∃!(g::god, Φ, ω=g.Ω)
    _, _, μ̃, _, ρ̃, _ = dxdy(g)
    ṫ = t(ω.Ο[ω] + 1)
    ρt = (one(T) - ṫ) * ○ # todo create non-eternally
    μt = ṫ + ρt
    μ = SA[μt, μ̃[2:end]...]
    ρ = SA[ρt, ρ̃[2:end]...]
    ϵ = ∃(ω, g.ẑero.d, μ, ρ, g.ẑero.∂, Φ)
    ∃!(ϵ, ω)
end

trivial(ϵ) = ϵ isa 𝕋 || ϵ.Φ === ○̂
function ∃̇(g::god, ω=g.Ω)
    dx, dy, μ, d, ρ, _ = dxdy(g)
    ϵ = ∃(ω, g.ẑero.d, μ, ρ, g.ẑero.∂, g.ẑero.Φ)
    ϵ = β(ϵ, ω, ω)
    istrivial = trivial(ϵ)
    i = fill(istrivial ? 0 : 1, g.♯..., GL_N - 1)
    Φ̃Φ̃ = []
    !istrivial && push!(Φ̃Φ̃, ϵ.Φ)
    owners!(g, ϵ, i, Φ̃Φ̃, 0, dx, dy, ω, istrivial)
    ΦΦ = ΦTuple(ntuple(i -> Φ̃Φ̃[i], length(Φ̃Φ̃)))
    hasdepth = !iszero(d[end])
    Π = hasdepth ? project3d! : project2d!
    out = project(g, ΦΦ, Π, i, dx, dy)
    if hasdepth
        ϕ = ○
        if !trivial(g.f̂ocus)
            ϵ = X(g.f̂ocus, g.∇̄, ω)
            ϕ = ϵ.Φ(g.f̂ocus.μ)
        end
        out .= out .+ ϕ
    end
    out
end
function owners!(g, ϵ, i, ΦΦ, ∇, dx, dy, ω, istrivial)
    if 0 < ∇ && ϵ isa ∃ && !istrivial
        intersects = pyramid_box_intersection!(
            i, length(ΦΦ) + 1,
            g.ẑero.μ, g.f̂ocus.μ,
            dx, dy,
            ϵ.μ .- ϵ.ρ, ϵ.μ .+ ϵ.ρ,
            g.♯..., GL_N - 1)
        intersects || return
        push!(ΦΦ, ϵ.Φ)
    end
    ∇ == g.∇̄ && return
    for ϵ̃ = ω.ϵ̃[ϵ]
        owners!(g, ϵ̃, i, ΦΦ, ∇ + 1, dx, dy, ω, trivial(ϵ̃))
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
function project(g, ΦΦ, Π, i, dx, dy)
    out = KernelAbstractions.zeros(GPU_BACKEND, T, g.♯...)
    i̇ = KernelAbstractions.allocate(GPU_BACKEND, UInt16, size(i))
    copyto!(i̇, i)
    Base.invokelatest() do
        Π(GPU_BACKEND, GPU_BACKEND_WORKGROUPSIZE)(
            out,
            ΦΦ, i̇, g.ẑero.μ, g.f̂ocus.μ, dx, dy, g.♯...,
            ndrange=g.♯
        )
    end
    KernelAbstractions.synchronize(GPU_BACKEND)
    Array(out)
end
# todo does x actually belong to ϵ
@kernel function project2d!(out, ΦΦ, i, ẑero, f̂ocus, dx, dy, nx, ny)
    (ix, iy) = @index(Global, NTuple)
    iϕ = i[ix, iy]
    if iszero(iϕ)
        out[ix, iy] = zero(T)
    else
        ĩx = T(ix - 1) / T(nx - 1)
        ĩy = T(iy - 1) / T(ny - 1)
        x = ẑero .+ ĩx * dx .+ ĩy * dy
        out[ix, iy] = Φ̇(ΦΦ, iϕ, x)
    end
end
# todo does x actually belong to ϵ
@kernel function project3d!(out, ΦΦ, i, ẑero, f̂ocus, dx, dy, nx, ny)
    ϕ = zero(T)
    (ix, iy) = @index(Global, NTuple)
    d = f̂ocus .- ẑero
    for iz = 1:GL_N-1
        iϕ = i[ix, iy, iz]
        if iszero(iϕ)
            ϕ += ○ * GL_WEIGHTS[iz]
            continue
        end
        t = GL_NODES[iz]
        ĩx = T(ix - 1) / T(nx - 1)
        ĩy = T(iy - 1) / T(ny - 1)
        zt = ẑero .+ t * d
        x = zt .+ ĩx * dx .+ ĩy * dy
        ϕ += Φ̇(ΦΦ, iϕ, x) * GL_WEIGHTS[iz]
    end
    out[ix, iy] = one(T) - exp(-ϕ)
end




function step(g::god, dt̂=one(T))
    if g.∂t₀
        ṫ = t()
        g.ẑero.μ[1] == ṫ && return g, false
        μ = SVector(ntuple(i -> i == 1 ? ṫ : g.ẑero.μ[i], length(g.ẑero.μ)))
    else
        δμ = g.f̂ocus.μ .- g.ẑero.μ
        all(d -> iszero(d), δμ) && return g, false
        α = clamp(g.v * dt̂, zero(T), one(T))
        μ = g.ẑero.μ .+ α .* δμ
    end
    g = move(g, μ)
    # g.ẑero.Φ !== ○̂ && ∃!(g.ẑero)
    g, true
end
jerk(g::god, δ) = accelerate(g, g.v * exp(δ))
accelerate(g::god, δ) = speed(g, iszero(g.v) ? δ : g.v * exp(δ))
speed(g::god, v) = god(g.ẑero, g.f̂ocus, g.∂t₀, clamp(T(v), zero(T), one(T)), g.ρ, g.Ω, g.⚷, g.♯, g.∇̄, g.norm)
# stop(g::god) = speed(g, zero(T))
# stoptime(g::god) = god(g.ẑero, g.f̂ocus, g.∂t₀, SA[zero(T), g.v[2:end]...], g.ρ, g.Ω, g.⚷, g.♯, g.∇̄, g.norm)
scale(g::god, δ) = begin
    # scale(g::god, δ, ω=Ω) = begin
    # ϵ = -(g.f̂ocus, g.ẑero, ω)
    # N = length(ϵ.μ)
    # ône = SVector(min.(ϵ.μ .+ ρ, ones(T, N)))
    # ẑero = SVector(max.(zeros(T, N), ϵ.μ .- ρ))
    # move(g, ône) # todo could be parallel
    # move(g, ẑero) # todo could be parallel
    ρ = max.(min.(g.ρ * exp(δ), ○), zero(T))
    god(g.ẑero, g.f̂ocus, g.∂t₀, g.v, ρ, g.Ω, g.⚷, g.♯, g.∇̄, g.norm)
end
move(g::god, ẑeroμ) = begin
    @show "move2", ẑeroμ
    god(
        ∃(g.ẑero.ϵ̂, g.ẑero.d, ẑeroμ, g.ẑero.ρ, g.ẑero.∂, g.ẑero.Φ),
        # ∃(g.f̂ocus.ϵ̂, g.f̂ocus.d, SA[ẑeroμ[1], g.f̂ocus.μ[2:end]...], g.f̂ocus.ρ, g.f̂ocus.∂, g.f̂ocus.Φ),
        g.f̂ocus,
        g.∂t₀, g.v, g.ρ, g.Ω, g.⚷, g.♯, g.∇̄, g.norm
    )
end
focus(g::god, ôneμ) =
    god(
        # ∃(g.ẑero.ϵ̂, g.ẑero.d, SA[ôneμ[1], g.ẑero.μ[2:end]...], g.ẑero.ρ, g.ẑero.∂, g.ẑero.Φ),
        g.ẑero,
        ∃(g.f̂ocus.ϵ̂, g.f̂ocus.d, ôneμ, g.f̂ocus.ρ, g.f̂ocus.∂, g.f̂ocus.Φ),
        g.∂t₀, g.v, g.ρ, g.Ω, g.⚷, g.♯, g.∇̄, g.norm
    )
move(g::god, d, δ) = begin
    @show "move1", d, δ
    move(g, SVector(ntuple(i -> begin
        g.ẑero.d[i] == d && return δ(g.ẑero.μ[i], T(0.01))
        g.ẑero.μ[i]
    end, length(g.ẑero.μ))))
end
focus(g::god, d, δ) = focus(g, SVector(ntuple(i -> begin
        g.f̂ocus.d[i] == d && return δ(g.f̂ocus.μ[i], T(0.01))
        g.f̂ocus.μ[i]
    end, length(g.f̂ocus.μ))))
focusup(g, d) = focus(g, d, +)
focusdown(g, d) = focus(g, d, -)
moveup(g, d) = move(g, d, +)
movedown(g, d) = move(g, d, -)
jerkdown(g) = jerk(g, T(-0.01))
jerkup(g) = jerk(g, T(0.01))
scaledown(g) = scale(g, T(-0.01))
scaleup(g) = scale(g, T(0.01))
flatten(g, d) = begin
    ôneμ = SVector(ntuple(i -> begin
            g.f̂ocus.d[i] == d && return g.ẑero.μ[i]
            g.f̂ocus.μ[i]
        end, length(g.f̂ocus.μ)))
    god(
        g.ẑero,
        ∃(g.f̂ocus.ϵ̂, g.f̂ocus.d, ôneμ, g.f̂ocus.ρ, g.f̂ocus.∂, g.f̂ocus.Φ),
        g.∂t₀, g.v, g.ρ, g.Ω, g.⚷, g.♯, g.∇̄, g.norm
    )
end

const GL_N = 8
const GL_NODES_RAW = SVector{GL_N,T}(
    -0.9602898564975363,
    -0.7966664774136267,
    -0.5255324099163290,
    -0.1834346424956498,
    0.1834346424956498,
    0.5255324099163290,
    0.7966664774136267,
    0.9602898564975363
)
const GL_NODES = ○ .+ GL_NODES_RAW ./ (GL_NODES_RAW[end] - GL_NODES_RAW[1]) # [0,1]
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
