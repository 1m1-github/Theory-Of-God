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
    N = length(d)+1
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

function octahedron(g::god)
    N = length(g.ẑero.d)
    d = g.ône.μ .- g.ẑero.μ
    D = g.norm(d)^2
    perm = sortperm(abs.(d[end-2:end]))
    i1, i2 = N - 3 + perm[1], N - 3 + perm[2]
    e1 = zeros(SVector{N,T})
    e1 = setindex(e1, one(T), i1)
    d1 = e1 - dot(e1, d) / D * d
    d1 /= g.norm(d1)
    e2 = zeros(SVector{N,T})
    e2 = setindex(e2, one(T), i2)
    d2 = e2 - dot(e2, d) / D * d
    d2 = d2 - dot(d2, d1) * d1
    d2 /= g.norm(d2)
    dx = cos(g.θ) * d1 + sin(g.θ) * d2
    dy = -sin(g.θ) * d1 + cos(g.θ) * d2
    dx *= 2 * g.ρ[1]
    dy *= 2 * g.ρ[2]
    μ = (g.ône.μ .+ g.ẑero.μ) * ○
    ρ = (abs.(dx) + abs.(dy) + abs.(d)) * ○
    dx, dy, d, μ, ρ, N
end

function ∃!(g::god, Φ, ω=g.Ω)
    g.ẑero.μ[1] < t(ω) && return
    _, _, _, μ, ρ, _ = octahedron(g)
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
        !owners!(g, ône, ϵ, i, ϵϵ, 0, dx, dy, nz, ω, istrivial)
        ΦΦ = ΦTuple(ntuple(i -> ϵϵ[i].Φ, length(ϵϵ)))
        μρϵϵ = map(ϵ -> μρΩ(ϵ), ϵϵ)
        ẑeros = SVector(ntuple(i -> μρϵϵ[i][1] .- μρϵϵ[i][2], length(μρϵϵ)))
        ônes = SVector(ntuple(i -> μρϵϵ[i][1] .+ μρϵϵ[i][2], length(μρϵϵ)))
        # ∂z = SVector(ntuple(i -> SVector(ntuple(j -> ϵϵ[i].∂[j][1], N)), length(ϵϵ)))
        # ∂o = SVector(ntuple(i -> SVector(ntuple(j -> ϵϵ[i].∂[j][2], N)), length(ϵϵ)))
        Π̂, Π, ôneϕ = if hasdepth
            z = @SVector zeros(T, N)
            ϵ = ∃(ω, g.ẑero.d, ône, z, g.ẑero.∂, ○̂)
            ϵ, found = X(ϵ, g.∇̄, ω)
            ϕ = !found || trivial(ϵ) ? ○ : ϵ.Φ(ône)
            project3d, project3d!, ϕ
        else
            project2d, project2d!, zero(T)
        end
        godẑero = μ .- (dx .+ dy) * ○
        godẑeroône = ône .- godẑero
        # Π̂(ΦΦ, Π, i, ẑeros, ônes, ∂z, ∂o, godẑero, godẑeroône, dx, dy, g.♯..., ôneϕ)
        Π̂(ΦΦ, Π, i, ẑeros, ônes, godẑero, godẑeroône, dx, dy, g.♯..., ôneϕ)
    catch e
        bt = catch_backtrace()
        showerror(stderr, e, bt)
        fill(○, g.♯...)
    end
end
function owners!(g, ône, ϵ, i, ϵϵ, ∇, dx, dy, nz, ω, istrivial)
    found = false
    if 0 < ∇ && ϵ isa ∃ && !istrivial
        intersects = pyramid_box_intersection!(
            i, length(ϵϵ) + 1,
            g.ẑero.μ, ône,
            dx, dy,
            ϵ.μ .- ϵ.ρ, ϵ.μ .+ ϵ.ρ,
            g.♯..., nz)
        intersects || return found
        push!(ϵϵ, ϵ)
        found = true
    end
    ∇ == g.∇̄ && return found
    for ϵ̃ = ω.ϵ̃[ϵ]
        found &= owners!(g, ône, ϵ̃, i, ϵϵ, ∇ + 1, dx, dy, nz, ω, trivial(ϵ̃))
    end
    found
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
# Π̂(ΦΦ, Π, i, godẑero, ẑeros, ônes, godẑeroône, dx, dy, g.♯..., ôneϕ)
# nx, ny = g.♯[1], g.♯[2]
# function project2d(ΦΦ, Π, i, ẑeros, ônes, ∂z, ∂o, godẑero, godẑeroône, dx, dy, nx, ny, ôneϕ)
function project2d(ΦΦ, Π, i, ẑeros, ônes, godẑero, godẑeroône, dx, dy, nx, ny, ôneϕ)
    out = KernelAbstractions.zeros(GPU_BACKEND, T, nx, ny)
    i̇ = KernelAbstractions.allocate(GPU_BACKEND, UInt16, size(i))
    copyto!(i̇, i)
    Π(GPU_BACKEND, GPU_BACKEND_WORKGROUPSIZE)(
        out,
        # ΦΦ, i̇, ẑeros, ônes, ∂z, ∂o, godẑero, dx, dy,
        ΦΦ, i̇, ẑeros, ônes, godẑero, dx, dy,
        nx, ny,
        ndrange=(nx, ny)
    )
    KernelAbstractions.synchronize(GPU_BACKEND)
    Array(out)
end
function project3d(ΦΦ, Π, i, godẑero, ẑeros, ônes, godẑeroône, dx, dy, nx, ny, ôneϕ)
    out = KernelAbstractions.zeros(GPU_BACKEND, T, nx, ny)
    i̇ = KernelAbstractions.allocate(GPU_BACKEND, UInt16, size(i))
    copyto!(i̇, i)
    Π(GPU_BACKEND, GPU_BACKEND_WORKGROUPSIZE)(
        out,
        ΦΦ, i̇, godẑero, godẑeroône, dx, dy, ôneϕ,
        nx, ny, GL_N - 1, GL_NODES, GL_WEIGHTS,
        ndrange=(nx, ny)
    )
    KernelAbstractions.synchronize(GPU_BACKEND)
    Array(out)
end
# @kernel function project2d!(out, ΦΦ, i, ẑeros, ônes, ∂z, ∂o, godẑero, dx, dy, nx, ny)
@kernel function project2d!(out, ΦΦ, i, ẑeros, ônes, godẑero, dx, dy, nx, ny)
    (ix, iy) = @index(Global, NTuple)
    # ix, iy = 2,2
    iϕ = i[ix, iy]
    if iszero(iϕ)
        out[ix, iy] = ○
    else
        ĩx = T(ix - 1) / T(nx - 1)
        ĩy = T(iy - 1) / T(ny - 1)
        x = godẑero .+ ĩx * dx .+ ĩy * dy
        ẋ = Base.setindex(x, godẑero[1], 1)
        zlocal = ẑeros[iϕ]
        olocal = ônes[iϕ]
        good = true
        for k = 2:4
            # ∂zϵ = ∂z[iϕ][k]
            # ∂oϵ = ∂o[iϕ][k]
            if x[k] < zlocal[k] ||
               olocal[k] < x[k] ||
               olocal[k] == zlocal[k]
            # if x[k] < zlocal[k] || (x[k] == zlocal[k] && ∂zϵ[k]) ||
            #    olocal[k] < x[k] || (x[k] == olocal[k] && ∂oϵ[k]) ||
            #    olocal[k] == zlocal[k]
                out[ix, iy] = ○
                good = false
                break
            end
            ẋ = Base.setindex(ẋ, (x[k] - zlocal[k]) / (olocal[k] - zlocal[k]), k)
            # ẋ[k] = (x[k] .- zlocal[k]) ./ (olocal[k] .- zlocal[k]) # todo case olocal < zlocal
        end
        if good
            # ẋ = Base.setindex(, (x .- zlocal) ./ (olocal .- zlocal), 1)
            out[ix, iy] = Φ̇(ΦΦ, iϕ, ẋ) # todo clamp to [0,1]
        end
    end
end
# todo does x actually belong to ϵ
@kernel function project3d!(out, ΦΦ, i, godẑero, godẑeroône, dx, dy, ôneϕ, nx, ny, nz, gl_nodes, gl_weights)
    (ix, iy) = @index(Global, NTuple)
    ĩx = T(ix - 1) / T(nx - 1)
    ĩy = T(iy - 1) / T(ny - 1)
    ϕ = ôneϕ
    for iz = 1:nz
        iϕ = i[ix, iy, iz]
        if iszero(iϕ)
            ϕ += ○ * gl_weights[iz]
            continue
        end
        t = gl_nodes[iz]
        t̃ = one(T) - t
        z = godẑero .+ t * godẑeroône
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
rotate!(g::god, θ) = g.θ = θ
scale!(g::god, ρ) = g.ρ = ρ # todo handle too large god
scale!(g::god, i, δ) = scale!(g, ntuple(ĩ -> begin
        ĩ == i && return g.ρ[ĩ] + δ
        g.ρ[ĩ]
    end, length(g.ρ)))
move!(g::god, ẑeroμ) = begin
    try
        any(g.ône.μ[2:end] .< ẑeroμ[2:end]) && return # todo handle ône<ẑero
        g.ẑero = ∃(g.ẑero.ϵ̂, g.ẑero.d, ẑeroμ, g.ẑero.ρ, g.ẑero.∂, g.ẑero.Φ)
        g.ône = ∃(g.ône.ϵ̂, g.ône.d, SA[ẑeroμ[1], g.ône.μ[2:end]...], g.ône.ρ, g.ône.∂, g.ône.Φ)
    catch
    end
end
focus!(g::god, ôneμ) = begin
    try
        any(ôneμ[2:end] .< g.ẑero.μ[2:end]) && return # todo handle ône<ẑero
        g.ẑero = ∃(g.ẑero.ϵ̂, g.ẑero.d, SA[ôneμ[1], g.ône.μ[2:end]...], g.ẑero.ρ, g.ẑero.∂, g.ẑero.Φ)
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
