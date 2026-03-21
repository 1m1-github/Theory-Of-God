const GL_N = 8
# const GL_NODES = (
#     -0.9602898564975363,
#     -0.7966664774136267,
#     -0.5255324099163290,
#     -0.1834346424956498,
#     0.1834346424956498,
#     0.5255324099163290,
#     0.7966664774136267,
#     0.9602898564975363,
# )
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
const GL_WEIGHTS = (
    0.1012285362903763,
    0.2223810344533745,
    0.3137066458778873,
    0.3626837833783620,
    0.3626837833783620,
    0.3137066458778873,
    0.2223810344533745,
    0.1012285362903763,
)
Kmin(v, a, b, N, i) = Int(clamp(ceil(N[i] * (v - a[i]) / (b[i] - a[i]) + 1 / 2), 1, N[i]))
Kmax(v, a, b, N, i) = Int(clamp(floor(N[i] * (v - a[i]) / (b[i] - a[i]) + 1 / 2), 1, N[i]))
function kminmax(a, b, c, d, h, t, N, dim)
    sminx = max(-t * h[dim], c[dim] - b[dim] + t * h[dim])
    smaxx = min(t * h[dim], d[dim] - b[dim] + t * h[dim])
    Lx = b[dim] - t * h[dim] + sminx
    Ux = b[dim] - t * h[dim] + smaxx
    Kmin(Lx, a, b, N, dim), Kmax(Ux, a, b, N, dim)
end
function owners!(ϵ, i, ΦΦ, ∇=1)
    μ, ρ = μρΩ(ϵ)
    dimx, dimy, dimz = length(g.ône.d) - 2, length(g.ône.d) - 1, length(g.ône.d)
    N = size(i)
    a, b = g.ẑero.μ, g.ône.μ
    h = (b .- a) / 2
    c, d = μ .- ρ, μ .+ ρ # same as a,b ?
    tminz = (b[dimz] - d[dimz]) / h[dimz]
    tmaxz = (b[dimz] - c[dimz]) / h[dimz]
    tmaxz < tminz && return
    Lz = b[dimz] - tmaxz * h[dimz]
    Uz = b[dimz] - tminz * h[dimz]
    kminz = Kmin(Lz, a, b, N, length(N))
    kmaxz = Kmax(Uz, a, b, N, length(N))
    push!(ΦΦ, ϵ.Φ)
    i̇ = length(ΦΦ)
    # k = collect(kminz:kmaxz)[1]
    for k = kminz:kmaxz
        tk = (b[dimz] - GL_NODES[k]) / h[dimz]
        kminx, kmaxx = kminmax(a, b, c, d, h, tk, N, dimx)
        kminy, kmaxy = kminmax(a, b, c, d, h, tk, N, dimy)
        i[kminx:kmaxx, kminy:kmaxy, k] .= i̇
    end
    iszero(∇) && return
    # (iϵ̃, ϵ̃) = collect(enumerate(get(God.ϵ̃,ϵ,[])))[1]
    for (iϵ̃, ϵ̃) = enumerate(get(God.ϵ̃, ϵ, []))
        owners!(ϵ̃, i, ΦΦ, ∇ - 1)
    end
end

