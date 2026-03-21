# observe.jl — god observes Ω
#
# observer = ZERO = g.ẑero       (watches)
# light    = ○    = origin        (illuminates)
# focus    = ONE  = g.ône         (direction)
# screen   = 2D disk at ○, ⊥ to ○→ONE
# rays     = cone-beam from ○ outward through screen pixels
#
# Pipeline:
#   CPU: compute GL sample points, find owner ϵ (pretopology walk)
#   GPU: eval_Φ at each point via ΦSet dispatch + GL weighted sum
#        textures via atlas, analytic Φ via ΦFunc
#
# N-dimensional rays, 2D output.

using LinearAlgebra: normalize, norm, dot
using KernelAbstractions: adapt
# using Adapt: adapt

# ── GL8 on [0,1] ────────────────────────────────────────────────────────────

const GL_NODES = SA[
    T(0.019855071751231884),  T(0.101666761293186631),
    T(0.237233795041835507),  T(0.408282678752175098),
    T(0.591717321247824902),  T(0.762766204958164493),
    T(0.898333238706813369),  T(0.980144928248768116)]
const GL_WEIGHTS = SA[
    T(0.050614268145188129),  T(0.111190517226687235),
    T(0.156853322938943644),  T(0.181341891689180991),
    T(0.181341891689180991),  T(0.156853322938943644),
    T(0.111190517226687235),  T(0.050614268145188129)]

# ── Φ wrappers (bitstype, GPU-safe) ─────────────────────────────────────────

struct ΦFunc{F}
    f::F
end
@inline (Φ::ΦFunc)(x, atlas) = Φ.f(x)

struct ΦTex3
    H::Int32
    W::Int32
    offset::Int32
end
@inline function (Φ::ΦTex3)(x, atlas)
    px = clamp(unsafe_trunc(Int32, x[1] * Φ.W) + Int32(1), Int32(1), Φ.W)
    py = clamp(unsafe_trunc(Int32, x[2] * Φ.H) + Int32(1), Int32(1), Φ.H)
    atlas[py, Φ.offset + px]
end

# ── ΦSet + @generated dispatch ───────────────────────────────────────────────

struct ΦSet{Fs}
    fs::Fs
end
ΦSet(fs...) = ΦSet(fs)

@generated function eval_Φ(φ::ΦSet{Fs}, idx, x, atlas) where Fs
    N = length(Fs.parameters)
    branches = []
    for i in 1:N
        push!(branches, quote
            if idx == $i
                return φ.fs[$i](x, atlas)
            end
        end)
    end
    quote
        $(branches...)
        return ○
    end
end

# ── gpu_safe ─────────────────────────────────────────────────────────────────

# function gpu_safe(Φ, N)
#     Φ isa ΦTex3 && return true
#     try
#         @kernel gpu(Φ, x) = Φ(x)
#         x = KernelAbstractions.zeros(GPU_BACKEND, T, N)
#         gpu(GPU_BACKEND, GPU_BACKEND_WORKGROUPSIZE)(Φ, x, ndrange=1)
#         true
#     catch
#         false
#     end
# end

# ── Typst → texture ─────────────────────────────────────────────────────────

# using PNGFiles
# using Colors: red

# const TYPST_TEMPLATE(content) = """
# #set page(width: auto, height: auto, margin: (top: 5pt, bottom: 5pt, left: 5pt, right: 5pt))
# #set text(size: 20pt)
# $content
# """
# const TYPST_CACHE = Dict{UInt,Matrix{T}}()

# function typst_to_matrix(typst_code::String; dpi=300)
#     h = hash((typst_code, dpi))
#     haskey(TYPST_CACHE, h) && return TYPST_CACHE[h]
#     cmd = `typst compile - --format png --ppi $dpi -`
#     rgba = pipeline(IOBuffer(TYPST_TEMPLATE(typst_code)), cmd) |> read |> IOBuffer |> PNGFiles.load
#     mat = Matrix{T}(one(T) .- T.(getfield.(rgba, :r)))
#     TYPST_CACHE[h] = mat
#     mat
# end

# function Φ_typst(code::String; offset::Int32=Int32(0), dpi=300)
#     mat = typst_to_matrix(code; dpi)
#     H, W = size(mat)
#     ΦTex3(Int32(H), Int32(W), offset), mat
# end

# ── Atlas builder ────────────────────────────────────────────────────────────

function build_atlas(tex_pairs)
    # tex_pairs: Vector of (ΦTex3, Matrix{T})
    isempty(tex_pairs) && return zeros(T, 1, 1), ΦTex3[]
    H_max = maximum(size(m, 1) for (_, m) in tex_pairs)
    total_W = sum(size(m, 2) for (_, m) in tex_pairs)
    atlas = fill(○, H_max, total_W)
    updated = ΦTex3[]
    col = Int32(0)
    for (φ, mat) in tex_pairs
        h, w = size(mat)
        atlas[1:h, col+1:col+w] .= mat
        push!(updated, ΦTex3(φ.H, φ.W, col))
        col += Int32(w)
    end
    atlas, updated
end

# ── Orthonormal basis: 2 vectors ⊥ v̂ in ℝᴺ ──────────────────────────────────

function screen_basis(v̂::SVector{N,T}) where {N}
    axes = sortperm(SVector(ntuple(i -> abs(v̂[i]), N)))
    e₁ = SVector(ntuple(i -> i == axes[1] ? one(T) : zero(T), N))
    û = normalize(e₁ - dot(e₁, v̂) * v̂)
    e₂ = SVector(ntuple(i -> i == axes[2] ? one(T) : zero(T), N))
    ŵ = normalize(e₂ - dot(e₂, v̂) * v̂ - dot(e₂, û) * û)
    û, ŵ
end

# ── CPU: compute sample points + find owners ─────────────────────────────────

function cpu_points_and_owners(g::god, W::Int, H::Int)
    ○_μ = (g.ẑero.μ .+ g.ône.μ) ./ 2
    v = g.ône.μ .- ○_μ
    r = norm(v)
    r < eps(T) && return nothing, nothing, nothing, r
    v̂ = v ./ r
    û, ŵ = screen_basis(v̂)
    N_dims = length(g.ẑero.d)

    points = Array{T}(undef, N_dims, 8, W * H)
    owners = zeros(UInt32, 8, W * H)

    owner_lock = ReentrantLock()
    owner_map = Dict{UInt,UInt32}()
    owner_list = ∃[]

    Threads.@threads for idx in 1:(W * H)
        py, px = divrem(idx - 1, W) .+ (1, 1)
        ξ₁ = T(2px - 1) / T(W) - one(T)
        ξ₂ = T(2py - 1) / T(H) - one(T)
        d̂ = normalize(v̂ .+ g.ρ .* (ξ₁ .* û .+ ξ₂ .* ŵ))

        for j in 1:8
            p = ○_μ .+ (GL_NODES[j] * r) .* d̂
            μ = SVector(ntuple(i -> i ≤ length(p) ? p[i] : ○, N_dims))
            for i in 1:N_dims
                points[i, j, idx] = μ[i]
            end

            x = ∃(g.ẑero.ϵ̂, g.ẑero.d, μ,
                  SVector(ntuple(_ -> zero(T), N_dims)),
                  SVector(ntuple(_ -> (true, true), N_dims)), ○̂)
            ϵ, found = X(x, g.∇)
            if found && !(ϵ isa 𝕋)
                h = hash(ϵ)
                lock(owner_lock)
                if !haskey(owner_map, h)
                    push!(owner_list, ϵ)
                    owner_map[h] = UInt32(length(owner_list))
                end
                owners[j, idx] = owner_map[h]
                unlock(owner_lock)
            end
        end
    end
    points, owners, owner_list, r
end

# ── Build ΦSet + atlas from owner list ───────────────────────────────────────

function build_φset(owner_list)
    tex_pairs = Tuple{ΦTex3, Matrix{T}}[]
    analytic_fs = []
    owner_types = Symbol[]  # :tex or :func per owner

    for ϵ in owner_list
        if ϵ.Φ isa ΦTex3
            push!(tex_pairs, (ϵ.Φ, TYPST_CACHE[hash(ϵ.Φ)]))  # TODO: better mat lookup
            push!(owner_types, :tex)
        else
            push!(analytic_fs, ΦFunc(ϵ.Φ))
            push!(owner_types, :func)
        end
    end

    atlas, updated_tex = build_atlas(tex_pairs)

    # build final tuple: replace tex entries with updated offsets
    ti, fi = 0, 0
    fs = map(owner_types) do t
        if t == :tex
            ti += 1
            updated_tex[ti]
        else
            fi += 1
            analytic_fs[fi]
        end
    end

    ΦSet(Tuple(fs)), atlas
end

# ── GPU kernel: eval Φ + GL weighted sum ─────────────────────────────────────

@kernel function κ_observe!(img, @Const(pts), @Const(own), @Const(gl_w),
                            @Const(r), @Const(N_dims), φ::ΦSet, @Const(atlas))
    idx = @index(Global)
    W = size(img, 2)
    py, px = divrem(idx - 1, W) .+ (1, 1)

    val = zero(T)
    for j in 1:8
        ow = own[j, idx]
        if iszero(ow)
            val += gl_w[j] * ○
        else
            x = ntuple(i -> pts[i, j, idx], N_dims)
            val += gl_w[j] * eval_Φ(φ, ow, x, atlas)
        end
    end
    img[py, px] = val * r
end

# ── Observe ──────────────────────────────────────────────────────────────────

function observe(g::god, W::Int, H::Int; backend=GPU_BACKEND)
    points, owners, owner_list, r = cpu_points_and_owners(g, W, H)
    r < eps(T) && return fill(○, H, W)

    φs, atlas = build_φset(owner_list)

    img_gpu   = KernelAbstractions.zeros(backend, T, H, W)
    pts_gpu   = adapt(backend, points)
    own_gpu   = adapt(backend, owners)
    gl_gpu    = adapt(backend, Array(GL_WEIGHTS))
    atlas_gpu = adapt(backend, atlas)
    N_dims    = Int32(length(g.ẑero.d))

    Base.invokelatest() do
        κ_observe!(backend, 256)(img_gpu, pts_gpu, own_gpu, gl_gpu,
                                  T(r), N_dims, φs, atlas_gpu, ndrange=W*H)
    end
    KernelAbstractions.synchronize(backend)
    Array(img_gpu)
end

# ── CPU fallback ─────────────────────────────────────────────────────────────

function observe_cpu(g::god, W::Int, H::Int)
    points, owners, owner_list, r = cpu_points_and_owners(g, W, H)
    r < eps(T) && return fill(○, H, W)

    φs, atlas = build_φset(owner_list)
    N_dims = length(g.ẑero.d)
    img = fill(○, H, W)

    Threads.@threads for idx in 1:(W * H)
        py, px = divrem(idx - 1, W) .+ (1, 1)
        val = zero(T)
        for j in 1:8
            ow = owners[j, idx]
            if iszero(ow)
                val += GL_WEIGHTS[j] * ○
            else
                x = ntuple(i -> points[i, j, idx], N_dims)
                val += GL_WEIGHTS[j] * eval_Φ(φs, ow, x, atlas)
            end
        end
        img[py, px] = val * r
    end
    img
end

# ── PPM ──────────────────────────────────────────────────────────────────────

function save_ppm(img::Matrix, path::String)
    H, W = size(img)
    open(path, "w") do io
        println(io, "P3\n$W $H\n255")
        for py in 1:H
            for px in 1:W
                c = round(Int, clamp(img[py, px], zero(eltype(img)), one(eltype(img))) * 255)
                print(io, c, ' ', c, ' ', c, ' ')
            end
            println(io)
        end
    end
end

function observe!(g::god, path::String, W::Int=512, H::Int=512; kwargs...)
    img = observe(g, W, H; kwargs...)
    save_ppm(img, path)
    img
end