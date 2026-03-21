using PNGFiles

# ── Typst Φ: text → pixel matrix → lookup ───────────────────────────────────

const DPI = 144
const TYPST_TEMPLATE(content) = """
#set page(width: auto, height: auto, margin: (top: 5pt, bottom: 5pt, left: 5pt, right: 5pt))
#set text(size: 20pt)
$content
"""

const TYPST_CACHE = Dict{UInt,Matrix{T}}()

function typst_to_matrix(typst_code::String)
    h = hash(typst_code)
    haskey(TYPST_CACHE, h) && return TYPST_CACHE[h]
    cmd = `typst compile - --format png --ppi $DPI -`
    rgba = pipeline(IOBuffer(TYPST_TEMPLATE(typst_code)), cmd) |> read |> IOBuffer |> PNGFiles.load
    # convert to T ∈ [0,1]: use luminance, invert (black text → high Φ)
    # mat = Matrix{T}(one(T) .- T.(Float64.(red.(rgba))))
    mat = Matrix{T}(one(T) .- T.(getfield.(rgba, :r)))
    TYPST_CACHE[h] = mat
    mat
end

function Φ_typst(mat::Matrix{T})
    H, W = size(mat)
    (x, y) -> begin
        px = clamp(round(Int, x * (W - 1)) + 1, 1, W)
        py = clamp(round(Int, y * (H - 1)) + 1, 1, H)
        mat[py, px]
    end
end

struct ΦTex
    tex::AbstractMatrix{T}
    H::Int32
    W::Int32
end

@inline function (Φ::ΦTex)(x, y)
    px = clamp(unsafe_trunc(Int32, x * Φ.W) + Int32(1), Int32(1), Φ.W)
    py = clamp(unsafe_trunc(Int32, y * Φ.H) + Int32(1), Int32(1), Φ.H)
    Φ.tex[py, px]
end
Φ_typst(code::String, backend=GPU_BACKEND) = 
    ΦTex(adapt(backend, typst_to_matrix(code))..., Int32.(size(typst_to_matrix(code)))...)
function Φ_typst(code::String, backend=GPU_BACKEND)
    mat = typst_to_matrix(code)
    H, W = size(mat)
    ΦTex(adapt(backend, mat), Int32(H), Int32(W))
end

struct ΦTex
    tex::Matrix{T}  # CPU matrix — always safe
    H::Int32
    W::Int32
end

@inline function (Φ::ΦTex)(x, y)
    px = clamp(unsafe_trunc(Int32, x * Φ.W) + Int32(1), Int32(1), Φ.W)
    py = clamp(unsafe_trunc(Int32, y * Φ.H) + Int32(1), Int32(1), Φ.H)
    Φ.tex[py, px]
end

# function Φ_typst(code::String)
#     mat = typst_to_matrix(code)
#     H, W = size(mat)
#     ΦTex(mat, Int32(H), Int32(W))
# end
function Φ_typst(code::String; offset::Int32=Int32(0))
    mat = typst_to_matrix(code)
    H, W = size(mat)
    ΦTex3(Int32(H), Int32(W), offset), mat
end