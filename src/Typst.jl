const DPI = 300
# const DPI = 144
const TYPST_TEMPLATE(content) = """
#set page(width: auto, height: auto, margin: (top: 5pt, bottom: 5pt, left: 5pt, right: 5pt))
#set text(size: 10pt)
$content
"""

const TYPST_CACHE = Dict{UInt,AbstractMatrix{T}}()
# collect(values(TYPST_CACHE))[1]
struct ΦTypst{M}
    mat::M
    H::Int32
    W::Int32
end
function (f::ΦTypst)(x)
    px = unsafe_trunc(UInt32, x[2] * (f.W - 1)) + 1
    py = unsafe_trunc(UInt32, (1 - x[3]) * (f.H - 1)) + 1
    f.mat[py, px]
end
function Adapt.adapt_structure(to, f::ΦTypst)
    ΦTypst(Adapt.adapt(to, f.mat), f.H, f.W)
end
function Adapt.adapt_structure(to, Φ::ΦTuple)
    ΦTuple(map(ϕ -> Adapt.adapt(to, ϕ), Φ.ϕ))
end

function typst(code)
    mat = typst_to_matrix(code)
    H, W = size(mat)
    ΦTypst(mat, Int32(H), Int32(W))
end

function typst_to_matrix(typst_code)
    h = hash(typst_code)
    haskey(TYPST_CACHE, h) && return TYPST_CACHE[h]
    cmd = `typst compile - --format png --ppi $DPI -`
    rgba = pipeline(IOBuffer(TYPST_TEMPLATE(typst_code)), cmd) |> read |> IOBuffer |> PNGFiles.load
    mat = KernelAbstractions.allocate(GPU_BACKEND, T, size(rgba)...)
    copyto!(mat, rgb2c.(rgba))
    TYPST_CACHE[h] = mat
end

# typst_code="A"
# mat = typst_to_matrix("A")
# φ_hi = Φ_typst(mat)
# gpu_safe(φ_hi, 4)
