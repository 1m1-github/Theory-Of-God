using StaticArrays

"""
Tighten an interval [lo, hi] given constraint L ≤ s*coeff + s_other*coeff_other ≤ U,
relaxing over s_other ∈ [-1,1]. Returns (lo, hi, empty).
"""
@inline function tighten(lo::Float64, hi::Float64, coeff::Float64, abs_other::Float64, L::Float64, U::Float64)
    if coeff > 0
        lo = max(lo, (L - abs_other) / coeff)
        hi = min(hi, (U + abs_other) / coeff)
    elseif coeff < 0
        lo = max(lo, (U + abs_other) / coeff)
        hi = min(hi, (L - abs_other) / coeff)
    else
        (L > abs_other || U < -abs_other) && return (lo, hi, true)
    end
    return (lo, hi, lo > hi)
end

"""
Convert continuous s ∈ [-1,1] bounds to 1-based grid indices on an axis of size ns.
s_i = (2i - 1 - ns) / (ns - 1)  ⟹  i = (s*(ns-1) + ns + 1) / 2
Returns (ilo, ihi) or (0, 0) if empty.
"""
@inline function s_to_indices(slo::Float64, shi::Float64, ns::Int)
    h = (ns - 1) * 0.5
    mid = (ns + 1) * 0.5
    ilo = clamp(ceil( Int, slo * h + mid), 1, ns)
    ihi = clamp(floor(Int, shi * h + mid), 1, ns)
    ilo > ihi && return (0, 0)
    return (ilo, ihi)
end

"""
    pyramid_box_intersection(z, o, ex, ey, wx, wy, a, b, nx, ny, nz)

Per pyramid slice k=1:nz, find (i_min, i_max, j_min, j_max) index ranges of grid
points inside axis-aligned box [a,b]. Empty slices → (0,0,0,0). Cost: O(nz·n).
"""
# i=i
# ĩ=length(ΦΦ)
# z=g.ẑero.μ
# o=g.ône.μ
# ex=
# ey=
# wx=
# wy=
# a=ϵ.μ .- ϵ.ρ
# b=ϵ.μ .+ ϵ.ρ
# nx=g.♯[1]
# ny=g.♯[2]
# nz=GL_N
function pyramid_box_intersection(
    i,ĩ,
    z::SVector{N,T}, o::SVector{N,T},
    ex,ey,wx,wy,a,b,nx,ny,nz
    # ex::SVector{N}, ey::SVector{N},
    # wx::Float64, wy::Float64,
    # a::SVector{N}, b::SVector{N},
    # nx::Int, ny::Int, nz::Int
) where {N}
    c = (z + o) * 0.5
    doc = o - c

    # out = Matrix{Int}(undef, 4, nz)
    intersects = false

    @inbounds for k in 1:nz
        t = k / (nz + 1)
        omt = 1 - t
        wxk = wx * omt
        wyk = wy * omt

        si_lo, si_hi = -1.0, 1.0
        sj_lo, sj_hi = -1.0, 1.0
        empty = false

        for m in 1:N
            ck = c[m] + t * doc[m]
            α = wxk * ex[m]
            β = wyk * ey[m]
            L = a[m] - ck
            U = b[m] - ck

            si_lo, si_hi, empty = tighten(si_lo, si_hi, α, abs(β), L, U)
            empty && break
            sj_lo, sj_hi, empty = tighten(sj_lo, sj_hi, β, abs(α), L, U)
            empty && break
        end

        if empty
            # out[1,k] = 0; out[2,k] = 0; out[3,k] = 0; out[4,k] = 0
        else
            ilo, ihi = s_to_indices(si_lo, si_hi, nx)
            jlo, jhi = s_to_indices(sj_lo, sj_hi, ny)
            if ilo == 0 || jlo == 0
                # out[1,k] = 0; out[2,k] = 0; out[3,k] = 0; out[4,k] = 0
            else
                # out[1,k] = ilo; out[2,k] = ihi; out[3,k] = jlo; out[4,k] = jhi
                i[ilo:ihi,jlo:jhi,k] .= ĩ
                intersects = true
            end
        end
    end

    intersects
end

# wx = W
# wy = H
# nx, ny, nz = 2000, 2000, 8
# z = SVector(0.0, 0.0, 0.0)
# o = SVector(1.0, 1.0, 1.0)
# a = SVector(0.75, 0.75, 0.75)
# b = SVector(1.0, 1.0, 1.0)
# @time out = pyramid_box_intersection(z, o, ex, ey, wx, wy, a, b, nx, ny, nz)

