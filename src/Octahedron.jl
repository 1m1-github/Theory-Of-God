"""
Tighten an interval [lo, hi] given constraint L ≤ s*coeff + s_other*coeff_other ≤ U,
relaxing over s_other ∈ [-1,1]. Returns (lo, hi, empty).
"""
@inline function tighten(lo, hi, coeff, abs_other, L, U)
    if zero(T) < coeff
        lo = max(lo, (L - abs_other) / coeff)
        hi = min(hi, (U + abs_other) / coeff)
    elseif coeff < zero(T)
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
@inline function s_to_indices(slo, shi, ns)
    h = (ns - 1) * ○
    mid = (ns + 1) * ○
    ilo = clamp(ceil(Int, slo * h + mid), 1, ns)
    ihi = clamp(floor(Int, shi * h + mid), 1, ns)
    ilo > ihi && return (0, 0)
    return (ilo, ihi)
end

"""
    pyramid_box_intersection(z, o, ex, ey, rx, ry, a, b, nx, ny, nz)

Per pyramid slice k=1:nz, find (i_min, i_max, j_min, j_max) index ranges of grid
points inside axis-aligned box [a,b]. Empty slices → (0,0,0,0). Cost: O(nz·n).
"""
# i=i
# ĩ=length(ΦΦ)
# z=g.ẑero.μ
# o=g.f̂ocus.μ
# ex=
# ey=
# rx=
# ry=
# a=ϵ.μ .- ϵ.ρ
# b=ϵ.μ .+ ϵ.ρ
# nx=g.♯[1]
# ny=g.♯[2]
# nz=GL_N
function pyramid_box_intersection!(
    i, ĩ, z, o, dx, dy, a, b, nx, ny, nz
)
    N = length(z)
    c = (z + o) * ○
    d = o - c
    intersects = false
    @inbounds for k = 1:nz
        t = GL_NODES[k]
        si_lo, si_hi = -one(T), one(T)
        sj_lo, sj_hi = -one(T), one(T)
        empty = false
        for m = 1:N
            ck = c[m] + t * d[m]
            α = dx[m] * ○
            β = dy[m] * ○
            L = a[m] - ck
            U = b[m] - ck
            si_lo, si_hi, empty = tighten(si_lo, si_hi, α, abs(β), L, U)
            empty && break
            sj_lo, sj_hi, empty = tighten(sj_lo, sj_hi, β, abs(α), L, U)
            empty && break
        end
        empty && continue
        ilo, ihi = s_to_indices(si_lo, si_hi, nx)
        jlo, jhi = s_to_indices(sj_lo, sj_hi, ny)
        (iszero(ilo) || iszero(jlo)) && continue
        i[ilo:ihi, jlo:jhi, k] .= ĩ
        intersects = true
    end
    intersects
end

# rx = W
# ry = H
# nx, ny, nz = 2000, 2000, 8
# z = SVector(0.0, 0.0, 0.0)
# o = SVector(1.0, 1.0, 1.0)
# a = SVector(0.75, 0.75, 0.75)
# b = SVector(1.0, 1.0, 1.0)
# @time out = pyramid_box_intersection(z, o, ex, ey, rx, ry, a, b, nx, ny, nz)

