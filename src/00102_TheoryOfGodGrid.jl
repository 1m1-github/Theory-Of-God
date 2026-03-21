# ♯=g.♯
# ♯ = (0,3,3,3,0)
function ∃̇(♯::NTuple{N,Int}, ϵ::∃{N,T}, God::𝕋{T}) where {N,T<:Real}
    # ks = ceil.(Int, log2.(♯))
    # ks = ♯
    # d = @SVector (!).(iszero.(♯))
    ♯̂ = [!iszero(g) for g = ♯]
    ♯̇ = ♯[♯̂]
    Ṅ = length(♯̇)
    Ns = ntuple(i -> (1 << ♯̇[i]) + 1, Ṅ)
    grid = fill(○(T), Ns...)

    center = CartesianIndex(ntuple(i -> (Ns[i] + 1) >> 1, Ṅ))
    # μ = MVector(ϵ.μ)
    # i̇ = 1
    # for i = 1:N
    #     iszero(♯[i]) && continue
    #     μ[i] = ϵ.μ[i] + ϵ.ρ[i]*(2*T(center[i̇])/T(Ns[i̇])-1)
    #     i̇ +=1
    # end
    ρ = @SVector zeros(T,N)
    x = ∃{N,T}(ϵ, "", ϵ.d, ϵ.μ, ρ, ntuple(_->(true, true), N), ϵ.Φ)
    # grid[center] = ϵ.Φ(center, T[], CartesianIndex{N}[])
    _, grid[center], _ = ∃̇(x, ϵ, God, T[], CartesianIndex{N}[])

    for b in 0:(1<<N)-1
        ci = CartesianIndex(ntuple(d -> 1 + (Ns[d] - 1) * ((b >> (d - 1)) & 1), N))
        grid[ci] = ϵ.Φ(ci, [grid[center]], [center])
    end

    for ℓ in 1:maximum(ks)
        ss = ntuple(d -> ℓ ≤ ks[d] ? 1 << (ks[d] - ℓ) : 0, N)
        strides = ntuple(d -> max(ss[d], 1), N)
        ci = CartesianIndices(ntuple(d -> 1:strides[d]:Ns[d], N))

        @threads for pt in ci
            isnan(grid[pt]) || continue

            odd_mask = ntuple(d -> ss[d] > 0 && isodd((pt[d] - 1) ÷ ss[d]), N)
            nodd = count(odd_mask)
            nodd == 0 && continue

            bit_pos = MVector{N,Int}(undef)
            j = 0
            for d in 1:N
                if odd_mask[d]
                    bit_pos[d] = j
                    j += 1
                end
            end

            pvals = MVector{1 << N,T}(undef)
            pcoords = MVector{1 << N,CartesianIndex{N}}(undef)
            np = 0
            for b in 0:(1<<nodd)-1
                valid = true
                coords = MVector{N,Int}(undef)
                for d in 1:N
                    if odd_mask[d]
                        coords[d] = pt[d] + ss[d] * (2 * ((b >> bit_pos[d]) & 1) - 1)
                        if coords[d] < 1 || coords[d] > Ns[d]
                            valid = false
                            break
                        end
                    else
                        coords[d] = pt[d]
                    end
                end
                if valid
                    np += 1
                    pvals[np] = grid[coords...]
                    pcoords[np] = CartesianIndex(Tuple(coords))
                end
            end
            grid[pt] = ϵ.Φ(pt, view(pvals, 1:np), view(pcoords, 1:np))
        end
    end
    grid
end
