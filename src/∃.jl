const ○ = one(T) / (one(T) + one(T))
const ○̂ = x -> ○
abstract type ∀ end
struct ∃{N,F,P<:∀} <: ∀
    ϵ̂::P
    d::SVector{N,T}
    μ::SVector{N,T}
    ρ::SVector{N,T}
    ∂::SVector{N,Tuple{Bool,Bool}}
    Φ::F
    h::UInt
    function ∃(ϵ̂::∀, d::SVector{N,T}, μ::SVector{N,T}, ρ::SVector{N,T}, ∂::SVector{N,Tuple{Bool,Bool}}, Φ::F) where {N,F}
        @assert all(zero(T) .≤ d .≤ one(T))
        @assert all(zero(T) .≤ μ .- ρ .≤ μ .+ ρ .≤ one(T))
        @assert gpu_safe(Φ, N)
        p = sortperm(d)
        d, μ, ρ = map(x -> x[p], (d, μ, ρ))
        ∂ = SVector(ntuple(i -> ∂[p[i]], N))
        h = hash(d, hash(μ, hash(ρ, hash(∂, hash(ϵ̂)))))
        new{N,F,typeof(ϵ̂)}(ϵ̂, d, μ, ρ, ∂, Φ, h) # todo wrap Φ in Φ̂ for global to local map and boundary check
    end
end
# function Φ̂(Φ)
#     xΩ -> begin
#     end
# end
Base.hash(ϵ::∃, h::UInt) = hash(ϵ.h, h)
struct 𝕋 <: ∀
    ϵ̃::Dict{∀,Vector{∃}}
    Ο::Dict{∀,Int}
    L::ReentrantLock
    s::Ref{Int}
    function 𝕋()
        ϵ̃ = Dict{∀,Vector{∃}}()
        Ο = Dict{∀,Int}()
        Ω = new(ϵ̃, Ο, ReentrantLock(), Ref(1))
        Ω.ϵ̃[Ω] = ∃[]
        Ω.Ο[Ω] = Ω.s[]
        Ω
    end
end
Base.hash(::𝕋, h::UInt) = hash(:Ω, h)
t(Ο::Int) = one(T) - one(T) / (one(T) + T(log(Ο)))
t(ω::∀) = t(ω.Ο[ω])
# δ(ϵ, ϵ)
# function δ(ϵ::∃, ϵ̂::∃)
#     nϵ̂ = length(ϵ.ϵ̂.d)
#     SVector(ntuple(length(ϵ.d)) do i
#         dᵢ = ϵ.d[i]
#         iϵ̂ = searchsortedfirst(ϵ.ϵ̂.d, dᵢ)
#         iϵ̂ ≤ nϵ̂ && !iszero(ϵ.ϵ̂.ρ[iϵ̂]) && return ϵ̂.d[iϵ̂]
#         dᵢ
#     end)
# end
# function Base.copy!(ϵ::∃, ϵ̂::∃, Ω::𝕋)
#     ϵϵ̃ = Ω.ϵ̃[ϵ]
#     d = δ(ϵ, ϵ̂)
#     ϵ̃ = ∃(ϵ̂, d, ϵ.μ, ϵ.ρ, ϵ.∂, ϵ.Φ)
#     ∃!(ϵ̃, Ω, ϵ̂)
#     for ϵ̃̃ = ϵϵ̃
#         copy!(ϵ̃̃, ϵ̃, Ω)
#     end
#     ϵ̃
# end
ρ̂(ρ, ρ̂) = 2 .* ρ̂ .* ρ
ρ̂(ϵ::∃) = ρ̂(ϵ.ρ, ϵ.ϵ̂.ρ)
μ̂(μ, μ̂, ρ̂) = μ̂ .- ρ̂ .+ 2 .* ρ̂ .* μ
μ̂(ϵ::∃) = μ̂(ϵ.μ, ϵ.ϵ̂.μ, ϵ.ϵ̂.ρ)
μ̃(ϵ::∃, μ) = (μ .- (ϵ.μ .- ϵ.ρ)) ./ ϵ.ρ ./ 2
ρ̃(ϵ::∃, ρ) = ρ ./ ϵ.ρ ./ 2
μρΩ(ϵ::∃) = begin
    ϵ.ϵ̂ isa 𝕋 && return ϵ.μ, ϵ.ρ
    μ̂, ρ̂ = μρΩ(ϵ.ϵ̂)
    μ = μ̂ .- ρ̂ .+ 2 .* ρ̂ .* ϵ.μ
    ρ = 2 .* ρ̂ .* ϵ.ρ
    μ, ρ
end
function μρ(ϵ::∃, d)
    i = searchsortedfirst(ϵ.d, d)
    N = length(ϵ.d)
    if i ≤ N && ϵ.d[i] == d
        return ϵ.μ[i], ϵ.ρ[i], ϵ.∂[i]
    end
    d₀ = d₁ = ϵ.d[1]
    dₙ = ϵ.d[N]
    μ₀ = μ₁ = ○
    if d < d₀
        d₀ = zero(T)
        μ₁ = ϵ.μ[1]
    elseif dₙ < d
        d₀, d₁ = dₙ, one(T)
        μ₀ = ϵ.μ[N]
    else
        d₀, d₁ = ϵ.d[i-1], ϵ.d[i]
        μ₀, μ₁ = ϵ.μ[i-1], ϵ.μ[i]
    end
    d = (d - d₀) / (d₁ - d₀)
    μ₀ + (μ₁ - μ₀) * d, zero(T), (true, true)
end
function ∂(x::∃, ::𝕋)
    zeroₓ = x.μ .- x.ρ
    any(==(zero(T)), zeroₓ) && return true
    oneₓ = x.μ .+ x.ρ
    any(==(one(T)), oneₓ)
end
function ∂(x::∃, ϵ::∃)
    zeroμ, oneμ = ϵ.μ .- ϵ.ρ, ϵ.μ .+ ϵ.ρ
    Threads.@threads for i = eachindex(ϵ.d)
        d = ϵ.d[i]
        iszero(ϵ.ρ[i]) && continue
        μₓ, _ = μρ(x, d)
        (μₓ == zeroμ[i] || μₓ == oneμ[i]) && return true
    end
    false
end
function Base.:(⊆)(zero₁, one₁, ∂₁::Tuple{Bool,Bool}, zero₂, one₂, ∂₂::Tuple{Bool,Bool})
    żero = zero₂ < zero₁ || (zero₂ == zero₁ && (!∂₁[1] || ∂₂[1]))
    ȯne = one₁ < one₂ || (one₁ == one₂ && (!∂₁[2] || ∂₂[2]))
    żero && ȯne
end
function ⪽(ϵ₁::∃, ϵ₂::∃)
    found = Threads.Atomic{Bool}(false)
    result = Threads.Atomic{Bool}(true)
    Threads.@threads for i₂ in eachindex(ϵ₂.d)
        d₂ = ϵ₂.d[i₂]
        ρ₂ = ϵ₂.ρ[i₂]
        μ₂ = ϵ₂.μ[i₂]
        if iszero(ρ₂)
            Threads.atomic_or!(found, true)
            μ₁, _, _ = μρ(ϵ₁, d₂)
            if μ₁ != μ₂
                Threads.atomic_and!(result, false)
            end
            continue
        end
        Threads.atomic_or!(found, true)
        result[] || continue  # skip work if already false
        μ₁, ρ₁, ∂₁ = μρ(ϵ₁, d₂)
        zero₁, one₁ = μ₁ - ρ₁, μ₁ + ρ₁
        zero₂, one₂ = μ₂ - ρ₂, μ₂ + ρ₂
        if !⊆(zero₁, one₁, ∂₁, zero₂, one₂, ϵ₂.∂[i₂])
            Threads.atomic_and!(result, false)
        end
    end
    found[] && result[]
end
function α(ϵ)
    αϵ = Set{∃}([ϵ])
    p = ϵ.ϵ̂
    while p isa ∃
        push!(αϵ, p)
        p = p.ϵ̂
    end
    αϵ
end
function α(ϵ₁::∃, ϵ₂::∃, ω)
    αϵ₁ = α(ϵ₁)
    for ϵ = αϵ₁
        ϵ == ϵ₂ && return ϵ₂
    end
    # ϵ₂ ∈ αϵ₁ && return ϵ₂
    ϵ̂ = ϵ₂.ϵ̂
    while ϵ̂ isa ∃
        for ϵ = αϵ₁
            ϵ == ϵ̂ && return ϵ̂
        end
        # ϵ̂ ∈ αϵ₁ && return ϵ̂
        ϵ̂ = ϵ̂.ϵ̂
    end
    ω
end
function ℼ(ϵ::∃)
    ϵ̂ = ϵ.ϵ̂
    ϵ̂ isa 𝕋 && return ϵ
    ϵ̂̂ = ∃(ϵ̂.ϵ̂, ϵ.d, μ̂(ϵ), ρ̂(ϵ), ϵ.∂, ϵ.Φ)
    ℼ(ϵ̂̂)
end
ℼ(ω::𝕋, ::Any) = ω
ℼ(ϵ, ::𝕋) = ℼ(ϵ)
function ℼ(ϵ₁::∃, ϵ₂::∃)
    ○̂̂ = SVector(ntuple(○̂, length(ϵ₁.d)))
    ϵ₁ == ϵ₂ && return ∃(ϵ₂, ϵ₁.d, ○̂̂, ○̂̂, ϵ₁.∂, ϵ₁.Φ)
    ϵ₁.ϵ̂ == ϵ₂ && return ϵ₁
    if ϵ₂ ∈ α(ϵ₁)
        return ℼ(∃(ϵ₁.ϵ̂.ϵ̂, ϵ₁.d, μ̂(ϵ₁), ρ̂(ϵ₁), ϵ₁.∂, ϵ₁.Φ), ϵ₂)
    end
    ϵ₁Ω = ℼ(ϵ₁)
    ϵ₂Ω = ℼ(ϵ₂)
    μ = μ̃(ϵ₂Ω, ϵ₁Ω.μ)
    ρ = ρ̃(ϵ₂Ω, ϵ₁Ω.ρ)
    ∃(ϵ₂, ϵ₁.d, μ, ρ, ϵ₁.∂, ϵ₁.Φ)
end
⫉(ϵ, ::𝕋) = true
function ⫉(ϵ₁::∃, ϵ₂::∃, ω)
    ϵ₁.ϵ̂ == ϵ₂.ϵ̂ && return ϵ₁ ⪽ ϵ₂
    ϵ̂ = α(ϵ₁, ϵ₂, ω)
    ℼ(ϵ₁, ϵ̂) ⪽ ℼ(ϵ₂, ϵ̂)
end
Base.:(==)(ϵ₁::∃, ϵ₂::∃) = ϵ₁.d == ϵ₂.d && ϵ₁.μ == ϵ₂.μ && ϵ₁.ρ == ϵ₂.ρ && ϵ₁.∂ == ϵ₂.∂
Base.:(==)(::∃, ::𝕋) = false
function β(ϵ₁::∃, ϵ₂::∀, ω)
    ϵ̃ = ω.ϵ̃[ϵ₂]
    ϵ̃₂ = filter(ϵ -> ⫉(ϵ₁, ϵ, ω), ϵ̃)
    isempty(ϵ̃₂) && return ϵ₂
    β(ϵ₁, only(ϵ̃₂), ω)
end
function Base.:∩(zero₁, one₁, ∂₁::Tuple{Bool,Bool}, zero₂, one₂, ∂₂::Tuple{Bool,Bool})
    żero = max(zero₁, zero₂)
    ȯne = min(one₁, one₂)
    żero < ȯne && return true
    żero ≠ ȯne && return false
    ∂₀₀ = zero₂ < zero₁ ? ∂₁[1] : (zero₁ < zero₂ ? ∂₂[1] : ∂₁[1] && ∂₂[1])
    ∂₀₁ = one₁ < one₂ ? ∂₁[2] : (one₂ < one₁ ? ∂₂[2] : ∂₁[2] && ∂₂[2])
    ∂₀₀ && ∂₀₁
end
Base.:∩(ϵ₁::∃, ::𝕋) = true
function Base.:∩(ϵ₁::∃, ϵ₂::∃, ω)
    ϵ₁ == ϵ₂ && return true
    if ϵ₁.ϵ̂ !== ϵ₂.ϵ̂
        ϵ̂ = α(ϵ₁, ϵ₂, ω)
        return ∩(ℼ(ϵ₁, ϵ̂), ℼ(ϵ₂, ϵ̂), ω)
    end
    d̂ = sort(ϵ₁.d ∪ ϵ₂.d)
    if !iszero(d̂[1])
        if !isone(d̂[end])
            d̂ = [zero(T), d̂..., one(T)]
        else
            d̂ = [zero(T), d̂...]
        end
    elseif !isone(d̂[end])
        d̂ = [d̂..., one(T)]
    end
    μ₁, ρ₁, ∂₁ = μρ(ϵ₁, zero(T))
    μ₂, ρ₂, ∂₂ = μρ(ϵ₂, zero(T))
    μ₁prev, μ₂prev = μ₁, μ₂
    for (i, d) = enumerate(d̂)
        if 1 < i
            μ₁, ρ₁, ∂₁ = μρ(ϵ₁, d)
            μ₂, ρ₂, ∂₂ = μρ(ϵ₂, d)
        end
        zero₁, one₁ = μ₁ - ρ₁, μ₁ + ρ₁
        zero₂, one₂ = μ₂ - ρ₂, μ₂ + ρ₂
        !∩(zero₁, one₁, ∂₁, zero₂, one₂, ∂₂) && return false
        i == 1 && continue
        μ₁prev < μ₂prev && μ₂ < μ₁ && return true
        μ₂prev < μ₁prev && μ₁ < μ₂ && return true
        μ₁prev, μ₂prev = μ₁, μ₂
    end
    ϵ̃ = get(ω.ϵ̃, ϵ₂, ∃[])
    isempty(ϵ̃) && return true
    ϵ̂ = ℼ(ϵ₁, ϵ₂)
    all(ϵ̃ -> ∩(ϵ̂, ϵ̃, ω), ϵ̃)
end
function ∃!(ϵ::∃, ω)
    lock(ω.L)
    ϵ̂ = β(ϵ, ω, ω)
    ϵ̃ = ω.ϵ̃[ϵ̂]
    any(ϵ̃ -> ∩(ϵ, ϵ̃, ω), ϵ̃) && (unlock(ω.L); return nothing)
    if ϵ̂ !== ϵ.ϵ̂
        ϵ = ∃(ϵ̂, ϵ.d, ϵ.μ, ϵ.ρ, ϵ.∂, ϵ.Φ)
    end
    while Sys.free_memory() < ω.s[] + sizeof(ϵ)
        rm!(ω)
    end
    ω.s[] += sizeof(ϵ)
    ω.Ο[ω] += 1
    ω.Ο[ϵ] = ω.Ο[ω]
    push!(ω.ϵ̃[ϵ̂], ϵ)
    ω.ϵ̃[ϵ] = ∃[]
    unlock(ω.L)
    ϵ
end
function rm!(ω::𝕋)
    ϵ̂̂ = argmin(ϵ -> ω.Ο[ϵ], filter(k -> k isa ∃, keys(ω.Ο)))
    ω.s[] -= sizeof(ϵ̂̂)
    filter!(ϵ -> ϵ !== ϵ̂̂, ω.ϵ̃[ϵ̂̂.ϵ̂])
    delete!(ω.ϵ̃, ϵ̂̂)
    delete!(ω.Ο, ϵ̂̂)
end
# function Base.:(-)(ϵ₁::∃, ϵ₂::∃, ω)
#     d̂ = sort!(ϵ₂.d ∪ ϵ₁.d)
#     N = length(d̂)
#     μ = MVector{N,T}(undef)
#     ρ = MVector{N,T}(undef)
#     ∂out = MVector{N,Tuple{Bool,Bool}}(undef)
#     Threads.@threads for i in eachindex(d̂)
#         d = d̂[i]
#         ϵ₂μ, ϵ₂ρ, ϵ₂∂ = μρ(ϵ₂, d)
#         ϵ₁μ, ϵ₁ρ, ϵ₁∂ = μρ(ϵ₁, d)
#         żero = ϵ₂μ - ϵ₂ρ
#         ȯne = ϵ₁μ + ϵ₁ρ
#         ρ[i] = abs(ȯne - żero) / 2
#         μ[i] = żero + ρ[i]
#         ∂out[i] = (ϵ₂∂[1], ϵ₁∂[2])
#     end
#     ϵ̂ = α(ϵ₁, ϵ₂, ω)
#     ∃(ϵ̂, SVector{N}(d̂), SVector{N}(μ), SVector{N}(ρ), SVector{N}(∂out), ϵ₁.Φ)
# end

function gpu_safe(Φ, N)
    try
        @kernel gpu(Φ, x) = Φ(x)
        x = KernelAbstractions.zeros(GPU_BACKEND, T, N)
        gpu(GPU_BACKEND, GPU_BACKEND_WORKGROUPSIZE)(Φ, x, ndrange=1)
        true
    catch e
        bt = catch_backtrace()
        showerror(stderr, e, bt)
        false
    end
end

function √(ϵ::∃)
    n = 0
    p = ϵ
    while !(p.ϵ̂ isa 𝕋)
        p = p.ϵ̂
        n += 1
    end
    p, n
end
# x=ϵ
function X(x::∃, ∇, ω)
    ϵ = β(x, ω, ω)
    ϵ === ω && return ω, true
    ∂(x, ϵ) && return ω, true
    ∩(x, ϵ, ω) && return ϵ, true # ?
    _, n = √(ϵ)
    ∇ < n && return ω, false
    ϵ̃ = ω.ϵ̃[ϵ]
    Threads.@threads for ϵ̃ = filter(ϵ̃ -> x ⫉ ϵ̃, ϵ̃)
        ∩(x, ϵ̃, ω) && return ϵ̃, true
        ϵ̂, found = X(x, ∇)
        found && return ϵ̂, true
    end
    ω, false
end
