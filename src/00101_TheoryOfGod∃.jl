using KernelAbstractions, StaticArrays

const ○ = one(T) / (one(T) + one(T))
# const ○̂ = _ -> ○
const ○̂ = (x...) -> ○
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
        @assert 1 ≤ N
        @assert all(zero(T) .≤ d .≤ one(T))
        @assert all(zero(T) .≤ μ .≤ one(T))
        @assert all(zero(T) .≤ ρ .≤ one(T))
        # @assert gpu_safe(Φ, N)
        p = sortperm(d)
        d, μ, ρ = map(x -> x[p], (d, μ, ρ))
        ∂ = SVector(ntuple(i -> ∂[p[i]], N))
        h = hash(d, hash(μ, hash(ρ, hash(∂, hash(ϵ̂)))))
        new{N,F,typeof(ϵ̂)}(ϵ̂, d, μ, ρ, ∂, Φ, h)
    end
end
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
t(ω::∀=Ω) = t(Ω.Ο[ω])
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
# μ(::𝕋) = 
# μ(ϵ::∃) = ϵ.μ
ρ̂(ρ, ρ̂) = 2 .* ρ̂ .* ρ
ρ̂(ϵ::∃) = ρ̂(ϵ.ρ, ϵ.ϵ̂.ρ)
# ρ̂(ϵ::∃) = 2 .* ϵ.ϵ̂.ρ .* ϵ.ρ
μ̂(μ, μ̂, ρ̂) = μ̂ .- ρ̂ .+ 2 .* ρ̂ .* μ
μ̂(ϵ::∃) = μ̂(ϵ.μ, ϵ.ϵ̂.μ, ϵ.ϵ̂.ρ)
# μ̂(ϵ::∃) = ϵ.ϵ̂.μ .- ϵ.ϵ̂.ρ .+ 2 .* ϵ.ϵ̂.ρ .* ϵ.μ
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
    # for (i, d) = enumerate(ϵ.d)
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
# (i₂, d₂) = collect(enumerate(ϵ₂.d))[1]
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
function α(ϵ₁::∃, ϵ₂::∃, ω=Ω)
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
ℼ(Ω::𝕋, ::Any) = Ω
ℼ(ϵ, ::𝕋) = ℼ(ϵ)
function ℼ(ϵ₁::∃, ϵ₂::∃)
    ○̂̂ = SVector(ntuple(○̂, length(ϵ₁.d)))
    ϵ₁ == ϵ₂ && return ∃(ϵ₂, ϵ₁.d, ○̂̂, ○̂̂, ϵ₁.∂, ϵ₁.Φ)
    ϵ₁.ϵ̂ == ϵ₂ && return ϵ₁
    if ϵ₂ ∈ α(ϵ₁)
        # @show typeof.([ϵ₁, ϵ₂])
        # @show ϵ₂ isa 𝕋
        # @show ϵ₂ === ϵ₁
        # @show ϵ₂ == ϵ₁
        # @show length(α(ϵ₁)), first(α(ϵ₁)) === Ω, typeof(first(α(ϵ₁))), first(α(ϵ₁))==ϵ₁, first(α(ϵ₁))==ϵ₂
        # @show typeof(ϵ₁.ϵ̂)
        # ϵ₂ isa 𝕋 && return ℼ(∃(ϵ₁.ϵ̂, ϵ₁.d, μ̂(ϵ₁), ρ̂(ϵ₁), ϵ₁.∂, ϵ₁.Φ))
        return ℼ(∃(ϵ₁.ϵ̂.ϵ̂, ϵ₁.d, μ̂(ϵ₁), ρ̂(ϵ₁), ϵ₁.∂, ϵ₁.Φ), ϵ₂)
    end
    ϵ₁Ω = ℼ(ϵ₁)
    ϵ₂Ω = ℼ(ϵ₂)
    μ = μ̃(ϵ₂Ω, ϵ₁Ω.μ)
    ρ = ρ̃(ϵ₂Ω, ϵ₁Ω.ρ)
    ∃(ϵ₂, ϵ₁.d, μ, ρ, ϵ₁.∂, ϵ₁.Φ)
end
⫉(ϵ, ::𝕋) = true
# ϵ₁ , ϵ₂ = ϵ₁,ϵ̃
function ⫉(ϵ₁::∃, ϵ₂::∃, ω=Ω)
    ϵ₁.ϵ̂ == ϵ₂.ϵ̂ && return ϵ₁ ⪽ ϵ₂
    ϵ̂ = α(ϵ₁, ϵ₂, ω)
    ℼ(ϵ₁, ϵ̂) ⪽ ℼ(ϵ₂, ϵ̂)
end
# ϵ₁=x
# ϵ₁=ϵ
# ϵ₁=-(g.ône, g.ẑero,ω)
# ϵ₁=view
# ϵ̂ = β(ϵ, Ω)
# ϵ₂=Ω
# ϵ₂=ω
# ϵ̃ = filter(ϵ̃ -> ϵ̃ ≠ ϵ₁, ϵ̃)[1]
# ϵ₁ ⫉ ϵ̃
# Ω.ϵ̃[ϵ₁]
# hash(ϵ)
# hash(ϵ̃[1])
# ϵ₁==ϵ̃[1]
# ϵ₂=ϵ̃[1]
Base.:(==)(ϵ₁::∃, ϵ₂::∃) = ϵ₁.d == ϵ₂.d && ϵ₁.μ == ϵ₂.μ && ϵ₁.ρ == ϵ₂.ρ && ϵ₁.∂ == ϵ₂.∂
Base.:(==)(::∃, ::𝕋) = false
function β(ϵ₁::∃, ϵ₂::∀, ω=Ω)
    ϵ̃ = ω.ϵ̃[ϵ₂]
    # ϵ̃₂ = filter(ϵ -> ϵ ≠ ϵ₁ && ⫉(ϵ₁, ϵ, ω), ϵ̃)
    ϵ̃₂ = filter(ϵ -> ⫉(ϵ₁, ϵ, ω), ϵ̃)
    isempty(ϵ̃₂) && return ω
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
# ϵ₁, ϵ₂=ϵ,ϵ̂
function Base.:∩(ϵ₁::∃, ϵ₂::∃, ω=Ω)
    ϵ₁ == ϵ₂ && return true
    if ϵ₁.ϵ̂ !== ϵ₂.ϵ̂
        ϵ̂ = α(ϵ₁, ϵ₂, ω)
        return ∩(ℼ(ϵ₁, ϵ̂),  ℼ(ϵ₂, ϵ̂), ω)
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
    all(ϵ̃ -> ∩(ϵ̂, ϵ̃,ω), ϵ̃)
end
# function √(ϵ::∃)
#     n = 0
#     p = ϵ
#     while !(p.ϵ̂ isa 𝕋)
#         p = p.ϵ̂
#         n += 1
#     end
#     p, n
# end
# # x=xϵ
# # ϵ=β(x, Ω, Ω)
# # ϵ.Φ(1)
# function X(x::∃, ∇)
#     ϵ = β(x, Ω)
#     ϵ === Ω && return Ω, true
#     ∂(x, ϵ) && return Ω, true
#     x ∩ ϵ && return ϵ, true # ?
#     _, n = √(ϵ)
#     ∇ < n && return Ω, false
#     ϵ̃ = Ω.ϵ̃[ϵ]
#     Threads.@threads for ϵ̃ = filter(ϵ̃ -> x ⫉ ϵ̃, ϵ̃)
#         x ∩ ϵ̃ && return ϵ̃, true
#         ϵ̂, found = X(x, ∇)
#         found && return ϵ̂, true
#     end
#     Ω, false
# end
# # x, ϵ, ∇ = xϵ, ϵ̂, typemax(UInt)
# # x, ϵ, ∇ = x, only(ϵ̃), ∇-1
# Φ(::Any, ::𝕋) = ○
# # Φ(x::∃, ϵ::∃) = ϵ.Φ(x.μ)
# Φ(x::∃, ϵ::∃) = ϵ.Φ(μ̃(ϵ, x.μ))
# function Φ(x::∃, ϵ::∀, ∇)
#     iszero(∇) && return Φ(x, ϵ)
#     ∂(x, ϵ) && return ○
#     ϵ̃ = filter(ϵ̃ -> x ⫉ ϵ̃, Ω.ϵ̃[ϵ])
#     isempty(ϵ̃) && return Φ(x, ϵ)
#     Φ(x, only(ϵ̃), ∇ - 1)
# end
# ϵ=ϵ̂
function ∃!(ϵ::∃, ω=Ω)
    lock(ω.L)
    ϵ̂ = β(ϵ, ω, ω)
    ϵ̃ = ω.ϵ̃[ϵ̂]
    any(ϵ̃ -> ∩(ϵ, ϵ̃,ω), ϵ̃) && (unlock(ω.L); return nothing)
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
# ϵ₁, ϵ₂ = ône, ẑero
# ϵ₁, ϵ₂ = g.ône, g.ẑero
# (i, d) = collect(enumerate(d̂))[1]
function Base.:(-)(ϵ₁::∃, ϵ₂::∃, ω=Ω)
    d̂ = sort!(ϵ₂.d ∪ ϵ₁.d)
    N = length(d̂)
    μ = MVector{N,T}(undef)
    ρ = MVector{N,T}(undef)
    ∂out = MVector{N,Tuple{Bool,Bool}}(undef)
    Threads.@threads for i in eachindex(d̂)
        d = d̂[i]
        ϵ₂μ, ϵ₂ρ, ϵ₂∂ = μρ(ϵ₂, d)
        ϵ₁μ, ϵ₁ρ, ϵ₁∂ = μρ(ϵ₁, d)
        żero = ϵ₂μ - ϵ₂ρ
        ȯne = ϵ₁μ + ϵ₁ρ
        ρ[i] = abs(ȯne - żero) / 2
        μ[i] = żero + ρ[i]
        ∂out[i] = (ϵ₂∂[1], ϵ₁∂[2])
    end
    ϵ̂ = α(ϵ₁, ϵ₂, ω)
    ∃(ϵ̂, SVector{N}(d̂), SVector{N}(μ), SVector{N}(ρ), SVector{N}(∂out), ϵ₁.Φ)
end

function gpu_safe(Φ, N)
    try
        @kernel gpu(Φ, x) = Φ(x)
        x = KernelAbstractions.zeros(GPU_BACKEND, T, N)
        gpu(GPU_BACKEND, GPU_BACKEND_WORKGROUPSIZE)(Φ, x, ndrange=1)
        true
    catch
        false
    end
end

# import Pkg
# Pkg.activate("tog")
# using KernelAbstractions
# T=Float32
# N=3
# using Metal
# const GPU_BACKEND = MetalBackend()
# const GPU_BACKEND_WORKGROUPSIZE = 2^2^3
# const ○ = one(T) / (one(T) + one(T))
# const ○̂ = _ -> ○
# const Φ=(x)->○
# const Φ=_->○
# const Φ=○̂
# gpu_safe(Φ, N)
# gpu_safe(○̂,1)

# function gpu_safe(Φ, N, backend=GPU_BACKEND)
#     try
#         @kernel function _gpu_safe_test(out, @Const(f))
#             I = @index(Global)
#             x = ntuple(_ -> ○, N)
#             out[I] = f(x...)
#         end
#         out = KernelAbstractions.zeros(backend, T, 1)
#         Base.invokelatest() do
#             _gpu_safe_test(backend, 64)(out, Φ, ndrange=1)
#         end
#         KernelAbstractions.synchronize(backend)
#         true
#     catch
#         false
#     end
# end
# gpu_safe(○̂, 2)

# @kernel function _test_hardcoded(out, @Const(f))
#     I = @index(Global)
#     out[I] = f(○)
# end
# out = KernelAbstractions.zeros(GPU_BACKEND, T, 1)
# _test_hardcoded(GPU_BACKEND, 64)(out, ○̂, ndrange=1)
# KernelAbstractions.synchronize(GPU_BACKEND)
# println("OK: ", Array(out))

# @generated function gpu_safe(Φ, ::Val{N}, backend=GPU_BACKEND) where N
#     quote
#         try
#             @kernel function _gpu_safe_test(out, @Const(f))
#                 I = @index(Global)
#                 x = ntuple(_ -> ○, $N)
#                 out[I] = f(x...)
#             end
#             out = KernelAbstractions.zeros(backend, T, 1)
#             Base.invokelatest() do
#                 _gpu_safe_test(backend, 64)(out, Φ, ndrange=1)
#             end
#             KernelAbstractions.synchronize(backend)
#             true
#         catch
#             false
#         end
#     end
# end
# gpu_safe(○̂, Val(2))

# function gpu_safe(Φ, N)
#     try
#         @kernel gpu(Φ, x) = Φ(x)
#         x = KernelAbstractions.zeros(GPU_BACKEND, T, N)
#         gpu(GPU_BACKEND, GPU_BACKEND_WORKGROUPSIZE)(Φ, x, ndrange=1)
#         true
#     catch
#         false
#     end
# end
# struct ΦFunc{F}
#     f::F
# end
# @inline (Φ::ΦFunc)(x, atlas) = Φ.f(x)
# gpu_safe(○̂, 2)


# struct ΦFunc{F}
#     f::F
# end
# @inline (Φ::ΦFunc)(x, atlas) = Φ.f(x...)

# struct ΦTex3
#     H::Int32
#     W::Int32
#     offset::Int32
# end
# @inline function (Φ::ΦTex3)(x, atlas)
#     px = clamp(unsafe_trunc(Int32, x[1] * Φ.W) + Int32(1), Int32(1), Φ.W)
#     py = clamp(unsafe_trunc(Int32, x[2] * Φ.H) + Int32(1), Int32(1), Φ.H)
#     atlas[py, Φ.offset + px]
# end

# φ_const = ΦFunc(○̂)
# φs = ΦSet((φ_const, φ_hi))

# @kernel function _test_both(out, φ::ΦSet, @Const(atlas))
#     I = @index(Global)
#     x = (T(0.5), T(0.5))
#     out[1] = eval_Φ(φ, UInt32(1), x, atlas)
#     out[2] = eval_Φ(φ, UInt32(2), x, atlas)
# end

# out = KernelAbstractions.zeros(GPU_BACKEND, T, 2)
# atlas_gpu = adapt(GPU_BACKEND, mat_hi)
# Base.invokelatest() do
#     _test_both(GPU_BACKEND, 64)(out, φs, atlas_gpu, ndrange=1)
# end
# KernelAbstractions.synchronize(GPU_BACKEND)
# println("const: ", Array(out)[1], "  tex: ", Array(out)[2])

# @generated function eval_Φ(φ::ΦSet{Fs}, idx, x, atlas) where Fs
#     N = length(Fs.parameters)
#     branches = []
#     for i in 1:N
#         push!(branches, quote
#             if idx == $i
#                 return φ.fs[$i](x, atlas)
#             end
#         end)
#     end
#     quote
#         $(branches...)
#         return ○
#     end
# end
# Base.invokelatest() do
#     _test_both(GPU_BACKEND, 64)(out, φs, atlas_gpu, ndrange=1)
# end
# KernelAbstractions.synchronize(GPU_BACKEND)
# println("const: ", Array(out)[1], "  tex: ", Array(out)[2])