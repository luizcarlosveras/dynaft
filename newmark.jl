using LinearAlgebra

function newmarkβ(
    u₀::Vector,
    v₀::Vector,
    Δt::Float64,
    t⁰::Float64,
    tᶠ::Float64,
    γ::Float64,
    β::Float64,
    𝛼,ϐ,
    m::Matrix,
    k::Matrix,
    f::Function,
)::AbstractArray
    """
    (1) Constant average acceleration method (γ = 1/2 , β = 1/4)
    (2) Linear acceleration method (γ = 1/2 , β = 1/6)
    k = k[udof,udof]
    m = m[udof,udof]
    f :: Function
    """
    # Eigenvectors and eigenvalues
    ω², Φ = eigen(k, m)
    # Natural frequency 
    ωₙ = real.(sqrt.(Complex.(ω²)))
    # Rayleigh damping
    # 𝛼 = 2ξ * ωₙ(1) * ωₙ(2) / (ωₙ(1) + ωₙ(2))
    # ϐ = 2ξ / (ωₙ(1) + ωₙ(2))
    c = 𝛼 * m + ϐ * k
    # Modal matrix
    md = Φ' * m * Φ
    cd = Φ' * c * Φ
    kd = Φ' * k * Φ
    # F = Φ' * f
    # Initial calculations
    t = t⁰:Δt:tᶠ
    q₀ = zeros(size(Φ)[2])
    ∂q₀ = zeros(size(Φ)[2])
    for i = 1:size(Φ)[2]
        q₀[i] = Φ[:,i]' * m * u₀ / (Φ[:,i]' * m * Φ[:,i])
        ∂q₀[i] = Φ[:,i]' * m * v₀ / (Φ[:,i]' * m * Φ[:,i])
    end
    ∂²q₀ = md \ (Φ' * f(t⁰) - cd * v₀ - kd * u₀)
    a₁ = (1 / (β * Δt^2)) * md + (γ / (β * Δt)) * cd
    a₂ = (1 / (β * Δt)) * md + (γ / β - 1) * cd
    a₃ = (1 / (2β) - 1) * md + Δt * (γ / (2β) - 1) * cd
    𝓚 = kd + a₁

    u = zeros(length(udof), length(t))
    q = zeros(length(udof), length(t))
    ∂q = zeros(length(udof), length(t))
    ∂²q = zeros(length(udof), length(t))

    u[:, 1] = u₀
    q[:, 1] = q₀
    ∂q[:, 1] = ∂q₀
    ∂²q[:, 1] = ∂²q₀

    for i = 1:length(t)-1
        Fᵢ₊₁ = Φ' * f(t[i]) + a₁ * q[:,i] + a₂ * ∂q[:,i] + a₃ * ∂²q[:,i]
        q[:, i+1] = 𝓚 \ Fᵢ₊₁
        ∂q[:, i+1] =
            γ / (β * Δt) * (q[:, i+1] - q[:, i]) +
            (1 - γ / β) * ∂q[:,i] +
            Δt * (1 - γ / (2β)) * ∂²q[:,i]
        ∂²q[:, i+1] =
            1 / (β * Δt^2) * (q[:, i+1] - q[:, i]) - 1 / (β * Δt) * ∂q[:,i] -
            (1 / (2β) - 1) * ∂²q[:,i]
        u[:,i+1] = Φ * q[:,i]
    end

    return u
end #function
