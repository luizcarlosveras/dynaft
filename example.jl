using DelimitedFiles

include("dynatruss.jl")
include("dynaframe.jl")
include("newmark.jl")

# ------------- #
# Truss example #
# ------------- #

# Example 12.15
# Adapted from Mechanical Vibrations by Singiresu S. Rao - 5ed

# Input data
ρ = 0.3 # lb/in³
E = 30e6 # psi or lb/in²
A₁ = 2.0 # in²
A₂ = 1.0 # in²
Iz₁ = 4 * A₁^4 / π^3 # in⁴
Iz₂ = 4 * A₂^4 / π^3 # in⁴
I₀ = 0.0 # in⁴

# Inertia neglected
propⁿ = [
    ρ E A₁ I₀
    ρ E A₁ I₀
    ρ E A₂ I₀
    ρ E A₂ I₀
]
# consistent mass
propᶜ = [
    ρ E A₁ Iz₁
    ρ E A₁ Iz₁
    ρ E A₂ Iz₂
    ρ E A₂ Iz₂
]
coord = [
    0.0 0.0
    100.0 0.0
    50.0 25.0
    200.0 100.0
] # X Y -> in
cnt = [
    1 3
    2 3
    3 4
    2 4
] # bar: n1 n2
support = [
    1 1 1
    2 1 1
] # node x_direction y_direction -> 0 is free 1 is fixed
load = [
    4 0.0 -1e3; # lb
]# node x_value y_value

# --------------Static analysis-------------- #
ndof, kdof, udof = dofTruss(coord, support)
K = stiffnessTruss(propⁿ, coord, cnt, ndof)
f = forceTruss(ndof, load)
u, d = solveTruss(ndof, udof, K, f)
σ = stressTruss(propⁿ, coord, cnt, d)
plotlabels = ["Truss - Example 12.15 from Rao", "Length (in)", "Height (in)"]
# plotTruss(plotlabels,coord,cnt,load,d,σ)
# savefig("truss_static.svg")
# ------------------------------------------- #



# --------------Dynamic analysis-------------- #

# Parameters for dynamic problem
u₀ = [0.0; 0.0; 0.0; 0.0]
v₀ = [0.0; 0.0; 0.0; 0.0]

Δt, t⁰, tᶠ = 0.01, 0.0, 10.0
# Constant average acceleration method
γ, β = 1 / 2, 1 / 4
# Rayleigh's coeffcients
𝛼, ϐ = 1.5, 2.5
# Force functions for each d.o.f. of the free nodes
ωᶠ = 2.0
pₘ = [1; -10; 1; -10]
p(tᵢ) =
    pₘ + [
        cos(tᵢ) * sin(ωᶠ * tᵢ)
        -sin(tᵢ) * sin(ωᶠ * tᵢ)
        10cos(tᵢ) * sin(ωᶠ * tᵢ)
        -10sin(tᵢ) * sin(ωᶠ * tᵢ)
    ]
# Stiffnes matrix
ndof, kdof, udof = dofTruss(coord, support)
k = stiffnessTruss(propⁿ, coord, cnt, ndof)
# Lumped Mass Matrix
begin
    masstype = "lumped"
    m1 = massTruss(masstype, propⁿ, coord, cnt, ndof)
    # Eigenvectors and eigenvalues
    ω², Φ = eigen(k[udof, udof], m1[udof, udof])
    # Natural frequency 
    ωₙ = real.(sqrt.(Complex.(ω²)))
    println("Natural frequency and Eigenvectors\nfor Lumped Mass Matrix")
    writedlm(stdout, ωₙ)
    writedlm(stdout, Φ)
    # Displacements for Lumped Mass Matrix
    u1 = newmarkβ(u₀, v₀, Δt, t⁰, tᶠ, γ, β, 𝛼, ϐ, m1[udof, udof], k[udof, udof], p)
end

# Mass Matrix With Rotational Inertia 
begin
    masstype = "consistent"
    m2 = massTruss(masstype, propᶜ, coord, cnt, ndof)
    # Eigenvectors and eigenvalues
    ω², Φ = eigen(k[udof, udof], m2[udof, udof])
    # Natural frequency 
    ωₙ = real.(sqrt.(Complex.(ω²)))
    println(
        "Natural frequency and Eigenvectors\nfor Consistent Mass Matrix\nWith Rotational Inertia",
    )
    writedlm(stdout, ωₙ)
    writedlm(stdout, Φ)
    # Displacements for Mass Matrix With Rotational Inertia 
    u2 = newmarkβ(u₀, v₀, Δt, t⁰, tᶠ, γ, β, 𝛼, ϐ, m2[udof, udof], k[udof, udof], p)
end

# Mass Matrix Without Rotational Inertia
begin
    masstype = "consistent"
    m3 = massTruss(masstype, propⁿ, coord, cnt, ndof)
    # Eigenvectors and eigenvalues
    ω², Φ = eigen(k[udof, udof], m3[udof, udof])
    # Natural frequency 
    ωₙ = real.(sqrt.(Complex.(ω²)))
    println(
        "Natural frequency and Eigenvectors\nfor consistent Mass Matrix\nWithout Rotational Inertia",
    )
    writedlm(stdout, ωₙ)
    writedlm(stdout, Φ)
    # Displacements for Mass Matrix Without Rotational Inertia 
    u3 = newmarkβ(u₀, v₀, Δt, t⁰, tᶠ, γ, β, 𝛼, ϐ, m3[udof, udof], k[udof, udof], p)
end

t = t⁰:Δt:tᶠ
plt1 = plot(
    t,
    [u1[1, :] u2[1, :] u3[1, :]],
    label = ["Lumped" "Consistent with\nrotational inertia" "Consistent without\nrotational inertia"],
    legend = :outertop,
    lw = [3 2 1],
    line = [:solid :dash :dash],
    linecolor = [:red :green :blue],
)
plt2 = plot(
    t,
    [u1[2, :] u2[2, :] u3[2, :]],
    legend = false,
    lw = [3 2 1],
    line = [:solid :dash :dash],
    linecolor = [:red :green :blue],
)
plt3 = plot(
    t,
    [u1[3, :] u2[3, :] u3[3, :]],
    legend = false,
    lw = [3 2 1],
    line = [:solid :dash :dash],
    linecolor = [:red :green :blue],
)
plt4 = plot(
    t,
    [u1[4, :] u2[4, :] u3[4, :]],
    legend = false,
    lw = [3 2 1],
    line = [:solid :dash :dash],
    linecolor = [:red :green :blue],
)
plot(plt1, plt2, plt3, plt4, layout = (2, 2))
savefig("displacements_dynamics.svg")


# ------------- #
# Frame example #
# ------------- #
