using Plots
pyplot()

include("dynatruss.jl")
include("dynaframe.jl")
include("newmark.jl")

# ------------- #
# Truss example #
# ------------- #

ρ = 78 # kN/m³
E = 1.2e7 # kN/m²
A = 3e-2 # m²
Iz = 4e-3 # m⁴

coord =[
    0.0   0.0;
    0.0   4.0;
    4.0   0.0;
] # X Y -> in
cnt = [
    1 2;
    2 3;
    3 1;
] # bar: n1 n2
support = [
    1 1 1;
    3 0 1
] # node x_direction y_direction -> 0 is free 1 is fixed
load = [
    4 0.0 -1e3; # lb
]# node x_value y_value
prop = [
    ρ E A Iz;
    ρ E A Iz;
    ρ E A Iz;
]

# Matriz de rigidez
ndof, kdof, udof = dofTruss(coord, support)
k = stiffnessTruss(prop,coord,cnt,ndof)
# Parâmetros para problema dinâmico
u₀ = [0.0; 0.0; 0.0]
v₀ = [0.0; 0.0; 0.0]
Δt,t⁰,tᶠ = 0.01,0.0,5.0
γ, β = 1/2, 1/4
𝛼, ϐ = 0.4, 0.7
ωᶠ = 5.0
f(tᵢ) = [
    2  + cos(tᵢ)*sin(ωᶠ*tᵢ);
    -5 + cos(tᵢ)*sin(ωᶠ*tᵢ);
    7  + sin(tᵢ)*sin(ωᶠ*tᵢ)
]
# Massa concentrada
masstype = "lumped"
m1 = massTruss(masstype,prop,coord,cnt,ndof)
u1 = newmarkβ(u₀,v₀,Δt,t⁰,tᶠ,γ,β,𝛼,ϐ,m1[udof,udof],k[udof,udof],f)
# Matriz de Massa Com Inércia Rotária
masstype = "consistent"
m2 = massTruss(masstype,prop,coord,cnt,ndof)
u2 = newmarkβ(u₀,v₀,Δt,t⁰,tᶠ,γ,β,𝛼,ϐ,m2[udof,udof],k[udof,udof],f)
# Matriz de Massa Sem Inércia Rotária
masstype = "consistent"
Iz = 0.0
prop = [
    ρ E A Iz;
    ρ E A Iz;
    ρ E A Iz;
]
m3 = massTruss(masstype,prop,coord,cnt,ndof)
u3 = newmarkβ(u₀,v₀,Δt,t⁰,tᶠ,γ,β,𝛼,ϐ,m3[udof,udof],k[udof,udof],f)

t = t⁰:Δt:tᶠ

plot(t,u1[1,:],line=(3,:solid,:red))
plot!(t,u2[1,:],line=(2,:dash,:green))
plot!(t,u3[1,:],line=(:dash,:blue))