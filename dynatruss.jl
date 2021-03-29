using LinearAlgebra
using Plots
pyplot()

function dofTruss(coord::AbstractArray, support::AbstractArray)::Tuple
    ndof = 2size(coord)[1]
    kdof = []
    for i = 1:size(support)[1]
        node = support[i, 1]
        u, v = support[i, 2:end]
        if u == 1 && v == 1 # x-fixed & y-fixed
            push!(kdof, 2node - 1, 2node)
        else
            if u == 1 # x-fixed & y-free
                push!(kdof, 2node - 1)
            elseif v == 1 # x-free & y-fixed
                push!(kdof, 2node)
            end
        end
    end
    udof = setdiff((1:ndof)', kdof)
    return ndof, kdof, udof
end

function stiffnessTrussBar(EA::Float64, x₁, y₁, x₂, y₂)::Matrix{Float64}
    # Stiffnes matrix
    L = sqrt((x₁ - x₂)^2 + (y₁ - y₂)^2)
    c = (x₂ - x₁) / L
    s = (y₂ - y₁) / L
    CS = [
        +c +s
        -s +c
    ]
    # Transformation matrix
    T = [
        CS zeros(2, 2)
        zeros(2, 2) CS
    ]

    Kl = [
        +EA/L 0.0 -EA/L 0.0
        0.0 0.0 0.0 0.0
        -EA/L 0.0 +EA/L 0.0
        0.0 0.0 0.0 0.0
    ]

    return T' * Kl * T
end

function stiffnessTruss(
    prop::AbstractArray,
    coord::AbstractArray,
    cnt::AbstractArray,
    ndof::Int,
)::Matrix{Float64}
    Kg = zeros(ndof, ndof)
    nbar = size(cnt)[1]
    for i = 1:nbar
        n₁ = cnt[i, 1]
        n₂ = cnt[i, 2]
        dof = [2n₁ - 1; 2n₁; 2n₂ - 1; 2n₂]
        x₁ = coord[n₁, 1]
        y₁ = coord[n₁, 2]
        x₂ = coord[n₂, 1]
        y₂ = coord[n₂, 2]
        E = prop[i, 2]
        A = prop[i, 3]
        Kg[dof, dof] += stiffnessTrussBar(E * A, x₁, y₁, x₂, y₂)
    end
    return Kg
end

function forceTruss(ndof::Number, f::AbstractArray)::Vector{Float64}
    fg = zeros(ndof)
    for i = 1:size(f)[1]
        node = Int(f[i, 1])
        fg[2node-1] += f[i, 2]
        fg[2node-0] += f[i, 3]
    end
    return fg
end

"""
Solve truss for Static analysis
"""
function solveTruss(
    ndof::Number,
    udof::AbstractArray,
    K::Matrix{Float64},
    F::Vector{Float64},
)::Tuple
    d = Vector{Float64}(undef, ndof)
    u = K[udof, udof] \ f[udof]
    d[udof] = u
    return u, d
end

function stressTruss(
    prop::AbstractArray,
    coord::AbstractArray,
    cnt::AbstractArray,
    d::Vector,
)::Vector{Float64}
    nbar = size(cnt)[1]
    σₓₓ = zeros(nbar)
    for i = 1:nbar
        n₁ = cnt[i, 1]
        n₂ = cnt[i, 2]
        dof = [2n₁ - 1; 2n₁; 2n₂ - 1; 2n₂]
        x₁ = coord[n₁, 1]
        y₁ = coord[n₁, 2]
        x₂ = coord[n₂, 1]
        y₂ = coord[n₂, 2]
        L = sqrt((x₁ - x₂)^2 + (y₁ - y₂)^2)
        c = (x₂ - x₁) / L
        s = (y₂ - y₁) / L
        cs = [-c -s c s]
        E = prop[i, 2]
        σₓₓ[i] = (E / L) * dot(cs, d[dof])
    end
    return σₓₓ
end

function massTrussBar(
    masstype::String,
    ρ::Float64,
    A::Float64,
    Iz::Float64,
    x₁,
    y₁,
    x₂,
    y₂,
)::Matrix{Float64}
    m = []
    L = sqrt((x₁ - x₂)^2 + (y₁ - y₂)^2)
    c = (x₂ - x₁) / L
    s = (y₂ - y₁) / L
    CS = [
        +c +s
        -s +c
    ]
    # Transformation matrix
    T = [
        CS zeros(2, 2)
        zeros(2, 2) CS
    ]
    if masstype == "lumped"
        m = (A * ρ * L / 2) * Matrix(I, 4, 4)
    elseif masstype == "consistent"
        m = [
            ρ*A*L/3 0.0 ρ*A*L/6 0.0
            0.0 ρ*(A*L^2+3Iz)/3L 0.0 ρ*(A*L^2-6Iz)/6L
            ρ*A*L/6 0.0 ρ*A*L/6 0.0
            0.0 ρ*(A*L^2-6Iz)/6L 0.0 ρ*(A*L^2+3Iz)/3L
        ]
    end
    return m
end

function massTruss(
    masstype::String,
    prop::AbstractArray,
    coord::AbstractArray,
    cnt::AbstractArray,
    ndof::Int,
)::Matrix{Float64}
    Mg = zeros(ndof, ndof)
    nbar = size(cnt)[1]
    for i = 1:nbar
        n₁ = cnt[i, 1]
        n₂ = cnt[i, 2]
        dof = [2n₁ - 1; 2n₁; 2n₂ - 1; 2n₂]
        x₁ = coord[n₁, 1]
        y₁ = coord[n₁, 2]
        x₂ = coord[n₂, 1]
        y₂ = coord[n₂, 2]
        ρ = prop[i, 1]
        A = prop[i, 3]
        Iz = prop[i, 4]
        Mg[dof, dof] += massTrussBar(masstype, ρ, A, Iz, x₁, y₁, x₂, y₂)
    end
    return Mg
end

function plotTruss(
    plotlabels::AbstractArray,
    coord::AbstractArray,
    cnt::AbstractArray,
    f::AbstractArray,
    δ::Vector,
    σ::Vector,
)
    plot(fontfamily = "serif", framestyle = :box, legend=false)
    nbar = size(cnt)[1]
    # undeformed
    for i = 1:size(cnt)[1]
        barx = [coord[cnt[i, 1], 1], coord[cnt[i, 2], 1]]
        bary = [coord[cnt[i, 1], 2], coord[cnt[i, 2], 2]]
        plot!(barx, bary, line = (2), c = :black)
    end

    # deformed
    scale = 100
    nnodes = size(coord)[1]
    δᵤ = δ[1:2:2nnodes-1]
    δᵥ = δ[2:2:2nnodes]
    newcoordx = coord[:, 1] + δᵤ * scale
    newcoordy = coord[:, 2] + δᵥ * scale
    stresscolor = [:red, :blue, :gray]
    for i = 1:nbar
        barx = [newcoordx[cnt[i, 1]], newcoordx[cnt[i, 2]]]
        bary = [newcoordy[cnt[i, 1]], newcoordy[cnt[i, 2]]]
        if σ[i] > 1e-6
            plot!(
                barx,
                bary,
                line = (:dash, 2),
                c = stresscolor[1],
            )
        elseif σ[i] < -1e-6
            plot!(
                barx,
                bary,
                line = (:dash, 2),
                c = stresscolor[2],
            )
        else
            plot!(
                barx,
                bary,
                line = (:dash, 2),
                c = stresscolor[3],
            )
        end
    end

    for i = 1:length(newcoordx)
        annotate!([(
            newcoordx[i],
            newcoordy[i],
            Plots.text("$i", 10, :magenta, :bottom),
        )])
    end
    for i = 1:size(f)[1]
        node = Int(f[i, 1])
        if f[i, 2] != 0.0
            ax = [
                newcoordx[node]
                newcoordx[node] + 0.1 * maximum(coord[:, 1]) * sign(f[i, 2])
            ]
            ay = [
                newcoordy[node]
                newcoordy[node]
            ]
            plot!(ax, ay, arrow = :closed, c = :magenta)
        elseif f[i, 3] != 0.0
            ax = [
                newcoordx[node]
                newcoordx[node]
            ]
            ay = [
                newcoordy[node]
                newcoordy[node] + 0.1 * maximum(coord[:, 2]) * sign(f[i, 3])
            ]
            plot!(ax, ay, arrow = :closed, c = :magenta)
        end
    end
    scatter!(
        coord[:, 1],
        coord[:, 2],
        marker = (:black, :circle),
        label = "Undeformed",
    )
    scatter!(
        newcoordx,
        newcoordy,
        marker = (:magenta, :circle),
        label = "Deformed",
    )
    title!(plotlabels[1])
    xaxis!(plotlabels[2])
    yaxis!(plotlabels[3])

end
