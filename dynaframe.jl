using LinearAlgebra

function stiffnessFrameBar(
    EA::Float64,
    EI::Float64,
    bartype::String,
    x₁, y₁, x₂, y₂
)::Matrix{Float64}
    # Stiffnes matrix
    Kl = zeros(6, 6)
    L = sqrt((x₁ - x₂)^2 + (y₁ - y₂)^2)
    c = (x₂ - x₁) / L
    s = (y₂ - y₁) / L
    CS = [
        +c +s 0.0
        -s +c 0.0
        0.0 0.0 1.0
    ]
    # Transformation matrix
    T = [
        CS zeros(3, 3)
        zeros(3, 3) CS
    ]
    L² = L^2
    L³ = L^3
    # Define the stiffness matrix by bar type
    if bartype == "ff" # fixed-fixed
        Kl = [
            +EA/L 0.0 0.0 -EA/L 0.0 0.0
            0.0 +12EI/L³ +6EI/L² 0.0 -12EI/L³ +6EI/L²
            0.0 +6EI/L² +4EI/L 0.0 -6EI/L² +2EI/L
            -EA/L 0.0 0.0 +EA/L 0.0 0.0
            0.0 -12EI/L³ -6EI/L² 0.0 +12EI/L³ -6EI/L²
            0.0 +6EI/L² +2EI/L 0.0 -6EI/L² +4EI/L
        ]
    elseif bartype == "pf" # pinned-fixed
        Kl = [
            +EA/L 0.0 0.0 -EA/L 0.0 0.0
            0.0 +3EI/L³ 0.0 0.0 -3EI/L³ +3EI/L²
            0.0 0.0 0.0 0.0 0.0 0.0
            -EA/L 0.0 0.0 +EA/L 0.0 0.0
            0.0 -3EI/L³ 0.0 0.0 +3EI/L³ -3EI/L²
            0.0 +3EI/L² 0.0 0.0 -3EI/L² +3EI/L
        ]
    elseif bartype == "fp" # fixed-pinned
        Kl = [
            +EA/L 0.0 0.0 -EA/L 0.0 0.0
            0.0 +3EI/L³ +3EI/L² 0.0 -3EI/L³ 0.0
            0.0 +3EI/L² +3EI/L 0.0 -3EI/L² 0.0
            -EA/L 0.0 0.0 +EA/L 0.0 0.0
            0.0 -3EI/L³ -3EI/L² 0.0 +3EI/L³ 0.0
            0.0 0.0 0.0 0.0 0.0 0.0
        ]
    else
        println("Type not supported!")
    end

    return T' * Kl * T

end

function stiffnessFrame(
    E::Float64,
    A::Float64,
    Iz::Float64,
    bartype::AbstractArray,
    coord::AbstractArray,
    cnt::AbstractArray,
    ndof::Int
)::Matrix{Float64}
    Kg = zeros(ndof,ndof)
    nbar = size(cnt)[1]
    for i = 1:nbar
        n₁ = cnt[i,1]
        n₂ = cnt[i,2]
        dof = [3n₁-2; 3n₁-1; 3n₁; 3n₂-2; 3n₂-1; 3n₂]
        x₁ = coord[n₁,1]; y₁ = coord[n₁,2]
        x₂ = coord[n₂,1]; y₂ = coord[n₂,2]
        Kg[dof,dof] += stiffnessFrameBar(E*A, E*Iz, bartype[i], x₁, y₁, x₂, y₂)
    end
    return Kg
end

function dofFrame(coord::AbstractArray, support::AbstractArray)::Tuple
    ndof = 3size(coord)[1]
    kdof = []
    for i = 1:size(support)[1]
        u, v, ϕ = support[i,2:end]
        if u == 1 && v == 1 && ϕ == 1 # fixed
            node = support[i,1]
            push!(kdof, 3node-2, 3node-1, 3node)
        else
            if u == 1 # x-fixed & y-free
                push!(kdof, 3node-2)
            elseif v == 1 # x-free & y-fixed
                push!(kdof, 3node-1)
            else
                println("Boundary conditions are not ok!")
            end
        end
    end
    udof = setdiff((1:ndof)', kdof)
    return ndof, kdof, udof
end

function forceFrame(ndof::Number, f::AbstractArray)::Vector{Float64}
    Fg = zeros(ndof)
    for i = 1:length(f)
        node = f[i,1]
        Fg[3node-2] = f[i,2]
        Fg[3node-1] = f[i,3]
        Fg[3node-0] = f[i,4]
    end
    return Fg
end

function solveFrame(
    ndof::Number, 
    udof::AbstractArray, 
    f::Vector{Float64}, 
    K::Matrix{Float64}
)::Vector{Float64}
    """
    Static analysis
    """
    d = Vector{Float64}(undef, ndof)
    u = K[udof,udof]\f[udof]
    d[udof] = u
    return d
end

function massFrameBar(
    masstype::String,
    ρ::Float64,
    A::Float64,
    Iz::Float64,
    x₁, y₁, x₂, y₂
)::Matrix{Float64}
    Ml = []
    L = sqrt((x₁ - x₂)^2 + (y₁ - y₂)^2)
    c = (x₂ - x₁) / L
    s = (y₂ - y₁) / L
    CS = [
        +c +s 0.0
        -s +c 0.0
        0.0 0.0 1.0
    ]
    # Transformation matrix
    T = [
        CS zeros(3, 3)
        zeros(3, 3) CS
    ]
    L²=L^2
    if masstype == "lumped"
        Ml = (ρ * A * L / 2) * Diagonal([1,1,L²/16,1,1,L²/16])
    elseif masstype == "consistent"
        Ml = [
            ρ*A*L/3 0.0 0.0 ρ*A*L/6 0.0 0.0;
            0.0 ρ*(13A*L²+42Iz)/35L ρ*(11A*L²+21Iz)/210L 0.0 3ρ*(3A*L²-28Iz)/70L -ρ*(13A*L²-42Iz)/420L;
            0.0 ρ*(11A*L²+21Iz)/210L ρ*(A*L²+14Iz)/105L 0.0 ρ*(13A*L²-42Iz)/420L -ρ*(3A*L²+14Iz)/420L;
            ρ*A*L/6 0.0 0.0 ρ*A*L/3 0.0 0.0;
            0.0 3ρ*(3A*L²-28Iz)/70L ρ*(13A*L²-42Iz)/420L 0.0 ρ*(13A*L²+42Iz)/35L -ρ*(11A*L²+21Iz)/210L;
            0.0 -ρ*(13A*L²-42Iz)/420L -ρ*(3A*L²+14Iz)/420L 0.0 -ρ*(11A*L²+21Iz)/210L ρ*(A*L²+14Iz)/105L;
        ]
    end

    return T' * Ml * T

end

function massFrame(
    masstype::String,
    ρ::Float64,
    A::Float64,
    Iz::Float64,
    coord::AbstractArray,
    cnt::AbstractArray,
    ndof::Int
)::Matrix{Float64}
    Mg = zeros(ndof,ndof)
    nbar = size(bar)[1]
    for i = 1:nbar
        n₁ = cnt[i,1]
        n₂ = cnt[i,2]
        dof = [2n₁ 2n₁+1 2n₂ 2n₂+1]
        x₁ = coord[n₁,1]; y₁ = coord[n₁,2]
        x₂ = coord[n₂,1]; y₂ = coord[n₂,2]
        Mg[dof,dof] += massFrameBar(masstype, ρ, A, Iz, x₁, y₁, x₂, y₂)
    end
    return Mg
end
