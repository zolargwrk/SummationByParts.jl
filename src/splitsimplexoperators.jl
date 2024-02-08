using SummationByParts
using SummationByParts.Cubature
using SummationByParts.OrthoPoly

# This file gathers functions that are used to construct multidimensional 
# SBP operators by splitting simplices into tensor-product quad and hexs

"""
### SummationByParts.tensor_lgl_quad_nodes
"""
function tensor_lgl_quad_nodes(p::Int)
    qf = 2*(p+1)-3 
    cub, vtx = SummationByParts.Cubature.quadrature(qf, internal=false)
    x = sort(vec(SymCubatures.calcnodes(cub, vtx)))'
    y = x[1] .*ones(1,cub.numnodes)
    xy = vcat(x,y)
    for i = 1:cub.numnodes
        xy = hcat(xy, vcat(x,x[i].*ones(1,cub.numnodes)))
    end
    return xy[:,cub.numnodes+1:end] 
end

"""
### SummationByParts.tensor_lgl_hex_nodes
"""
function tensor_lgl_hex_nodes(p::Int)
    qf = 2*(p+1)-3 
    cub, vtx = SummationByParts.Cubature.quadrature(qf, internal=false)
    x = sort(vec(SymCubatures.calcnodes(cub, vtx)))'
    y = x[1] .*ones(1,cub.numnodes)
    z = x[1] .*ones(1,cub.numnodes)
    xyz = vcat(x,y,z)
    for i = 1:cub.numnodes
        xy = vcat(x, x[i].*ones(1,cub.numnodes))
        for j = 1:cub.numnodes
            xyz = hcat(xyz, vcat(xy,x[j].*ones(1,cub.numnodes)))
        end
    end
    xyz[[2,3],:] = xyz[[3,2],:]
    return xyz[:,cub.numnodes+1:end]
end

"""
### SummationByParts.square_quad_map
"""
function square_quad_map(xp::Array{T},quad_vert::Array{T}) where T
    xi = xp[1]
    eta = xp[2]
    psi = []
    push!(psi, 1/4*(1-xi)*(1-eta))
    push!(psi, 1/4*(1+xi)*(1-eta))
    push!(psi, 1/4*(1-xi)*(1+eta))
    push!(psi, 1/4*(1+xi)*(1+eta))
    
    x = zeros(2,1)
    for i=1:2
        for j=1:4
            x[i] += quad_vert[j,i]*psi[j]
        end
    end
    return x 
end

"""
### SummationByParts.cube_hex_map.jl 
"""
function cube_hex_map(xp::Array{T},hex_vert::Array{T}) where T
    xi = xp[1]
    eta = xp[2]
    zeta = xp[3]
    psi = []

    push!(psi, 1/8*(1-xi)*(1-eta)*(1-zeta))
    push!(psi, 1/8*(1+xi)*(1-eta)*(1-zeta))
    push!(psi, 1/8*(1-xi)*(1+eta)*(1-zeta))
    push!(psi, 1/8*(1+xi)*(1+eta)*(1-zeta))
    push!(psi, 1/8*(1-xi)*(1-eta)*(1+zeta))
    push!(psi, 1/8*(1+xi)*(1-eta)*(1+zeta))
    push!(psi, 1/8*(1-xi)*(1+eta)*(1+zeta))
    push!(psi, 1/8*(1+xi)*(1+eta)*(1+zeta))

    x = zeros(3,1)
    for i=1:3
        for j=1:8
            x[i] += hex_vert[j,i]*psi[j]
        end
    end

    return x 
end

"""
### SummationByParts.square_to_tri_map
"""
function square_to_tri_map(xi::Array{T}) where T
    quad_vert = get_quad_vert()
    n = size(xi,2)
    x = zeros(2,3*n)
    for i=0:2 
        for j=1:n
            xp = square_quad_map(xi[:,j], quad_vert[i+1])
            x[:,i*n+j] = xp 
        end
    end
    return x
end

function get_quad_vert()
    quad_vert = [[-1 -1; 0 -1; -1 0; -1/3 -1/3], 
                 [0 -1; 1 -1; -1/3 -1/3; 0 0], 
                 [-1 0; -1/3 -1/3; -1 1; 0 0]]
    return quad_vert
end

function get_hex_vert(;T=Float64)
    v1 = T[-1 -1 -1]
    v2 = T[1 -1 -1]
    v3 = T[-1 1 -1]
    v4 = T[-1 -1 1]
    v5 = T[0 0 -1]
    v6 = T[-1 0 -1]
    v7 = T[0 -1 -1]
    v8 = T[-1 0 0]
    v9 = T[-1 -1 0]
    v10 = T[0 -1 0]
    v11 = T[-1/3 -1/3 -1]
    v12 = T[-1/3 -1/3 -1/3]
    v13 = T[-1 -1/3 -1/3]
    v14 = T[-1/3 -1 -1/3]
    v15 = T[-1/2 -1/2 -1/2]
    hex_vert = [[v2; v5; v7; v11; v10; v12; v14; v15],
                     [v5; v3; v11; v6; v12; v8; v15; v13],
                     [v7; v11; v1; v6; v14; v15; v9; v13],
                     [v10; v4; v12; v8; v14; v9; v15; v13]]
    return hex_vert
end

# function get_hex_vert()
#     v1 = [-1 -1/sqrt(3) -sqrt(6)/6]
#     v2 = [1 -1/sqrt(3) -sqrt(6)/6]
#     v3 = [0 2/sqrt(3) -sqrt(6)/6]
#     v4 = [0 0 sqrt(6)/2]
#     v5 = [1/2 sqrt(3)/6 -sqrt(6)/6]
#     v6 = [-1/2 sqrt(3)/6 -sqrt(6)/6]
#     v7 = [0 -1/sqrt(3) -sqrt(6)/6]
#     v8 = [0 1/sqrt(3) sqrt(6)/6]
#     v9 = [-1/2 -1/(2*sqrt(3)) sqrt(6)/6]
#     v10 = [1/2 -1/(2*sqrt(3)) sqrt(6)/6]
#     v11 = [0 0 -sqrt(6)/6]
#     v12 = [1/3 sqrt(3)/9 sqrt(6)/18]
#     v13 = [-1/3 sqrt(3)/9 sqrt(6)/18]
#     v14 = [0 -2*sqrt(3)/9 sqrt(6)/18]
#     v15 = [0 0 0]
#     # hex_vert_equi = [[v2; v5; v11; v7; v14; v10; v12; v15],
#     #                  [v3; v6; v11; v5; v12; v8; v13; v15],
#     #                  [v1; v7; v11; v6; v13; v9; v14; v15],
#     #                  [v4; v8; v13; v9; v14; v10; v12; v15]]
#     # hex_vert_equi = [[v2; v5; v7; v11; v10; v12; v14; v15],
#     #                  [v5; v3; v11; v6; v12; v8; v15; v13],
#     #                  [v7; v11; v1; v6; v14; v15; v9; v13],
#     #                  [v10; v12; v14; v15; v9; v13; v4; v8]]
#     hex_vert_equi = [[v2; v5; v7; v11; v10; v12; v14; v15],
#                      [v5; v3; v11; v6; v12; v8; v15; v13],
#                      [v7; v11; v1; v6; v14; v15; v9; v13],
#                      [v10; v4; v12; v8; v14; v9; v15; v13]]
#     # hex_vert_equi = [[v7; v2; v11; v5; v14; v10; v15; v12],
#     #                  [v5; v3; v11; v6; v12; v8; v15; v13],
#     #                  [v7; v11; v1; v6; v14; v15; v9; v13],
#     #                  [v10; v4; v12; v8; v14; v9; v15; v13]]
#     hex_vert = []
#     for i = 1:4 
#         vert = Matrix((hex_vert_equi[i])')
#         push!(hex_vert, Matrix(equi_to_perp_tet_map(vert)'))
#     end
#     return hex_vert
# end

"""
### SummationByParts.cube_to_tet_map.jl 
"""
function cube_to_tet_map(xi::Array{T}) where T
    hex_vert = get_hex_vert()
    n = size(xi,2)
    x = zeros(3,4*n)
    for i=0:3 
        for j=1:n
            xp = cube_hex_map(xi[:,j], Matrix(hex_vert[i+1]))
            x[:,i*n+j] = xp 
        end
    end
    return x
end

"""
### SummationByParts.perp_to_equi_tri_map 
"""
function perp_to_equi_tri_map(x::Array{T}) where T
    vtx = T[-1 -1/sqrt(3); 1 -1/sqrt(3); 0 2/sqrt(3)]
    xequi = zeros(size(x))
    for i = 1:2
        xequi[i,:] .= -0.5 .* (x[1,:].+x[2,:])*vtx[1,i] .+
                    0.5 .*(x[1,:].+1.0)*vtx[2,i] .+
                    0.5 .*(x[2,:].+1.0)*vtx[3,i]
    end
    return xequi
end

"""
### SummationByParts.perp_to_equi_tri_map
"""
function perp_to_equi_tet_map(x::Array{T}) where T
    v1 = [-1,-1/sqrt(3), -1/sqrt(6)]; v2 = [ 1,-1/sqrt(3), -1/sqrt(6)];
    v3 = [ 0, 2/sqrt(3), -1/sqrt(6)]; v4 = [ 0, 0/sqrt(3),  3/sqrt(6)];
    xequi = zeros(size(x))
    for i = 1:3
        xequi[i,:] .= - 0.5 .* (x[1,:].+x[2,:].+x[3,:].+1.0).*v1[i] .+ 
                        0.5 .*(x[1,:].+1.0).*v2[i] .+ 
                        0.5 .*(x[2,:].+1.0).*v3[i] .+ 
                        0.5 .*(x[3,:].+1.0).*v4[i]
    end
    return xequi
end

"""
### SummationByParts.perp_to_equi_tri_map
"""
function equi_to_perp_tet_map(x::Array{T}) where T
    v1 = [-1,-1/sqrt(3), -1/sqrt(6)]; v2 = [ 1,-1/sqrt(3), -1/sqrt(6)];
    v3 = [ 0, 2/sqrt(3), -1/sqrt(6)]; v4 = [ 0, 0/sqrt(3),  3/sqrt(6)];

    rhs = x - 0.5 .* (v2'.+v3'.+v4' .-v1')' *ones(1,size(x,2))
    A = [0.5.*(v2 .- v1)'; 0.5 .*(v3 .- v1)'; 0.5 .* (v4 .- v1)']'
    xperp = A\rhs
    return xperp
end

"""
### SummationByParts.metric_tri
"""
function metric_tri!(xp::Array{T},quad_vert::Array{T},dxi::SubArray{T},dx::SubArray{T},Jac::SubArray{T}) where T
    ξ = xp[1]
    η = xp[2]

    ∂Ψ∂ξ = []
    ∂Ψ∂η = []
    push!(∂Ψ∂ξ, -1/4*(1-η))
    push!(∂Ψ∂ξ, 1/4*(1-η))
    push!(∂Ψ∂ξ, -1/4*(1+η))
    push!(∂Ψ∂ξ, 1/4*(1+η))

    push!(∂Ψ∂η, -1/4*(1-ξ))
    push!(∂Ψ∂η, -1/4*(1+ξ))
    push!(∂Ψ∂η, 1/4*(1-ξ))
    push!(∂Ψ∂η, 1/4*(1+ξ))

    ∂x∂ξ = 0.0
    ∂x∂η = 0.0
    ∂y∂ξ = 0.0
    ∂y∂η = 0.0

    for j=1:4
        ∂x∂ξ += quad_vert[j,1]*∂Ψ∂ξ[j]
        ∂x∂η += quad_vert[j,1]*∂Ψ∂η[j]
        ∂y∂ξ += quad_vert[j,2]*∂Ψ∂ξ[j]
        ∂y∂η += quad_vert[j,2]*∂Ψ∂η[j]
    end

    J = ∂x∂ξ*∂y∂η - ∂x∂η*∂y∂ξ
    ∂ξ∂x = ∂y∂η/J 
    ∂ξ∂y = -∂x∂η/J 
    ∂η∂x = -∂y∂ξ/J 
    ∂η∂y = ∂x∂ξ/J 

    Jac[1,1]=J
    dxi[1,1]=∂x∂ξ
    dxi[2,1]=∂x∂η
    dxi[3,1]=∂y∂ξ
    dxi[4,1]=∂y∂η

    dx[1,1]=∂ξ∂x
    dx[2,1]=∂ξ∂y
    dx[3,1]=∂η∂x
    dx[4,1]=∂η∂y
end

"""
### SummationByParts.metric_tet
"""
function metric_tet!(xp::Array{T},hex_vert::Array{T},dxi::SubArray{T},dx::SubArray{T},Jac::SubArray{T}) where T
    ξ = xp[1]
    η = xp[2]
    ζ = xp[3]
    ∂Ψ∂ξ = []
    ∂Ψ∂η = []
    ∂Ψ∂ζ = []

    push!(∂Ψ∂ξ, -1/8*(1-η)*(1-ζ))
    push!(∂Ψ∂ξ, 1/8*(1-η)*(1-ζ))
    push!(∂Ψ∂ξ, -1/8*(1+η)*(1-ζ))
    push!(∂Ψ∂ξ, 1/8*(1+η)*(1-ζ))
    push!(∂Ψ∂ξ, -1/8*(1-η)*(1+ζ))
    push!(∂Ψ∂ξ, 1/8*(1-η)*(1+ζ))
    push!(∂Ψ∂ξ, -1/8*(1+η)*(1+ζ))
    push!(∂Ψ∂ξ, 1/8*(1+η)*(1+ζ))

    push!(∂Ψ∂η, -1/8*(1-ξ)*(1-ζ))
    push!(∂Ψ∂η, -1/8*(1+ξ)*(1-ζ))
    push!(∂Ψ∂η, 1/8*(1-ξ)*(1-ζ))
    push!(∂Ψ∂η, 1/8*(1+ξ)*(1-ζ))
    push!(∂Ψ∂η, -1/8*(1-ξ)*(1+ζ))
    push!(∂Ψ∂η, -1/8*(1+ξ)*(1+ζ))
    push!(∂Ψ∂η, 1/8*(1-ξ)*(1+ζ))
    push!(∂Ψ∂η, 1/8*(1+ξ)*(1+ζ))

    push!(∂Ψ∂ζ, -1/8*(1-ξ)*(1-η))
    push!(∂Ψ∂ζ, -1/8*(1+ξ)*(1-η))
    push!(∂Ψ∂ζ, -1/8*(1-ξ)*(1+η))
    push!(∂Ψ∂ζ, -1/8*(1+ξ)*(1+η))
    push!(∂Ψ∂ζ, 1/8*(1-ξ)*(1-η))
    push!(∂Ψ∂ζ, 1/8*(1+ξ)*(1-η))
    push!(∂Ψ∂ζ, 1/8*(1-ξ)*(1+η))
    push!(∂Ψ∂ζ, 1/8*(1+ξ)*(1+η))

    ∂x∂ξ = 0.0
    ∂x∂η = 0.0
    ∂x∂ζ = 0.0
    ∂y∂ξ = 0.0
    ∂y∂η = 0.0
    ∂y∂ζ = 0.0 
    ∂z∂ξ = 0.0 
    ∂z∂η = 0.0
    ∂z∂ζ = 0.0

    for j=1:8
        ∂x∂ξ += hex_vert[j,1]*∂Ψ∂ξ[j]
        ∂x∂η += hex_vert[j,1]*∂Ψ∂η[j]
        ∂x∂ζ += hex_vert[j,1]*∂Ψ∂ζ[j]

        ∂y∂ξ += hex_vert[j,2]*∂Ψ∂ξ[j]
        ∂y∂η += hex_vert[j,2]*∂Ψ∂η[j]
        ∂y∂ζ += hex_vert[j,2]*∂Ψ∂ζ[j]

        ∂z∂ξ += hex_vert[j,3]*∂Ψ∂ξ[j]
        ∂z∂η += hex_vert[j,3]*∂Ψ∂η[j]
        ∂z∂ζ += hex_vert[j,3]*∂Ψ∂ζ[j]
    end

    J = (∂x∂ξ*∂y∂η*∂z∂ζ + ∂x∂η*∂y∂ζ*∂z∂ξ + ∂x∂ζ*∂y∂ξ*∂z∂η - 
            ∂x∂ζ*∂y∂η*∂z∂ξ - ∂x∂η*∂y∂ξ*∂z∂ζ - ∂x∂ξ*∂y∂ζ*∂z∂η)
    
    ∂ξ∂x = 1/J * (∂y∂η*∂z∂ζ - ∂y∂ζ*∂z∂η) 
    ∂η∂x = 1/J * (∂y∂ζ*∂z∂ξ - ∂y∂ξ*∂z∂ζ)
    ∂ζ∂x = 1/J * (∂y∂ξ*∂z∂η - ∂y∂η*∂z∂ξ)
    
    ∂ξ∂y = 1/J * (∂x∂ζ*∂z∂η - ∂x∂η*∂z∂ζ)
    ∂η∂y = 1/J * (∂x∂ξ*∂z∂ζ - ∂x∂ζ*∂z∂ξ)
    ∂ζ∂y = 1/J * (∂x∂η*∂z∂ξ - ∂x∂ξ*∂z∂η)

    ∂ξ∂z = 1/J * (∂x∂η*∂y∂ζ - ∂x∂ζ*∂y∂η)
    ∂η∂z = 1/J * (∂x∂ζ*∂y∂ξ - ∂x∂ξ*∂y∂ζ)
    ∂ζ∂z = 1/J * (∂x∂ξ*∂y∂η - ∂x∂η*∂y∂ξ)

    Jac[1,1]=J
    dxi[1,1]=∂x∂ξ
    dxi[2,1]=∂x∂η
    dxi[3,1]=∂x∂ζ
    dxi[4,1]=∂y∂ξ
    dxi[5,1]=∂y∂η
    dxi[6,1]=∂z∂ζ
    dxi[7,1]=∂z∂ξ
    dxi[8,1]=∂z∂η
    dxi[9,1]=∂z∂ζ

    dx[1,1]=∂ξ∂x
    dx[2,1]=∂ξ∂y
    dx[3,1]=∂ξ∂z
    dx[4,1]=∂η∂x
    dx[5,1]=∂η∂y
    dx[6,1]=∂η∂z
    dx[7,1]=∂ζ∂x
    dx[8,1]=∂ζ∂y
    dx[9,1]=∂ζ∂z
end

"""
### SummationByParts.tensor_operators
"""
function tensor_operators(p::Int, dim::Int; T=Float64)
    n = p+1
    oper_lgl = getLineSegSBPLobbato(degree=p)
    perm = sortperm(vec(SymCubatures.calcnodes(oper_lgl.cub, oper_lgl.vtx)))
    Q1 = oper_lgl.Q[:,:,1][:,perm][perm,:]
    H1 = diagm(vec(oper_lgl.w)[perm])
    D1 = inv(H1)*Q1[:,:,1]
    E1 = Q1 + Q1'
    In = I(n)
    tR = zeros(T,(n,1)); tR[end] = 1.0
    tL = zeros(T,(n,1)); tL[1] = 1.0

    H = zeros(T, (n^dim,n^dim))
    Q = zeros(T, (n^dim,n^dim,dim))
    D = zeros(T, (n^dim,n^dim,dim))
    E = zeros(T, (n^dim,n^dim,dim))
    R = zeros(T, (n^(dim-1),n^dim,2*dim))

    if dim==2
        H[:,:] = kron(H1, H1)
        Q[:,:,1] = kron(H1, Q1)
        Q[:,:,2] = kron(Q1, H1)
        D[:,:,1] = kron(In, D1)
        D[:,:,2] = kron(D1, In)
        E[:,:,1] = kron(In, E1)
        E[:,:,2] = kron(E1, In)
        R[:,:,1] = kron(In, tL') #left facet (x=-1)
        R[:,:,2] = kron(In, tR') #right facet (x=1)
        R[:,:,3] = kron(tL', In) #bottom facet (y=-1)
        R[:,:,4] = kron(tR', In) #top facet (y=1)
    elseif dim==3 
        H[:,:] = kron(H1,kron(H1, H1))
        Q[:,:,1] = kron(H1, kron(H1, Q1))
        Q[:,:,2] = kron(H1, kron(Q1, H1))
        Q[:,:,3] = kron(Q1, kron(H1, H1))
        D[:,:,1] = kron(In, kron(In, D1))
        D[:,:,2] = kron(In, kron(D1, In))
        D[:,:,3] = kron(D1, kron(In, In))
        E[:,:,1] = kron(In, kron(In, E1))
        E[:,:,2] = kron(In, kron(E1, In))
        E[:,:,3] = kron(E1, kron(In, In))
        R[:,:,1] = kron(In,kron(In, tL')) #left facet (x=-1)
        R[:,:,2] = kron(In,kron(In, tR')) #right facet (x=1)
        R[:,:,3] = kron(In,kron(tL', In)) #back facet (y=-1)
        R[:,:,4] = kron(In,kron(tR', In)) #front facet (y=1)
        R[:,:,5] = kron(tL',kron(In, In)) #bottom facet (z=-1) 
        R[:,:,6] = kron(tR',kron(In, In)) #top facet (z=1)
    end
    return H, Q, D, E, R
end

"""
### SummationByParts.square_normals_facet_nodes
"""
function square_normals_facet_nodes(x::Array{T}) where T
   
    dim=2
    n = size(x,2)
    nf = convert(Int, sqrt(n))
    N = zeros(T, (dim,nf,4)) # normal vector for each facet
    facet_node_idx = zeros(Int, (4,nf))

    N[1,:,1] .= -1.0
    N[1,:,2] .= 1.0
    N[2,:,3] .= -1.0
    N[2,:,4] .= 1.0

    facet_node_idx[1,:] = 1:nf:n
    facet_node_idx[2,:] = nf:nf:n
    facet_node_idx[3,:] = 1:nf 
    facet_node_idx[4,:] = (n+1-nf):n

    return N, facet_node_idx
end

"""
### SummationByParts.cube_normals_facet_nodes
"""
function cube_normals_facet_nodes(x::Array{T}) where T
   
    dim=3
    n = size(x,2)
    n1 = convert(Int, round(n^(1/3)))
    nf = convert(Int, round(n^(2/3)))
    N = zeros(T, (dim,nf,6)) # normal vector for each facet
    facet_node_idx = zeros(Int, (6,nf))

    N[1,:,1] .= -1.0
    N[1,:,2] .= 1.0
    N[2,:,3] .= -1.0
    N[2,:,4] .= 1.0
    N[3,:,5] .= -1.0
    N[3,:,6] .= 1.0

    facet_node_idx[1,:] = 1:n1:n
    facet_node_idx[2,:] = n1:n1:n
    facet_node_idx[3,:] = collect(Iterators.flatten([1:n1...] .+ (i-1)*nf for i in 1:n1))'
    facet_node_idx[4,:] = collect(Iterators.flatten([nf+1-n1:nf...] .+ (i-1)*nf for i in 1:n1))'
    facet_node_idx[5,:] = 1:nf 
    facet_node_idx[6,:] = n+1-nf:n

    return N, facet_node_idx
end

"""
### SummationByParts.map_tensor_operators_to_tri
"""
function map_tensor_operators_to_tri(p::Int; T=Float64)
    dim = 2
    xs = tensor_lgl_quad_nodes(p) #nodes on square
    xt = square_to_tri_map(xs)    #nodes on triangle 
    n = size(xs,2)
    nf = convert(Int, sqrt(n))

    qf = 2*(p+1)-3 
    cub_lgl, vtx_lgl = SummationByParts.Cubature.quadrature(qf, internal=false)
    perm = sortperm(vec(SymCubatures.calcnodes(cub_lgl, vtx_lgl)))
    B = SymCubatures.calcweights(cub_lgl)[perm]

    quad_vert = get_quad_vert()
    Nhat, facet_node_idx = square_normals_facet_nodes(xs)
    
    dxis = []
    dxs = []
    Js = []
    Ns = []
    for i = 1:3
        dxi = zeros(4,n)
        dx = zeros(4,n)
        J = zeros(1,n)
        for j=1:n
            metric_tri!(xs[:,j],quad_vert[i],view(dxi,:,j),view(dx,:,j),view(J,:,j))
        end
        push!(dxis, dxi)
        push!(dxs, dx)
        push!(Js, J)

        N = zeros(T,(dim,nf,4))
        for k=1:4
            N[1,:,k] = J[facet_node_idx[k,:]].*(dx[1,facet_node_idx[k,:]].*Nhat[1,:,k] .+ dx[3,facet_node_idx[k,:]].*Nhat[2,:,k])
            N[2,:,k] = J[facet_node_idx[k,:]].*(dx[2,facet_node_idx[k,:]].*Nhat[1,:,k] .+ dx[4,facet_node_idx[k,:]].*Nhat[2,:,k])
        end
        push!(Ns, N)
    end

    Hhat, Qhat, Dhat, Ehat, Rhat = tensor_operators(p, dim, T=T)
    Es = []
    for k=1:3
        E = zeros(T, (n,n,dim))
        for i=1:dim
            for j=1:4
                E[:,:,i] += Rhat[:,:,j]'*diagm(Ns[k][i,:,j].*B)*Rhat[:,:,j]
            end
        end
        push!(Es,E)
    end

    Hs = []
    Qs = []
    Ds = []
    for i=1:3
        Q = zeros(T, (n,n,dim))
        E = Es[i]
        D = zeros(T, (n,n,dim))
        H = diagm(vec(Js[i]))*Hhat 
        push!(Hs, H)
        Sx = 0.5*(diagm(vec(Js[i]).*dxs[i][1,:]) * Qhat[:,:,1] + diagm(vec(Js[i]).*dxs[i][3,:]) * Qhat[:,:,2]) - 
             0.5*(Qhat[:,:,1]' * diagm(vec(Js[i]).*dxs[i][1,:]) + Qhat[:,:,2]' * diagm(vec(Js[i]).*dxs[i][3,:]))
        Sy = 0.5*(diagm(vec(Js[i]).*dxs[i][2,:]) * Qhat[:,:,1] + diagm(vec(Js[i]).*dxs[i][4,:]) * Qhat[:,:,2]) - 
             0.5*(Qhat[:,:,1]' * diagm(vec(Js[i]).*dxs[i][2,:]) + Qhat[:,:,2]' * diagm(vec(Js[i]).*dxs[i][4,:]))
        Q[:,:,1] = Sx + 0.5.*E[:,:,1]
        Q[:,:,2] = Sy + 0.5.*E[:,:,2]
        push!(Qs, Q)
        D[:,:,1] = inv(H)*Q[:,:,1]
        D[:,:,2] = inv(H)*Q[:,:,2]
        # D[:,:,1] = diagm(vec(dxs[i][1,:]))*Dhat[:,:,1] + diagm(vec(dxs[i][3,:]))*Dhat[:,:,2]
        # D[:,:,2] = diagm(vec(dxs[i][2,:]))*Dhat[:,:,1] + diagm(vec(dxs[i][4,:]))*Dhat[:,:,2]
        push!(Ds, D)
    end 

    return Hs,Qs,Ds,Es
end

"""
### SummationByParts.map_tensor_operators_to_tet
"""
function map_tensor_operators_to_tet(p::Int; T=Float64)
    dim = 3
    xc = tensor_lgl_hex_nodes(p)   #nodes on cube
    xt = cube_to_tet_map(xc)       #nodes on tetrahedron 
    n = size(xc,2)
    n1 = round(n^(1/3))
    nf = convert(Int, round(n1^2))

    qf = 2*(p+1)-3 
    cub_lgl, vtx_lgl = SummationByParts.Cubature.quadrature(qf, internal=false)
    perm = sortperm(vec(SymCubatures.calcnodes(cub_lgl, vtx_lgl)))
    B = SymCubatures.calcweights(cub_lgl)[perm]
    B = vec(kron(B,B))

    hex_vert = get_hex_vert()
    Nhat, facet_node_idx = cube_normals_facet_nodes(xc) 

    dxis = []
    dxs = []
    Js = []
    Ns = []
    for i = 1:4
        dxi = zeros(9,n)
        dx = zeros(9,n)
        J = zeros(1,n)
        N = zeros(3,n)
        for j=1:n
            metric_tet!(xc[:,j],hex_vert[i],view(dxi,:,j),view(dx,:,j),view(J,:,j))
        end
        push!(dxis, dxi)
        push!(dxs, dx)
        push!(Js, J)

        N = zeros(T,(dim,nf,2*dim))
        for k=1:6
            for id=1:dim
                N[id,:,k] = J[facet_node_idx[k,:]].*(dx[id,facet_node_idx[k,:]].*Nhat[1,:,k] .+ 
                                                     dx[dim+id,facet_node_idx[k,:]].*Nhat[2,:,k] .+
                                                     dx[2*dim+id,facet_node_idx[k,:]].*Nhat[3,:,k])
            end
        end
        push!(Ns, N)
    end

    Hhat, Qhat, Dhat, Ehat, Rhat = tensor_operators(p, dim, T=T)
    Es = []
    for k=1:dim+1
        E = zeros(T, (n,n,dim))
        for i=1:dim
            for j=1:2*dim
                E[:,:,i] += Rhat[:,:,j]'*diagm(Ns[k][i,:,j].*B)*Rhat[:,:,j]
            end
        end
        push!(Es,E)
    end

    Hs = []
    Qs = []
    Ds = []
    for i=1:dim+1
        Q = zeros(T, (n,n,dim))
        E = Es[i]
        D = zeros(T, (n,n,dim))
        H = diagm(vec(Js[i]))*Hhat 
        push!(Hs, H)
        Sx = 0.5*(diagm(vec(Js[i]).*dxs[i][1,:]) * Qhat[:,:,1] + diagm(vec(Js[i]).*dxs[i][1+dim,:]) * Qhat[:,:,2] + diagm(vec(Js[i]).*dxs[i][1+2*dim,:]) * Qhat[:,:,3]) - 
             0.5*(Qhat[:,:,1]' * diagm(vec(Js[i]).*dxs[i][1,:]) + Qhat[:,:,2]' * diagm(vec(Js[i]).*dxs[i][1+dim,:]) + Qhat[:,:,3]' * diagm(vec(Js[i]).*dxs[i][1+2*dim,:]))
        Sy = 0.5*(diagm(vec(Js[i]).*dxs[i][2,:]) * Qhat[:,:,1] + diagm(vec(Js[i]).*dxs[i][2+dim,:]) * Qhat[:,:,2] + diagm(vec(Js[i]).*dxs[i][2+2*dim,:]) * Qhat[:,:,3]) - 
             0.5*(Qhat[:,:,1]' * diagm(vec(Js[i]).*dxs[i][2,:]) + Qhat[:,:,2]' * diagm(vec(Js[i]).*dxs[i][2+dim,:]) + Qhat[:,:,3]' * diagm(vec(Js[i]).*dxs[i][2+2*dim,:]))
        Sz = 0.5*(diagm(vec(Js[i]).*dxs[i][3,:]) * Qhat[:,:,1] + diagm(vec(Js[i]).*dxs[i][3+dim,:]) * Qhat[:,:,2] + diagm(vec(Js[i]).*dxs[i][3+2*dim,:]) * Qhat[:,:,3]) - 
             0.5*(Qhat[:,:,1]' * diagm(vec(Js[i]).*dxs[i][3,:]) + Qhat[:,:,2]' * diagm(vec(Js[i]).*dxs[i][3+dim,:]) + Qhat[:,:,3]' * diagm(vec(Js[i]).*dxs[i][3+2*dim,:]))
        
        Q[:,:,1] = Sx + 0.5.*E[:,:,1]
        Q[:,:,2] = Sy + 0.5.*E[:,:,2]
        Q[:,:,3] = Sz + 0.5.*E[:,:,3]
        push!(Qs, Q)
        D[:,:,1] = inv(H)*Q[:,:,1]
        D[:,:,2] = inv(H)*Q[:,:,2]
        D[:,:,3] = inv(H)*Q[:,:,3]
        # D[:,:,1] = diagm(vec(dxs[i][1,:]))*Dhat[:,:,1] + diagm(vec(dxs[i][4,:]))*Dhat[:,:,2] + diagm(vec(dxs[i][7,:]))*Dhat[:,:,3]
        # D[:,:,2] = diagm(vec(dxs[i][2,:]))*Dhat[:,:,1] + diagm(vec(dxs[i][5,:]))*Dhat[:,:,2] + diagm(vec(dxs[i][8,:]))*Dhat[:,:,3]
        # D[:,:,3] = diagm(vec(dxs[i][3,:]))*Dhat[:,:,1] + diagm(vec(dxs[i][6,:]))*Dhat[:,:,2] + diagm(vec(dxs[i][9,:]))*Dhat[:,:,3]
        push!(Ds, D)
    end 

    return Hs,Qs,Ds,Es
end