using SummationByParts
using SummationByParts.Cubature
using SummationByParts.OrthoPoly
using SparseArrays

# This file gathers functions that are used to construct multidimensional 
# SBP operators by splitting simplices into tensor-product quad and hexs

"""
### SummationByParts.tensor_lgl_quad_nodes

Computes tensor-product nodes on a quadrilateral 
**Inputs**
* `p`: degree of the operator 
* `opertype`: the type of 1d operator used to construct the SST-SBP method (lgl or csbp)
* `n1d`: number of nodes in the 1D operator 

**Outputs** 
* `xy`: node coordinates
* `w`: weights of the 1D operator 
"""
function tensor_quad_nodes(p::Int;opertype::String="lgl", n1d::Int=-1)
    if opertype=="lgl"
        qf = 2*(p+1)-3 
        cub, vtx = SummationByParts.Cubature.quadrature(qf, internal=false)
    elseif opertype=="csbp"
        qf = 2p
        if n1d==-1
            n1d=4*p 
        end
        cub, vtx = SummationByParts.Cubature.getLineCubatureGregory(qf, n1d)
    else 
        error("Operator not implemented. Should choose between 'lgl' and 'csbp'.")
    end
    x = vec(SymCubatures.calcnodes(cub, vtx))
    perm = sortperm(x)
    x = x[perm]'
    y = x[1] .*ones(1,cub.numnodes)
    xy = vcat(x,y)
    for i = 1:cub.numnodes
        xy = hcat(xy, vcat(x,x[i].*ones(1,cub.numnodes)))
    end
    w1 = (vec(SymCubatures.calcweights(cub)))[perm]
    return xy[:,cub.numnodes+1:end], w1
end

"""
### SummationByParts.tensor_lgl_hex_nodes

Computes tensor-product nodes on a hexahedron
**Inputs**
* `p`: degree of the operator 
* `opertype`: the type of 1d operator used to construct the SST-SBP method (lgl or csbp)
* `n1d`: number of nodes in the 1D operator 

**Outputs** 
* `xyz`: node coordinates
* `w`: weights of the 1D operator 
"""
function tensor_hex_nodes(p::Int; opertype::String="lgl", n1d::Int=-1)
    # qf = 2*(p+1)-3 
    # cub, vtx = SummationByParts.Cubature.quadrature(qf, internal=false)
    if opertype=="lgl"
        qf = 2*(p+1)-3 
        cub, vtx = SummationByParts.Cubature.quadrature(qf, internal=false)
    elseif opertype=="csbp"
        qf = 2p
        if n1d==-1
            n1d=4*p 
        end
        cub, vtx = SummationByParts.Cubature.getLineCubatureGregory(qf, n1d)
    else 
        error("Operator not implemented. Should choose between 'lgl' and 'csbp'.")
    end
    # x = sort(vec(SymCubatures.calcnodes(cub, vtx)))'
    x = vec(SymCubatures.calcnodes(cub, vtx))
    perm = sortperm(x)
    x = x[perm]'
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
    w1 = (vec(SymCubatures.calcweights(cub)))[perm]
    return xyz[:,cub.numnodes+1:end], w1
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
function tensor_operators(p::Int, dim::Int; opertype::String="lgl", n1d::Int=-1, T=Float64)
    # n = p+1
    # oper_lgl = getLineSegSBPLobbato(degree=p)
    # perm = sortperm(vec(SymCubatures.calcnodes(oper_lgl.cub, oper_lgl.vtx)))
    # Q1 = oper_lgl.Q[:,:,1][:,perm][perm,:]
    # H1 = diagm(vec(oper_lgl.w)[perm])

    if opertype=="lgl"
        n = p+1
        oper_lgl = getLineSegSBPLobbato(degree=p)
        perm = sortperm(vec(SymCubatures.calcnodes(oper_lgl.cub, oper_lgl.vtx)))
        Q1 = oper_lgl.Q[:,:,1][:,perm][perm,:]
        H1 = diagm(vec(oper_lgl.w)[perm])
        D1 = inv(H1)*Q1[:,:,1]
        E1 = Q1 + Q1'
    elseif opertype=="csbp"
        if n1d==-1
            n1d=4*p 
        end
        n=n1d
        H1,Q1,E1,D1 = build_csbp_operators(p, n1d)
    else 
        error("Operator not implemented. Should choose between 'lgl' and 'csbp'.")
    end

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
### SummationByParts.normals_square
"""
function normals_square(x::Array{T}) where T
   
    dim=2
    n = size(x,2)
    nf = convert(Int, sqrt(n))
    N = zeros(T, (dim,nf,4)) # normal vector for each facet

    N[1,:,1] .= -1.0
    N[1,:,2] .= 1.0
    N[2,:,3] .= -1.0
    N[2,:,4] .= 1.0

    return N
end

"""
### SummationByParts.facet_nodes_square
"""
function facet_nodes_square(x::Array{T}) where T
    n = size(x,2)
    nf = convert(Int, sqrt(n))
    facet_node_idx = zeros(Int, (4,nf))

    facet_node_idx[1,:] = 1:nf:n
    facet_node_idx[2,:] = nf:nf:n
    facet_node_idx[3,:] = 1:nf 
    facet_node_idx[4,:] = (n+1-nf):n

    return facet_node_idx
end

"""
### SummationByParts.normals_cube
"""
function normals_cube(x::Array{T}) where T
   
    dim=3
    n = size(x,2)
    nf = convert(Int, round(n^(2/3)))
    N = zeros(T, (dim,nf,6)) # normal vector for each facet

    N[1,:,1] .= -1.0
    N[1,:,2] .= 1.0
    N[2,:,3] .= -1.0
    N[2,:,4] .= 1.0
    N[3,:,5] .= -1.0
    N[3,:,6] .= 1.0

    return N
end

"""
### SummationByParts.cube_normals_facet_nodes
"""
function facet_nodes_cube(x::Array{T}) where T
 
    n = size(x,2)
    n1 = convert(Int, round(n^(1/3)))
    nf = convert(Int, round(n^(2/3)))
    facet_node_idx = zeros(Int, (6,nf))

    facet_node_idx[1,:] = 1:n1:n
    facet_node_idx[2,:] = n1:n1:n
    facet_node_idx[3,:] = collect(Iterators.flatten([1:n1...] .+ (i-1)*nf for i in 1:n1))'
    facet_node_idx[4,:] = collect(Iterators.flatten([nf+1-n1:nf...] .+ (i-1)*nf for i in 1:n1))'
    facet_node_idx[5,:] = 1:nf 
    facet_node_idx[6,:] = n+1-nf:n

    return facet_node_idx
end

"""
### SummationByParts.map_tensor_operators_to_tri
"""
function map_tensor_operators_to_tri(p::Int; opertype::String="lgl", n1d::Int=-1, T=Float64)
    dim = 2
    xs, B = tensor_quad_nodes(p,opertype=opertype,n1d=n1d) #nodes on square
    n = size(xs,2)
    nf = convert(Int, sqrt(n))

    # qf = 2*(p+1)-3 
    # cub_lgl, vtx_lgl = SummationByParts.Cubature.quadrature(qf, internal=false)
    # perm = sortperm(vec(SymCubatures.calcnodes(cub_lgl, vtx_lgl)))
    # B = SymCubatures.calcweights(cub_lgl)[perm]
    
    quad_vert = get_quad_vert()
    Nhat = normals_square(xs)
    facet_node_idx = facet_nodes_square(xs)
    
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

    Hhat, Qhat, Dhat, Ehat, Rhat = tensor_operators(p, dim, opertype=opertype, n1d=n1d, T=T)
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

    return Hs,Qs,Ds,Es,Ns
end

"""
### SummationByParts.map_tensor_operators_to_tet
"""
function map_tensor_operators_to_tet(p::Int; opertype::String="lgl", n1d::Int=-1, T=Float64)
    dim = 3
    xc, B = tensor_hex_nodes(p, opertype=opertype, n1d=n1d)   #nodes on cube
    n = size(xc,2)
    n1 = round(n^(1/3))
    nf = convert(Int, round(n1^2))

    # qf = 2*(p+1)-3 
    # cub_lgl, vtx_lgl = SummationByParts.Cubature.quadrature(qf, internal=false)
    # perm = sortperm(vec(SymCubatures.calcnodes(cub_lgl, vtx_lgl)))
    # B = SymCubatures.calcweights(cub_lgl)[perm]
    B = vec(kron(B,B))

    hex_vert = get_hex_vert()
    Nhat = normals_cube(xc) 
    facet_node_idx = facet_nodes_cube(xc) 

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

    Hhat, Qhat, Dhat, Ehat, Rhat = tensor_operators(p, dim, opertype=opertype, n1d=n1d, T=T)
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

    return Hs,Qs,Ds,Es,Ns
end

"""
### SummationByParts.global_node_index_tri
"""
function global_node_index_tri(p::Int;opertype::String="lgl",n1d::Int=-1, T=Float64)
    xs,_ = tensor_quad_nodes(p,opertype=opertype,n1d=n1d)
    xt = square_to_tri_map(xs)
    n = size(xs,2)

    facet_node_idx = facet_nodes_square(xs)
    xg = copy(xt)
    remove_idx = collect(Iterators.flatten([n.+facet_node_idx[1,:], (2*n).+facet_node_idx[2,:], (2*n).+facet_node_idx[3,:]]))
    xg = xg[:, filter(x -> !(x in remove_idx), 1:size(xg, 2))]

    loc_glob_idx = []
    for k = 1:3
        x = zeros(Int,(2,n))
        x[1,:] = 1:n
        for i=1:n
            xgidx = xg .- xt[:,(k-1)*n+i]
            col_norm = [norm(xgidx[:, j]) for j in 1:size(xg,2)]
            ig = argmin(col_norm)
            x[2,i] = ig
        end
        push!(loc_glob_idx, x)
    end
    return xg, loc_glob_idx
end

"""
### SummationByParts.global_node_index_tet
"""
function global_node_index_tet(p::Int; opertype::String="lgl",n1d::Int=-1, T=Float64)
    xh,_ = tensor_hex_nodes(p,opertype=opertype,n1d=n1d)
    xt = cube_to_tet_map(xh)
    n = size(xh,2)

    facet_node_idx = facet_nodes_cube(xh)
    xg = copy(xt)
    remove_idx = collect(Iterators.flatten([n.+facet_node_idx[1,:], 
                                            (2*n).+facet_node_idx[2,:],
                                            (2*n).+facet_node_idx[3,:],
                                            (3*n).+facet_node_idx[1,:],
                                            (3*n).+facet_node_idx[4,:],
                                            (3*n).+facet_node_idx[6,:]]))
    xg = xg[:, filter(x -> !(x in remove_idx), 1:size(xg, 2))]

    loc_glob_idx = []
    for k = 1:4
        x = zeros(Int,(2,n))
        x[1,:] = 1:n
        for i=1:n
            xgidx = xg .- xt[:,(k-1)*n+i]
            col_norm = [norm(xgidx[:, j]) for j in 1:size(xg,2)]
            ig = argmin(col_norm)
            x[2,i] = ig
        end
        push!(loc_glob_idx, x)
    end
    return xg, loc_glob_idx
end

"""
### SummationByParts.construct_zmatrix
"""
function construct_zmatrix(glob_idx::Array{T},i::Int,j::Int,nglob::Int) where T
    ihat = glob_idx[2,i]
    jhat = glob_idx[2,j]
    ei = sparse(I,nglob,nglob)[:, ihat]
    ej = sparse(I,nglob,nglob)[:, jhat]
    Z = ei*ej'
    return Z
end

"""
### SummationByParts.construct_split_operator_tri
"""
function construct_split_operator_tri(p::Int; opertype::String="lgl", n1d::Int=-1, T=Float64)
    Hs,Qs,Ds,Es,_ = map_tensor_operators_to_tri(p, opertype=opertype, n1d=n1d)
    xg, loc_glob_idx = global_node_index_tri(p, opertype=opertype, n1d=n1d)
    n = size(Hs[1],1)
    nglob = size(xg,2) 
    dim = 2

    H = spzeros(nglob,nglob)
    Q = [spzeros(nglob,nglob),spzeros(nglob,nglob)]
    D = [spzeros(nglob,nglob),spzeros(nglob,nglob)]
    E = [spzeros(nglob,nglob),spzeros(nglob,nglob)]
    for k=1:dim+1
        for i=1:n
            for j=1:n 
                glob_idx = loc_glob_idx[k]
                Z = construct_zmatrix(glob_idx,i,j,nglob)
                if Hs[k][i,j]!=0.0
                    H[:,:] += (Hs[k][i,j]*Z)
                end
                for id=1:dim
                    if Qs[k][i,j,id]!=0.0
                        Q[id] += (Qs[k][i,j,id]*Z)
                    end
                    if Es[k][i,j,id]!=0.0
                        E[id] += (Es[k][i,j,id]*Z)
                    end
                end
            end
        end
    end
    
    for id=1:dim 
        D[id] = inv(Matrix(H))*Q[id]
    end
    Q = [Matrix(m) for m in Q]
    D = [Matrix(m) for m in D]
    E = [Matrix(m) for m in E]
    return Matrix(H),Q,D,E
end

"""
### SummationByParts.construct_split_operator_tet
"""
function construct_split_operator_tet(p::Int; opertype::String="lgl", n1d::Int=-1, T=Float64)
    Hs,Qs,Ds,Es,_ = map_tensor_operators_to_tet(p, opertype=opertype, n1d=n1d)
    xg, loc_glob_idx = global_node_index_tet(p, opertype=opertype, n1d=n1d)
    n = size(Hs[1],1)
    nglob = size(xg,2) 
    dim = 3

    H = spzeros(nglob,nglob)
    Q = [spzeros(nglob,nglob),spzeros(nglob,nglob),spzeros(nglob,nglob)]
    D = [spzeros(nglob,nglob),spzeros(nglob,nglob),spzeros(nglob,nglob)]
    E = [spzeros(nglob,nglob),spzeros(nglob,nglob),spzeros(nglob,nglob)]
    for k=1:dim+1
        for i=1:n
            for j=1:n 
                glob_idx = loc_glob_idx[k]
                Z = construct_zmatrix(glob_idx,i,j,nglob)
                if Hs[k][i,j]!=0.0
                    H[:,:] += (Hs[k][i,j]*Z)
                end
                for id=1:dim
                    if Qs[k][i,j,id]!=0.0
                        Q[id] += (Qs[k][i,j,id]*Z)
                    end
                    if Es[k][i,j,id]!=0.0 
                        E[id] += (Es[k][i,j,id]*Z)
                    end
                end
            end
        end
    end
    
    for id=1:dim 
        D[id] = inv(Matrix(H))*Q[id]
    end
    Q = [Matrix(m) for m in Q]
    D = [Matrix(m) for m in D]
    E = [Matrix(m) for m in E]
    return Matrix(H),Q,D,E
end

"""
### SummationByParts.global_node_index_tri_facet
"""
function global_node_index_tri_facet(p::Int;opertype::String="lgl",n1d::Int=-1, T=Float64)
    dim = 2
    xs,_ = tensor_quad_nodes(p, opertype=opertype, n1d=n1d)
    xt = square_to_tri_map(xs)
    xg,_ = global_node_index_tri(p, opertype=opertype, n1d=n1d)
    n = size(xs,2)
    n1 = convert(Int, round(n^(1/dim)))
    nf = n1^(dim-1)

    facet_node_idx = facet_nodes_square(xs)
    xf = zeros(T, (dim, dim*nf-1, dim+1))
    keep_idx = []
    push!(keep_idx, collect(Iterators.flatten([(n).+facet_node_idx[2,:], ((2*n).+facet_node_idx[4,1:end-1])])))
    push!(keep_idx, collect(Iterators.flatten([((2*n).+facet_node_idx[1,:]), (facet_node_idx[1,1:end-1])])))
    push!(keep_idx, collect(Iterators.flatten([facet_node_idx[3,:], (n).+facet_node_idx[3,2:end]])))
    for k=1:3
        xf[:,:,k] = xt[:,keep_idx[k]]
    end
    xfall= zeros(T,(dim,nf,dim*(dim+1)))
    xfall[:,:,1] = xt[:,(n).+facet_node_idx[2,:]]
    xfall[:,:,2] = xt[:,((2*n).+facet_node_idx[4,:])]
    xfall[:,:,3] = xt[:,((2*n).+facet_node_idx[1,:])]
    xfall[:,:,4] = xt[:,(facet_node_idx[1,:])]
    xfall[:,:,5] = xt[:,facet_node_idx[3,:]]
    xfall[:,:,6] = xt[:,(n).+facet_node_idx[3,:]]

    loc_glob_facet_idx = []
    for k = 1:(dim*(dim+1))
        x = zeros(Int,(3,nf))
        x[1,:] = 1:nf
        for i=1:nf
            xfidx = xf[:,:,convert(Int,ceil(k/dim))] .- xfall[:,i,k] 
            col_norm = [norm(xfidx[:, j]) for j in 1:size(xf,2)]
            ig = argmin(col_norm)
            x[2,i] = ig
        end
        for i=1:nf
            xgidx = xg .- xfall[:,i,k]
            col_norm = [norm(xgidx[:, j]) for j in 1:size(xg,2)]
            ig = argmin(col_norm)
            x[3,i] = ig
        end
        push!(loc_glob_facet_idx, x)
    end
    return xf, loc_glob_facet_idx
end

"""
### SummationByParts.global_node_index_tet_facet
"""
function global_node_index_tet_facet(p::Int;opertype::String="lgl",n1d::Int=-1, T=Float64)
    dim = 3
    xh,_= tensor_hex_nodes(p, opertype=opertype, n1d=n1d)
    xt = cube_to_tet_map(xh)
    xg,_ = global_node_index_tet(p, opertype=opertype, n1d=n1d)
    n = size(xh,2)
    n1 = convert(Int, round(n^(1/dim)))
    nf = n1^(dim-1)

    facet_node_idx = facet_nodes_cube(xh)
    xf = zeros(T, (dim, dim*nf-3*n1+1, dim+1))
    keep_idx = []
    push!(keep_idx, collect(Iterators.flatten([facet_node_idx[3,:],
                                               (n.+facet_node_idx[3,:]), 
                                               ((3*n).+facet_node_idx[5,:])])))
    push!(keep_idx, collect(Iterators.flatten([(n.+facet_node_idx[2,:]), 
                                                (2*n).+(facet_node_idx[4,:]),
                                                (3*n).+(facet_node_idx[2,:])])))
    push!(keep_idx, collect(Iterators.flatten([(facet_node_idx[1,:]), 
                                                (2*n).+(facet_node_idx[1,:]),
                                                (3*n).+(facet_node_idx[3,:])])))
    push!(keep_idx, collect(Iterators.flatten([(facet_node_idx[5,:]), 
                                                (n).+(facet_node_idx[5,:]),
                                                (2*n).+(facet_node_idx[5,:])])))

    unique_idx = [[keep_idx[1][1]],[keep_idx[2][1]],[keep_idx[3][1]],[keep_idx[4][1]]]
    for k=1:dim+1
        xtf = xt[:,keep_idx[k]]
        for i=2:nf*dim
            xtf_temp = xtf .- xtf[:,i]
            col_norm = [norm(xtf_temp[:, j]) for j in 1:size(xtf,2)]
            idx = argmin(col_norm[1:(i-1)])
            if norm(xtf_temp[:,idx]) > 1e-14
                push!(unique_idx[k],keep_idx[k][i])
            end
        end
    end

    for k=1:dim+1
        xf[:,:,k] = xt[:,unique_idx[k]]
    end
    xfall= zeros(T,(dim,nf,dim*(dim+1)))
    xfall[:,:,1] = xt[:,facet_node_idx[3,:]]
    xfall[:,:,2] = xt[:,(n.+facet_node_idx[3,:])]
    xfall[:,:,3] = xt[:,((3*n).+facet_node_idx[5,:])]
    xfall[:,:,4] = xt[:,(n.+facet_node_idx[2,:])]
    xfall[:,:,5] = xt[:,(2*n).+(facet_node_idx[4,:])]
    xfall[:,:,6] = xt[:,(3*n).+(facet_node_idx[2,:])]
    xfall[:,:,7] = xt[:,(facet_node_idx[1,:])]
    xfall[:,:,8] = xt[:,(2*n).+(facet_node_idx[1,:])]
    xfall[:,:,9] = xt[:,(3*n).+(facet_node_idx[3,:])]
    xfall[:,:,10] = xt[:,(facet_node_idx[5,:])]
    xfall[:,:,11] = xt[:,(n).+(facet_node_idx[5,:])]
    xfall[:,:,12] = xt[:,(2*n).+(facet_node_idx[5,:])]

    loc_glob_facet_idx = []
    for k = 1:(dim*(dim+1))
        x = zeros(Int,(3,nf))
        x[1,:] = 1:nf
        for i=1:nf
            xfidx = xf[:,:,convert(Int,ceil(k/dim))] .- xfall[:,i,k] 
            col_norm = [norm(xfidx[:, j]) for j in 1:size(xf,2)]
            ig = argmin(col_norm)
            x[2,i] = ig
        end
        for i=1:nf
            xgidx = xg .- xfall[:,i,k]
            col_norm = [norm(xgidx[:, j]) for j in 1:size(xg,2)]
            ig = argmin(col_norm)
            x[3,i] = ig
        end
        push!(loc_glob_facet_idx, x)
    end
    return xf, loc_glob_facet_idx
end

"""
### SummationByParts.construct_split_facet_operator_tri
"""
function construct_split_facet_operator_tri(p::Int; opertype::String="lgl", n1d::Int=-1, T=Float64)
    dim = 2 
    # qf = 2*(p+1)-3 
    # cub_lgl, vtx_lgl = SummationByParts.Cubature.quadrature(qf, internal=false)
    # perm = sortperm(vec(SymCubatures.calcnodes(cub_lgl, vtx_lgl)))
    # Bhat = SymCubatures.calcweights(cub_lgl)[perm]
    _, Bhat = tensor_quad_nodes(p,opertype=opertype,n1d=n1d)

    _,_,_,_,Nhat = map_tensor_operators_to_tri(p, opertype=opertype, n1d=n1d)
    nf = length(Bhat)
    n = (dim+1)*nf^dim - (dim+dim^(dim-2))*nf^(dim-1) + (dim-2)*((dim+1)*nf - 2) + 1
    # N_idx = [2 2; 3 1; 1 3] #first column contains element number, and second column contains facet number of the element

    xf, loc_glob_idx = global_node_index_tri_facet(p, opertype=opertype, n1d=n1d)
    nglob = size(xf,2) 
    B = zeros(T,(nglob,nglob,dim+1))
    for k=1:dim+1
        for i=1:dim
            glob_idx = loc_glob_idx[(k-1)*dim+i]
            for j=1:nf
                Z = construct_zmatrix(glob_idx,j,j,nglob)
                B[:,:,k] += (Bhat[j]*Z)
            end
        end
    end
    
    # N = ones(T, (dim, size(xf,2),dim+1))
    # for k=1:dim+1
    #     for i=1:dim 
    #         N[i,:,k] = Nhat[N_idx[k,1]][i,1,N_idx[k,2]] * N[i,:,k]
    #     end
    # end

    N_idx = [[2 2; 3 4],[3 1; 1 1],[1 3; 2 3]]
    N = zeros(T, (dim, nglob,dim+1))
    for id=1:dim
        for k=1:dim+1
            jj = dim*(k-1)
            idx_vec=(collect(Iterators.flatten([loc_glob_idx[jj+1][2,:],loc_glob_idx[jj+2][2,:]])))
            idx = [findfirst(isequal(num), idx_vec) for num in 1:nglob]
            N[id,:,k] = collect(Iterators.flatten([Nhat[N_idx[k][1,1]][id,:,N_idx[k][1,2]],
                                                   Nhat[N_idx[k][2,1]][id,:,N_idx[k][2,2]]]))[idx]
        end
    end

    R = zeros(T, (nglob,n,dim+1))
    for k=1:dim+1 
        jj = dim*(k-1)
        loc_idx = unique(collect(Iterators.flatten([loc_glob_idx[jj+1][2,:],loc_glob_idx[jj+2][2,:]])))
        glob_idx = unique(collect(Iterators.flatten([loc_glob_idx[jj+1][3,:],loc_glob_idx[jj+2][3,:]])))
        for i=1:nglob
            R[loc_idx[i],glob_idx[i],k] = 1.0
        end 
    end

    E = zeros(T, (n,n,dim))
    for i=1:dim 
        for k=1:dim+1
            E[:,:,i] += R[:,:,k]'*diagm(N[i,:,k])*B[:,:,k]*R[:,:,k]
        end
    end
    return B, N, R, E
end

"""
### SummationByParts.construct_split_facet_operator_tet
"""
function construct_split_facet_operator_tet(p::Int; opertype::String="lgl", n1d::Int=-1, T=Float64)
    dim = 3 
    # qf = 2*(p+1)-3 
    # cub_lgl, vtx_lgl = SummationByParts.Cubature.quadrature(qf, internal=false)
    # perm = sortperm(vec(SymCubatures.calcnodes(cub_lgl, vtx_lgl)))
    # B1 = SymCubatures.calcweights(cub_lgl)[perm]
    _, B1 = tensor_quad_nodes(p,opertype=opertype,n1d=n1d)

    Bhat = kron(diagm(B1),diagm(B1))
    _,_,_,_,Nhat = map_tensor_operators_to_tet(p,opertype=opertype,n1d=n1d)
    n1 = length(B1)
    nf = n1^(dim-1)
    n = (dim+1)*n1^dim - (dim+dim^(dim-2))*n1^(dim-1) + (dim-2)*((dim+1)*n1 - 2) + 1
    # N_idx = [1 3; 2 2; 3 1; 1 5] #first column contains element number, and second column contains facet number of the element
    # N_idx = [[1 3; 2 3; 4 5],[2 2; 3 4; 4 2],[1 1; 3 1; 4 3],[1 5; 2 5; 3 5]]
    N_idx = [[1 3; 2 3; 4 5],[2 2; 3 4; 4 2],[1 1; 3 1; 4 3],[1 5; 2 5; 3 5]]
    
    xf, loc_glob_idx = global_node_index_tet_facet(p,opertype=opertype,n1d=n1d)
    nglob = size(xf,2) 
    B = zeros(T,(nglob,nglob,dim+1))
    N = zeros(T, (dim, nglob,dim+1))
    for k=1:dim+1
        for i=1:dim
            glob_idx = loc_glob_idx[(k-1)*dim+i]
            for j=1:nf
                Z = construct_zmatrix(glob_idx,j,j,nglob)
                B[:,:,k] += (Bhat[j,j]*Z)
            end
        end
    end
    
    for id=1:dim
        for k=1:dim+1
            jj = dim*(k-1)
            idx_vec=(collect(Iterators.flatten([loc_glob_idx[jj+1][2,:],loc_glob_idx[jj+2][2,:],loc_glob_idx[jj+3][2,:]])))
            idx = [findfirst(isequal(num), idx_vec) for num in 1:nglob]
            N[id,:,k] = collect(Iterators.flatten([Nhat[N_idx[k][1,1]][id,:,N_idx[k][1,2]],
                                                   Nhat[N_idx[k][2,1]][id,:,N_idx[k][2,2]],
                                                   Nhat[N_idx[k][3,1]][id,:,N_idx[k][3,2]]]))[idx]
        end
    end
    # N = ones(T, (dim, size(xf,2),dim+1))
    # for k=1:dim+1
    #     for i=1:dim 
    #         N[i,:,k] = Nhat[N_idx[k,1]][i,:,N_idx[k,2]] * N[i,:,k]
    #     end
    # end

    R = zeros(T, (nglob,n,dim+1))
    for k=1:dim+1 
        jj = dim*(k-1)
        loc_idx = unique(collect(Iterators.flatten([loc_glob_idx[jj+1][2,:],loc_glob_idx[jj+2][2,:],loc_glob_idx[jj+3][2,:]])))
        glob_idx = unique(collect(Iterators.flatten([loc_glob_idx[jj+1][3,:],loc_glob_idx[jj+2][3,:],loc_glob_idx[jj+3][3,:]])))
        for i=1:nglob
            R[loc_idx[i],glob_idx[i],k] = 1.0
        end 
    end

    # xg, _ = global_node_index_tet(p)
    E = zeros(T, (n,n,dim))
    for i=1:dim 
        for k=1:dim+1
            E[:,:,i] += R[:,:,k]'*diagm(N[i,:,k])*B[:,:,k]*R[:,:,k]
        end
    end
    return B, N, R, E
end

"""
### SummationByParts.stat_sst_sbp

Returns DOF and numbor of nonzero elements for SST SBP operators

**Inputs**
* `p`: degree of the operator 
* `dim`: dimesion (2d or 3d)
* `opertype`: the type of 1d operator used to construct the SST-SBP method (lgl or csbp)
* `n1d`: number of nodes in the 1D operator 

**Outputs** 
* `dof`: degrees of freedom in the element (number of nodes)
* `nnz`: number of nonzero elements in the derivative matrix
"""
function stat_sst_sbp(p::Int,dim::Int; opertype::String="lgl", n1d::Int=-1)
    if opertype=="lgl"
        n1=p+1
        nn = n1
    elseif opertype=="csbp"
        if n1d==-1
            n1d=4*p 
        end
        n1=n1d 
        nn = 2*p
    else 
        error("Operator not implemented. Should choose between 'lgl' and 'csbp'.")
    end
    # n1 = p+1
    dof = (dim+1)*n1^dim - (dim+dim^(dim-2))*n1^(dim-1) + (dim-2)*((dim+1)*n1 - 2) + 1
    # nnz = (dim*p + 1)*dof
    if opertype=="lgl"
        nnz = (nn+(nn-1)*(dim-1))*dof
    elseif opertype=="csbp"
        nnz = (nn+(nn-1)*(dim-1))*dof #(4*p*(dim+1)*(2*dim)) + (nn+(nn-1)*(dim-1))*(dof-(4*p*(dim+1)*(2*dim)))
    end
    return dof, nnz
end

"""
### SummationByParts.csbp_interior_operators
"""
function csbp_interior_operators(p::Int; T=Float64)
    s = convert(Int,p/2) #+ mod(p,2)
    m = 2*s+1
    r = s+2 #convert(Int, ceil(s/2))+2
    Qint = zeros(T,(1,m))
    if p <= 2
        Qint[r:end] = T[1/2]
    elseif p<=4
        Qint[r:end] = T[8,-1]./12
    elseif p<=6
        Qint[r:end] = T[45,-9,1]./60
    elseif p<=8
        Qint[r:end] = T[672,-168,32,-3]./840
    elseif p<=10
        Qint[r:end] = T[2100,-600,150,-25,2]./2520
    end

    Qint[1:r-2] = -reverse(Qint[r:end])
    return Qint
end
"""
### SummationByParts.build_csbp_operators
"""
function build_csbp_operators(p::Int,n::Int; T=Float64)
    @assert(n >= 2*(2*p))
    q = 2*p

    cub,vtx = SummationByParts.Cubature.getLineCubatureGregory(q,n)
    x = SymCubatures.calcnodes(cub,vtx)
    w = SymCubatures.calcweights(cub)
    perm = sortperm(vec(x))
    x = x[perm]
    w = w[perm]

    s = p
    Q = zeros(T, (n,n))
    Qint = csbp_interior_operators(2*p)
    for i=1:n 
        if i>s && i<=n-s
            Q[i,i-s:i+s] = Qint 
        end
    end

    V, Vdx = OrthoPoly.vandermonde(p, x)
    r = 2*s
    Qhat = Q[1:r,1:r]
    Hhat = diagm(w[1:r])
    Ehat = zeros(T,(n,n))[1:r,1:r]; Ehat[1,1]=-1.0; 
    Dhat = (diagm(1.0./w)*Q)[1:r,r+1:r+s]
    Vhat = V[1:r,:]
    Vdxhat = Vdx[1:r,:]
    Ahat = (Vdxhat - Dhat*V[r+1:r+s,:])

    Qhat = pocs(Qhat,Hhat,Ehat,Vhat,Ahat)
    Q[1:r,1:r] = Qhat 
    Q[end-r+1:end,end-r+1:end] = -reverse(Qhat)
    H = diagm(w)
    E = zeros(T,(n,n)); E[1,1]=-1.0; E[end,end]=1.0
    D = inv(H)*Q
    return H,Q,E,D
end

"""
### SummationByParts.pocs
"""
function pocs(Q::Array{T,2},H::Array{T,2},E::Array{T,2},V::Array{T,2},A::Array{T,2}) where T
    tol = 4e-14
    err1 = 1.0
    err2 = 1.0
  
    while ((err1 > tol) || (err2 > tol))
        Q = 0.5.*(E + Q - Q')
        Q = Q + (H*A - Q*V)*pinv(V) 
        
        err1 = norm(Q + Q' - E)
        err2 = norm(Q*V - H*A)
        # println("err1: ", err1, "   ", "err2: ", err2)
    end
  
    return Q
end