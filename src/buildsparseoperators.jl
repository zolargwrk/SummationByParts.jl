using Interpolations
using SummationByParts
using SummationByParts.Cubature
using SummationByParts.OrthoPoly
using DataStructures
#using SymPy

# This file gathers together functions used to build the SBP operators

"""
### SummationByParts.sparse_stencil

Identify the stencil of sparse SBP operators

**Inputs**

* `cub`: symmetric cubature rule
* `p`: maximum total degree for the Proriol polynomials

**Outputs**

* `stencil_lines`: a matrix containing the stencil in each direction for each node

"""

function sparse_stencil(cub::TriSymCub{T}, p::Int) where {T}

    # use equilateral triangle
    vtx = T[0 0; 1 0; 1/2 sqrt(3/4)]

    # set number of lines and create empty stecil table
    numlines = cub.numedge*3 + 3 
    midedges = 0
    if cub.midedges
        numlines += 3
        midedges = 1
    end
    centroid=0
    if cub.centroid
        numlines +=2
        centroid = 1
    end

    lines = zeros(Int64,numlines, cub.numnodes)
    if cub.centroid
        lines[numlines-1,cub.numnodes] = cub.numnodes
        lines[numlines,cub.numnodes] = cub.numnodes
    end
    xy = SymCubatures.calcnodes(cub, vtx);

    # add vertices and S21 nodes (each node in these groups will have its own line)
    lines[diagind(lines)[1:numlines-2*centroid]] .= 1:numlines-2*centroid
    # add stencil for the vertices
    lines[1,2]=2; lines[2,3]=3; lines[3,4]=1;
    # add nodes in the numedge group to lines of the vertices
    j = 3+midedges*3+cub.numS21*3+1
    for m = 1:cub.numedge
        for i = 1:3
            for k=1:2
                lines[i, j] = j
                j+=1
            end
        end
    end
    # add nodes in the midedge groups to lines of the vertices
    if cub.centroid
        lines[1,4]=4
        lines[2,5]=5
        lines[3,6]=6
    end

    # add connection between S21 and numedge nodes
    # first get the node ordering on line 1
    idx_l1 = lines[1,findall(x -> x != 0,lines[1,:])[3:end]]
    idx_l1 = idx_l1[sortperm(xy[1,idx_l1])]
    # get node ordering on S21 line
    lst_S21 = []
    for i = (3+3*midedges+1):(3+3*midedges+3*cub.numS21)
        if abs(xy[1,i]-0.5)< 1e-6
            push!(lst_S21,i)
        end
    end
    idx_sort_S21 = sortperm(xy[2,lst_S21],rev=true)
    idx_S21 = lst_S21[idx_sort_S21]

    # get S21 nodes that connect with the edge nodes (top of centroid) and that do not (bottom)
    idx_top_S21 = idx_S21[1:cub.numedge]
    idx_bot_S21 = idx_S21[cub.numedge+1:end]

    # add S21 nodes that connect to the edge nodes to the stencil table
    k=1
    for i=1:length(idx_top_S21)
        for j=1:3
            lines[3+3*midedges+j+3*(i-1), 3+3*midedges+j+3*(i-1)] = idx_top_S21[i]+j-1
        end
        k+=1
    end

    # connect the S21 node stencil to the respective edge nodes
    m=1
    for i = (3+3*midedges+1):3:(3+3*midedges+length(idx_top_S21)*3)
        j=idx_l1[m]
        lines[i,j] = j
        lines[i,j+3] = j+3
        for k=1:2
            lines[i+k,j+3-k] = j+3-k
            lines[i+k,j+3+3-k] = j+3+3-k
        end
        m+=1
    end

    # add connection between S21 nodes 
    k=1
    for i = (3+3*midedges+1):3:(3+3*midedges+length(idx_top_S21)*3)
        idx_temp = collect(idx_top_S21[k]:idx_top_S21[k]+2)
        for j=0:2
            circshift!(idx_temp,-1)
            lines[i+j,idx_temp[1]] = idx_temp[1]
        end
        k+=1
    end
    
    # find and add the nodes in the S111 group and S21(below centroid) using polynomial fit
    for i=(4+3*midedges):numlines-2*centroid
        x = xy[1,lines[i,findall(a->a!=0, lines[i,:])]]
        y = xy[2,lines[i,findall(a->a!=0, lines[i,:])]]
        idx_sort = sortperm(x)
        x = x[idx_sort]
        y = y[idx_sort]
        # println("nodes = ", lines[i,findall(a->a!=0, lines[i,:])])
        # println("x = ", x)
        # println("y = ", y)
        interp = interpolate(x,y,FritschCarlsonMonotonicInterpolation())
        lines_temp=[]
        dist_temp=[]

        # add S111 nodes
        for j=(3+3*midedges+cub.numS21*3+6*cub.numedge+1):cub.numnodes
            if (xy[1,j]>=x[1] && xy[1,j]<=x[end])
                dist = norm(Array([xy[1,j];interp(xy[1,j])])-xy[:,j])
                # println("node: ", j, "  fit= ", interp(xy[1,j]), "   y = ", xy[1,j], "  dist = ", dist)
                push!(dist_temp,dist)
                push!(lines_temp,j)
            end
        end

        # add S21 nodes at the bottom of the centroid
        for j in idx_bot_S21
            for k=0:2
                if (xy[1,j+k]>=x[1] && xy[1,j+k]<=x[end])
                    dist = norm(Array([xy[1,j+k];interp(xy[1,j+k])])-xy[:,j+k])
                    push!(dist_temp,dist)
                    push!(lines_temp,j+k)
                end
            end
        end

        sort_dist = sortperm(dist_temp)
        lines_temp = lines_temp[sort_dist][1:p+2-length(x)]
        lines[i,lines_temp] = lines_temp
    end

    # add centroid connections
    if cub.centroid
        lines[numlines-1,4]=4
        lines[numlines-1,6]=6
        lines[numlines,4]=4
        lines[numlines,5]=5
    end
    if cub.midedges
        lines[4,3]=3
        lines[5,1]=1
        lines[6,2]=2
    end

    # find and add the midedge groups
    if cub.midedges
        for i=4:3+3*midedges
            x = xy[1,lines[i,findall(a->a!=0, lines[i,:])]]
            y = xy[2,lines[i,findall(a->a!=0, lines[i,:])]]
            idx_sort = sortperm(x)
            x = x[idx_sort]
            y = y[idx_sort]
            interp = interpolate(x,y,FritschCarlsonMonotonicInterpolation())
            lines_temp=[]
            dist_temp=[]

            # add all nodes on the midlines
            for j=7:cub.numnodes
                if (xy[1,j]>=x[1] && xy[1,j]<=x[end])
                    dist = norm(Array([xy[1,j];interp(xy[1,j])])-xy[:,j])
                    push!(dist_temp,dist)
                    push!(lines_temp,j)
                end
            end

            sort_dist = sortperm(dist_temp)
            lines_temp = lines_temp[sort_dist][1:p+2-length(x)]
            lines[i,lines_temp] = lines_temp
        end
    end

    # find and add the centroid groups
    if cub.centroid
        lines[numlines-1,idx_bot_S21] = idx_bot_S21
        lines[numlines-1,idx_bot_S21.+2] = idx_bot_S21.+2

        lines[numlines,idx_bot_S21] = idx_bot_S21
        lines[numlines,idx_bot_S21.+1] = idx_bot_S21.+1
    end
    
    # stencil_lines = lines
    stencil_lines = zeros(Int64, numlines, p+2)
    for i = 1:numlines
        stencil_lines[i,:] = lines[i,findall(x -> x != 0,lines[i,:])]
    end
    return stencil_lines
end

"""
### SummationByParts.node_sparse_stencil

Identify the stencil of sparse SBP operators for each node

**Inputs**

* `cub`: symmetric cubature rule
* `p`: maximum total degree for the Proriol polynomials

**Outputs**

* `node_stencil`: a cube containing the stencil in each direction for each node

"""
function node_sparse_stencil(cub::TriSymCub{T}, p::Int) where {T}

    # vtx = T[0 0; 1 0; 1/2 sqrt(3/4)]
    vtx = T[-1 -1; 1 -1; -1 1]
    xy = SymCubatures.calcnodes(cub, vtx)
    lines = sparse_stencil(cub, p)
    node_stencil = zeros(Int64,cub.numnodes,p+2,2)

    for i = 1:cub.numnodes
        idx = findall(x -> i in x, eachrow(lines))
        line1 = lines[idx[1],:]
        line2 = lines[idx[2],:]

        line1_x = line1[sortperm(xy[1,line1])]
        line2_x = line2[sortperm(xy[1,line2])]
        line1_y = line1[sortperm(xy[2,line1])]
        line2_y = line2[sortperm(xy[2,line2])]
        line1 = line1_x
        line2 = line2_y

        ylonger = (abs((xy[2,line1_y[end]])-(xy[2,line1_y[1]])) - abs((xy[1,line1_x[end]])-(xy[1,line1_x[1]])) )>1e-12
        xlonger = (abs((xy[1,line2_x[end]])-(xy[1,line2_x[1]])) - abs((xy[2,line2_y[end]])-(xy[2,line2_y[1]])) )>1e-12
        if ylonger || xlonger
            line1 = line2_x
            line2 = line1_y
        end
        if ylonger 
            line1 = line2_x
            line2 = line1_y
        end
        node_stencil[i,:,:] = Array([line1;line2])'
    end

    ddxi = zeros(Float64,2,cub.numnodes)
    ddeta = zeros(Float64,2,cub.numnodes)
    jac = zeros(Float64,1,cub.numnodes)
    oper_lgl = getLineSegSBPLobbato(degree=p+1)
    Q = oper_lgl.Q[:,:,1]
    H = diagm(oper_lgl.w)
    D = inv(H)*Q[:,:,1]

    x1 = vec(SymCubatures.calcnodes(oper_lgl.cub, oper_lgl.vtx))
    x1_sort = sort(x1)
    perm = []
    for i=1:length(x1)
        push!(perm, findall(a->x1[i] in a, x1_sort)[1])
    end
    for i=1:cub.numnodes
        for j=1:2
            node_stencil[i,:,j] = node_stencil[i,perm,j]
        end
    end

    D2 = zeros(Float64,cub.numnodes,cub.numnodes,2)
    for i = 1:2
        for j=1:cub.numnodes
            D2[j,node_stencil[j,:,i],i] = D[findall(a->j in a, node_stencil[j,:,i])[1],:]
        end
    end
    # for i = 1:2
    #     ddxi[i,:] = D2[:,:,1] * (xy'[:,i])
    #     ddeta[i,:] = D2[:,:,2] * (xy'[:,i])
    # end

    # x1 = SymCubatures.calcnodes(oper_lgl.cub, oper_lgl.vtx)
    # y1 = 0 .*x1
    # V, Vdx, Vdy = OrthoPoly.vandermonde(p, reshape(x1,(length(x1),1)), reshape(y1,(length(x1),1)))
    # Vd = [Vdx,Vdy]
    for i = 1:2
        ddxi[i,:] = D2[:,:,1] * (xy'[:,i])
        ddeta[i,:] = D2[:,:,2] * (xy'[:,i])
        # for j = 1:cub.numnodes
        #     c1 = pinv(V)*xy[i,node_stencil[j,:,1]]
        #     ddxi[i,j] = (Vd[i]*c1)[findall(a->j in a,node_stencil[j,:,1])][1]
        #     c2 = pinv(V)*xy[i,node_stencil[j,:,2]]
        #     ddeta[i,j] = (Vd[i]*c2)[findall(a->j in a,node_stencil[j,:,2])][1]
        # end
    end

    for i = 1:cub.numnodes
        jac[1,i] = ddxi[1,i]*ddeta[2,i] - ddeta[1,i]*ddxi[2,i]
        if jac[1,i] <= 0.0
            # xtemp = node_stencil[i,:,1]
            # node_stencil[i,:,1] = node_stencil[i,:,2]
            # node_stencil[i,:,2] = xtemp
            node_stencil[i,1:2,1] = node_stencil[i,:,1][2:-1:1]
            node_stencil[i,3:end,1] = node_stencil[i,:,1][end:-1:3]
        end
    end

    return node_stencil
end

# function node_sparse_stencil(cub::TriSymCub{T}, p::Int) where {T}

#     vtx = T[0 0; 1 0; 1/2 sqrt(3/4)]
#     xy = SymCubatures.calcnodes(cub, vtx)
#     lines = sparse_stencil(cub, p)
#     node_stencil = zeros(Int64,cub.numnodes,p+2,2)

#     for i = 1:cub.numnodes
#         idx = findall(x -> i in x, eachrow(lines))
#         line1 = lines[idx[1],:]
#         line2 = lines[idx[2],:]

#         line1_x = line1[sortperm(xy[1,line1])]
#         line2_x = line2[sortperm(xy[1,line2])]

#         x = line1_x
#         y = line2_x
#         # if xy[1,x[1]] > xy[1,y[1]]
#         #     x = line2_x
#         #     y = line1_x
#         # end
#         # if (xy[2,x[1]] < 1e-10 && xy[1,x[1]] > 1e-10)
#         #     xtemp = copy(x)
#         #     x = y
#         #     y = xtemp
#         # end

#         node_stencil[i,:,:] = Array([x;y])'
#     end

#     ddxi = zeros(Float64,2,cub.numnodes)
#     ddeta = zeros(Float64,2,cub.numnodes)
#     jac = zeros(Float64,1,cub.numnodes)
#     oper_lgl = getLineSegSBPLobbato(degree=p+1)
#     Q = oper_lgl.Q[:,:,1]
#     H = diagm(oper_lgl.w)
#     D = inv(H)*Q[:,:,1]

#     x1 = vec(SymCubatures.calcnodes(oper_lgl.cub, oper_lgl.vtx))
#     x1_sort = sort(x1)
#     perm = []
#     for i=1:length(x1)
#         push!(perm, findall(a->x1[i] in a, x1_sort)[1])
#     end
#     for i=1:cub.numnodes
#         for j=1:2
#             node_stencil[i,:,j] = node_stencil[i,perm,j]
#         end
#     end

#     Dxi = zeros(Float64,cub.numnodes,cub.numnodes,2)
#     for i = 1:2
#         for j=1:cub.numnodes
#             Dxi[j,node_stencil[j,:,i],i] = D[findall(a->j in a, node_stencil[j,:,i])[1],:]
#         end
#     end
#     for i = 1:2
#         ddxi[i,:] = Dxi[:,:,1] * (xy'[:,i])
#         ddeta[i,:] = Dxi[:,:,2] * (xy'[:,i])
#     end

#     for i = 1:cub.numnodes
#         jac[1,i] = ddxi[1,i]*ddeta[2,i] - ddeta[1,i]*ddxi[2,i]
#         if jac[1,i] <= 0.0
#             xtemp = node_stencil[i,:,1]
#             node_stencil[i,:,1] = node_stencil[i,:,2]
#             node_stencil[i,:,2] = xtemp
#         end
#     end

#     # for i=1:cub.numnodes
#     #     if i==9 || i==8 || i==6 || i==2
#     #         xtemp = node_stencil[i,:,1]
#     #         node_stencil[i,:,1] = node_stencil[i,:,2]
#     #         node_stencil[i,:,2] = xtemp
#     #         # node_stencil[i,1:2,2] = node_stencil[i,:,2][2:-1:1]
#     #         node_stencil[i,1:2,1] = node_stencil[i,:,1][2:-1:1]
#     #         # node_stencil[i,3:end,2] = node_stencil[i,:,2][end:-1:3]
#     #         node_stencil[i,3:end,1] = node_stencil[i,:,1][end:-1:3]
#     #     end
#     # end
   
#     return node_stencil
# end

"""
### SummationByParts.sparse_metric_terms

Computes the metric terms at each node for sparse SBP operators. 
The metric terms are approximated using 1D SBP operator.

**Inputs**

* `cub`: symmetric cubature rule
* `p`: maximum total degree for the Proriol polynomials

**Outputs**

* `ddxi` : contains a vector of [dx/dξ; dy/dξ]'
* `ddeta`: contains a vector of [dx/dη; dy/dη]'
* `dxid` : contains a vector of [dξ/dx; dξ/dy]'
* `detad`: contains a vector of [dη/dx; dη/dy]'
* `jac`  : contains the determinant of the metric Jacobian at each node
"""

function sparse_metric_terms(cub::TriSymCub{T}, vtx::Array{T,2}, p::Int) where {T}
    
    ddxi = zeros(Float64,2,cub.numnodes)
    dxid = zeros(Float64,2,cub.numnodes)
    ddeta = zeros(Float64,2,cub.numnodes)
    detad = zeros(Float64,2,cub.numnodes)
    jac = zeros(Float64,1,cub.numnodes)
    node_lines = node_sparse_stencil(cub, p)
    xy = SymCubatures.calcnodes(cub, vtx)

    oper_lgl = getLineSegSBPLobbato(degree=p+1)
    Q = oper_lgl.Q[:,:,1]
    H = diagm(oper_lgl.w)
    D = inv(H)*Q[:,:,1]
    x1 = SymCubatures.calcnodes(oper_lgl.cub, oper_lgl.vtx)
    y1 = 0 .*x1

    # for i = 1:cub.numnodes
    #     for j = 1:2
    #         ddxi[j,i] = ((D * reshape(xy[j,node_lines[i,:,1]],(p+2,1)))[findall(a->i in a, node_lines[i,:,1])])[1]
    #         ddeta[j,i] = ((D * reshape(xy[j,node_lines[i,:,2]],(p+2,1)))[findall(a->i in a, node_lines[i,:,2])])[1]
    #     end
    # end

    D2 = zeros(Float64,cub.numnodes,cub.numnodes,2)
    for i = 1:2
        for j=1:cub.numnodes
            D2[j,node_lines[j,:,i],i] = D[findall(a->j in a, node_lines[j,:,i])[1],:]
            # V, dV = OrthoPoly.vandermonde_monomial(p+1, xy[i, node_lines[j,:,i]])
            # dd = dV*inv(V)
            # D2[j,node_lines[j,:,i],i] = dd[findall(a->j in a, node_lines[j,:,i])[1],:]
        end
    end
    # V, Vdx, Vdy = OrthoPoly.vandermonde(p, reshape(x1,(length(x1),1)), reshape(y1,(length(x1),1)))
    # Vd = [Vdx,Vdy]
    for i = 1:2
        ddxi[i,:] = D2[:,:,1] * ((xy')[:,i])
        ddeta[i,:] = D2[:,:,2] * ((xy')[:,i])
        # for j = 1:cub.numnodes
        #     c1 = pinv(V)*xy[i,node_lines[j,:,1]]
        #     ddxi[i,j] = (Vd[i]*c1)[findall(a->j in a,node_lines[j,:,1])][1]
        #     c2 = pinv(V)*xy[i,node_lines[j,:,2]]
        #     ddeta[i,j] = (Vd[i]*c2)[findall(a->j in a,node_lines[j,:,2])][1]
        # end
    end

    for i = 1:cub.numnodes
        jac[1,i] = ddxi[1,i]*ddeta[2,i] - ddeta[1,i]*ddxi[2,i]

        if jac[1,i] <= 0.0
            @warn "There is a negative metric Jacobian term."
        end
    end

    for i = 1:cub.numnodes
        dxid[1,i] = 1.0/jac[1,i] * ddeta[2,i] 
        dxid[2,i] = -1.0/jac[1,i] * ddeta[1,i] 
        detad[1,i] = -1.0/jac[1,i] * ddxi[2,i]
        detad[2,i] = 1.0/jac[1,i] * ddxi[1,i]
    end
    # dxdxi = ddxi[1,:]
    # dydxi = ddxi[2,:]
    # dxdeta = ddeta[1,:]
    # dydeta = ddeta[2,:]
    # dxidx = dxid[1,:]
    # dxidy = dxid[2,:]
    # detadx = detad[1,:]
    # detady = detad[2,:]
    return ddxi, ddeta, dxid, detad, jac
end


"""
### SummationByParts.sparse_sbp_operator

Constructs a sparse SBP operator

**Inputs**

* `cub`: symmetric cubature rule
* `p`: maximum total degree for the Proriol polynomials

**Outputs**

* `D`: the strong derivative matrix
* `Q`: the weak derivative matrix
* `H`: the norm matrix
* `E`: the boundary integration oprator
"""
function sparse_sbp_operator(cub::TriSymCub{T}, vtx::Array{T,2}, p::Int) where {T}

    ddxi, ddeta, dxid, detad, jac = sparse_metric_terms(cub, vtx, p)
    node_lines = node_sparse_stencil(cub, p)
    oper_lgl = getLineSegSBPLobbato(degree=p+1)
    Q1 = oper_lgl.Q[:,:,1]
    H1 = diagm(oper_lgl.w)
    D1 = inv(H1)*Q1[:,:,1]

    H2 = zeros(T,cub.numnodes,cub.numnodes,2)
    for i=1:cub.numnodes
        ix = findall(a->i in a, node_lines[i,:,1])[1]
        iy = findall(a->i in a, node_lines[i,:,2])[1]
        H2[i,i,1] = H1[ix,ix] 
        H2[i,i,2] = H1[iy,iy]
    end

    Q2 = zeros(T,cub.numnodes,cub.numnodes,2)
    for i = 1:2
        for j=1:cub.numnodes
            Q2[j,node_lines[j,:,i],i] = Q1[findall(a->j in a, node_lines[j,:,i])[1],:]
        end
    end

    D2 = zeros(T,cub.numnodes,cub.numnodes,2)
    for i = 1:2
        for j=1:cub.numnodes
            D2[j,node_lines[j,:,i],i] = D1[findall(a->j in a, node_lines[j,:,i])[1],:]
        end
    end

    # dxid, detad = accurate_metric_terms_frobenius(p, D2)

    w, Qdense, E = SummationByParts.buildoperators(cub, vtx, p,vertices=true)
    H = diagm(w) #.*diagm(jac[1,:])
    Q = zeros(T,cub.numnodes,cub.numnodes,2)
    D = zeros(T,cub.numnodes,cub.numnodes,2)
    S = zeros(T,cub.numnodes,cub.numnodes,2)

    # S[:,:,1] = 0.5.*(diagm(ddeta[2,:])*Qxi[:,:,1] - Qxi[:,:,1]' * diagm(ddeta[2,:])) + 
    #            0.5.*(-diagm(ddxi[2,:])*Qxi[:,:,2] + Qxi[:,:,2]' * diagm(ddxi[2,:]))
    # S[:,:,2] = 0.5.*(-diagm(ddeta[1,:])*Qxi[:,:,1] + Qxi[:,:,1]' * diagm(ddeta[1,:])) + 
    #            0.5.*(diagm(ddxi[1,:])*Qxi[:,:,2] - Qxi[:,:,2]' * diagm(ddxi[1,:]))
    
    S[:,:,1] = 0.5.*H*(diagm(dxid[1,:])*D2[:,:,1] + diagm(detad[1,:])*D2[:,:,2]) -
               0.5.*(D2[:,:,1]'*diagm(dxid[1,:]) + D2[:,:,2]'*diagm(detad[1,:]))*H
    S[:,:,2] = 0.5.*H*(diagm(dxid[2,:])*D2[:,:,1] + diagm(detad[2,:])*D2[:,:,2]) -
               0.5.*(D2[:,:,1]'*diagm(dxid[2,:]) + D2[:,:,2]'*diagm(detad[2,:]))*H

    for i = 1:2
        Q[:,:,i] = S[:,:,i] + 0.5 .* E[:,:,i]
        D[:,:,i] = inv(H)*Q[:,:,i]
    end

    D[:,:,1] = diagm(dxid[1,:])*D2[:,:,1] + diagm(detad[1,:])*D2[:,:,2]
    D[:,:,2] = diagm(dxid[2,:])*D2[:,:,1] + diagm(detad[2,:])*D2[:,:,2]
    Q[:,:,1] = H*D[:,:,1]
    Q[:,:,2] = H*D[:,:,2]

    return D, Q, H, E, S, D2, Q2, H2, ddxi, ddeta, dxid, detad, jac, D1
    # return D, Q, H, E, S
end

function test_accuracy(cub::TriSymCub{T},p::Int) where{T}
    # T = Float64
    vtx = T[-1 -1; 1 -1; -1 1];
    # q=2*p-1
    # cub,vtx = SummationByParts.getTriCubatureDiagE(q, T, vertices=true)
    xy = SymCubatures.calcnodes(cub, vtx)
    x = (xy')[:,1]
    y = (xy')[:,2]
    D, Q, H, E, S = SummationByParts.sparse_sbp_operator(cub, vtx, p)
    Dx = D[:,:,1]
    Dy = D[:,:,2]

    # oper_tri = SummationByParts.getTriSBPDiagE(degree=p,vertices=true,quad_degree=q)
    # Qx = oper_tri.Q[:,:,1];
    # Qy = oper_tri.Q[:,:,2];
    # H = diagm(oper_tri.w);
    # Ex = oper_tri.E[:,:,1];
    # Ey = oper_tri.E[:,:,2];
    # Dx = inv(H)*Qx;
    # Dy = inv(H)*Qy;

    errs_val = []
    errs_x = []
    errs_y = []
    # h_vec = [1,1/4,1/16,1/32,1/64,1/128]
    h_vec = [4.0^-i for i in 0:5]
    ω = 2*π
    for h in h_vec
        f = sin.(ω*h*x).*sin.(ω*h*y)
        fx = ω*h*cos.(ω*h*x).*sin.(ω*h*y)
        fy = ω*h*sin.(ω*h*x).*cos.(ω*h*y)

        fx_num = Dx*f
        fy_num = Dy*f

        push!(errs_val, (fx .- fx_num))
        err_x = sqrt((fx .- fx_num)' * H *(fx .- fx_num))
        err_y = sqrt((fy .- fy_num)' * H *(fy .- fy_num))
        # err_x = norm(fx .- fx_num)
        # err_y = norm(fy .- fy_num)
        push!(errs_x, err_x)
        push!(errs_y, err_y)
    end

    rate_x = []
    rate_y = []
    for i = 1:length(h_vec)-1
        push!(rate_x, log10(errs_x[i+1]/errs_x[i])/log10(h_vec[i+1]/h_vec[i]))
        push!(rate_y, log10(errs_y[i+1]/errs_y[i])/log10(h_vec[i+1]/h_vec[i]))
    end

    return errs_x, errs_y, rate_x, rate_y, h_vec, errs_val
end

function test_accuracy(oper::TriSBP{T},p::Int) where{T}
    # T = Float64
    vtx = T[-1 -1; 1 -1; -1 1];
    # q=2*p-1
    # cub,vtx = SummationByParts.getTriCubatureDiagE(q, T, vertices=true)
    xy = AsymCubatures.calcnodes(oper.cub)
    x = (xy')[:,1]
    y = (xy')[:,2]
    H = diagm(oper.w)
    Dx = inv(H)*oper.Q[:,:,1]
    Dy = inv(H)*oper.Q[:,:,2]


    # oper_tri = SummationByParts.getTriSBPDiagE(degree=p,vertices=true,quad_degree=q)
    # Qx = oper_tri.Q[:,:,1];
    # Qy = oper_tri.Q[:,:,2];
    # H = diagm(oper_tri.w);
    # Ex = oper_tri.E[:,:,1];
    # Ey = oper_tri.E[:,:,2];
    # Dx = inv(H)*Qx;
    # Dy = inv(H)*Qy;

    errs_val = []
    errs_x = []
    errs_y = []
    # h_vec = [1,1/4,1/16,1/32,1/64,1/128]
    h_vec = [4.0^-i for i in 0:5]
    ω = 2*π
    for h in h_vec
        f = sin.(ω*h*x).*sin.(ω*h*y)
        fx = ω*h*cos.(ω*h*x).*sin.(ω*h*y)
        fy = ω*h*sin.(ω*h*x).*cos.(ω*h*y)

        fx_num = Dx*f
        fy_num = Dy*f

        push!(errs_val, (fx .- fx_num))
        err_x = sqrt((fx .- fx_num)' * H *(fx .- fx_num))
        err_y = sqrt((fy .- fy_num)' * H *(fy .- fy_num))
        # err_x = norm(fx .- fx_num)
        # err_y = norm(fy .- fy_num)
        push!(errs_x, err_x)
        push!(errs_y, err_y)
    end

    rate_x = []
    rate_y = []
    for i = 1:length(h_vec)-1
        push!(rate_x, log10(errs_x[i+1]/errs_x[i])/log10(h_vec[i+1]/h_vec[i]))
        push!(rate_y, log10(errs_y[i+1]/errs_y[i])/log10(h_vec[i+1]/h_vec[i]))
    end

    return errs_x, errs_y, rate_x, rate_y, h_vec, errs_val
end

function test_accuracy(Q::Array{T,2},w::Array{T,1},x::Array{T,1},y::Array{T,1};refine=4) where{T}

    H = diagm(w)
    # Dx = inv(H)*Q[:,:,1]
    Dx = inv(H)*Q
    
    errs_val = []
    errs_x = []

    h_vec = [4.0^-i for i in 0:refine]
    ω = 2*π
    for h in h_vec
        f = sin.(ω*h*x).*sin.(ω*h*y)
        fx = ω*h*cos.(ω*h*x).*sin.(ω*h*y)

        fx_num = Dx*f

        push!(errs_val, (fx .- fx_num))
        err_x = sqrt((fx .- fx_num)' * H *(fx .- fx_num))
        # err_x = norm(fx .- fx_num)
        push!(errs_x, err_x)
    end

    rate_x = []
    for i = 1:length(h_vec)-1
        push!(rate_x, log10(errs_x[i+1]/errs_x[i])/log10(h_vec[i+1]/h_vec[i]))
    end

    return errs_x, rate_x, h_vec, errs_val
end

function accurate_metric_terms_newton(p::Int)
    tol = 1e-13
    err1 = 1.0
    T = Float64

    q=2*p
    cub,vtx = SummationByParts.getTriCubatureDiagE(q, T, vertices=true)
    # cub,vtx = Cubature.getTriCubatureDiagE(q, T, vertices=true)
    xy = SymCubatures.calcnodes(cub, vtx)
    x = (xy')[:,1]
    y = (xy')[:,2]
    V, Vdx, Vdy = OrthoPoly.vandermonde(p+1, x, y)
    D, Q, H, E, S, D2, Q2, H2 = SummationByParts.sparse_sbp_operator(cub, vtx, p)
    xi_x = diagm(dxid[1,:])
    eta_x = diagm(detad[1,:])
    D_xi = D2[:,:,1]
    D_eta = D2[:,:,2]

    # while ((err1 > tol))
    #     xi_x = diag((Vdx*pinv(V) - diagm(eta_x)*D_eta)*pinv(D_xi))
    #     eta_x = diag((Vdx*pinv(V) - diagm(xi_x)*D_xi)*pinv(D_eta))
        
    #     err1 = norm((diagm(xi_x)*D_xi + diagm(eta_x)*D_eta)*V - Vdx)
    #     println("err1: ", err1)
    # end
    N1 = ones(cub.numnodes,1)
    m = Array([xi_x; eta_x])*N1
    d = Array([D_xi; D_eta])
    J = V'*d'
    Jinv = pinv(J)
    for i = 1:10
        F = V'*d'*m - Vdx'*N1
        m = m - Jinv*F
        println(norm(F))
    end

    xi_x_opt = m[1:cub.numnodes,1]
    eta_x_opt = m[cub.numnodes+1:2*cub.numnodes,1]
    err1 = N1'*((diagm(xi_x_opt)*D_xi + diagm(eta_x_opt)*D_eta)*V - Vdx)
    err2 = ((diagm(xi_x_opt)*D_xi + diagm(eta_x_opt)*D_eta)*V - Vdx)
    writedlm(stdout, round.(err1, sigdigits=6))
    writedlm(stdout, round.(err2, sigdigits=6))
    return m, xi_x, eta_x
end

function accurate_metric_terms_frobenius(p::Int, D2::Array)
    T = Float64

    q=2*p
    cub,vtx = Cubature.getTriCubatureDiagE(q, T, vertices=true)
    xy = SymCubatures.calcnodes(cub, vtx)
    x = (xy')[:,1]
    y = (xy')[:,2]
    V, Vdx, Vdy = OrthoPoly.vandermonde(p, x, y)
    dxid = zeros(Float64,2,cub.numnodes)
    detad = zeros(Float64,2,cub.numnodes)

    N1 = ones(cub.numnodes,1)
    Nv = ones(size(V)[2],1)

    # m = Array([xi_x; eta_x])*N1
    a = D2[:,:,1]*V
    b = D2[:,:,2]*V
    
    Ax = Array([diagm(vec((a.*a)*Nv)) diagm(vec((a.*b)*Nv)); diagm(vec((b.*a)*Nv)) diagm(vec((b.*b)*Nv))])
    Ay = Array([diagm(vec((a.*a)*Nv)) diagm(vec((a.*b)*Nv)); diagm(vec((b.*a)*Nv)) diagm(vec((b.*b)*Nv))])
    Cx = Array([(a.*Vdx)*Nv; (b.*Vdx)*Nv])
    Cy = Array([(a.*Vdy)*Nv; (b.*Vdy)*Nv])

    # Ax = Array([diagm(vec((a.*a)*Nv)) diagm(vec((a.*b)*Nv)); diagm(vec((b.*a)*Nv)) diagm(vec((b.*b)*Nv));
    #            diagm(vec((Vdx.*a)*Nv)) diagm(vec((Vdx.*b)*Nv))])
    # Ay = Array([diagm(vec((a.*a)*Nv)) diagm(vec((a.*b)*Nv)); diagm(vec((b.*a)*Nv)) diagm(vec((b.*b)*Nv));
    #            diagm(vec((Vdy.*a)*Nv)) diagm(vec((Vdy.*b)*Nv))])
    # Cx = Array([(a.*Vdx)*Nv; (b.*Vdx)*Nv; (Vdx.*Vdx)*Nv])
    # Cy = Array([(a.*Vdy)*Nv; (b.*Vdy)*Nv; (Vdx.*Vdx)*Nv])

    m_opt_x = pinv(Ax)*Cx
    m_opt_y = pinv(Ay)*Cy
    dxid[1,:] = m_opt_x[1:cub.numnodes,1]
    detad[1,:] = m_opt_x[cub.numnodes+1:2*cub.numnodes,1]
    dxid[2,:] = m_opt_y[1:cub.numnodes,1]
    detad[2,:] = m_opt_y[cub.numnodes+1:2*cub.numnodes,1]

    # err1 = N1'*((diagm(dxid[1,:])*D2[:,:,1] + diagm(detad[1,:])*D2[:,:,2])*V - Vdx)*Nv
    # err2 = N1'*((diagm(dxid[2,:])*D2[:,:,1] + diagm(detad[2,:])*D2[:,:,2])*V - Vdy)*Nv
    # writedlm(stdout, round.(err1, sigdigits=6))
    # writedlm(stdout, round.(err2, sigdigits=6))

    return dxid, detad
end

function test_accuracy_D2(cub::TriSymCub{T},p::Int) where{T}
    # T = Float64
    # q=2*p-1
    # cub,vtx = Cubature.getTriCubatureDiagE(q, T, vertices=true)    
    vtx = T[-1 -1; 1 -1; -1 1]
    D, Q, H, E, S, D2, Q2, H2 = SummationByParts.sparse_sbp_operator(cub, vtx, p)
    node_lines = SummationByParts.node_sparse_stencil(cub, p)

    oper_lgl = getLineSegSBPLobbato(degree=p+1)
    x1 = vec(SymCubatures.calcnodes(oper_lgl.cub, oper_lgl.vtx))
    V, Vdx = OrthoPoly.vandermonde_monomial(p+1, x1)

    V1 = zeros(T, cub.numnodes, size(V)[2], cub.numnodes)
    Vdx1 = zeros(T, cub.numnodes, size(V)[2], cub.numnodes)
    V2 = zeros(T, cub.numnodes, size(V)[2], cub.numnodes)
    Vdx2 = zeros(T, cub.numnodes, size(V)[2], cub.numnodes)
    for i = 1:cub.numnodes
        for j = 1:size(V)[2]
            V1[node_lines[i,:,1], j, i] = V[:,j]
            Vdx1[node_lines[i,:,1], j, i] = Vdx[:,j] 
            V2[node_lines[i,:,2], j, i] = V[:,j]
            Vdx2[node_lines[i,:,2], j, i] = Vdx[:,j] 
        end
    end
    
    err_x = 0.0
    err_y = 0.0
    for i = 1:cub.numnodes
        err_x += norm((D2[:,:,1]*V1[:,:,i] - Vdx1[:,:,i])[i,:])
        err_y += norm((D2[:,:,2]*V2[:,:,i] - Vdx2[:,:,i])[i,:])
        # println((D2[:,:,1]*V1[:,:,i] - Vdx1[:,:,i])[i,:])
        # println((D2[:,:,2]*V2[:,:,i] - Vdx2[:,:,i])[i,:])
    end
    return err_x, err_y
end 

#------------------------------------------------------

function compute_mutual_coherence(A)
    M = 0.0
    for i = 1:size(A,2)-1
        for j = i+1:size(A,2)
            m = abs(A[:,i]'*A[:,j])/(norm(A[:,i])*norm(A[:,j]))
            if m > M
                M = m 
            end
        end
    end
    return M
end

"""
### SummationByParts.closest_node_stencil

Identify the closest k nodes to each node in the element

**Inputs**

* `cub`: symmetric cubature rule
* `k`: maximum number of nodes

**Outputs**

* `node_stencil`: a cube containing the stencil in each direction for each node

"""
function closest_node_stencil(cub::Union{TriSymCub{T},TriAsymCub{T}}, k::Int) where {T}

    vtx = T[-1 -1; 1 -1; -1 1]
    if typeof(cub)==TriSymCub{T}
        xy = SymCubatures.calcnodes(cub, vtx)'
    elseif typeof(cub)==TriAsymCub{T}
        xy = AsymCubatures.calcnodes(cub)'
    end
    closest_stencil = zeros(Int64,cub.numnodes,k)
    function dist(point1, point2)
        return norm(point1 .- point2)
    end
    for j = 1:cub.numnodes
        d = [dist(xy[i, :], xy[j,:]) for i in 1:cub.numnodes]
        closest_stencil[j,:] = sortperm(d)[1:k]
    end
    return closest_stencil
end

function closest_node_stencil(cub::TetSymCub{T}, k::Int) where {T}

    vtx = T[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
    xyz = SymCubatures.calcnodes(cub, vtx)'
    closest_stencil = zeros(Int64,cub.numnodes,k)
    function dist(point1, point2)
        return norm(point1 .- point2)
    end
    for j = 1:cub.numnodes
        d = [dist(xyz[i, :], xyz[j,:]) for i in 1:cub.numnodes]
        closest_stencil[j,:] = sortperm(d)[1:k]
    end
    return closest_stencil
end

"""
### SummationByParts.init_tri_nodes

Initialize nodal points for sparse sbp operators

**Inputs**

* `q`: maximum total degree for the quadrature rule

**In/Outs**

* `cub`: symmetric cubature rule

"""
function init_tri_nodes(p::Int,q::Int; ne::Int=-1, T=Float64)

    # use equilateral triangle
    node_tol = 1e-14
    vtx = T[-1 -1/sqrt(3); 1 -1/sqrt(3); 0 2/sqrt(3)] 
    if ne==-1
        ne=p+2
    end
    qf = 2*ne-3
    cub_lgl, _ = SummationByParts.Cubature.quadrature(qf, internal=false)
    xedge = cub_lgl.params
    wlgl = cub_lgl.weights

    # vertices = true
    # midedges = convert(Bool,mod(p,2))
    # numedge = convert(Int,(p-mod(p,2))/2)
    # numS21 = convert(Int,(1+mod(p,2))*numedge)
    # numS111 = convert(Int,1/2*(numedge^2 - numedge))
    # centroid = convert(Bool,mod(p,2))

    vertices = true
    midedges = convert(Bool,mod(ne-2,2))
    numedge = convert(Int,(ne-2-mod(ne-2,2))/2)
    numS21 = convert(Int,(1+mod(ne-2,2))*numedge)
    numS111 = convert(Int,1/2*(numedge^2 - numedge))
    centroid = convert(Bool,mod(ne-2,2))

    cub = SymCubatures.TriSymCub{T}(vertices=vertices,
                                    midedges=midedges,
                                    numS21=numS21,
                                    numedge=numedge,
                                    numS111=numS111,
                                    centroid = centroid)
    
    # set the edge parameters 
    init_tri_cub(cub, xedge)

    # get interior nodes
    interior_nodes, weights = get_interior_points_tri(p, ne=ne)

    # get barycentric coordinates of the interior nodes
    params=[]
    weight_vert=[]
    weight_mid=[]
    weight_cent=[]
    if interior_nodes!=[]
        param_S21, param_S111, weight_S21, weight_S111= get_bary_coord_interior(interior_nodes, weights, vtx)
        params = convert.(T,collect(Iterators.flatten([2.0 .*param_S21,cub_lgl.params,2.0 .*param_S111])))
        J = 1/2*sqrt(3)/3
        if vertices
            weight_vert=wlgl[1]^2
        end
        if midedges
            weight_mid=wlgl[end]*wlgl[end-convert(Int,centroid)]
        end
        if centroid
            weight_cent=J*wlgl[end]^2
        end
        weight_edge= wlgl[2:end-convert(Int,centroid)].*wlgl[1]
        weights = convert.(T,collect(Iterators.flatten([J.*weight_vert,J.*weight_mid, J.*weight_S21, 
                                                        J.*weight_S111, J.*weight_edge, J.*weight_cent])))
    end
    if params != []
        SymCubatures.setparams!(cub, params)
        SymCubatures.setweights!(cub, weights)
    end
    # xy = SymCubatures.calcnodes(cub,vtx)
    # SummationByParts.plot_tri_nodes(x=xy, vtx=vtx, q=q, n=cub.numnodes,write_title=true,label_nodes=false)

    xinit = convert.(T,collect(Iterators.flatten([cub.params,cub.weights])))
    mask = SymCubatures.getInternalParamMask(cub)
    append!(mask, (cub.numparams+1):(cub.numparams+cub.numweights))  
    Cubature.solvecubature!(cub, q, mask, tol=1e-14, hist=true, verbose=true, xinit=xinit, delta1=1e-4, delta2=1e-4)
    println("\n", cub.params,"\n")
    println(cub.weights,"\n")
    return cub
end

function init_tri_nodes_omega(p::Int, q::Int; ne::Int=-1, T=Float64)

    # use equilateral triangle
    node_tol = 1e-14
    vtx = T[-1 -1/sqrt(3); 1 -1/sqrt(3); 0 2/sqrt(3)] 
    if ne==-1
        ne=p+1
    end
    qf = 2*ne-1
    cub_lg, _ = SummationByParts.Cubature.quadrature(qf, internal=true)
    xedge = cub_lg.params
    wlg = cub_lg.weights
    
    nn = convert(Int,(ne-0-mod(ne-0,2))/2)

    vertices = false
    midedges = false #convert(Bool,mod(ne-2,2))
    numedge = 0 #convert(Int,(ne-2-mod(ne-2,2))/2)
    numS21 = convert(Int,(1+mod(ne-0,2))*nn)
    numS111 = convert(Int,1/2*(nn^2 - nn))
    centroid = convert(Bool,mod(ne-0,2))

    cub = SymCubatures.TriSymCub{T}(vertices=vertices,
                                    midedges=midedges,
                                    numS21=numS21,
                                    numedge=numedge,
                                    numS111=numS111,
                                    centroid = centroid)
    
    # set the edge parameters 
    init_tri_cub(cub, xedge)

    # get interior nodes
    interior_nodes, weights = get_interior_points_tri_omega(p)

    # get barycentric coordinates of the interior nodes
    params=[]
    weight_vert=[]
    weight_mid=[]
    weight_cent=[]
    if interior_nodes!=[]
        param_S21, param_S111, weight_S21, weight_S111= get_bary_coord_interior(interior_nodes, weights, vtx)
        params = convert.(T,collect(Iterators.flatten([2.0 .*param_S21, 2.0 .*param_S111])))
        J = 1 #1/2*sqrt(3)/2
        if centroid
            weight_cent=J*wlg[end]^2
        end
        weights = convert.(T,collect(Iterators.flatten([J.*weight_S21, 
                                                        J.*weight_S111, 
                                                        J.*weight_cent])))
    end
    if params != []
        SymCubatures.setparams!(cub, params)
        SymCubatures.setweights!(cub, weights)
    end
    # xy = SymCubatures.calcnodes(cub,vtx)
    # SummationByParts.plot_tri_nodes(x=xy, vtx=vtx, q=q, n=cub.numnodes,write_title=true,label_nodes=false)
    # println(params)
    # println(weights)
    cub, _ = SummationByParts.Cubature.getTriCubatureOmegaLG(q-1)

    xinit = convert.(T,collect(Iterators.flatten([cub.params,cub.weights])))
    mask = Int[]
    append!(mask, 1:cub.numparams+cub.numweights) 
    # Cubature.solvecubature!(cub, q, mask, tol=5e-14, hist=true, verbose=true, xinit=xinit, delta1=1e-4, delta2=1e-4)
    res = Cubature.solvecubaturelma!(cub, q, mask, tol=5e-14, maxiter=1000, hist=true, verbose=true, xinit=xinit)
    println("\n", cub.params,"\n")
    println(cub.weights,"\n")
    return cub, res
end

function init_tet_nodes(p::Int,q::Int; ne::Int=-1, T=Float64)

    # use equilateral triangle
    node_tol = 1e-14
    vtx = T[1 -sqrt(3)/3 -sqrt(6)/6;
            0 2*sqrt(3)/3 -sqrt(6)/6;
            -1 -sqrt(3)/3 -sqrt(6)/6;
            0 0 sqrt(6)/2]
    if ne==-1
        ne=p+2
    end
    qf = 2*ne-3
    cub_lgl, _ = SummationByParts.Cubature.quadrature(qf, internal=false)
    cub_tri,_ = SummationByParts.Cubature.getTriCubatureDiagE(q, T, vertices=true)
    xedge = cub_tri.params
    wlgl = cub_lgl.weights

    # get interior nodes
    interior_nodes, weights = get_interior_points_tet(p)

    # get barycentric coordinates of the interior nodes
    params=[]
    weight_vert=[]
    weight_mid=[]
    weight_cent=[]
    if interior_nodes!=[]
        param_S31, param_S22, param_S211, param_S1111, param_centroid, weight_S31, weight_S22, weight_S211, weight_S1111, weight_centroid= get_bary_coord_interior_tet(interior_nodes, vec(weights), vtx)        
        
        weight_vert=wlgl[1]^3      
        if convert(Bool,mod(p,2))
            weight_mid=wlgl[convert(Int, floor(p/2)+1)]^2*wlgl[end]
            # weight_cent=wlgl[convert(Int, floor(p/2)+1)]^3
        end
        # weight_edge= wlgl[2:end-convert(Int,centroid)].*wlgl[1]^2
        # weights = convert.(T,collect(Iterators.flatten([J.*weight_vert, weight_S31, J.*weight_mid, weight_S22, 
        #                                                 J.*weight_S211, 
        #                                                 J.*weight_S1111, J.*weight_edge, J.*weight_cent])))
    end

    vertices = true
    midedges = convert(Bool,mod(p,2))
    numedge = convert(Int,(p-mod(p,2))/2)
    numfaceS21 = convert(Int,(1+mod(p,2))*numedge)
    numfaceS111 = convert(Int,1/2*(numedge^2 - numedge))
    facecentroid = convert(Bool,mod(p,2))
    numS31 = length(param_S31)
    numS22 = length(param_S22)
    numS211 = convert(Int, length(param_S211)/2)
    numS1111 = convert(Int, length(param_S1111)/3)
    centroid = convert(Bool,mod(p,2))

    cub = SymCubatures.TetSymCub{T}(vertices=vertices,
                                    midedges=midedges,
                                    numfaceS21=numfaceS21,
                                    numedge=numedge,
                                    numfaceS111=numfaceS111,
                                    facecentroid = facecentroid,
                                    numS31=numS31,
                                    numS22=numS22,
                                    numS211=numS211,
                                    numS1111=numS1111,
                                    centroid=centroid)
    
    # set the edge parameters 
    init_tet_cub(cub, xedge)
    faceS21_params = cub_tri.params[1 : cub_tri.numS21]
    faceS21_w = wlgl[1].*cub_tri.weights[convert(Int,cub.vertices)+convert(Int,cub.midedges)+1 : convert(Int,cub.vertices)+convert(Int,cub.midedges)+cub_tri.numS21]
    edge_params = cub_tri.params[cub_tri.numS21+1 : cub_tri.numS21+cub_tri.numedge]
    edge_w = wlgl[1].*cub_tri.weights[convert(Int,cub.vertices)+convert(Int,cub.midedges)+cub.numfaceS21+1 : convert(Int,cub.vertices)+convert(Int,cub.midedges)+cub.numfaceS21+cub.numedge]
    faceS111_params = cub_tri.params[cub_tri.numS21+cub_tri.numedge+1 : cub_tri.numS21+cub_tri.numedge+2*cub_tri.numS111]
    faceS111_w =  wlgl[1].*cub_tri.weights[convert(Int,cub.vertices)+convert(Int,cub.midedges)+cub.numfaceS21+cub.numedge+1 : convert(Int,cub.vertices)+convert(Int,cub.midedges)+cub.numfaceS21+cub.numedge+cub.numfaceS111]
    facecentroid_w = [] 
    if mod(p,2)==1
        facecentroid_w=wlgl[1].*cub_tri.weights[end]
    end
    params = convert.(T,collect(Iterators.flatten([3.0 .*param_S31,2.0 .*param_S22, faceS21_params, 
                                                   edge_params, 2.0*param_S211, faceS111_params, 2.0 .*param_S1111])))
    
    # p=7 works with J=1/800                                               
    J = 1/4000
    Jf = 1/10                    
    weights = J.*convert.(T,collect(Iterators.flatten([weight_vert, weight_S31, weight_mid, weight_S22, 
                                                    Jf*faceS21_w, Jf*edge_w, weight_S211, Jf*faceS111_w, 
                                                    Jf*facecentroid_w, weight_S1111, weight_centroid])))
                           
    println(params,"\n")
    println(weights)                                                
    if params != []
        SymCubatures.setparams!(cub, params)
        SymCubatures.setweights!(cub, weights)
    end
    # xy = SymCubatures.calcnodes(cub,vtx)
    # SummationByParts.plotly_tet_nodes(q=q, x=xy, vtx=vtx)
    # SummationByParts.plotly_tet_nodes(q=q, x=interior_nodes, vtx=vtx)

    # n=mod(p,2)
    # for i=0:convert(Int,floor(p/2))
    #     ne = convert(Int, (p-2*i-mod(p,2))/2)
    #     n += 4*(3*ne^2 + 3*mod(p,2)*ne + mod(p,2)) + 6*(2*ne + mod(p,2)) + 4
    # end
    # println("\n","n = ", n)
    # println("numnodes = ", cub.numnodes)

    xinit = convert.(T,collect(Iterators.flatten([cub.params,cub.weights])))
    mask = SymCubatures.getInternalParamMask(cub)
    append!(mask, (cub.numparams+1):(cub.numparams+cub.numweights))  
    Cubature.solvecubature!(cub, q, mask, tol=5e-14, hist=true, verbose=true, xinit=xinit, delta1=1e-4, delta2=1e-4)
    println("\n", cub.params,"\n")
    println(cub.weights,"\n")
    return cub
end

# function init_tet_nodes_omega(p::Int,q::Int; ne::Int=-1, T=Float64)

#     # use equilateral triangle
#     node_tol = 1e-14
#     vtx = T[1 -sqrt(3)/3 -sqrt(6)/6;
#             0 2*sqrt(3)/3 -sqrt(6)/6;
#             -1 -sqrt(3)/3 -sqrt(6)/6;
#             0 0 sqrt(6)/2]
#     if ne==-1
#         ne=p+2
#     end
#     qf = 2*ne-3
#     cub_lgl, _ = SummationByParts.Cubature.quadrature(qf, internal=false)
#     cub_tri,_ = SummationByParts.Cubature.getTriCubatureDiagE(q, T, vertices=true)
#     xedge = cub_tri.params
#     wlgl = cub_lgl.weights

#     # get interior nodes
#     interior_nodes, weights = get_interior_points_tet(p, ne=ne)

#     # get barycentric coordinates of the interior nodes
#     params=[]
#     weight_vert=[]
#     weight_mid=[]
#     weight_cent=[]
#     if interior_nodes!=[]
#         param_S31, param_S22, param_S211, param_S1111, param_centroid, weight_S31, weight_S22, weight_S211, weight_S1111, weight_centroid= get_bary_coord_interior_tet(interior_nodes, vec(weights), vtx)        
        
#         weight_vert=wlgl[1]^3      
#         if convert(Bool,mod(p,2))
#             weight_mid=wlgl[convert(Int, floor(p/2)+1)]^2*wlgl[end]
#             # weight_cent=wlgl[convert(Int, floor(p/2)+1)]^3
#         end
#         # weight_edge= wlgl[2:end-convert(Int,centroid)].*wlgl[1]^2
#         # weights = convert.(T,collect(Iterators.flatten([J.*weight_vert, weight_S31, J.*weight_mid, weight_S22, 
#         #                                                 J.*weight_S211, 
#         #                                                 J.*weight_S1111, J.*weight_edge, J.*weight_cent])))
#     end

#     vertices = false #true
#     midedges = false #convert(Bool,mod(p,2))
#     numedge = 0 #convert(Int,(p-mod(p,2))/2)
#     numfaceS21 = 0 #convert(Int,(1+mod(p,2))*numedge)
#     numfaceS111 = 0 #convert(Int,1/2*(numedge^2 - numedge))
#     facecentroid = false #convert(Bool,mod(p,2))
#     numS31 = length(param_S31)
#     numS22 = length(param_S22)
#     numS211 = convert(Int, length(param_S211)/2)
#     numS1111 = convert(Int, length(param_S1111)/3)
#     centroid = convert(Bool,mod(p,2))

#     cub = SymCubatures.TetSymCub{T}(vertices=vertices,
#                                     midedges=midedges,
#                                     numfaceS21=numfaceS21,
#                                     numedge=numedge,
#                                     numfaceS111=numfaceS111,
#                                     facecentroid = facecentroid,
#                                     numS31=numS31,
#                                     numS22=numS22,
#                                     numS211=numS211,
#                                     numS1111=numS1111,
#                                     centroid=centroid)
    
#     # set the edge parameters 
#     init_tet_cub(cub, xedge)
#     faceS21_params = cub_tri.params[1 : cub_tri.numS21]
#     faceS21_w = wlgl[1].*cub_tri.weights[convert(Int,cub.vertices)+convert(Int,cub.midedges)+1 : convert(Int,cub.vertices)+convert(Int,cub.midedges)+cub_tri.numS21]
#     edge_params = cub_tri.params[cub_tri.numS21+1 : cub_tri.numS21+cub_tri.numedge]
#     edge_w = wlgl[1].*cub_tri.weights[convert(Int,cub.vertices)+convert(Int,cub.midedges)+cub.numfaceS21+1 : convert(Int,cub.vertices)+convert(Int,cub.midedges)+cub.numfaceS21+cub.numedge]
#     faceS111_params = cub_tri.params[cub_tri.numS21+cub_tri.numedge+1 : cub_tri.numS21+cub_tri.numedge+2*cub_tri.numS111]
#     faceS111_w =  wlgl[1].*cub_tri.weights[convert(Int,cub.vertices)+convert(Int,cub.midedges)+cub.numfaceS21+cub.numedge+1 : convert(Int,cub.vertices)+convert(Int,cub.midedges)+cub.numfaceS21+cub.numedge+cub.numfaceS111]
#     facecentroid_w = [] 
#     if mod(p,2)==1
#         facecentroid_w=wlgl[1].*cub_tri.weights[end]
#     end
#     params = convert.(T,collect(Iterators.flatten([3.0 .*param_S31,2.0 .*param_S22,  
#                                                    2.0*param_S211, 2.0 .*param_S1111])))
    
#     # p=7 works with J=1/800                                               
#     J = 1/4000
#     Jf = 1/10                    
#     weights = J.*convert.(T,collect(Iterators.flatten([weight_S31,  weight_S22, 
#                                                        weight_S211, 
#                                                        weight_S1111, weight_centroid])))
                           
#     println(params,"\n")
#     println(weights)                                                
#     if params != []
#         SymCubatures.setparams!(cub, params)
#         SymCubatures.setweights!(cub, weights)
#     end
#     # xy = SymCubatures.calcnodes(cub,vtx)
#     # SummationByParts.plotly_tet_nodes(q=q, x=xy, vtx=vtx)

#     xinit = convert.(T,collect(Iterators.flatten([cub.params,cub.weights])))
#     mask = Int[]
#     append!(mask, 1 : (cub.numparams+cub.numweights))  
#     Cubature.solvecubature!(cub, q, mask, tol=5e-14, hist=true, verbose=true, xinit=xinit, delta1=1e-4, delta2=1e-4)
#     println("\n", cub.params,"\n")
#     println(cub.weights,"\n")
#     return cub
# end

function init_tet_nodes_omega(p::Int,q::Int; ne::Int=-1, T=Float64)

    # use equilateral triangle
    node_tol = 1e-14
    vtx = T[1 -sqrt(3)/3 -sqrt(6)/6;
            0 2*sqrt(3)/3 -sqrt(6)/6;
            -1 -sqrt(3)/3 -sqrt(6)/6;
            0 0 sqrt(6)/2]
    if ne==-1
        ne=p+1
    end
    qf = 2*ne-1
    cub_lg, _ = SummationByParts.Cubature.quadrature(qf, internal=true)

    nn = convert(Int,(ne-0-mod(ne-0,2))/2)
    cub_tri = SymCubatures.TriSymCub{T}(vertices=false,
                                    midedges= false ,
                                    numS21=convert(Int,(1+mod(ne-0,2))*nn),
                                    numedge=0,
                                    numS111=convert(Int,1/2*(nn^2 - nn)),
                                    centroid = convert(Bool,mod(ne-0,2)))
    cub_tri = SummationByParts.init_tri_nodes_omega(p, q)
    init_tri_cub(cub_tri, cub_lg.params)
    xedge = cub_tri.params
    wlg = cub_lg.weights

    # get interior nodes
    interior_nodes, weights = get_interior_points_tet_omega(p)

    # get barycentric coordinates of the interior nodes
    params=[]
    weight_cent=[]
    if interior_nodes!=[]
        param_S31, param_S22, param_S211, param_S1111, param_centroid, weight_S31, weight_S22, weight_S211, weight_S1111, weight_centroid= get_bary_coord_interior_tet(interior_nodes, vec(weights), vtx)        

        if convert(Bool,mod(ne,2))
            weight_cent=wlg[convert(Int, floor(p/2)+1)]^3
        end
        J=1
        weights = convert.(T,collect(Iterators.flatten([weight_S31, weight_S22, 
                                                        J.*weight_S211, 
                                                        J.*weight_S1111, J.*weight_cent])))
    end

    vertices = false
    midedges = false #convert(Bool,mod(p,2))
    numedge = 0 #convert(Int,(p-mod(p,2))/2)
    numfaceS21 = 0 #convert(Int,(1+mod(p,2))*nn)
    numfaceS111 = 0 #convert(Int,1/2*(nn^2 - nn))
    facecentroid = false #convert(Bool,mod(p,2))
    numS31 = length(param_S31)
    numS22 = length(param_S22)
    numS211 = convert(Int, length(param_S211)/2)
    numS1111 = convert(Int, length(param_S1111)/3)
    centroid = convert(Bool,mod(ne,2))

    cub = SymCubatures.TetSymCub{T}(vertices=vertices,
                                    midedges=midedges,
                                    numfaceS21=numfaceS21,
                                    numedge=numedge,
                                    numfaceS111=numfaceS111,
                                    facecentroid = facecentroid,
                                    numS31=numS31,
                                    numS22=numS22,
                                    numS211=numS211,
                                    numS1111=numS1111,
                                    centroid=centroid)
    
    # set the edge parameters 
    init_tet_cub(cub, xedge)
    params = convert.(T,collect(Iterators.flatten([3.0 .*param_S31,2.0 .*param_S22, 
                                                   2.0*param_S211, 2.0 .*param_S1111])))
    
    # p=7 works with J=1/800
    Jc = 1/1000                                           
    J = 1/10          
    weights = J.*convert.(T,collect(Iterators.flatten([J*weight_S31,
                                                       J*weight_S22, 
                                                       J*weight_S211, 
                                                       J*weight_S1111, 
                                                       Jc*weight_cent])))
                           
    println(params,"\n")
    println(weights)                                                
    if params != []
        SymCubatures.setparams!(cub, params)
        SymCubatures.setweights!(cub, weights)
    end
    # xy = SymCubatures.calcnodes(cub,vtx)
    # SummationByParts.plotly_tet_nodes(q=q, x=xy, vtx=vtx)
    # SummationByParts.plotly_tet_nodes(q=q, x=interior_nodes, vtx=vtx)

    xinit = convert.(T,collect(Iterators.flatten([cub.params,cub.weights])))
    mask = Int[]
    append!(mask, 1 : cub.numparams+cub.numweights) 
    Cubature.solvecubature!(cub, q, mask, tol=5e-14, hist=true, verbose=true, xinit=xinit, delta1=1e-6, delta2=1e-6)
    println("\n", cub.params,"\n")
    println(cub.weights,"\n")
    return cub
end

# function init_tri_nodes(p::Int,q::Int; T=Float64)

#     # use equilateral triangle
#     node_tol = 1e-14
#     # vtx = T[0 0; 1 0; 1/2 sqrt(3/4)]
#     # vtx = T[-1 -1; 1 -1; 0 sqrt(3)-1]
#     vtx = T[-1 -1/sqrt(3); 1 -1/sqrt(3); 0 2/sqrt(3)] #T[-1 1-sqrt(3); 1 1-sqrt(3); 0 1] 
#     # p = convert(Int64,ceil(q/2))
#     # qf = q+mod(q,2)+1 
#     qf = 2*(p+2)-3
#     cub_lgl, _ = SummationByParts.Cubature.quadrature(qf, internal=false)
#     xedge = cub_lgl.params

#     vertices = true
#     midedges = convert(Bool,mod(p,2))
#     numedge = convert(Int,(p-mod(p,2))/2)
#     numS21 = convert(Int,(1+mod(p,2))*numedge)
#     numS111 = convert(Int,1/2*(numedge^2 - numedge))
#     centroid = convert(Bool,mod(p,2))

#     cub = SymCubatures.TriSymCub{T}(vertices=vertices,
#                                     midedges=midedges,
#                                     numS21=numS21,
#                                     numedge=numedge,
#                                     numS111=numS111,
#                                     centroid = centroid)
    
#     # set the edge parameters 
#     init_tri_cub(cub, xedge)

#     # find indices of the edge nodes of each interior line
#     numlines = cub.numedge*3 + 3*convert(Int,cub.centroid)
#     line_node_idx = zeros(Int64, 2, numlines)
#     numedge_idx0 = 3*convert(Int,vertices)+3*convert(Int,midedges)+3*numS21
#     # numedge_idx1 = numedge_idx0 + 6*numedge
#     k = numedge_idx0
#     for i=0:cub.numedge-1
#         for j=1:3
#             line_node_idx[1,i*3+j]=k+j
#             line_node_idx[2,i*3+j]=k+j+3
#         end
#         k+=6
#     end
#     if midedges
#         line_node_idx[1,cub.numedge*3+1:end] = Int[4, 5, 6]
#         line_node_idx[2,cub.numedge*3+1:end] = Int[3, 1, 2]
#     end

#     # find the intersection of each interior line and add to the nodes
#     xy = SymCubatures.calcnodes(cub, vtx)

#     num_int_nodes = convert(Int,centroid)+numS21*3+numS111*6
#     interior_nodes = zeros(T, (2,num_int_nodes))
#     k = 1
#     for i=1:numlines-1
#         for j=i+1:numlines
#             int_node = get_intersection(xy[:,line_node_idx[1,i]],xy[:,line_node_idx[2,i]], 
#                                         xy[:,line_node_idx[1,j]],xy[:,line_node_idx[2,j]],
#                                         vtx=vtx, node_tol=node_tol)
#             if norm(int_node)!= 0.0 
#                 if k <= num_int_nodes
#                     if sum(map(col -> norm(col), eachcol(interior_nodes[:, 1:k].-int_node)) .>= node_tol)==k
#                         interior_nodes[:,k]=int_node
#                         k+=1
#                     end
#                 end
#             end
#         end
#     end

#     # get barycentric coordinates of the interior nodes
#     param_S21, param_S111, param_centroid = get_bary_coord_interior(interior_nodes, vtx)
#     params = convert.(T,collect(Iterators.flatten([2.0 .*param_S21,cub_lgl.params,2.0 .*param_S111])))
#     if params != []
#         SymCubatures.setparams!(cub, params)
#         redistribute_nodes(cub,p)
#     end
#     xy = SymCubatures.calcnodes(cub,vtx)
#     SummationByParts.plot_tri_nodes(x=xy, vtx=vtx, q=q, n=cub.numnodes,write_title=true,label_nodes=true)

#     # xy_bent = bend_lines(cub, p, vtx)
#     # curve the lines 
#     # xy_int_bent = warp(interior_nodes, p, vtx, alpha=0.0)
#     # param_S21, param_S111, param_centroid = get_bary_coord_interior(xy_int_bent, vtx)
#     # params = convert.(T,collect(Iterators.flatten([2.0 .*param_S21,cub_lgl.params,2.0 .*param_S111])))
#     # if params != []
#     #     SymCubatures.setparams!(cub, params)
#     # end
#     # xy_bent = SymCubatures.calcnodes(cub,vtx)

#     # SummationByParts.plot_tri_nodes(x=xy_bent, vtx=vtx, q=q, n=cub.numnodes,write_title=true,label_nodes=true)

#     xinit = convert.(T,collect(Iterators.flatten([cub.params,cub.weights])))
#     mask = SymCubatures.getInternalParamMask(cub)
#     append!(mask, (cub.numparams+1):(cub.numparams+cub.numweights))  
#     Cubature.solvecubature!(cub, q, mask, tol=1e-14, hist=true, verbose=true, xinit=xinit, delta1=1e-4, delta2=1e-4)
    
#     return cub
# end

function get_bary_coord_interior(interior_nodes::Array{T,2},weights::Array{T,1},vtx::Array{T,2};node_tol::T=1e-14) where {T}
    # get barycentric coordinates of the interior nodes
    vtx1 = vcat(vtx',ones(T,(1,3)))
    interior_nodes1 = vcat(interior_nodes, ones(T,(1,size(interior_nodes,2))))
    bary_coord = vtx1\interior_nodes1
    
    param_S21=[]
    param_S111=[]
    weight_S21=[]
    weight_S111=[]

    for i=1:size(bary_coord,2)
        if ((abs(bary_coord[1,i]-bary_coord[2,i])<node_tol)||(abs(bary_coord[1,i]-bary_coord[3,i])<node_tol))
            push!(param_S21, bary_coord[1,i])
            push!(weight_S21, weights[i])
        elseif (abs(bary_coord[1,i] + bary_coord[2,i] - 1.0) <= node_tol)
            push!(param_S111, bary_coord[1,i])
            push!(param_S111, bary_coord[2,i])
            push!(weight_S111, weights[i])
        elseif (abs(bary_coord[1,i] + bary_coord[3,i] - 1.0) <= node_tol)
            push!(param_S111, bary_coord[1,i])
            push!(param_S111, bary_coord[3,i])
            push!(weight_S111, weights[i])
        else 
            push!(param_S111, bary_coord[2,i])
            push!(param_S111, bary_coord[3,i])
            push!(weight_S111, weights[i])
        end
    end

    return param_S21, param_S111, weight_S21, weight_S111
end

function get_bary_coord_interior_tet(interior_nodes::Array{T,2},weights::Array{T,1},vtx::Array{T,2};node_tol::T=1e-14) where {T}
    # get barycentric coordinates of the interior nodes
    vtx1 = vcat(vtx',ones(T,(1,4)))
    interior_nodes1 = vcat(interior_nodes, ones(T,(1,size(interior_nodes,2))))
    bary_coord = vtx1\interior_nodes1
    
    param_S31=[]
    param_S22=[]
    param_S211=[]
    param_S1111 =[]
    param_centroid =[]

    weight_S31=[]
    weight_S22=[]
    weight_S211=[]
    weight_S1111=[]
    weight_centroid=[]

    for i=1:size(bary_coord,2)
        if (any(sum(isapprox.(bary_coord[:, i], bary_coord[1, i], atol=node_tol)) .== 3) || 
            any(sum(isapprox.(bary_coord[:, i], bary_coord[2, i], atol=node_tol)) .== 3))
            push!(param_S31, bary_coord[1,i])
            push!(weight_S31, weights[i])
        elseif (any(sum(isapprox.(bary_coord[:, i], bary_coord[1, i], atol=node_tol)) .== 2) &&
                any(sum(isapprox.(bary_coord[:, i], bary_coord[2, i], atol=node_tol)) .== 2) &&
                any(sum(isapprox.(bary_coord[:, i], bary_coord[3, i], atol=node_tol)) .== 2))
            push!(param_S22, bary_coord[1,i])
            push!(weight_S22, weights[i])
        elseif (any(sum(isapprox.(bary_coord[:, i], bary_coord[1, i], atol=node_tol)) .== 2) ||
                any(sum(isapprox.(bary_coord[:, i], bary_coord[2, i], atol=node_tol)) .== 2) ||
                any(sum(isapprox.(bary_coord[:, i], bary_coord[3, i], atol=node_tol)) .== 2))
            idx = []
            for j=1:3
                for k=j+1:4
                    if isapprox(bary_coord[j,i], bary_coord[k, i], atol=node_tol)
                        append!(idx, j)
                        append!(idx, k)
                    end
                end
            end
            alpha_idx = idx[1]
            beta_idx = collect(setdiff(Set(1:4),idx))[1]
            push!(param_S211, bary_coord[alpha_idx,i])
            push!(param_S211, bary_coord[beta_idx,i])
            push!(weight_S211, weights[i])
        elseif any(sum(isapprox.(bary_coord[:, i], bary_coord[1, i], atol=node_tol)) .== 4)
            push!(param_centroid, bary_coord[1,i])
            push!(weight_centroid, weights[i])
        else 
            push!(param_S1111, bary_coord[1,i])
            push!(param_S1111, bary_coord[2,i])
            push!(param_S1111, bary_coord[3,i])
            push!(weight_S1111, weights[i])
        end
    end

    return  param_S31, param_S22, param_S211, param_S1111, param_centroid, weight_S31, weight_S22, weight_S211, weight_S1111, weight_centroid 
end

# function get_bary_coord_interior(interior_nodes::Array{T,2},vtx::Array{T,2}; node_tol::T=1e-14) where {T}
#     # get barycentric coordinates of the interior nodes
#     vtx1 = vcat(vtx',ones(T,(1,3)))
#     interior_nodes1 = vcat(interior_nodes, ones(T,(1,size(interior_nodes,2))))
#     bary_coord = vtx1\interior_nodes1
#     for i in 1:size(bary_coord)[2]
#         bary_coord[:, i] = sort(bary_coord[:, i])
#     end
#     int_S21 = []
#     int_S111 = []
#     int_centroid = []
#     for i = 1:size(bary_coord)[2]
#         idx = sort(findall(map(col -> norm(col), eachcol(bary_coord .- bary_coord[:,i])) .<= node_tol))
#         if length(idx)==1
#             if idx[1] ∉ int_centroid
#                 push!(int_centroid,idx[1])
#             end
#         elseif length(idx)==3
#             if idx[1] ∉ int_S21
#                 push!(int_S21,idx[1])
#             end
#         elseif length(idx)==6
#             if idx[1] ∉ int_S111
#                 push!(int_S111,idx[1])
#             end
#         else
#             error("There is no matching symmetry group.")
#         end
#     end

#     param_S21=[]
#     param_S111=[]
#     param_centroid=[]
    
#     for i=1:length(int_S21)
#         if (abs(bary_coord[:,int_S21[i]][1] - bary_coord[:,int_S21[i]][2]) <= node_tol) || 
#            (abs(bary_coord[:,int_S21[i]][1] - bary_coord[:,int_S21[i]][3]) <= node_tol)
#             push!(param_S21, bary_coord[:,int_S21[i]][1])
#         else
#             push!(param_S21, bary_coord[:,int_S21[i]][2])
#         end
#     end
    
#     for i=1:length(int_S111)
#         if (abs(bary_coord[:,int_S111[i]][1] + bary_coord[:,int_S111[i]][2] - 1.0) <= node_tol)
#             push!(param_S111, bary_coord[:,int_S111[i]][1])
#             push!(param_S111, bary_coord[:,int_S111[i]][2])
#         elseif (abs(bary_coord[:,int_S111[i]][1] + bary_coord[:,int_S111[i]][3] - 1.0) <= node_tol)
#             push!(param_S111, bary_coord[:,int_S111[i]][1])
#             push!(param_S111, bary_coord[:,int_S111[i]][3])
#         else 
#             push!(param_S111, bary_coord[:,int_S111[i]][2])
#             push!(param_S111, bary_coord[:,int_S111[i]][3])
#         end
#     end

#     if int_centroid != []
#         push!(param_centroid, bary_coord[1,int_centroid[1]])
#     end

#     return param_S21, param_S111, param_centroid
# end

function init_tri_cub(cub::TriSymCub{T}, xedge::Array{T,1}) where {T}
    xinit = [] 
    xinit_sym_group=[]
    if cub.vertices==true
        push!(xinit_sym_group,"vertices")
    end
    if cub.midedges==true
        push!(xinit_sym_group,"midedges")
    end
    if cub.numS21 != 0
        push!(xinit_sym_group,"numS21")
    end
    if cub.numedge != 0
        push!(xinit_sym_group,"numedge")
    end
    if cub.numS111 != 0
        push!(xinit_sym_group,"numS111")
    end
    if cub.centroid==true
        push!(xinit_sym_group,"centroid")
    end

    xedge_sym_group=[]
    if cub.vertices==true
        push!(xedge_sym_group,"vertices")
    end
    if cub.midedges==true
        push!(xedge_sym_group,"midedges")
    end
    if cub.numedge != 0
        push!(xedge_sym_group,"numedge")
    end
    # compute the number of parameters and weights
    numparams = cub.numparams
    numweights = cub.numweights

    xinit = 0.1 .* ones(numparams+numweights, 1)
    sym_group = ["vertices", "midedges", "numS21", "numedge", "numS111", "centroid"]
    param_dict = OrderedDict{String, Vector{Float64}}()
    weight_dict = OrderedDict{String, Vector{Float64}}()

    # sort to match symmetry groups in xinit
    p_loc = [0]
    w_loc = [numparams]
    for s in xinit_sym_group
        if s == sym_group[1]
            push!(w_loc, w_loc[end]+1)
        elseif s ==sym_group[2]
            push!(w_loc, w_loc[end]+1)
        elseif s ==sym_group[3]
            push!(p_loc, p_loc[end]+cub.numS21)
            push!(w_loc, w_loc[end]+cub.numS21)        
        elseif s ==sym_group[4]
            push!(p_loc, p_loc[end]+cub.numedge)
            push!(w_loc, w_loc[end]+cub.numedge)
        elseif s==sym_group[5]
            push!(p_loc, p_loc[end]+2*cub.numS111)
            push!(w_loc, w_loc[end]+cub.numS111)
        elseif s==sym_group[6]
            push!(w_loc, w_loc[end]+1)
        end
    end

    p_cnt = 1
    for i=eachindex(xinit_sym_group)
        if (xinit_sym_group[i] ∉ ["vertices", "midedges", "centroid"])
            param_dict[xinit_sym_group[i]] = xinit[p_loc[p_cnt]+1:p_loc[p_cnt+1]]
            p_cnt += 1
        end
    end

    for i=2:length(xinit_sym_group)+1
        weight_dict[xinit_sym_group[i-1]] = xinit[w_loc[i-1]+1:w_loc[i]]
    end

    param_sorted = OrderedDict(key => param_dict[key] for key in sym_group if haskey(param_dict, key))
    weight_sorted = OrderedDict(key => weight_dict[key] for key in sym_group if haskey(weight_dict, key))
    xinit_param = values(param_sorted)
    xinit_weight = values(weight_sorted)
    xinit = collect(Iterators.flatten(collect(Iterators.flatten(hcat(xinit_param, xinit_weight)))))

    # sort to match symmetry groups in xedge
    p_loc = [0]
    for s in xedge_sym_group
        if s ==sym_group[4]
            push!(p_loc, p_loc[end]+cub.numedge)
        end
    end

    param_dict_edge = OrderedDict{String, Vector{Float64}}()
    p_cnt = 1
    for i=eachindex(xedge_sym_group)
        if (xedge_sym_group[i] ∉ ["vertices", "midedges", "centroid"])
            param_dict_edge[xedge_sym_group[i]] = xedge[p_loc[p_cnt]+1:p_loc[p_cnt+1]]
            p_cnt += 1
        end
    end

    for i=eachindex(xedge_sym_group)
        if (xedge_sym_group[i] ∉ ["vertices", "midedges", "centroid"])
            param_dict[xedge_sym_group[i]] = param_dict_edge[xedge_sym_group[i]]
        end
    end

    param_sorted = OrderedDict(key => param_dict[key] for key in sym_group if haskey(param_dict, key))
    weight_sorted = OrderedDict(key => weight_dict[key] for key in sym_group if haskey(weight_dict, key))
    xinit_param = values(param_sorted)
    xinit_weight = values(weight_sorted)
    # xinit = collect(Iterators.flatten(collect(Iterators.flatten(hcat(xinit_param, xinit_weight)))))

    cub.params = collect(Iterators.flatten(xinit_param))
    cub.weights = collect(Iterators.flatten(xinit_weight))

    return
end

function init_tet_cub(cub::TetSymCub{T}, xedge::Array{T,1}) where {T}
    xinit = [] 
    xinit_sym_group=[]
    if xinit_sym_group==[]
        if cub.vertices==true
            push!(xinit_sym_group,"vertices")
        end
        if cub.numS31!=0
            push!(xinit_sym_group,"numS31")
        end
        if cub.midedges==true
            push!(xinit_sym_group,"midedges")
        end
        if cub.numS22 != 0
            push!(xinit_sym_group,"numS22")
        end
        if cub.numfaceS21 != 0
            push!(xinit_sym_group,"numfaceS21")
        end
        if cub.numedge != 0
            push!(xinit_sym_group,"numedge")
        end
        if cub.numS211 != 0
            push!(xinit_sym_group,"numS211")
        end
        if cub.numfaceS111 != 0
            push!(xinit_sym_group,"numfaceS111")
        end
        if cub.facecentroid != 0
            push!(xinit_sym_group,"facecentroid")
        end
        if cub.numS1111 != 0
            push!(xinit_sym_group,"numS1111")
        end
        if cub.centroid==true
            push!(xinit_sym_group,"centroid")
        end
    end

    xedge_sym_group=[]
    if xedge_sym_group==[]
        if cub.vertices==true
            push!(xedge_sym_group,"vertices")
        end
        if cub.midedges==true
            push!(xedge_sym_group,"midedges")
        end
        if cub.numfaceS21 != 0
            push!(xedge_sym_group,"numfaceS21")
        end
        if cub.numedge != 0
            push!(xedge_sym_group,"numedge")
        end
        if cub.numfaceS111 != 0
            push!(xedge_sym_group,"numfaceS111")
        end
        if cub.facecentroid != 0
            push!(xedge_sym_group,"facecentroid")
        end
    end

    # compute the number of parameters and weights
    numparams = cub.numparams
    numweights = cub.numweights

    xinit = 0.1 .* ones(numparams+numweights, 1)
    sym_group =  ["vertices", "numS31", "midedges", "numS22", "numfaceS21", "numedge", 
                 "numS211", "numfaceS111", "facecentroid", "numS1111", "centroid"]
    param_dict = OrderedDict{String, Vector{Float64}}()
    weight_dict = OrderedDict{String, Vector{Float64}}()

    # sort to match symmetry groups in xinit
    p_loc = [0]
    w_loc = [numparams]
    for s in xinit_sym_group
        if s == sym_group[1]
            push!(w_loc, w_loc[end]+1)
        elseif s ==sym_group[2]
            push!(p_loc, p_loc[end]+cub.numS31)
            push!(w_loc, w_loc[end]+cub.numS31)
        elseif s ==sym_group[3]
            push!(w_loc, w_loc[end]+1)        
        elseif s ==sym_group[4]
            push!(p_loc, p_loc[end]+cub.numS22)
            push!(w_loc, w_loc[end]+cub.numS22)
        elseif s==sym_group[5]
            push!(p_loc, p_loc[end]+cub.numfaceS21)
            push!(w_loc, w_loc[end]+cub.numfaceS21)
        elseif s==sym_group[6]
            push!(p_loc, p_loc[end]+cub.numedge)
            push!(w_loc, w_loc[end]+cub.numedge)
        elseif s==sym_group[7]
            push!(p_loc, p_loc[end]+2*cub.numS211)
            push!(w_loc, w_loc[end]+cub.numS211)
        elseif s==sym_group[8]
            push!(p_loc, p_loc[end]+2*cub.numfaceS111)
            push!(w_loc, w_loc[end]+cub.numfaceS111)
        elseif s==sym_group[9]
            push!(w_loc, w_loc[end]+1)
        elseif s==sym_group[10]
            push!(p_loc, p_loc[end]+3*cub.numS1111)
            push!(w_loc, w_loc[end]+cub.numS1111)
        elseif s==sym_group[11]
            push!(w_loc, w_loc[end]+1)
        end
    end

    p_cnt = 1
    for i=eachindex(xinit_sym_group)
        if (xinit_sym_group[i] ∉ ["vertices", "midedges", "facecentroid", "centroid"])
            param_dict[xinit_sym_group[i]] = xinit[p_loc[p_cnt]+1:p_loc[p_cnt+1]]
            p_cnt += 1
        end
    end

    for i=2:length(xinit_sym_group)+1
        weight_dict[xinit_sym_group[i-1]] = xinit[w_loc[i-1]+1:w_loc[i]]
    end

    param_sorted = OrderedDict(key => param_dict[key] for key in sym_group if haskey(param_dict, key))
    weight_sorted = OrderedDict(key => weight_dict[key] for key in sym_group if haskey(weight_dict, key))
    xinit_param = values(param_sorted)
    xinit_weight = values(weight_sorted)
    xinit = collect(Iterators.flatten(collect(Iterators.flatten(hcat(xinit_param, xinit_weight)))))

    # sort to match symmetry groups in xedge
    p_loc = [0]
    for s in xedge_sym_group
        if s ==sym_group[5]
            push!(p_loc, p_loc[end]+cub.numfaceS21)        
        elseif s ==sym_group[6]
            push!(p_loc, p_loc[end]+cub.numedge)
        elseif s==sym_group[8]
            push!(p_loc, p_loc[end]+2*cub.numfaceS111)
        end
    end

    param_dict_edge = OrderedDict{String, Vector{Float64}}()
    p_cnt = 1
    for i=eachindex(xedge_sym_group)
        if (xedge_sym_group[i] ∉ ["vertices", "midedges", "facecentroid","centroid"])
            param_dict_edge[xedge_sym_group[i]] = xedge[p_loc[p_cnt]+1:p_loc[p_cnt+1]]
            p_cnt += 1
        end
    end

    for i=eachindex(xedge_sym_group)
        if (xedge_sym_group[i] ∉ ["vertices", "midedges","facecentroid","centroid"])
            param_dict[xedge_sym_group[i]] = param_dict_edge[xedge_sym_group[i]]
        end
    end

    param_sorted = OrderedDict(key => param_dict[key] for key in sym_group if haskey(param_dict, key))
    weight_sorted = OrderedDict(key => weight_dict[key] for key in sym_group if haskey(weight_dict, key))
    xinit_param = values(param_sorted)
    xinit_weight = values(weight_sorted)
    # xinit = collect(Iterators.flatten(collect(Iterators.flatten(hcat(xinit_param, xinit_weight)))))

    cub.params = collect(Iterators.flatten(xinit_param))
    cub.weights = collect(Iterators.flatten(xinit_weight))

    return
end

function get_intersection(node1::Array{T,1}, node2::Array{T,1}, 
                          node3::Array{T,1}, node4::Array{T,1};
                          vtx::Array{T,2}=T[-1 -1; 1 -1; 0 sqrt(3)],
                          node_tol::T=1e-14) where {T}

    vtx_equilateral = [0 0; 1 0; 1/2 sqrt(3/4)]
    vtx_equilateral_big = [-1 -1/sqrt(3); 1 -1/sqrt(3); 0 2/sqrt(3)]#[-1 1-sqrt(3); 1 1-sqrt(3); 0 1]

    node = zeros(T,(2,1))
    if abs((node2[1]-node1[1])) <= node_tol
        node[1] = node1[1]
        node[2] = ((node4[2]-node3[2])/(node4[1]-node3[1]))*(node[1]-node3[1]) + node3[2]
    elseif abs((node4[1]-node3[1])) <= node_tol
        node[1] = node3[1]
        node[2] = ((node2[2]-node1[2])/(node2[1]-node1[1]))*(node[1]-node1[1]) + node1[2]
    else
        m1 = (node2[2]-node1[2])/(node2[1]-node1[1])
        m2 = (node4[2]-node3[2])/(node4[1]-node3[1])
    
        if m1 == m2
            return node
        else
            node[1] = 1/(m1-m2)*(m1*node1[1] - m2*node3[1] + node3[2] - node1[2])
            node[2] = m1*(node[1]-node1[1]) + node1[2]
        end
        
        if norm(vtx .- vtx_equilateral) <= 1e-13
            if ((node[2] > -sqrt(3)*(node[1]-1)) || (node[2] > sqrt(3)*node[1]) || node[2]<0 || node[1] < 0 || node[1]>1)
                return 0.0 .*node
            end
        elseif norm(vtx .- vtx_equilateral_big) <= 1e-13
            # if ((node[2] > sqrt(3)*(1-node[1])-1) || (node[2] > sqrt(3)*(node[1]+1)-1) || node[2]<-1 || node[1] < -1 || node[1]>1)
            if ((node[2] > -sqrt(3)*node[1]+sqrt(3)-1/sqrt(3)) || 
                (node[2] > sqrt(3)*node[1]+sqrt(3)-1/sqrt(3)) || node[2]< -1/sqrt(3) || 
                node[1] < -1 || node[1]>1)
                return 0.0 .*node
            end
        end
    end
    return node
end 

function redistribute_nodes(cub::TriSymCub{T}, p::Int) where {T}
    # Assumes equilateral triangle with centroid at (0,0)
    qf = 2*(p+2)-3 
    cub_lgl, vtx_lgl = SummationByParts.Cubature.quadrature(qf, internal=false)
    dist = sort(vec(SymCubatures.calcnodes(cub_lgl, Matrix(reshape([0.0;1],2,1)))'))

    param_S21 = sort(vec(cub.params[1:cub.numS21])./2)
    numedge = cub.numedge
    param_S111 = vec(cub.params[length(param_S21)+numedge+1 : end])./2
    n_above_cent = cub.numedge
    n_below_cent = 0
    if mod(p,2)==1
        n_below_cent=numedge
    end
    # dist = range(0.0,0.8,numedge+1)
    
    k=2
    for i=1:n_above_cent
        # if param_S21[i] > 1/3
        param_S21[i] = (1-0.2)*param_S21[i] #dist[k] * 3/4#*1/3
        # end
        k+=1
    end

    if n_below_cent > 0
        # k=n_above_cent
        for i=n_above_cent+1:length(param_S21)
            param_S21[i] = (1+0.8*(0.5-param_S21[i]))*param_S21[i] #dist[k] * #* 1/3 #*1/3
            # k-=1
        end
    end

    m=1
    for i=1:convert(Int,length(param_S111)/2)
        d = 0.1*param_S111[m+1]
        param_S111[m] = param_S111[m]*(1-d)
        param_S111[m+1] = param_S111[m+1]*(1+d)
        m+=2
    end

    params = cub.params
    params[1:length(param_S21)] = param_S21.*2
    params[length(param_S21)+numedge+1 : end] = param_S111.*2
    params = convert.(T,collect(Iterators.flatten(params)))
    SymCubatures.setparams!(cub,params)
end

function bend_lines(cub::TriSymCub{T}, p::Int,vtx::Array{T,2}) where {T}
    q = 2*p
    lines = sparse_stencil(cub, p)
    xy = SymCubatures.calcnodes(cub,vtx)
    xy_bent = zeros(T,size(xy))
    numlines = size(lines,1)
    nperline = size(lines,2)
    cent = [0.5,1/3]
    xh = 0.0
    yh = 0.0
    h = 0.0
    for i = 4:numlines
        line = xy[:,lines[i,:]]
        n1 = xy[:,lines[i,argmin(line[1,:])]]
        n2 = xy[:,lines[i,argmax(line[1,:])]]
        m = (n2[2]-n1[2])/(n2[1]-n1[1])
        m_orth = -1/m
        xc = 0.5*(n1[1]+n2[2])
        yc = m*(xc - n1[1]) + n1[2]

        if (m > 1e-13 || m < 1e-13)
            if cent[2] < yc
                h = 1/q
            else
                h = sqrt((cent[2] - yc)^2 + (cent[1] - xc)^2) 
                h = (1+1/q)*h
            end
        else
            h = 1/q
        end

        if (m > 1e-13) #&& (n1[2] > cent[2] || n2[2] >cent[2]) 
            xh = xc - sqrt(h^2/(1+m_orth^2)) 
            yh = m_orth*(xh-xc)+yc
        elseif (m < -1e-13)
            xh = xc + sqrt(h^2/(1+m_orth^2)) 
            yh = m_orth*(xh-xc)+yc
        else
            xh = xc
            yh = yc - h*yc
        end

        A = zeros(T, (3,3))
        xx = Array([n1[1];n2[1];xh])
        yy = Array([n1[2];n2[2];yh])
        A[:,1] = xx.^2
        A[:,2] = xx.^1
        A[:,3] = xx.^0
        aa = A\yy
        
        # for j = 1:nperline
        #     x1 = -(aa[2]-m_orth) + sqrt((aa[2]-m_orth)^2 - 4*aa[1]*(aa[3]-xy[1,lines[i,j]]/m - xy[2,lines[i,j]]))
        #     x2 = -(aa[2]-m_orth) - sqrt((aa[2]-m_orth)^2 - 4*aa[1]*(aa[3]-xy[1,lines[i,j]]/m - xy[2,lines[i,j]]))
        #     if (m > 0) && (n1[2] > cent[2] || n2[2] >cent[2]) 
        #         xy_bent[1,lines[i,j]] = minimum([x1,x2])
        #         xy_bent[2,lines[i,j]] = aa[1]*xy_bent[1,lines[i,j]]^2 + aa[2]*xy_bent[1,lines[i,j]] + aa[3]
        #     elseif (m < 0) && (n1[2] > cent[2] || n2[2] >cent[2]) 
        #         xy_bent[1,lines[i,j]] = maximum([x1,x2])
        #         xy_bent[2,lines[i,j]] = aa[1]*xy_bent[1,lines[i,j]]^2 + aa[2]*xy_bent[1,lines[i,j]] + aa[3]
        #     end
        # end
        for j = 1:nperline
            if m > 1e-13
                
            end
        end
    end
    return xy_bent
end

function warp(interior_nodes::Array{T,2}, p::Int, vtx::Array{T,2}; alpha::T=5/3) where {T}

    vtx1 = vcat(vtx',ones(T,(1,3)))
    interior_nodes1 = vcat(interior_nodes, ones(T,(1,size(interior_nodes,2))))
    bary_coord = vtx1\interior_nodes1

    L1 = bary_coord'[:,1]
    L2 = bary_coord'[:,2]
    L3 = bary_coord'[:,3]

    b1 = 1 .* L2 .* L3
    b2 = 1 .* L1 .* L3
    b3 = 1 .* L1 .* L2

    warpf1 = warp_factor(p, vec(L3 .- L2))
    warpf2 = warp_factor(p, vec(L1 .- L3))
    warpf3 = warp_factor(p, vec(L2 .- L1))

    warp1 = b1 .* warpf1 .* (1 .+ (alpha .* L1).^2)
    warp2 = b2 .* warpf2 .* (1 .+ (alpha .* L2).^2)
    warp3 = b3 .* warpf3 .* (1 .+ (alpha .* L3).^2)

    xy = zeros(T, (2, length(L1)))
    xy[1,:] = -L2 .+ L3
    xy[2,:] = (-L2.-L3.+2 .*L1)./(sqrt(3))

    xy[1,:] = reshape(xy[1,:],(1,length(xy[1,:])))  .+ (1 .* warp1 .+ cos(2*pi/3).*warp2 .+ cos(4*pi/3).*warp3)'
    xy[2,:] = reshape(xy[2,:],(1,length(xy[2,:])))  .+ (0 .* warp1 .+ sin(2*pi/3).*warp2 .+ sin(4*pi/3).*warp3)'

    return xy
end

function warp_factor(p::Int, r::Array{T,1}) where{T}
    N = p+1
    # q = 2*p
    qf = 2*(p+2)-3 #q+mod(q,2)+1 
    cub_lgl, vtx_lgl = SummationByParts.Cubature.quadrature(qf, internal=false)
    LGLr = sort(vec(SymCubatures.calcnodes(cub_lgl, vtx_lgl)'))
    # get uniform nodes 
    if rem(N+1,2) == 0 
        centroid=false
    else
        centroid=true
    end
    numedge = div(N-1,2)
    cub_eq = LineSymCub{T}(vertices=true, centroid=centroid,numedge=numedge)
    alpha = zeros(T, (numedge))
    dx = 1.0/(N)
    for i = 1:numedge
    alpha[i] = i*dx 
    end    
    SymCubatures.setparams!(cub_eq, alpha)

    # cub_eq, vtx_eq = SummationByParts.Cubature.quadratureUniform(N, N+1, internal=false)
    req = sort(vec(SymCubatures.calcnodes(cub_eq, vtx_lgl)'))

    Veq,_ = OrthoPoly.vandermonde(N,collect(Iterators.flatten(req)))
    Nr = length(r)
    Pmat = zeros(N+1, Nr)
    for i=1:N+1
        Pmat[i,:] = OrthoPoly.jacobipoly(vec(r), 0.0, 0.0, i-1)'
    end
    Lmat = Veq'\Pmat

    warp = Lmat'*(LGLr .- req)

    zerof = (abs.(r) .< 1-1e-10)
    sf = 1 .- (zerof .* r).^2
    warp = warp./sf .+ warp.*(zerof .- 1)

    return warp
end


#--------------------------------
function init_tri_uni_nodes(p::Int, q::Int; T=Float64)

    # use equilateral triangle
    node_tol = 1e-14
    vtx = T[-1 -1/sqrt(3); 1 -1/sqrt(3); 0 2/sqrt(3)]
    # p = convert(Int64,ceil(q/2))
    # qf = q+mod(q,2)+1 
    qf = 2*(p+2)-3
    cub_lgl, _ = SummationByParts.Cubature.quadrature(qf, internal=false)
    N = convert(Int,(p+2)*(p+3)/2);

    xy = uni_warp_nodes(p, alpha=5.0/3.0)

    # find indices of the edge nodes
    ne = Int[]
    for i = 1:N
        if (abs(xy[2,i] - (-sqrt(3)*xy[1,i]+sqrt(3)-1/sqrt(3))) < node_tol || 
            abs(xy[2,i] - (sqrt(3)*xy[1,i]+sqrt(3)-1/sqrt(3))) < node_tol ||
            abs(xy[2,i] - (-1/sqrt(3))) < node_tol)
           append!(ne,i) 
        end
    end
    
    # check if centroid node is present
    nc=0
    for i = 1:N
        if (abs(xy[1,i])<=node_tol && abs(xy[2,i])<=node_tol)
            nc += 1
        end
        if nc > 1
            error("More than one centroid node identified.")
        end
    end

    int_node_idx = vec(setdiff(range(1,N), ne))
    interior_nodes = xy[:,int_node_idx]

    # get barycentric coordinates of the interior nodes
    param_S21, param_S111, param_centroid = get_bary_coord_interior(interior_nodes, vtx)

    vertices = true
    midedges = convert(Bool,mod(p,2))
    numedge = convert(Int,(p-mod(p,2))/2)
    numS21 = length(param_S21) #convert(Int,(1+mod(p,2))*numedge)
    numS111 = convert(Int,1/2*length(param_S111)) #convert(Int,1/2*(numedge^2 - numedge))
    centroid = convert(Bool,nc) #convert(Bool,mod(p+1,2))

    cub = SymCubatures.TriSymCub{T}(vertices=vertices,
                                    midedges=midedges,
                                    numS21=numS21,
                                    numedge=numedge,
                                    numS111=numS111,
                                    centroid = centroid)
    params = convert.(T,collect(Iterators.flatten([2.0 .*param_S21,cub_lgl.params,2.0 .*param_S111])))
    if params != []
        SymCubatures.setparams!(cub, params)
    end
    xy = SymCubatures.calcnodes(cub,vtx)
    # SummationByParts.plot_tri_nodes(x=xy, vtx=vtx, q=q, n=cub.numnodes,write_title=true,label_nodes=true)

    xinit = convert.(T,collect(Iterators.flatten([cub.params,cub.weights])))
    mask = SymCubatures.getInternalParamMask(cub)
    append!(mask, (cub.numparams+1):(cub.numparams+cub.numweights))  
    Cubature.solvecubature!(cub, q, mask, tol=1e-14, hist=true, verbose=true, xinit=xinit, delta1=1e-4, delta2=1e-4)
    
    return cub
end

function uni_warp_nodes(p::Int; alpha::T=5/3) where {T}
  
    N=p+1
    Np = convert(Int,(N+1)*(N+2)/2);

    L1 = zeros(T,(Np,1)); L2 = zeros(T,(Np,1)); L3 = zeros(T,(Np,1));
    sk = 1;
    for n=1:N+1
    for m=1:N+2-n
        L1[sk] = (n-1)/N; L3[sk] = (m-1)/N;
        sk = sk+1;
    end
    end
    L2 = 1.0 .-L1 .-L3;
    x = -L2 .+L3; y = (-L2.-L3 .+2 .*L1)./sqrt(3.0);
    xy = hcat(x,y)'

    b1 = 4 .* L2 .* L3
    b2 = 4 .* L1 .* L3
    b3 = 4 .* L1 .* L2

    warpf1 = warp_factor(p, vec(L3 .- L2))
    warpf2 = warp_factor(p, vec(L1 .- L3))
    warpf3 = warp_factor(p, vec(L2 .- L1))

    warp1 = b1 .* warpf1 .* (1 .+ (alpha .* L1).^2)
    warp2 = b2 .* warpf2 .* (1 .+ (alpha .* L2).^2)
    warp3 = b3 .* warpf3 .* (1 .+ (alpha .* L3).^2)

    xy = zeros(T, (2, length(L1)))
    xy[1,:] = -L2 .+ L3
    xy[2,:] = (-L2.-L3.+2 .*L1)./(sqrt(3))

    xy[1,:] = reshape(xy[1,:],(1,length(xy[1,:])))  .+ (1 .* warp1 .+ cos(2*pi/3).*warp2 .+ cos(4*pi/3).*warp3)'
    xy[2,:] = reshape(xy[2,:],(1,length(xy[2,:])))  .+ (0 .* warp1 .+ sin(2*pi/3).*warp2 .+ sin(4*pi/3).*warp3)'

    return xy
end

function get_cut_rect_to_equi_tri_mapping(p::Int; T=Float64)

    qf = 2*(p+2)-3
    cub_lgl, vtx_lgl = SummationByParts.Cubature.quadrature(qf, internal=false)
    xlgl = sort(vec(SymCubatures.calcnodes(cub_lgl, vtx_lgl)'))
    a = xlgl[convert(Int, ceil(length(xlgl)/2)-mod(p,2))]
    b = xlgl[convert(Int, ceil(length(xlgl)/2)+1)]

    xy_rect = Array(T[b 1 1 -1 -1 xlgl[2] xlgl[3] xlgl[end] xlgl[end] 0; 
                     -1 -1 1 1 b xlgl[end] xlgl[end] xlgl[2] xlgl[3] 0])'
    xy_tri = Array(T[b 1 0 -1 a (-1+1/2*(xlgl[2]+1)) (-1+1/2*(xlgl[3]+1)) (1-1/2*(1-xlgl[end-1])) (1-1/2*(1-xlgl[end-2])) 0; 
                  -1/sqrt(3) -1/sqrt(3) 2/sqrt(3) -1/sqrt(3) -1/sqrt(3) (-1/sqrt(3)+sqrt(3)/2*(xlgl[2]+1)) (-1/sqrt(3)+sqrt(3)/2*(xlgl[3]+1)) (-1/sqrt(3)+sqrt(3)/2*(1-xlgl[end-1])) (-1/sqrt(3)+sqrt(3)/2*(1-xlgl[end-2])) 0])'

    # xy_rect = Array(T[b 1 1 -1 -1 ; 
    #                  -1 -1 1 1 b ])'
    # xy_tri = Array(T[b 1 0 -1 a ; 
    #                 -1/sqrt(3) -1/sqrt(3) 2/sqrt(3) -1/sqrt(3) -1/sqrt(3)])'

    # v = Array([xy_rect[:,1].^2; xy_rect[:,2].^2; xy_rect[:,1]; xy_rect[:,2]; ones(T,(5,1))])'
    # v = reshape(v,(5,5))

    # v = Array([xy_rect[:,1].^2; xy_rect[:,2].^2; xy_rect[:,1]; xy_rect[:,2]; ones(T,(5,1))])'
    # v = reshape(v,(5,5))

    v,_ = OrthoPoly.vandermonde(2, xy_rect[:,1],xy_rect[:,2])
    v2 = []
    idx= Int[]
    for i=1:size(v,1)
        push!(v2, v[i,:])
        if rank(hcat(v2...)') < size(hcat(v2...)',1)
            pop!(v2)
        else
            append!(idx,i)
        end
    end
    v=v[idx,:]
    # v2 = copy(v[1:6,1:6])
    # for i = 7:size(xy_rect,1) 
    #     if rank(v2)<6
    #         v2[6,:] = v2[i]
    #     end
    # end
    # v,_ = OrthoPoly.vandermonde(2, xy_rect[:,1],xy_rect[:,2])
    # coef = pinv(v)*xy_tri
    coef = v\xy_tri[idx,:]

    return coef
end

function generate_cut_rect(p::Int; T=Float64)
    node_tol=1e-14
    n = p+2
    qf = 2*n-3
    cub_lgl, vtx_lgl = SummationByParts.Cubature.quadrature(qf, internal=false)
    perm = sortperm(vec(SymCubatures.calcnodes(cub_lgl, vtx_lgl)))
    xlgl = SymCubatures.calcnodes(cub_lgl, vtx_lgl)[perm]
    b = xlgl[convert(Int, ceil(length(xlgl)/2)+1)]
    centroid = convert(Int, cub_lgl.centroid)
    # wlgl = SymCubatures.calcweights(cub_lgl)[perm]

    xgrid = repeat(xlgl, outer=length(xlgl))
    ygrid = repeat(xlgl, inner=length(xlgl))
    rect = hcat(xgrid, ygrid)'

    ncut = convert(Int, n^2 - 1/4*(n^2 - centroid*(2*n-1)))
    cut_rect = zeros(T, (2,ncut))
    k = 1
    for i = 1:length(xgrid)
        if (rect[1,i] >= (b-node_tol) || rect[2,i] >= (b-node_tol)) 
            cut_rect[:,k] = rect[:,i]
            k+=1
        end
    end

    return cut_rect
end

function map_cut_rec_to_equi_tri(p::Int; T=Float64)

    c = get_cut_rect_to_equi_tri_mapping(p)
    xy_rect = generate_cut_rect(p)
    n = size(xy_rect,1)

    # v = Array([xy_rect[:,1].^2; xy_rect[:,2].^2; xy_rect[:,1]; xy_rect[:,2]; ones(T,(n,1))])'
    # v = reshape(v,(size(xy_rect,1),5))
    v,_ = OrthoPoly.vandermonde(2, xy_rect[1,:],xy_rect[2,:])

    xy_tri = (v*c)'

    q = 2*p
    vtx = T[-1 -1/sqrt(3); 1 -1/sqrt(3); 0 2/sqrt(3)]
    SummationByParts.plot_tri_nodes(x=xy_tri, vtx=vtx, q=q, n=n,write_title=true,label_nodes=true)
    return xy_tri
end

function get_interior_points_tri(p::Int; ne::Int=-1, T=Float64, node_tol=1e-14) 

    vtx = [-1 -1/sqrt(3); 1 -1/sqrt(3); 0 2/sqrt(3)]
    if ne==-1
        ne=p+2
    end
    points = []
    weights = T[]
    numlines = convert(Int, floor((ne-2)/2))
    facet1_lgl, wlgl_f1 = get_lgl_nodes_2d(p, Array(T[1,-1/sqrt(3)]), Array(T[-1,-1/sqrt(3)]), ne=ne, half=false)
    facet2_lgl, wlgl_f2 = get_lgl_nodes_2d(p, Array(T[1,-1/sqrt(3)]), Array(T[0,2/sqrt(3)]), ne=ne, half=false)
    a = 0/(2*p)
    centerline_1, wlgl_c1 = get_lgl_nodes_2d(p, Array(T[0,-1/sqrt(3)]), Array(T[0, a]), ne=ne, half=true)
    d=sqrt((1/2)^2+(sqrt(3)/2 - 1/sqrt(3))^2)
    c=a+d
    x1= 1/2 - c*cos(pi/6)
    y1= sqrt(3)/2-1/sqrt(3) - c*sin(pi/6)
    centerline_2, wlgl_c2 = get_lgl_nodes_2d(p, Array(T[1/2,sqrt(3)/2-1/sqrt(3)]), Array(T[x1,y1]), ne=ne, half=true)

    for i = 2:numlines+1
        for j = i:numlines+1
            point = get_intersection(facet2_lgl[:,i], centerline_1[:,i], 
                            facet1_lgl[:,j], centerline_2[:,j], vtx=vtx, node_tol=node_tol) 
            append!(points, [point])
            append!(weights, wlgl_f2[i]*wlgl_c2[j])
        end
    end

    if mod(ne-2,2)==1
        for i = 2:convert(Int,(ne-2-1)/2)+1
            push!(points,centerline_1[:,i])
            append!(weights, wlgl_c1[end]*wlgl_c1[i])
        end
    end
    points = hcat(points...)
    return points, weights
end 

function get_interior_points_tri_omega(p::Int; ne::Int=-1, T=Float64, node_tol=1e-14) 

    vtx = [-1 -1/sqrt(3); 1 -1/sqrt(3); 0 2/sqrt(3)]
    if ne==-1
        ne=p+1
    end
    points = []
    weights = T[]
    numlines = convert(Int, floor((ne-0)/2))
    facet1_lg, wlg_f1 = get_lg_nodes_2d(p, Array(T[1,-1/sqrt(3)]), Array(T[-1,-1/sqrt(3)]), ne=ne, half=false)
    facet2_lg, wlg_f2 = get_lg_nodes_2d(p, Array(T[1,-1/sqrt(3)]), Array(T[0,2/sqrt(3)]), ne=ne, half=false)
    a = 0/(2*p)
    centerline_1, wlg_c1 = get_lg_nodes_2d(p, Array(T[0,-1/sqrt(3)]), Array(T[0, a]), ne=ne, half=true)
    d=sqrt((1/2)^2+(sqrt(3)/2 - 1/sqrt(3))^2)
    c=a+d
    x1= 1/2 - c*cos(pi/6)
    y1= sqrt(3)/2-1/sqrt(3) - c*sin(pi/6)
    centerline_2, wlg_c2 = get_lg_nodes_2d(p, Array(T[1/2,sqrt(3)/2-1/sqrt(3)]), Array(T[x1,y1]), ne=ne, half=true)

    for i = 1:numlines
        for j = i:numlines
            point = get_intersection(facet2_lg[:,i], centerline_1[:,i], 
                            facet1_lg[:,j], centerline_2[:,j], vtx=vtx, node_tol=node_tol) 
            append!(points, [point])
            append!(weights, wlg_f2[i]*wlg_c2[j])
        end
    end

    if mod(ne,2)==1
        for i = 1:convert(Int,(ne-1)/2)
            push!(points,centerline_1[:,i])
            append!(weights, wlg_c1[end]*wlg_c1[i])
        end
    end
    points = hcat(points...)
    return points, weights
end 

function get_interior_points_tet(p::Int; ne::Int=-1, T=Float64, node_tol=1e-14) 

    vtx = [1 -sqrt(3)/3 -sqrt(6)/6;
            0 2*sqrt(3)/3 -sqrt(6)/6;
            -1 -sqrt(3)/3 -sqrt(6)/6;
            0 0 sqrt(6)/2]
    if ne==-1
        ne=p+2
    end
    points = []
    weights = []
    numlines = convert(Int, floor(p/2))
    vert_hex = [[1/2, sqrt(3)/6, -sqrt(6)/6],
                [0, 2/sqrt(3), -sqrt(6)/6],
                [0, 1/sqrt(3), sqrt(6)/6],
                [1/3, sqrt(3)/9, sqrt(6)/18],
                [0, 0, -sqrt(6)/6],
                [-1/2, sqrt(3)/6, -sqrt(6)/6],
                [-1/3, sqrt(3)/9, sqrt(6)/18],
                [0, 0, 0]]
    a =0
    line1, w1 = get_lgl_nodes_3d(p, Array(vert_hex[2]), Array(vert_hex[1]), half=true)
    line2, w2 = get_lgl_nodes_3d(p, Array(vert_hex[2]), Array(vert_hex[3].+a), half=true)
    line3, w3 = get_lgl_nodes_3d(p, Array(vert_hex[3]), Array(vert_hex[4]), half=true)
    line4, w4 = get_lgl_nodes_3d(p, Array(vert_hex[1]), Array(vert_hex[4].+a), half=true)
    line5, w5 = get_lgl_nodes_3d(p, Array(vert_hex[6]), Array(vert_hex[5]), half=true)
    line6, w6 = get_lgl_nodes_3d(p, Array(vert_hex[6]), Array(vert_hex[7].+a), half=true)
    line7, w7 = get_lgl_nodes_3d(p, Array(vert_hex[7]), Array(vert_hex[8]), half=true)
    line8, w8 = get_lgl_nodes_3d(p, Array(vert_hex[5]), Array(vert_hex[8].+a), half=true)
    line9, w9 = get_lgl_nodes_3d(p, Array(vert_hex[1]), Array(vert_hex[5]), half=true)
    line10, w10 = get_lgl_nodes_3d(p, Array(vert_hex[4]), Array(vert_hex[8]), half=true)
    line11, w11 = get_lgl_nodes_3d(p, Array(vert_hex[2]), Array(vert_hex[6]), half=true)
    line12, w12 = get_lgl_nodes_3d(p, Array(vert_hex[3]), Array(vert_hex[7]), half=true)


    n = convert(Int,(ne-mod(p,2))/2) + mod(p,2)
    # for i = 2:n
    #     line26, w26 = get_lgl_nodes_3d(p, Array(line2[i]), Array(line6[i]), half=true)
    #     line48, w48 = get_lgl_nodes_3d(p, Array(line4[i]), Array(line8[i]), half=true)
    #     for j = 2:n
    #         line2648, w2648  =  get_lgl_nodes_3d(p, Array(line26[j]), Array(line48[j]), half=true)
    #         push!(points, line2648[2:n])
    #         push!(weights, (w2648[2:n].^2 .*reverse(w2648[2:n]))[:])
    #     end
    # end
    # points = hcat(hcat(points...)...)
    # weights = vcat(weights...)'

    for i = 2:n
        line26, w26 = get_lgl_nodes_3d(p, Array(line2[i]), Array(line6[i].+a), half=true)
        line48, w48 = get_lgl_nodes_3d(p, Array(line4[i]), Array(line8[i].+a), half=true)
        for j = 2+(i-2):n
            line2648, w2648  =  get_lgl_nodes_3d(p, Array(line26[j]), Array(line48[j]), half=true)
            for k = 2+(j-2):n
                push!(points, line2648[k])
                push!(weights, (w2648[k]^2 *reverse(w2648)[k]))
            end
        end
    end
    points = hcat(points...)
    weights = hcat(weights...)
    return points, weights
end 

function get_interior_points_tet_omega(p::Int; ne::Int=-1, T=Float64, node_tol=1e-14) 

    vtx = [1 -sqrt(3)/3 -sqrt(6)/6;
            0 2*sqrt(3)/3 -sqrt(6)/6;
            -1 -sqrt(3)/3 -sqrt(6)/6;
            0 0 sqrt(6)/2]
    if ne==-1
        ne=p+1
    end
    points = []
    weights = []
    numlines = convert(Int, floor(p/2))
    vert_hex = [[1/2, sqrt(3)/6, -sqrt(6)/6],
                [0, 2/sqrt(3), -sqrt(6)/6],
                [0, 1/sqrt(3), sqrt(6)/6],
                [1/3, sqrt(3)/9, sqrt(6)/18],
                [0, 0, -sqrt(6)/6],
                [-1/2, sqrt(3)/6, -sqrt(6)/6],
                [-1/3, sqrt(3)/9, sqrt(6)/18],
                [0, 0, 0]]
    a =0
   
    line2, w2 = get_lg_nodes_3d(p, Array(vert_hex[2]), Array(vert_hex[3].+a), half=true)
    line4, w4 = get_lg_nodes_3d(p, Array(vert_hex[1]), Array(vert_hex[4].+a), half=true)
    line6, w6 = get_lg_nodes_3d(p, Array(vert_hex[6]), Array(vert_hex[7].+a), half=true)
    line8, w8 = get_lg_nodes_3d(p, Array(vert_hex[5]), Array(vert_hex[8].+a), half=true)

    # n = convert(Int,(ne-mod(p,2))/2) + mod(p,2)
    n = convert(Int,(ne-mod(p+1,2))/2) + mod(p+1,2)
    for i = 1:n
        line26, w26 = get_lg_nodes_3d(p, Array(line2[i]), Array(line6[i].+a), half=true)
        line48, w48 = get_lg_nodes_3d(p, Array(line4[i]), Array(line8[i].+a), half=true)
        for j = 2+(i-2):n
            line2648, w2648  =  get_lg_nodes_3d(p, Array(line26[j]), Array(line48[j]), half=true)
            for k = 2+(j-2):n
                push!(points, line2648[k])
                push!(weights, (w2648[k]^2 *reverse(w2648)[k]))
            end
        end
    end
    points = hcat(points...)
    weights = hcat(weights...)
    return points, weights
end 

function get_interior_points_parabola_tri(p::Int; node_tol=1e-14) 

    points = []
    numparab = convert(Int, floor(p/2))
    parab = SummationByParts.get_parabola_tri(p)
    
    # rotate parabolas and compute intersections using SymPy
    x,y=symbols("x,y",real=true)
    t = 2*pi/3

    for i = 1:numparab
        k1 = parab[1,i]
        a1 = parab[2,i]
        eq1 = y - k1 - a1*x^2
        for j = 1:numparab
            k2 = parab[1,j]
            a2 = parab[2,j]
            eq2 = x*sin(t)+y*cos(t)-k2 - a2*(x*cos(t)-y*sin(t))^2
            sol = solve([eq1,eq2])
            sols=[[convert(Float64,sol[1][x]),convert(Float64,sol[1][y])],
                  [convert(Float64,sol[2][x]),convert(Float64,sol[2][y])]]
            for k=1:2
                if ((sols[k][1] - 1) < node_tol && sols[k][1] + 1 > -node_tol && 
                    (sols[k][2] + 1/sqrt(3)) > -node_tol &&
                    (sols[k][2]+sqrt(3)*sols[k][1] - sqrt(3) + 1/sqrt(3)) < node_tol && 
                    (sols[k][2]-sqrt(3)*sols[k][1] - sqrt(3) + 1/sqrt(3)) < node_tol)

                    push!(points,sols[k])
                    #rotate point by 120 and 240 degrees add them to the list
                    push!(points,[sols[k][1]*cos(2*pi/3)-sols[k][2]*sin(2*pi/3), 
                                  sols[k][1]*sin(2*pi/3)+sols[k][2]*cos(2*pi/3)])
                    push!(points,[sols[k][1]*cos(4*pi/3)-sols[k][2]*sin(4*pi/3), 
                                  sols[k][1]*sin(4*pi/3)+sols[k][2]*cos(4*pi/3)])
                end
            end
        end
    end

    if mod(p,2)==1
        centerline1 = get_lgl_nodes_2d(p, Array(Float64[0,-1/sqrt(3)]), Array(Float64[0, 1/p]), half=false)
        for i = 2:convert(Int,(p-1)/2)+1
            push!(points,centerline1[:,i])
            #rotate point by 120 and 240 degrees add them to the list
            push!(points,[centerline1[1,i]*cos(2*pi/3)-centerline1[2,i]*sin(2*pi/3), 
                          centerline1[1,i]*sin(2*pi/3)+centerline1[2,i]*cos(2*pi/3)])
            push!(points,[centerline1[1,i]*cos(4*pi/3)-centerline1[2,i]*sin(4*pi/3), 
                          centerline1[1,i]*sin(4*pi/3)+centerline1[2,i]*cos(4*pi/3)])
        end
        push!(points,[0.0,0.0])
    end
    points = hcat(points...)
    return points
end 

# function get_intersection_parabola_tri(p1::Array{T,1}, p2::Array{T,1}; node_tol::Float64=1e-14) where {T}

#     # vtx = [-1 -1/sqrt(3); 1 -1/sqrt(3); 0 2/sqrt(3)]

#     a = p1[1] - p2[1]
#     b = p1[2] - p2[2]
#     c = p1[3] - p2[3]

#     point = []
#     if (a==0 && b==0 && c==0)
#         point = []
#     elseif (a!=0 && (b^2 - 4*a*c)>node_tol)
#         x1 = (-b + sqrt(b^2 - 4*a*c))/(2*a)
#         x2 = (-b - sqrt(b^2 - 4*a*c))/(2*a)
#         y1 = p1[1]*x1^2 + p1[2]*x1 + p1[3]
#         y2 = p1[1]*x2^2 + p1[2]*x2 + p1[3]
#         x = [x1,x2]
#         y = [y1,y2]
#         for i = 1:2
#             if ((x[i] - 1) < node_tol && x[i] + 1 > -node_tol && (y1 + 1/sqrt(3)) > -node_tol &&
#                 (y[i]+sqrt(3)*x[i] - sqrt(3) + 1/sqrt(3)) < node_tol && 
#                 (y[i]-sqrt(3)*x[i] - sqrt(3) + 1/sqrt(3)) < node_tol)
#                 push!(point,[x[i],y[i]])
#             end
#         end
#     end
        
#     return point
# end 

function get_lgl_nodes_2d(p::Int, node1::Array{T,1}, node2::Array{T,1}; half::Bool=false, ne::Int=-1, node_tol=1e-14) where {T}

    if ne==-1
        ne=p+2
    end
    qf = 2*ne-3

    d = sqrt((node2[2]-node1[2])^2 + (node2[1]-node1[1])^2)    
    x1 = node1[1]; y1 = node1[2]; 
    x2 = node2[1]; y2 = node2[2];
    m = 0
    if abs(x2 - x1) > node_tol
        m = (y2 - y1)/(x2 - x1)
    end

    if half
        if abs(x2-x1) < node_tol
            if y2 - y1 > node_tol
                y2 = y2 + (y2-y1)
            else
                y2 = y2 - (y1 - y2)
            end 
        elseif x2 - x1 > node_tol
            x2 = x2 + sqrt(4*d^2/(m^2+1)) 
            y2 = m*(x2-x1) + y1
        elseif x1 - x2 > node_tol
            x2 = x2 - sqrt(4*d^2/(m^2+1)) 
            y2 = m*(x2-x1) + y1
        end
        d = 2*d
    end

    cub_lgl, vtx_lgl = SummationByParts.Cubature.quadrature(qf, internal=false)
    vtx_lgl[1]=0; vtx_lgl[2]=d; 
    perm = sortperm(vec(SymCubatures.calcnodes(cub_lgl, vtx_lgl)))
    xlgl = vec(SymCubatures.calcnodes(cub_lgl, vtx_lgl))[perm]
    wlgl = (sqrt(3)/3) .*(vec(SymCubatures.calcweights(cub_lgl))[perm])

    xlgl_2d = zeros(T, (2,length(xlgl)))
    for i = 1:length(xlgl)
        if x2 - x1 > node_tol
            xlgl_2d[1,i] = x1 + sqrt(xlgl[i]^2/(m^2+1))
            xlgl_2d[2,i] = m*(xlgl_2d[1,i] - x1) + y1
        elseif x1 - x2 > node_tol
            xlgl_2d[1,i] = x1 - sqrt(xlgl[i]^2/(m^2+1))
            xlgl_2d[2,i] = m*(xlgl_2d[1,i] - x1) + y1
        elseif abs(x2-x1) < node_tol
            xlgl_2d[1,i] = x1
            xlgl_2d[2,i] = y1+xlgl[i]
        end
    end

    if half 
        xlgl_2d = xlgl_2d[:,1:convert(Int, ceil(size(xlgl_2d,2)/2))]
        wlgl = wlgl[:,1:convert(Int, ceil(size(wlgl,2)/2))]
    end
    return xlgl_2d, wlgl
end

function get_lg_nodes_2d(p::Int, node1::Array{T,1}, node2::Array{T,1}; half::Bool=false, ne::Int=-1, node_tol=1e-14) where {T}

    if ne==-1
        ne=p+1
    end
    qf = 2*ne-1

    d = sqrt((node2[2]-node1[2])^2 + (node2[1]-node1[1])^2)    
    x1 = node1[1]; y1 = node1[2]; 
    x2 = node2[1]; y2 = node2[2];
    m = 0
    if abs(x2 - x1) > node_tol
        m = (y2 - y1)/(x2 - x1)
    end

    if half
        if abs(x2-x1) < node_tol
            if y2 - y1 > node_tol
                y2 = y2 + (y2-y1)
            else
                y2 = y2 - (y1 - y2)
            end 
        elseif x2 - x1 > node_tol
            x2 = x2 + sqrt(4*d^2/(m^2+1)) 
            y2 = m*(x2-x1) + y1
        elseif x1 - x2 > node_tol
            x2 = x2 - sqrt(4*d^2/(m^2+1)) 
            y2 = m*(x2-x1) + y1
        end
        d = 2*d
    end

    cub_lg, vtx_lg = SummationByParts.Cubature.quadrature(qf, internal=true)
    vtx_lg[1]=0; vtx_lg[2]=d; 
    perm = sortperm(vec(SymCubatures.calcnodes(cub_lg, vtx_lg)))
    xlg = vec(SymCubatures.calcnodes(cub_lg, vtx_lg))[perm]
    wlg = (sqrt(3)/3) .*(vec(SymCubatures.calcweights(cub_lg))[perm])

    xlg_2d = zeros(T, (2,length(xlg)))
    for i = 1:length(xlg)
        if x2 - x1 > node_tol
            xlg_2d[1,i] = x1 + sqrt(xlg[i]^2/(m^2+1))
            xlg_2d[2,i] = m*(xlg_2d[1,i] - x1) + y1
        elseif x1 - x2 > node_tol
            xlg_2d[1,i] = x1 - sqrt(xlg[i]^2/(m^2+1))
            xlg_2d[2,i] = m*(xlg_2d[1,i] - x1) + y1
        elseif abs(x2-x1) < node_tol
            xlg_2d[1,i] = x1
            xlg_2d[2,i] = y1+xlg[i]
        end
    end

    if half 
        xlg_2d = xlg_2d[:,1:convert(Int, ceil(size(xlg_2d,2)/2))]
        wlg = wlg[:,1:convert(Int, ceil(size(wlg,2)/2))]
    end
    return xlg_2d, wlg
end

function get_lgl_nodes_3d(p::Int, node1::Array{T,1}, node2::Array{T,1}; ne::Int=-1, half::Bool=false) where {T}

    if ne==-1
        ne=p+2
    end
    qf = 2*ne-3

    cub_lgl, vtx_lgl = SummationByParts.Cubature.quadrature(qf, internal=false)
    vtx_lgl[1]=0; vtx_lgl[2]=1; 
    perm = sortperm(vec(SymCubatures.calcnodes(cub_lgl, vtx_lgl)))
    xlgl = vec(SymCubatures.calcnodes(cub_lgl, vtx_lgl))[perm]
    # scale the weights by the volume divided by four (since there are four hexahedra)
    # also the weight is based on [-1,1], so should be divided by 2
    wlgl = (sqrt(2)/6*1/2) .*(vec(SymCubatures.calcweights(cub_lgl))[perm])

    # get the direction vector
    m = node2 .- node1
    xlgl_3d = []
    if half
        for i=1:convert(Int,ceil(length(xlgl)/2))
            push!(xlgl_3d, node1 .+ (2*xlgl[i]).*m)
        end
    else
        for i=1:length(xlgl)
            push!(xlgl_3d, node1 .+ xlgl[i].*m)
        end
    end

    if half 
        wlgl = wlgl[:,1:convert(Int, ceil(size(wlgl,2)/2))]
    end

    return xlgl_3d, wlgl
end

function get_lg_nodes_3d(p::Int, node1::Array{T,1}, node2::Array{T,1}; ne::Int=-1, half::Bool=false) where {T}

    if ne==-1
        ne=p+1
    end
    qf = 2*ne-1

    cub_lg, vtx_lg = SummationByParts.Cubature.quadrature(qf, internal=true)
    vtx_lg[1]=0; vtx_lg[2]=1; 
    perm = sortperm(vec(SymCubatures.calcnodes(cub_lg, vtx_lg)))
    xlg = vec(SymCubatures.calcnodes(cub_lg, vtx_lg))[perm]
    # scale the weights by the volume divided by four (since there are four hexahedra)
    # also the weight is based on [-1,1], so should be divided by 2
    wlg = (sqrt(2)/6*1/2) .*(vec(SymCubatures.calcweights(cub_lg))[perm])

    # get the direction vector
    m = node2 .- node1
    xlg_3d = []
    if half
        for i=1:convert(Int,ceil(length(xlg)/2))
            push!(xlg_3d, node1 .+ (2*xlg[i]).*m)
        end
    else
        for i=1:length(xlg)
            push!(xlg_3d, node1 .+ xlg[i].*m)
        end
    end

    if half 
        wlg = wlg[:,1:convert(Int, ceil(size(wlg,2)/2))]
    end

    return xlg_3d, wlg
end

# function get_parabola_tri(p::Int; T=Float64)

#     n = p+2
#     nnode = convert(Int, floor(p/2))
#     numparab = 3*nnode
#     parab = zeros(T, (numparab,2,3))
#     parab_coef = zeros(T, (3,numparab))

#     facet1_lgl = get_lgl_nodes_2d(p, Array(T[-1,-1/sqrt(3)]), Array(T[1,-1/sqrt(3)]), half=false)
#     facet2_lgl = get_lgl_nodes_2d(p, Array(T[1,-1/sqrt(3)]), Array(T[0,2/sqrt(3)]), half=false)
#     facet3_lgl = get_lgl_nodes_2d(p, Array(T[0,2/sqrt(3)]), Array(T[-1,-1/sqrt(3)]), half=false)

#     centerline_1 = get_lgl_nodes_2d(p, Array(T[0,-1/sqrt(3)]), Array(T[0, 0]), half=true)
#     centerline_2 = get_lgl_nodes_2d(p, Array(T[1/2,sqrt(3)/2-1/sqrt(3)]), Array(T[0, 0]), half=true)
#     centerline_3 = get_lgl_nodes_2d(p, Array(T[-1/2,sqrt(3)/2-1/sqrt(3)]), Array(T[0, 0]), half=true)

#     for i = 1:nnode
#         parab[i,:,1] = facet1_lgl[:,i+1]
#         parab[i,:,2] = centerline_3[:,i+1]
#         parab[i,:,3] = facet2_lgl[:,n-i]

#         parab[nnode+i,:,1] = facet1_lgl[:,n-i]
#         parab[nnode+i,:,2] = centerline_2[:,i+1]
#         parab[nnode+i,:,3] = facet3_lgl[:,i+1]

#         parab[2*nnode+i,:,1] = facet2_lgl[:,i+1]
#         parab[2*nnode+i,:,2] = centerline_1[:,i+1]
#         parab[2*nnode+i,:,3] = facet3_lgl[:,n-i]
#     end

#     for i=1:numparab
#         x = parab[i,1,:]
#         y = parab[i,2,:]
#         v = hcat([x.^2, x, ones(3,1)]...)
#         parab_coef[:,i] = v\y
#     end

#     return parab, parab_coef
# end

function get_parabola_tri(p::Int; T=Float64)

    numparab = convert(Int, floor(p/2))
    parab = zeros(T, (2,numparab)) # contains [k;a] in y-k = a*x^2

    facet2_lgl = get_lgl_nodes_2d(p, Array(T[1,-1/sqrt(3)]), Array(T[0,2/sqrt(3)]), half=false)
    centerline_1 = get_lgl_nodes_2d(p, Array(T[0,-1/sqrt(3)]), Array(T[0, 1/p]), half=true)

    for i = 1:numparab
        parab[1,i] = centerline_1[2,i+1]
        parab[2,i] = (facet2_lgl[2,i+1] - parab[1,i])/(facet2_lgl[1,i+1]^2)
    end

    return parab
end



function init_tri_staggered_nodes(p::Int,q::Int; ne::Int=-1, T=Float64)

    # use equilateral triangle
    node_tol = 1e-14
    vtx = T[-1 -1/sqrt(3); 1 -1/sqrt(3); 0 2/sqrt(3)] 
    if ne==-1
        ne=p+2
    end
    qf = 2*ne-3
    cub_lgl, _ = SummationByParts.Cubature.quadrature(qf, internal=false)
    xedge = cub_lgl.params
    wlgl = cub_lgl.weights

    # get interior nodes
    cub_omega, _ = SummationByParts.Cubature.getTriCubatureOmega(q-1, T)
    param_S21 = cub_omega.params[1:cub_omega.numS21]
    param_S111 = cub_omega.params[cub_omega.numS21+1:cub_omega.numS21+2*cub_omega.numS111]
    weight_S21 = cub_omega.weights[1:cub_omega.numS21]
    weight_S111 = cub_omega.weights[cub_omega.numS21+1:cub_omega.numS21+cub_omega.numS111]

    add_S21 = []
    add_S21_w = []
    add_S111 = []
    add_S111_w = []
    add_centroid = false

    vertices = true
    midedges = convert(Bool,mod(ne-2,2))
    numedge = convert(Int,(ne-2-mod(ne-2,2))/2)
    numS21 = cub_omega.numS21 + length(add_S21)
    numS111 = cub_omega.numS111 + convert(Int, length(add_S111)/2)
    centroid = (cub_omega.centroid || add_centroid)


    cub = SymCubatures.TriSymCub{T}(vertices=vertices,
                                    midedges=midedges,
                                    numS21=numS21,
                                    numedge=numedge,
                                    numS111=numS111,
                                    centroid = centroid)
    
    # set the edge parameters 
    init_tri_cub(cub, xedge)

    # get barycentric coordinates of the interior nodes
    params=[]
    weight_vert=[]
    weight_mid=[]
    weight_cent=[]

    params = convert.(T,collect(Iterators.flatten([param_S21,add_S21,cub_lgl.params,param_S111,add_S111])))
    J = 1 #1/2*sqrt(3)/3
    if vertices
        weight_vert=wlgl[1]^2
    end
    if midedges
        weight_mid=wlgl[end]*wlgl[end-convert(Int,centroid)]
    end
    if centroid
        weight_cent=cub_omega.weights[end]
    end
    weight_edge= wlgl[2:end-convert(Int,cub_lgl.centroid)].*wlgl[1]
    weights = convert.(T,collect(Iterators.flatten([J.*weight_vert,J.*weight_mid, weight_S21, add_S21_w,
                                                    weight_S111, add_S111_w, J.*weight_edge, weight_cent])))

    if params != []
        SymCubatures.setparams!(cub, params)
        SymCubatures.setweights!(cub, weights)
    end
    # xy = SymCubatures.calcnodes(cub,vtx)
    # SummationByParts.plot_tri_nodes(x=xy, vtx=vtx, q=q, n=cub.numnodes,write_title=true,label_nodes=false)

    xinit = convert.(T,collect(Iterators.flatten([cub.params,cub.weights])))
    # mask = SymCubatures.getInternalParamMask(cub)
    mask = collect(Iterators.flatten([length(param_S21)+1:cub.numS21, 
                                      cub.numS21+cub.numedge+length(param_S111)+1: cub.numS21+cub.numedge+2*cub.numS111]))
    append!(mask, (cub.numparams+1):(cub.numparams+cub.numweights)) 
    convert(Array{Int},mask)
    Cubature.solvecubature!(cub, q, mask, tol=1e-14, hist=true, verbose=true, xinit=xinit, delta1=1e-3, delta2=1e-4)
    println("\n", cub.params,"\n")
    println(cub.weights,"\n")
    return cub
end

function init_tet_staggered_nodes(p::Int,q::Int; ne::Int=-1, T=Float64)

    # use equilateral tet
    vtx = T[1 -sqrt(3)/3 -sqrt(6)/6;
            0 2*sqrt(3)/3 -sqrt(6)/6;
            -1 -sqrt(3)/3 -sqrt(6)/6;
            0 0 sqrt(6)/2]
            
    if ne==-1
        ne=p+2
    end
    qf = 2*ne-3
    cub_lgl, _ = SummationByParts.Cubature.quadrature(qf, internal=false)
    # cub_tri,_ = SummationByParts.Cubature.getTriCubatureDiagE(q, T, vertices=true)
    # cub_tri,_ = SummationByParts.Cubature.getTriCubatureForTetFaceDiagE(q, faceopertype=:Omega)
    cub_tri,_ = SummationByParts.Cubature.getTriCubatureForTetFaceDiagE(q, faceopertype=:DiagE)
    xedge = cub_tri.params
    wlgl = cub_lgl.weights

    # get interior nodes
    cub_omega, _ = SummationByParts.Cubature.getTetCubatureOmega(q-1, T)
    param_S31 = cub_omega.params[1:cub_omega.numS31]
    weight_S31 = cub_omega.weights[1:cub_omega.numS31]
    nn=cub_omega.numS31
    nw=cub_omega.numS31
    param_S22 = cub_omega.params[nn+1:nn+cub_omega.numS22]
    weight_S22 = cub_omega.weights[nw+1:nw+cub_omega.numS22]
    nn+=cub_omega.numS22
    nw+=cub_omega.numS22
    param_S211 = cub_omega.params[nn+1:nn+2*cub_omega.numS211]
    weight_S211 = cub_omega.weights[nw+1:nw+cub_omega.numS211]
    nn+=2*cub_omega.numS211
    nw+=cub_omega.numS211
    param_S1111 = cub_omega.params[nn+1:nn+3*cub_omega.numS1111]
    weight_S1111 = cub_omega.weights[nw+1:nw+cub_omega.numS1111]

    add_S31 = []#[0.2128785129992084, 0.9673676540009212, 0.5162909518831329, 0.8593193682825483]
    add_S31_w = []#[0.00012587350281596975, 0.012244256930795668, 0.010856385828906114, 0.03437297074454886,]
    add_S22 = []#[0.8923932601020271, 0.733783076383442,]
    add_S22_w = []#[0.019933209306358283, 0.0004056776588533049]
    add_S211 = []#[0.12845378081353867, 1.2838828878296134, 0.4076479469911398, 0.07616274419943524, 0.7628782111929711, 0.1180398395722123]
    add_S211_w = []#[0.01563787368302707,0.018637092172174203, 0.0025937227508652174]
    add_S1111 = []
    add_S1111_w = []
    add_centroid = true

    numedge = cub_tri.numedge #convert(Int,(p-mod(p,2))/2)
    vertices = cub_tri.vertices #true
    midedges = cub_tri.midedges #convert(Bool,mod(p,2))
    facecentroid = cub_tri.centroid #convert(Bool, mod(p,2))
    numfaceS21 = cub_tri.numS21 #convert(Int,(1+mod(p,2))*numedge)
    numfaceS111 = cub_tri.numS111 #convert(Int,1/2*(numedge^2 - numedge))
    numS31 = cub_omega.numS31 + length(add_S31)
    numS22 = cub_omega.numS22 + length(add_S22)
    numS211= cub_omega.numS211 + convert(Int, length(add_S211)/2)
    numS1111 = cub_omega.numS1111 + convert(Int, length(add_S1111)/3)
    centroid = (cub_omega.centroid || add_centroid)

    cub = SymCubatures.TetSymCub{T}(vertices=vertices,
                                    numS31=numS31,
                                    midedges=midedges,
                                    numS22=numS22,
                                    numfaceS21=numfaceS21,
                                    numedge=numedge,
                                    numS211=numS211,
                                    numfaceS111=numfaceS111,
                                    facecentroid=facecentroid,
                                    numS1111=numS1111,
                                    centroid=centroid)
    
    # set the edge parameters 
    init_tet_cub(cub, xedge)
    faceS21_params = cub_tri.params[1 : cub_tri.numS21]
    faceS21_w = wlgl[1].*cub_tri.weights[convert(Int,cub.vertices)+convert(Int,cub.midedges)+1 : convert(Int,cub.vertices)+convert(Int,cub.midedges)+cub_tri.numS21]
    edge_params = cub_tri.params[cub_tri.numS21+1 : cub_tri.numS21+cub_tri.numedge]
    edge_w = wlgl[1].*cub_tri.weights[convert(Int,cub.vertices)+convert(Int,cub.midedges)+cub.numfaceS21+1 : convert(Int,cub.vertices)+convert(Int,cub.midedges)+cub.numfaceS21+cub.numedge]
    faceS111_params = cub_tri.params[cub_tri.numS21+cub_tri.numedge+1 : cub_tri.numS21+cub_tri.numedge+2*cub_tri.numS111]
    faceS111_w =  wlgl[1].*cub_tri.weights[convert(Int,cub.vertices)+convert(Int,cub.midedges)+cub.numfaceS21+cub.numedge+1 : convert(Int,cub.vertices)+convert(Int,cub.midedges)+cub.numfaceS21+cub.numedge+cub.numfaceS111]
    facecentroid_w = [] 
    if facecentroid #mod(p,2)==1
        facecentroid_w=wlgl[1].*cub_tri.weights[end]
    end
    # get barycentric coordinates of the interior nodes
    params=[]
    weight_vert=[]
    weight_mid=[]
    weight_cent=[]

    params = convert.(T,collect(Iterators.flatten([param_S31,add_S31,param_S22,add_S22,faceS21_params,
                                                   edge_params, 
                                                   param_S211,add_S211,
                                                   faceS111_params, param_S1111,add_S1111])))
    # params = convert.(T,collect(Iterators.flatten([add_S31,add_S22,faceS21_params,
    #                                                edge_params, 
    #                                                add_S211,
    #                                                faceS111_params, param_S1111,add_S1111])))
    J = 1
    if vertices
        weight_vert=[0.00012587350281596975] #wlgl[1]^3
    end
    if midedges
        weight_mid= [0.0004056776588533049] #wlgl[end]^2*wlgl[end-convert(Int,centroid)]
    end
    if centroid
        weight_cent=cub_omega.weights[end]
    end
    weights = convert.(T,collect(Iterators.flatten([J.*weight_vert, weight_S31, add_S31_w, J.*weight_mid, 
                                                    weight_S22, add_S22_w, faceS21_w, edge_w, weight_S211, add_S211_w,
                                                    faceS111_w, facecentroid_w,
                                                    weight_S1111, add_S1111_w, weight_cent])))

    # weights = convert.(T,collect(Iterators.flatten([J.*weight_vert,  add_S31_w, J.*weight_mid, 
    #                                              add_S22_w, faceS21_w, edge_w,  add_S211_w,
    #                                             faceS111_w, facecentroid_w,
    #                                             weight_S1111, add_S1111_w, weight_cent])))
    if params != []
        SymCubatures.setparams!(cub, params)
        SymCubatures.setweights!(cub, weights)
    end
    # xy = SymCubatures.calcnodes(cub,vtx)
    # SummationByParts.plot_tri_nodes(x=xy, vtx=vtx, q=q, n=cub.numnodes,write_title=true,label_nodes=false)

    xinit = convert.(T,collect(Iterators.flatten([cub.params,cub.weights])))
    # mask = SymCubatures.getInternalParamMask(cub)
    mask = collect(Iterators.flatten([length(param_S31)+1 : cub.numS31, 
                                      cub.numS31+length(param_S22)+1 : cub.numS31+cub.numS22,
                                      cub.numS31+cub.numS22+cub.numfaceS21+cub.numedge+length(param_S211)+1 : cub.numS31+cub.numS22+cub.numfaceS21+cub.numedge+2*cub.numS211,
                                      cub.numS31+cub.numS22+cub.numfaceS21+cub.numedge+2*cub.numS211+2*cub.numfaceS111+length(param_S1111)+1 : cub.numS31+cub.numS22+cub.numfaceS21+cub.numedge+2*cub.numS211+2*cub.numfaceS111+3*cub.numS1111]))
    append!(mask, (cub.numparams+1):(cub.numparams+cub.numweights)) 
    convert(Array{Int},mask)
    Cubature.solvecubature!(cub, q, mask, tol=1e-14, hist=true, verbose=true, xinit=xinit, delta1=1e-5, delta2=1e-3)
    println("\n", cub.params,"\n")
    println(cub.weights,"\n")
    return cub
end

# function eliminate_nodes(cub::TriSymCub{T}, p::Int, q::Int) where {T}
#     vtx = T[-1 -1; 1 -1; -1 1]

#     n = cub.numweights
#     j=1
#     # j = convert(Int, floor(n/2))
#     res_min = -1.0
#     # while (j <= convert(Int, ceil(n/3)) && cub.numnodes>=1)
#     # while (j >0 && cub.numnodes>=1)
#     while (j <= n && cub.numnodes>=1)
#         idx = []
#         param_idx = []
#         weight_idx = []
#         sym_group = String[]
#         iig = 0
#         iig_idx = []
#         ptr = 0
#         paramptr = 0
#         weightptr = 0
#         if cub.vertices
#             append!(idx, ptr+1)
#             ptr += 3
#             append!(weight_idx, weightptr+1)
#             weightptr +=1
#             push!(sym_group, "vertices")
#         end
#         if cub.midedges
#             append!(idx, ptr+1)
#             ptr += 3
#             append!(weight_idx, weightptr+1)
#             weightptr +=1
#             push!(sym_group, "midedges")
#         end
#         for i=1:cub.numS21
#             append!(idx, ptr+1)
#             ptr += 3
#             append!(param_idx, [paramptr+1])
#             paramptr += 1
#             append!(weight_idx, weightptr+1)
#             weightptr +=1
#             push!(sym_group, "S21")
#             iig += 1
#             push!(iig_idx, iig)
#         end
#         for i = 1:cub.numedge
#             append!(idx, ptr+1)
#             ptr += 6
#             append!(param_idx, [paramptr+1])
#             paramptr += 1
#             append!(weight_idx, weightptr+1)
#             weightptr +=1
#             push!(sym_group, "edge")
#             iig += 1
#             push!(iig_idx, iig)
#         end
#         for i = 1:cub.numS111
#             append!(idx, ptr+1)
#             ptr += 6
#             append!(param_idx, [paramptr+1:paramptr+2])
#             paramptr += 2
#             append!(weight_idx, weightptr+1)
#             weightptr +=1
#             push!(sym_group, "S111", "S111")
#             iig += 1
#             push!(iig_idx, iig, iig)
#         end
#         if cub.centroid
#             append!(idx, ptr+1)
#             ptr += 1
#             append!(weight_idx, weightptr+1)
#             weightptr +=1
#             push!(sym_group, "centroid")
#         end
#         xy = SymCubatures.calcnodes(cub, vtx)
#         xy_sym = xy[:,idx]

#         # compute the Vandermonde matrix
#         V, _,_= OrthoPoly.vandermonde(p, xy_sym[1,:],xy_sym[2,:],compute_grad=false)
#         # V, _,_,_ = OrthoPoly.vandermonde_monomial(p, xy_sym[1,:],xy_sym[2,:],compute_grad=false,compute_integ=false)
#         # V, _,_,_= OrthoPoly.vandermonde_arnoldi(p, xy_sym[1,:],xy_sym[2,:],compute_grad=false)
#         V2 = V.^2*ones(size(V,2),1)
#         w = SymCubatures.calcweights(cub)
#         s = sortperm(vec(diagm(w[idx])*V2),rev=false)
#         # s = sortperm(vec(diagm(w[idx])*V2),rev=true)
#         # s = sortperm(vec(V2),rev=false)
#         # s = sortperm(vec(V2),rev=true)

#         # ss = sym_group[s]
#         # nk = convert(Int,floor(n)) 
#         # ssnk = ss[1:nk]

#         # ind_S111 = findall(x -> x == "S111", ssnk)
#         # if length(ind_S111) !=0
#         #     s[1:length(ind_S111)]=s[ind_S111]
#         # end
#         # ind_S21 = findall(x -> x == "S21", ssnk)
#         # if length(ind_S21) !=0
#         #     s[length(ind_S111)+1:length(ind_S111)+length(ind_S21)] = s[ind_S21]
#         # end
#         # ind_cent= findall(x -> x == "centroid", ssnk)
#         # if length(ind_cent) !=0
#         #     s[length(ind_S111)+length(ind_S111)+1: length(ind_S111)+length(ind_S111)+length(ind_S21)]=s[ind_cent]
#         # end

#         # ind_cent= findall(x -> x == "centroid", ssnk)
#         # if length(ind_cent) !=0
#         #     s[1]=s[ind_cent]
#         # end
#         # ind_S111 = findall(x -> x == "S111", ssnk)
#         # if length(ind_S111) !=0
#         #     s[length(ind_cent)+1:length(ind_cent)+length(ind_S111)]=s[ind_S111]
#         # end
#         # ind_S21 = findall(x -> x == "S21", ssnk)
#         # if length(ind_S21) !=0
#         #     s[length(ind_cent)+length(ind_S111)+1:length(ind_cent)+length(ind_S111)+length(ind_S21)] = s[ind_S21]
#         # end

#         ii = s[j]
#         sg=splice!(sym_group, ii)
#         params = copy(cub.params)
#         if sg ∉ ["vertices", "midedges", "centroid"]
#             splice!(params, param_idx[iig_idx[ii]]) 
#         end
#         weights = copy(cub.weights)
#         splice!(weights, weight_idx[ii]) 

#         cub2 = SymCubatures.TriSymCub{T}(vertices=convert(Bool,count(x -> x == "vertices", sym_group)),
#                                         midedges=convert(Bool,count(x -> x == "midedges", sym_group)),
#                                         numS21=count(x -> x == "S21", sym_group),
#                                         numedge=count(x -> x == "edge", sym_group),
#                                         numS111=convert(Int,floor(count(x -> x == "S111", sym_group)/2)),
#                                         centroid =convert(Bool,count(x -> x == "centroid", sym_group)))
#         SymCubatures.setparams!(cub2, params)
#         SymCubatures.setweights!(cub2, weights)
        
#         if cub2.numnodes>=1
#             xinit = convert.(T,collect(Iterators.flatten([cub2.params,cub2.weights])))
#             mask = Int[]
#             append!(mask, 1:cub2.numparams+cub2.numweights) 
#             # res = Cubature.solvecubature!(cub2, q, mask, tol=5e-14, hist=true, verbose=true,  xinit=xinit, delta1=1e-4, delta2=1e-4)
#             res = Cubature.solvecubaturelma!(cub2, q, mask, tol=5e-14, hist=true, verbose=true,  maxiter=30, xinit=xinit, delta1=1e-1, delta2=1e-1)
#             println("\n", cub2.params,"\n")
#             println(cub2.weights,"\n")
#             println("j = ", j, ":  n = ", n, ":  numnodes = ", cub.numnodes)

#             if res < 1e-1 && res > 5e-3
#                 xinit = convert.(T,collect(Iterators.flatten([cub2.params,cub2.weights])))
#                 # xinit = (1.0 .-rand(size(xinit))*0.05).*xinit
#                 # res = Cubature.solvecubature!(cub2, q, mask, tol=5e-14, hist=true, verbose=true,  xinit=xinit, delta1=1e-4, delta2=1e-4)
#                 res = Cubature.solvecubaturelma!(cub2, q, mask, tol=5e-14, hist=true, verbose=true,  maxiter=50, xinit=xinit, delta1=1e-4, delta2=1e-4)
#             end
#             if res <= 1e-2 && res > 5e-14
#                 xinit = convert.(T,collect(Iterators.flatten([cub2.params,cub2.weights])))
#                 res = Cubature.solvecubaturelma!(cub2, q, mask, tol=5e-11, hist=true, verbose=true,  maxiter=100, xinit=xinit, delta1=1e-4, delta2=1e-4)
#             end

#             j+=1
#             # j-=1
#             if res < 1e-10
#                 n = cub2.numweights 
#                 j = 1
#                 # j=convert(Int, floor(n/2))
#                 cub = cub2
#                 res_min = res
#             end
#         else
#             break 
#         end
#     end
#     println("--------------------------")
#     if res_min==-1.0
#         print("No new solution was found.")
#         println("\n", cub.params,"\n")
#         println(cub.weights,"\n")
#     else
#         println("res norm = ", res_min)
#         println("\n", cub.params,"\n")
#         println(cub.weights,"\n")
#     end
#     return cub, res_min
# end 

function eliminate_nodes(cub::TriSymCub{T}, p::Int, q::Int) where {T}
    vtx = T[-1 -1; 1 -1; -1 1]
 
    res_min = -1.0
    nu = 1e1
    for k = 1:4
        n = cub.numweights
        j=n
        while (j >= 1 && cub.numnodes>=1)
            idx = []
            param_idx = []
            weight_idx = []
            sym_group = String[]
            iig = 0
            iig_idx = []
            ptr = 0
            paramptr = 0
            weightptr = 0
            if cub.vertices
                append!(idx, ptr+1)
                ptr += 3
                append!(weight_idx, weightptr+1)
                weightptr +=1
                push!(sym_group, "vertices")
            end
            if cub.midedges
                append!(idx, ptr+1)
                ptr += 3
                append!(weight_idx, weightptr+1)
                weightptr +=1
                push!(sym_group, "midedges")
            end
            for i=1:cub.numS21
                append!(idx, ptr+1)
                ptr += 3
                append!(param_idx, [paramptr+1])
                paramptr += 1
                append!(weight_idx, weightptr+1)
                weightptr +=1
                push!(sym_group, "S21")
                iig += 1
                push!(iig_idx, iig)
            end
            for i = 1:cub.numedge
                append!(idx, ptr+1)
                ptr += 6
                append!(param_idx, [paramptr+1])
                paramptr += 1
                append!(weight_idx, weightptr+1)
                weightptr +=1
                push!(sym_group, "edge")
                iig += 1
                push!(iig_idx, iig)
            end
            for i = 1:cub.numS111
                append!(idx, ptr+1)
                ptr += 6
                append!(param_idx, [paramptr+1:paramptr+2])
                paramptr += 2
                append!(weight_idx, weightptr+1)
                weightptr +=1
                push!(sym_group, "S111")
                iig += 1
                push!(iig_idx, iig, iig)
            end
            if cub.centroid
                append!(idx, ptr+1)
                ptr += 1
                append!(weight_idx, weightptr+1)
                weightptr +=1
                push!(sym_group, "centroid")
            end
            
            w = SymCubatures.calcweights(cub)
            # xy = SymCubatures.calcnodes(cub, vtx)
            # ssf = []
            # for m=1:cub.numnodes
            #     xy_m = xy[:, setdiff(1:end, m)]
            #     mindist = convert(T, Inf)
            #     if cub.numnodes > 1
            #         mindist=calc_mindist_to_facet(xy_m)
            #     end 
            #     # ww = w[idx][setdiff(1:end,m),:] 
            #     push!(ssf, mindist)
            # end
            # ss = ssf[idx]
            # rev_sort = false
            
            if (k<=2 && cub.numnodes>1)
                nmin,sym_min = minnodes_bound_omega_tri(q)
                sym_cnt = 0
                ii=j -convert(Int,cub.centroid)
                sym_cnt=cub.numweights
                # if cub.centroid && sym_min[1]==0 && sym_cnt-convert(Int,cub.centroid) < j <= sym_cnt
                #     ii = j 
                #     println("was here ------- centroid")
                # elseif sym_cnt- convert(Int,cub.centroid) < j <= sym_cnt 
                #     j = sym_cnt - convert(Int,cub.centroid)
                # end
                sym_cnt -= convert(Int,cub.centroid)

                if cub.numS111 > sym_min[3] && sym_cnt-cub.numS111 < j <= sym_cnt
                    s = sortperm(cub.weights[sym_cnt+1-cub.numS111:sym_cnt],rev=true)
                    # s = sortperm(ss[sym_cnt+1+convert(Int,cub.centroid)-cub.numS111:sym_cnt+convert(Int,cub.centroid)],rev=rev_sort)
                    ii = s[j-(sym_cnt-cub.numS111)]+(sym_cnt-cub.numS111)
                    # ii = j 
                    println("was here ------- S111")
                elseif sym_cnt-cub.numS111 < j <= sym_cnt
                    j = sym_cnt-cub.numS111
                end
                sym_cnt -= cub.numS111
                if cub.numS21 > sym_min[2] && sym_cnt-cub.numS21 < j <= sym_cnt
                    s = sortperm(cub.weights[sym_cnt+1-cub.numS21:sym_cnt],rev=true)
                    # s = sortperm(ss[sym_cnt+1+convert(Int,cub.centroid)-cub.numS21:sym_cnt+convert(Int,cub.centroid)],rev=rev_sort)
                    ii = s[j-(sym_cnt-cub.numS21)]+(sym_cnt-cub.numS21)
                    # ii = j 
                    println("was here ------- S21")
                elseif sym_cnt-cub.numS21 < j <= sym_cnt
                    j = sym_cnt-cub.numS21
                end
                sym_cnt -= cub.numS21

                if j <= 0
                    break
                end 
            else 
                # s = sortperm(cub.weights,rev=true)
                # ii = s[j]
                ii = sortperm(vec(w[idx]),rev=true)[j]
                # ii = sortperm(vec(ss), rev=rev_sort)[j]
            end

            sg=splice!(sym_group, ii)
            params = copy(cub.params)
            if sg ∉ ["vertices", "midedges", "centroid"]
                splice!(params, param_idx[ii]) 
            end
            weights = copy(cub.weights)
            if sg=="centroid" && k<=3
                # don't splice
                push!(sym_group,"centroid")
            else
                splice!(weights, weight_idx[ii]) 
            end
            cub2 = SymCubatures.TriSymCub{T}(vertices=convert(Bool,count(x -> x == "vertices", sym_group)),
                                        midedges=convert(Bool,count(x -> x == "midedges", sym_group)),
                                        numS21=count(x -> x == "S21", sym_group),
                                        numedge=count(x -> x == "edge", sym_group),
                                        numS111=convert(Int,floor(count(x -> x == "S111", sym_group))),
                                        centroid =convert(Bool,count(x -> x == "centroid", sym_group)))
            # println(params)
            # println(cub2.params)
            SymCubatures.setparams!(cub2, params)
            SymCubatures.setweights!(cub2, weights)

            if cub2.numnodes>=1
                xinit = convert.(T,collect(Iterators.flatten([cub2.params,cub2.weights])))
                mask = Int[]
                append!(mask, 1:cub2.numparams+cub2.numweights) 
                # res = Cubature.solvecubature!(cub2, q, mask, tol=5e-14, hist=true, verbose=true,  xinit=xinit, delta1=1e-6, delta2=1e-6)
                res = Cubature.solvecubaturelma!(cub2, q, mask, tol=5e-14, hist=true, verbose=true,  maxiter=10, xinit=xinit, nu=nu)

                res_old=copy(res)
                if res > 5e-14 && res < 1e4
                    xinit = convert.(T,collect(Iterators.flatten([cub2.params,cub2.weights])))
                    res = Cubature.solvecubaturelma!(cub2, q, mask, tol=5e-12, hist=true, verbose=true,  maxiter=10, xinit=xinit, nu=nu)
                end

                # kk = 0
                for i=1:1000
                    if ((res_old-res)/res_old>1e-2 && res>5e-5 && res <1e4)
                        res_old=res
                        xinit = convert.(T,collect(Iterators.flatten([cub2.params,cub2.weights])))
                        res = Cubature.solvecubaturelma!(cub2, q, mask, tol=5e-14, hist=true, verbose=true,  maxiter=10, xinit=xinit, nu=nu)
                    end
                    if (res<5e-5 && res>1e-8)
                        xinit = convert.(T,collect(Iterators.flatten([cub2.params,cub2.weights])))
                        res = Cubature.solvecubaturelma!(cub2, q, mask, tol=5e-14, hist=true, verbose=true,  maxiter=10, xinit=xinit, nu=nu)
                    end
                end
                
                j-=1
                if res <= 1e-8
                    n = cub2.numweights 
                    cub = cub2
                    res_min = res
                end
                println("\n", cub2.params,"\n")
                println(cub2.weights,"\n")
                println("sym_group : ", "[", convert(Int, cub2.centroid),",",cub2.numS21,",",cub2.numS111,"]", ",   numnodes = ",cub2.numnodes)
                # println(cub2,"\n")
                println("k = ", k, ":  j = ", j, ":  n = ", n, ":  numnodes = ", cub.numnodes)
            else
                break 
            end
        end
    end
    println("--------------------------")
    if res_min==-1.0
        println("\n", cub.params,"\n")
        println(cub.weights,"\n")
        println("sym_group : ", "[", convert(Int, cub.centroid),",",cub.numS21,",",cub.numS111,"]", ",   numnodes = ",cub.numnodes)
        # println(cub,"\n")
        print("\n","No new solution was found.","\n")
    else
        println("res norm = ", res_min)
        println("\n", cub.params,"\n")
        println(cub.weights,"\n")
        println("sym_group : ", "[", convert(Int, cub.centroid),",",cub.numS21,",",cub.numS111,"]", ",   numnodes = ",cub.numnodes)
        # println(cub,"\n")
    end
    return cub, res_min
end 

function eliminate_nodes(cub::TetSymCub{T}, p::Int, q::Int) where {T}
    vtx = T[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
 
    res_min = -1.0
    nu = 1e0 #[1e-2,1e0,1e-2,1e0]
    for k = 1:4
        n = cub.numweights
        j=n
        while (j >= 1 && cub.numnodes>=1)
            idx = []
            param_idx = []
            weight_idx = []
            sym_group = String[]
            iig = 0
            iig_idx = []
            ptr = 0
            paramptr = 0
            weightptr = 0
            if cub.vertices
                append!(idx, ptr+1)
                ptr += 4
                append!(weight_idx, weightptr+1)
                weightptr +=1
                push!(sym_group, "vertices")
            end
            for i=1:cub.numS31
                append!(idx, ptr+1)
                ptr += 4
                append!(param_idx, [paramptr+1])
                paramptr += 1
                append!(weight_idx, weightptr+1)
                weightptr +=1
                push!(sym_group, "S31")
                iig += 1
                push!(iig_idx, iig)
            end
            if cub.midedges
                append!(idx, ptr+1)
                ptr += 6
                append!(weight_idx, weightptr+1)
                weightptr +=1
                push!(sym_group, "midedges")
            end
            for i=1:cub.numS22
                append!(idx, ptr+1)
                ptr += 6
                append!(param_idx, [paramptr+1])
                paramptr += 1
                append!(weight_idx, weightptr+1)
                weightptr +=1
                push!(sym_group, "S22")
                iig += 1
                push!(iig_idx, iig)
            end
            for i=1:cub.numfaceS21
                append!(idx, ptr+1)
                ptr += 3
                append!(param_idx, [paramptr+1])
                paramptr += 1
                append!(weight_idx, weightptr+1)
                weightptr +=1
                push!(sym_group, "faceS21")
                iig += 1
                push!(iig_idx, iig)
            end
            for i = 1:cub.numedge
                append!(idx, ptr+1)
                ptr += 6
                append!(param_idx, [paramptr+1])
                paramptr += 1
                append!(weight_idx, weightptr+1)
                weightptr +=1
                push!(sym_group, "edge")
                iig += 1
                push!(iig_idx, iig)
            end
            for i=1:cub.numS211
                append!(idx, ptr+1)
                ptr += 12
                append!(param_idx, [paramptr+1:paramptr+2])
                paramptr += 2
                append!(weight_idx, weightptr+1)
                weightptr +=1
                push!(sym_group, "S211")
                iig += 1
                push!(iig_idx, iig, iig)
            end
            for i = 1:cub.numfaceS111
                append!(idx, ptr+1)
                ptr += 6
                append!(param_idx, [paramptr+1:paramptr+2])
                paramptr += 2
                append!(weight_idx, weightptr+1)
                weightptr +=1
                push!(sym_group, "faceS111")
                iig += 1
                push!(iig_idx, iig, iig)
            end
            if cub.facecentroid
                append!(idx, ptr+1)
                ptr += 1
                append!(weight_idx, weightptr+1)
                weightptr +=1
                push!(sym_group, "facecentroid")
            end
            for i = 1:cub.numS1111
                append!(idx, ptr+1)
                ptr += 24
                append!(param_idx, [paramptr+1:paramptr+3])
                paramptr += 3
                append!(weight_idx, weightptr+1)
                weightptr +=1
                push!(sym_group, "S1111")
                iig += 1
                push!(iig_idx, iig, iig, iig)
            end
            if cub.centroid
                append!(idx, ptr+1)
                ptr += 1
                append!(weight_idx, weightptr+1)
                weightptr +=1
                push!(sym_group, "centroid")
            end
            xy = SymCubatures.calcnodes(cub, vtx)
            xy_sym = xy[:,idx]

            # compute the Vandermonde matrix
            # V, _,_,_= OrthoPoly.vandermonde(p, xy_sym[1,:],xy_sym[2,:],xy_sym[3,:],compute_grad=false)
            # V, _,_,_ = OrthoPoly.vandermonde_monomial(p, xy_sym[1,:],xy_sym[2,:],compute_grad=false,compute_integ=false)
            # V, _,_,_= OrthoPoly.vandermonde_arnoldi(p, xy_sym[1,:],xy_sym[2,:],compute_grad=false)
            # V2 = V.^2*ones(size(V,2),1)
            # V2 = abs.(V)*ones(size(V,2),1)
            # V2 = V2./maximum(V2)
            w = SymCubatures.calcweights(cub)
            # w = w./maximum(w)
            # ss = vec(diagm(w[idx])*V2)
            # ss = vec(V2)
            # ss = vec(cub.weights)
            # ss = w[idx]

            # ss = []
            # for m=1:cub.numweights
            #     num_eq = convert(Int, (q+1)*(q+2)/2)
            #     F = zeros(T, (num_eq) )
            #     F[1] = -2.0/sqrt(2.0)
            #     xy = xy_sym[:, setdiff(1:end, m)]
            #     ww = w[idx][setdiff(1:end,m),:] 
            #     V, _,_= OrthoPoly.vandermonde(q, xy[1,:], xy[2,:], compute_grad=false)
            #     F = V'*ww + F
            #     push!(ss, norm(F))
            # end
            # println(ss)
            # ss =vec(diagm(w[idx])*ss)

            # ss = []
            # for m=1:cub.numweights 
            #     xy = xy_sym[:, setdiff(1:end, m)]
            #     mindist = convert(T, Inf)
            #     if cub.numnodes > 1
            #         for i = 1:size(xy,2)
            #             for j = i+1:size(xy,2)
            #                 mindist = min(mindist, norm(xy[:,i] - xy[:,j],2))
            #             end
            #         end   
            #     end 
            #     # ww = w[idx][setdiff(1:end,m),:] 
            #     push!(ss, mindist)
            # end    
            # println(ss)
            # ss =vec(diagm(w[idx])*ss)

            # ssd = []
            # for m=1:cub.numnodes
            #     xy_m = xy[:, setdiff(1:end, m)]
            #     mindist = convert(T, Inf)
            #     if cub.numnodes > 1
            #         for im = 1:size(xy_m,2)
            #             for jm = im+1:size(xy_m,2)
            #                 mindist = min(mindist, norm(xy_m[:,im] - xy_m[:,jm],2))
            #             end
            #         end   
            #     end 
            #     # ww = w[idx][setdiff(1:end,m),:] 
            #     push!(ssd, mindist)
            # end
            # ss = ssd[idx]
            # # ss = w[idx]./(min.(w[idx]))+ss./(min.(ss))
            # ss =vec(diagm(w[idx])*(ss))

            ssf = []
            for m=1:cub.numnodes
                xy_m = xy[:, setdiff(1:end, m)]
                mindist = convert(T, Inf)
                if cub.numnodes > 1
                    mindist=calc_mindist_to_facet(xy_m)
                end 
                # ww = w[idx][setdiff(1:end,m),:] 
                push!(ssf, mindist)
            end
            ss = ssf[idx]
            # ss = ssf[idx].*w[idx]
            # ss = (ssd./minimum(ssd)) .* (ssf./minimum(ssf))
            
            # w = SymCubatures.calcweights(cub)
            # ss = sortperm(vec(diagm(w[idx])*V2),rev=false)
            # ss = sortperm(vec(diagm(w[idx])*V2),rev=true)
            # ss = sortperm(vec(V2),rev=false)
            # ss = sortperm(vec(V2),rev=true)
            # s = sortperm(cub.weights,rev=false)
            # ss = sortperm(vec(diagm(w[idx])*V2),rev=true)
            # s = cub.weights

            # ss = sym_group[s]
            # nk = convert(Int,floor(cub.numweights)) 
            # ssnk = ss[1:nk]
            # ng =0
            # ind_S1111 = findall(x -> x == "S1111", ssnk)
            # if length(ind_S1111) !=0
            #     ng+=length(ind_S1111)
            #     s[1:ng]=s[ind_S1111]
            # end
            # ind_S211 = findall(x -> x == "S211", ssnk)
            # if length(ind_S211) !=0
            #     s[ng+1:ng+length(ind_S211)]=s[ind_S211]
            #     ng+=length(ind_S211)
            # end
            # ind_S22 = findall(x -> x == "S22", ssnk)
            # if length(ind_S22) !=0
            #     s[ng+1:ng+length(ind_S22)]=s[ind_S22]
            #     ng+=length(ind_S22)
            # end
            # ind_S31 = findall(x -> x == "S31", ssnk)
            # if length(ind_S31) !=0
            #     s[ng+1:ng+length(ind_S31)]=s[ind_S31]
            #     ng+=length(ind_S31)
            # end
            
            # ind_cent= findall(x -> x == "centroid", ssnk)
            # if length(ind_cent) !=0
            #     s[ng+1:ng+length(ind_cent)]=s[ind_cent]
            #     ng+=length(ind_cent)
            # end
            rev_sort = false
            if (k<=2 && cub.numnodes>1)
                nmin,sym_min = minnodes_bound_omega_tet(q)
                sym_cnt = 0
                ii=j - convert(Int, cub.centroid)
                sym_cnt=cub.numweights
                # if cub.centroid && sym_min[1]==0 && sym_cnt-convert(Int,cub.centroid) < j <= sym_cnt
                #     ii = j 
                #     println("was here ------- centroid")
                # elseif sym_cnt- convert(Int,cub.centroid) < j <= sym_cnt 
                #     j = sym_cnt - convert(Int,cub.centroid)
                # end
                sym_cnt -= convert(Int,cub.centroid)

                if cub.numS1111 > sym_min[5] && sym_cnt-cub.numS1111 < j <= sym_cnt
                    # s = sortperm(cub.weights[sym_cnt+1+convert(Int,cub.centroid)-cub.numS1111:sym_cnt+convert(Int,cub.centroid)],rev=true)
                    s = sortperm(ss[sym_cnt+1+convert(Int,cub.centroid)-cub.numS1111:sym_cnt+convert(Int,cub.centroid)],rev=rev_sort)
                    # s = sortperm(ss[sym_cnt+1-cub.numS1111:sym_cnt],rev=true)
                    ii = s[j-(sym_cnt-cub.numS1111)]+(sym_cnt-cub.numS1111)
                    println("was here ------- S1111")
                elseif sym_cnt-cub.numS1111 < j <= sym_cnt
                    j = sym_cnt-cub.numS1111
                end
                sym_cnt -= cub.numS1111

                if cub.numS211 > sym_min[4] && sym_cnt-cub.numS211 < j <= sym_cnt
                    # s = sortperm(cub.weights[sym_cnt+1+convert(Int,cub.centroid)-cub.numS211:sym_cnt+convert(Int,cub.centroid)],rev=true)
                    s = sortperm(ss[sym_cnt+1+convert(Int,cub.centroid)-cub.numS211:sym_cnt+convert(Int,cub.centroid)],rev=rev_sort)
                    # s = sortperm(ss[sym_cnt+1-cub.numS211:sym_cnt],rev=true)
                    ii = s[j-(sym_cnt-cub.numS211)]+(sym_cnt-cub.numS211)
                    println("was here ------- S211")
                elseif sym_cnt-cub.numS211 < j <= sym_cnt
                    j = sym_cnt-cub.numS211
                end
                sym_cnt -= cub.numS211
                if cub.numS22 > sym_min[3] && sym_cnt-cub.numS22 < j <= sym_cnt
                    # s = sortperm(cub.weights[sym_cnt+1+convert(Int,cub.centroid)-cub.numS22:sym_cnt+convert(Int,cub.centroid)],rev=true)
                    s = sortperm(ss[sym_cnt+1+convert(Int,cub.centroid)-cub.numS22:sym_cnt+convert(Int,cub.centroid)],rev=rev_sort)
                    # s = sortperm(ss[sym_cnt+1-cub.numS22:sym_cnt],rev=true)
                    ii = s[j-(sym_cnt-cub.numS22)]+(sym_cnt-cub.numS22)
                    # pp = cub.params[cub.numS31+1:cub.numS31+cub.numS22]
                    # ii = argmax(norm.(pp.-1/4))
                    println("was here ------- S22")
                elseif sym_cnt-cub.numS22 < j <= sym_cnt
                    j = sym_cnt-cub.numS22
                end
                sym_cnt -= cub.numS22
                if cub.numS31 > sym_min[2] && sym_cnt-cub.numS31 < j <= sym_cnt
                    # s = sortperm(cub.weights[sym_cnt+1+convert(Int,cub.centroid)-cub.numS31:sym_cnt+convert(Int,cub.centroid)],rev=true)
                    s = sortperm(ss[sym_cnt+1+convert(Int,cub.centroid)-cub.numS31:sym_cnt+convert(Int,cub.centroid)],rev=rev_sort)
                    # s = sortperm(ss[sym_cnt+1-cub.numS31:sym_cnt],rev=true)
                    ii = s[j-(sym_cnt-cub.numS31)]+(sym_cnt-cub.numS31)
                    # pp = cub.params[1:cub.numS31]
                    # ii = argmax(norm.(pp.-1/4))
                    println("was here ------- S31")
                elseif sym_cnt-cub.numS31 < j <= sym_cnt
                    j = sym_cnt-cub.numS31
                end
                sym_cnt -= cub.numS31

                if j <= 0
                    break
                end 
            else 
                # s = sortperm(cub.weights,rev=true)
                # V, _,_,_= OrthoPoly.vandermonde(p, xy_sym[1,:],xy_sym[2,:],xy_sym[3,:],compute_grad=false)
                # V2 = abs.(V)*ones(size(V,2),1)
                # s = sortperm(vec(V2),rev=true)
                # V2 = V.^2*ones(size(V,2),1)
                # w = SymCubatures.calcweights(cub)
                # s = sortperm(vec(diagm(w[idx])*V2),rev=true)
                # s = 1:length(cub.weights)
                ii = sortperm(vec(ss), rev=rev_sort)[j]
                # ii = j
                # ii = sortperm(vec(w[idx]),rev=true)[j]
            end

            # ii = s[j]
            sg=splice!(sym_group, ii)
            params = copy(cub.params)
            if sg ∉ ["vertices", "midedges", "centroid", "facecentroid"]
                # splice!(params, param_idx[iig_idx[ii]]) 
                splice!(params, param_idx[ii]) 
            end
            weights = copy(cub.weights)
            if sg=="centroid" && k<=2
                # don't splice
                push!(sym_group,"centroid")
            else
                splice!(weights, weight_idx[ii]) 
            end
            cub2 = SymCubatures.TetSymCub{T}(vertices=convert(Bool,count(x -> x == "vertices", sym_group)),
                                            midedges=convert(Bool,count(x -> x == "midedges", sym_group)),
                                            numfaceS21=count(x -> x == "faceS21", sym_group),
                                            numedge=count(x -> x == "edge", sym_group),
                                            numfaceS111=convert(Int,floor(count(x -> x == "faceS111", sym_group))),
                                            facecentroid =convert(Bool,count(x -> x == "facecentroid", sym_group)), 
                                            numS31=count(x ->x == "S31", sym_group),
                                            numS22=count(x ->x == "S22", sym_group),
                                            numS211=convert(Int,floor(count(x ->x == "S211", sym_group))),
                                            numS1111=convert(Int,floor(count(x ->x == "S1111", sym_group))),
                                            centroid =convert(Bool,count(x -> x == "centroid", sym_group)))
    
            SymCubatures.setparams!(cub2, params)
            SymCubatures.setweights!(cub2, weights)

            if cub2.numnodes>=1
                xinit = convert.(T,collect(Iterators.flatten([cub2.params,cub2.weights])))
                mask = Int[]
                append!(mask, 1:cub2.numparams+cub2.numweights) 
                # res = Cubature.solvecubature!(cub2, q, mask, tol=5e-14, hist=true, verbose=true,  xinit=xinit, delta1=1e-4, delta2=1e-4)
                res = Cubature.solvecubaturelma!(cub2, q, mask, tol=5e-14, hist=true, verbose=true,  maxiter=10, xinit=xinit, nu=nu)

                res_old=copy(res)
                if res > 5e-14 && res < 1e4
                    xinit = convert.(T,collect(Iterators.flatten([cub2.params,cub2.weights])))
                    res = Cubature.solvecubaturelma!(cub2, q, mask, tol=5e-14, hist=true, verbose=true,  maxiter=10, xinit=xinit, nu=nu)
                end
                res_ratio=res/res_old 
                println("res_here = ", res)
                # res_old=res
                for i=1:100
                    # # if ((res_ratio <= 9.8e-1 && res > 1e-2) || (res_ratio <= 9.99e-1 && res < 1e-2 && res > 1e-10))
                    # if ((res_old-res)/res_old>1e-2 && res>1e-8 && res <1e4)
                    #     res_old=res
                    #     xinit = convert.(T,collect(Iterators.flatten([cub2.params,cub2.weights])))
                    #     res = Cubature.solvecubaturelma!(cub2, q, mask, tol=5e-14, hist=true, verbose=true,  maxiter=10, xinit=xinit, delta1=1e-4, delta2=1e-4)
                    #     # res_ratio=res/res_old 
                    #     # res_old=res
                    # end

                    if ((res_old-res)/res_old>1e-2 && res>5e-4 && res <1e4)
                        res_old=res
                        xinit = convert.(T,collect(Iterators.flatten([cub2.params,cub2.weights])))
                        res = Cubature.solvecubaturelma!(cub2, q, mask, tol=5e-14, hist=true, verbose=true,  maxiter=10, xinit=xinit, nu=nu)
                    end
                    
                    if (res<5e-4 && res>1e-8)
                        xinit = convert.(T,collect(Iterators.flatten([cub2.params,cub2.weights])))
                        res = Cubature.solvecubaturelma!(cub2, q, mask, tol=5e-14, hist=true, verbose=true,  maxiter=10, xinit=xinit, nu=nu)
                    end
                end
                
                j-=1
                if res <= 1e-8
                    n = cub2.numweights 
                    # j = 1
                    # j += 1
                    cub = cub2
                    res_min = res
                end
                println("\n", cub2.params,"\n")
                println(cub2.weights,"\n")
                println("sym_groups : ", "[", convert(Int, cub2.centroid),",",cub2.numS31,",",cub2.numS22,",",cub2.numS211,",",cub2.numS1111,"]", ",   numnodes = ",cub2.numnodes)
                # println(cub2,"\n")
                println("k = ", k, ":  j = ", j, ":  n = ", n, ":  numnodes = ", cub.numnodes)
            else
                break 
            end
        end
    end
    println("--------------------------")
    if res_min==-1.0
        println("\n", cub.params,"\n")
        println(cub.weights,"\n")
        println("sym_groups : ", "[", convert(Int, cub.centroid),",",cub.numS31,",",cub.numS22,",",cub.numS211,",",cub.numS1111,"]", ",   numnodes = ",cub.numnodes)
        # println(cub,"\n")
        print("\n","No new solution was found.","\n")
    else
        println("res norm = ", res_min)
        println("\n", cub.params,"\n")
        println(cub.weights,"\n")
        # println(cub,"\n")
        println("sym_groups : ", "[", convert(Int, cub.centroid),",",cub.numS31,",",cub.numS22,",",cub.numS211,",",cub.numS1111,"]", ",   numnodes = ",cub.numnodes)
    end
    return cub, res_min
end 

function combine_tet_cubs(cuba::TetSymCub{T}, cubb::TetSymCub{T}) where {T}

    cubs = [cuba, cubb]

    w_vertices=[]
    w_S31= []
    w_midedges= []
    w_S22= []
    w_faceS21= []
    w_edge= []
    w_S211= []
    w_faceS111= []
    w_facecentroid= []
    w_S1111= []
    w_centroid= []

    p_S31= []
    p_S22= []
    p_faceS21= []
    p_edge= []
    p_S211= []
    p_faceS111= []
    p_S1111= []

    sym_group=[]

    for cub in cubs
        w_ptr=0
        p_ptr=0
        ws = cub.weights
        ps = cub.params
        if (cub.vertices && w_vertices==[])
            append!(w_vertices, ws[1])
            w_ptr+=1
            push!(sym_group, "vertices")
        end
        for i=1:cub.numS31
            append!(w_S31,ws[w_ptr+1])
            w_ptr += 1
            append!(p_S31,ps[p_ptr+1])
            p_ptr += 1
            push!(sym_group, "S31")
        end
        if (cub.midedges && w_midedges==[])
            append!(w_midedges,ws[w_ptr+1])
            w_ptr += 1
            push!(sym_group, "midedges")
        end
        for i=1:cub.numS22
            append!(w_S22,ws[w_ptr+1])
            w_ptr += 1
            append!(p_S22,ps[p_ptr+1])
            p_ptr += 1
            push!(sym_group, "S22")
        end
        for i=1:cub.numfaceS21
            append!(w_faceS21,ws[w_ptr+1])
            w_ptr += 1
            append!(p_faceS21,ps[p_ptr+1])
            p_ptr += 1
            push!(sym_group, "faceS21")
        end
        for i = 1:cub.numedge
            append!(w_edge,ws[w_ptr+1])
            w_ptr += 1
            append!(p_edge,ps[p_ptr+1])
            p_ptr += 1
            push!(sym_group, "edge")
        end
        for i=1:cub.numS211
            append!(w_S211,ws[w_ptr+1])
            w_ptr += 1
            append!(p_S211,ps[p_ptr+1:p_ptr+2])
            p_ptr += 2
            push!(sym_group, "S211", "S211")
        end
        for i = 1:cub.numfaceS111
            append!(w_faceS111,ws[w_ptr+1])
            w_ptr += 1
            append!(p_faceS111,ps[p_ptr+1:p_ptr+2])
            p_ptr += 2
            push!(sym_group, "faceS111", "faceS111")
        end
        if (cub.facecentroid && w_facecentroid==[])
            append!(w_facecentroid,ws[w_ptr+1])
            w_ptr += 1
            push!(sym_group, "facecentroid")
        end
        for i = 1:cub.numS1111
            append!(w_S1111,ws[w_ptr+1])
            w_ptr += 1
            append!(p_S1111,ps[p_ptr+1:p_ptr+3])
            p_ptr += 3
            push!(sym_group, "S1111", "S1111", "S1111")
        end
        if (cub.centroid && w_centroid==[])
            append!(w_centroid,ws[w_ptr+1])
            w_ptr += 1
            push!(sym_group, "centroid")
        end
    end
    
    weights = convert.(T,collect(Iterators.flatten([w_vertices,w_S31,w_midedges,w_S22,w_faceS21,w_edge,w_S211,w_faceS111,w_facecentroid,w_S1111,w_centroid])))
    params = convert.(T,collect(Iterators.flatten([p_S31,p_S22,p_faceS21,p_edge,p_S211,p_faceS111,p_S1111]))) 
    
    cub2 = SymCubatures.TetSymCub{T}(vertices=convert(Bool,count(x -> x == "vertices", sym_group)),
                                    midedges=convert(Bool,count(x -> x == "midedges", sym_group)),
                                    numfaceS21=count(x -> x == "faceS21", sym_group),
                                    numedge=count(x -> x == "edge", sym_group),
                                    numfaceS111=convert(Int,floor(count(x -> x == "faceS111", sym_group)/2)),
                                    facecentroid =convert(Bool,count(x -> x == "facecentroid", sym_group)), 
                                    numS31=count(x ->x == "S31", sym_group),
                                    numS22=count(x ->x == "S22", sym_group),
                                    numS211=convert(Int,floor(count(x ->x == "S211", sym_group)/2)),
                                    numS1111=convert(Int,floor(count(x ->x == "S1111", sym_group)/3)),
                                    centroid =convert(Bool,count(x -> x == "centroid", sym_group)))
    SymCubatures.setparams!(cub2, params)
    SymCubatures.setweights!(cub2, weights)

    return cub2
end

function minnodes_bound_omega_tet(d::Int)
    d12 = d - 12
    mp3 = Int(round(((d12+4)^3 + 3*(d12+4)^2 - 9*(d12+4)*mod(d12+4,2))/144))*iverson(d12,0)
    me = Int(round(((d+4)^3 + 3*(d+4)^2 - 9*(d+4)*mod(d+4,2))/144))*iverson(d,0)
    mp2 = Int(round((d+3)^2/12))*iverson(d,0)
    mp1 = Int(floor((d+2)/2))*iverson(d,0)
    mp0 = iverson(d,0)
    m4 = mp3
    m3 = Int(floor((d/2-2)^2))*iverson(d,6)
    m2 = Int(floor(d/2-1))*iverson(d,4)
    m1 = (d-2)*iverson(d,2)
    m12 = iverson(d,2)
    m0 = 1

    n4 = Int(ceil(m4/4))
    n3 = Int(ceil((m4+m3 - 4*n4)/3))
    n2 = Int(ceil((m4+m3+m2-3*n3-4*n4)/2))
    n1 = Int(floor((me-2*n2-3*n3-4*n4)/2))
    n0 = me - 2*n1 - 2*n2 - 3*n3 - 4*n4

    n = n0 + 4*n1 + 6*n2 + 12*n3 + 24*n4

    # println("m0=", m0, ": m12=", m12, ": m1=", m1, ": m2=", m2, ": m3=", m3, ": m4=", m4, ": me=", me)
    # println("omega_sym_min: ", "[", n0,",",n1,",",n2,",",n3,",",n4,"]", ",   numnodes = ",n )
    sym_group = [n0,n1,n2,n3,n4]
    return n, sym_group
end

function iverson(d::Int, a::Int)
    if d >= a 
        return 1
    else
        return 0
    end
end

function numnodes_omega_lg_tri(q::Int)
    p = Int(ceil(q/2))
    ne=p+1
    numedge = convert(Int,(p+1-mod(p+1,2))/2)
    numS21= n1 = convert(Int,(1+mod(ne,2))*numedge)
    numS111= n2 = convert(Int,1/2*(numedge^2 - numedge))
    centroid = n0 = convert(Int,mod(ne,2))
    
    n = centroid + 3*numS21 + 6*numS111
    println("omega_sym_lg : ", "[", n0,",",n1,",",n2,"]", ",   numnodes = ",n )
    
    nmin,sym_min = minnodes_bound_omega_tri(q)
    println("omega_sym_min: ", "[", sym_min[1],",",sym_min[2],",",sym_min[3],"]", ",   numnodes = ",nmin )
    r= n/nmin
    println("ratio = ", round(r,digits=4))
    return n
end

function numnodes_omega_lg_tet(q::Int)
    p = Int(ceil(q/2))
    ne=p+1
    numedge=Int((ne-mod(ne,2))/2)
    numS31=n1=Int((1+mod(ne,2))*numedge)
    numS22 =n2= Int(mod(ne,2)*numedge)
    numS211=n3= Int((1+2*mod(ne,2))/(1+mod(ne,2))*(numedge*(numedge-1)))
    numS1111=n4=Int(1/6*((numedge-1)^3-(numedge-1)))
    centroid=n0= mod(ne,2)

    n = centroid + 4*numS31 + 6*numS22 + 12*numS211 + 24*numS1111
    println("omega_sym_lg : ", "[", n0,",",n1,",",n2,",",n3,",",n4,"]", ",   numnodes = ",n )
    
    nmin,sym_min = minnodes_bound_omega_tet(q)
    println("omega_sym_min: ", "[", sym_min[1],",",sym_min[2],",",sym_min[3],",",sym_min[4],",",sym_min[5],"]", ",   numnodes = ",nmin )
    r= n/nmin
    println("ratio = ", round(r,digits=4))
    return n
end

"""
### minnodes_bound_omega_tri

A function that finds the minimum possible number and corresponding symmetry
group for symmetric quadrature rule on the triangle. This is based on the work 
of Lyness and Jespersen, Moderate Degree Symmetric Quadrature Rules for the Triangle (1975)

**Inputs**
* `q`: degree of the quadrature rule 

**Outputs**
* `n`: minimum number of nodes

"""
function minnodes_bound_omega_tri(d::Int)
    alphas = [3, -4, -1, 0, -1, -4]
    alpha_d = alphas[mod(d,6)+1]

    function E(d, alpha_d)
        return convert(Int, 1/12*((d+3)^2 + alpha_d))
    end

    n2=0
    if d>=6
        n2 = convert(Int,floor((E(d-6,alpha_d) + 2)/3))
    end
    n1 = convert(Int,floor((E(d,alpha_d)-3*n2)/2))
    n0 = 1
    if n0+2*n1+3*n2 > E(d,alpha_d)
        n0 = 0
    end

    n = n0 + 3*n1 + 6*n2 
    sym_group = [n0,n1,n2]

    return n, sym_group
end

"""
SummationByParts.sym_group_for_elimination

Node elimination algorithm based on Slobodkins and Tausch (2022)

**Inpouts**
* `cub`: cubature data
* `q`: the degree of the quadrature rule

**Outputs**
* `z`: unsorted step size for the linearization of the objective function
"""
function sym_group_for_elimination(cub::SymCub{T}, q::Int) where{T}
    Jac = SymCubatures.calcjacobian(cub)
    _, dF = Cubature.cubatureresidual(cub, q, compute_grad=true)
    J = dF*Jac
    Q = nullspace(J')
    w = cub.weights
    z = zeros(1,length(w))

    for k = 1:length(w)
        yk = Q[:,k] .*(w[k]/(norm(Q[:,k])^2))
        z[k] = (Q'*yk)[k]
    end
    println(z)
    return z
end

function calc_mindist_to_facet(x::Array{T}) where T
    #assume vertex is vtx= T[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
    mindist = minimum([minimum(abs.(x.-1)),minimum(abs.((x'*ones(size(x,1),1)).+1))])
    return mindist
end