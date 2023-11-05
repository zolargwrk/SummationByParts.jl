using Interpolations
using SummationByParts.Cubature

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

    vtx = T[0 0; 1 0; 1/2 sqrt(3/4)]
    xy = SymCubatures.calcnodes(cub, vtx)
    lines = SummationByParts.sparse_stencil(cub, p)
    node_stencil = zeros(Int64,cub.numnodes,p+2,2)

    for i = 1:cub.numnodes
        idx = findall(x -> i in x, eachrow(lines))
        line1 = lines[idx[1],:]
        line2 = lines[idx[2],:]

        line1_x = line1[sortperm(xy[1,line1])]
        line2_x = line2[sortperm(xy[1,line2])]

        x = line1_x
        y = line2_x
        if xy[1,x[1]] > xy[1,y[1]]
            x = line2_x
            y = line1_x
        end
        if (xy[2,x[1]] < 1e-10 && xy[1,x[1]] > 1e-10)
            xtemp = copy(x)
            x = y
            y = xtemp
        end

        node_stencil[i,:,:] = Array([x;y])'
    end

    return node_stencil
end


"""
### SummationByParts.sparse_metric_terms

Computes the metric terms at each node for sparse SBP operators. 
The metric terms are approximated using 1D SBP operator.

**Inputs**

* `cub`: symmetric cubature rule
* `p`: maximum total degree for the Proriol polynomials

**Outputs**

* `metric_terms`: a matrix containing the stencil in each direction for each node

"""

function sparse_metric_terms(cub::TriSymCub{T}, vtx::Array{T,2}, p::Int) where {T}
    
    ddxi = zeros(Float64,2,cub.numnodes)
    dxid = zeros(Float64,2,cub.numnodes)
    ddeta = zeros(Float64,2,cub.numnodes)
    detad = zeros(Float64,2,cub.numnodes)
    jac = zeros(Float64,1,cub.numnodes)
    node_lines = SummationByParts.node_sparse_stencil(cub, p)
    xy = SymCubatures.calcnodes(cub, vtx)

    oper_lgl = getLineSegSBPLobbato(degree=p+1)
    Q = oper_lgl.Q[:,:,1]
    H = diagm(oper_lgl.w)
    D = inv(H)*Q[:,:,1]
    # writedlm(stdout, round.(D,digits=4))


    for i = 1:cub.numnodes
        for j = 1:2
            ddxi[j,i] = ((D * reshape(xy[j,node_lines[i,:,1]],(p+2,1)))[findall(a->i in a, node_lines[i,:,1])])[1]
            ddeta[j,i] = ((D * reshape(xy[j,node_lines[i,:,2]],(p+2,1)))[findall(a->i in a, node_lines[i,:,2])])[1]
        end
    end

    for i = 1:cub.numnodes
        jac[1,i] = ddxi[1,i]*ddeta[2,i] - ddeta[1,i]*ddxi[2,i]
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

