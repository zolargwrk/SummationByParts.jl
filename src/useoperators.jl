# This file gathers together functions related to using the SBP operators

@doc """
### SummationByParts.calcnodes

This function returns the node coordinates for an SBP operator.  The nodes are
ordered first by vertex, then by edge (with nodes ordered in sequence along the
directed edge), then by face (if appropriate), and then finally by volumn nodes.
This function assumes the element mapping is linear, i.e. edges are lines.

**Inputs**

* `sbp`: an SBP operator
* `vtx`: the vertices that define the element

**Outputs**

* `x`: the node coordinates; 1st dimension is the coordinate, the second the node

"""->
function calcnodes{T}(sbp::TriSBP{T}, vtx::Array{T})
  perm = SummationByParts.getnodepermutation(sbp.cub, sbp.degree)
  x = zeros(T, (2, sbp.numnodes))
  x[1,:], x[2,:] = SymCubatures.calcnodes(sbp.cub, vtx)
  return x[:,perm]
end

function calcnodes{T}(sbp::TetSBP{T}, vtx::Array{T})
  perm = SummationByParts.getnodepermutation(sbp.cub, sbp.degree)
  x = zeros(T, (3, sbp.numnodes))
  x[1,:], x[2,:], x[3,:] = SymCubatures.calcnodes(sbp.cub, vtx)
  return x[:,perm]
end

@doc """
### SummationByParts.weakdifferentiate!

Applies the SBP Q matrix operator to data in `u` and stores the result in `res`.
Different methods are available depending on the rank of `u`:

* For *scalar* fields, it is assumed that `u` is a rank-2 array, with the first
dimension for the local-node index, and the second dimension for the element
index.
* For *vector* fields, `u` is a rank-3 array, with the first dimension for the
index of the vector field, the second dimension for the local-node index, and
the third dimension for the element index.

Naturally, the number of entries in the dimension of `u` (and `res`)
corresponding to the nodes must be equal to the number of nodes in the SBP
operator sbp.

**Inputs**

* `sbp`: an SBP operator type
* `di`: direction index of the operator that is desired (di=1 for Qx, etc)
* `u`: the array that the operator is applied to

**In/Outs**

* `res`: where the result of applying Q[:,:,di] to u is stored

"""->
function weakdifferentiate!{T}(sbp::SBPOperator{T}, di::Int,
                               u::AbstractArray{T,2}, res::AbstractArray{T,2})
  @assert( sbp.numnodes == size(u,1) && sbp.numnodes == size(res,1) )
  @assert( length(u) == length(res) )
  @assert( di > 0 && di <= size(sbp.Q,3) )
  for elem = 1:size(u,2)
    for i = 1:sbp.numnodes
      for j = 1:sbp.numnodes
        res[i,elem] += sbp.Q[i,j,di]*u[j,elem]
      end
    end
  end
end

function weakdifferentiate!{T}(sbp::SBPOperator{T}, di::Int,
                               u::AbstractArray{T,3}, res::AbstractArray{T,3})
  @assert( sbp.numnodes == size(u,2) && sbp.numnodes == size(res,2) )
  @assert( length(u) == length(res) )
  @assert( di > 0 && di <= size(sbp.Q,3) )
  for elem = 1:size(u,3)
    for i = 1:sbp.numnodes
      for j = 1:sbp.numnodes
        for field = 1:size(u,1)
          res[field,i,elem] += sbp.Q[i,j,di]*u[field,j,elem]
        end
      end
    end
  end
end

@doc """
### SummationByParts.differentiate!

Applies the SBP differentiation matrix operator, D, to data in `u` and stores
the result in `res`.  Different methods are available depending on the rank of
`u`:

* For *scalar* fields, it is assumed that `u` is a rank-2 array, with the first
dimension for the local-node index, and the second dimension for the element
index.
* For *vector* fields, `u` is a rank-3 array, with the first dimension for the
index of the vector field, the second dimension for the local-node index, and
the third dimension for the element index.

Naturally, the number of entries in the dimension of `u` (and `res`)
corresponding to the nodes must be equal to the number of nodes in the SBP
operator sbp.

**Inputs**

* `sbp`: an SBP operator type
* `di`: direction index of the operator that is desired (di=1 for Dx, etc)
* `u`: the array that the operator is applied to

**In/Outs**

* `res`: where the result of applying inv(H)*Q[:,:,di] to u is stored

"""->
function differentiate!{T}(sbp::SBPOperator{T}, di::Int,
                                u::AbstractArray{T,2}, res::AbstractArray{T,2})
  @assert( sbp.numnodes == size(u,1) && sbp.numnodes == size(res,1) )
  @assert( length(u) == length(res) )
  @assert( di > 0 && di <= size(sbp.Q,3) )
  Hinv = 1./sbp.w
  for elem = 1:size(u,2)
    for i = 1:sbp.numnodes
      for j = 1:sbp.numnodes
        res[i,elem] += sbp.Q[i,j,di]*u[j,elem]
      end
      res[i,elem] *= Hinv[i]
    end
  end
end

function differentiate!{T}(sbp::SBPOperator{T}, di::Int,
                                u::AbstractArray{T,3}, res::AbstractArray{T,3})
  @assert( sbp.numnodes == size(u,2) && sbp.numnodes == size(res,2) )
  @assert( length(u) == length(res) )
  @assert( di > 0 && di <= size(sbp.Q,3) )
  Hinv = 1./sbp.w
  for elem = 1:size(u,3)
    for i = 1:sbp.numnodes
      for j = 1:sbp.numnodes
        for field = 1:size(u,1)
          res[field,i,elem] += sbp.Q[i,j,di]*u[field,j,elem]
        end
      end
      for field = 1:size(u,1)
        res[field,i,elem] *= Hinv[i]
      end
    end
  end
end

@doc """
### SummationByParts.volumeintegrate!

Applies the SBP mass matrix operator, H, to data in `u` and stores
the result in `res`.  Different methods are available depending on the rank of
`u`:

* For *scalar* fields, it is assumed that `u` is a rank-2 array, with the first
dimension for the local-node index, and the second dimension for the element
index.
* For *vector* fields, `u` is a rank-3 array, with the first dimension for the
index of the vector field, the second dimension for the local-node index, and
the third dimension for the element index.

Naturally, the number of entries in the dimension of `u` (and `res`)
corresponding to the nodes must be equal to the number of nodes in the SBP
operator sbp.

**Inputs**

* `sbp`: an SBP operator type
* `u`: the array that the operator is applied to

**In/Outs**

* `res`: where the result of applying H to u is stored

"""->
function volumeintegrate!{T}(sbp::SBPOperator{T}, u::AbstractArray{T,2},
                            res::AbstractArray{T,2})
  @assert( sbp.numnodes == size(u,1) && sbp.numnodes == size(res,1) )
  @assert( length(u) == length(res) )
  for elem = 1:size(u,2)
    for i = 1:sbp.numnodes
      res[i,elem] += sbp.w[i]*u[i,elem]
    end
  end
end

function volumeintegrate!{T}(sbp::SBPOperator{T}, u::AbstractArray{T,3},
                            res::AbstractArray{T,3})
  @assert( sbp.numnodes == size(u,2) && sbp.numnodes == size(res,2) )
  @assert( length(u) == length(res) )
  for elem = 1:size(u,3)
    for i = 1:sbp.numnodes
      for field = 1:size(u,1)
        res[field,i,elem] += sbp.w[i]*u[field,i,elem]
      end
    end
  end
end

@doc """
### SummationByParts.mappingjacobian!

Evaluates the (scaled) Jacobian of the mapping from reference coordinates to
physical coordinates, as well as the determinant of the Jacobian.  The values
returned in dxidx are scaled by the determinant, so they have the same units as
the boundary measure (i.e. length in 2D, or length^2 in 3D).  This scaling is
adopted, because conservation laws written in conservative form in the reference
frame use the scaled Jacobian.

**Inputs**

* `sbp`: an SBP operator type
* `x`: the physical coordinates; 1st dim = coord, 2nd dim = node, 3rd dim = elem

**In/Outs**

* `dxidx`: the scaled Jacobian of the mapping; 1st dim = ref coord, 2nd dim =
  phys coord, 3rd dim = node, 3rd dim = elem
* `jac`: the determinant of the Jacobian; 1st dim = node, 2nd dim = elem

"""->
function mappingjacobian!{T}(sbp::TriSBP{T}, x::AbstractArray{T,3},
                             dxidx::AbstractArray{T,4}, jac::AbstractArray{T,2})
  @assert( sbp.numnodes == size(x,2) && sbp.numnodes == size(dxidx,3) )
  @assert( size(x,3) == size(dxidx,4) )
  @assert( size(x,1) == 2 && size(dxidx,1) == 2 && size(dxidx,2) == 2 )
  fill!(dxidx, zero(T))
  dxdxi = zeros(T, (2,sbp.numnodes,size(x,3)))
  # compute d(x,y)/dxi and set deta/dx and deta/dy
  differentiate!(sbp, 1, x, dxdxi)
  dxidx[2,1,:,:] = -dxdxi[2,:,:]
  dxidx[2,2,:,:] = dxdxi[1,:,:]
  # compute d(x,y)/deta and set dxi/dx and dxi/dy
  fill!(dxdxi, zero(T))
  differentiate!(sbp, 2, x, dxdxi)
  dxidx[1,2,:,:] = -dxdxi[1,:,:]
  dxidx[1,1,:,:] = dxdxi[2,:,:]
  # compute the determinant of the Jacobian
  for elem = 1:size(x,3)
    for i = 1:sbp.numnodes
      jac[i,elem] = dxidx[1,1,i,elem]*dxidx[2,2,i,elem] - 
      dxidx[1,2,i,elem]*dxidx[2,1,i,elem]
    end
  end
  # check for negative jac here?
end

function mappingjacobian!{T}(sbp::TetSBP{T}, x::AbstractArray{T,3},
                             dxidx::AbstractArray{T,4}, jac::AbstractArray{T,2})
  @assert( sbp.numnodes == size(x,2) && sbp.numnodes == size(dxidx,3) )
  @assert( size(x,3) == size(dxidx,4) )
  @assert( size(x,1) == 3 && size(dxidx,1) == 3 && size(dxidx,2) == 3 )
  numelem = size(x,3)
  dxdxi = zeros(T, (3,sbp.numnodes,numelem,3))
  # calculate the derivative of the coordinates with respect to (xi,eta,zeta)
  # using the SBP operator
  for di = 1:3
    differentiate!(sbp, di, x, sub(dxdxi,:,:,:,di)) 
  end
  fill!(dxidx, zero(T))
  # calculate the metrics: the outer loop calculates the derivatives of
  # coordinates-times-coordinate-derivatives, for example, d/deta( z * d(y)/d
  # zeta), and this contribution is added to the relevant metric.
  for di = 1:3
    it1 = mod(di,3)+1
    it2 = mod(di+1,3)+1
    for elem = 1:numelem
      for i = 1:sbp.numnodes
        for j = 1:sbp.numnodes
          # the 0.5 factor in the SBP entry accounts for the averaging over the
          # two possible analytical formulations
          coeff = 0.5*sbp.Q[i,j,di]
          # this second set of indices (di2, it21, it22) denote the x-, y-,
          # z-coordinate indices, while the first set (di, it1, it2) denotes the xi-,
          # eta-, zeta-indices
          for di2 = 1:3
            it21 = mod(di2,3)+1
            it22 = mod(di2+1,3)+1
            dxidx[it1,di2,i,elem] += 
            coeff*(x[it22,j,elem]*dxdxi[it21,j,elem,it2] -
                   x[it21,j,elem]*dxdxi[it22,j,elem,it2])            
            dxidx[it2,di2,i,elem] += 
            coeff*(x[it21,j,elem]*dxdxi[it22,j,elem,it1] -
                   x[it22,j,elem]*dxdxi[it21,j,elem,it1])
          end
        end
      end
    end
  end
  # scale metrics by norm and calculate the determinant of the mapping
  for elem = 1:numelem
    for i = 1:sbp.numnodes
      dxidx[:,:,i,elem] /= sbp.w[i]
      jac[i,elem] = 1.0/(
                         dxdxi[1,i,elem,1]*dxdxi[2,i,elem,2]*dxdxi[3,i,elem,3] +
                         dxdxi[1,i,elem,2]*dxdxi[2,i,elem,3]*dxdxi[3,i,elem,1] +
                         dxdxi[1,i,elem,3]*dxdxi[2,i,elem,1]*dxdxi[3,i,elem,2] -
                         dxdxi[1,i,elem,1]*dxdxi[2,i,elem,3]*dxdxi[3,i,elem,2] -
                         dxdxi[1,i,elem,2]*dxdxi[2,i,elem,1]*dxdxi[3,i,elem,3] -
                         dxdxi[1,i,elem,3]*dxdxi[2,i,elem,2]*dxdxi[3,i,elem,1])
    end
  end
  # check for negative jac here?
end