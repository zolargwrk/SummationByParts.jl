# This file gathers together functions related to using the SBP operators

@doc """
### SummationByParts.Boundary

Used to identify boundary faces in a finite-element grid.

**Fields**

* `element` : index of the element to which the boundary face belongs
* `face` : the face index of the boundary (local index to the element)

**Example**

To mark face 2 of element 7 to be a boundary face, use `Boundary(7,2)`

"""->
immutable Boundary
  element::UInt64
  face::UInt8
end

@doc """
### SummationByParts.Interface

Used to identify interfaces between elements in a finite-element grid.

**Fields**

* `elementL` : index of the so-called left element in the pair
* `elementR` : index of the so-called right element in the pair
* `faceL` : the face index of the interface with respect to the left element
* `faceR` : the face index of the interface with respect to the right element

**Example**

Consider an interface between elements 2 and 5.  Suppose the interface is on
face 1 of element 2 and face 3 of element 5.  This can be indicated as
`Interface(2,5,1,3)`

"""->
immutable Interface
  elementL::UInt64
  elementR::UInt64
  faceL::UInt8
  faceR::UInt8
end

@doc """
### SummationByParts.getnbrnodeindex

Returns the face-node index on `face.faceR` equivalent to index `i` on
`face.faceL`.

**Inputs**

* `sbp`: an SBP operator
* `face`: an element interface
* `i`: face-node index on `face.faceL`

**Returns**

* `j`: face-node index on `face.faceR`

"""->
function getnbrnodeindex{T}(sbp::TriSBP{T}, face::Interface, i::Int)
  return [2; 1; sbp.numfacenodes:-1:3][i]
end

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

Applies the SBP stiffness matrix (or its transpose) to data in `u` and **adds**
the result to `res`.  Different methods are available depending on the rank of
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
* `di`: direction index of the operator that is desired (di=1 for Qx, etc)
* `u`: the array that the operator is applied to
* `trans` (optional): if true, the transpose operation is applied

**In/Outs**

* `res`: where the result of applying Q[:,:,di] to u is stored

"""->
function weakdifferentiate!{T}(sbp::SBPOperator{T}, di::Int,
                               u::AbstractArray{T,2}, res::AbstractArray{T,2};
                               trans::Bool=false)
  @assert( sbp.numnodes == size(u,1) && sbp.numnodes == size(res,1) )
  @assert( length(u) == length(res) )
  @assert( di > 0 && di <= size(sbp.Q,3) )
  if trans # apply transposed Q
    for elem = 1:size(u,2)
      for i = 1:sbp.numnodes
        for j = 1:sbp.numnodes
          res[i,elem] += sbp.Q[j,i,di]*u[j,elem]
        end
      end
    end
  else # apply Q
    for elem = 1:size(u,2)
      for i = 1:sbp.numnodes
        for j = 1:sbp.numnodes
          res[i,elem] += sbp.Q[i,j,di]*u[j,elem]
        end
      end
    end
  end
end

function weakdifferentiate!{T}(sbp::SBPOperator{T}, di::Int,
                               u::AbstractArray{T,3}, res::AbstractArray{T,3};
                               trans::Bool=false)
  @assert( sbp.numnodes == size(u,2) && sbp.numnodes == size(res,2) )
  @assert( length(u) == length(res) )
  @assert( di > 0 && di <= size(sbp.Q,3) )
  if trans # apply transposed Q
    for elem = 1:size(u,3)
      for i = 1:sbp.numnodes
        for j = 1:sbp.numnodes
          for field = 1:size(u,1)
            res[field,i,elem] += sbp.Q[j,i,di]*u[field,j,elem]
          end
        end
      end
    end
  else # apply Q
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
end

@doc """
### SummationByParts.differentiate!

Applies the SBP differentiation matrix operator, D, to data in `u` and **adds**
the result to `res`.  Different methods are available depending on the rank of
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
### SummationByParts.directionaldifferentiate

Performs a directional derivative (in reference space) at a given node.  In
constrast with other differentiate routines, the input field `u` is for **a
single element**, not a collection of elements.

**Inputs**

* `sbp`: an SBP operator type
* `dir`: a direction vector for the directional derivative
* `u`: the field that is being differentiated (either a scalar or vector)
* `i`: index of the node at which the derivative is desired

**Returns**

* `Ddir`: derivative of `u` in direction `dir`

"""->
function directionaldifferentiate{T}(sbp::SBPOperator{T}, dir::Array{T,1}, 
                                     u::AbstractArray{T,1}, i::Int)
  @assert( size(sbp.Q, 3) == size(dir,1) )
  Ddir = zero(T)
  for di = 1:size(sbp.Q, 3)
    tmp = zero(T)
    for j = 1:sbp.numnodes
      tmp += sbp.Q[i,j,di]*u[j]
    end
    Ddir += dir[di]*tmp/sbp.w[i]
  end
  return Ddir
end

function directionaldifferentiate{T}(sbp::SBPOperator{T}, dir::Array{T,1}, 
                                     u::AbstractArray{T,2}, i::Int)
  @assert( size(sbp.Q, 3) == size(dir,1) )
  Ddir = zeros(T, size(u,1))
  for di = 1:size(sbp.Q, 3)
    tmp = zeros(Ddir)
    for j = 1:sbp.numnodes
      tmp[:] += sbp.Q[i,j,di]*u[:,j]
    end
    Ddir[:] += dir[di]*tmp[:]/sbp.w[i]
  end
  return Ddir
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
### SummationByParts.boundaryintegrate!

Integrates a numerical flux over a boundary using appropriate mass matrices
defined on the element faces.  Different methods are available depending on the
rank of `u`:

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
* `bndryfaces`: list of boundary faces stored as an array of `Boundary`s
* `u`: the array of data that the boundary integral depends on
* `x` (optional): Cartesian coordinates stored in (coord,node,element) format
* `dξdx`: Jacobian of the element mapping (as output from `mappingjacobian!`)
* `bndryflux`: function to compute the numerical flux over the boundary

**In/Outs**

* `res`: where the result of the integration is stored

"""->
function boundaryintegrate!{T}(sbp::SBPOperator{T}, bndryfaces::Array{Boundary},
                               u::AbstractArray{T,2}, dξdx::AbstractArray{T,4},
                               bndryflux::Function, res::AbstractArray{T,2})
  @assert( sbp.numnodes == size(u,1) == size(res,1) == size(dξdx,3) )
  @assert( size(dξdx,4) == size(u,2) == size(res,2) )
  @assert( length(u) == length(res) )
  for bndry in bndryfaces
    for i = 1:sbp.numfacenodes
      # j = element-local index for ith node on face 
      j = sbp.facenodes[i, bndry.face]
      flux = bndryflux(u[j,bndry.element], dξdx[:,:,j,bndry.element],
                       sbp.facenormal[:,bndry.face])
      for i2 = 1:sbp.numfacenodes
        j2 = sbp.facenodes[i2, bndry.face]
        res[j2,bndry.element] += sbp.wface[i2,i]*flux
      end
    end
  end
end

function boundaryintegrate!{T}(sbp::SBPOperator{T}, bndryfaces::Array{Boundary},
                               u::AbstractArray{T,3}, dξdx::AbstractArray{T,4},
                               bndryflux::Function, res::AbstractArray{T,3})
  @assert( sbp.numnodes == size(u,2) == size(res,2) == size(dξdx,3) )
  @assert( size(dξdx,4) == size(u,3) == size(res,3) )
  @assert( length(u) == length(res) )
  flux = zeros(T, (size(u,1)))
  for bndry in bndryfaces
    for i = 1:sbp.numfacenodes
      # j = element-local index for ith node on face 
      j = sbp.facenodes[i, bndry.face]
      flux = bndryflux(u[:,j,bndry.element], dξdx[:,:,j,bndry.element],
                       sbp.facenormal[:,bndry.face])
      for i2 = 1:sbp.numfacenodes
        j2 = sbp.facenodes[i2, bndry.face]
        res[:,j2,bndry.element] += sbp.wface[i2,i]*flux
      end
    end
  end
end

function boundaryintegrate!{T}(sbp::SBPOperator{T}, bndryfaces::Array{Boundary},
                               u::AbstractArray{T,2}, x::AbstractArray{T,3},
                               dξdx::AbstractArray{T,4}, bndryflux::Function,
                               res::AbstractArray{T,2})
  @assert( sbp.numnodes == size(u,1) == size(res,1) == size(dξdx,3) == size(x,2) )
  @assert( size(dξdx,4) == size(u,2) == size(res,2) == size(x,3) )
  @assert( length(u) == length(res) )
  for bndry in bndryfaces
    for i = 1:sbp.numfacenodes
      # j = element-local index for ith node on face 
      j = sbp.facenodes[i, bndry.face]
      flux = bndryflux(u[j,bndry.element], x[:,j,bndry.element], 
                       dξdx[:,:,j,bndry.element], sbp.facenormal[:,bndry.face])
      for i2 = 1:sbp.numfacenodes
        j2 = sbp.facenodes[i2, bndry.face]
        res[j2,bndry.element] += sbp.wface[i2,i]*flux
      end
    end
  end
end

function boundaryintegrate!{T}(sbp::SBPOperator{T}, bndryfaces::Array{Boundary},
                               u::AbstractArray{T,3}, x::AbstractArray{T,3},
                               dξdx::AbstractArray{T,4}, bndryflux::Function,
                               res::AbstractArray{T,3})
  @assert( sbp.numnodes == size(u,2) == size(res,2) == size(dξdx,3) == size(x,2) )
  @assert( size(dξdx,4) == size(u,3) == size(res,3) == size(x,3) )
  @assert( length(u) == length(res) )
  flux = zeros(T, (size(u,1)))
  for bndry in bndryfaces
    for i = 1:sbp.numfacenodes
      # j = element-local index for ith node on face 
      j = sbp.facenodes[i, bndry.face]
      flux = bndryflux(u[:,j,bndry.element], x[:,j,bndry.element], 
                       dξdx[:,:,j,bndry.element], sbp.facenormal[:,bndry.face])
      for i2 = 1:sbp.numfacenodes
        j2 = sbp.facenodes[i2, bndry.face]
        res[:,j2,bndry.element] += sbp.wface[i2,i]*flux
      end
    end
  end
end

@doc """
### SummationByParts.mappingjacobian!

Evaluates the (scaled) Jacobian of the mapping from reference coordinates to
physical coordinates, as well as the determinant of the Jacobian.  The values
returned in dξdx are scaled by the determinant, so they have the same units as
the boundary measure (i.e. length in 2D, or length^2 in 3D).  This scaling is
adopted, because conservation laws written in conservative form in the reference
frame use the scaled Jacobian.

**Inputs**

* `sbp`: an SBP operator type
* `x`: the physical coordinates; 1st dim = coord, 2nd dim = node, 3rd dim = elem

**In/Outs**

* `dξdx`: the scaled Jacobian of the mapping; 1st dim = ref coord, 2nd dim =
  phys coord, 3rd dim = node, 3rd dim = elem
* `jac`: the determinant of the Jacobian; 1st dim = node, 2nd dim = elem

"""->
function mappingjacobian!{T}(sbp::TriSBP{T}, x::AbstractArray{T,3},
                             dξdx::AbstractArray{T,4}, jac::AbstractArray{T,2})
  @assert( sbp.numnodes == size(x,2) && sbp.numnodes == size(dξdx,3) )
  @assert( size(x,3) == size(dξdx,4) )
  @assert( size(x,1) == 2 && size(dξdx,1) == 2 && size(dξdx,2) == 2 )
  fill!(dξdx, zero(T))
  dxdξ = zeros(T, (2,sbp.numnodes,size(x,3)))
  # compute d(x,y)/dxi and set deta/dx and deta/dy
  differentiate!(sbp, 1, x, dxdξ)
  dξdx[2,1,:,:] = -dxdξ[2,:,:]
  dξdx[2,2,:,:] = dxdξ[1,:,:]
  # compute d(x,y)/deta and set dxi/dx and dxi/dy
  fill!(dxdξ, zero(T))
  differentiate!(sbp, 2, x, dxdξ)
  dξdx[1,2,:,:] = -dxdξ[1,:,:]
  dξdx[1,1,:,:] = dxdξ[2,:,:]
  # compute the determinant of the Jacobian
  for elem = 1:size(x,3)
    for i = 1:sbp.numnodes
      jac[i,elem] = 1.0/(dξdx[1,1,i,elem]*dξdx[2,2,i,elem] - 
      dξdx[1,2,i,elem]*dξdx[2,1,i,elem])
    end
  end
  # check for negative jac here?
end

function mappingjacobian!{T}(sbp::TetSBP{T}, x::AbstractArray{T,3},
                             dξdx::AbstractArray{T,4}, jac::AbstractArray{T,2})
  @assert( sbp.numnodes == size(x,2) && sbp.numnodes == size(dξdx,3) )
  @assert( size(x,3) == size(dξdx,4) )
  @assert( size(x,1) == 3 && size(dξdx,1) == 3 && size(dξdx,2) == 3 )
  numelem = size(x,3)
  dxdξ = zeros(T, (3,sbp.numnodes,numelem,3))
  # calculate the derivative of the coordinates with respect to (xi,eta,zeta)
  # using the SBP operator
  for di = 1:3
    differentiate!(sbp, di, x, sub(dxdξ,:,:,:,di)) 
  end
  fill!(dξdx, zero(T))
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
            dξdx[it1,di2,i,elem] += 
            coeff*(x[it22,j,elem]*dxdξ[it21,j,elem,it2] -
                   x[it21,j,elem]*dxdξ[it22,j,elem,it2])            
            dξdx[it2,di2,i,elem] += 
            coeff*(x[it21,j,elem]*dxdξ[it22,j,elem,it1] -
                   x[it22,j,elem]*dxdξ[it21,j,elem,it1])
          end
        end
      end
    end
  end
  # scale metrics by norm and calculate the determinant of the mapping
  for elem = 1:numelem
    for i = 1:sbp.numnodes
      dξdx[:,:,i,elem] /= sbp.w[i]
      jac[i,elem] = 1.0/(
                         dxdξ[1,i,elem,1]*dxdξ[2,i,elem,2]*dxdξ[3,i,elem,3] +
                         dxdξ[1,i,elem,2]*dxdξ[2,i,elem,3]*dxdξ[3,i,elem,1] +
                         dxdξ[1,i,elem,3]*dxdξ[2,i,elem,1]*dxdξ[3,i,elem,2] -
                         dxdξ[1,i,elem,1]*dxdξ[2,i,elem,3]*dxdξ[3,i,elem,2] -
                         dxdξ[1,i,elem,2]*dxdξ[2,i,elem,1]*dxdξ[3,i,elem,3] -
                         dxdξ[1,i,elem,3]*dxdξ[2,i,elem,2]*dxdξ[3,i,elem,1])
    end
  end
  # check for negative jac here?
end

@doc """
### SummationByParts.edgestabilize!

Adds edge-stabilization (see Burman and Hansbo, doi:10.1016/j.cma.2003.12.032)
to a given residual.  Different methods are available depending on the rank of
`u`:

* For *scalar* fields, it is assumed that `u` is a rank-2 array, with the first
dimension for the local-node index, and the second dimension for the element
index.
* For *vector* fields, `u` is a rank-3 array, with the first dimension for the
index of the vector field, the second dimension for the local-node index, and
the third dimension for the element index.

Naturally, the number of entries in the dimension of `u` (and `res`)
corresponding to the nodes must be equal to the number of nodes in the SBP
operator `sbp`.

**Inputs**

* `sbp`: an SBP operator type
* `ifaces`: list of element interfaces stored as an array of `Interface`s
* `u`: the array of solution data
* `x`: Cartesian coordinates stored in (coord,node,element) format
* `dξdx`: scaled Jacobian of the mapping (as output from `mappingjacobian!`)
* `jac`: determinant of the Jacobian
* `α`: array of transformation terms (see below)
* `stabscale`: function to compute the edge-stabilization scaling (see below)

**In/Outs**

* `res`: where the result of the integration is stored

**Details**

The array `α` is used to compute the directional derivative normal to the faces.
For a 2-dimensional problem, it can be computed as follows:

      for k = 1:mesh.numelem
        for i = 1:sbp.numnodes
          for di1 = 1:2
            for di2 = 1:2
              α[di1,di2,i,k] = (dξdx[di1,1,i,k].*dξdx[di2,1,i,k] + 
                                dξdx[di1,2,i,k].*dξdx[di2,2,i,k])*jac[i,k]
            end
          end
        end
      end

The function `stabscale` has the signature

  function stabscale(u, dξdx, nrm)

where `u` is the solution at a node, `dξdx` is the (scaled) Jacobian at the same
node, and `nrm` is the normal to the edge in reference space.  `stabscale`
should return the necessary scaling to ensure the edge-stabilization has the
desired asymptotic convergence rate.

"""->
function edgestabilize!{T}(sbp::SBPOperator{T}, ifaces::Array{Interface},
                           u::AbstractArray{T,2}, x::AbstractArray{T,3},
                           dξdx::AbstractArray{T,4}, jac::AbstractArray{T,2},
                           α::AbstractArray{T,4}, stabscale::Function,
                           res::AbstractArray{T,2})
  @assert( sbp.numnodes == size(u,1) == size(res,1) == size(dξdx,3) == size(x,2) 
          == size(α,3) )
  @assert( size(dξdx,4) == size(α,4) == size(u,2) == size(res,2) == size(x,3) )
  @assert( length(u) == length(res) )
  dim = size(sbp.Q, 3)
  Dn = zero(T)
  dirL = zeros(3)
  dirR = zeros(3)
  tmpL = zero(T)
  tmpR = zero(T)
  for face in ifaces
    EDn = zeros(T, (sbp.numfacenodes) )
    for i = 1:sbp.numfacenodes
      # iL = element-local index for ith node on left element face
      # iR = element-local index for ith node on right element face
      iL = sbp.facenodes[i, face.faceL]
      iR = sbp.facenodes[getnbrnodeindex(sbp, face, i), face.faceR]
      # apply the normal-derivative difference operator along the face
      dir = α[:,:,iL,face.elementL]*sbp.facenormal[:,face.faceL]
      Dn = directionaldifferentiate(sbp, dir, u[:,face.elementL], iL)
      dir = α[:,:,iR,face.elementR]*sbp.facenormal[:,face.faceR]
      Dn += directionaldifferentiate(sbp, dir, u[:,face.elementR], iR)
      # get differential area element: need 1/ds for each Dn term (here and loop
      # below)to get unit normal, and then need ds for integration, so net
      # result is 1/ds
      nrm = sbp.facenormal[:,face.faceL]
      dξ = dξdx[:,:,iL,face.elementL]
      ds = norm(nrm[1]*dξ[1,:] + nrm[2]*dξ[2,:])
      # apply the scaling function
      Dn *= stabscale(u[iL,face.elementL], dξdx[:,:,iL,face.elementL],
                      sbp.facenormal[:,face.faceL])./ds # note that u[iL] = u[iR]
      # add the face-mass matrix contribution
      for j = 1:sbp.numfacenodes
        EDn[j] += sbp.wface[j,i]*Dn
      end
    end
    # here we use hand-coded reverse-mode to apply the transposed
    # normal-derivative difference operator
    for i = 1:sbp.numfacenodes
      iL = sbp.facenodes[i, face.faceL]
      iR = sbp.facenodes[getnbrnodeindex(sbp, face, i), face.faceR]
      dirL = α[:,:,iL,face.elementL]*sbp.facenormal[:,face.faceL]
      dirR = α[:,:,iR,face.elementR]*sbp.facenormal[:,face.faceR]
      for di = 1:size(sbp.Q, 3)
        tmpL = dirL[di]*EDn[i]/sbp.w[iL]
        tmpR = dirR[di]*EDn[i]/sbp.w[iR]
        for j = 1:sbp.numnodes
          res[j,face.elementL] += sbp.Q[iL,j,di]*tmpL
          res[j,face.elementR] += sbp.Q[iR,j,di]*tmpR
        end
      end
    end
  end
end

function edgestabilize!{T}(sbp::SBPOperator{T}, ifaces::Array{Interface},
                           u::AbstractArray{T,3}, x::AbstractArray{T,3},
                           dξdx::AbstractArray{T,4}, jac::AbstractArray{T,2},
                           α::AbstractArray{T,4}, stabscale::Function,
                           res::AbstractArray{T,3})
  @assert( sbp.numnodes == size(u,2) == size(res,2) == size(dξdx,3) == size(x,2) 
          == size(α,3) )
  @assert( size(dξdx,4) == size(α,4) == size(u,3) == size(res,3) == size(x,3) )
  @assert( length(u) == length(res) )
  Dn = zeros(T, (size(u,1)) )
  tmpL = zeros(Dn)
  tmpR = zeros(Dn)
  for face in ifaces
    EDn = zeros(T, (size(u,1), sbp.numfacenodes) )
    for i = 1:sbp.numfacenodes
      # iL = element-local index for ith node on left element face
      # iR = element-local index for ith node on right element face
      iL = sbp.facenodes[i, face.faceL]
      iR = sbp.facenodes[getnbrnodeindex(sbp, face, i), face.faceR]
      # apply the normal-derivative difference operator along the face
      dir = α[:,:,iL,face.elementL]*sbp.facenormal[:,face.faceL]
      Dn[:] = directionaldifferentiate(sbp, dir, u[:,:,face.elementL], iL)
      dir = α[:,:,iR,face.elementR]*sbp.facenormal[:,face.faceR]
      Dn[:] += directionaldifferentiate(sbp, dir, u[:,:,face.elementR], iR)
      # get differential area element: need 1/ds for each Dn term (here and loop
      # below) to get unit normals, and then need ds for integration, so net
      # result is 1/ds
      nrm = sbp.facenormal[:,face.faceL]
      dξ = dξdx[:,:,iL,face.elementL]
      ds = norm(nrm[1]*dξ[1,:] + nrm[2]*dξ[2,:])
      # apply the scaling function
      Dn[:] *= stabscale(u[:,iL,face.elementL], dξdx[:,:,iL,face.elementL],
                         sbp.facenormal[:,face.faceL])./ds # note that u[iL] = u[iR]
      # add the face-mass matrix contribution
      for j = 1:sbp.numfacenodes
        EDn[:,j] += sbp.wface[j,i]*Dn[:]
      end
    end
    for i = 1:sbp.numfacenodes
      iL = sbp.facenodes[i, face.faceL]
      iR = sbp.facenodes[getnbrnodeindex(sbp, face, i), face.faceR]
      dirL = α[:,:,iL,face.elementL]*sbp.facenormal[:,face.faceL]
      dirR = α[:,:,iR,face.elementR]*sbp.facenormal[:,face.faceR]
      # here we use hand-coded reverse-mode to apply the transposed
      # normal-derivative difference operator
      for di = 1:size(sbp.Q, 3)
        tmpL[:] = dirL[di]*EDn[:,i]/sbp.w[iL]
        tmpR[:] = dirR[di]*EDn[:,i]/sbp.w[iR]
        for j = 1:sbp.numnodes
          res[:,j,face.elementL] += sbp.Q[iL,j,di]*tmpL[:]
          res[:,j,face.elementR] += sbp.Q[iR,j,di]*tmpR[:]
        end
      end
    end
  end
end
