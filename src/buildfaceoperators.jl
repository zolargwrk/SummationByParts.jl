# This file gathers together functions used to build the SBP face operators
using LinearAlgebra

"""
### SummationByParts.buildfacereconstruction

Builds a matrix operator that can reconstruct a field from a set of volume nodes
to a set of face nodes.  The reconstruction operator is only constructed for one
face, but a permutation array is returned that allows the same operator to be
applied on all the faces.

**Inputs**

* `facecub`: symmetric cubature rule for the face
* `cub`: symmetric cubature rule for the volume
* `vtx`: vertices of the right simplex
* `d`: maximum total degree for the Proriol polynomials

**Returns**

* `R`: the volume- to face-node reconstruction operator
* `perm`: a permutation array that allows `R` to be applied on all faces

"""
function buildfacereconstruction(facecub::PointSymCub{T}, cub::LineSymCub{T},
                                    vtx::Array{T,2}, d::Int) where {T}
  # first, decide whether or not to use volume nodes or just face nodes
  if SymCubatures.getnumfacenodes(cub) == 1
    perm = SymCubatures.getfacebasedpermutation(cub, faceonly=true)
  else
    perm = SymCubatures.getfacebasedpermutation(cub, faceonly=false)
  end
  # evaluate the basis at the volume and face cubature points
  N = convert(Int, d+1 )
  Pv = zeros(T, (size(perm,1),N) )  
  Pf = zeros(T, (facecub.numnodes,N) ) 
  xv = SymCubatures.calcnodes(cub, vtx)
  xf = SymCubatures.calcnodes(facecub, vtx[[1],:])
  ptr = 1
  for i = 0:d
    Pv[:,ptr] = OrthoPoly.jacobipoly(vec(xv[1,perm[:,1]]), 0.0, 0.0, i)
    Pf[:,ptr] = OrthoPoly.jacobipoly(vec(xf[1,:]), 0.0, 0.0, i)
    ptr += 1
  end
  if SymCubatures.getnumfacenodes(cub) == 1
    A = kron(Pv',I(facecub.numnodes))
    b = vec(Pf)
    R = zeros(facecub.numnodes*size(perm,1))
    SummationByParts.calcSparseSolution!(A, b, R)
    R = reshape(R, (facecub.numnodes,size(perm,1)))
  else
    R = Matrix((pinv(Pv')*Pf')')
  end
  return R, perm
end

function buildfacereconstruction(facecub::LineSymCub{T}, cub::TriSymCub{T},
                                    vtx::Array{T,2}, d::Int; faceopertype::Symbol=:Omega) where {T}
  # first, decide whether or not to use volume nodes or just face nodes
  if SymCubatures.getnumfacenodes(cub) >= (d+1)
    perm = SymCubatures.getfacebasedpermutation(cub, faceonly=true)
  else
    perm = SymCubatures.getfacebasedpermutation(cub, faceonly=false)
  end
  # evaluate the basis at the volume and face cubature points
  N = convert(Int, (d+1)*(d+2)/2 )
  Pv = zeros(T, (size(perm,1),N) )  
  Pf = zeros(T, (facecub.numnodes,N) ) 
  xv = SymCubatures.calcnodes(cub, vtx)
  xf = SymCubatures.calcnodes(facecub, vtx[[1;2],:])
  ptr = 1
  for r = 0:d
    for j = 0:r
      i = r-j
      Pv[:,ptr] = OrthoPoly.proriolpoly(vec(xv[1,perm[:,1]]),
                                        vec(xv[2,perm[:,1]]), i, j)
      Pf[:,ptr] = OrthoPoly.proriolpoly(vec(xf[1,:]), vec(xf[2,:]), i, j)
      ptr += 1
    end
  end
  #R = Pf/Pv
  #R = (pinv(Pv')*Pf')'
  # notes on sparse solution:
  # A = kron(eye(size(perm,1)),Pv)
  # b = vec(Pf)
  # use basis pursuit on A x = b
  # find dominant entries and resolve for R
  #----------
  # below is commented for test, uncomment if needed
  # if SymCubatures.getnumfacenodes(cub) >= (d+1)
  #   A = kron(Pv', I(facecub.numnodes))
  #   #println("size(A) = ",size(A))
  #   #println("rank(A) = ",rank(A))
  #   b = vec(Pf)
  #   #println("size(b) = ",size(b))
  #   R = zeros(facecub.numnodes*size(perm,1))
  #   SummationByParts.calcSparseSolution!(A, b, R)
  #   R = reshape(R, (facecub.numnodes,size(perm,1)))
  # else
  #   R = (pinv(Pv')*Pf')'
  # end
  #----------
  if faceopertype == :DiagE
    R = 1.0*I(facecub.numnodes)
  else
    R = (pinv(Pv')*Pf')'
  end
  #----------

  return R, perm
end

function buildfacereconstruction(facecub::TriSymCub{T}, cub::TetSymCub{T},
                                    vtx::Array{T,2}, d::Int; faceopertype::Symbol=:Omega) where {T}
  # first, decide whether or not to use volume nodes or just face nodes
  if SymCubatures.getnumfacenodes(cub) >= (d+1)
    perm = SymCubatures.getfacebasedpermutation(cub, faceonly=true)
  else
    perm = SymCubatures.getfacebasedpermutation(cub, faceonly=false)
  end
  # evaluate the basis at the volume and face cubature points
  N = convert(Int, (d+1)*(d+2)*(d+3)/6 )
  Pv = zeros(T, (size(perm,1),N) )  
  Pf = zeros(T, (facecub.numnodes,N) ) 
  xv = SymCubatures.calcnodes(cub, vtx)
  xf = SymCubatures.calcnodes(facecub, vtx[[1;2;3],:])
  ptr = 1
  for r = 0:d
    for k = 0:r
      for j = 0:r-k
        i = r-j-k
        Pv[:,ptr] = OrthoPoly.proriolpoly(vec(xv[1,perm[:,1]]),
                                          vec(xv[2,perm[:,1]]),
                                          vec(xv[3,perm[:,1]]),
                                          i, j, k)
        Pf[:,ptr] = OrthoPoly.proriolpoly(vec(xf[1,:]), vec(xf[2,:]),
                                          vec(xf[3,:]), i, j, k)
        ptr += 1
      end
    end
  end
  #R = Pf/Pv
  # R = (pinv(Pv')*Pf')'
  if faceopertype == :DiagE
    R = 1.0*I(facecub.numnodes)
  else
    R = (pinv(Pv')*Pf')'
  end

  return R, perm
end

"""
### SummationByParts.buildfacederivatives

Builds matrix operators that can differentiate a polynomial field of degree `d`
from a set of volume nodes to a set of face nodes.  The derivative operators are
only constructed for one face, but a permutation array is returned that allows
the same operators to be applied on all the faces.

**Note**: the derivative operators are for the tangential and normal directions,
  and do not necessarily correspond to the directions ξ and η (and ζ) used for
  the volume derivatives.  These face derivatives are intended for edge
  stabilization.

**Inputs**

* `facecub`: symmetric cubature rule for the face
* `cub`: symmetric cubature rule for the volume
* `vtx`: vertices of the right simplex
* `d`: maximum total degree for the Proriol polynomials

**Returns**

* `D`: derivative operators in [face node, vol node, direction] format
* `perm`: a permutation array that allows `D` to be applied on all faces

"""
function buildfacederivatives(facecub::PointSymCub{T}, cub::LineSymCub{T},
                                 vtx::Array{T,2}, d::Int) where {T}
  perm = SymCubatures.getfacebasedpermutation(cub)
  # evaluate the basis at the volume and face cubature points
  N = d+1
  Pv = zeros(T, (cub.numnodes,N) )  
  dPdx = zeros(T, (facecub.numnodes,N) )
  xv = SymCubatures.calcnodes(cub, vtx)
  xf = SymCubatures.calcnodes(facecub, vtx[[1;],:])
  ptr = 1
  for i = 0:d
    Pv[:,ptr] = OrthoPoly.jacobipoly(vec(xv[1,perm[:,1]]), 0.0, 0.0, i)
    dPdx[:,ptr] = OrthoPoly.diffjacobipoly(vec(xf[1,:]), 0.0, 0.0, i)
    ptr += 1
  end
  A = pinv(Pv')
  D = zeros(cub.numnodes, facecub.numnodes, 1)
  D[:,:,1] = A*dPdx'
  return D, perm
end

function buildfacederivatives(facecub::LineSymCub{T}, cub::TriSymCub{T},
                                 vtx::Array{T,2}, d::Int) where {T}
  perm = SymCubatures.getfacebasedpermutation(cub)
  # evaluate the basis at the volume and face cubature points
  N = convert(Int, (d+1)*(d+2)/2 )
  Pv = zeros(T, (cub.numnodes,N) )  
  dPdx = zeros(T, (facecub.numnodes,N) )
  dPdy = zeros(T, (facecub.numnodes,N) )
  xv = SymCubatures.calcnodes(cub, vtx)
  xf = SymCubatures.calcnodes(facecub, vtx[[1;2],:])
  ptr = 1
  for r = 0:d
    for j = 0:r
      i = r-j
      Pv[:,ptr] = OrthoPoly.proriolpoly(vec(xv[1,perm[:,1]]),
                                        vec(xv[2,perm[:,1]]), i, j)
      dPdx[:,ptr], dPdy[:,ptr] = 
      OrthoPoly.diffproriolpoly(vec(xf[1,:]), vec(xf[2,:]), i, j)
      ptr += 1
    end
  end
  A = pinv(Pv')
  D = zeros(cub.numnodes, facecub.numnodes, 2)
  D[:,:,1] = A*dPdx'
  D[:,:,2] = A*dPdy'
  return D, perm
end
