# SBP abstract and concrete type definitions

"""
### SBP.AbstractSBP

`AbstractSBP` is a parametric abstract type that defines summation-by-parts
finite-difference operators.

"""
abstract type AbstractSBP{T<:Number} end

"""
### SBP.TriSBP

Defines diagonal-norm SBP first-derivative operators on a line segment.

**Fields**

* `degree` : maximum polynomial degree for which the derivatives are exact
* `numnodes` : number of nodes on line segment required for these operators
* `cub` : a symmetric cubature type for line segments (usually LG or LGL)
* `vtx` : vertices of the reference element in computational space
* `w` : cubature weights, i.e. the diagonal SBP norm, stored as an array
* `Q[:,:,1]` : discrete stiffness matrix operator
  """
struct LineSegSBP{T} <: AbstractSBP{T}
  degree::Int
  numnodes::Int
  cub::LineSymCub{T}
  vtx::Array{T,2}
  w::Array{T,1}
  Q::Array{T,3}
  E::Array{T,3}

  # inner constructor
  function LineSegSBP{T}(degree::Int, cub::LineSymCub{T}, vtx::Array{T,2},
                         w::Array{T,1}, Q::Array{T,3}, E::Array{T,3}) where T
    numnodes = cub.numnodes
    @assert( size(Q,1) == size(Q,2) == size(w,1) == numnodes )
    @assert( size(Q,3) == 1 )
    @assert( size(E,3) == 1 )
    new{T}(degree, numnodes, cub, vtx, w, Q, E)
  end
end


"""
### SBP.TriSBP

Defines diagonal-norm SBP first-derivative operators on a right-triangle.

**Fields**

* `degree` : maximum polynomial degree for which the derivatives are exact
* `numnodes` : number of nodes in the triangle required for these operators
* `cub` : a symmetric cubature type for triangles
* `vtx` : vertices of the reference element in computational space
* `w` : cubature weights, i.e. the diagonal SBP norm, stored as an array
* `Q[:,:,i]` : discrete stiffness matrix operator in ith coordinate direction

"""
struct TriSBP{T} <: AbstractSBP{T}
  degree::Int
  numnodes::Int
  cub::TriSymCub{T}
  vtx::Array{T,2}
  w::Array{T,1}
  Q::Array{T,3}
  E::Array{T,3}

  # inner constructor
  function TriSBP{T}(degree::Int, cub::TriSymCub{T}, vtx::Array{T,2},
                     w::Array{T,1}, Q::Array{T,3}, E::Array{T,3}) where T
    @assert( degree >= 1 && degree <= 10)
    numnodes = cub.numnodes
    @assert( size(Q,1) == size(Q,2) == size(w,1) == numnodes )
    @assert( size(Q,3) == 2 )
    @assert( size(E,3) == 2 )
    new{T}(degree, numnodes, cub, vtx, w, Q, E)
  end
end

"""
### SBP.SparseTriSBP

Defines diagonal-norm SBP first-derivative operators on a right-triangle using a
cubature rule that is greater than 2*p-1.  This provides additional flexiblity
in the SBP operator that is used to make a sparse S.

**Fields**

* `degree` : maximum polynomial degree for which the derivatives are exact
* `numnodes` : number of nodes in the triangle required for these operators
* `cub` : a symmetric cubature type for triangles
* `vtx` : vertices of the reference element in computational space
* `w` : cubature weights, i.e. the diagonal SBP norm, stored as an array
* `Q[:,:,i]` : discrete stiffness matrix operator in ith coordinate direction

"""
struct SparseTriSBP{T} <: AbstractSBP{T}
  degree::Int
  numnodes::Int
  cub::TriSymCub{T}
  vtx::Array{T,2}
  w::Array{T,1}
  Q::Array{T,3}

  function SparseTriSBP{T}(;degree::Int=1, faceorder::Array{Int,1}=[1;2;3], 
                           internal=false, cubdegree::Int=2*degree+1) where T
    @assert( degree >= 1 && degree <= 4 )
    if internal
      cub, vtx = getTriCubatureOmega(cubdegree, T)
    else
      cub, vtx = getTriCubatureGamma(cubdegree, T)
    end
    numnodes = cub.numnodes
    Q = zeros(T, (numnodes, numnodes, 2))
    w, Q = SummationByParts.buildsparseoperators(cub, vtx, degree)
    new{T}(degree, numnodes, cub, vtx, w, Q)
  end
end

"""
### SBP.TetSBP

Defines diagonal-norm SBP first-derivative operators on a right-tetrahedron.

**Fields**

* `degree` : maximum polynomial degree for which the derivatives are exact
* `numnodes` : number of nodes in the tetrahedron required for these operators
* `cub` : a symmetric cubature type for tetrahedra
* `vtx` : vertices of the reference element in computational space
* `w` : cubature weights, i.e. the diagonal SBP norm, stored as an array
* `Q[:,:,i]` : discrete stiffness matrix operator in ith coordinate direction

"""
struct TetSBP{T} <: AbstractSBP{T}
  degree::Int
  numnodes::Int
  cub::TetSymCub{T}
  vtx::Array{T,2}
  w::Array{T,1}
  Q::Array{T,3}
  E::Array{T,3}
  
  # inner constructor
  function TetSBP{T}(degree::Int, cub::TetSymCub{T}, vtx::Array{T,2},
                     w::Array{T,1}, Q::Array{T,3}, E::Array{T,3}) where T
    @assert( degree >= 1 && degree <= 5)
    numnodes = cub.numnodes
    @assert( size(Q,1) == size(Q,2) == size(w,1) == numnodes )
    @assert( size(Q,3) == 3 )
    @assert( size(E,3) == 3 )
    new{T}(degree, numnodes, cub, vtx, w, Q, E)
  end
end

"""
### SBP.SparseTetSBP

Defines diagonal-norm SBP first-derivative operators on a right-tetrahedron
using a cubature rule that is greater than 2*p-1.  This provides additional
flexiblity in the SBP operator that is used to make a sparse S.

**Fields**

* `degree` : maximum polynomial degree for which the derivatives are exact
* `numnodes` : number of nodes in the tetrahedron required for these operators
* `cub` : a symmetric cubature type for tetrahedra
* `vtx` : vertices of the reference element in computational space
* `w` : cubature weights, i.e. the diagonal SBP norm, stored as an array
* `Q[:,:,i]` : discrete stiffness matrix operator in ith coordinate direction

"""
struct SparseTetSBP{T} <: AbstractSBP{T}
  degree::Int
  numnodes::Int
  cub::TetSymCub{T}
  vtx::Array{T,2}
  w::Array{T,1}
  Q::Array{T,3}

  function SparseTetSBP{T}(;degree::Int=1, faceorder::Array{Int,1}=[1;2;3;4],
                           internal=false, cubdegree::Int=2*degree-1) where T
    @assert( degree >= 1 && degree <= 3 )
    if internal
      cub, vtx = getTetCubatureOmega(cubdegree, T)
    else
      cub, vtx = getTetCubatureGamma(cubdegree, T)
    end
    numnodes = cub.numnodes
    Q = zeros(T, (numnodes, numnodes, 3))
    w, Q = SummationByParts.buildsparseoperators(cub, vtx, degree)
    new{T}(degree, numnodes, cub, vtx, w, Q)
  end
end

"""
### SBP.AbstractFace

`AbstractFace` is a parametric abstract type that defines face-based data and
operations (e.g. volume-to-face reconstruction, face integration, etc) for
summation-by-parts finite-difference operators.

"""
abstract type AbstractFace{T<:Number} end

"""
### SBP.DenseFace

`DenseFace` is a parametric abstract type that defines face-based data and
operations (e.g. volume-to-face reconstruction, face integration, etc) for
summation-by-parts finite-difference operators.  This is a subtype for which
interpolation is a dense matrix.

"""
abstract type DenseFace{T} <: AbstractFace{T} end

"""
### SBP.LineSegFace

Defines a "face" between two LineSegSBP operators with the same cubature nodes.

**Fields**

* `degree` : face integration is exact for polys of degree 2*`degree`
* `numnodes` : number of cubature nodes (always 1)
* `stencilsize` : number of nodes in the reconstruction stencil
* `dstencilsize` : number of nodes in the derivative operator stencils
* `cub` : a symmetric cubature type for line-segment faces (i.e. points)
* `vtx` : the vertices of the face in reference space, [-1]
* `wface` : mass matrix (quadrature) for the face (always 1.0)
* `interp[:,:]` : volume-to-face-nodes reconstruction operator
* `perm[:,:]` : permutation for volume nodes so `interp` can be used on both sides
* `deriv[:,:]` : derivative operators for face-based coordinate system
* `dperm[:,:]` : permutation for volume nodes so `deriv` can be used on both sides
* `nbrperm[:,:]` : permutation for face nodes on neighbour element

"""
struct LineSegFace{T} <: DenseFace{T}
  degree::Int
  numnodes::Int
  stencilsize::Int
  dstencilsize::Int
  cub::PointSymCub{T}
  vtx::Array{T,2}
  wface::Array{T,1}
  normal::Array{T,2}
  interp::Array{T,2}
  perm::Array{Int,2}
  deriv::Array{T,3}
  dperm::Array{Int,2}
  nbrperm::Array{Int,2}

  # inner constructor
  function LineSegFace{T}(degree::Int, facecub::PointSymCub{T},
                          facevtx::Array{T,2}, interp::Array{T,2},
                          perm::Array{Int,2}, deriv::Array{T,3},
                          dperm::Array{Int,2}) where T
    @assert( degree >= 1 )
    numnodes = facecub.numnodes
    @assert( size(interp,2) == size(deriv,2) == numnodes )
    normal = T[-1; 1]'
    nbrperm = SymCubatures.getneighbourpermutation(facecub)
    wface = SymCubatures.calcweights(facecub)
    stencilsize = size(interp,1)
    dstencilsize = size(deriv,1)
    new{T}(degree, facecub.numnodes, stencilsize, dstencilsize, facecub, facevtx, 
        wface, normal, interp, perm, deriv, dperm, nbrperm)
  end
end

"""
### SBP.TriFace

Defines a face between two TriSBP operators with the same cubature nodes

**Fields**

* `degree` : face integration is exact for polys of degree 2*`degree`
* `numnodes` : number of cubature nodes
* `stencilsize` : number of nodes in the reconstruction stencil
* `dstencilsize` : number of nodes in the derivative operator stencils
* `cub` : a symmetric cubature type for triangle faces (i.e. edges)
* `vtx` : the vertices of the face in reference space, [-1,1]
* `wface` : mass matrix (quadrature) for the face
* `interp[:,:]` : volume-to-face-nodes reconstruction operator
* `perm[:,:]` : permutation for volume nodes so `interp` can be used on all sides
* `deriv[:,:]` : derivative operators for face-based coordinate system
* `dperm[:,:]` : permutation for volume nodes so `deriv` can be used on all sides
* `nbrperm[:,:]` : permutation for face nodes on neighbour element

"""
struct TriFace{T} <: DenseFace{T}
  degree::Int
  numnodes::Int
  stencilsize::Int
  dstencilsize::Int
  cub::LineSymCub{T}
  vtx::Array{T,2}
  wface::Array{T,1}
  normal::Array{T,2}
  interp::Array{T,2}
  perm::Array{Int,2}
  deriv::Array{T,3}
  dperm::Array{Int,2}
  nbrperm::Array{Int,2}

  # inner constructor
  function TriFace{T}(degree::Int, facecub::LineSymCub{T}, facevtx::Array{T,2},
                      interp::Array{T,2}, perm::Array{Int,2},
                      deriv::Array{T,3}, dperm::Array{Int,2}) where T
    @assert( degree >= 1 && degree <= 10 )
    numnodes = facecub.numnodes
    @assert( size(interp,2) == size(deriv,2) == numnodes )
    normal = T[0 -1; 1 1; -1 0]'
    #-----
    # normal = T[0 -1; sqrt(3)/2 1/2; -sqrt(3)/2 1/2]'
    # wface = SymCubatures.calcweights(facecub)./2.0
    #---
    nbrperm = SymCubatures.getneighbourpermutation(facecub)
    wface = SymCubatures.calcweights(facecub)
    stencilsize = size(interp,1)
    dstencilsize = size(deriv,1)
    new{T}(degree, facecub.numnodes, stencilsize, dstencilsize, facecub, facevtx, 
        wface, normal, interp, perm, deriv, dperm, nbrperm)
  end
end

"""
### SBP.TetFace

Defines a face between two TetSBP operators with the same cubature nodes

**Fields**

* `degree` : face integration is exact for polys of degree 2*`degree`
* `numnodes` : number of cubature nodes
* `stencilsize` : number of nodes in the reconstruction stencil
* `dstencilsize` : number of nodes in the derivative operator stencils
* `cub` : a symmetric cubature type for tetrahedral faces (i.e. triangles)
* `vtx` : the vertices of the face in the reference space of the face
* `wface` : mass matrix (quadrature) for the face
* `interp[:,:]` : volume-to-face-nodes reconstruction operator
* `perm[:,:]` : permutation for volume nodes so `interp` can be used on all sides
* `deriv[:,:]` : derivative operators for face-based coordinate system
* `dperm[:,:]` : permutation for volume nodes so `deriv` can be used on all sides
* `nbrperm[:,:]` : permutation for face nodes on neighbour element

"""
struct TetFace{T} <: DenseFace{T}
  degree::Int
  numnodes::Int
  stencilsize::Int
  #dstencilsize::Int
  cub::TriSymCub{T}
  vtx::Array{T,2}
  wface::Array{T,1}
  normal::Array{T,2}
  interp::Array{T,2}
  perm::Array{Int,2}
  #deriv::Array{T,3}
  #dperm::Array{Int,2}
  nbrperm::Array{Int,2}

  # inner constructor
  function TetFace{T}(degree::Int, facecub::TriSymCub{T}, facevtx::Array{T,2},
                      interp::Array{T,2}, perm::Array{Int,2}) where T
    @assert( degree >= 1 && degree <= 5 )
    numnodes = facecub.numnodes
    @assert( size(interp,2) == numnodes )
    normal = T[0 0 -1; 0 -1 0; 1 1 1; -1 0 0]'
    nbrperm = SymCubatures.getneighbourpermutation(facecub)
    wface = SymCubatures.calcweights(facecub)
    stencilsize = size(interp,1)
    new{T}(degree, numnodes, stencilsize, facecub, facevtx, wface, normal, interp,
        perm, nbrperm)
  end
end

"""
### SBP.SparseFace

`SparseFace` is a parametric abstract type that defines face-based data and
operations (e.g. volume-to-face reconstruction, face integration, etc) for
summation-by-parts finite-difference operators in the case where the
face-cubature nodes and volume nodes coincide (i.e. diagonal E operators).

"""
abstract type SparseFace{T} <: AbstractFace{T} end

"""
### SBP.TriSparseFace

Defines a face between two TriSBP operators with the same cubature nodes, in
which the face-cubature nodes and volume nodes coincide (i.e. diagonal E
operators).

**Fields**

* `degree` : face integration is exact for polys of degree 2*`degree`
* `numnodes` : number of cubature nodes
* `dstencilsize` : number of nodes in the derivative operator stencils
* `cub` : a symmetric cubature type for triangle faces (i.e. edges)
* `vtx` : the vertices of the face in reference space, [-1,1]
* `wface` : mass matrix (quadrature) for the face
* `perm[:,:]` : maps volume nodes to face nodes on each side
* `deriv[:,:]` : derivative operators for face-based coordinate system
* `dperm[:,:]` : permutation for volume nodes so `deriv` can be used on all sides
* `nbrperm[:,:]` : permutation for face nodes on neighbour element

"""
struct TriSparseFace{T} <: SparseFace{T}
  degree::Int
  numnodes::Int
  dstencilsize::Int
  cub::LineSymCub{T}
  vtx::Array{T,2}
  wface::Array{T,1}
  normal::Array{T,2}
  perm::Array{Int,2}
  deriv::Array{T,3}
  dperm::Array{Int,2}
  nbrperm::Array{Int,2}

  # inner constructor
  function TriSparseFace{T}(degree::Int, facecub::LineSymCub{T},
                            facevtx::Array{T,2}, perm::Array{Int,2},
                            deriv::Array{T,3}, dperm::Array{Int,2}) where T
    # @assert( degree >= 1 && degree <= 5 )
    numnodes = facecub.numnodes
    @assert( size(deriv,2) == numnodes )
    normal = T[0 -1; 1 1; -1 0]'
    nbrperm = SymCubatures.getneighbourpermutation(facecub)
    wface = SymCubatures.calcweights(facecub)
    dstencilsize = size(deriv,1)
    new{T}(degree, facecub.numnodes, dstencilsize, facecub, facevtx, wface, normal,
        perm, deriv, dperm, nbrperm)
  end
end

"""
### SBP.TetSparseFace

Defines a face between two TetSBP operators with the same cubature nodes, in
which the face-cubature nodes and volume nodes coincide (i.e. diagonal E
operators).

**Fields**

* `degree` : face integration is exact for polys of degree 2*`degree`
* `numnodes` : number of cubature nodes
* `dstencilsize` : number of nodes in the derivative operator stencils
* `cub` : a symmetric cubature type for tetrahedral faces (i.e. triangles)
* `vtx` : the vertices of the face in the reference space of the face
* `wface` : mass matrix (quadrature) for the face
* `perm[:,:]` : permutation for volume nodes to face nodes on each side
* `deriv[:,:]` : derivative operators for face-based coordinate system
* `dperm[:,:]` : permutation for volume nodes so `deriv` can be used on all sides
* `nbrperm[:,:]` : permutation for face nodes on neighbour element

"""
struct TetSparseFace{T} <: SparseFace{T}
  degree::Int
  numnodes::Int
  #dstencilsize::Int
  cub::TriSymCub{T}
  vtx::Array{T,2}
  wface::Array{T,1}
  normal::Array{T,2}
  perm::Array{Int,2}
  #deriv::Array{T,3}
  #dperm::Array{Int,2}
  nbrperm::Array{Int,2}

  # inner constructor
  function TetSparseFace{T}(degree::Int, facecub::TriSymCub{T},
                            facevtx::Array{T,2}, perm::Array{Int,2}) where T
    @assert( degree >= 1 && degree <= 5 )
    numnodes = facecub.numnodes
    normal = T[0 0 -1; 0 -1 0; 1 1 1; -1 0 0]'
    nbrperm = SymCubatures.getneighbourpermutation(facecub)
    wface = SymCubatures.calcweights(facecub)
    new{T}(degree, numnodes, facecub, facevtx, wface, normal, perm, nbrperm)
  end
end
