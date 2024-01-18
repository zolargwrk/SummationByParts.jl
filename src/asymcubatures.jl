module AsymCubatures
using LinearAlgebra
# types and methods for mapping between asymmetry nodes for cubatures
# on various domains

export AsymCub, LineAsymCub, TriAsymCub, TetAsymCub

"""
### AsymCubatures.AsymCub

`AsymCub` is an abstract type that defines cubatures for asymmetric
nodal distributions.  It is parameterized on `T` in order to allow for future
implementations of arbitrary precision types.  The parameterization also permits
the use of the complex-step method for verification.

""" abstract type AsymCub{T<:Number} end

"""
### AsymCubatures.LineAsymCub
  
Defines an asymmetric quadrature rule on the interval [-1,1].  Current choices are uniform,
Legendre-Gauss-Lobatto (LGL), or Legendre-Gauss (LG) rules.

**Fields**

* `numparams` : total number of nodal degrees of freedom
* `numweights` : total number of unique weights
* `numnodes` : total number of nodes
* `numedgenodes` : number of edge nodes along the added dimension
* `params` : the actual values of the orbit nodal locations
* `weights` : values of the unique weights

"""
mutable struct LineAsymCub{T} <: AsymCub{T}
  numparams::Int
  numweights::Int
  numnodes::Int
  numedgenodes::Int
  edgeparams::Array{T,1}
  params::Array{T,1}
  weights::Array{T,1}
    
  function LineAsymCub{T}(;numedgenodes::Int=0) where T
    @assert(numedgenodes >= 0)
    numnodes = numedgenodes
    # compute the number of degrees of freedom and unique weights
    numparams = numnodes
    numweights = numnodes
    # initialize parameter arrays
    @assert(numweights == numedgenodes)
    edgeparams = zeros(T,numedgenodes)
    params = zeros(T, numparams)
    weights = zeros(T, numweights)
    new{T}(numparams, numweights, numnodes, numedgenodes, edgeparams, params, weights)
  end
end

"""
### AsymCubatures.TriAsymCub

Used to define asymmetric cubature rules on the triangle.  The `params` array
determines the position of the nodes, and the `weights` array
determines the value of the weight for each symmetric orbit.  

**Fields**

* `numparams` : total number of nodal degrees of freedom
* `numweights` : total number of unique weights
* `numnodes` : total number of nodes
* `numedgenodes` : number of edge nodes along the added dimension
* `params` : the actual values of the orbit nodal locations
* `weights` : values of the unique weights

"""
mutable struct TriAsymCub{T} <: AsymCub{T}
  numparams::Int
  numweights::Int
  numnodes::Int
  numedgenodes::Int
  edgeparams::Array{T,1}
  params::Array{T,1}
  weights::Array{T,1}

  function TriAsymCub{T}(;numedgenodes::Int=0) where T
    @assert(numedgenodes >= 0)
    # calculate total number of nodes
    numnodes = convert(Int, (numedgenodes+1)*numedgenodes/2.0)
    # compute the number of degrees of freedom and weights
    numparams = numnodes*2
    numweights = numnodes
    # initialize parameter arrays
    edgeparams = zeros(T,numedgenodes)
    params = zeros(T, numparams)
    weights = zeros(T, numweights)
    new{T}(numparams, numweights, numnodes, numedgenodes, edgeparams, params, weights)
  end
end

"""
### AsymCubatures.TetAsymCub

Used to define asymmetric cubature rules on the tetrahedron.  The `params` array
determines the position of the nodes, and the `weights` array
determines the value of the weight for each symmetric orbit. 

**Fields**

* `numparams` : total number of nodal degrees of freedom
* `numweights` : total number of unique weights
* `numnodes` : total number of nodes
* `numedgenodes` : number of edge nodes along the added dimension
* `params` : the actual values of the orbit nodal parameters
* `weights` : values of the unique weights

"""
mutable struct TetAsymCub{T} <: AsymCub{T}
  numparams::Int
  numweights::Int
  numnodes::Int
  numedgenodes::Int
  edgeparams::Array{T,1}
  params::Array{T,1}
  weights::Array{T,1}

  function TetAsymCub{T}(;numedgenodes::Int=0) where T
    @assert(numedgenodes >= 0)
    # compute the total number of nodes
    numnodes = convert(Int, numedgenodes*(numedgenodes+1)*(numedgenodes+2)/6.0)
    # compute the number of degrees of freedom and weights
    numparams = numnodes*3
    numweights = numnodes
    # initialize parameter arrays
    edgeparams = zeros(T,numedgenodes)
    params = zeros(T, numparams)
    weights = zeros(T, numweights)
    new{T}(numparams, numweights, numnodes, numedgenodes, edgeparams, params, weights)
  end
end

"""
### AsymCubatures.getnumboundarynodes

Returns the number of (explicit) boundary nodes

*Notes*: if the parameter value for an internal orbit is such that the
corresponding node lies on the boundary, this node is **NOT** included in the
boundary-node count returned.

**Inputs**

* `cub`: symmetric cubature rule

**Outputs**

* `numboundary`: number of boundary nodes

"""

function getnumboundarynodes(cub::LineAsymCub{T}) where {T}
  return 2
end

function getnumboundarynodes(cub::TriAsymCub{T}) where {T}
  numboundary = cub.numedgenodes*3
  return numboundary
end

function getnumboundarynodes(cub::TetAsymCub{T}) where {T}
  numboundary = convert(Int, cub.numedgenodes*(cub.numedgenodes+1)/2.0)*4
  return numboundary
end

"""
### SymCubatures.getnumfacenodes

Returns the number of nodes on an individual face of the element.

**Inputs**

* `cub`: asymmetric cubature rule

**Outputs**

* `numfacenodes`: number of nodes on a face

"""

function getnumfacenodes(cub::LineAsymCub{T}) where {T}
  return 1
end

function getnumfacenodes(cub::TriAsymCub{T}) where {T}
  numfacenodes = cub.numedgenodes
  return numfacenodes
end

function getnumfacenodes(cub::TetAsymCub{T}) where {T}
  numfacenodes = convert(Int, cub.numedgenodes*(cub.numedgenodes+1)/2.0)
  return numfacenodes
end

"""
### AsymCubatures.setparams!

Sets the parameter vector with a uniform nodal locations in the element given
a nodal location on an edge of the element.

**Inputs**

* `edgeparams` : the nodal distribution on the edge of the element

**In/Outs**

* `cub`: asymmetric cubature rule whose nodal parameters are being updated

"""

function setparams!(cub::LineAsymCub{T}, edgeparams::Array{T}) where {T}
    @assert( length(edgeparams) == cub.numedgenodes)
    cub.edgeparams = sort(vec(edgeparams))
    cub.params = sort(vec(edgeparams))
end

function setparams!(cub::TriAsymCub{T}, edgeparams::Array{T}) where {T}
    @assert( length(edgeparams) == cub.numedgenodes)
    n = length(edgeparams)
    xedge = sort(vec(edgeparams))
    xy = zeros(2,cub.numnodes)
    idx = 1
    for i = 1:n
        for j = 1:n-(i-1)
            xy[1,idx] = xedge[j]
            xy[2,idx] = xedge[i]
            idx += 1
        end
    end
    cub.edgeparams = sort(vec(edgeparams))
    cub.params = vec(xy')
end

function setparams!(cub::TetAsymCub{T}, edgeparams::Array{T}) where {T}
    @assert( length(edgeparams) == cub.numedgenodes)
    n = length(edgeparams)
    xedge = sort(vec(edgeparams))
    xyz = zeros(3,cub.numnodes)
    idx = 1
    for i = 1:n
        for j = 1:n-(i-1)
            for k = 1:n-(j-1)-(i-1)
                xyz[1,idx] = xedge[k]
                xyz[2,idx] = xedge[j]
                xyz[3,idx] = xedge[i]
                idx += 1
            end
        end
    end
    cub.edgeparams = sort(vec(edgeparams))
    cub.params = vec(xyz')
end

"""
### AsymCubatures.setweights!

Sets a cubature's (unique) weights.

**Inputs**

* `weights`: cubature weights 

**In/Outs**

* `cub`: asymmetric cubature rule whose weights are being updated

"""
function setweights!(cub::AsymCub{T}, weights::Array{T}) where {T}
  @assert( length(weights) == cub.numweights )
  cub.weights = vec(weights)
end

"""
### AsymCubatures.calcnodes

Returns the cubature nodes 

**Inputs**

* `cub`: symmetric cubature rule
**Outputs**

* `x`: cubature's nodal coordinates 

"""
function calcnodes(cub::LineAsymCub) 
  @assert(cub.numparams >= 0)
  @assert(cub.numnodes >= 1)
  x = cub.params
  return x
end

function calcnodes(cub::TriAsymCub) 
  @assert(cub.numparams >= 0)
  @assert(cub.numnodes >= 1)
  x = reshape(cub.params, (cub.numnodes,2))'
  return x
end

function calcnodes(cub::TetAsymCub) 
  @assert(cub.numparams >= 0)
  @assert(cub.numnodes >= 1)
  x = reshape(cub.params, (cub.numnodes,3))'
  return x
end

"""
### SymCubatures.calcweights

Map the unique cubature weights to the weights of all nodes.

**Inputs**

* `cub`: symmetric cubature rule

**Outputs**

* `w`: cubature's weights at all nodes

"""

function calcweights(cub::LineAsymCub{T}) where {T}
  @assert(cub.numweights >= 0)
  @assert(cub.numnodes >= 1)
  w = cub.weights
  return w
end

function calcweights(cub::TriAsymCub{T}) where {T}
  @assert(cub.numweights >= 0)
  @assert(cub.numnodes >= 1)

  w = cub.weights
  return w
end

function calcweights(cub::TetAsymCub{T}) where {T}
  @assert(cub.numweights >= 0)
  @assert(cub.numnodes >= 1)

  w = cub.weights
  return w
end

"""
### AsymCubatures.calcjacobianofnodes

Returns the Jacobian of the nodes with respect the nodes themselves.

*Notes*: Jac stores all the x-coordinate Jacobians first, then y (then z)

**Inputs**

* `cub`: symmetric cubature rule
* `vtx`: vertices of the cubature domain

**Outputs**

* `Jac`: Jacobian of the mapping from node parameters to nodes

"""

function calcjacobianofnodes(cub::TriAsymCub{T}) where {T}
  @assert(cub.numparams >= 0)
  @assert(cub.numnodes >= 1)

  Jac = zeros(T, (2*cub.numnodes, cub.numparams) )
  Jac[1:cub.numnodes, 1:cub.numnodes] = I(cub.numnodes)
  Jac[cub.numnodes+1:2*cub.numnodes, 1:cub.numnodes] = I(cub.numnodes)
  return Jac
end

function calcjacobianofnodes(cub::TetAsymCub{T}) where {T}
  @assert(cub.numparams >= 0)
  @assert(cub.numnodes >= 1)
  
  Jac = zeros(T, (3*cub.numnodes, cub.numparams) )
  Jac[1:cub.numnodes, 1:cub.numnodes] = I(cub.numnodes)
  Jac[cub.numnodes+1:2*cub.numnodes, 1:cub.numnodes] = I(cub.numnodes)
  Jac[2*cub.numnodes+1:3*cub.numnodes, 1:cub.numnodes] = I(cub.numnodes)
  return Jac
end

"""
### AsymCubatures.calcjacobianofweights

Returns the Jacobian of the nodal weights with respect to the unique weights.

**Inputs**

* `cub`: symmetric cubature rule

**Outputs**

* `Jac`: Jacobian of the mapping from (unique) weights to nodal weights

"""
function calcjacobianofweights(cub::LineAsymCub{T}) where {T}
  @assert(cub.numweights >= 0)
  @assert(cub.numnodes >= 1)

  Jac = zeros(T, (cub.numnodes, cub.numweights) )
  Jac[1:cub.numnodes,1:cub.numweights] = I(cub.numnodes)
  return Jac
end

function calcjacobianofweights(cub::TriAsymCub{T}) where {T}
  @assert(cub.numweights >= 0)
  @assert(cub.numnodes >= 1)

  Jac = zeros(T, (cub.numnodes, cub.numweights) )
  Jac[1:cub.numnodes,1:cub.numweights] = I(cub.numnodes)
  return Jac
end

function calcjacobianofweights(cub::TetAsymCub{T}) where {T}
  @assert(cub.numweights >= 0)
  @assert(cub.numnodes >= 1)

  Jac = zeros(T, (cub.numnodes, cub.numweights) )
  Jac[1:cub.numnodes,1:cub.numweights] = I(cub.numnodes)
  return Jac
end

"""
### AsymCubatures.calcjacobian

Returns the Jacobian of the nodal coordinates and weights with respect to their
parameters.  

**Inputs**

* `cub`: symmetric cubature rule
* `vtx`: vertices of the cubature domain

**Outputs**

* `Jac`: Jacobian of the nodal coordinates and weights

"""
function calcjacobian(cub::TriAsymCub{T}) where {T}
  Jac = zeros(T, (3*cub.numnodes, cub.numparams + cub.numweights) )
  Jac[1:2*cub.numnodes, 1:cub.numparams] = AsymCubatures.calcjacobianofnodes(cub)
  Jac[2*cub.numnodes+1:end, cub.numparams+1:end] =AsymCubatures.calcjacobianofweights(cub)
  return Jac
end

function calcjacobian(cub::TetAsymCub{T}) where {T}
  Jac = zeros(T, (4*cub.numnodes, cub.numparams + cub.numweights) )
  Jac[1:3*cub.numnodes, 1:cub.numparams] = AsymCubatures.calcjacobianofnodes(cub)
  Jac[3*cub.numnodes+1:end, cub.numparams+1:end] = AsymCubatures.calcjacobianofweights(cub)
  return Jac
end

"""
### AsymCubatures.getfacenodeindices

Returns the indices of the nodes that lie on each face.  See getbndrynodeindices
for a method that returns a single array of boundary nodes.

**Inputs**

* `cub`: a symmetric cubature rule whose boundary-node indices are sought

**Outputs**

* `bndryindices`: indicies of nodes that lie on boundary; there is a separate
  column of indices for each edge/face.

"""
function getfacenodeindices(cub::TriAsymCub{T}) where {T}
  # get the number of nodes on one edge
  numedge = getnumfacenodes(cub)
  bndryindices = zeros(Int, (numedge,3) )
  # assumes a right triangle with vtx = [-1 -1; 1 -1; -1 1] 
  xy = AsymCubatures.calcnodes(cub)'
  x = xy[:,1]
  y = xy[:,2]
  # facet 1
  j = 1
  for i = 1:cub.numnodes
    if (abs(y[i]+1.0)<1e-13)
      bndryindices[j,1] = i
      j+=1
    end
  end
  # facet 2
  j = 1
  for i = 1:cub.numnodes
    if (abs(x[i] + y[i]-0.0) < 1e-13)
      bndryindices[j,2] = i
      j+=1
    end
  end
  # facet 3
  j = numedge
  for i = 1:cub.numnodes
    if (abs(x[i]+1.0)<1e-13)
      bndryindices[j,3] = i
      j-=1
    end
  end

  return bndryindices
end

"""
### AsymCubatures.getfacebasedpermutation

Returns a permutation of the volume nodes (or a subset of them) for each face,
such that the same face operator can be applied to all faces.  This is useful
for volume-to-face interpolation or differentiation.

**Inputs**

* `cub`: an asymmetric cubature rule for which a face-based permutation is sought
* `faceonly`: if true, only face nodes are used in the permutation.

**Outputs**

* `perm`: permutation of the volume nodes for each face

"""
function getfacebasedpermutation(cub::TriAsymCub{T}; faceonly::Bool=false) where {T}
  if faceonly
    perm = zeros(Int, (getnumfacenodes(cub), 3))
    perm = getfacenodeindices(cub)
  else
    error("Not implemented.")
  end
  return perm
end

end