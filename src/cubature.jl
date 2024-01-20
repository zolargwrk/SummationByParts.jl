module Cubature
# routines for constructing cubatures

using LinearAlgebra
using LeastSquaresOptim
using ..OrthoPoly
using ..SymCubatures, ..AsymCubatures
using ..Optimizer

export pointCubature
export quadrature, quadratureUniform
export getTriCubatureGamma, getTriCubatureOmega, getTriCubatureDiagE, getTriCubatureForTetFaceDiagE
export getTetCubatureGamma, getTetCubatureOmega, getTetCubatureDiagE

"""
### Cubature.cubatureresidual

This method computes the residuals, `F`, between a cubature, defined by `cub`,
and the true value of an integral.  Each residual corresponds with an orthogonal
polynomial on the simplex up to degree `q`.  The Jacobian, `dF`, of the
residual, with respect to the quadrature nodes and weights, is also returned.

**Inputs**

* `cub`: defines the nodes and weights of the cubature via symmetry orbits
* `q`: maximum degree of the othogonal polynomials used in the conditions
* `calc_grad`: indicates whether to calculate the gradients the orthogonal polynomials

**Outputs**

* `F`: the accuracy conditions for orthogonal polynomials up to degree q
* `dF`: derivative of F with respect to x, y, (z,) w, in that order

"""
function cubatureresidual(cub::LineSymCub{T}, q::Int; compute_grad=true) where {T}
  # compute the nodes and weights defined by cub
  vtx = reshape(T[-1; 1], (2,1))
  x = SymCubatures.calcnodes(cub, vtx)
  w = SymCubatures.calcweights(cub)
  num_eq = convert(Int, q+1)
  # loop over orthogonal polynomials of degree r <= q and form conditions
  F = zeros(T, (num_eq) )
  F[1] = -2.0/sqrt(2.0)
  dF = zeros(T, (num_eq, 2*cub.numnodes) )
  ptr = 1
  for r = 0:q
    P = OrthoPoly.jacobipoly(vec(x[1,:]), 0.0, 0.0, r)
    F[ptr] += (w'*P)[1]
    if compute_grad
      dPdx = OrthoPoly.diffjacobipoly(vec(x[1,:]), 0.0, 0.0, r)
      dF[ptr,:] = [w.*dPdx P]
    end
    ptr += 1
  end
  return F, dF 
end

# function cubatureresidual(cub::TriSymCub{T}, q::Int; compute_grad=true) where {T}
#   # compute the nodes and weights defined by cub
#   vtx = T[-1 -1; 1 -1; -1 1]
#   x = SymCubatures.calcnodes(cub, vtx)
#   w = SymCubatures.calcweights(cub)
#   num_eq = convert(Int, (q+1)*(q+2)/2)

#   # loop over orthogonal polynomials of degree r <= q and form conditions
#   F = zeros(T, (num_eq) )
#   F[1] = -2.0/sqrt(2.0)
#   dF = zeros(T, (num_eq, 3*cub.numnodes) )
#   ptr = 1
#   for r = 0:q
#     for j = 0:r
#       i = r-j
#       P = OrthoPoly.proriolpoly(vec(x[1,:]), vec(x[2,:]), i, j)
#       F[ptr] += (w'*P)[1]
#       if compute_grad
#         dPdx, dPdy = OrthoPoly.diffproriolpoly(vec(x[1,:]), vec(x[2,:]), i, j)
#         dF[ptr,:] = [w.*dPdx w.*dPdy P]
#       end
#       ptr += 1
#       #print("(i,j,k) = (",i,",",j,",",k,"): i+j+k = ",i+j+k,"\n")
#     end
#   end
#   return F, dF 
# end

function cubatureresidual(cub::TriSymCub{T}, q::Int; compute_grad=true) where {T}
  # compute the nodes and weights defined by cub
  vtx = T[-1 -1; 1 -1; -1 1]
  x = SymCubatures.calcnodes(cub, vtx)
  w = SymCubatures.calcweights(cub)
  num_eq = convert(Int, (q+1)*(q+2)/2)

  # loop over orthogonal polynomials of degree r <= q and form conditions
  F = zeros(T, (num_eq) )
  F[1] = -2.0/sqrt(2.0)
  dF = zeros(T, (num_eq, 3*cub.numnodes) )

  V, Vdx, Vdy = OrthoPoly.vandermonde(q, x[1,:], x[2,:], compute_grad=compute_grad)
  F = V'*w + F

  # P = Optimizer.preconditioner(Matrix(V'))
  # F = (P*V')*w + P*F

  # P = Optimizer.preconditioner_shannon(V)
  # F = (P'*V')*w + P'*F

  if compute_grad
    dF = [Vdx'*Diagonal(w) Vdy'*Diagonal(w) V']
    # dF = [P*Vdx'*Diagonal(w) P*Vdy'*Diagonal(w) P*V']
    # dF = [P'*Vdx'*Diagonal(w) P'*Vdy'*Diagonal(w) P'*V']
  end
  return F, dF 
end

# function cubatureresidual(cub::TriSymCub{T}, q::Int; compute_grad=true) where {T}
#   # compute the nodes and weights defined by cub
#   # vtx = T[-1 -1; 1 -1; -1 1]
#   vtx = T[0 0; 1 0; 0 1]
#   x = SymCubatures.calcnodes(cub, vtx)
#   w = SymCubatures.calcweights(cub)
#   num_eq = convert(Int, (q+1)*(q+2)/2)

#   # loop over orthogonal polynomials of degree r <= q and form conditions
#   F = zeros(T, (num_eq) )
#   # F[1] = -2.0/sqrt(2.0)
#   dF = zeros(T, (num_eq, 3*cub.numnodes) )

#   V, Vdx, Vdy, Vinteg = OrthoPoly.vandermonde_monomial(q, x[1,:], x[2,:], compute_grad=compute_grad, compute_integ=true)
#   F = V'*w - Vinteg
#   # P = Optimizer.preconditioner(Matrix(V'))
#   # F = (P*V')*w - P*Vinteg

#   if compute_grad
#     dF = [Vdx'*Diagonal(w) Vdy'*Diagonal(w) V']
#     # dF = [P*Vdx'*Diagonal(w) P*Vdy'*Diagonal(w) P*V']
#   end
#   return F, dF 
# end

# function cubatureresidual(cub::TetSymCub{T}, q::Int; compute_grad=true) where {T}
#   # compute the nodes and weights defined by cub
#   vtx = T[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
#   x = SymCubatures.calcnodes(cub, vtx)
#   w = SymCubatures.calcweights(cub)
#   num_eq = convert(Int, (q+1)*(q+2)*(q+3)/6)
#   # loop over orthogonal polynomials of degree r <= q and form conditions
#   F = zeros(T, (num_eq) )
#   F[1] = -2.0/sqrt(3.0)
#   dF = zeros(T, (num_eq, 4*cub.numnodes) )
#   ptr = 1
#   for r = 0:q
#     for k = 0:r
#       for j = 0:r-k
#         i = r-j-k
#         P = OrthoPoly.proriolpoly(vec(x[1,:]), vec(x[2,:]), vec(x[3,:]), i, j, k)
#         F[ptr] += (w'*P)[1]
#         if compute_grad
#           dPdx, dPdy, dPdz = OrthoPoly.diffproriolpoly(vec(x[1,:]), vec(x[2,:]),vec(x[3,:]), i, j, k)
#           dF[ptr,:] = [w.*dPdx w.*dPdy w.*dPdz P]
#         end
#         ptr += 1
#         #print("(i,j,k) = (",i,",",j,",",k,"): i+j+k = ",i+j+k,"\n")
#       end
#     end
#   end
#   return F, dF 
# end

function cubatureresidual(cub::TetSymCub{T}, q::Int; compute_grad=true) where {T}
  # compute the nodes and weights defined by cub
  vtx = T[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
  # vtx = T[0 0 0; 1 0 0; 0 1 0; 0 0 1]
  x = SymCubatures.calcnodes(cub, vtx)
  w = SymCubatures.calcweights(cub) 
  num_eq = convert(Int, (q+1)*(q+2)*(q+3)/6)
  # loop over orthogonal polynomials of degree r <= q and form conditions
  F = zeros(T, (num_eq) )
  F[1] = -2.0/sqrt(3.0)
  dF = zeros(T, (num_eq, 4*cub.numnodes) )

  V, Vdx, Vdy, Vdz = OrthoPoly.vandermonde(q, x[1,:], x[2,:], x[3,:], compute_grad=compute_grad)
  F = V'*w + F
  # F = V'*w - Vinteg
  # P = Optimizer.preconditioner(Matrix(V'))
  # F = (P*V')*w + P*F

  if compute_grad
    dF = [Vdx'*Diagonal(w) Vdy'*Diagonal(w) Vdz'*Diagonal(w) V']
    # dF = [P*Vdx'*Diagonal(w) P*Vdy'*Diagonal(w) P*V']
  end
  return F, dF 
end

"""
### Cubature.solvecubature!{SymCub{T}}
  
Attempts to solve for the nodes and weights of a cubature that is exact for
polynomials of degree r <= `q`.  The nodes and weights of the cubature are
defined by `cub`, which is a parametric abstract type (see symcubatures.jl).

**Inputs**

* `cub`: symmetric cubature rule
* `q`: maximum (desired) degree for which the cubature is exact
* `mask`: array of indicies of parameters and weights that are free
* `tol`: tolerance with which to solve the accuracy conditions
* `hist`: if true, print the residual-norm convergence history
* `xinit`: initial parameter guess
* `delta1`: purturbation constant that determines how much the parameters are perturbed
            in the case of stagnation with fmin <= 0.1
* `delta2`: purturbation constant that determines how much the parameters are perturbed
            in the case of stagnation with fmin > 0.1

**In/Outs**

* `cub`: on entry, defines the initial guess for the cubature nodes and weights.
  on exit, defines the nodes and weights that satisfy the desired accuracy.

"""
function solvecubature!(cub::SymCub{T}, q::Int, mask::AbstractArray{Int64,1};
  tol=10*eps(typeof(real(one(T)))),hist::Bool=false, verbose::Bool=false, xinit=[],
  delta1::Float64=1e-2, delta2::Float64=1e-2) where {T}
  @assert( length(mask) <= cub.numparams + cub.numweights )

  # Particle Swarm Optimization 
  n = cub.numparams + cub.numweights
  v = zeros(T, n)
  v[1:cub.numparams] = cub.params
  v[cub.numparams+1:end] = cub.weights
  # xinit=[]
  if xinit==[]
    xinit = copy(v)
  end
  nperturb = 0
  k=0
  v1 = copy(v)
  v2 = copy(v)
  
  if xinit!=[]
    SymCubatures.setparams!(cub, xinit[1:cub.numparams])
    SymCubatures.setweights!(cub, xinit[cub.numparams+1:end])
    xinit[1:cub.numparams] = cub.params
    xinit[cub.numparams+1:end] = cub.weights
  end

  nperturb_all = 0
  iter_pso = 0
  iter_lma = 0
  fmin,xinit,k = Optimizer.levenberg_marquardt(Cubature.cubatureresidual,cub,q,mask, xinit=xinit, maxiter=100, tol=tol, nu=100.0, verbose=1)
  iter_lma+=k
  # v = xinit
  # println(fmin)
  # println(xinit)

  fbest = 1.0
  xbest = copy(v)
  # verbose=true
  for i = 1:5000
    fmin1,v1,_,_,k_pso,nperturb = Optimizer.pso(Cubature.cubatureresidual, n, cub, q, mask,xinit=xinit, np=10,maxiter=800, 
                                            tol=tol, delta1=delta1, delta2=delta2, save_iter=false,verbose=verbose)
    fmin2,v2,k_lma = Optimizer.levenberg_marquardt(Cubature.cubatureresidual,cub,q,mask, xinit=v1, xL=1e-8, xR=2.0, maxiter=60, tol=tol, nu=1000.0, verbose=verbose)

    iter_pso += k_pso
    iter_lma += k_lma
    nperturb_all += nperturb

    # if (fmin2 < fmin1)
    #   xinit = v2
    # else
    #   xinit = v1
    # end
    use_v2 = false
    fmin3,_,_,_,_,_ = Optimizer.pso(Cubature.cubatureresidual, n, cub, q, mask,xinit=v2, np=2, maxiter=10, 
                                            tol=tol, delta1=delta1, delta2=delta2, save_iter=false,verbose=false)
    if (fmin3 < fmin1)
      use_v2 = true
    end

    if use_v2
      xinit = v2
    else
      xinit = v1
    end

    if (fbest < 1.0 && fmin2-fbest > 0.10 && fmin1-fbest > 0.10)
      xinit = xbest
    elseif use_v2 #(fmin2 < fmin1)
      xbest = v2
      fbest = fmin2
    else
      xbest = v1
      fbest = fmin1
    end

    if (fmin2<5e-14||fmin1<5e-14) && (minimum(v2)>0||minimum(v1)>0)
      if verbose
        println(v2)
      end
      break
    end
    if verbose
      println(v2)
    end
  end
  v = v2

  SymCubatures.setparams!(cub, v[1:cub.numparams])
  SymCubatures.setweights!(cub, v[cub.numparams+1:end])
  F, _ = Cubature.cubatureresidual(cub, q)
  res = norm(F)
  # hist ? print("\titer ",k, ":  nrestart ",nrestart,":  res norm = ",res,"\n") : nothing
  hist ? println("----------------------------------------------------------------------------------------------") : nothing
  hist ? println("iter_pso = ",iter_pso, ":  iter_lma = ",iter_lma, ":  nperturb_pso = ",nperturb_all,":  res norm = ",res) : nothing
  hist ? println("----------------------------------------------------------------------------------------------") : nothing
  # hist ? print("res norm = ",res,"\n") : nothing
  if res < tol
    return res
  end
  return res
end

# function solvecubaturelma!(cub::SymCub{T}, q::Int, mask::AbstractArray{Int64,1};
#   tol=10*eps(typeof(real(one(T)))),hist::Bool=false, verbose::Bool=false, xinit=[], maxiter::Int=100,
#   delta1::Float64=1e-2, delta2::Float64=1e-2) where {T}
#   @assert( length(mask) <= cub.numparams + cub.numweights )

#   # Particle Swarm Optimization 
#   n = cub.numparams + cub.numweights
#   v = zeros(T, n)
#   v[1:cub.numparams] = cub.params
#   v[cub.numparams+1:end] = cub.weights
#   # xinit=[]
#   if xinit==[]
#     xinit = copy(v)
#   end
#   nperturb = 0
#   k=0
#   v1 = copy(v)
#   v2 = copy(v)
  
#   if xinit!=[]
#     SymCubatures.setparams!(cub, xinit[1:cub.numparams])
#     SymCubatures.setweights!(cub, xinit[cub.numparams+1:end])
#     xinit[1:cub.numparams] = cub.params
#     xinit[cub.numparams+1:end] = cub.weights
#   end

#   nperturb_all = 0
#   iter_pso = 0
#   iter_lma = 0
#   fmin,xinit,k = Optimizer.levenberg_marquardt(Cubature.cubatureresidual,cub,q,mask, xinit=xinit, maxiter=100, tol=tol, nu=100.0, verbose=1)
#   iter_lma+=k
#   # v = xinit
#   # println(fmin)
#   # println(xinit)

#   fbest = 1.0
#   xbest = copy(v)
#   mask_weight=Int[]
#   append!(mask_weight, 1:cub.numweights)
#   # verbose=true
#   for i = 1:2 #1:5000
#     fmin1,vw,_,_,k_pso,nperturb = Optimizer.pso_weight(Cubature.cubatureresidual, cub.numweights, cub, q, mask_weight,
#                                                        xinit=xinit[cub.numparams+1:cub.numparams+cub.numweights],np=10,maxiter=800, 
#                                                        tol=tol, delta1=delta1, delta2=delta2, save_iter=false,verbose=verbose)
#     v1 = convert.(T,collect(Iterators.flatten([xinit[1:cub.numparams],vw])))
#     fmin2,v2,k_lma = Optimizer.levenberg_marquardt(Cubature.cubatureresidual,cub,q,mask, xinit=v1, xL=1e-8, xR=2.0, maxiter=30, tol=tol, nu=100000.0, verbose=verbose)

#     iter_pso += k_pso
#     iter_lma += k_lma
#     nperturb_all += nperturb

#     # if (fmin2 < fmin1)
#     #   xinit = v2
#     # else
#     #   xinit = v1
#     # end
#     use_v2 = false
#     fmin3,_,_,_,_,_ = Optimizer.pso(Cubature.cubatureresidual, n, cub, q, mask,xinit=v2, np=2, maxiter=10, 
#                                             tol=tol, delta1=delta1, delta2=delta2, save_iter=false,verbose=false)
#     if (fmin3 < fmin1)
#       use_v2 = true
#     end

#     if use_v2
#       xinit = v2
#     else
#       xinit = v1
#     end

#     if (fbest < 1.0 && fmin2-fbest > 0.10 && fmin1-fbest > 0.10)
#       xinit = xbest
#     elseif use_v2 #(fmin2 < fmin1)
#       xbest = v2
#       fbest = fmin2
#     else
#       xbest = v1
#       fbest = fmin1
#     end

#     if (fmin2<5e-14||fmin1<5e-14) && (minimum(v2)>0||minimum(v1)>0)
#       if verbose
#         println(v2)
#       end
#       break
#     end
#     if verbose
#       println(v2)
#     end
#   end
#   v = v2

#   SymCubatures.setparams!(cub, v[1:cub.numparams])
#   SymCubatures.setweights!(cub, v[cub.numparams+1:end])
#   F, _ = Cubature.cubatureresidual(cub, q)
#   res = norm(F)
#   # hist ? print("\titer ",k, ":  nrestart ",nrestart,":  res norm = ",res,"\n") : nothing
#   hist ? println("----------------------------------------------------------------------------------------------") : nothing
#   hist ? println("iter_pso = ",iter_pso, ":  iter_lma = ",iter_lma, ":  nperturb_pso = ",nperturb_all,":  res norm = ",res) : nothing
#   hist ? println("----------------------------------------------------------------------------------------------") : nothing
#   # hist ? print("res norm = ",res,"\n") : nothing
#   if res < tol
#     return res
#   end
#   return res
# end

function solvecubaturelma!(cub::SymCub{T}, q::Int, mask::AbstractArray{Int64,1};
  tol=10*eps(typeof(real(one(T)))),hist::Bool=false, maxiter::Int=100, verbose::Bool=false, xinit=[],
  delta1::Float64=1e-2, delta2::Float64=1e-2) where {T}
  @assert( length(mask) <= cub.numparams + cub.numweights )

  # Particle Swarm Optimization 
  n = cub.numparams + cub.numweights
  v = zeros(T, n)
  v[1:cub.numparams] = cub.params
  v[cub.numparams+1:end] = cub.weights
  # xinit=[]
  if xinit==[]
    xinit = copy(v)
  end
  nperturb = 0
  k=0
  # v1 = copy(v)
  # v2 = copy(v)
  
  if xinit!=[]
    SymCubatures.setparams!(cub, xinit[1:cub.numparams])
    SymCubatures.setweights!(cub, xinit[cub.numparams+1:end])
    xinit[1:cub.numparams] = cub.params
    xinit[cub.numparams+1:end] = cub.weights
  end

  nperturb_all = 0
  iter_pso = 0
  iter_lma = 0
  fmin,v,k = Optimizer.levenberg_marquardt(Cubature.cubatureresidual,cub,q,mask, xL=0.0, xR=2.0, xinit=xinit, maxiter=maxiter, tol=tol, nu=1e3, verbose=verbose)
  iter_lma+=k

  if (fmin<5e-14) && (minimum(v)>0.0)
    if verbose
      println(v)
    end
  end

  SymCubatures.setparams!(cub, v[1:cub.numparams])
  SymCubatures.setweights!(cub, v[cub.numparams+1:end])
  F, _ = Cubature.cubatureresidual(cub, q)
  res = norm(F)
  # hist ? print("\titer ",k, ":  nrestart ",nrestart,":  res norm = ",res,"\n") : nothing
  hist ? println("----------------------------------------------------------------------------------------------") : nothing
  hist ? println("iter_lma = ",iter_lma, ":  res norm = ",res) : nothing
  hist ? println("----------------------------------------------------------------------------------------------") : nothing
  # hist ? print("res norm = ",res,"\n") : nothing
  return res

end

"""
### Cubature.solvecubatureweights!{SymCub{T}}
  
Attempts to solve for the weights of a cubature that is exact for
polynomials of degree r <= `q`.  The weights (and nodes) of the cubature are
defined by `cub`, which is a parametric abstract type (see symcubatures.jl).

**Inputs**

* `q`: maximum (desired) degree for which the cubature is exact
* `tol`: tolerance with which to solve the accuracy conditions
* `hist`: if true, print the residual-norm convergence history

**In/Outs**

* `cub`: on entry, defines the initial guess for the cubature nodes and weights.
  on exit, defines the nodes and weights that satisfy the desired accuracy.

"""
function solvecubatureweights!(cub::SymCub{T}, q::Int;
                                  tol=10*eps(typeof(real(one(T)))),
                                  hist::Bool=false) where {T}
  Jac = SymCubatures.calcjacobianofweights(cub)

  # compute accuracy for initial guess 
  F, dF = Cubature.cubatureresidual(cub, q)
  res = norm(F)
  res0 = res
  res_old = res
  hist ? print("solvecubatureweights!:\n") : nothing
  hist ? print("\titer ",0,": res norm = ",res,"\n") : nothing
  if (res < tol)
    return
  end
  
  # Levenberg–Marquardt loop
  maxiter = 2000
  nu = 1000.0 #100.0
  v = zeros(T, (cub.numweights) )
  v = cub.weights
  for k = 1:maxiter
    JtJ = Jac'*(dF[:,end-cub.numnodes+1:end]'*
      dF[:,end-cub.numnodes+1:end])*Jac
    H = JtJ + nu*diagm(diag(JtJ))
    g = -Jac'*dF[:,end-cub.numnodes+1:end]'*F
    dv = H\g

    # update cubature definition and check for convergence
    v += dv
    SymCubatures.setweights!(cub, v)
    F, dF = Cubature.cubatureresidual(cub, q)
    res = norm(F)
    hist ? print("\titer ",k,": res norm = ",res,"\n") : nothing
    if res < tol
      #println("size(JtJ) = ",size(JtJ))
      #println("rank(JtJ) = ",rank(JtJ))
      return
    end

    # trust-region like update
    if res > res_old
      v -= dv
      SymCubatures.setweights!(cub, v)
      F, dF = Cubature.cubatureresidual(cub, q)
      nu *= 4.0
    else
      nu /= 2.0
      res_old = res
    end

  end
  error("solvecubatureweights failed to find solution in ",maxiter," iterations")
end

function solvecubatureweights!(cub::SymCub{T}, q::Int, mask::AbstractArray{Int64,1};
                              tol=10*eps(typeof(real(one(T)))),hist::Bool=false, verbose::Bool=false, xinit=[],
                              delta1::Float64=1e-2, delta2::Float64=1e-2) where {T}
  
  # Particle Swarm Optimization 
  n = cub.numweights
  v = cub.weights

  if xinit==[]
    xinit = copy(v)
  end
  nperturb = 0
  k=0
  v1 = copy(v)
  v2 = copy(v)
  
  if xinit!=[]
    SymCubatures.setweights!(cub, xinit)
    xinit = cub.weights
  end

  nperturb_all = 0
  iter_pso = 0
  iter_lma = 0
  fmin,xinit,k = Optimizer.levenberg_marquardt_weight(Cubature.cubatureresidual,cub,q,mask, xinit=xinit, maxiter=200, tol=tol, nu=10.0, verbose=0)
  iter_lma+=k
  # println(fmin)
  # println(xinit)

  fbest = 1.0
  xbest = copy(v)
  # verbose=true
  for i = 1:5000
    fmin1,v1,_,_,k_pso,nperturb = Optimizer.pso_weight(Cubature.cubatureresidual, n, cub, q, mask,xinit=xinit, np=10,maxiter=400, 
                                            tol=tol, delta1=delta1, delta2=delta2, save_iter=false,verbose=verbose)
    fmin2,v2,k_lma = Optimizer.levenberg_marquardt_weight(Cubature.cubatureresidual,cub,q,mask, xinit=v1, maxiter=200, tol=tol, nu=1000.0, verbose=verbose)

    iter_pso += k_pso
    iter_lma += k_lma
    nperturb_all += nperturb

    # if (fmin2 < fmin1)
    #   xinit = v2
    # else
    #   xinit = v1
    # end
    use_v2 = false
    fmin3,_,_,_,_,_ = Optimizer.pso_weight(Cubature.cubatureresidual, n, cub, q, mask,xinit=v2, np=10, maxiter=100, 
                                            tol=tol, delta1=delta1, delta2=delta2, save_iter=false,verbose=false)
    if (fmin3 < fmin1)
      use_v2 = true
    end

    if use_v2
      xinit = v2
    else
      xinit = v1
    end

    if (fbest < 1.0 && fmin2-fbest > 0.10 && fmin1-fbest > 0.10)
      xinit = xbest
    elseif use_v2 #(fmin2 < fmin1)
      xbest = v2
      fbest = fmin2
    else
      xbest = v1
      fbest = fmin1
    end

    if (fmin2<1e-12||fmin1<1e-12) && (minimum(v2)>0||minimum(v1)>0)
      if verbose
        println(v2)
      end
      break
    end
    if verbose
      println(v2)
    end
  end
  v = v2

  SymCubatures.setweights!(cub, v)
  F, _ = Cubature.cubatureresidual(cub, q)
  res = norm(F)
  # hist ? print("\titer ",k, ":  nrestart ",nrestart,":  res norm = ",res,"\n") : nothing
  hist ? println("----------------------------------------------------------------------------------------------") : nothing
  hist ? println("iter_pso = ",iter_pso, ":  iter_lma = ",iter_lma, ":  nperturb_pso = ",nperturb_all,":  res norm = ",res) : nothing
  hist ? println("----------------------------------------------------------------------------------------------") : nothing
  # hist ? print("res norm = ",res,"\n") : nothing
  if res < tol
    return
  end

  # Jac = SymCubatures.calcjacobianofweights(cub)

  # # compute accuracy for initial guess 
  # F, dF = fun(cub, q)
  # res = norm(F)
  # res0 = res
  # res_old = res
  # hist ? print("solvecubatureweights!:\n") : nothing
  # hist ? print("\titer ",0,": res norm = ",res,"\n") : nothing
  # if (res < tol)
  #   return
  # end
  
  # # Levenberg–Marquardt loop
  # alpha = 1.0
  # maxiter = 200
  # nu = 1000.0 #100.0
  # v = zeros(T, (cub.numweights) )
  # v = cub.weights
  # dv = zeros(size(v))
  # for k = 1:maxiter
  #   J = dF[:,end-cub.numnodes+1:end]*Jac
  #   JtJ = J'*J
  #   H = JtJ + nu*diagm(diag(JtJ))
  #   g = -J'*F
  #   # dv = H\g
  #   fill!(dv, zero(T))
  #   Hred = H[mask,mask]
  #   dv[mask] = pinv(Hred,1e-14)*g[mask]

  #   # update cubature definition and check for convergence
  #   v += dv
  #   SymCubatures.setweights!(cub, v)
  #   F, dF = fun(cub, q)
  #   res = norm(F)
  #   hist ? print("\titer ",k,": res norm = ",res,"\n") : nothing
  #   if res < tol
  #     #println("size(JtJ) = ",size(JtJ))
  #     #println("rank(JtJ) = ",rank(JtJ))
  #     return
  #   end

  #   # trust-region like update
  #   if res > res_old
  #     v -= dv
  #     SymCubatures.setweights!(cub, v)
  #     F, dF = fun(cub, q)
  #     nu *= 4.0
  #   else
  #     nu /= 2.0
  #     res_old = res
  #   end

  #   eps = 1e-4 
  #   for i=1:length(axes(v,1))
  #     if v[i]<1e-10
  #         v -= alpha*dv
  #         alpha = (eps - v[i])/(dv[i])
  #         v += alpha*dv
  #     end
  #   end
  #   alpha = 1.0

  # end
  # error("solvecubatureweights failed to find solution in ",maxiter," iterations")
end


"""
### Cubature.pointCubature

This returns a (trivial) point cubature and default vertex -1

**Inputs**

* `T`: the data type used to represent the cubature

**Outputs**

* `cub`: a symmetric cubature for point
* `vtx`: vertex, [-1]

"""
function pointCubature(T::Type=Float64)
  pt = PointSymCub{T}()
  SymCubatures.setweights!(pt, T[1;])
  vtx = reshape(T[-1;], (1,1))
  return pt, vtx
end

"""
### Cubature.quadrature{T}

This high-level function computes and returns a symmetric cubature of requested
accuracy on the interval [-1,1]

**Inputs**

* `q`: maximum degree of polynomial for which the cubature is exact
* `T`: the data type used to represent the cubature
* `internal`: if true, all nodes are strictly internal (default false)
* `tol`: tolerance with which to solve the cubature

**Outputs**

* `cub`: a symmetric cubature for the interval [-1,1]
* `vtx`: vertices, [-1,1]

"""
function quadrature(q::Int, T=Float64; internal::Bool=false)
  if internal
    # all nodes are internal (LG quadrature)
    N = div(q+2,2)
    x, w = OrthoPoly.lgnodes(N, T)
    alpha = zeros(T, (div(N,2)))
    weight = zeros(T, (div(N+1,2)))
    for i = 1:div(N,2)
      alpha[i] = (1 + x[i])/2
      weight[i] = w[i]
    end
    if rem(N,2) == 0
      centroid = false
    else
      centroid = true
      weight[div(N+1,2)] = w[div(N+1,2)]
    end
    quad = LineSymCub{T}(vertices=false, centroid=centroid,
                         numedge=div(N,2))
    SymCubatures.setparams!(quad, alpha)
    SymCubatures.setweights!(quad, weight)
  else
    # vertices are included (LGL quadrature)
    N = div(q+4,2)
    x, w = OrthoPoly.lglnodes(N-1, T)
    alpha = zeros(T, (div(N-2,2)))
    weight = zeros(T, (div(N+1,2)))
    weight[1] = w[1]
    for i = 1:div(N-2,2)
      alpha[i] = (1 - x[i+1])/2
      weight[i+1] = w[i+1]
    end
    if rem(N,2) == 0 
      centroid=false
    else
      centroid=true
      weight[div(N+1,2)] = w[div(N+1,2)]
    end
    quad = LineSymCub{T}(vertices=true, centroid=centroid,
                         numedge=div(N-2,2))
    SymCubatures.setparams!(quad, alpha)
    SymCubatures.setweights!(quad, weight)
  end
  vtx = reshape(T[-1; 1], (2,1))
  return quad, vtx
end

"""
### Cubature.quadratureUniform{T}

This high-level function computes and returns a uniform cubature of requested
accuracy on the interval [-1,1]

**Inputs**

* `q`: maximum degree of polynomial for which the cubature is exact
* `N`: number of nodes (N >= q+1)
* `T`: the data type used to represent the cubature

**Outputs**

* `cub`: a symmetric cubature for the interval [-1,1]
* `vtx`: vertices, [-1,1]

"""
function quadratureUniform(q::Int, N::Int, T=Float64; internal::Bool=false)
  @assert(N >= q+1)

  if rem(N,2) == 0 
    centroid=false
  else
    centroid=true
  end
  numedge = div(N-2,2)
  quad = LineSymCub{T}(vertices=true, centroid=centroid,
                       numedge=numedge)
  alpha = zeros(T, (numedge))
  dx = 1.0/(N-1)
  for i = 1:numedge
    alpha[i] = i*dx 
  end    
  SymCubatures.setparams!(quad, alpha)
  weight = (2.0/N)*ones(T, (numedge + 1 + centroid))
  SymCubatures.setweights!(quad, weight)
  solvecubatureweights!(quad, q)
  #mask = zeros(Int64, (0))
  #append!(mask, (quad.numparams+1):(quad.numparams+quad.numweights))
  #Cubature.solvecubature!(quad, q, mask, tol=1e-15)

  vtx = reshape(T[-1; 1], (2,1))
  return quad, vtx
end

"""
### Cubature.quadratureGregory{T}

This high-level function computes and returns a uniform Gregory-type quadrature 
rules of requested accuracy on the interval [-1,1]

**Inputs**

* `q`: maximum degree of polynomial for which the cubature is exact
* `N`: number of nodes (N >= q+1)
* `T`: the data type used to represent the cubature

**Outputs**

* `cub`: a symmetric cubature for the interval [-1,1]
* `vtx`: vertices, [-1,1]

"""
function getLineCubatureGregory(q::Int, N::Int; T=Float64)
  @assert(N >= maximum([2*(q-1),2]))

  if rem(N,2) == 0 
    centroid=false
  else
    centroid=true
  end
  numedge = div(N-2,2)
  cub = LineSymCub{T}(vertices=true, centroid=centroid, numedge=numedge)
  alpha = zeros(T, (numedge))
  dx = 1.0/(N-1)
  for i = 1:numedge
    alpha[i] = i*dx 
  end    
  SymCubatures.setparams!(cub, alpha)
  weight = (2.0/(N-1))*ones(T, (numedge + 1 + centroid))
  SymCubatures.setweights!(cub, weight)

  if q <= 1
    cub.weights = (2.0/N)*ones(T, (numedge + 1 + centroid))
  elseif q <= 2
    cub.weights[1] = (2.0/(N-1))*1.0/2.0
  elseif q <= 3
    cub.weights[1:2] = (2.0/(N-1)).* T[5/12, 13/12]
  elseif q <= 4
    cub.weights[1:3] = (2.0/(N-1)).* T[3/8, 7/6, 23/24]
  elseif q <= 5
    cub.weights[1:4] = (2.0/(N-1)).* T[251/720, 299/240, 211/240, 739/720]
  elseif q <= 6
    cub.weights[1:5] = (2.0/(N-1)).* T[95/288, 317/240, 23/30, 793/720, 157/160]
  elseif q <= 7
    cub.weights[1:6] = (2.0/(N-1)).* T[19087/60480, 84199/60480, 18869/30240, 37621/30240, 55031/60480, 61343/60480]
  elseif q <= 8
    cub.weights[1:7] = (2.0/(N-1)).* T[5257/17280, 22081/15120, 54851/120960, 103/70, 89437/120960, 16367/15120, 23917/24192]
  elseif q <= 9
    cub.weights[1:8] = (2.0/(N-1)).* T[1070017/3628800, 5537111/3628800, 103613/403200, 261115/145152, 298951/725760, 515677/403200, 3349879/3628800, 3662753/3628800]
  end
  vtx = reshape(T[-1; 1], (2,1))
  return cub, vtx
end

"""
### Cubature.getTriCubatureGregory{T}

Returns a cubature rule and vertices for the SBP DiagE operators on triangles;
these are cubatures that have nodes on the boundary that correspond with 
Gregory quadrature rules, which then leads to diagonal E.

**Inputs**

* `q`: maximum degree of polynomial for which the cubature is exact
* `T`: the data type used to represent the cubature
* `vertices`: if true then vertices are included

**Outputs**

* `cub`: an asymmetric cubature for the right triangle
* `vtx`: vertices for the right triangle

"""
function getTriCubatureGregory(q::Int, N::Int; T=Float64)
  quad, vtx_1d = Cubature.getLineCubatureGregory(q, N)
  nodes1d = SymCubatures.calcnodes(quad,vtx_1d)
  perm = sortperm(vec(nodes1d))
  nodes1d = nodes1d[perm]
  weights1d = (SymCubatures.calcweights(quad))[perm]
  weights2d = zeros(T,convert(Int,(N+1)*N/2))
  nbnd = q-1

  cub = AsymCubatures.TriAsymCub{T}(numedgenodes=N) 
  AsymCubatures.setparams!(cub, nodes1d)

  node = 1
  for j=1:N
    for i=1:(N-j+1)
      dist = N-j+2-i
      # node = convert(Int64,(N+1)*(j-1) - j*(j-1)/2 + i)
      if ((i <= nbnd) && (j <= nbnd)) 
        weights2d[node] = weights1d[i]*weights1d[j]
      elseif ((j <= nbnd) && (dist <= nbnd))
        weights2d[node] = weights1d[j]*weights1d[dist]
      elseif ((i <= nbnd) && (dist <= nbnd))
        weights2d[node] = weights1d[i]*weights1d[dist]
      elseif ((i > nbnd) && (dist <= nbnd))            
        weights2d[node] = weights1d[maximum([i,j])]*weights1d[dist]
      else
        weights2d[node] = weights1d[i]*weights1d[j]
      end
      node += 1
    end
  end

  AsymCubatures.setweights!(cub, weights2d)
  vtx = T[-1 -1; 1 -1; -1 1]

  return cub, vtx
end

"""
### Cubature.getTetCubatureGregory{T}

Returns a cubature rule and vertices for the SBP DiagE operators on tetrahedra;
these are cubatures that have nodes on the boundary that correspond with 
Gregory quadrature rules, which then leads to diagonal E.

**Inputs**

* `q`: maximum degree of polynomial for which the cubature is exact
* `T`: the data type used to represent the cubature
* `vertices`: if true then vertices are included

**Outputs**

* `cub`: an asymmetric cubature for the right triangle
* `vtx`: vertices for the right triangle

"""
function getTetCubatureGregory(q::Int, N::Int; T=Float64)
  quad, vtx_1d = Cubature.getLineCubatureGregory(q, N)
  nodes1d = SymCubatures.calcnodes(quad,vtx_1d)
  perm = sortperm(vec(nodes1d))
  nodes1d = nodes1d[perm]
  weights1d = (SymCubatures.calcweights(quad))[perm]
  weights3d = zeros(T,convert(Int,(N+1)*(N+2)*N/6))
  nbnd = q-1

  cub = AsymCubatures.TetAsymCub{T}(numedgenodes=N) 
  AsymCubatures.setparams!(cub, nodes1d)

  node = 1
  for k=1:N
    for j=1:N-(k-1)
      for i=1:N-(j-1)-(k-1)
        dist = N-i-j-k+3
        if ((i <= nbnd) && (j <= nbnd) && (k <= nbnd)) 
          weights3d[node] = weights1d[i]*weights1d[j]*weights1d[k]
        elseif ((i <= nbnd) && (j <= nbnd) && (dist <= nbnd))
          weights3d[node] = weights1d[i]*weights1d[j]*weights1d[dist]
        elseif ((j <= nbnd) && (k <= nbnd) && (dist <= nbnd))
          weights3d[node] = weights1d[dist]*weights1d[j]*weights1d[k]
        elseif ((i <= nbnd) && (k <= nbnd) && (dist <= nbnd))
          weights3d[node] = weights1d[i]*weights1d[dist]*weights1d[k]
        elseif ((k <= nbnd) && (i >nbnd) && (dist <= nbnd))
          weights3d[node] = weights1d[maximum([i,j])]*weights1d[dist]*weights1d[k]
        elseif ((i <= nbnd) && (j >nbnd) && (dist <= nbnd))
          weights3d[node] = weights1d[i]*weights1d[maximum([j,k])]*weights1d[dist]
        elseif ((j <= nbnd) && (k >nbnd) && (dist <= nbnd))
          weights3d[node] = weights1d[maximum([i,k])]*weights1d[j]*weights1d[dist]
        elseif ((i > nbnd) && (dist <= nbnd))            
          weights3d[node] = weights1d[maximum([i,j])]*weights1d[maximum([i,k])]*weights1d[dist]
        else
          weights3d[node] = weights1d[i]*weights1d[j]*weights1d[k]
        end
        node += 1
      end
    end
  end
  AsymCubatures.setweights!(cub, weights3d)
  vtx = T[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]

  return cub, vtx
end

"""
### Cubature.tricubature{T}

Deprecated; this function will be removed in the future

**Inputs**

* `q`: maximum degree of polynomial for which the cubature is exact
* `T`: the data type used to represent the cubature

**Outputs**

* `cub`: a symmetric cubature for the right triangle
* `vtx`: vertices for the right triangle

"""
function tricubature(q::Int, T=Float64)
  return getTriCubatureGamma(q, T)
end

"""
### Cubature.getTriCubatureGamma{T}

Returns a cubature rule and vertices for the SBP Gamma operators on triangles;
these are operators with p+1 nodes on each face, where, typically, p =
(`q`+1)/2.

**Inputs**

* `q`: maximum degree of polynomial for which the cubature is exact
* `T`: the data type used to represent the cubature
* `tol`: tolerance with which to solve the cubature

**Outputs**

* `cub`: a symmetric cubature for the right triangle
* `vtx`: vertices for the right triangle

"""
function getTriCubatureGamma(q::Int, T=Float64;
                             tol=10*eps(typeof(real(one(T)))))
  @assert( q >= 1 && q <= 9 && mod(q,2) == 1 )
  cub_degree = q
  mask = zeros(Int64, (0))
  if q <= 1
    # P1 (vertices only); 2nd order cubature
    cub = SymCubatures.TriSymCub{T}(vertices=true) 
    SymCubatures.setweights!(cub, T[2/3])
    cub_degree = 1
  elseif q <= 3
    # P2 + 1 bubble node; 4th order cubature
    cub = SymCubatures.TriSymCub{T}(midedges=true, centroid=true)
    SymCubatures.setweights!(cub, T[9/10, 1/10, 4/15])
    cub_degree = 3
  elseif q <= 5
    # P3 + 3 bubble nodes; 6th order cubature
    cub = SymCubatures.TriSymCub{T}(numedge=1, numS21=1)
    SymCubatures.setweights!(cub, T[0.02974582604964118,0.4415541156808217,
                                    0.09768336246810204])
    SymCubatures.setparams!(cub, T[0.41469035132718185,0.29346955590904017])
    cub_degree = 5
  elseif q <= 7
    # P4 + 6 bubble nodes; 8th order cubature
    cub = SymCubatures.TriSymCub{T}(midedges=true, numedge=1, numS21=2)
    SymCubatures.setweights!(cub, T[0.012698412698412695,0.05079365079365077,
                                    0.2023354595827503,0.3151248578775673,
                                    0.04285714285714284])
    SymCubatures.setparams!(cub, T[0.2615831876594899,0.8495279234516212,
                                   0.2113248654051872])
    cub_degree = 7
  elseif q <= 9
    # P5 + 10 bubble nodes; 10th order cubature
    cub = SymCubatures.TriSymCub{T}(numedge=2, centroid=false, numS21=2,
                                    numS111=1)
    SymCubatures.setweights!(cub, T[0.005060870857201095,0.1656605522739661,
                                    0.11124762548151651,0.02601835402818015,
                                    0.022007571591738946,0.1443228834070724])
    SymCubatures.setparams!(cub, T[0.5264875797340474,0.2020312767621901,
                                   0.3647863788168577,0.12582399442561498,
                                   0.17313630713608186,0.6472196801547492])
    cub_degree = 9
  elseif q <= 13
    # P7; 14th order cubature
    cub = SymCubatures.TriSymCub{T}(numedge=3, centroid=false, numS21=3,
                                    numS111=3)
    SymCubatures.setweights!(cub, T[0.0013434826332230758,0.07754170837489897,
                                    0.0392937103109862,0.08132263825474927,
                                    0.009447719133224907,0.011139320053630379,
                                    0.007018823640551441,0.055270427044833675,
                                    0.06131670240670696,0.08938957126745721])
    SymCubatures.setparams!(cub, T[0.5749380303918797,0.11974220085793316,
                                   0.33343014155598055,0.20437477140696825,
                                   0.39432452613154606,0.06536359154212519,
                                   0.38875869465672286,0.10195028749023816,
                                   1.1491810826793598,0.09609353164480232,
                                   0.6657786329998556,1.0308822535578346])
    cub_degree = 13
  else
    error("polynomial degree must be <= 9 (presently)\n")
  end
  mask = 1:(cub.numparams+cub.numweights)
  vtx = T[-1 -1; 1 -1; -1 1]
  Cubature.solvecubature!(cub, cub_degree, mask, tol=tol)
  return cub, vtx
end

"""
### Cubature.getTriCubatureOmega{T}

Returns a cubature rule and vertices for the SBP Omega operators on triangles;
these are cubatures that are analogous to Gauss-Legendre in 1D, and they are
strictly internal to the triangle. 

**Inputs**

* `q`: maximum degree of polynomial for which the cubature is exact
* `T`: the data type used to represent the cubature
* `tol`: tolerance with which to solve the cubature

**Outputs**

* `cub`: a symmetric cubature for the right triangle
* `vtx`: vertices for the right triangle

"""
function getTriCubatureOmega(q::Int, T=Float64;
                             tol=10*eps(typeof(real(one(T)))))
  cub_degree = q
  mask = zeros(Int64, (0))
  if q <= 1
    # P1; 3rd order cubature
    cub = SymCubatures.TriSymCub{T}(vertices=false, numS21=1)      
    SymCubatures.setweights!(cub, T[2/3])
    SymCubatures.setparams!(cub, T[4/5])
    cub_degree = 1
  elseif q <= 2
    # P1; 3rd order cubature
    cub = SymCubatures.TriSymCub{T}(vertices=false, numS21=1)      
    SymCubatures.setweights!(cub, T[2/3])
    SymCubatures.setparams!(cub, T[1/3])
    cub_degree = 2
  elseif q <= 3
    # P2; 4th order cubature
    cub = SymCubatures.TriSymCub{T}(vertices=false, numS21=2, numS111=0)
    SymCubatures.setparams!(cub, T[0.32232557541453105, 0.9443569106417902])
    SymCubatures.setweights!(cub, T[0.5467137619817344, 0.11995290468493017])
    cub_degree = 3
  elseif q <= 4
    # P2; 5th order cubature
    cub = SymCubatures.TriSymCub{T}(vertices=false, numS21=2)
    SymCubatures.setweights!(cub, T[0.44676317935602283;
                                    0.2199034873106437])
    SymCubatures.setparams!(cub, T[0.8918969818319298;
                                   0.18315242701954149])
    cub_degree = 4
  elseif q <= 5
    # P3; 6th order cubature
    cub = SymCubatures.TriSymCub{T}(vertices=false, centroid=true,
                                    numS21=1, numS111=1)
    SymCubatures.setweights!(cub, T[0.11550472674301035;
                                    0.20924480696331949;
                                    0.39801697799105223])
    SymCubatures.setparams!(cub, T[0.13862330627662678;
                                   0.14215944055500324;
                                   0.6226442585632832])
    cub_degree = 5
  elseif q <= 6
    # P3; 7th order cubature
    cub = SymCubatures.TriSymCub{T}(vertices=false, centroid=false,
                                    numS21=2, numS111=1)
    SymCubatures.setweights!(cub, T[0.1016898127404136;
                                    0.23357255145275863;
                                    0.16570215123674722])
    SymCubatures.setparams!(cub, T[0.1261780289830045;
                                   0.49857349034182097;
                                   0.10629009968963399;
                                   0.6207049020675687])
    cub_degree = 6
  elseif q <= 7
    # P4; 8th order cubature
    cub = SymCubatures.TriSymCub{T}(vertices=false, numS21=3, numS111=1)
    SymCubatures.setweights!(cub, T[0.045386157905236965;
                                    0.1458284149509071;
                                    0.2543369199180239;
                                    0.11055758694624956])
    SymCubatures.setparams!(cub, T[0.08433122881886425;
                                   0.9485893782350207;
                                   0.4841719475189572;
                                   0.4231241172761849;
                                   0.0959626827429292])
    cub_degree = 7
    # JEH The commented version below has much worse conditioning/CFL limit;
    # leaving it here for reference, and to avoid
    # SymCubatures.setweights!(cub, T[0.10482661091570668;
    #                                 0.2253930198733382;
    #                                 0.057547518977195254;
    #                                 0.13944975845021326])
    # SymCubatures.setparams!(cub, T[0.1290634461434249;
    #                                0.4731163893279408;
    #                                0.8413468069012109;
    #                                0.08815074437486997;
    #                                0.624003943088726])
  elseif q <= 8
    # P4; 9th order cubature
    # cub = SymCubatures.TriSymCub{T}(vertices=false, centroid=true,
    #                                 numS21=3, numS111=1)
    # SymCubatures.setweights!(cub, T[0.20643474106943666;
    #                                 0.19018326853457126;
    #                                 0.06491699524639659;
    #                                 0.05446062834886685;
    #                                 0.2886312153555818])
    # SymCubatures.setparams!(cub, T[0.3411386155035168;
    #                                0.9185851765854467;
    #                                0.10109445663405793;
    #                                0.5262256592692812;
    #                                0.016789554819910586])      
    cub = SymCubatures.TriSymCub{T}(vertices=false, numS21=4, numS111=1)
    SymCubatures.setweights!(cub, T[0.05415768359243055;
                                    0.11606844251884181;
                                    0.16305903205856434;
                                    0.18042409969012593;
                                    0.07647870440335205])
    SymCubatures.setparams!(cub, T[0.09199464041197784;
                                   0.9547402313677645;
                                   0.353676230440559;
                                   0.8037465193903021;
                                   0.4612850345245523;
                                   0.06105167151116454])
    cub_degree = 8
    tol = 6e-14
  elseif q <= 9 
    # 21 nodes
    cub = SymCubatures.TriSymCub{T}(vertices=false, numS21=3, numS111=2, centroid=false)
    SymCubatures.setparams!(cub, T[0.09037801956875349, 0.9630396695666223, 0.8072079596358801, 0.43658014194276157, 0.27398240252980893, 0.06084872345763977, 0.4441263310746356])
    SymCubatures.setweights!(cub, T[0.05198714206463918, 0.10323440513804151, 0.18816014691671232, 0.09093907609523856, 0.07070341017839829])

    #19 nodes (not enough nodes to construct degree 4 operator)    
    # cub = SymCubatures.TriSymCub{T}(vertices=false, numS21=4, numS111=1, centroid=true)
    # SymCubatures.setparams!(cub, T[0.3764070712380649, 0.9793650383974742, 0.08945902678890531, 0.8741791829858713, 0.07367682410947267, 0.4439259783215317])
    # SymCubatures.setweights!(cub, T[0.1592954778544202, 0.06266940045428016, 0.05115535131739605, 0.15565508200954808, 0.08656707875457871, 0.19427159256559468])
    
    cub_degree = 9
    tol = 5e-15
  elseif q <= 10 # 25 nodes     
    cub = SymCubatures.TriSymCub{T}(vertices=false, numS21=2, numS111=3, centroid=true)
    SymCubatures.setparams!(cub, T[0.06411074643388708, 0.28432220211312814, 0.296265771567641, 1.0601082378546884, 0.32740346747436544, 0.05673533067987688, 0.05923977897745965, 0.7382935636556222])
    SymCubatures.setweights!(cub, T[0.02670593762629918, 0.09191592720948939, 0.12780981279284825, 0.05059551541457666, 0.0683692963259188, 0.16348665829257195])
    cub_degree = 10
    tol = 6e-14
  elseif q <= 11 # 28 nodes     
    # cub = SymCubatures.TriSymCub{T}(vertices=false, numS21=5, numS111=2, centroid=true)
    # SymCubatures.setparams!(cub, T[0.05233095719752731, 0.9894031463490898, 0.18960563931430668, 0.8777841974682057, 0.41501935044082044, 0.2832849334674954, 0.0010190676116574875, 0.5558104967282995, 1.3543678890763444])
    # SymCubatures.setweights!(cub, T[0.017711439632666408, 0.03746416080447641, 0.07620451986110559, 0.13844887139337872, 0.14413384886902336, 0.015060319238336077, 0.08200532556201282, 0.17571760951595808])
    #34 nodes
    # cub = SymCubatures.TriSymCub{T}(vertices=false, numS21=1, numS111=5, centroid=true)
    # SymCubatures.setparams!(cub, T[0.06517612308759263, 0.04519866864934196, 0.31848143875936, 0.0370316853857241, 0.7045798860690692, 0.09129861405067456, 0.7916177719149112, 0.34312389714411, 0.19429850362582882, 0.30641879844403425, 1.0637524989723475])
    # SymCubatures.setweights!(cub, T[0.027082523005573095, 0.03745121993811072, 0.0344496896152223, 0.03779638743243422, 0.05794564876237383, 0.12592411810881812, 0.15735004784152562])
    # 30 nodes
    cub = SymCubatures.TriSymCub{T}(vertices=false, numS21=4, numS111=3, centroid=false)
    SymCubatures.setparams!(cub, T[0.0653339850268845, 0.7848925385554962, 0.4576439501834159, 0.2748895760200856, 0.016564338513558364, 0.7635155543630708, 0.6922762730907215, 0.1566169791427997, 0.33206144090886697, 0.0532930008473203])
    SymCubatures.setweights!(cub, T[0.027582256878190624, 0.11392363359624305, 0.09159074449758878, 0.08739254511710136, 0.026587894037230206, 0.09818904821579671, 0.04831180103574457])

    cub_degree = 11
    tol = 6e-14
  elseif q <= 12 # 33 nodes     
    cub = SymCubatures.TriSymCub{T}(vertices=false, numS21=5, numS111=3, centroid=false)
    SymCubatures.setparams!(cub, T[0.5424207700242318, 0.04263470090642063, 0.25515229108317217, 0.8794487845889205, 0.9764347795476097, 0.051468101096660494, 0.2325038318151941, 0.04567666444451398, 0.5626511619798789, 0.5514265393710289, 0.23068698906939594])
    SymCubatures.setweights!(cub, T[0.12571644843577023, 0.012332522103117988, 0.0695922258614181, 0.08738508907607664, 0.051462132880910685, 0.03463246221731784, 0.04471354640460689, 0.08074311553276181])
    cub_degree = 12
    tol = 6e-14
  elseif q <= 13 # 37 nodes     
    cub = SymCubatures.TriSymCub{T}(vertices=false, numS21=5, numS111=3, centroid=false)
    SymCubatures.setparams!(cub, T[0.04301936221768631, 0.44274457258366556, 0.9781538929050788, 0.8538828285196006, 0.3271948021357017, 0.17579096606439407, 0.22184408560692612, 0.04874037380218756, 0.13602448710841317, 0.6168835217842368, 0.010252778204764737, 0.5450316355468587])
    SymCubatures.setweights!(cub, T[0.012104674207078276, 0.1165569702384, 0.04798880385778925, 0.11120393506090662, 0.04835807962318797, 0.02993080221033121, 0.06928255228169665, 0.019181362007086523, 0.1359200731736634])
    cub_degree = 13
    tol = 5e-15
  else
    error("polynomial degree must be <= 20 (presently)\n")
  end
  mask = 1:(cub.numparams+cub.numweights)
  vtx = T[-1 -1; 1 -1; -1 1]
  Cubature.solvecubature!(cub, cub_degree, mask, tol=tol)
  return cub, vtx
end

"""
### Cubature.getTriCubatureForTetFaceDiagE{T}

Returns a cubature rule and vertices for facets of the SBP DiagE operators 
on tetrahedra; these should not be used for 2D problems as they do not 
satisfy the accruacy requirements along their edge.

**Inputs**

* `q`: maximum degree of polynomial for which the cubature is exact
* `T`: the data type used to represent the cubature
* `faceopertype`: the operator type on the facets of the tetrahedron
* `tol`: tolerance with which to solve the cubature

**Outputs**

* `cub`: a symmetric cubature for the right triangle
* `vtx`: vertices for the right triangle

"""
function getTriCubatureForTetFaceDiagE(q::Int, T=Float64; faceopertype::Symbol=:DiagE,
  tol=10*eps(typeof(real(one(T)))))
  cub_degree = q
  mask = zeros(Int64, (0))

  if faceopertype == :Omega
    cub, vtx = getTriCubatureOmega(q)
    tol=1e-14
  else
    if q<=2 #3 nodes
      cub = SymCubatures.TriSymCub{T}(vertices=false, midedges=true, centroid=false)
      SymCubatures.setweights!(cub, T[0.6666666666666666])
      cub_degree = 2
      tol = 5e-15

      #4 nodes
      # cub = SymCubatures.TriSymCub{T}(vertices=true, midedges=false, centroid=true)
      # SymCubatures.setweights!(cub, T[0.16666666666666624, 1.4999999999999973])
      # cub_degree = 2
      # tol = 5e-15
    elseif q<=4 #9 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = true,
                                      centroid = false,
                                      numedge = 0,
                                      numS21 = 1,
                                      numS111 = 0)
      SymCubatures.setparams!(cub, T[0.3771609693928902])
      SymCubatures.setweights!(cub, T[0.0410802706918665, 0.123950997736534, 0.5016353982382646])
      cub_degree = 4
      tol = 5e-15

      #7 nodes, but this leads to tet element with 26 nodes while the above gives 23 nodes
      # cub = SymCubatures.TriSymCub{T}(vertices = false,
      #                                 midedges = true,
      #                                 centroid = true,
      #                                 numedge = 0,
      #                                 numS21 = 1,
      #                                 numS111 = 0)
      # SymCubatures.setparams!(cub, T[0.22222222222222215])
      # SymCubatures.setweights!(cub, T[0.1523809523809521, 0.28928571428571387, 0.6749999999999992])
      # cub_degree = 4
      # tol = 5e-15
    elseif q<=6 #15 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = false,
                                      centroid = false,
                                      numedge = 1,
                                      numS21 = 2,
                                      numS111 = 0)
      SymCubatures.setparams!(cub, T[0.8506802519794945, 0.23722737279318576, 0.3077459416259917])
      SymCubatures.setweights!(cub, T[0.01426071861440897, 0.3303589772911334, 0.20376930605390392, 0.059138832353610636])
      cub_degree = 6
      tol = 5e-15
      
      # omega type facet nodes
      # cub = SymCubatures.TriSymCub{T}(vertices=false, centroid=false,
      #                                 numS21=2, numS111=1)
      # SymCubatures.setweights!(cub, T[0.1016898127404136;
      #                                 0.23357255145275863;
      #                                 0.16570215123674722])
      # SymCubatures.setparams!(cub, T[0.1261780289830045;
      #                                0.49857349034182097;
      #                                0.10629009968963399;
      #                                0.6207049020675687])
      # cub_degree = 6
    elseif q<=8 #22 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = true,
                                      centroid = true,
                                      numedge = 1,
                                      numS21 = 1,
                                      numS111 = 1)
      SymCubatures.setparams!(cub, T[0.16099183834007516, 0.800367892880542, 0.6058255660767269, 0.21518364356973504])
      SymCubatures.setweights!(cub, T[0.006036623735435305, 0.040748592504715506, 0.09377442432581949, 0.02798160585550063, 0.19106553442197458, 0.2640382366372378])
      cub_degree = 8
      tol = 7e-15

      # omega type facet nodes
      # cub = SymCubatures.TriSymCub{T}(vertices=false, centroid=true,
      # numS21=3, numS111=1)
      # SymCubatures.setweights!(cub, T[0.20643474106943666;
      #   0.19018326853457126;
      #   0.06491699524639659;
      #   0.05446062834886685;
      #   0.2886312153555818])
      # SymCubatures.setparams!(cub, T[0.3411386155035168;
      #   0.9185851765854467;
      #   0.10109445663405793;
      #   0.5262256592692812;
      #   0.016789554819910586])      
      # cub_degree = 8
      # tol = 6e-14
    elseif q<=10 #28 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = true,
                                      centroid = true,
                                      numedge = 1,
                                      numS21 = 3,
                                      numS111 = 1)
      SymCubatures.setparams!(cub, T[0.40849548761177584, 0.1693510205976773, 0.895541492424858, 0.9192858150590572, 0.5199457984883094, 0.07662224108118755])
      SymCubatures.setweights!(cub, T[0.0017341485435839296, 0.02676086117258247, 0.1584314816374832, 0.07114324960955178, 0.15198175500182232, 0.013615409089319942, 0.08154977731366873, 0.1988543936869974])
      cub_degree = 10
      tol = 5e-15
    else
      error("polynomial degree must be <= 10 (presently)\n")
    end
  end
  mask = SymCubatures.getInternalParamMask(cub)
  append!(mask, (cub.numparams+1):(cub.numparams+cub.numweights))
  vtx = T[-1 -1; 1 -1; -1 1]
  Cubature.solvecubature!(cub, cub_degree, mask, tol=tol)
  return cub, vtx
end

"""
### Cubature.getTriCubatureDiagE{T}

Returns a cubature rule and vertices for the SBP DiagE operators on triangles;
these are cubatures that have nodes on the boundary that correspond with LG or
LGL quadrature rules, which then leads to diagonal E.

**Inputs**

* `q`: maximum degree of polynomial for which the cubature is exact
* `T`: the data type used to represent the cubature
* `vertices`: if true then vertices are included
* `tol`: tolerance with which to solve the cubature

**Outputs**

* `cub`: a symmetric cubature for the right triangle
* `vtx`: vertices for the right triangle

"""
function getTriCubatureDiagE(q::Int, T=Float64; vertices::Bool=true,
                             tol=10*eps(typeof(real(one(T)))))
  cub_degree = q
  mask = zeros(Int64, (0))
  if vertices
    # include the vertices in the cubature
    if q<=1 #6 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = true,
                                      centroid = false,
                                      numedge = 0,
                                      numS21 = 0,
                                      numS111 = 0)                
      SymCubatures.setweights!(cub, T[0.5525545450892143; 0.11411212157745271])
      # SymCubatures.setweights!(cub, T[1/6; 1/2])
      # SymCubatures.setweights!(cub, T[1/8; 13/24])
      cub_degree = 1
      tol = tol

    elseif q<=2 #7 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = true,
                                      centroid = true,
                                      numedge = 0,
                                      numS21 = 0,
                                      numS111 = 0)

      SymCubatures.setweights!(cub, T[0.04761904761904726, 0.47619047619047483, 0.4285714285714286])
      cub_degree = 2
      tol = 5e-15

      # SymCubatures.setweights!(cub, T[9/10, 1/10, 4/15])
      # cub_degree = 3
      # tol = 1e-14
    elseif q <=3 #10 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = false,
                                      centroid = true,
                                      numedge = 1,
                                      numS21 = 0,
                                      numS111 = 0)

      SymCubatures.setparams!(cub, T[0.7236067977499789])
      SymCubatures.setweights!(cub, T[0.03333333333333329, 0.16666666666666605, 0.8999999999999969])
      cub_degree = 3
      tol=1e-14
    elseif q <= 4
      #12 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = false,
                                      numS21 = 1,
                                      numedge = 1,
                                      numS111 = 0,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.4257087142236166, 0.7236067977499789])
      SymCubatures.setweights!(cub, T[0.025044506019598876, 0.42703175864395354, 0.10729520100155711])
      # cub = SymCubatures.TriSymCub{T}(vertices=true, numedge=1,
      #                                 midedges=false,
      #                                 numS21=1, numS111=0,
      #                                 centroid=false)
      # SymCubatures.setparams!(cub, T[0.4257087142201423;
      #                                0.5*(1 + sqrt(1/5))])
      # SymCubatures.setweights!(cub, T[0.02504450602156441;
      #                                 0.4270317586471588;
      #                                 0.10729520099967835])
      cub_degree = 4
      tol = 1e-14
    elseif q <=5 #15 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = true,
                                      centroid = false,
                                      numedge = 1,
                                      numS21 = 1,
                                      numS111 = 0)
      SymCubatures.setparams!(cub, T[0.41469035132718185; 0.8273268353539885])
      SymCubatures.setweights!(cub, T[0.014698618394803228, 0.09752600361864236, 0.44155411568082115, 0.056443964486199594])
      # SymCubatures.setweights!(cub, T[0.02323751046092639; 0.07171714179335946; 0.4415541156808216; 0.06507894936577964])
      cub_degree = 5
      tol = 1e-14
    elseif q <= 6 #18 nodes
      # cub = SymCubatures.TriSymCub{T}(vertices = true,
      #                                 midedges = true,
      #                                 numS21 = 2,
      #                                 numedge = 1,
      #                                 numS111 = 0,
      #                                 centroid = false)
      # SymCubatures.setparams!(cub, T[0.8487720503437628, 0.28401016819355557, 0.8273268353539885])
      # SymCubatures.setweights!(cub, T[0.009130264572198617, 0.06201621057208993, 0.29870677879935964, 0.20607267198227744, 0.045370370370370325])
      # cub_degree = 6
      # tol = 1e-14
      #18 nodes
      # cub = SymCubatures.TriSymCub{T}(vertices=true, numedge=1,
      #                                 midedges=true,
      #                                 numS21=2, numS111=0, centroid=false)
      # SymCubatures.setparams!(cub, T[0.8487720503426771;
      #                                0.28401016818370567;
      #                                0.5*(1 + sqrt(3/7))])
      # SymCubatures.setweights!(cub, T[0.00913026457472031;
      #                                 0.06201621056804736;
      #                                 0.2987067788024998;
      #                                 0.20607267197855683;
      #                                 0.04537037036896131])
      # tol = 1e-14
      #19 nodes
      # cub = SymCubatures.TriSymCub{T}(vertices=true, 
      #                                 midedges=true, 
      #                                 numS21=2, 
      #                                 numedge=1, 
      #                                 numS111=0, 
      #                                 centroid=true)
      # SymCubatures.setparams!(cub, T[0.2841256618085499, 0.849330059385702, 0.8273268353539885])
      # SymCubatures.setweights!(cub, T[0.009139015771102226, 0.06184185447856355, 0.20623085642293718, 0.2972006564063465, 0.04537037037037026, 0.0045406285409288994])
      # cub_degree = 6
      # tol = 1e-14
      cub = SymCubatures.TriSymCub{T}(vertices=true, 
                                      midedges=true, 
                                      numS21=2, 
                                      numedge=1, 
                                      numS111=0, 
                                      centroid=true)
      SymCubatures.setparams!(cub, T[0.2891640824232583, 0.8750401312391274, 0.8273268353539885])
      SymCubatures.setweights!(cub, T[0.009519647256771995, 0.05213852810721973, 0.21326078563097062, 0.2444355957641196, 0.04537037037037035, 0.16971410750053226])
      cub_degree = 6
      tol = 1e-14
    elseif q <=7 #24 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = false,
                                      centroid = false,
                                      numedge = 2,
                                      numS21 = 1,
                                      numS111 = 1)
      SymCubatures.setparams!(cub, T[0.8461370386526059, 0.8825276619647324, 0.642615758240322, 0.33879488493771764, 0.19830545894321575])
      SymCubatures.setweights!(cub, T[0.00753161345765886, 0.30480530834098646, 0.02108627120396598, 0.045504525339300515, 0.1105740758907442])
      cub_degree = 7
      tol=5e-15
    elseif q <= 8 # 27 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = false,
                                      numS21 = 2,
                                      numedge = 2,
                                      numS111 = 1,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.5306627609684195, 0.20735501628561037, 0.8825276619647324, 0.6426157582403226, 0.17654792120316234, 0.6492809445301031])
      SymCubatures.setweights!(cub, T[0.004361575619973876, 0.15984019211140704, 0.1135887626556459, 0.022078990549428416, 0.02786777148585481, 0.14449130610453667])
      cub_degree = 8
      tol = 1e-14
      # cub = SymCubatures.TriSymCub{T}(vertices=true, numedge=2,
      #                                 numS21=2, numS111=1)
      # SymCubatures.setparams!(cub, T[0.20735501628574252;
      #                                0.5306627609507977;
      #                                0.5*(1 + sqrt(1/3 - 2*sqrt(7)/21));
      #                                0.5*(1 + sqrt(1/3 + 2*sqrt(7)/21));
      #                                0.6492809445444747;
      #                                1.1741711342683223])
      # SymCubatures.setweights!(cub, T[0.004361575620524937;
      #                                 0.11358876265867929;
      #                                 0.15984019213585915;
      #                                 0.027867771483031517;
      #                                 0.02207899054885172;
      #                                 0.14449130609379232])  
      # cub_degree = 8

      # tol = 1e-14
    elseif q <=9 #34 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = true,
                                      centroid = false,
                                      numedge = 2,
                                      numS21 = 3,
                                      numS111 = 1)
      SymCubatures.setparams!(cub, T[0.5110998040469292, 0.9089231881407306, 0.19994428870176986, 0.9151119481392835, 0.7344243967353571, 0.1465093205286518, 0.5715830463482526])
      SymCubatures.setweights!(cub, T[0.0019241266200829663, 0.022987519027599365, 0.1866822178508963, 0.09270056884965712, 0.0985056912129899, 0.016540728394130698, 0.01779352639593344, 0.09759901676265635])
      cub_degree = 9
      tol=5e-15
    elseif q <= 10 #36 nodes
      # cub = SymCubatures.TriSymCub{T}(vertices = true,
      #                                 midedges = true,
      #                                 numS21 = 2,
      #                                 numedge = 2,
      #                                 numS111 = 2,
      #                                 centroid = false)
      # SymCubatures.setparams!(cub, T[0.1025947577289524, 0.514543490059506, 0.9151119481392835, 0.7344243967353571, 0.16739689557448953, 0.736084860943316, 0.37663587800802634, 0.15391188146518084])
      # SymCubatures.setweights!(cub, T[0.001135617009180763, 0.022982971732153294, 0.03238151607768464, 0.1823790802107045, 0.009620012578487587, 0.019610669249206004, 0.09949978404801857, 0.08516327494275952])
      # cub_degree = 10
      # tol = 1e-14
      # # 34 nodes (this is better than what we have in the paper)
      # cub = SymCubatures.TriSymCub{T}(vertices = true,
      #                                 midedges = true,
      #                                 numS21 = 3,
      #                                 numedge = 2,
      #                                 numS111 = 1,
      #                                 centroid = true)
      # SymCubatures.setparams!(cub, T[0.8963473087865563, 0.4216026177435242, 0.1768201173272758, 0.9358700742548033, 0.7958500907165711, 0.09048622650978146, 1.38127245795577])
      # SymCubatures.setweights!(cub, T[0.0001400555146967039, 0.027547415560250886, 0.1445482011568161, 0.15428412182476955, 0.07658455873454234, 0.012503331996950234, 0.00770237566957249, 0.0802768737939981, 0.18779145286364823])
      # cub_degree = 10
      # tol = 1e-14

      # 37 nodes for sparse
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = true,
                                      numS21 = 4,
                                      numedge = 2,
                                      numS111 = 1,
                                      centroid = true)
      SymCubatures.setparams!(cub, T[0.4127652030196255, 0.1577404405372036, 0.9383113442093055, 0.828770971730687, 0.9151119481392835, 0.7344243967353571, 0.5020748121582398, 1.3787314834940385])
      SymCubatures.setweights!(cub, T[0.002253032562263613, 0.01822260666165858, 0.12932769175775988, 0.06624360894290958, 0.08375983980529911, 0.10840144360935694, 0.012240466552266514, 0.013479866725941795, 0.08349773113726795, 0.12006694348939934])
      cub_degree = 10
      tol = 1e-14
    elseif q <=11 #40 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = false,
                                      centroid = true,
                                      numedge = 3,
                                      numS21 = 4,
                                      numS111 = 1)
      SymCubatures.setparams!(cub, T[0.40821896354263243, 0.981013735325053, 0.16748106658466438, 0.8690323743735683, 0.9358700742548033, 0.7958500907165711, 0.6046496089512394, 0.10247525036525074, 0.5346096921202054])
      SymCubatures.setweights!(cub, T[0.0008860287032078429, 0.1399375485517579, 0.04956216973314194, 0.07404852901772899, 0.13638424140745345, 0.010480045148169446, 0.011623548319685177, 0.0012118355050860282, 0.08102106300279119, 0.1715254959057385])
      cub_degree = 11
      tol=5e-15
      # 43 nodes
      # cub = SymCubatures.TriSymCub{T}(vertices = true,
      #                                 midedges = false,
      #                                 centroid = true,
      #                                 numedge = 3,
      #                                 numS21 = 3,
      #                                 numS111 = 2)
      # SymCubatures.setparams!(cub, T[0.8571396116924792, 0.42137206563390533, 0.15416411447545023, 0.9358700742548033, 0.7958500907165711, 0.6046496089512394, 0.8193910152685325, 0.07542901871217762, 0.1301469442971666, 0.4770052452777402])
      # SymCubatures.setweights!(cub, T[0.0011717118884880214, 0.12473132707644294, 0.12730130899704278, 0.06274770142440542, 0.009242885505864128, 0.013103502174005579, 0.008618295447807097, 0.04030615822887277, 0.07902105209405838, 0.15039249113721462])
      # cub_degree = 11
      # tol=5e-15
      # 43 nodes
      # cub = SymCubatures.TriSymCub{T}(vertices = true,
      #                                 midedges = false,
      #                                 centroid = true,
      #                                 numedge = 3,
      #                                 numS21 = 5,
      #                                 numS111 = 1)
      # SymCubatures.setparams!(cub, T[0.8690570506162132, 0.4099529922719056, 0.1674893211310927, 0.9810521073173698, 0.40735472249904886, 0.9358700742548033, 0.7958500907165711, 0.6046496089512394, 0.10243650581409067, 1.3628948408784771])
      # SymCubatures.setweights!(cub, T[0.0008858048750959239, 0.13640002400783297, 0.04612209214143705, 0.07405649995653696, 0.049554179974936816, 0.0938448125463987, 0.010480865652606294, 0.011622127797623866, 0.0011870772410632817, 0.08102044324683055, 0.17154667586454087])
      # cub_degree = 11
      # tol=1e-14
    elseif q <= 12 # 48 nodes
      # cub = SymCubatures.TriSymCub{T}(vertices = true,
      #                                 midedges = false,
      #                                 numS21 = 3,
      #                                 numedge = 3,
      #                                 numS111 = 3,
      #                                 centroid = false)
      # SymCubatures.setparams!(cub, T[0.3416879952650103, 0.11946629335996972, 0.5747680456926054, 0.9358700742548033, 0.7958500907165711, 0.6046496089512394, 0.38458694989344594, 0.1066659705465824, 0.6718411050524655, 0.30250431498537617, 0.7519048745331144, 0.09573484698587167])
      # SymCubatures.setweights!(cub, T[0.0012716011685488162, 0.08154407642608534, 0.03881033382161809, 0.07780339023363679, 0.006972694460371076, 0.010101834216810933, 0.011019799229776561, 0.05674688135939982, 0.08690355893881288, 0.06187386430321755])
      # cub_degree = 12
      # tol = 1e-14
      
      #48 nodes for sparse
      # cub = SymCubatures.TriSymCub{T}(vertices=true, numedge=3,
      #                                 numS21=3, numS111=3)
      # SymCubatures.setparams!(cub, T[0.5747680454804257; 0.34168799534894295;
      #                                0.1194662933444393; 0.6046496089512394;
      #                                0.7958500907165711; 0.9358700742548033;
      #                                0.38458694986468334; 0.10666597060509266;
      #                                0.3025043145989276; 0.6718411052897879;
      #                                0.09573484681006857;0.7519048745731671])
      # SymCubatures.setweights!(cub, T[0.001271601169161372; 0.07780339049594198;
      #                                 0.08154407642794283; 0.03881033381664769;
      #                                 0.01101979920649565; 0.010101834223603934;
      #                                 0.006972694458714173; 0.056746881379720164;
      #                                 0.08690355889320943; 0.06187386421672949])  
      # cub_degree = 12
      # tol = 1e-14

      cub = SymCubatures.TriSymCub{T}(vertices=true, numedge=3, numS21=3, numS111=3)
      SymCubatures.setparams!(cub, T[0.11946629335997022, 0.3416879952650081, 0.5747680456926114, 0.9358700742548033, 0.7958500907165711, 0.6046496089512394, 1.5087470795599724, 0.10666597054658086, 1.15236027848101, 0.09573484698587671, 1.025654579962154, 0.30250431498538716])
      SymCubatures.setweights!(cub, T[0.0012716011685487963, 0.03881033382161821, 0.08154407642608524, 0.07780339023362927, 0.006972694460371118, 0.010101834216810752, 0.011019799229777215, 0.05674688135939928, 0.06187386430322001, 0.0869035589388141])  
      cub_degree = 12
      tol = 1e-14
    elseif q <= 13 #55 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = true,
                                      centroid = true,
                                      numedge = 3,
                                      numS21 = 4,
                                      numS111 = 3)
      SymCubatures.setparams!(cub, T[0.3218395626601594, 0.4957183097055999, 0.9721805710603578, 0.10876332278702848, 0.94987899770573, 0.8385931397553689, 0.6815587319130891, 0.9943269168205083, 0.7196871268867262, 0.11795721403434152, 0.6250454566805667, 0.0943997083693319, 0.3344326521708524])
      SymCubatures.setweights!(cub, T[0.0007222803300461352, 0.0017482849407310453, 0.07989770555605481, 0.07525421818189708, 0.04552726413680475, 0.030615775227096503, 0.005199952283163629, 0.006832142397590008, 0.01150160780019001, 0.07355921265284922, 0.05745315943375207, 0.043511497615255, 0.11035798178531019])
      cub_degree = 13
      # # #61 nodes
      # cub = SymCubatures.TriSymCub{T}(vertices = true,
      #                                 midedges = true,
      #                                 centroid = true,
      #                                 numedge = 3,
      #                                 numS21 = 6,
      #                                 numS111 = 3)
      # SymCubatures.setparams!(cub, T[0.4757066927127499, 0.8640750851214627, 0.848071695696415, 0.09588932222312066, 0.9730573610986599, 0.2970876410785625, 0.94987899770573, 0.8385931397553689, 0.6815587319130891, 0.08772329380398057, 0.303910195608789, 0.583928734667523, 0.08412812390769837, 0.6214646072392227, 0.18555290129142354])
      # SymCubatures.setweights!(cub, T[0.0007495859766713502, 0.0004006124634786362, 0.09741666818077659, 0.018352298042798858, 0.0765136751062513, 0.024240908916246806, 0.050029973029268995, 0.06791617164641447, 0.004523675585682937, 0.00614216438323085, 0.01039041823865825, 0.038711997393405265, 0.031386407110491864, 0.05529527857696424, 0.11444067218367934])
      # cub_degree = 13
      tol=5e-15
    elseif q <= 14 # 57 nodes
      # cub = SymCubatures.TriSymCub{T}(vertices = true,
      #                                 midedges = true,
      #                                 numS21 = 5,
      #                                 numedge = 3,
      #                                 numS111 = 3,
      #                                 centroid = false)
      # SymCubatures.setparams!(cub, T[0.09123182474143035, 0.29679838270099523, 0.5601853746013142, 0.9547456706858821, 0.8556466743175067, 0.94987899770573, 0.8385931397553689, 0.6815587319130891, 0.09449105513358476, 0.29883210637720314, 0.28045625536995294, 0.5626405122741717, 0.08624428771540175, 0.6026786176032012])
      # SymCubatures.setweights!(cub, T[0.0007538778164594649, 0.009414839808665405, 0.023018612111280527, 0.058120471428513544, 0.09454606798115053, 0.05373738699005786, 0.0727513309290953, 0.004178640246830719, 0.0072314035856924815, 0.008184710625980713, 0.040253997904260874, 0.06900626592215424, 0.04830702151580295])
      # cub_degree = 14
      # tol = 1e-14
      #61 nodes for sparse
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = true,
                                      numS21 = 6,
                                      numedge = 3,
                                      numS111 = 3,
                                      centroid = true)
      SymCubatures.setparams!(cub, T[0.9540069927604912, 0.09149922166977983, 0.5272375124699329, 0.8540637650586674, 0.29746375857253127, 0.7370347878007151, 0.94987899770573, 0.8385931397553689, 0.6815587319130891, 0.08118879249452095, 0.6045976879703241, 0.2993309303393121, 0.09452178262313014, 0.5659292908830563, 0.2642822989570226])
      SymCubatures.setweights!(cub, T[0.0007544183587686153, 0.009618517048991415, 0.054301359425760844, 0.02311383313585573, 0.07479536196607821, 0.0724047539507651, 0.058530504993888254, 0.023085894583599966, 0.004194700437036379, 0.007216774672018443, 0.007658144705588096, 0.04571774280785687, 0.040342997894998796, 0.06569686953205542, 0.025222689311551826])
      cub_degree = 14
      tol = 1e-14
    elseif q <= 15 #69 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = false,
                                      centroid = false,
                                      numedge = 4,
                                      numS21 = 6,
                                      numS111 = 4)
      SymCubatures.setparams!(cub, T[0.453639886926089, 0.757300794432156, 0.9833658031003171, 0.2673140133674513, 0.09859868374369533, 0.9313655727780035, 0.9597669540832294, 0.8693869325527526, 0.7389624749052223, 0.5826394788331934, 0.05491312392961762, 0.6373866699608978, 0.9962538574472307, 0.7237541325698926, 0.19371136334839004, 0.543875090517202, 0.31621286243906677, 0.07653701224015233])
      SymCubatures.setweights!(cub, T[0.00041188664492458155, 0.07744252741964325, 0.07978483547877466, 0.022442598130765025, 0.055748397969630784, 0.02612364913831872, 0.047345838243333084, 0.00378591167322386, 0.00494107871623112, 0.005183734332323488, 0.002737379819615485, 0.03279819517102539, 0.041923249456950765, 0.05253038629189004, 0.03478353135937816])
      # #78 nodes
      # cub = SymCubatures.TriSymCub{T}(vertices = true,
      #                                 midedges = false,
      #                                 centroid = false,
      #                                 numedge = 4,
      #                                 numS21 = 7,
      #                                 numS111 = 5)
      # SymCubatures.setparams!(cub, T[0.4077750904687781, 0.7603167059782701, 0.09923314655025814, 0.5392503479960836, 0.7336297585443722, 0.8840881622927267, 0.9647641607531148, 0.9597669540832294, 0.8693869325527526, 0.7389624749052223, 0.5826394788331934, 0.5906777979138421, 0.22187382664358296, 0.4240838353085315, 0.9561037136275863, 0.32312514688527244, 0.1922339643728221, 0.3059644597113809, 0.06164364302370151, 0.61117935281228, 0.06810139957582922])
      # SymCubatures.setweights!(cub, T[0.00038560167936041434, 0.051828320319575245, 0.04897897079965172, 0.025976977006482195, 0.013736493536953337, 0.008950065388866491, 0.060619930406823286, 0.04266756023023606, 0.0038903087243377906, 0.003952860653539124, 0.005383423214263016, 0.006285453618608518, 0.0547038473668017, 0.028976866620338368, 0.039091444929817185, 0.026295287780152696, 0.03818188074150067])
      # SymCubatures.setparams!(cub, T[0.926073884304195, 0.6324214560155823, 0.2732070265874975, 0.9800392449549591, 0.9275624142028438, 0.09803698045583892, 0.4678434864713421, 0.9597669540832294, 0.8693869325527526, 0.7389624749052223, 0.5826394788331934, 0.6281478611396281, 0.05988739351469942, 0.5609548008743171, 0.2126550295547661, 0.31652535690496086, 0.08109010595545985, 0.8143850982384193, 0.7648461377280185, 0.932544753465994, 0.30649918776292634])
      # SymCubatures.setweights!(cub, T[0.00043899617260247656, 0.02832593621352893, 0.0284355895772454, 0.05708633451373208, 0.023936844515619497, 0.022118700922270935, 0.02609816009431366, 0.07746633875563608, 0.0036965993872807804, 0.005287955943337329, 0.005332396620054094, 0.00394995932661379, 0.03487470485376105, 0.06094440138731721, 0.03641363868298402, 0.024295609746207345, 0.026584617003303227])
      tol=5e-15
    elseif q <= 16 # 72 nodes
      # cub = SymCubatures.TriSymCub{T}(vertices = true,
      #                                 midedges = false,
      #                                 numS21 = 5,
      #                                 numedge = 4,
      #                                 numS111 = 5,
      #                                 centroid = false)
      # SymCubatures.setparams!(cub, T[0.06522886751875544, 0.22580096183071185, 0.7596981294914541, 0.4709436867165254, 0.9626832123175285, 0.9597669540832294, 0.8693869325527526, 0.7389624749052223, 0.5826394788331934, 0.6837234484596031, 0.07332394938311322, 0.0696195668842623, 0.212933356067202, 0.23253598278247095, 0.45142512187590106, 0.07099634477078105, 0.4271102931595644, 0.7291817098161258, 0.24027678681212955])
      # SymCubatures.setweights!(cub, T[0.00046398080499042475, 0.011683737917833769, 0.03748058638601257, 0.08029394321068659, 0.07225832122458112, 0.0350040593727091, 0.002438057356674435, 0.004350159403927255, 0.005775046495650795, 0.006715480702033995, 0.03300813082784586, 0.021283575846091102, 0.05161492134502783, 0.028186978295003252, 0.061368668602672])
      # cub_degree = 16
      # tol = 1e-14
      #75 nodes
      # cub = SymCubatures.TriSymCub{Float64}(vertices=true, numedge=4,
      #                                       numS21=4, numS111=6)
      # SymCubatures.setparams!(cub, T[0.0768946752469594; 0.6109907336234316; 0.4179369130153705;
      #         0.23221732669622028; 0.5826394788331936; 0.7389624749052223;
      #         0.8693869325527526; 0.9597669540832294; 0.930330150896981;
      #         0.6679157686119799; 0.8075697058065031; 1.1286554029515519;
      #         0.07116357646128006; 0.25084609708056116; 0.20705942876211147;
      #         0.46901560437791967; 0.06341212586405608; 0.504346239131436;
      #         0.7497163430497378; 1.0452430326898021])
      # SymCubatures.setweights!(cub, T[0.0004934174763938973; 0.016280379328615396; 0.03461636692031435;
      #         0.0544553219641883; 0.04016179760262383; 0.005790378245871778;
      #         0.005099447551845056; 0.004407182316786121; 0.002824162676005338;
      #         0.05536224506473959; 0.03347658309313654; 0.025582552500936707;
      #         0.047069550256364515; 0.029816087172401348; 0.050901502809179516])
      # cub_degree = 16
      # tol = 1e-14
      cub = SymCubatures.TriSymCub{Float64}(vertices=true, numedge=4, numS21=4, numS111=6)
      SymCubatures.setparams!(cub, T[0.07695246604559251, 0.23155195440690438, 0.4214410348442252, 0.5984805733893661, 0.9597669540832294, 0.8693869325527526, 0.7389624749052223, 0.5826394788331934, 1.6780802584948735, 0.070868673270381, 1.4310000152356053, 0.06475029714620971, 1.1307663966420336, 0.059871502734567955, 1.3223892173002267, 0.21047585794328214, 1.0527058565418306, 0.19351851169866596, 0.9390386879493979, 0.38284849575318963])
      SymCubatures.setweights!(cub, T[0.0004936028234827967, 0.016304902766397903, 0.04002920177118074, 0.05424925987527771, 0.04433748243927064, 0.002825931524738883, 0.004388548114996859, 0.005219115996897478, 0.00542370270562261, 0.025500278386860246, 0.030360271992401475, 0.0315343758493793, 0.047411491325143024, 0.04896454269626057, 0.053997849903228125])
      cub_degree = 16
      tol = 1e-14
    elseif q<=17 #78 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = true,
                                      centroid = false,
                                      numedge = 4,
                                      numS21 = 6,
                                      numS111 = 5)
      SymCubatures.setparams!(cub, T[0.9609712932934811, 0.041766468223703515, 0.28302191880327227, 0.760185750106773, 0.46292528148197715, 0.173093973905005, 0.9670007152040296, 0.8922417368315723, 0.7826176634981026, 0.6478790677934697, 0.1610470535559693, 0.04528225579552883, 0.7226288131709319, 0.2400281193275206, 0.36438411202861126, 0.07009663480762927, 0.06638724199110738, 0.6541614471267511, 0.20944043731114892, 0.46300947657489105])
      SymCubatures.setweights!(cub, T[0.00028908921124315597, 0.006749278476735112, 0.04010182067917665, 0.006389003512011476, 0.031350458689112376, 0.08284141010739084, 0.07608792396452779, 0.028383320820692025, 0.0008879784296625767, 0.0025339572186700047, 0.005105107386578676, 0.004945455787338109, 0.013263503268280243, 0.06280733563151443, 0.028782576708273472, 0.0346732468195029, 0.04423801935306822])
      #81 nodes
      # cub = SymCubatures.TriSymCub{T}(vertices = true,
      #                                 midedges = true,
      #                                 centroid = false,
      #                                 numedge = 4,
      #                                 numS21 = 7,
      #                                 numS111 = 5)
      # SymCubatures.setparams!(cub, T[0.9593342591248476, 0.055020860803980914, 0.28526621062594437, 0.8387490974325792, 0.7575470957605024, 0.4643330654690022, 0.16444412944990114, 0.9670007152040296, 0.8922417368315723, 0.7826176634981026, 0.6478790677934697, 0.18880867265986118, 0.035555254484675405, 0.7184784536560056, 0.23691797352429475, 0.3741094228156754, 0.07603572398114046, 0.062436028204803366, 0.6623692442916642, 0.20752925425214702, 0.4787260323140681])
      # SymCubatures.setweights!(cub, T[0.00039358697853621074, 0.0071635142083468235, 0.04019813807428178, 0.009282438655081281, 0.03869786195664633, 0.011258459094176387, 0.07803814033007511, 0.07672709475642626, 0.02944041155775245, 0.0013436036287747311, 0.0013462041905516386, 0.0055617515345457525, 0.004323946528663438, 0.01238657052502277, 0.05821065906057333, 0.029593585323296128, 0.03325710394842388, 0.04171008578782039])
      cub_degree = 17
      tol = 1e-14
    elseif q <= 18 # 93 nodes
      # cub = SymCubatures.TriSymCub{T}(vertices = true,
      #                                 midedges = true,
      #                                 numS21 = 5,
      #                                 numedge = 4,
      #                                 numS111 = 8,
      #                                 centroid = false)
      # SymCubatures.setparams!(cub, T[0.9707518160195144, 0.20641902821770264, 0.05909234695936303, 0.7110519489326181, 0.39622899824483554, 0.9670007152040296, 0.8922417368315723, 0.7826176634981026, 0.6478790677934697, 0.3906058825830485, 0.19331653245290645, 0.8008189404907055, 0.17109449910705912, 0.24603845950256573, 0.9650206949829071, 0.20142808520129168, 0.5868483166437959, 0.4088076515219679, 0.6463428792626565, 0.4093195068370922, 0.05981030752842776, 0.1973673755323838, 0.06444759127052886, 0.673801874205518, 0.06033239562749046])
      # SymCubatures.setweights!(cub, T[0.00031737159521531795, 0.00491592992346147, 0.02867493867862121, 0.029107246299861025, 0.009795965466384334, 0.031357172956706196, 0.05321207918099946, 0.0017927521866819998, 0.0033404931412304374, 0.004059862113462083, 0.004795116727932584, 0.03247349913877226, 0.022030499945312042, 0.0192471893187778, 0.03594121829355437, 0.05970283389294487, 0.02407056230858391, 0.018810685588525493, 0.028378268626930925])
      # cub_degree = 18
      # tol = 1e-14

      #91 for sparse 
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = true,
                                      numS21 = 8,
                                      numedge = 4,
                                      numS111 = 6,
                                      centroid = true)
      SymCubatures.setparams!(cub, T[0.06262408792230445, 0.19777933032042158, 0.37882457759149224, 0.5522678470067186, 0.7417773969871918, 0.8384879881534788, 0.9192489794614319, 0.9751369967510601, 0.9670007152040296, 0.8922417368315723, 0.7826176634981026, 0.6478790677934697, 0.3969247802147992, 0.19029720152483445, 0.17280786654060148, 0.6444759719518841, 0.34728044663382485, 0.6010071331698775, 0.4182247647898829, 0.05843657476537568, 0.20577249022713864, 0.06088553891723073, 0.6831211550516773, 0.05282478749715142])
      SymCubatures.setweights!(cub, T[0.0003297824717475295, 0.004124392067225074, 0.01087901972194179, 0.029045091983495044, 0.04316393078212389, 0.03924781094944372, 0.034073175803525244, 0.04207173428673133, 0.0383678047460057, 0.02433697237065314, 0.0018849300430363557, 0.0031445820532631672, 0.003988033707787848, 0.004163737246541059, 0.036726442736361095, 0.0391750110870636, 0.043924272917320734, 0.023405456462375726, 0.018160144500102976, 0.02473059563066814, 0.007261616144198268])
      cub_degree = 18
      tol = 1e-14
    elseif q<=19 #96 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = false,
                                      centroid = false,
                                      numedge = 5,
                                      numS21 = 5,
                                      numS111 = 8)
      SymCubatures.setparams!(cub, T[0.584118160705067, 0.05594736023975348, 0.35206685470664284, 0.8193250037475288, 0.9106302368288487, 0.9724496361114411, 0.9096396608220033, 0.8164380765159303, 0.6997654704826745, 0.5682764664274638, 0.36516822326597015, 0.05300184234489419, 0.39755390108675603, 0.17428704397839728, 0.05879006206078798, 0.1815645896189941, 0.17519871987759317, 0.2119731854287182, 0.5725515218226456, 0.35807920419385625, 0.054591434526018896, 0.8428999975969531, 0.17752640292448835, 0.6418878905694481, 0.5915035416317636, 0.05407771456049964])
      SymCubatures.setweights!(cub, T[0.0002316646038971307, 0.058241783304694345, 0.008602602423026792, 0.04069174280772659, 0.05226677627664586, 0.042505709104633006, 0.0014109503198198763, 0.00256960026599889, 0.00309831866176889, 0.003741487391936303, 0.00407918767311633, 0.018331624199230855, 0.034095002356800554, 0.015212277721266073, 0.015179420870583436, 0.04887476434442788, 0.023431653507726604, 0.04029884352862801, 0.021740063231717795])
      
      #99 nodes
      # cub = SymCubatures.TriSymCub{T}(vertices = true,
      #                                 midedges = false,
      #                                 centroid = false,
      #                                 numedge = 5,
      #                                 numS21 = 6,
      #                                 numS111 = 8)
      # SymCubatures.setparams!(cub, T[0.3358980598350944, 0.9789706742159475, 0.9289390131831053, 0.7440240615283973, 0.06501518566521776, 0.19142118573089564, 0.9724496361114411, 0.9096396608220033, 0.8164380765159303, 0.6997654704826745, 0.5682764664274638, 0.1583978333661437, 0.3969064666813339, 0.042212410750339456, 0.6904634855900522, 0.42850251199567174, 0.04755503506287861, 0.9879382272825736, 0.46840576395362366, 0.9475849438753418, 0.29884780482324186, 0.6542434414374048, 0.14303581650779246, 0.5459090121405945, 0.29780285727201, 0.21229041305026342, 0.05731545867768344])
      # SymCubatures.setweights!(cub, T[0.00021393564288845607, 0.03852805201142253, 0.020863457513064957, 0.03579097887221835, 0.055167566276037854, 0.01162155336267164, 0.02853211993973326, 0.001663406453602736, 0.0025528949149080204, 0.002946876800429615, 0.0029333115372134344, 0.003080267714905865, 0.03252630590418761, 0.019886495371279594, 0.019158063559635667, 0.028048303639912196, 0.036293762745006924, 0.034722432073318175, 0.03654424141681496, 0.01761813939309992])
      cub_degree = 19
      tol = 5e-15
    elseif q <= 20 # 103 nodes
      # cub = SymCubatures.TriSymCub{T}(vertices = true,
      #                                 midedges = false,
      #                                 numS21 = 5,
      #                                 numedge = 5,
      #                                 numS111 = 9,
      #                                 centroid = true)
      # SymCubatures.setparams!(cub, T[0.5706569630167525, 0.045476786416401474, 0.3386020213960199, 0.8173378833539444, 0.921990277619815, 0.9724496361114411, 0.9096396608220033, 0.8164380765159303, 0.6997654704826745, 0.5682764664274638, 0.33734350166550275, 0.051994927730595886, 0.39122487789147525, 0.17029307283717016, 0.0555956144787866, 0.15691138175383615, 0.6048213057399955, 0.16379237699452232, 0.1617644050989192, 0.21607595265283405, 0.56000817972725, 0.34672951075674346, 0.053470945381353656, 0.8336640057658877, 0.19072411914891063, 0.7807442886366605, 0.5689238612251344, 0.05028584743589604])
      # SymCubatures.setweights!(cub, T[0.00020287039679960927, 0.05480860557627294, 0.0059385235222309644, 0.03842886096227663, 0.05547473001542249, 0.0169887933138356, 0.0011773441836911722, 0.0024488282936603527, 0.0030066460137970824, 0.0034778130522919964, 0.004015524716182081, 0.018032136009343752, 0.030456546263333304, 0.013428745082940518, 0.030846728655919257, 0.015872166479982706, 0.0501088680789255, 0.02445659700889463, 0.025248076325159272, 0.02109467351877063, 0.02244868654213111])
      # cub_degree = 20
      # tol = 1e-14

      #108 nodes for sparse sbp
      # cub = SymCubatures.TriSymCub{T}(vertices = true,
      #                                 midedges = false,
      #                                 numS21 = 5,
      #                                 numedge = 5,
      #                                 numS111 = 10,
      #                                 centroid = false)
      # SymCubatures.setparams!(cub, T[0.05294490986531675, 0.16773289449842865, 0.31622279360326333, 0.4879667756898061, 0.6513493837496125, 0.9724496361114411, 0.9096396608220033, 0.8164380765159303, 0.6997654704826745, 0.5682764664274638, 0.051195998723422086, 0.17426772855417563, 0.15559719406273465, 0.3410916700392985, 0.047462163492437585, 0.3558632464690134, 0.3015188695888322, 0.5149248126378921, 0.1485561447226546, 0.5580518026197324, 0.045282999147861325, 0.5845980081189597, 0.47766585164678715, 0.6759879691102543, 0.2984203296477869, 0.7360991066041834, 0.14866968526150143, 0.8002558630925941, 0.04562318536898432, 0.8428283926597013])
      # SymCubatures.setweights!(cub, T[0.0002292509138976338, 0.007783575266946965, 0.021386288997735533, 0.03282340926195889, 0.03737509863710259, 0.01232514825527602, 0.0013325681607189903, 0.0022223561249509146, 0.002768972623204779, 0.0031274574625577323, 0.0034137825308797613, 0.013023905490374535, 0.02669382193873107, 0.016372911768879273, 0.036527379608411736, 0.0301333208225669, 0.01857476823754868, 0.032598392548470634, 0.03828529775483018, 0.032083926453061136, 0.02021308614168822])
      # cub_degree = 20
      # tol = 5e-14

      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = false,
                                      numS21 = 5,
                                      numedge = 5,
                                      numS111 = 10,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.05310181454886258, 0.16669909638072203, 0.31890625932107397, 0.4769062635958667, 0.6122291577074322, 0.9724496361114411, 0.9096396608220033, 0.8164380765159303, 0.6997654704826745, 0.5682764664274638, 1.7744542331519155, 0.05080126431977686, 1.5955445904517396, 0.04801899242996552, 1.3704216313092366, 0.043836148777747434, 1.114531215089879, 0.03979079119920485, 1.5034294955202192, 0.15742672682529332, 1.2985603401515144, 0.14404642861985076, 1.062482496970096, 0.13060240593213884, 1.186451789763432, 0.29251018581533667, 0.9830892364157693, 0.2651997011457629, 0.8833988726389188, 0.43207240512538564])
      SymCubatures.setweights!(cub, T[0.0002294787365987622, 0.007828268985620133, 0.021195156339203238, 0.03320663292184484, 0.03691909632788907, 0.029040086649929967, 0.0013362233467532185, 0.002205307837565235, 0.002803342974009274, 0.003023258919430133, 0.0029678334323056895, 0.012953389146554724, 0.01656716962783203, 0.01799963185402161, 0.017696124066043997, 0.02688605319231424, 0.029377664531438242, 0.029006832813656704, 0.03575222698099652, 0.035239251856487463, 0.035309662773381226])
      cub_degree = 20
      tol = 5e-14

    elseif q <= 22
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = true,
                                      numS21 = 10,
                                      numedge = 5,
                                      numS111 = 10,
                                      centroid = true)
      SymCubatures.setparams!(cub, T[0.0452671188791389, 0.14251445923625328, 0.28073020926730075, 0.42127731315299516, 0.558591685525777, 0.9808784475256531, 0.941755226640668, 0.8924521179816877, 0.8344191924303693, 0.7628077735811749, 0.976654923321082, 0.9231737823259362, 0.8430942345408787, 0.741454910545668, 0.62464346505312, 1.80711191161222, 0.04340190963100611, 1.652396972843757, 0.04108975716835016, 1.453269828612064, 0.04044688560715337, 1.2280021473835083, 0.033797064982277754, 1.5726788533048874, 0.13666382887880732, 1.390247092801606, 0.12988285384097822, 1.179582780981497, 0.11621906093358203, 1.2785482904192165, 0.2569247204326685, 1.0899451542969714, 0.24050646929898697, 0.9970083540636144, 0.3906008774317397])
      SymCubatures.setweights!(cub, T[0.00016426637898436108, 0.0027265597361666068, 0.005705761130911468, 0.015544790593692227, 0.027151260034983787, 0.03113815956460653, 0.029416408312620403, 0.015222191200137928, 0.02197699035855633, 0.022810347487063024, 0.025560959844254653, 0.029261240672412963, 0.0009654224796097691, 0.0016101346598994335, 0.0020596785887459393, 0.0024838239054362035, 0.002187871336388239, 0.009517122805401587, 0.012400968731461239, 0.014439412299382478, 0.013900677732350162, 0.020563237305462546, 0.022675243107811598, 0.02455287706085823, 0.028010465565942614, 0.02908858712160794, 0.0300164778576649, 0.03313119070869239])
      cub_degree = 22
      tol = 5e-14
    elseif q <= 24
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = false,
                                      numS21 = 6,
                                      numedge = 6,
                                      numS111 = 15,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.03880020746539761, 0.12469808500950672, 0.24546909601117012, 0.3824152837478771, 0.5133193869937424, 0.6213749471902272, 0.9799675226336304, 0.9339005269151737, 0.8644342995456631, 0.7753197014643236, 0.6713620066713564, 0.5581659344418519, 1.8338443095729768, 0.03780922263195278, 1.6994790270415137, 0.03640951307210373, 1.526293098464035, 0.03429184874986117, 1.3234587981554498, 0.031190632422757072, 1.1015360584910252, 0.027579295831113855, 1.6242840516199397, 0.12000396109006645, 1.4624546933179234, 0.11317147001885593, 1.2735879053298942, 0.10322501854397217, 1.0653951368868222, 0.09170711961572527, 1.3626513277391037, 0.2316819648099439, 1.194470962768185, 0.21217173127961658, 1.0086952027493412, 0.189774133359083, 1.0911070815110906, 0.35183079053218175, 0.9340684516928085, 0.3169984087390115, 0.8466357861808691, 0.46557399075669426])
      SymCubatures.setweights!(cub, T[0.00012069840580176264, 0.004198664879865943, 0.012064586107175031, 0.020632562484240925, 0.026184605800877482, 0.026332465393874897, 0.020583953115869035, 0.0007103036654031911, 0.0012087324073259897, 0.0015971192723511888, 0.0018297705187362292, 0.0018697848709738776, 0.0017448419361082413, 0.007158296860976668, 0.00951187093723075, 0.010972580906638338, 0.011288927275732316, 0.010602982512680697, 0.015896277087086595, 0.018314096666680088, 0.018976980128332852, 0.01808891876434126, 0.023474607350174835, 0.024286758745376905, 0.02345506450457874, 0.026695959163965766, 0.02578249016436156, 0.024808201500424714])
      cub_degree = 24
      tol = 5e-14
    elseif q <= 26
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = true,
                                      numS21 = 12,
                                      numedge = 6,
                                      numS111 = 15,
                                      centroid = true)
      SymCubatures.setparams!(cub, T[0.03393923611308314, 0.10873649956316811, 0.21482988698982655, 0.3398694855014765, 0.4633755619477792, 0.5765169904802263, 0.9866066659382519, 0.957310865314053, 0.9167523374350691, 0.8686175130656929, 0.8131598018284998, 0.7482527988848714, 0.9826229632519192, 0.9425410221114882, 0.8817598449759076, 0.8031266027349229, 0.7103190273568363, 0.6076769776818971, 1.8547247311426134, 0.03286746714508377, 1.7366818435509916, 0.031520503014385394, 1.5832211302981727, 0.029980899742798908, 1.4010761443451814, 0.028752652405317998, 1.2005374709609065, 0.025115089079989942, 1.6716587080221865, 0.10434769233046208, 1.5270795860587258, 0.09920116448130761, 1.354249854084208, 0.09440983064095608, 1.1644201761813233, 0.08468894753424339, 1.4380112429828442, 0.20422577449383864, 1.280888624803641, 0.19215819989412208, 1.1065564041141398, 0.17677134397464672, 1.1837996453343822, 0.3161530200288055, 1.0291162067290274, 0.29455381556320365, 0.9443063768751075, 0.4286992022860749])
      SymCubatures.setweights!(cub, T[9.069238190165345e-5, 0.0016179539521252906, 0.0032155845831953237, 0.009250891464077097, 0.01620659319694568, 0.021767107373752377, 0.023263766820218256, 0.02152640541290746, 0.009449584793961761, 0.014900000545125563, 0.017691308186816863, 0.018806596685280074, 0.01979181103155395, 0.021097197260307152, 0.0005392245829151473, 0.0009152304716258027, 0.0012124371481207756, 0.0014189727760605922, 0.0015564220064443746, 0.0014406670450625991, 0.005466568965680323, 0.007271477201919237, 0.008524845972745438, 0.009307932744515836, 0.008892154836773757, 0.012287607557515269, 0.014406327743146069, 0.015519882011292993, 0.015562805983626506, 0.018878712663747736, 0.01981605153748861, 0.020451592586636647, 0.022316281638671734, 0.02226915826729402, 0.022126730269780764, 0.022887014869109382])
      cub_degree = 26
      tol = 5e-14
    elseif q <= 28
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = false,
                                      numS21 = 7,
                                      numedge = 7,
                                      numS111 = 21,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.0297684342510763, 0.0955990655973306, 0.19291754638860226, 0.3071372298707168, 0.42368986691738036, 0.5320960490073281, 0.6246702357765312, 0.984784023135109, 0.949600266546736, 0.8960041459309076, 0.8261943514412465, 0.7430297109435688, 0.6499152344503816, 0.5506631367609747, 1.872342032026377, 0.028871097803943917, 1.7675811074108279, 0.028158691401833948, 1.6306934015947248, 0.02719784703705251, 1.4677125664505128, 0.025321800402231868, 1.2846879400321705, 0.02347760065222769, 1.0891749012037502, 0.022023861231761973, 1.7093628889667856, 0.09344851689576034, 1.5801489061325036, 0.08990447329553126, 1.4262344614567835, 0.08409466147453422, 1.2515984975183934, 0.07779747113489333, 1.0624982098653115, 0.07221835901759747, 1.4931838100842858, 0.18443446862241503, 1.3509572181334975, 0.1729159926865194, 1.191532901609541, 0.1600935095819077, 1.0192254255841438, 0.14693939098728045, 1.2587688670743364, 0.2877527514030317, 1.1161287832133815, 0.2667227445499398, 0.963965368480241, 0.2427509005204612, 1.0307078228125615, 0.39285142735909817, 0.8993858298065905, 0.35757238978691713, 0.826436434580063, 0.4877366783995886])
      SymCubatures.setweights!(cub, T[6.94309976584969e-5, 0.0024780002921293448, 0.007172658254169629, 0.013254435516353235, 0.018012182717792716, 0.020052958218257165, 0.019171216903482207, 0.01687828214819848, 0.000414121332183576, 0.0007067697110965416, 0.0009570685938221425, 0.0011484008367847046, 0.0012290067714772204, 0.0012459481505339895, 0.001223584842576201, 0.004233379821013063, 0.005763255725052125, 0.0069045787546334475, 0.007438572872977677, 0.007516612526806014, 0.007316609968826768, 0.009766599872356634, 0.011620907020761661, 0.012683811680144989, 0.012905579101697942, 0.012404455913928735, 0.01546497039525781, 0.016594719193546714, 0.016757909866986224, 0.015839460335576226, 0.019156203279994795, 0.019186851712129985, 0.01799924291577809, 0.019989076700119, 0.01927921509390719, 0.019041837819343238])
      cub_degree = 28
      tol = 5e-14
    elseif q <= 30
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = true,
                                      numS21 = 14,
                                      numedge = 7,
                                      numS111 = 21,
                                      centroid = true)
      SymCubatures.setparams!(cub, T[0.026263712933993107, 0.08517773895617918, 0.17279615006913418, 0.27589945905594404, 0.38841435055824713, 0.4974390314519614, 0.5848434660310842, 0.9888569686145411, 0.9661320178091347, 0.939352201678907, 0.9019197972335102, 0.8491616753371906, 0.7917391940597742, 0.73363931848329, 0.9865660883157092, 0.9554399979577868, 0.9078481256108851, 0.8455144903138423, 0.7706926996650507, 0.6860872167827385, 0.5947559867591588, 1.8870221434136396, 0.02570366713733081, 1.793983492536147, 0.025119958538073715, 1.671668422695357, 0.02432310058038001, 1.5249412598428131, 0.02255867065673869, 1.3578052693232254, 0.021693571868118818, 1.1783277514824677, 0.019013639177892776, 1.7405684419522667, 0.0834737046514834, 1.6245714394071344, 0.08045190568031999, 1.48581264880407, 0.07531756667200452, 1.3259706249257068, 0.07132176062782467, 1.152384095100876, 0.06529596804242158, 1.544801702685509, 0.1652266243762263, 1.4137371720619094, 0.15610982592370423, 1.2659378667488168, 0.14566531636274646, 1.1041469369299484, 0.13811994187867957, 1.329023725898816, 0.26259117857381364, 1.1969481364030816, 0.24335464892304964, 1.0461028692407564, 0.23226408985909297, 1.108279879880839, 0.3623302386259286, 0.9830644472062909, 0.3357220260332416, 0.9068645543524504, 0.4456950580906604])
      SymCubatures.setweights!(cub, T[5.405044786504608e-5, 0.001192880437707203, 0.0019312314810255172, 0.005711073426261259, 0.010761257315153134, 0.014756398945140059, 0.0174157684168835, 0.01864510095651114, 0.017963788108972353, 0.006825296184896629, 0.00955374782283768, 0.010427681661980178, 0.015375884367431562, 0.01793569875819078, 0.014889786796888176, 0.01574799038977883, 0.0003225828816021466, 0.0005572258909639724, 0.0007594058158235478, 0.0009194823497034421, 0.0009865450398853174, 0.0010564102265268403, 0.0009516861253625811, 0.003338566145875818, 0.004578739868326248, 0.005538718036858022, 0.006014310823223363, 0.006334876293383455, 0.006020570758806448, 0.007827674900678702, 0.009358253728476271, 0.010361418602777934, 0.010715785443984336, 0.010965016992870157, 0.012640541607047326, 0.013979569962459408, 0.013996330146940432, 0.014536435857847646, 0.01637096964300807, 0.016709472707133105, 0.01633162255120627, 0.018049100713842133, 0.01549868127417652, 0.01636315816154465, 0.015938178145424926])
      cub_degree = 30
      tol = 5e-14
    elseif q <= 32
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = false,
                                      numS21 = 8,
                                      numedge = 8,
                                      numS111 = 28,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.02321871305376844, 0.07661618717846833, 0.1559791421368959, 0.2521738623055426, 0.3545692287061555, 0.45767816195407424, 0.549463769469837, 0.6306162242106939, 0.9880527787060993, 0.9603245926737669, 0.917796767609045, 0.8618396646416213, 0.7942524171593308, 0.717207518456062, 0.6331813264391405, 0.544874546742326, 1.8996496108572756, 0.023092166887234995, 1.81676414176096, 0.022715838377297376, 1.7074073590386896, 0.02190016787386545, 1.5748173470551243, 0.02064368665731146, 1.422688000649648, 0.019585055825585337, 1.2561048718941046, 0.01861624970654637, 1.0808477405119028, 0.0172820077297201, 1.766051495211983, 0.07538386482064902, 1.6609851555541308, 0.07267709679349796, 1.5346215483282108, 0.06852877906576661, 1.388707690087819, 0.06506031014508493, 1.2284430082860829, 0.061407475271634704, 1.0591012767277843, 0.05706415698149338, 1.5879766214823252, 0.15025750903992952, 1.4691418539032406, 0.14168359490021026, 1.3318603782356464, 0.13452633663066355, 1.1827648341624186, 0.12591884275049853, 1.0247834740282877, 0.117066670587145, 1.3868287877740464, 0.23790955375732145, 1.2603790725236477, 0.225354316669128, 1.1248510950616952, 0.209547690202714, 0.9807770173315405, 0.1944139815353357, 1.1785158683308932, 0.33441715076647405, 1.057898055579627, 0.3099719835521217, 0.9289685153380944, 0.2868439182948843, 0.9825174624912336, 0.42421743735884065, 0.870007736902752, 0.3930708327461422, 0.8057659998170382, 0.5102884711042318])
      SymCubatures.setweights!(cub, T[4.268165075517583e-5, 0.001511180896766916, 0.004639497480618228, 0.008767080714884294, 0.012645153794488592, 0.015290856591938302, 0.015800848754198134, 0.015020041334806205, 0.012642495139653596, 0.00025363281773869924, 0.0004462611165941801, 0.0006149328457454014, 0.0007431700207732406, 0.0008199857920184081, 0.000866651267053437, 0.0008845297431302931, 0.0008479789298755368, 0.0026627010371127985, 0.0036907517806568247, 0.0044879434346930406, 0.004976043679418772, 0.0052754219943729926, 0.005358653756577263, 0.005148561081254316, 0.006386540567391889, 0.007715605263605968, 0.008537788031568457, 0.009068569882674191, 0.009113223323352715, 0.008799404660587205, 0.010526267347682484, 0.011547663258573355, 0.012111889546242195, 0.011942573836421836, 0.011496055006004618, 0.0138721095378295, 0.014316784715312027, 0.013939070959970315, 0.013258350321481118, 0.015647895583310877, 0.015223884592364727, 0.01440793565364018, 0.015511248676234423, 0.014973429781573644, 0.01467990531144163])
      cub_degree = 32
      tol = 5e-14
    elseif q <= 34
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = true,
                                      numS21 = 16,
                                      numedge = 8,
                                      numS111 = 28,
                                      centroid = true)
      SymCubatures.setparams!(cub, T[0.02087996297483963, 0.0687657959397549, 0.1400942075337306, 0.22904399887821697, 0.3244246918508021, 0.4225035034208078, 0.5164179478612008, 0.5968962009656194, 0.9909754186855135, 0.9714893060655033, 0.9467140834714922, 0.9187810654688159, 0.8798110837026452, 0.8314878536112061, 0.7787816679235553, 0.726123093532526, 0.98930588311104, 0.9644507640762932, 0.926230288898323, 0.8757471012763065, 0.8144540686326103, 0.7441146428403568, 0.6667524239122493, 0.5845930117046407, 1.9097809206961698, 0.020701947099993952, 1.8350983922783348, 0.020317844321620204, 1.7360431634277984, 0.019794122468429938, 1.6155468191667306, 0.018796019451068163, 1.4767619143526631, 0.01759796701910264, 1.3233048593378425, 0.01693455896928361, 1.1611521286328597, 0.014501061604756487, 1.7899081894668152, 0.06752946481275546, 1.6942579383447123, 0.06579642432776667, 1.5792257729511674, 0.062356861707573014, 1.4461675962696445, 0.058864788836341886, 1.2981244140858863, 0.055880025454212334, 1.1411760531336008, 0.050078261953227406, 1.6277875980656276, 0.13631496767508472, 1.5187331981980707, 0.12891768995190417, 1.3918465629393817, 0.12246242971781919, 1.2525337540363704, 0.11516622190496596, 1.105052441785466, 0.10633720288126963, 1.441741186638359, 0.21687345876510444, 1.3245514845186277, 0.20616666989104557, 1.1961373214124558, 0.1940723697547783, 1.059271392227268, 0.18054317027492742, 1.2453601669765948, 0.30691184233155966, 1.1274913457938105, 0.2904222097435548, 1.0053609064657212, 0.2690501205016952, 1.0512083267736851, 0.39935883483120954, 0.9434781388674804, 0.3693737003918547, 0.8785179611793953, 0.4785347431708478])
      SymCubatures.setweights!(cub, T[3.4179821416342116e-5, 0.0008534354553188183, 0.0012226458876656932, 0.003749346858448903, 0.0071470429085852445, 0.010517413107750343, 0.01318313983195964, 0.014265826959449116, 0.01346103442025923, 0.013048297319309163, 0.005010578107677077, 0.007752751435923993, 0.008051980644171388, 0.009939078236524944, 0.012795256673525275, 0.013527636875234963, 0.012809075067295467, 0.012042100560377649, 0.00020419654721598002, 0.0003587830321044669, 0.0004949337983951107, 0.0006073518535294149, 0.0006797163208568804, 0.0007110744450001241, 0.000745757999800162, 0.0006456370941810816, 0.002150867325788182, 0.0029812262722087772, 0.0036757406705631625, 0.0041224060587171036, 0.004360236976258753, 0.004501227974834951, 0.004144125325497679, 0.005172199839721902, 0.00634069796090139, 0.007068495787636817, 0.007589217011283014, 0.007712606390084903, 0.007607124682732092, 0.008719393631837252, 0.009676407229202679, 0.010335420225646269, 0.010374591279870584, 0.010322642808118154, 0.011688261389681903, 0.01226643422024623, 0.01254998858768061, 0.012172518636603926, 0.013657886642867099, 0.013843129114929999, 0.013165742417160842, 0.01416152298271302, 0.013577062655542899, 0.01318213207943946, 0.01236699587420319])
      cub_degree = 34
      tol = 5e-14
    elseif q <= 36
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = false,
                                      numS21 = 9,
                                      numedge = 9,
                                      numS111 = 36,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.01894176483035853, 0.061816073159419206, 0.12643600829156484, 0.20850063414635683, 0.30321120631787235, 0.39991144153174507, 0.49159217052238263, 0.5809793286938126, 0.6447965478902975, 0.9903718524469571, 0.9679672494063327, 0.933438989044975, 0.8876841304760279, 0.8318882011451556, 0.7674964320159432, 0.6961765918569547, 0.6197758529614933, 0.5402729686194109, 1.9183549033456793, 0.018571681116591108, 1.8504983815398175, 0.018308977800245974, 1.7604201504525863, 0.017942212872515937, 1.65040574936511, 0.017427672026863218, 1.5233685705884894, 0.016588533230919168, 1.3821944273495466, 0.015932908007251292, 1.2307658587225636, 0.015168906392015537, 1.0734132411859731, 0.013578807143490055, 1.81057839628979, 0.060899702444417045, 1.7232545175922793, 0.059663272186598844, 1.616578872729244, 0.05803280037805141, 1.4937559918624357, 0.05525690315823585, 1.3566074432759891, 0.052838950513637024, 1.2090095362007585, 0.05043856759487502, 1.0561123392325074, 0.045723136769648634, 1.6631059751082655, 0.12380169064276847, 1.5608943750162068, 0.12056324959425499, 1.444653013184646, 0.1149084436767297, 1.3148827077957512, 0.10925376682930046, 1.174047537156922, 0.10452589518367211, 1.0283930220217068, 0.09626773427103735, 1.485580142216857, 0.20313929586660695, 1.3766726349187943, 0.1938618040283219, 1.2569224171974605, 0.18323485028526507, 1.1259205826542678, 0.17533447001195032, 0.9903943984236928, 0.16416980204736373, 1.2922312602203085, 0.2895955678281579, 1.183090376406563, 0.27254445495898655, 1.0640427143587934, 0.2598230255879276, 0.9421013632971927, 0.24690690455392922, 1.102522816613522, 0.3761951569474201, 0.9968975134323241, 0.3563094643728103, 0.8881042879362973, 0.3402668446718541, 0.9239301286091602, 0.46542992626003427, 0.833587851232978, 0.4386134304285687, 0.7691290963988814, 0.5396788226945799])
      SymCubatures.setweights!(cub, T[2.770001564797166e-5, 0.0010062785005788472, 0.003046414801394318, 0.005838889420348297, 0.008849952199635313, 0.011552282381041936, 0.01312600347961231, 0.013277550238860935, 0.012390473018105084, 0.0066812450428201066, 0.00016681678492628578, 0.00029017288916928274, 0.0004034773141323907, 0.0005001595437516682, 0.0005745073685150242, 0.0006174654622842023, 0.00064691262211419, 0.0006490071080928749, 0.0005927803865428596, 0.0017521888546774128, 0.0024392620787933293, 0.0030288402997879333, 0.0034906215859451736, 0.0037583792392536803, 0.003926616345764845, 0.003955344730089208, 0.0036600726311437006, 0.004231686084821131, 0.00525088641261227, 0.006068789540069106, 0.006549202818373323, 0.006793860462687638, 0.006879471092431658, 0.006527552099696962, 0.007201742155965866, 0.008287834370056296, 0.008932994504270536, 0.009167613093281653, 0.009303716292930937, 0.009087932396937793, 0.01013115661387792, 0.01084944066606852, 0.010940437680141307, 0.01096086115422666, 0.010979595266702844, 0.01231649221266526, 0.012266730351759916, 0.01176760154159052, 0.011664636255788965, 0.01313073456410588, 0.012311206086884597, 0.011342216418399734, 0.012653838025683163, 0.010568382858040258, 0.00881570251925599])
      cub_degree = 36
      tol = 5e-14
    elseif q <= 38
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = true,
                                      numS21 = 18,
                                      numedge = 9,
                                      numS111 = 36,
                                      centroid = true)
      SymCubatures.setparams!(cub, T[0.017130461689212957, 0.056255950615688645, 0.11525556122196583, 0.19049292290426412, 0.27500927170437267, 0.3653216631347806, 0.4503942845194703, 0.5321603379829839, 0.6086339024072717, 0.9922198977344079, 0.9749681231366449, 0.9530078521863322, 0.9322140444198274, 0.9012549058521676, 0.8601544588818287, 0.8151184876160785, 0.7693430258227226, 0.7242809798928879, 0.991286148302274, 0.9709881484798728, 0.9396473776617953, 0.8980009630388561, 0.8470255130311116, 0.7879159801309153, 0.7220578916395011, 0.6509949282543824, 0.5763927579010928, 1.9260189539740704, 0.016896739142096075, 1.8644642012883639, 0.01664716288365133, 1.7824850475909926, 0.016319288532510072, 1.6820682824023627, 0.015732724086128032, 1.565528115609385, 0.014831963088829664, 1.4349615951543966, 0.014442479192257869, 1.2942265142051905, 0.013410364916528264, 1.1467622601757383, 0.01130838277257212, 1.8275301869038898, 0.05542781403266054, 1.7477965296845848, 0.05429789300619424, 1.6504309064587088, 0.052423969153442926, 1.5377802922450499, 0.04946146802257783, 1.4103479520687638, 0.04787249153033048, 1.2735832834487555, 0.044574685113981015, 1.1301650703594774, 0.03987713700916072, 1.6925471728489834, 0.11283137996688535, 1.5993905725086583, 0.108990480925924, 1.4926436496439215, 0.10316314830002635, 1.3716640739680208, 0.09875046761353254, 1.241060626805168, 0.09291914012525138, 1.1031663098324538, 0.08599390045395436, 1.5306779043076437, 0.18376066231745639, 1.429427928705733, 0.1748138133304228, 1.31680852752447, 0.16477160937559418, 1.192992036221837, 0.1576918914627673, 1.065209376259242, 0.14683841726772687, 1.3561890451272893, 0.26334303820584115, 1.2555865233160035, 0.24539602517020204, 1.1379877320912346, 0.23726684337456486, 1.0206205763182443, 0.21950500341316162, 1.1805578360991846, 0.3402977014644074, 1.075202109917414, 0.32746881419190976, 0.9715641532305442, 0.30356432332239647, 1.0074227794717605, 0.4242583994648633, 0.9116694974531468, 0.39767029102624096, 0.8542668216872868, 0.5004385346933509])
      SymCubatures.setweights!(cub, T[2.267301415850696e-5, 0.0006599580813047184, 0.0008235738055169466, 0.0025259054926754974, 0.004880419963247925, 0.007485951635056811, 0.009617412283551008, 0.011317563085990448, 0.012480876917240607, 0.011309310744603252, 0.01052717936147283, 0.0039323917502098, 0.006284544056575955, 0.006134881071190988, 0.006743523150868109, 0.009860047177015097, 0.011088606891134804, 0.010373966543882259, 0.009572989996311196, 0.00932691952703692, 0.00013653982231953874, 0.00023932212098579198, 0.0003331984500004814, 0.00041461880301125277, 0.000474478525026942, 0.0005079070480752788, 0.0005431590659173278, 0.0005362888087780533, 0.00044738071881541695, 0.0014446868865527227, 0.0020158778229340864, 0.002512169133071432, 0.0028872480094793767, 0.0030957864079450325, 0.0032998695431059183, 0.0032668141982831185, 0.002968048435635268, 0.003517961041659069, 0.004374542398960253, 0.005037203322298377, 0.005423189938344538, 0.00570842892732402, 0.0057035398590174055, 0.005595529038486804, 0.006041252064128942, 0.006910850778325506, 0.007470490913490806, 0.007663029602061697, 0.007900708770930551, 0.00785136346520369, 0.008521283815278257, 0.009278464139737529, 0.009201484804863025, 0.009620385907361846, 0.009220768197779347, 0.010537428491295757, 0.010533082700003795, 0.010846869112144773, 0.010207554221491798, 0.011534235089219342, 0.010939956655926321, 0.010931615627955301, 0.011413758381172642, 0.011101895662857044, 0.010944988084680246, 0.010162387481254524])
      cub_degree = 38
      tol = 5e-14
    elseif q <= 40
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = false,
                                      numS21 = 10,
                                      numedge = 10,
                                      numS111 = 45,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.015520195916893246, 0.051361835022348444, 0.10586161437431742, 0.17545300117935672, 0.25619669042302295, 0.34284594418975295, 0.4293872002621678, 0.5124420709366497, 0.5900107185979331, 0.6481200613115444, 0.992076219228823, 0.9736021419996144, 0.9450311450954523, 0.9069744638059606, 0.8602436199806011, 0.8058347191421296, 0.7449074375949511, 0.6787603550694598, 0.6088032925796425, 0.5365272700054492, 1.9328122124888099, 0.015422948202503575, 1.87676324921528, 0.015270747393298109, 1.8019422285286835, 0.014990275739266935, 1.7098170420035164, 0.014588484060144732, 1.602255837288058, 0.014090639605248229, 1.4815519124855931, 0.013457665634837926, 1.3501924770536882, 0.012837971455718829, 1.2110827019324186, 0.012142854715495081, 1.0673633108293215, 0.011137468970963366, 1.8422260290816894, 0.05086522748851043, 1.7691905317471874, 0.04991204833554682, 1.679553136565619, 0.048595097268011346, 1.5750337031979533, 0.04695383891511884, 1.4578292464741107, 0.04480631014918908, 1.3298091807397527, 0.042749167202004806, 1.1939138360347499, 0.04046860208607936, 1.0533422877744056, 0.03740904514137683, 1.7171768822425573, 0.10380633717448347, 1.6302624641466859, 0.10107928806113173, 1.529519425842615, 0.09767473292024843, 1.4174823173481104, 0.09311246760488218, 1.2951201915526183, 0.08879797208386385, 1.1651674058219412, 0.0842628238297553, 1.0306068277786176, 0.07859062752831351, 1.5659628061625543, 0.17083588671670435, 1.469889065097346, 0.16507458674346523, 1.3639545366717825, 0.15726725255304355, 1.2486795194846063, 0.14963869803640817, 1.1261406821085884, 0.1425712911966415, 0.9996827145834297, 0.13399948780146256, 1.3973616702595404, 0.2474868638716685, 1.2985419024092064, 0.23589110420097514, 1.1918136001921487, 0.22362199224160173, 1.077658513621494, 0.2138944472861949, 0.9611926851714845, 0.2023080066841032, 1.223718742192, 0.32739888034830167, 1.1273516141147384, 0.3095991276917263, 1.022908057505157, 0.29617571170406565, 0.9167704298899507, 0.28157897118524744, 1.0546955259769366, 0.4065760213185983, 0.9626075774216891, 0.38721344016843695, 0.8683483205798233, 0.36893446854103373, 0.8932110316253533, 0.48541232469064305, 0.8147171532779959, 0.4603246894147091, 0.7573104809598555, 0.5543869544469994])
      SymCubatures.setweights!(cub, T[1.8733529912621262e-5, 0.000676555638551697, 0.002107571670047606, 0.004134549177571797, 0.006383872897881562, 0.008498526748266772, 0.010027908945458948, 0.010841576456269942, 0.011101729830317963, 0.009956121376098925, 0.005088036584553657, 0.00011247915181118212, 0.00019892682227794998, 0.0002788529732188051, 0.00034831964216604206, 0.0004042122058318268, 0.00044512877436219485, 0.0004686258328935599, 0.00047882017377938983, 0.00047324308153912244, 0.000441773562256894, 0.0011975838070135236, 0.0016844417316453032, 0.0021105675129485008, 0.0024580337018506888, 0.0027147496854323374, 0.0028612333284024306, 0.0029278657933091185, 0.002897513057483757, 0.002727443277679162, 0.0029543829455717886, 0.0036871844072078313, 0.0042888703864718694, 0.004737652578352103, 0.004990352568067134, 0.005114765042966848, 0.005081938457541204, 0.004859955362031306, 0.005141573426722237, 0.00595330607801518, 0.006539242918431697, 0.006845608868783824, 0.006977655426671063, 0.006970381174762545, 0.006763214925896961, 0.007377398627382866, 0.008077010005250247, 0.008414945995250592, 0.008455012335309077, 0.008467422227418659, 0.008274814231094654, 0.00926651275353102, 0.009654188756017663, 0.009560472827671269, 0.009463069481406305, 0.00924447366834803, 0.010442682270972937, 0.010327449027846826, 0.00997676340721409, 0.009672189358248803, 0.010788826770672999, 0.010033328161951571, 0.009530420375103106, 0.01016550779054263, 0.008445069902179421, 0.007142285249059427])
      cub_degree = 40
      tol = 5e-14
    elseif q <= 42
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = true,
                                      numS21 = 20,
                                      numedge = 10,
                                      numS111 = 45,
                                      centroid = true)
      SymCubatures.setparams!(cub, T[0.014233148926712185, 0.04701338775399012, 0.09660872767377096, 0.16127739638660843, 0.2347205649806202, 0.3154121440631567, 0.39687572043788605, 0.47629443397932136, 0.5426729083242992, 0.610553369219444, 0.9782055115268008, 0.9933093126526155, 0.9313496589712144, 0.9561125838596619, 0.9080510683102087, 0.8790481420356077, 0.8426952349487032, 0.8061034777255904, 0.7707527118783836, 0.7249585097471642, 0.9927635779393662, 0.9758789778553552, 0.9497292790201726, 0.9148255483256429, 0.8718475205860303, 0.821631822230068, 0.7651558855684221, 0.7035189689572374, 0.6379207744728965, 0.5696381020203342, 1.9383935422120617, 0.01411354137541273, 1.8869741155620188, 0.013892950830894674, 1.8180212720305167, 0.013722732302447265, 1.733002590928326, 0.013348440287731347, 1.6334548072825574, 0.012838728642781055, 1.521135037434328, 0.012466651053483763, 1.3988798181635114, 0.011471068726915166, 1.2683900788357, 0.011053468580701496, 1.133351415222857, 0.009789353457658351, 1.855786173558956, 0.04632672880520466, 1.7885711855993294, 0.04575989542416093, 1.7064873069897608, 0.04441762949436742, 1.610335019722311, 0.04296946924799931, 1.5018517945256211, 0.04132684149182165, 1.3830585199448528, 0.038849394002806965, 1.2556027071161464, 0.03645698825209156, 1.1213124263168772, 0.034386249512717064, 1.7408586981297864, 0.09532205117777173, 1.6609313916139088, 0.09226413495832157, 1.5664431184139582, 0.08968732566841942, 1.4613615677537795, 0.08533197455630548, 1.3460986220982665, 0.08159582624325751, 1.2248938358405959, 0.0756776214142127, 1.0958374549223058, 0.07380953552826391, 1.6018130878414616, 0.1559896937375575, 1.5120163684900338, 0.15199440846089113, 1.4143801428844034, 0.14408203499012362, 1.3048795669278968, 0.1387076748173073, 1.1899544215741675, 0.12975405290336522, 1.0662247208323798, 0.12550018391789214, 1.4453855145420302, 0.22784980289664855, 1.352637407046122, 0.21667300745334464, 1.2496213909434069, 0.20731652557044283, 1.1410205993910338, 0.19770619799201775, 1.0285814095200863, 0.18645905570607216, 1.2843752699585265, 0.30217104860432603, 1.1918376406693225, 0.28695580521922914, 1.0885329253379485, 0.2762867683506422, 0.9860252231017388, 0.25805968175106847, 1.1240620660590723, 0.3773669775413808, 1.035786574850825, 0.35843544327897964, 0.9389109900731019, 0.3414161915534071, 0.9739975623234071, 0.4421894880946545, 0.884727697265363, 0.4277183366353174, 0.8399486237870114, 0.5134701565353947])
      SymCubatures.setweights!(cub, T[1.5619144975969795e-5, 0.0005133331562655253, 0.0005691905905626049, 0.0017667475053100875, 0.003470654011584678, 0.005416650529285685, 0.007293460230385215, 0.008480348427245613, 0.009363971349533947, 0.009876762911873392, 0.00953090686987899, 0.009096133555257336, 0.005223033024299842, 0.0030833316492791625, 0.006142913336803549, 0.006454432558952336, 0.006005389460082897, 0.007785741836859099, 0.008517441484072198, 0.007086709191050727, 0.007752594892733474, 0.009635371782259289, 9.4211095713958e-5, 0.00016643956634208033, 0.0002323186899573745, 0.00029271268619258066, 0.00034104240467025255, 0.00037449156832683353, 0.00040488451678151704, 0.0003963854512917873, 0.0004092977414984883, 0.00035544218791497745, 0.0010062525600777216, 0.0014097116006511936, 0.0017814769920649282, 0.002077299058190709, 0.002298359342157026, 0.002460474653175375, 0.0024663441642099916, 0.002457136677666414, 0.002328316888555884, 0.002468765864369683, 0.0031037459974505707, 0.0035922831328284825, 0.0040118553223337855, 0.004236919586696688, 0.00442082404565826, 0.004290099140617136, 0.004435900340243366, 0.004362830782311943, 0.005043848257685337, 0.005651047155591732, 0.005864047673208063, 0.006139753665498885, 0.00596892640609939, 0.006090773075605744, 0.0062326114858510615, 0.006885244988684392, 0.007202980639374879, 0.007465236753222948, 0.007638588325325451, 0.0072559139447266035, 0.007968152874777005, 0.008478371620378176, 0.008398115191202428, 0.008745197660573458, 0.007967174056272841, 0.009111188215937486, 0.009188709943191748, 0.0091824982352256, 0.008954571165541868, 0.009535152723616649, 0.008574536872563748, 0.009508210476947769, 0.008889107391972897, 0.008337395084674924, 0.008534061539701987, 0.010234362615932248])
      cub_degree = 42
      tol = 5e-14
    elseif q <= 44
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = false,
                                      numS21 = 11,
                                      numedge = 11,
                                      numS111 = 55,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.013095897394543955, 0.043017094342087635, 0.08907497059301889, 0.1490344265403236, 0.2195232051017919, 0.2958351936562362, 0.37517021807125356, 0.45345538925577883, 0.5290398793213522, 0.5986058599251294, 0.6481290680149059, 0.9933652767525805, 0.9778741104649431, 0.9538528375567532, 0.9217320350774361, 0.8820852412102467, 0.8356200526320643, 0.7831656789896477, 0.7256581866071613, 0.6641238066877555, 0.5996606266954163, 0.5334189968686143, 1.9434109278357383, 0.0128930425234367, 1.895995393501014, 0.012816651787807158, 1.832549921759642, 0.012676782788731843, 1.754255195686827, 0.012354613812116454, 1.662378722834007, 0.011936533742460722, 1.558445359174444, 0.011528007438317786, 1.4442670765959433, 0.011176469698951023, 1.3220108869332803, 0.01068096201552778, 1.1938638760113487, 0.010016662811489628, 1.061982210956519, 0.009452998559355243, 1.8674988903120613, 0.042731153609825656, 1.805364085094794, 0.04225646522326377, 1.7289269304075328, 0.041230176763217365, 1.6392626512024728, 0.03984428843122967, 1.5375978109990696, 0.038434275809777256, 1.4255366751954346, 0.037247021642905234, 1.3055921217920914, 0.03562087846972901, 1.179986319515297, 0.03333675245819139, 1.0501923832817048, 0.0316367218637269, 1.7610118289470487, 0.08801974335601451, 1.6867285413107613, 0.08598747282575896, 1.6001883590962602, 0.0831134816707293, 1.5022352502936511, 0.08008565669432069, 1.3940252661754777, 0.0775233894588859, 1.2783627501678592, 0.07420944833297323, 1.1574842640435608, 0.06942938741817942, 1.0315035188763377, 0.06617935597931238, 1.6300910228168468, 0.14578968327032946, 1.547212690406877, 0.14091762799595228, 1.4537398429792807, 0.13566669115435814, 1.3505486410792749, 0.13102471798315016, 1.240648670625957, 0.12546855807065588, 1.1262925397145391, 0.11783617276973499, 1.0062000649449052, 0.1123555083867259, 1.481958741104801, 0.21216322972418483, 1.3929838067980165, 0.204077555800164, 1.2951718411645479, 0.19642437815371983, 1.1919937946486683, 0.1877895363828565, 1.0852729686084612, 0.17792971472789618, 0.9742528965147907, 0.1690995785228743, 1.325394610372731, 0.2844370648664548, 1.2342145966561928, 0.2728215754016323, 1.1380536048400207, 0.2599290173023054, 1.0374549655815244, 0.24878676438659975, 0.9368200130436194, 0.23544459791073136, 1.1667120488094747, 0.3590236773605816, 1.079294301783339, 0.3412178661841566, 0.9856722480813209, 0.3283089924309305, 0.8947980594157001, 0.31065773731505464, 1.0148684806023258, 0.43159798403333594, 0.9317052600119387, 0.4141514762099363, 0.8503816975260998, 0.3930140855442751, 0.8669483052763327, 0.5034636907338625, 0.7997343992828871, 0.47759802058293266, 0.7467613785074341, 0.5629920026429184])
      SymCubatures.setweights!(cub, T[1.3127387389360765e-5, 0.0004818370837918223, 0.0014857469048168963, 0.0029466607772683493, 0.004653996952482066, 0.00641132680858248, 0.007859615378077393, 0.008834700126159833, 0.009101413656291295, 0.009400148903858579, 0.008277753492774086, 0.004435609789777861, 7.949056446085082e-5, 0.00013948295766722936, 0.00019707335887574369, 0.0002491537438768069, 0.0002911630717816652, 0.0003231129020692911, 0.00034712400468129, 0.0003643399361880555, 0.00036841100245558014, 0.00035890626395456294, 0.000343798454817005, 0.000846157997558198, 0.0011959671022030197, 0.0015140275011545218, 0.0017739996840068404, 0.0019718730848074916, 0.0021184089836480536, 0.002225821190301804, 0.002254616058103938, 0.002192749979403082, 0.0021142835344030407, 0.002097103664216032, 0.0026535946395624295, 0.0031172207885047792, 0.003469500378131763, 0.0037235647082263503, 0.003910526674248321, 0.003968751618401844, 0.0038523249492846442, 0.003749941194247354, 0.0037104268799080554, 0.004351148884993413, 0.004826732559949885, 0.005160222127219566, 0.0053940910550221155, 0.0054818246470354605, 0.005350047490902341, 0.005226390632049418, 0.005454649016585168, 0.006026244583241074, 0.006406875661843602, 0.006628104519493036, 0.006698633688162521, 0.0066332424426098855, 0.006434599564604805, 0.0070854782710881945, 0.007498964825518566, 0.007640463554177111, 0.007564193917665245, 0.007583651487490952, 0.0072481361250443885, 0.008330752700493915, 0.008456505801299004, 0.008265344291584533, 0.008223819233606064, 0.007742507707503516, 0.008958664141921365, 0.008771328590273627, 0.008469324690807549, 0.007981178727476753, 0.009039249282115193, 0.008331245794695373, 0.0079432637168756, 0.008402730641758612, 0.006696328222547606, 0.005553509233893231])
      cub_degree = 44
      tol = 5e-14
    elseif q <= 46
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = true,
                                      numS21 = 22,
                                      numedge = 11,
                                      numS111 = 55,
                                      centroid = true)
      SymCubatures.setparams!(cub, T[0.011904952370168838, 0.039852593950582965, 0.0829884103203517, 0.13886479515659572, 0.20345452443843495, 0.2750070556578638, 0.3510896549488877, 0.4264000857463327, 0.5017666950377873, 0.5672767303405978, 0.6149817097746308, 0.9416432188136856, 0.9940360895298631, 0.9804064679778927, 0.9603623059929579, 0.924721803432006, 0.8943053348875658, 0.8597755177796159, 0.8253448367344777, 0.7874122526260056, 0.7481822425226893, 0.7044105581737128, 0.9938949724657469, 0.9796320691262672, 0.9574913853673113, 0.9278382329176582, 0.8911598296203584, 0.8480585244075671, 0.7992420736399966, 0.7455120574094392, 0.6877507289296136, 0.6269065320844383, 0.5639785297415535, 1.9482810448153494, 0.0119494253376906, 1.904909938290526, 0.011928506309555763, 1.8467041326419955, 0.011780809384832329, 1.774513843553716, 0.01146120765862152, 1.6893569641819162, 0.011036495674415267, 1.5924800602242417, 0.010659224994593852, 1.4855633102009433, 0.010244198757427108, 1.370328093515278, 0.009991217171331682, 1.2488065160038515, 0.009804233967438245, 1.1245123051457566, 0.0068457937961872144, 1.8770930741357128, 0.03978160134251367, 1.8197023489792952, 0.03928698311392778, 1.7492380029747525, 0.03820707351904685, 1.6664749400690755, 0.036846646471199763, 1.5722470588182353, 0.03558276505382042, 1.4681020687528012, 0.03421270855921786, 1.3550792467373145, 0.033340821297086654, 1.2356798976036139, 0.03215441923522565, 1.1138899975460455, 0.02691399572182326, 1.777273045250314, 0.0819029570854946, 1.7080898930808894, 0.07958690630372657, 1.6269309131951906, 0.07683592691010777, 1.5347681772331696, 0.07418026411279749, 1.4335611842193727, 0.07131739974269342, 1.324016659534221, 0.06944633746933634, 1.2102363592116052, 0.06573346687848165, 1.092863914939362, 0.06093075880481106, 1.6558509350897084, 0.1349590442807187, 1.5779560289774455, 0.13044463335254972, 1.489424690261929, 0.12599041960636, 1.3922962193057102, 0.12104038008855451, 1.2868714768715421, 0.11744815989416263, 1.1782611989710725, 0.11061928977707565, 1.0639357976670099, 0.10679814197075947, 1.5193766149936263, 0.19687713478429916, 1.4355327252939407, 0.19030265744190009, 1.3442737635335391, 0.18265187032163094, 1.245272607119078, 0.1765248615631583, 1.1429233465937254, 0.167669651175569, 1.0354494703678374, 0.1602879905731187, 1.3707219784401807, 0.2657789976578297, 1.2842893923298435, 0.2546885012742395, 1.1896570695292359, 0.2453052289960689, 1.0921770078382658, 0.23446569729743225, 0.9953823736783091, 0.2193899177849963, 1.2190709805026294, 0.3363512918702627, 1.132809995790961, 0.324272475596237, 1.0445739490148254, 0.3099573914320439, 0.9564519798671332, 0.2903977073634726, 1.0664385642019094, 0.41078480936330186, 0.9861428636729425, 0.38938865410318946, 0.9040673039435533, 0.37122225855678065, 0.9276870685460854, 0.4762168172604875, 0.8555901928047401, 0.4557775812529678, 0.8035020389655912, 0.5322772030169447])
      SymCubatures.setweights!(cub, T[1.1101611688514401e-5, 0.00042023214035481364, 0.0003985971451932693, 0.0012740969740858436, 0.0025662945094323424, 0.004064704833765581, 0.005531833507147237, 0.0069003492633715365, 0.007865329401393744, 0.008465010521951706, 0.008028527547337922, 0.0076244999108218536, 0.0059854263894260195, 0.0030493484095230824, 0.00254799827059452, 0.004346749049190449, 0.005336167197023197, 0.005921268930492443, 0.0071049822098239864, 0.007118992558291747, 0.006410094364508517, 0.0068150097531801445, 0.005534150241936276, 0.006617122025084275, 6.648134094147061e-5, 0.00011910794247860838, 0.00016912074022665014, 0.00021378284713331156, 0.00025007210207056766, 0.0002772515048544471, 0.00029892160409153635, 0.0003122383086737853, 0.00032436153133969626, 0.00033488017357522295, 0.00019528426210252325, 0.0007149936083282077, 0.0010191961854064354, 0.0012947835637481362, 0.0015214742879965533, 0.0016955148237978567, 0.0018316509603062828, 0.0019181379138901702, 0.0019903909063751593, 0.0020226830089409824, 0.0016464755241484177, 0.0018074102024248553, 0.002282279035961409, 0.0026661200463418233, 0.0029701455836976727, 0.0032075913084645626, 0.0033669872223430135, 0.0035035996668574855, 0.0034600313095041607, 0.0034486477578626446, 0.003236634682930242, 0.003775381829180317, 0.00419775859221513, 0.004506698655770465, 0.004689209810875999, 0.00482070799571421, 0.004644534747263025, 0.004960264905737076, 0.004737360404172615, 0.005263477181257924, 0.005651333682872205, 0.0058614579228295684, 0.005956208107206291, 0.005880570399919043, 0.0059219744368001715, 0.006152413221196505, 0.006584645573024671, 0.0067749804544887755, 0.006821967062724842, 0.006921214386951353, 0.006394785767904748, 0.007404694220755323, 0.007646679196483499, 0.007704808712540679, 0.007589382489289738, 0.006814719768021105, 0.0080574573024807, 0.008091792886971585, 0.007751076916619296, 0.0077983150063185925, 0.008523157342995925, 0.007924150993348393, 0.007794698677943712, 0.007828629433817453, 0.006850315110201509, 0.006106357537806797, 0.004769755560064076])
      cub_degree = 46
      tol = 5e-14
    elseif q <= 48
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = false,
                                      numS21 = 12,
                                      numedge = 12,
                                      numS111 = 66,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.011119467889261261, 0.03665811945604501, 0.07642878524650094, 0.12828261003304436, 0.18905837280791843, 0.25709026807170055, 0.3290307417855613, 0.40336778436094683, 0.47476617247194147, 0.5429455345240304, 0.6058697975943952, 0.6493613068957415, 0.9943637061557378, 0.9811889373838587, 0.9607177734087786, 0.9332621619795618, 0.8992385915537187, 0.8591629081813326, 0.8136426497461584, 0.7633678710149392, 0.7091006935331234, 0.6516637564296264, 0.5919277476350274, 0.5307982058909598, 1.95189859440565, 0.01097658116099749, 1.9114643377750045, 0.010980334720813344, 1.8572875759091876, 0.010867196116201355, 1.7901894039420845, 0.010565154148322126, 1.7109537483342225, 0.01027562187303299, 1.6207730863640524, 0.00995841220110371, 1.5210415339105365, 0.009530516355868443, 1.4130929978939826, 0.009235710489974793, 1.2987190621588334, 0.008842486096954464, 1.1795962118454495, 0.00845438788457881, 1.0575573819556179, 0.008093447008205752, 1.8867602716890373, 0.03662684845395589, 1.8334531874414257, 0.03625893873871482, 1.7677789456130042, 0.035278569354970556, 1.69007336626751, 0.03431530358245048, 1.6015845360060845, 0.03324888789321421, 1.5038061513930057, 0.03183763151339249, 1.3975054190572203, 0.030908062108643445, 1.2851224428323913, 0.02942383141935864, 1.1676996044530366, 0.028332978733737465, 1.0475331285822191, 0.026898810397137103, 1.7946329087307231, 0.07562530871783346, 1.7305854441872717, 0.07366035284800768, 1.6548951575658746, 0.07165032559032392, 1.56891515455282, 0.0694043019979584, 1.4742982269266964, 0.06651184246404333, 1.3709942914217061, 0.06469568230464036, 1.26260927868362, 0.061284578396489604, 1.1483338527269216, 0.05936312729846246, 1.0318346686332043, 0.055931607483866136, 1.6814677158173759, 0.12516203117379518, 1.6085173439114335, 0.1217169137312298, 1.5258071465423024, 0.11789809075212407, 1.4350567978100672, 0.11311159385949959, 1.3353982539709535, 0.10999687714407783, 1.2317872696689134, 0.1041508096122119, 1.121793537203718, 0.10089964074233762, 1.010545904157881, 0.09479874397125267, 1.5519413594269993, 0.18380652955612933, 1.4726967317511719, 0.17796768932720183, 1.3861454661907777, 0.1710589848144764, 1.2915696894507496, 0.16576630677274123, 1.1936235157669703, 0.1578133303762729, 1.0896852811396518, 0.15200041280328527, 0.9842473206945903, 0.1432020843749508, 1.410817873071676, 0.24865490392538594, 1.3278104524233427, 0.2395592841341166, 1.2381884697992513, 0.230481153863, 1.1447090415465841, 0.22119029653591937, 1.0492001395769393, 0.2114724472778892, 0.9521707783601906, 0.20061579960891388, 1.2645173155179505, 0.31778584834049756, 1.1820172979885823, 0.3040216601970652, 1.092288103499682, 0.2928793826728261, 1.0030807268457225, 0.27856665612212417, 0.9141045379228541, 0.26625355889287794, 1.1194565875275264, 0.3858742544072065, 1.040014307708822, 0.3712882896867678, 0.9595801191170028, 0.3535663233589538, 0.8758566624979641, 0.33917730467247104, 0.9783476201829914, 0.45415445511721475, 0.9080091688882351, 0.4323149561470506, 0.8341239645562638, 0.41572396578236787, 0.8488921011871938, 0.5162185097244729, 0.7908126193025644, 0.48955249301396486, 0.69402085505317, 0.5649016344389982])
      SymCubatures.setweights!(cub, T[9.467629934036714e-6, 0.0003475658941123335, 0.0010820564195515088, 0.002178937491421468, 0.0034795334674424926, 0.004842003997642271, 0.006109601952193788, 0.00704782958020841, 0.007624162229197525, 0.007962820101563673, 0.007683465020482708, 0.007316496374129097, 0.0040683445595667755, 5.733974727686064e-5, 0.00010103352104201348, 0.00014404476532493772, 0.00018277943638650944, 0.00021407446803176788, 0.00024050312222235334, 0.00026084129760777213, 0.0002723456307638426, 0.0002818198430368332, 0.00028387155106901376, 0.0002788869905551724, 0.00027238923837056146, 0.0006128903467566493, 0.0008735074545395374, 0.0011106843685432792, 0.0013040552777922323, 0.001467636877280983, 0.0015944667796046878, 0.001668459148011683, 0.0017313567973100293, 0.0017361349326978411, 0.0017202439843011761, 0.0016660197049720104, 0.0015397615186090075, 0.0019581357692882756, 0.0023026346841536005, 0.002592243737543934, 0.0028151864636366567, 0.002950014277926946, 0.0030723899298758647, 0.003051910093442544, 0.0030589133950389134, 0.0029227404733863168, 0.0027585907097379837, 0.003241809662907659, 0.003637645501191212, 0.003938000996939949, 0.004121881631329918, 0.0042894998472234845, 0.004246100234082178, 0.004273529178514525, 0.004061697810084652, 0.0040919639913471915, 0.004573097825860406, 0.004938899040600418, 0.005175806805507545, 0.0053446356459945, 0.005335147963520614, 0.005285208251604212, 0.005057330145350648, 0.005418297166425593, 0.005818758802509795, 0.006085979362346366, 0.006150661940914018, 0.006247715579416631, 0.006076864538896629, 0.005918014755396861, 0.006575773788370595, 0.006932073648151604, 0.0068582971985326644, 0.006913796225602649, 0.006536308474902952, 0.0064681437171873035, 0.0074411745205399974, 0.007476386807519159, 0.007444564238995607, 0.006969480903176569, 0.006605338601409084, 0.0076524986638299795, 0.007362309123766692, 0.007151480838961931, 0.006896692424857007, 0.007507339553594687, 0.0068334350528061975, 0.006419728657278589, 0.007162428955640784, 0.005206399748014299, 0.004617086817369184])
      cub_degree = 48
      tol = 5e-14
    elseif q <= 50
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = true,
                                      numS21 = 24,
                                      numedge = 12,
                                      numS111 = 66,
                                      centroid = true)
      SymCubatures.setparams!(cub, T[0.010317128425468431, 0.03408426486550986, 0.07082151117757383, 0.11861052838130542, 0.1750016963994337, 0.23885225665041623, 0.3065998661697263, 0.375758190185502, 0.44392226193785383, 0.5084602975306254, 0.5706177093665489, 0.617597868177169, 0.9949406487714625, 0.9833246801982458, 0.9675906662382323, 0.9614631052379022, 0.9407492839317382, 0.9129976047509656, 0.8828494384761305, 0.8521361436156975, 0.8206956604453594, 0.784620033084231, 0.7470964097571796, 0.7063596178106568, 0.9947804818642754, 0.9825742012254095, 0.9635917293625579, 0.9381010431072612, 0.9064602434479061, 0.8691135749232299, 0.8265853318484048, 0.7794725304712806, 0.7284365378070412, 0.6741937909945144, 0.6175057415514591, 0.5591681669492605, 1.955335236325035, 0.010207122524116667, 1.917819827533068, 0.010154138723891258, 1.8674493098010911, 0.01002352486955699, 1.8049408411117627, 0.00976244706716222, 1.73104117871429, 0.009481397969342244, 1.646681642586269, 0.009286670764823477, 1.553179068907099, 0.008974932224594594, 1.4518305324641818, 0.008513942140714445, 1.3436948047567707, 0.008441978894758064, 1.2309539129083646, 0.0076944714004082305, 1.1154450227318788, 0.005698405143023764, 1.8948642063588015, 0.03390175822635092, 1.8453782974795938, 0.03345147646922794, 1.7841341431506381, 0.03260860515867212, 1.7116339031041863, 0.0316819586881003, 1.6285798025363603, 0.03097999152191126, 1.5366066161229144, 0.029962856306171518, 1.4371126572412596, 0.02845204442554342, 1.3302665394453743, 0.02806774586290369, 1.2199503382265187, 0.02549310278862716, 1.1061793316368003, 0.0225056916191099, 1.809796581219209, 0.06982443737856277, 1.7500155403492186, 0.06811923894037257, 1.6793686707896436, 0.06623708594687433, 1.5984495451774479, 0.06462654624626273, 1.5090440161638787, 0.0625279228507765, 1.4125980009286467, 0.05952633518980003, 1.3084249575278484, 0.058130945813481244, 1.201597241257891, 0.0532986714900251, 1.0895229345466915, 0.05100769198061558, 1.7050619073165363, 0.1158053625052347, 1.636774231074243, 0.1127281764858152, 1.5586786606209218, 0.1096911612693891, 1.4724439035689354, 0.1061424252738812, 1.3798392793871326, 0.10137008410154746, 1.2802613135266026, 0.0976618451757949, 1.1772406184106081, 0.0915917357013412, 1.0680290685322047, 0.08921896987921049, 1.584619043820077, 0.17057509148576416, 1.510183271900166, 0.1655506181329295, 1.4275570417017813, 0.16013625246913368, 1.338745037187639, 0.1532779527369395, 1.244083263271472, 0.14604550659885895, 1.1442176277784128, 0.1397574038615874, 1.0407829692455532, 0.13401351671190034, 1.452058042456292, 0.23139015366294188, 1.3741196876338295, 0.22365239394793002, 1.291012066418219, 0.21425788642165203, 1.203118282932698, 0.2039367714141411, 1.1088796630299793, 0.19606692171281462, 1.0127620342233394, 0.18500540352398506, 1.3140009416135918, 0.2958395173150682, 1.2357567195304788, 0.2830521851200057, 1.1516872614478613, 0.2707532690109715, 1.06443853290276, 0.2574649618682823, 0.9752755475806139, 0.244380464047942, 1.1779866333480444, 0.35992825828179786, 1.1009861086788812, 0.34617205307080967, 1.0235321725434754, 0.32673658892091983, 0.9384402882872149, 0.31285079774583713, 1.0407927614509036, 0.4250114141170961, 0.9691584437281432, 0.4025873340894361, 0.892587363661434, 0.38555185874211345, 0.9171248534435187, 0.4869386093797817, 0.8521537085621156, 0.46042215723508895, 0.8043128980897319, 0.5343944748937508])
      SymCubatures.setweights!(cub, T[8.116671431254852e-6, 0.0003293428299921317, 0.00029930386300423313, 0.0009352337692633848, 0.0018781667657277819, 0.002995488158239551, 0.0041621053811475295, 0.005310774244266686, 0.006221854812218463, 0.006715528398147508, 0.007302054506262997, 0.006712278833133041, 0.006424948500154451, 0.0059996558473221955, 0.001998814222570394, 0.003420306809491281, 0.0026868627863105, 0.0024497459913375564, 0.005586695725104741, 0.0062478819931096885, 0.006022267077863082, 0.005904951778023037, 0.0058627736390862275, 0.006311283994494153, 0.00562204474212676, 0.006065489909139349, 4.9269862960777685e-5, 8.707840729347656e-5, 0.00012352435477229576, 0.00015664639185127179, 0.00018411759498593425, 0.00020702832582241657, 0.00022774595990227573, 0.00024079374434481565, 0.0002450313009694429, 0.00025650356620319137, 0.00024337686191510435, 0.00014893716768506982, 0.0005292261978613284, 0.0007512734640535887, 0.0009535258647910984, 0.0011232118889133024, 0.0012645666950834428, 0.0013905151508894784, 0.0014740484952512883, 0.0015020857170040351, 0.0015679676470621564, 0.001478513034688286, 0.0012736293644521727, 0.0013268063579868754, 0.0016823025020816442, 0.0019851343151075446, 0.002238602310181626, 0.0024551582376269596, 0.002605055188453562, 0.0026623050698167494, 0.0027458077735128877, 0.00259403346749828, 0.0026676890839772405, 0.0023736256043617294, 0.0027967911859927484, 0.0031501707437048623, 0.0034324065053480336, 0.003637938744175916, 0.0037385731927140615, 0.0037705457194186936, 0.003699967101642545, 0.0038762966386522517, 0.003530763044440027, 0.003981690643411783, 0.004310901946380905, 0.0045538879844447076, 0.004686993109894903, 0.004615047502682983, 0.004748671543335188, 0.0047020230144745965, 0.0047023903339686824, 0.005079485919877876, 0.005362344781406777, 0.005521710359975425, 0.00543207214012247, 0.005570518565282292, 0.0051163882890995875, 0.005720409299363746, 0.0059778097484784526, 0.00607905024798127, 0.00613104730573785, 0.006011414009728255, 0.0057413712865332376, 0.006512710262512335, 0.006641671830294419, 0.006787643972577052, 0.006277636609836412, 0.006281854557588678, 0.006869349398357839, 0.006929207692987254, 0.0067327528827480195, 0.0067503118948556325, 0.00704838679798439, 0.007036911862058815, 0.006424306121648763, 0.006966536565821554, 0.006084910779034753, 0.005846487302541092, 0.005471107802224679])
      cub_degree = 50
      tol = 5e-14
    elseif q <= 52
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = false,
                                      numS21 = 13,
                                      numedge = 13,
                                      numS111 = 78,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.00949525804860511, 0.03178634682982848, 0.06621618283408648, 0.11136441304721016, 0.1650445566610981, 0.22552551993911316, 0.29069698705000924, 0.35876388439553775, 0.4265535971002155, 0.4917850583113431, 0.5531285870766821, 0.6101307143607642, 0.651921892158718, 0.9951527013092271, 0.9838121429285657, 0.9661625835607792, 0.9424355086056514, 0.912940485028169, 0.8780620970027848, 0.8382550644647866, 0.7940383449185878, 0.7459883769657897, 0.6947315687881814, 0.6409361333108012, 0.5853033776540022, 0.5285585608467565, 1.9587511893015488, 0.009517042409715917, 1.9241059367626043, 0.009494980425366041, 1.8774903140515626, 0.009394947026362644, 1.8194415473144836, 0.009209166550113172, 1.750618902771926, 0.008965965929039714, 1.6718266889655926, 0.008709160875547859, 1.5840499497076719, 0.008448183305746341, 1.4884064076089032, 0.008189846443630121, 1.3862123943298215, 0.00784016379105014, 1.2787634444887457, 0.007472643295741368, 1.1674847847440686, 0.007124427114682028, 1.0538502384317967, 0.006849487548362965, 1.9018520318567336, 0.031695902714484615, 1.85557181205915, 0.031360306194271104, 1.7982157435219468, 0.03074008535665102, 1.7304062569518959, 0.02993440255981651, 1.6528502652045653, 0.029089790051147812, 1.566557826320941, 0.02817308610254698, 1.4724593710745748, 0.027383668325285045, 1.3722292935330542, 0.026150924812230886, 1.266683618660731, 0.025005034210172864, 1.1572467405500355, 0.023854901187115234, 1.0452532924031508, 0.02287563345773851, 1.821915318041066, 0.06551062456694705, 1.7656419664899288, 0.06421728049208841, 1.6992995245823634, 0.06255663187213907, 1.6235083811451643, 0.06080979058825037, 1.539408508143152, 0.058803284807235345, 1.4474242504860677, 0.05727261141382469, 1.350030489068599, 0.05463175688460539, 1.2472340467425858, 0.052351576491801985, 1.1406872833335153, 0.050029846002648606, 1.0316366643150774, 0.047832094397880644, 1.7225257705114096, 0.10917985302568294, 1.6580866871681272, 0.10641413150604136, 1.5845512658099368, 0.10343760311448717, 1.503254749606294, 0.09995462148391383, 1.4141478257423614, 0.09735551083691955, 1.3205284540220221, 0.09303264264854984, 1.221612555578173, 0.08918649932083207, 1.118836125286758, 0.08539344485025151, 1.01351109530042, 0.08145499954310628, 1.6079754861770568, 0.16098461730888627, 1.5373192129636772, 0.15643681544776888, 1.459123293480321, 0.1512166340584492, 1.373227160101755, 0.14689493601806675, 1.2830768233792575, 0.14089193786432805, 1.1884440025692204, 0.13506204684884157, 1.0905556714454785, 0.1293463464920354, 0.9904944925617503, 0.12337940677705385, 1.4817785625665343, 0.21911329617829522, 1.4073871181356148, 0.21201774648442243, 1.3261900482991738, 0.20524142546349225, 1.2404527382563337, 0.1975187196802033, 1.15062309986416, 0.1897173484381411, 1.0581576847717467, 0.18131678660409864, 0.9634534252093442, 0.17308108245910392, 1.348251461879078, 0.2815326490083148, 1.2718508120245604, 0.27171620496110716, 1.1908110624151134, 0.26152933439459203, 1.106034805918239, 0.2520785842821584, 1.0204630118508908, 0.24082803789770862, 0.932673309522442, 0.2297578429243045, 1.2114782540116367, 0.34573687864379493, 1.135654723639584, 0.3320195506532412, 1.0551517363079974, 0.3199659617620703, 0.9751066540066046, 0.3068872256498823, 0.8969438852471131, 0.292495569687757, 1.078187016738345, 0.40914126419750646, 1.005419011860457, 0.3924140427198347, 0.929546912374982, 0.37815428003174967, 0.8583524920983869, 0.36097188457589957, 0.9524713292236696, 0.4699468279237398, 0.8865758365908832, 0.45145728005620306, 0.8205943040509457, 0.4333886134305382, 0.8345157215110051, 0.5269117169627492, 0.7809740550099639, 0.503168722226905, 0.7352298477824402, 0.57384286619585])
      SymCubatures.setweights!(cub, T[6.994962237307251e-6, 0.0002536461057183018, 0.0008139266646389697, 0.001642644805847096, 0.0026503149441682527, 0.003715909187886374, 0.004757496523993951, 0.005657147367038863, 0.006376743403110127, 0.006723000348764388, 0.006723904883733086, 0.006661089335775839, 0.0063556444571070715, 0.0033084938852923815, 4.211022891107712e-5, 7.545271523351331e-5, 0.00010744115375245904, 0.00013669080617844152, 0.00016199548789004898, 0.00018292973630068773, 0.00019986110687306236, 0.00021308914115599257, 0.00022199388721413959, 0.0002251538851174764, 0.00022316415361982422, 0.00021859611732755774, 0.0002132668788291394, 0.0004548989028419482, 0.0006489428512891986, 0.0008280883737981209, 0.0009846157647324357, 0.0011154345738636812, 0.0012224815937351306, 0.0013043072552116565, 0.0013657589065991932, 0.0013837180673480505, 0.0013779683322676063, 0.0013503149994099706, 0.0013142149417103845, 0.001157799374108548, 0.0014738158930768918, 0.0017482292897469131, 0.0019772607550943557, 0.0021648874925964457, 0.0023004730276085994, 0.0024161747563669245, 0.002438372852545071, 0.0024405754440570966, 0.002398040839661893, 0.002325244411767343, 0.002087294601346198, 0.002471452123014578, 0.0027911572681865093, 0.0030492442595375356, 0.003229041302266415, 0.003389908014173067, 0.0034231635334637193, 0.0034246190940541298, 0.0033742595961919896, 0.0032563187915006456, 0.0031368222110432563, 0.0035418283237459904, 0.0038575164368777276, 0.004075533524363364, 0.004241483436736261, 0.004307216167928729, 0.004297862740850322, 0.004244423048025661, 0.004101239106281145, 0.004202379982516455, 0.00457536570139127, 0.0048482027263555385, 0.004994801437034231, 0.005091035284594073, 0.005058634693860585, 0.0049383713810226635, 0.004788771201297008, 0.005176361003828357, 0.005481009379671213, 0.005596460214899668, 0.0056850766601914005, 0.005699385611136264, 0.005520247487240062, 0.0053327410559146264, 0.006008148966208984, 0.006136330881892544, 0.006118574888050726, 0.006084396927532283, 0.005923253988019839, 0.005690157528295197, 0.00655082087551223, 0.006527356669484292, 0.006309896069358967, 0.006093974906400766, 0.005691001034937031, 0.006710804081085576, 0.006374690786246178, 0.006088044504019631, 0.00558941874944241, 0.006401002671808006, 0.005671042583332715, 0.005427972383435156, 0.0059378294383445026, 0.004479931504981677, 0.003992614058907151])
      cub_degree = 52
      tol = 5e-14
    elseif q <= 54
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = true,
                                      numS21 = 26,
                                      numedge = 13,
                                      numS111 = 78,
                                      centroid = true)
      SymCubatures.setparams!(cub, T[0.008936634794849808, 0.029518629262917817, 0.06151247058556487, 0.10365547108361463, 0.1533081811004403, 0.2102802435823432, 0.271517318726108, 0.3351555925413538, 0.40054107901851815, 0.464095576146433, 0.5213080039781178, 0.5761995289985111, 0.6258282247890377, 0.9849945166693987, 0.9954525460307236, 0.9671488897980254, 0.9701337563607053, 0.9478139314618067, 0.9237874378705382, 0.8985147753371993, 0.8716204761243544, 0.8437141669311823, 0.8122117161875361, 0.7770012702058653, 0.7398629458973487, 0.7080858660373063, 0.9954864941342849, 0.9849229036439682, 0.9684713592604912, 0.9463328599880441, 0.9187763681408934, 0.8861364486032344, 0.848809330678184, 0.8072481261017165, 0.7619573371859845, 0.7134867358567472, 0.6624246914209555, 0.6093910291421305, 0.555029506697796, 1.9613052421465447, 0.008834302465163544, 1.9287346069700284, 0.00880524728834118, 1.8849052968707654, 0.008739809754215367, 1.8304332289868401, 0.008520333284186707, 1.765843841585007, 0.008267027271532798, 1.6917510280798354, 0.00816186960389974, 1.6091863931069668, 0.00800316274534257, 1.5191685303865896, 0.0077015391876706, 1.4227062402904915, 0.007322401542953523, 1.3210908017716512, 0.006665400232373193, 1.2149146608889894, 0.00699107297397932, 1.1075431053041787, 0.004695628194939424, 1.908832893437691, 0.029415863224207094, 1.86564424090354, 0.029181958385296188, 1.8122036336036893, 0.02847449287967893, 1.7487672851519525, 0.027652529809171333, 1.6756101331687039, 0.027247421340192708, 1.594163321574639, 0.02669617751040634, 1.505660256367769, 0.02572360791835193, 1.4110392585403888, 0.02443949273982199, 1.3113670411828129, 0.022743503866892965, 1.206621762221501, 0.022400656292026214, 1.0999541906415407, 0.01958616449599274, 1.8343452835230465, 0.06096250731195918, 1.7820470815943195, 0.059522482296807795, 1.7200086364139053, 0.05790480696953999, 1.6482061145693667, 0.05689910196270043, 1.5682460092520945, 0.05566985116525726, 1.4814742036074657, 0.05366661760698872, 1.3887652757373723, 0.05101800241887257, 1.2909505146020108, 0.048259322413239165, 1.189119023706006, 0.04554039466784551, 1.0841949019611499, 0.04486370649469845, 1.741973193682388, 0.10128126138280494, 1.6817685018172805, 0.09874436180196827, 1.612395470106324, 0.09668107141052038, 1.535405560029677, 0.09451932065777212, 1.4526553797452257, 0.09114399980295453, 1.3642130122107419, 0.08712285128855322, 1.27032680494972, 0.08290093734091118, 1.1719268796371698, 0.07813291976024679, 1.067434136671434, 0.07785018728086907, 1.635280012187895, 0.14977872687188193, 1.5684022247382958, 0.14607404430719337, 1.4931576942188254, 0.14257353946346876, 1.4123062427812751, 0.13724102203792557, 1.325651846885923, 0.13189696683412191, 1.2355733178134247, 0.12528611978128038, 1.1410216403594264, 0.12044177032115522, 1.04416477392495, 0.11464165270263361, 1.5160843371251855, 0.20450447790788676, 1.444857141533289, 0.19921712695352908, 1.368797387827877, 0.19190073517064, 1.286412571139805, 0.1850672321474733, 1.200434934973735, 0.17613928258332617, 1.108412277384031, 0.1708075223191229, 1.0189461184819282, 0.15819970090626045, 1.3902046968039274, 0.26359907606690813, 1.3180186550841886, 0.2545510852723101, 1.241268833719637, 0.24508851850796023, 1.160657914697629, 0.23533299706483798, 1.075647668043917, 0.22612259986735142, 0.9886235751700969, 0.21319476837619045, 1.261958153125164, 0.3246748851746998, 1.1903366625420155, 0.31142832712063606, 1.1125361968456458, 0.300598177549036, 1.035562460846933, 0.2855795051282515, 0.9540507951226508, 0.274660818885039, 1.1336097697211285, 0.38460985572878503, 1.0631514192601637, 0.36907245893315505, 0.9893342301002254, 0.35248068891456097, 0.9133995711735049, 0.3363014443352938, 1.0112353578310416, 0.4404450392192418, 0.9432308282697722, 0.4261272321815249, 0.8794546284622144, 0.402522003628898, 0.8960166833103033, 0.49844326274896045, 0.8348572764802761, 0.47376323772750595, 0.7965672213260155, 0.5467340779778493])
      SymCubatures.setweights!(cub, T[6.066578269434501e-6, 0.0002744648476459219, 0.0002246494434711907, 0.0007027104390879714, 0.0014224703054383095, 0.002299972844442826, 0.003229739201981281, 0.004195091189675818, 0.00498099927303127, 0.005499351156906855, 0.005969628831700068, 0.006347000706756289, 0.00637426766987938, 0.005591466541842817, 0.005502651752048112, 0.0028891786326648177, 0.0016720229237619816, 0.0020110735374265945, 0.002063345894611114, 0.004666065113061425, 0.00490196663699584, 0.004948822621050313, 0.005181237956230799, 0.004837744683420618, 0.005466863410652697, 0.005537743181172141, 0.0047611232485706935, 0.004498128763077223, 3.690687415611264e-5, 6.524373010278725e-5, 9.289700733128489e-5, 0.00011875463625942165, 0.0001401447326368385, 0.00015801345428330662, 0.00017606544465064908, 0.00019008993208530578, 0.00019732376517748152, 0.0001993207526245431, 0.00018709316836189717, 0.00021008870501991047, 9.971699697589793e-5, 0.000397298350439537, 0.0005658637212936313, 0.0007235744658508771, 0.0008555697928034794, 0.0009659898582037899, 0.0010753751412457256, 0.0011621172225193298, 0.0012099467900879195, 0.001222176438654234, 0.0011791144314003137, 0.0012329575904819429, 0.0010277738571349735, 0.0010007740472793153, 0.0012789281669984833, 0.0015153494418741647, 0.0017157534818915418, 0.0019032299259115412, 0.0020512799691201596, 0.002137170671423938, 0.002156247835427318, 0.0021524535997702746, 0.0020490957271269363, 0.0022219926118824487, 0.001811281840603601, 0.0021432823268266966, 0.002431419578753789, 0.0026809769103515106, 0.002885837468689379, 0.0030101690947625143, 0.003055131006257567, 0.0030743789332224217, 0.002883144132600974, 0.0031666822503002047, 0.0027206850586329696, 0.0030869139312968276, 0.003361349232301003, 0.003589920690826795, 0.0037200002532280104, 0.0038367707982892877, 0.0038492966841592335, 0.003891956293620111, 0.0037384400245665706, 0.003685396942603902, 0.004019319807667004, 0.004304929565031352, 0.004458351959358785, 0.004583709393152876, 0.004499989467412575, 0.0046558203801681484, 0.003933192431578766, 0.004547170383959634, 0.0047959786781129445, 0.004993790652030696, 0.005121526984956532, 0.005194052441919575, 0.005130376701240762, 0.004826168753195443, 0.0052317995182362504, 0.005465899624118684, 0.005430111700217583, 0.005618909417301179, 0.0053051058086077945, 0.0056698141944962295, 0.005825302433393776, 0.0058852516845512275, 0.00586330461184867, 0.005527485727360928, 0.005609146038416579, 0.006140935613801044, 0.005827922328929832, 0.006126062457874812, 0.005265514814333989, 0.005866572503572728, 0.005701227559093005, 0.005543200968189821, 0.005641931343243777, 0.00539701856133372, 0.005146159821072999, 0.005145777124416784])
      cub_degree = 54
      tol = 5e-14
    elseif q <= 56
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = false,
                                      numS21 = 14,
                                      numedge = 14,
                                      numS111 = 91,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.008273529655134502, 0.02767779152454998, 0.0577201590418012, 0.09736059314479527, 0.14497523704840842, 0.19929939221682266, 0.2579451070659714, 0.32009364614220975, 0.38346560430976595, 0.4471385169592702, 0.5056457227664508, 0.5624746752190101, 0.6148487972133834, 0.6529017487048017, 0.9957869714202501, 0.9859233015831346, 0.9705523904755285, 0.9498496090996384, 0.9240497435900991, 0.8934451786187736, 0.8583826993185426, 0.8192595879037792, 0.7765191300475265, 0.7306455950841203, 0.6821587502112245, 0.6316079718597869, 0.5795660213129252, 0.5266225552427434, 1.9640604174882317, 0.008282333232636897, 1.9338352153252059, 0.008263825584371907, 1.8931018541470348, 0.00819224498375019, 1.8422781705651523, 0.008057735060472358, 1.7818732622337567, 0.007874365463783116, 1.7124925386764047, 0.007676145319027634, 1.6348757380197438, 0.007472679823457755, 1.5498655308388984, 0.007268520360934068, 1.4584566340106357, 0.007006225343221441, 1.3616655649285798, 0.006714926734799525, 1.2605924863530158, 0.006442615298732943, 1.1564513553454945, 0.006145180119064605, 1.0504609527986593, 0.005837048685941709, 1.9144820401403582, 0.027603387971697587, 1.8739775881725758, 0.027360227393468644, 1.8236079995361711, 0.026914022480609555, 1.7638877085375975, 0.02630185407057656, 1.6953504500637409, 0.02565354937128413, 1.6187969382369762, 0.024934255253403564, 1.5349601082593556, 0.024296134922056137, 1.4450818983208378, 0.023390263451214414, 1.349916566548346, 0.02246355334470107, 1.250439643884566, 0.021575824741192427, 1.1478668268692844, 0.020537545812941226, 1.043152697709654, 0.0196337271686886, 1.8445469412713724, 0.057204356616617484, 1.7950100301576, 0.05627811592377364, 1.7364527213818324, 0.055016052837921964, 1.669290494628982, 0.053677324968326875, 1.5944114452708849, 0.05209660742652303, 1.5121319551559234, 0.05082578687729815, 1.4241838225487968, 0.048901984432891676, 1.3309463391569496, 0.0470589172362201, 1.2336031396130511, 0.045212527267538556, 1.1334972592095323, 0.043037641456296204, 1.031242791051601, 0.04132632941562825, 1.7568604837024187, 0.09578646104068093, 1.699718532189719, 0.0937054342277059, 1.6343120187467257, 0.0913996895385416, 1.561672721002621, 0.0886570032741212, 1.4818206100964064, 0.08648301391376775, 1.396992185542577, 0.08331140629242133, 1.3069670264321855, 0.0803294609263235, 1.2128854212755054, 0.07710799283514337, 1.1156091103521508, 0.07362328400733019, 1.0157859392779294, 0.07065380171352345, 1.6547608621043293, 0.14198873076850352, 1.591674229799606, 0.13838514755602746, 1.5214187561403758, 0.13428410336978988, 1.4439474679499995, 0.13071486564376328, 1.3615464641286634, 0.12618391744074087, 1.2743515000397414, 0.12183467074972287, 1.1841855692102967, 0.11669972614007931, 1.0910296628288703, 0.11194863717863557, 0.9958132827058498, 0.10713826676964172, 1.5409579234002653, 0.19405094304144505, 1.4734188847314096, 0.18846650245562013, 1.3994177048041103, 0.1829855858455153, 1.320473930630761, 0.17705497741080867, 1.2371201529751545, 0.17118762434644874, 1.151426078363718, 0.16378800074034525, 1.062449960786634, 0.15744430948180993, 0.9722466724036608, 0.15032809093165506, 1.420026721048694, 0.2507176011039852, 1.349899677948564, 0.24294166891204433, 1.2742958472236448, 0.23528157721909926, 1.1948107744307752, 0.22748329552825702, 1.1137228996155646, 0.21815157872113025, 1.030301348540546, 0.2092651265073536, 0.94546915748307, 0.19982088691610922, 1.2941938706406613, 0.30999231910882125, 1.223052261344584, 0.30019836855344323, 1.14848451348847, 0.28950968048752534, 1.0710317915436751, 0.2791547601554792, 0.9938210145676761, 0.26686762332704045, 0.9153549525661338, 0.25507781625256526, 1.1666042376155277, 0.3711104346646404, 1.0978936147652716, 0.35622515639880326, 1.0238155154589883, 0.34485421910313846, 0.9519274420741594, 0.3299203183271917, 0.8812076353315272, 0.3155034038910355, 1.0441378790390883, 0.4279565827808882, 0.9769533147350372, 0.41338251963938166, 0.9098912832880053, 0.39799977576702933, 0.8454328226345293, 0.3809620376185938, 0.9289233241463704, 0.4851890667540783, 0.8680390261269431, 0.46746730585114177, 0.8092440103228008, 0.4494755070369711, 0.8194611659934062, 0.5373466112954692, 0.7703566036248477, 0.5146069252745412, 0.7289932571201156, 0.5802585120329661])
      SymCubatures.setweights!(cub, T[5.282711326452953e-6, 0.00019262322819311594, 0.0006181102063689229, 0.0012525702929889633, 0.00203824752477019, 0.002891256360201973, 0.003776579397442453, 0.004536631617426465, 0.0051907117009995915, 0.005673104534981359, 0.005875778144053051, 0.00587018175306527, 0.00580962815502939, 0.005484866926128754, 0.0028491332897460932, 3.1892766096225695e-5, 5.712681047063228e-5, 8.148404933303753e-5, 0.00010410778211540227, 0.00012414934621423008, 0.0001412208223358388, 0.00015546527653842104, 0.00016715537536456264, 0.00017579960117691877, 0.0001805752274534454, 0.0001814858401056411, 0.00018032160488116523, 0.00017631581406909167, 0.00016876473116595643, 0.000345277632590447, 0.0004931103338501017, 0.0006313371762708257, 0.0007549171303819756, 0.0008608665187562369, 0.0009505487159105603, 0.0010226396013238923, 0.00108011586598519, 0.0011100288031825744, 0.0011193653252430176, 0.0011142205626353964, 0.0010866810659382223, 0.001047817959971678, 0.000881059087839201, 0.001126228662015531, 0.001344868082132332, 0.001531518754521616, 0.0016897635447431316, 0.0018101679485425213, 0.0019138365474920165, 0.0019609522318980934, 0.001983757490425913, 0.001978842125953309, 0.0019281532705190806, 0.001880287906076681, 0.0015982071522338626, 0.0019051890953809434, 0.0021680322799119255, 0.002388436000539703, 0.002553432783741912, 0.00270014345108524, 0.0027689651858738493, 0.0028054308840540073, 0.0027892198620765372, 0.00272299104374955, 0.0026528003808407538, 0.002428007552272216, 0.002764614844452662, 0.0030332951641081417, 0.003234638581709465, 0.0033934697978967488, 0.00348738774531342, 0.003542287017941019, 0.0035174717667371493, 0.003479848598215802, 0.00336945983073106, 0.0032997191058163, 0.0036153465133297517, 0.0038710814920687917, 0.004039116939043227, 0.004162640241687986, 0.004210100881611787, 0.0041324331060890295, 0.004091840071121457, 0.003937880522310097, 0.004132334858556508, 0.00441868840017971, 0.004573206828606683, 0.004702154361132697, 0.004751324158413318, 0.004685439373918698, 0.004579325139088671, 0.0043997108528432746, 0.004863947448978945, 0.005054965018608284, 0.005173722581522091, 0.005156829672773772, 0.005131938566688306, 0.004919402116235496, 0.004763003740228883, 0.005408273384979739, 0.005522710985231223, 0.005428184464613638, 0.005454927339298679, 0.005167008987827663, 0.004971913873747036, 0.005790929729363306, 0.005644767357393278, 0.005567506705943306, 0.0053297575718334956, 0.004936430272912872, 0.005762435069243016, 0.005392922536162968, 0.005267114226867358, 0.00490347153537737, 0.0055130266591474205, 0.004893971709156552, 0.004652502667494252, 0.005194007405431523, 0.0036265343755212485, 0.0033308751531865295])
      cub_degree = 56
      tol = 5e-14
    elseif q <= 58
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = true,
                                      numS21 = 28,
                                      numedge = 14,
                                      numS111 = 91,
                                      centroid = true)
      SymCubatures.setparams!(cub, T[0.007791744392251025, 0.02584954272709248, 0.05404252143304533, 0.09111515199153958, 0.1353022121852652, 0.18680734571029725, 0.24238258352840367, 0.30019622573201243, 0.36161205210001296, 0.4214235735744895, 0.4800419772440027, 0.5347447656927256, 0.5813443016544089, 0.6287720059825438, 0.9674364795066475, 0.9960228693886543, 0.9747045155336808, 0.9869819549108823, 0.9522056812439704, 0.9310184574068863, 0.9090144418110911, 0.8871252370993976, 0.862108454688492, 0.833189623752548, 0.8013791889183032, 0.7675188099202325, 0.7354544234499601, 0.7050774440416341, 0.9960584221732405, 0.9868274679078682, 0.9724345851040197, 0.9530334757206349, 0.9288299976487278, 0.9000807715962316, 0.8670905681545376, 0.8302091013057628, 0.789827328604001, 0.7463733095494163, 0.7003076691402809, 0.6521187156364363, 0.6023172646237622, 0.5514312243803411, 1.9662148816542544, 0.00773513197352199, 1.9377351855733163, 0.007725186757158534, 1.8993570782914295, 0.007661521467342815, 1.8515013647542484, 0.007502640190083773, 1.7945831929559184, 0.007331478090705328, 1.7291291408124811, 0.007227973148119096, 1.6559129668596992, 0.0070638957694912005, 1.575722951736622, 0.0067638248317294985, 1.4891659620692552, 0.006619531956210375, 1.3973790175327758, 0.006261894428445822, 1.3012282203143755, 0.005874564891316181, 1.2013937293116448, 0.006242504274090594, 1.1010111006424794, 0.003569466829089728, 1.9200816566378103, 0.025822266898747125, 1.8822126220508362, 0.025593905012106326, 1.8351385479766196, 0.025078195703521645, 1.779063579846171, 0.024520322961345504, 1.7143334045990861, 0.02414346418466404, 1.6419627337242468, 0.023596691389495307, 1.5629848823886279, 0.02259624524303214, 1.4774174956614252, 0.022145572679117043, 1.3872885855587802, 0.020920948754112975, 1.292617697771869, 0.019954905171347303, 1.1939866559086052, 0.01984090600796206, 1.0944419778480925, 0.01653359900900349, 1.854413482381937, 0.0535116154696635, 1.8080623016980113, 0.05245120429580059, 1.752928457425704, 0.05134092766032, 1.6893454963538934, 0.05047223779043293, 1.6184895663143648, 0.04928045618906106, 1.5413956776841666, 0.04728321501604546, 1.4571643164696855, 0.046285376660132695, 1.368692170998722, 0.04386452525995837, 1.2753134841137739, 0.04209368926514558, 1.1790480789141955, 0.039942648705353974, 1.0807979988741396, 0.03869102053181708, 1.7726391490336084, 0.08932638094791645, 1.7186868657746663, 0.08756657636594402, 1.6564921706884648, 0.08589993001967594, 1.5873874546235054, 0.08372718765148629, 1.5126234441303203, 0.08064829127244802, 1.4314690283881817, 0.07861925189679776, 1.346523438924314, 0.07526467849936337, 1.2573779280547424, 0.07171612834658964, 1.1639921778944908, 0.06822550335485852, 1.0659377036828874, 0.06758994314866829, 1.6772942978741188, 0.13291174570905687, 1.6176318143307153, 0.13008921218842337, 1.5511451020539277, 0.12657307062993556, 1.4782205266067454, 0.122478356735255, 1.398927839683802, 0.11863226971447613, 1.3147612415884933, 0.11455205288472997, 1.2285272208863725, 0.10814606760824161, 1.137533446693557, 0.10499010753672135, 1.0459429204949813, 0.09979064204760948, 1.5689206138430638, 0.18234281754567544, 1.5050850948780314, 0.17708722878104657, 1.4350012924172821, 0.17203923031912924, 1.3604286759745285, 0.1661563152865207, 1.2806487873261772, 0.1610019621829817, 1.1988452411243016, 0.1525227688849302, 1.11049933315382, 0.1485798568244077, 1.0240393550964784, 0.1379174055137364, 1.4552171029579708, 0.23517557965850372, 1.388031413810047, 0.22890878153084623, 1.3167324443904018, 0.2210677176132759, 1.2412012413223048, 0.21332615872300958, 1.1635393446682405, 0.2044213861000466, 1.081741659168664, 0.1971715288822305, 0.9989338211060592, 0.1855111654399098, 1.3363815135226265, 0.2924506828969786, 1.2691484104775255, 0.28278841151234907, 1.1978244883703741, 0.27175140264046777, 1.1215279190524, 0.2616402464471571, 1.0437805819201522, 0.2504317438769851, 0.9658446767670523, 0.23876362306718, 1.2144324761822896, 0.3495174061797348, 1.1489520589022757, 0.3367229978743541, 1.081357615362776, 0.3229076433295217, 1.008524379061975, 0.3112962917522816, 0.9351544482951188, 0.29544085153378297, 1.0951032042631543, 0.40677322250861203, 1.0324425584277181, 0.38760445264398047, 0.9635855110751916, 0.3753259108364754, 0.9002734860225334, 0.354681647361528, 0.9840420995314307, 0.45965363871247683, 0.9236009206469115, 0.4407937910961164, 0.8621775370432926, 0.41938796483108837, 0.880013710807643, 0.5077979633811022, 0.8242197383826382, 0.48718022143095696, 0.7877625210857278, 0.5574635690843643])
      SymCubatures.setweights!(cub, T[4.624619483966628e-6, 0.0002248501725093247, 0.00017086409597539314, 0.0005392006737992557, 0.0011013983656237205, 0.0017898553269139095, 0.002526535050736664, 0.0033530792511240115, 0.004014136109917459, 0.004567739168534267, 0.005062542019452174, 0.005312641831068645, 0.005182862307355278, 0.005243989155957571, 0.004703114333841408, 0.004444184608558975, 0.0017217614430284243, 0.0013643744578540804, 0.0021549932030650756, 0.0023068356973078405, 0.0037779303318208448, 0.004107697574275374, 0.00401022782170233, 0.0039432033160277016, 0.0045504113321050655, 0.004854114849054312, 0.004868636922171954, 0.004581473570571962, 0.004132454580345047, 0.004076024245390253, 2.8099242884639098e-5, 4.9940107429011534e-5, 7.133560452346626e-5, 9.129734759735466e-5, 0.00010851690074698334, 0.00012362074264844896, 0.0001380194733946463, 0.0001490630778818499, 0.00015492560122240203, 0.00016173876818205958, 0.00016131592880847185, 0.00015552731672209728, 0.00017629105153413777, 4.3415841488893754e-5, 0.0003037513550322621, 0.0004343735188157398, 0.0005562720514512761, 0.0006624079419915032, 0.000755505767477079, 0.0008431708417992076, 0.0009125494827762496, 0.000949274530282474, 0.0009944979884446892, 0.0009908814699283567, 0.000977009078774421, 0.0010299154326020741, 0.0008173602054397562, 0.0007708553313770876, 0.000985660370695239, 0.0011752063904035513, 0.0013435674152032816, 0.0014974072830590847, 0.0016193377464146675, 0.0016851976780211868, 0.0017642694033731044, 0.0017555684815890149, 0.001769563716189198, 0.0016814337378087615, 0.001822831525825254, 0.001404249037810091, 0.0016701699727620254, 0.0019071400414874976, 0.002112808464907076, 0.0022760228941168716, 0.0023862327563319136, 0.0024884268079534934, 0.0025172703097416185, 0.0025038197691295925, 0.0023503235611900225, 0.00260772853574966, 0.0021298675655035725, 0.0024415171951963108, 0.0026919363493185617, 0.0028790074741281715, 0.003022400828805239, 0.003089421930446522, 0.003186594417988153, 0.003098543722148068, 0.0031833125094889113, 0.003071513246797846, 0.002903680210801975, 0.0031920447264062873, 0.003422567948957008, 0.0036492638484536136, 0.003729319897728112, 0.0038483249919330155, 0.0037165489928178918, 0.00381847301122625, 0.003250753536611758, 0.0036668694370556, 0.0038962050649933495, 0.004109081853609266, 0.004187411363862816, 0.004249991847887383, 0.004331447876235667, 0.004245735504399163, 0.004002628556636109, 0.004286453076251913, 0.00453204010652774, 0.0046637778772446025, 0.004578281397406255, 0.004702169808588006, 0.0044769500677399755, 0.004669374802831543, 0.00482620057621344, 0.005001694468687796, 0.005058056840264928, 0.005034252289686144, 0.004858963466172247, 0.004706900112186112, 0.005137060691789382, 0.0052117697572990635, 0.004948223476116001, 0.005190666734392235, 0.004803208290722894, 0.005426912548441108, 0.0053015985605568815, 0.004962372499583978, 0.004770257487033898, 0.005332531618905396, 0.004875951954545258, 0.004900840705477848, 0.0047240970974496335, 0.004693494442357779, 0.004612416359019425, 0.004507837149601399])
      cub_degree = 58
      tol = 5e-14
    elseif q <= 60
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = false,
                                      numS21 = 15,
                                      numedge = 15,
                                      numS111 = 105,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.007251670779451793, 0.02434218305518954, 0.05084065424725807, 0.08585170705300875, 0.1281169072736125, 0.1765494979315885, 0.23009066350079405, 0.28707589128941857, 0.34583020768970724, 0.4058148749378949, 0.4641109877537756, 0.5188810888777332, 0.5703877900437582, 0.619162005345566, 0.6533038439026205, 0.9963044669863808, 0.9876473452413546, 0.9741424192086162, 0.9559249695318659, 0.9331762380063378, 0.9061223658887212, 0.8750322469683374, 0.8402148778077754, 0.8020162935742106, 0.7608161314407826, 0.7170238586009234, 0.6710747032694407, 0.6234253294251026, 0.5745492984068238, 0.5249323625232967, 1.968473172973586, 0.00727975098993256, 1.9419441617949182, 0.0072722784546967595, 1.9061597421301892, 0.007210228513060478, 1.8614303274091117, 0.0070944895778347165, 1.808130364050401, 0.006945674517194969, 1.7467136426477199, 0.0067890368707273245, 1.6777327182769735, 0.00663565887379255, 1.6018335278119626, 0.006483632741925444, 1.5197613759834225, 0.006314396906950535, 1.4323393303498226, 0.0061090674131063635, 1.34045194600753, 0.005862441715661695, 1.2450339370261194, 0.005590208399198204, 1.147050986925407, 0.005339551644309833, 1.0474896558325695, 0.005174477247179303, 1.924721702388307, 0.024299105224569736, 1.8889706598588119, 0.024088918732567074, 1.8444150082810362, 0.023703088393705295, 1.7914282098407261, 0.023205213045386235, 1.730446570793267, 0.02268938366214746, 1.6620783980185239, 0.02216483277476832, 1.5869971303539883, 0.021655390040784275, 1.5060232169782355, 0.02108553739562781, 1.4199836249048394, 0.020399905929270813, 1.3296720421632997, 0.019594372770598535, 1.2358721987094983, 0.018705799359251474, 1.1393537051049292, 0.017879566152551072, 1.0409657464923077, 0.017337811334381587, 1.8629865251560132, 0.05039837648125345, 1.8191366708948893, 0.04960384589556019, 1.7670937662615565, 0.04857367431086085, 1.7071737960226092, 0.04751236249232246, 1.6399904607046232, 0.04638697839040638, 1.56610257514842, 0.045303455495493944, 1.4863902802988522, 0.044096087893792144, 1.4017421830959231, 0.042666462379508094, 1.313032800809671, 0.0410083728343135, 1.2211153684063396, 0.03919511275516262, 1.12669897658165, 0.037527587360429604, 1.0304849683312072, 0.03637243465790513, 1.7854014135913758, 0.08452596531122049, 1.734604215671282, 0.08280919780829588, 1.676161875396856, 0.08102208225759176, 1.6107965640955573, 0.07907495187023639, 1.5390054226393328, 0.07717982913480016, 1.4616814298178085, 0.07511236836787255, 1.3795893788700122, 0.07274406996924941, 1.2934640928437329, 0.06996815917776868, 1.2039165761571773, 0.06694064505492647, 1.111468241870924, 0.06423973793839437, 1.0170391512961658, 0.062070819790631705, 1.6944748661432625, 0.12560405952500273, 1.637726134203402, 0.12289972902677776, 1.57422536950507, 0.11991577934378123, 1.5044350828577777, 0.1169110174902077, 1.4292737008194083, 0.11371667847045322, 1.349623441832499, 0.11024852718437265, 1.2665710666160435, 0.1060933314213191, 1.180473554305898, 0.10157673299379538, 1.0912011349818822, 0.09777655330862117, 0.9999510463139305, 0.0940881139712794, 1.5921395432747465, 0.1727451367479026, 1.5308111528249044, 0.16852109434613752, 1.463452006451261, 0.16410124376592689, 1.3908581922875587, 0.15950132698629038, 1.3138074690229438, 0.1547982665292155, 1.233904450128909, 0.14903943161475153, 1.1514870547396403, 0.14276303419835312, 1.0660581245730791, 0.13767081660997596, 0.9792973979297642, 0.13203814167909783, 1.4811917629014173, 0.22442481045493803, 1.416738352498075, 0.21831840626306476, 1.3470509219407147, 0.21203557639188345, 1.2728835179201732, 0.20586380456525624, 1.1963954276692141, 0.19837105628978483, 1.1179042772558452, 0.19027477047824856, 1.037044810583227, 0.18324193655446686, 0.9554083305447656, 0.17559276560677897, 1.365070028264125, 0.2791077239764376, 1.2991308458260256, 0.2709981122522994, 1.2287506497595375, 0.26296253234717165, 1.1560547129828116, 0.2536008608785604, 1.0813721789712514, 0.24402790214042613, 1.0057815182498073, 0.2340200696115761, 0.9289682179535144, 0.22446929309482527, 1.2464627366582475, 0.33587100878031306, 1.180629280075315, 0.325461795571433, 1.1124019604994493, 0.31372505440929876, 1.0410746253355545, 0.30325489061165445, 0.9710700167980668, 0.2899004207673365, 0.8998328735347181, 0.27826274059286504, 1.1272088081978067, 0.3925328303452188, 1.06398520805875, 0.37765786612529145, 0.9964255847284751, 0.36595117143316547, 0.9314967033702795, 0.3505063589403834, 0.8671566306193936, 0.3365705785534148, 1.0146055843195467, 0.44614193794545387, 0.9539087284993859, 0.4309167672250851, 0.8928461929833194, 0.4155179280181685, 0.833647570677911, 0.39926172516811864, 0.9077996586244426, 0.49796627140181693, 0.8525146664640595, 0.48049542254946953, 0.7991046554497919, 0.46365467621862533, 0.8076852798558873, 0.5468158676146623, 0.7642305515446226, 0.5232070911938174, 0.691753394307883, 0.5839735754237688])
      SymCubatures.setweights!(cub, T[4.0628106204281e-6, 0.00014800661548916252, 0.00047889955138299853, 0.0009741098126060645, 0.0015915229305339923, 0.002277940948962854, 0.0029958531798422276, 0.0036823952689203666, 0.0042723377117095335, 0.004717943755678808, 0.005074354473995535, 0.005118074830859075, 0.005119882931333268, 0.005007176683625881, 0.004775323405586535, 0.0025661822257579125, 2.4523120132730433e-5, 4.407774574977598e-5, 6.303450446002676e-5, 8.068097198223918e-5, 9.646876585316746e-5, 0.0001102388147639833, 0.0001221090843555717, 0.00013226228593411982, 0.00014052390606839892, 0.0001464973692796303, 0.00014964519727235883, 0.0001497145037326732, 0.0001472015941879087, 0.00014352168200230288, 0.00014045332444631806, 0.0002662290487262932, 0.0003809911506360988, 0.0004886011611573817, 0.0005856746859956674, 0.0006709761171579266, 0.0007454162405921528, 0.00080906036679887, 0.0008618682747639518, 0.000900601914111442, 0.0009218922183487332, 0.0009245739699469028, 0.0009105603988314195, 0.0008880079563809486, 0.0008696867069577512, 0.0006841539251046939, 0.000876507629222672, 0.0010494372732856934, 0.0012002766703237042, 0.0013313126142907251, 0.001439945370886491, 0.0015293299589951776, 0.0015934500061351952, 0.001628459251156322, 0.0016346929658800759, 0.0016143213718336005, 0.0015794195395413826, 0.0015494414584540847, 0.0012456959125761277, 0.0014900216181536132, 0.0017034451466088763, 0.0018892192228193532, 0.00204178355082251, 0.002166853661944911, 0.0022566944670181738, 0.0023062491411220405, 0.0023122733173954117, 0.002280872057345306, 0.0022339476169244425, 0.0021801264540173685, 0.001902936704016913, 0.0021745018059092427, 0.002406213339144581, 0.00259288810960378, 0.0027397136008973566, 0.0028478173428654915, 0.002918754766853116, 0.002935918861893629, 0.0029076228283322464, 0.00286839415216409, 0.0027787240617576456, 0.002608282504936532, 0.00288644186308862, 0.0031109343060708524, 0.0032781347699730513, 0.0033958040634508597, 0.0034742851785505404, 0.003482575626551562, 0.0034431706694692597, 0.0034134871035667157, 0.0032912259312913947, 0.0033187609469007315, 0.0035760070231343963, 0.0037585560586780114, 0.0038821747030544796, 0.003967137709763172, 0.003969136559997059, 0.003908694578190016, 0.003837493881642498, 0.0036959261785667764, 0.003970591523697284, 0.0041763201196738705, 0.00430990113276311, 0.0043800382609877975, 0.004371778580298081, 0.004310013080735603, 0.0041393478459741595, 0.00399856902040058, 0.004491262767069942, 0.004642964393571003, 0.004695772592784067, 0.004672500520439521, 0.004640818759364855, 0.004383267628223129, 0.004239369574094451, 0.004887701650671983, 0.00492393291782514, 0.004863632391142119, 0.004856208888767718, 0.004584235906633665, 0.004379025432136933, 0.005094369318702593, 0.004993407084021285, 0.004855899646231922, 0.0046735915285663916, 0.004324127203779936, 0.005034953961829448, 0.004671062451367132, 0.004585933985953013, 0.004314346834210719, 0.004774692504609587, 0.004201940516171697, 0.00393350865278463, 0.004548525023234443, 0.003051008776668955, 0.0028519407766949677])
      cub_degree = 60
      tol = 5e-13
    else
      error("polynomial degree must be <= 60 (presently)\n")
    end

  else
    # do not include vertices in the cubature
    if q <= 1 #6 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = false,
                                      numS21 = 0,
                                      numedge = 1,
                                      numS111 = 0,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.21132486540518702])
      SymCubatures.setweights!(cub, T[0.3333333333333332])
      cub_degree = 1
      tol = tol
      # # Taken from Chen-Shu, 2017
      # cub = SymCubatures.TriSymCub{T}(vertices=false, centroid=false, midedges=false, 
      #                                 numedge=1, numS21=0, numS111=0)
      # SymCubatures.setparams!(cub,T[0.5 - sqrt(3)/6])
      # SymCubatures.setweights!(cub, T[1/3])
      # cub_degree = 1
    elseif q <= 2 # 7nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = false,
                                      numS21 = 0,
                                      numedge = 1,
                                      numS111 = 0,
                                      centroid = true)
      SymCubatures.setparams!(cub, T[0.21132486540518702])
      SymCubatures.setweights!(cub, T[0.16666666666666638, 0.9999999999999993])
      cub_degree = 2
      tol = 1e-14
      # cub = SymCubatures.TriSymCub{T}(vertices=false, numedge=1,
      #                                 centroid=true)
      # SymCubatures.setweights!(cub, T[1/6; 1.0])
      # SymCubatures.setparams!(cub, T[0.5*(1 + 1/sqrt(3))])
      # cub_degree = 2
    elseif q <=3 #10 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = true,
                                      numS21 = 0,
                                      numedge = 1,
                                      numS111 = 0,
                                      centroid = true)
      SymCubatures.setparams!(cub, T[0.1127016653792583])
      SymCubatures.setweights!(cub, T[0.1999999999999999, 0.08333333333333327, 0.8999999999999997])
      cub_degree = 3
      tol = 1e-14
      # # Taken from Chen-Shu, 2017
      # cub = SymCubatures.TriSymCub{T}(vertices=false, centroid=true, midedges=true, 
      #                     numedge=1, numS21=0, numS111=0)
      # SymCubatures.setparams!(cub,T[0.5 - sqrt(15)/10])
      # SymCubatures.setweights!(cub, T[1/5; 1/12; 9/10])
      # cub_degree = 3
    elseif q <= 4 #12 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = true,
                                      numS21 = 1,
                                      numedge = 1,
                                      numS111 = 0,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.40936128314152403, 0.1127016653792583])
      SymCubatures.setweights!(cub, T[0.1168096321294691, 0.44903690872888175, 0.05041006290415779])
      cub_degree = 4
      tol = 1e-14
      # cub = SymCubatures.TriSymCub{T}(vertices=false, numedge=1,
      #                                 midedges=true, numS21=1, numS111=0,
      #                                 centroid=false)
      # SymCubatures.setparams!(cub, T[0.4093612831422611;
      #                                0.5*(1 + sqrt(3/5))])
      # SymCubatures.setweights!(cub, T[0.11680963211922607;
      #                                 0.449036908703613;
      #                                 0.050410062902983145])
      # cub_degree = 4
      # tol = 5e-15
    elseif q <=5 #18 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = false,
                                      numS21 = 0,
                                      numedge = 2,
                                      numS111 = 1,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.06943184420297371, 0.33000947820757187, 0.457538013666131, 0.37414775838255504])
      SymCubatures.setweights!(cub, T[0.03019802974513124, 0.08091308136597997, 0.22222222222222215])
      cub_degree = 5
      tol = 1e-14
      # # Taken from Chen-Shu, 2017
      # cub = SymCubatures.TriSymCub{T}(vertices=false, centroid=false, midedges=false, numedge=2,
      #                   numS21=0, numS111=1)
      # SymCubatures.setparams!(cub,[0.330009478207572; 
      #                             0.0694318442029737; 
      #                             0.3741477583825542; 
      #                             1.168314227951314])
      # SymCubatures.setweights!(cub, T[0.08091308136598; 
      #                                 0.0301980297451312; 
      #                                 0.2222222222222224])
      # cub_degree = 5
      # tol = 5e-15
    elseif q <= 6 #21 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = false,
                                      numS21 = 1,
                                      numedge = 2,
                                      numS111 = 1,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.8478748148896262, 0.06943184420297371, 0.33000947820757187, 0.31541668315764937, 0.21896598857216268])
      SymCubatures.setweights!(cub, T[0.30870082526638737, 0.019432110768375376, 0.053997899611202924, 0.10555291032056137])
      cub_degree = 6
      tol = 1e-14
      # cub = SymCubatures.TriSymCub{T}(vertices=false, numedge=2,
      #                                 midedges=false,
      #                                 numS21=1, numS111=1, centroid=false)
      # SymCubatures.setparams!(cub, T[0.8478748148895112; 
      #                                0.5*(1 + sqrt(3/7 - 2/7*sqrt(6/5)));
      #                                0.5*(1 + sqrt(3/7 + 2/7*sqrt(6/5)));
      #                                0.3154166831592224;
      #                                0.2189659885706149])
      # SymCubatures.setweights!(cub, T[0.30870082526604714;
      #                                 0.0539978996110132;
      #                                 0.01943211076834772;
      #                                 0.10555291032094156])
      # cub_degree = 6
      # tol = 3e-15
    elseif q <= 7 #22 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = true,
                                      numS21 = 2,
                                      numedge = 2,
                                      numS111 = 0,
                                      centroid = true)
      SymCubatures.setparams!(cub, T[0.8768479048816373, 0.27886746283090746, 0.04691007703066802, 0.23076534494715845])
      SymCubatures.setweights!(cub, T[0.03707416966789951, 0.24947346457954658, 0.21085865924168884, 0.013202630162003234, 0.04106091936085787, 0.18219982239542837])
      cub_degree = 7
      tol = 1e-14
      # # Taken from Chen-Shu, 2017
      # cub = SymCubatures.TriSymCub{T}(vertices=false, centroid=true, midedges=true, numedge=2,
      # numS21=2, numS111=0)
      # SymCubatures.setparams!(cub, T[0.876847904881637;
      #       0.2788674628309072; 
      #       0.230765344947159; 
      #       0.046910077030668])
      # SymCubatures.setweights!(cub, T[0.03707416966789956; 
      #       0.2494734645795472; 
      #       0.2108586592416888; 
      #       0.041060919360858; 
      #       0.01320263016200324;
      #       0.1821998223954268])
      # cub_degree = 7
    elseif q <= 8 #28 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = true,
                                      numS21 = 2,
                                      numedge = 2,
                                      numS111 = 1,
                                      centroid = true)
      SymCubatures.setparams!(cub, T[0.18847440090927203, 0.4284337169945454, 0.04691007703066802, 0.23076534494715845, 0.16905856872047995, 0.6829175075115606])
      SymCubatures.setweights!(cub, T[0.021995779196974322, 0.106502996693885, 0.11909567136143454, 0.009155740220708433, 0.029322113201317515, 0.13315455373827445, 0.227422215281317])
      cub_degree = 8
      tol = 1e-14
      # cub = SymCubatures.TriSymCub{T}(vertices=false, numedge=2,
      #                                 midedges=true,
      #                                 numS21=0, numS111=2, centroid=true)
      # SymCubatures.setparams!(cub, T[0.5*(1 +(1/3)*sqrt(5 - 2*sqrt(10/7)));
      #                                0.5*(1 +(1/3)*sqrt(5 + 2*sqrt(10/7)));
      #                                0.22099843842186342;
      #                                1.1631912073287645;
      #                                0.24591837943530243;
      #                                0.12484120739275503])
      # SymCubatures.setweights!(cub, T[0.03946097492219484;
      #                                 0.022750835455651795;
      #                                 0.008818681441998117;
      #                                 0.18379350807584505;
      #                                 0.055866234615190226;
      #                                 0.2542415177018296])
      # cub_degree = 8
      # tol=1e-14
    elseif q <= 9 #34 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = false,
                                      numS21 = 1,
                                      numedge = 3,
                                      numS111 = 2,
                                      centroid = true)
      SymCubatures.setparams!(cub, T[0.4490679047936404, 0.033765242898423975, 0.16939530676686776, 0.38069040695840156, 0.13350211698690437, 0.2551152473178156, 0.6591253866697547, 0.18678645216668188])
      SymCubatures.setweights!(cub, T[0.0952257519287601, 0.006658503841344434, 0.016313739546142424, 0.025165807591009056, 0.06038167775800903, 0.14235832969148426, 0.20905439364578426])
      cub_degree = 9
      tol = 1e-14
      #36 nodes
      # cub = SymCubatures.TriSymCub{T}(vertices = false,
      #                                 midedges = false,
      #                                 numS21 = 2,
      #                                 numedge = 3,
      #                                 numS111 = 2,
      #                                 centroid = false)
      # SymCubatures.setparams!(cub, T[0.5267710694668858, 0.22581307560584146, 0.033765242898423975, 0.16939530676686776, 0.38069040695840156, 0.6556337684088949, 0.1729852501295953, 0.23384053251428352, 0.07010081712018432])
      # SymCubatures.setweights!(cub, T[0.16521286119089668, 0.09190991876959623, 0.006155310695930864, 0.015651771073189685, 0.02308837408072088, 0.13970715125683908, 0.020169336246406353])
      # cub_degree = 9
      # tol = 1e-14
    elseif q <= 10 #39 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = false,
                                      numS21 = 1,
                                      numedge = 3,
                                      numS111 = 3,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.07910043691664528, 0.033765242898423975, 0.16939530676686776, 0.38069040695840156, 0.7326534517622117, 0.16178998830038974, 0.15435874951213008, 0.3673737616984141, 0.47137911926303555, 0.5572798360390397])
      SymCubatures.setweights!(cub, T[0.03388444217225276, 0.0008514906758374474, 0.015608482275332266, 0.021266529714946317, 0.09830017774056278, 0.08553487986793913, 0.09482955197258906])
      cub_degree = 10
      tol = 1e-14
    elseif q <= 11 #42 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = true,
                                      numS21 = 5,
                                      numedge = 3,
                                      numS111 = 1,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.35406264500305373, 0.5349706832314904, 0.14620092105924026, 0.9344133285246975, 0.8608802427946676, 0.025446043828620757, 0.12923440720030277, 0.2970774243113014, 0.11288494321746362, 0.5189980362860592])
      SymCubatures.setweights!(cub, T[0.019568317008650336, 0.09027010441205281, 0.12668034231431224, 0.06200006585228387, 0.06976602848889582, 0.0831816128947378, 0.003678928929350306, 0.012519064458309727, 0.009408240262222357, 0.08199386419798452])
      cub_degree = 11
      tol = 1e-14
    elseif q <= 12 #49 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = true,
                                      numS21 = 3,
                                      numedge = 3,
                                      numS111 = 3,
                                      centroid = true)
      SymCubatures.setparams!(cub, T[0.944532850510785, 0.061319596747426944, 0.3670399210288128, 0.025446043828620757, 0.12923440720030277, 0.2970774243113014, 0.9643317137065957, 0.6791559176167298, 0.2800879551160052, 0.1108945564420906, 0.11759149845575509, 0.5898408487390895])
      SymCubatures.setweights!(cub, T[0.013033842870623723, 0.06517521728204109, 0.018835815821139735, 0.09104423407752954, 0.0007806066251635613, 0.008577761805886172, 0.012522944584712063, 0.08556960015635551, 0.04947754382659775, 0.0651662877443358, 0.10316420138769361])
      cub_degree = 12
      tol = 1e-14

    elseif q <= 13 #54 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = false,
                                      numS21 = 4,
                                      numedge = 4,
                                      numS111 = 3,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.08465074046995401, 0.9461615506412667, 0.674767791012365, 0.323267179957891, 0.019855071751231856, 0.10166676129318664, 0.2372337950418355, 0.4082826787521751, 0.5887945858454915, 0.10444119220651123, 0.6453371751137852, 0.34019090826916276, 0.09961977082750291, 0.29042446777707714])
      SymCubatures.setweights!(cub, T[0.021591599231498757, 0.06539932723237385, 0.041145418196861044, 0.07713054222104607, 0.0017474343859017836, 0.0059951229414363285, 0.008892253927912585, 0.010775615503604317, 0.05805264312502276, 0.10317694130545241, 0.0420598787031132])
      cub_degree = 13
      tol = 1e-14
    elseif q <= 14 #60 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = false,
                                      numS21 = 6,
                                      numedge = 4,
                                      numS111 = 3,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.3438598694407348, 0.08667709181254661, 0.23438830525749135, 0.5520607432679493, 0.9587852382179086, 0.8565327499763229, 0.019855071751231856, 0.10166676129318664, 0.2372337950418355, 0.4082826787521751, 0.09143726717923181, 0.31045214698586704, 0.24432313634373679, 0.5801444096258322, 0.07777721618947092, 0.6080084479074117])
      SymCubatures.setweights!(cub, T[0.054542527743094635, 0.023228939031591043, 0.02513762015408967, 0.10503584963092831, 0.05101963837324333, 0.07707180717717142, 0.0017455842845511217, 0.006113152456609967, 0.007165616322632858, 0.007968949263384184, 0.03911755531627651, 0.06148029583288167, 0.04172398880193796])
      cub_degree = 14
      tol = 1e-14
    elseif q <= 15 #69 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = true,
                                      numS21 = 4,
                                      numedge = 4,
                                      numS111 = 5,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.07246516124875677, 0.9597333330363812, 0.5578271063863942, 0.25839266543174033, 0.015919880246186957, 0.08198444633668206, 0.19331428364970477, 0.33787328829809554, 0.08945010025247659, 0.7665909769206356, 0.890366171837056, 0.8251504415876004, 0.27617698719951544, 0.5241882285519328, 0.24339513574480046, 0.08051715832948862, 0.4824242225832316, 0.08250669708663373])
      SymCubatures.setweights!(cub, T[0.007975211738078417, 0.01555795199453867, 0.0142287699046545, 0.09767682471952172, 0.04995668055544628, 0.0011984999191526156, 0.004057653525909466, 0.0057771570575940985, 0.007510437420733494, 0.041736293516403324, 0.043023900670943645, 0.07279841199275705, 0.027493765777003604, 0.03703949399671621])
      cub_degree = 15
      tol = 1e-14
    elseif q <= 16 #72 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = true,
                                      numS21 = 5,
                                      numedge = 4,
                                      numS111 = 5,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.04538479805133307, 0.22611045578916183, 0.7597184979885856, 0.4715660213460911, 0.962698662025747, 0.015919880246186957, 0.08198444633668206, 0.19331428364970477, 0.33787328829809554, 0.6642487752933809, 0.07317173680920826, 0.06769959788003582, 0.18359591467251551, 0.23363169305593603, 0.45309708381258046, 0.07154380405671455, 0.3955820133312661, 0.7299488157852585, 0.2402324353785272])
      SymCubatures.setweights!(cub, T[0.0068099595522426995, 0.007902251972551224, 0.03822679813612106, 0.08015956854713617, 0.07195779765864842, 0.03751715583277405, 0.0006453194613776787, 0.0033040868832377697, 0.0051468228469984285, 0.006303392104459153, 0.035057823925930626, 0.019961382833761403, 0.051626832885843146, 0.028924233456013058, 0.06107667308597509])
      cub_degree = 16
      tol = 1e-14

    elseif q <= 17 #81 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = false,
                                      numS21 = 5,
                                      numedge = 5,
                                      numS111 = 6,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.9673783153130373, 0.2285847353426557, 0.05840202475924722, 0.6899436810572654, 0.42617465029864315, 0.013046735741414128, 0.06746831665550773, 0.16029521585048778, 0.2833023029353764, 0.4255628305091844, 0.45997260718861244, 0.21436618225556944, 0.21119339160494885, 0.7433245013095352, 0.41869075536293493, 0.6667211762963206, 0.41159337035640003, 0.06703091561356278, 0.2024719492716725, 0.07034605000725505, 0.673480779947899, 0.06531116254401612])
      SymCubatures.setweights!(cub, T[0.03235009501507383, 0.03894783481393893, 0.010304291352971073, 0.023353742881921363, 0.05363694764402083, 0.000780347117914888, 0.0028284869642989496, 0.004156670349616389, 0.004888093779533072, 0.00532868081300272, 0.047922726227699664, 0.05405338591315981, 0.056272829027598656, 0.026314129806583212, 0.020943399709124497, 0.030548127770838585])
      cub_degree = 17
      tol = 1e-14
    elseif q <= 18 #93 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = false,
                                      numS21 = 5,
                                      numedge = 5,
                                      numS111 = 8,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.9733910151194378, 0.21626534843071651, 0.05150550298873318, 0.736133255549125, 0.37909135549530876, 0.013046735741414128, 0.06746831665550773, 0.16029521585048778, 0.2833023029353764, 0.4255628305091844, 0.3896301312622026, 0.1746624859673528, 0.7730831662450132, 0.08870182260636436, 0.21955949136089495, 0.978215349126951, 0.19128424042417297, 0.5706040919579012, 0.3950618272652804, 0.6260607152152675, 0.39506027304380237, 0.05823520227247069, 0.18823228992001764, 0.06686354946549783, 0.628350432052669, 0.04462834824245546])
      SymCubatures.setweights!(cub, T[0.021617628605668425, 0.030859684307638195, 0.008492717071468023, 0.048302368344678256, 0.05248145434612612, 0.0006756820364337081, 0.0026372964003534435, 0.003926222398471446, 0.003254151514351532, 0.00491785483632988, 0.025721171124103152, 0.020968675195521817, 0.035338100044051406, 0.037915704809379495, 0.057440515681043525, 0.021551248401012964, 0.019465345373257455, 0.0186444391812341])
      cub_degree = 18
      tol = 5e-13
    elseif q <= 19 #96 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = true,
                                      numS21 = 5,
                                      numedge = 5,
                                      numS111 = 8,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.5856059294178624, 0.047740969681407466, 0.36491553785928776, 0.818481134242735, 0.9099765582336337, 0.010885670926971514, 0.05646870011595234, 0.13492399721297532, 0.2404519353965941, 0.36522842202382755, 0.344515016089839, 0.05634768668558025, 0.3892623612407014, 0.1816493745951281, 0.0584137152546641, 0.16730180032945222, 0.5835945178396305, 0.3615699266691124, 0.19685560078378486, 0.1865954834733815, 0.05507891514214126, 0.8345201463254764, 0.17940998372119654, 0.635984563672555, 0.5715152111527353, 0.05498532291050753])
      SymCubatures.setweights!(cub, T[0.004144191713042269, 0.05671201687778216, 0.006952954178298502, 0.04182915174250011, 0.049236730600897276, 0.04362659548797906, 0.0005315273623156564, 0.0019591053495906494, 0.002952187933786907, 0.0035504848806637655, 0.003981570473101632, 0.018993243389823684, 0.03526459339555751, 0.014596159925672434, 0.047289681869773156, 0.014035668547812305, 0.025039659732181335, 0.041217429530185055, 0.022671200642619533])
      cub_degree = 19
      tol = 1e-14
    elseif q <= 20 #103 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = true,
                                      numS21 = 5,
                                      numedge = 5,
                                      numS111 = 9,
                                      centroid = true)
      SymCubatures.setparams!(cub, T[0.5436910855480633, 0.040950582871744044, 0.33887371496698127, 0.8109732389249305, 0.9624177194107271, 0.010885670926971514, 0.05646870011595234, 0.13492399721297532, 0.2404519353965941, 0.36522842202382755, 0.3334872800151244, 0.05388954156501192, 0.3895234197816185, 0.1718277785489534, 0.05577562117440405, 0.15389464625504432, 0.599253467648353, 0.14153780032196547, 0.1643849285186809, 0.21519353620251885, 0.5588940962149358, 0.3266182856250745, 0.04660262406085601, 0.801240166454181, 0.18939648231227538, 0.7890574015599582, 0.5525675668379594, 0.04329675777029787])
      SymCubatures.setweights!(cub, T[0.0043014134434171, 0.051486376804974164, 0.005540917165528219, 0.03694293372504983, 0.05781723615931413, 0.016819298793165075, 0.00044140660013591816, 0.0018162275341354346, 0.002944300406209069, 0.0028771277952727914, 0.003198212554961351, 0.0180917026903521, 0.029958760377674705, 0.013668758583780152, 0.02893920761789538, 0.01556903102368063, 0.04880488827336611, 0.01915230387979439, 0.035450671026311734, 0.017313930955951226, 0.05191629580852894])
      cub_degree = 20
      tol = 5e-14
    else
      error("polynomial degree must be <= 20 (presently)\n")
    end
  end
  mask = SymCubatures.getInternalParamMask(cub)
  append!(mask, (cub.numparams+1):(cub.numparams+cub.numweights))
  vtx = T[-1 -1; 1 -1; -1 1]
  Cubature.solvecubature!(cub, cub_degree, mask, tol=tol)
  return cub, vtx
end

function getTriCubatureOmegaLG(q::Int, T=Float64; tol=10*eps(typeof(real(one(T)))))
  @assert(mod(q,2)==0, "Quadrature degree must be even.")

  cub_degree = q
  mask = zeros(Int64, (0))
  p=convert(Int, q/2)
  ne = p+1
  numedge = convert(Int,(p+1-mod(p+1,2))/2)
  
  cub = SymCubatures.TriSymCub{T}(vertices=false,
                                  midedges=false,
                                  numS21=convert(Int,(1+mod(ne,2))*numedge),
                                  numedge=0,
                                  numS111=convert(Int,1/2*(numedge^2 - numedge)),
                                  centroid = convert(Bool,mod(ne,2)))
  if q <=2 # 3 nodes
    SymCubatures.setparams!(cub, T[0.33333333333333337])
    SymCubatures.setweights!(cub, T[0.6666666666666665])
    cub_degree = 2
    tol = 5e-14
  elseif q <=4 # 7 nodes
      SymCubatures.setparams!(cub, T[0.20233996166486962, 0.9396471009775607])
      SymCubatures.setweights!(cub, T[0.2514666037829547, 0.2664855267617265, 0.4461436083659291])
      cub_degree = 4
      tol = 5e-14
  elseif q <=6 # 12 nodes
    SymCubatures.setparams!(cub, T[0.12617802898300437, 0.4985734903418216, 1.2730049982427976, 0.10629009968963431])
    SymCubatures.setweights!(cub, T[0.10168981274041346, 0.2335725514527577, 0.16570215123674778])
    cub_degree = 6
    tol = 5e-14
  elseif q <=8 # 19 nodes
    SymCubatures.setparams!(cub, T[0.08637645651445221, 0.3898737393738146, 0.9750459783078007, 0.8632061170560124, 1.4849013997942493, 0.08090764459358865])
    SymCubatures.setweights!(cub, T[0.04814167840637482, 0.15337274330018216, 0.07488530992046893, 0.14871838633460427, 0.09350540993913128, 0.1636131864803212])
    cub_degree = 8
    tol = 5e-14
  elseif q <=10 # 27 nodes
    SymCubatures.setparams!(cub, T[0.06383921995140682, 0.29382965196340055, 0.5589751589490284, 1.6169889378735758, 0.05874192024953637, 1.2102133325028497, 0.050592079568966254, 1.0800239104248088, 0.2526484925513121])
    SymCubatures.setweights!(cub, T[0.026398989779595514, 0.09669613038232432, 0.10501730726845547, 0.051947647265061916, 0.05913975371921684, 0.10818971863386688])
    cub_degree = 10
    tol = 5e-14
  elseif q <=12 # 37 nodes
    # SymCubatures.setparams!(cub, T[0.04720134460354175, 0.23940087545921981, 0.4717396028948488, 0.9839248598439716, 0.9208772362624151, 0.8167789030246517, 1.708226268471533, 0.04785641668725236, 1.3872596590693262, 0.039467585843392945, 1.2547279674840568, 0.2011979341707161])
    # SymCubatures.setweights!(cub, T[0.01456411960169562, 0.0645983991965432, 0.08973969660331668, 0.034311756526720034, 0.0630737622175107, 0.08087271975980655, 0.03262024154086622, 0.03815621206243732, 0.07433302627360634, 0.08786175902176188])
   
    cub = SymCubatures.TriSymCub{T}(vertices=false,
                                  midedges=false,
                                  numS21=5,
                                  numedge=0,
                                  numS111=3,
                                  centroid=true)
    SymCubatures.setparams!(cub, T[0.04291787752125799, 0.257036207401034, 0.5107799738099695, 0.9754548840787666, 0.871924074925403, 1.7147294450904589, 0.05187389718509316, 1.394674530442853, 0.04184678032999987, 1.2337408348133263, 0.21104185636719278])
    SymCubatures.setweights!(cub, T[0.012473078997049074, 0.07047448976185372, 0.10838609608150095, 0.053535271338014945, 0.09312328568570066, 0.03494255432571696, 0.04099673260786038, 0.0754882765029034, 0.0774579537887578])
    cub_degree = 12
    tol = 5e-14
  elseif q <=14 # 48 nodes
    SymCubatures.setparams!(cub, T[0.03824798363517264, 0.18604036881658711, 0.39309856232165447, 0.589075871239251, 1.7667150057905754, 0.03632927991114931, 1.5038165491081694, 0.03301317443117578, 1.1676100984889524, 0.02943954702794549, 1.3954169939084178, 0.16887292886222252, 1.0942508605597852, 0.15041792985495467, 0.9722814163452018, 0.3486609061637871])
    SymCubatures.setweights!(cub, T[0.009533061380317499, 0.0412887986214927, 0.06518636222300395, 0.05691112344412955, 0.020056896904514798, 0.026015192788788614, 0.027030089794648958, 0.052892878141109795, 0.054926280081376344, 0.06595232278842289])
    cub_degree = 14
    tol = 5e-14
  elseif q <=16 # 61 nodes
    SymCubatures.setparams!(cub, T[0.030117843214746778, 0.1556519026128595, 0.33372894626399824, 0.5208327510323139, 0.9890255213475688, 0.9447249558670946, 0.8764734157647733, 0.7901961478694208, 1.8128573835969648, 0.030461705222347018, 1.5985309983228466, 0.027169587707121436, 1.3125710849128402, 0.025434956518319913, 1.4950653732615455, 0.1403170523845155, 1.2359624790668944, 0.1288560915430722, 1.1179100967131026, 0.29833558379257136])
    SymCubatures.setweights!(cub, T[0.005943216716219028, 0.02888666490266314, 0.05015588399565907, 0.054477578051096036, 0.018565524076026574, 0.036455789130213793, 0.042646342119967456, 0.048931161420656545, 0.013542772329513833, 0.017837779632168958, 0.020202917344582575, 0.03764374684567826, 0.04037637064367518, 0.051293710447503346, 0.05642973530376209])
    cub_degree = 16
    tol = 5e-14
  elseif q <=18 # 75 nodes
    SymCubatures.setparams!(cub, T[0.025392077582729995, 0.12671923598165075, 0.2824551127151226, 0.45648224953853817, 0.6068157196477899, 1.8437839991009135, 0.02449299060864029, 1.6625961864313492, 0.02300942825305127, 1.4204459548248227, 0.02126241445853266, 1.1388466147621379, 0.019051890198018086, 1.5794034488483435, 0.11898074296413326, 1.3550729275428062, 0.10983678254064705, 1.09259345380257, 0.09857546905948249, 1.2407139809876682, 0.26009443752060496, 1.0124012589321478, 0.2339918832282794, 0.9051975049525322, 0.4111820731149451])
    SymCubatures.setweights!(cub, T[0.004215716173423838, 0.019689247210426628, 0.03671250259566127, 0.0435961633359648, 0.03488587526584107, 0.009160369464051342, 0.01270834945600676, 0.014519084091668306, 0.01433753420583174, 0.02711834383018822, 0.030900586011567848, 0.030750146151876323, 0.04103413402702512, 0.0409529718843363, 0.04230206192012251])
    cub_degree = 18
    tol = 5e-14
  elseif q <=20 # 91 nodes
    SymCubatures.setparams!(cub, T[0.020963366959704022, 0.10795676521683514, 0.24446191111329224, 0.4026106291098426, 0.5488598343301395, 0.9926352500587178, 0.9619952856157791, 0.9108348444105016, 0.8464688567298789, 0.768520547539854, 1.8698414047026417, 0.020867645762484472, 1.7172709979730068, 0.019728258745342897, 1.5094938689177169, 0.018373496150538085, 1.2614829988201983, 0.01698644959278518, 1.6403118664047645, 0.10233054116516908, 1.4463473743742652, 0.09504597143168426, 1.2141791028254918, 0.08747594753486108, 1.3386197062718102, 0.22637839842626514, 1.1336247919862927, 0.20603105827470233, 1.0275096038401648, 0.3614368826752523])
    SymCubatures.setweights!(cub, T[0.0028807083331741866, 0.014321895363501704, 0.02815573726738288, 0.03649868822992064, 0.03613731776986962, 0.010306404493305683, 0.022057812947204034, 0.028850643940112424, 0.030507920649083954, 0.03378023183492541, 0.006516328548844917, 0.009254341880951164, 0.01088819790239361, 0.011403012815993614, 0.020113686591566204, 0.023352634166911843, 0.024204216097467007, 0.032168628659847, 0.032144271388721855, 0.03536925764896547, 0.03702046330458218])
    cub_degree = 20
    tol = 5e-14
  elseif q <=22 # 108 nodes
    SymCubatures.setparams!(cub, T[0.01794951619962759, 0.0917829334508462, 0.21065329229044205, 0.3533989792029848, 0.493294539143401, 0.6157127577585276, 1.8886819841434717, 0.017679795003491382, 1.7573871794640863, 0.016812592391618062, 1.5770295625872561, 0.01601409172483065, 1.3594857110449263, 0.01448717461409048, 1.1185520526457524, 0.012654004942552944, 1.6924522743511388, 0.08759264031454209, 1.5225060741926355, 0.08299107939024877, 1.317951272513367, 0.07550530552922186, 1.0888611903033274, 0.06635867196602277, 1.4237794487157651, 0.19790404483727267, 1.238887286812796, 0.18082124234920402, 1.032983177731578, 0.16098670513109353, 1.1355936848254982, 0.32312375037445623, 0.9568092778613884, 0.29173410363465735, 0.8666487186090581, 0.4488774154616421])
    SymCubatures.setweights!(cub, T[0.0021122781049839373, 0.010435853456364467, 0.02146775343852233, 0.029168477528704623, 0.030594524955223325, 0.025097419342539, 0.004739071418119735, 0.006806884729587624, 0.008269906198968913, 0.008639396776380516, 0.008062917558623582, 0.014925907511225106, 0.01786823770072817, 0.018951498736209646, 0.018123906234834047, 0.02517114291026744, 0.026378715984731883, 0.02571168624203832, 0.030284009690206507, 0.029911358343193233, 0.030050539885049736])
    cub_degree = 6
    tol = 5e-14
  elseif q <=24 # 127 nodes
    SymCubatures.setparams!(cub, T[0.015311224419790784, 0.07913277823095072, 0.18678441434204623, 0.31508302605222216, 0.4486750354591625, 0.5696770764824054, 0.9946959966276336, 0.9727783744644215, 0.9349379634957916, 0.882584429099999, 0.8208740954421498, 0.7522930732202445, 1.9049010917214353, 0.01516469223932067, 1.791734057295909, 0.015084917541981982, 1.6362867555719716, 0.013737190392591024, 1.444247587504066, 0.013056019514933697, 1.2267491952310887, 0.011561556966449489, 1.7308692345706498, 0.0781897205626327, 1.5829372014822098, 0.07185826322016423, 1.3989565011809344, 0.0675280252151803, 1.1926591915134788, 0.06088185427873798, 1.4916094963016353, 0.1739660869354208, 1.325384081990777, 0.16138797907293181, 1.1362622540139842, 0.1481032835889141, 1.2241805365962013, 0.28985704453680955, 1.0577834183247574, 0.26649135108352884, 0.9656723565870593, 0.40729335261167304])
    SymCubatures.setweights!(cub, T[0.0015373333986622828, 0.00785364121321484, 0.01663201023451273, 0.024372785574901353, 0.027587172959587007, 0.02574171216232596, 0.006353887790398897, 0.013583917838379566, 0.019556483810405742, 0.023029935353030558, 0.022673761467626653, 0.023387688177378938, 0.0034790248766018127, 0.005245802161283398, 0.006200483190627099, 0.006919888088994088, 0.006763086671263192, 0.011640082957440501, 0.013913221312176495, 0.01507802963746926, 0.01509331401876735, 0.02000496012244715, 0.021166035560835062, 0.021590553993757686, 0.025492399355196538, 0.024693634518574804, 0.025673880775368113, 0.025342626613911946])
    cub_degree = 24
    tol = 5e-14
  elseif q <=26 # 147 nodes
    SymCubatures.setparams!(cub, T[0.013346900919935074, 0.06939106610127004, 0.1615022661602191, 0.2788194105116561, 0.40373806676488133, 0.5219723356127394, 0.6231818947891458, 1.9169150795826768, 0.013302562795141856, 1.818104920313233, 0.012854300065964483, 1.6806923418646658, 0.012198489859669176, 1.5103729676559772, 0.011602842230753684, 1.3152794365554963, 0.01054460105809997, 1.1040389503693362, 0.009507136331565582, 1.7654493545477143, 0.06696300816749344, 1.6328499945438377, 0.06366203935828145, 1.4687519535425586, 0.06031148690686611, 1.2821484507711145, 0.055108602567997676, 1.0803275203641778, 0.04960675897649916, 1.5541327787381636, 0.15398943824326666, 1.402659039846475, 0.14509309348668187, 1.2295901694977003, 0.1336297725429019, 1.0412277251288813, 0.12010089398234516, 1.306595042191614, 0.2610801276104397, 1.1515527825636767, 0.2416530345418009, 0.984451320711026, 0.21834640830372123, 1.059547751226172, 0.3733483826043285, 0.9144570543574598, 0.34040145296969954, 0.8358821169204574, 0.4792737156258915])
    SymCubatures.setweights!(cub, T[0.001169057429760084, 0.006038179403032311, 0.012802726387537508, 0.01943962348369659, 0.022833165173914385, 0.02221608508472214, 0.018329367940057177, 0.0026684808070790724, 0.003930893623179237, 0.004860293995318182, 0.005485371520465311, 0.0055522869929046335, 0.005272852445069182, 0.008822651296459805, 0.010882957714261303, 0.012123968936954787, 0.012349486996597181, 0.011652702379806772, 0.015781491176125133, 0.017375596672941123, 0.017974734858227606, 0.017130978607230417, 0.021042500771456658, 0.02139383404077275, 0.020717448346548772, 0.02286714157032748, 0.02232601885795241, 0.021707539272295455])
    cub_degree = 26
    tol = 5e-14
  elseif q <=28 # 169 nodes
    SymCubatures.setparams!(cub, T[0.011716806964887383, 0.060677495913924594, 0.1440099178923534, 0.250847675082151, 0.3667386069108808, 0.4828178938148006, 0.5845696982964287, 0.9961179993474005, 0.9795186950956482, 0.9501160481685136, 0.9089933651213439, 0.857815411569713, 0.8028348395471911, 0.7442437821118563, 1.9271267926853528, 0.011599231440311563, 1.8399436933230913, 0.011473707135939646, 1.7185198470288054, 0.0108915045904305, 1.5668295761565023, 0.010330504364747627, 1.3910616735936816, 0.009392447040974113, 1.1979994365892728, 0.008621039211915585, 1.793145822731474, 0.05975823688065656, 1.675419443709008, 0.05700002822810855, 1.5289727291599016, 0.05372296053175259, 1.3600915725794935, 0.04946680618636102, 1.1745626723798273, 0.04511044450554109, 1.6012250805409602, 0.13821632206795104, 1.4653431656393416, 0.1295082416629256, 1.306927332333613, 0.12069384116239218, 1.1334160264111408, 0.1101354094058535, 1.3753702631306723, 0.23482231593257402, 1.2319453867783619, 0.21961229128487028, 1.0749462123232953, 0.2024976542562706, 1.1401793584033422, 0.3414509095543681, 1.0015028442756833, 0.31650111931722863, 0.9200660639223832, 0.4429869822425094])
    SymCubatures.setweights!(cub, T[0.0009009578771334429, 0.004646636566006364, 0.01019656127317631, 0.015895255812201484, 0.019894678115797638, 0.020871592824661854, 0.019161279405285045, 0.004067562880425828, 0.009224401242591496, 0.013536977238227658, 0.016795376531586175, 0.018036222013664255, 0.01610310966645167, 0.017418837613092167, 0.0020455138017436802, 0.003091041967327358, 0.0038501620400637837, 0.004372550949784034, 0.0044946004290830404, 0.004402835155270262, 0.00695221872976019, 0.008699527747900072, 0.009715289388782378, 0.010162255665511431, 0.009874289407605846, 0.012773366159742651, 0.014163090831944648, 0.01496894296141418, 0.01472865022763936, 0.017561239635569346, 0.018114235509647154, 0.01822829468628484, 0.020134072142523305, 0.01933487097811864, 0.01894165544964919, 0.02009942962690391])
    cub_degree = 28
    tol = 5e-14
  elseif q <=30 # 192 nodes
    SymCubatures.setparams!(cub, T[0.010339788456256407, 0.053898419453375326, 0.12810568967376337, 0.22413301045229955, 0.3324501401307675, 0.44194648095227895, 0.5426025316993733, 0.6288522782884479, 1.9355291534491872, 0.01030901530586667, 1.8582173915009468, 0.010119209205447439, 1.7497772551066502, 0.009762553565002285, 1.613862699209418, 0.009154230536998981, 1.4547829925673499, 0.008767542433637794, 1.27869791595414, 0.00806047878595933, 1.091788508027231, 0.007393933057319797, 1.8164485632693244, 0.05288655748515248, 1.7110570474293019, 0.050979190402608184, 1.5798142764330367, 0.0479482275704201, 1.4255487282967814, 0.045699338452326724, 1.2553118604284121, 0.042216854206935724, 1.0742319622865966, 0.038571839809832295, 1.6436822175466377, 0.1233374463975945, 1.519768575161266, 0.11649339189017176, 1.3742069485737125, 0.11027108017149774, 1.2138083173506566, 0.10244740387987994, 1.0432874489695694, 0.0933016528778228, 1.4377461303656685, 0.21257301346955987, 1.304642302914858, 0.19998813491623058, 1.1570190096167923, 0.18638147499127838, 1.0003636316691822, 0.16982915955767489, 1.2167489911453562, 0.3116596154757736, 1.0853658420377015, 0.29018308491514283, 0.9459966978938363, 0.26583027856667946, 1.0036914041872222, 0.4101166460361616, 0.8824498034474251, 0.37807739049643657, 0.8132618176972721, 0.5019246025849355])
    SymCubatures.setweights!(cub, T[0.000702220500263789, 0.0036635760586336756, 0.008195989137872497, 0.012856146510561892, 0.01661894070905062, 0.018091907507522004, 0.017036283625664544, 0.013976658412315934, 0.0016097745550546477, 0.002427163854655689, 0.0030831816644048596, 0.0034915175971017673, 0.003793272180286967, 0.0037827478240638843, 0.003603412214324012, 0.005492179937471774, 0.006933785607039577, 0.007873815853959569, 0.008460044073785832, 0.008507777821028835, 0.0080585599843871, 0.010282145762146986, 0.011690863166053713, 0.01235957287042385, 0.012492574750444383, 0.011828484657067699, 0.014637988957650972, 0.015327120486555539, 0.015420857586478948, 0.014753438195535987, 0.01729113101726095, 0.017027645070498528, 0.016505211714807577, 0.01762804079753721, 0.017019543107004806, 0.01638062079535912])
    cub_degree = 30
    tol = 5e-14
  elseif q <=32 # 217 nodes
    SymCubatures.setparams!(cub, T[0.009205719729844901, 0.04801375915397127, 0.11482522226359032, 0.20297525964568752, 0.3029379482038123, 0.40783482716763764, 0.5080079568899072, 0.5946147321064766, 0.9968900796419312, 0.9836541401042558, 0.9602627800653321, 0.9277153668106405, 0.8868767414224639, 0.8390103743315466, 0.7881143748709221, 0.7334406987845693, 1.9425875784217277, 0.00917130883550347, 1.8735600909388948, 0.009056987484021261, 1.7764598159122518, 0.00879184278954369, 1.6541392007738034, 0.008313990965931684, 1.5100488759123964, 0.007904907319391608, 1.349009853220651, 0.007256205291989251, 1.1759868181933166, 0.006759929348902141, 1.8360039110285078, 0.04735672315691457, 1.7410646192497414, 0.04597565098141886, 1.6223856337185119, 0.04352734338728842, 1.4824182467993607, 0.041324273981863203, 1.326463381428876, 0.03816463935214324, 1.1583295104530684, 0.035356194148873046, 1.6793354225140038, 0.11147274803959478, 1.5667120607220373, 0.10576161520065076, 1.4338221149728556, 0.10016923689339916, 1.286004540667027, 0.09323268844874721, 1.1267340247593274, 0.08603601476679355, 1.4894052419105017, 0.19313696307045156, 1.3664824172201555, 0.18253814059409076, 1.2294421718591961, 0.17087967301475362, 1.0820954432455299, 0.1580124165756617, 1.2827687707157112, 0.2861329322192584, 1.1588039409645212, 0.26811997539555427, 1.0254359231031414, 0.24938790772532793, 1.0755553508404063, 0.3812857090309576, 0.9586867382650691, 0.3549836397925539, 0.8870737997549597, 0.46911118619078135])
    SymCubatures.setweights!(cub, T[0.0005567205228851702, 0.0029164500118277355, 0.0066137955254072745, 0.010657921592627427, 0.014084873753065627, 0.016118800825126966, 0.016107137746153112, 0.014580783754360437, 0.002879743827561943, 0.006538711150467082, 0.009655426565589516, 0.012026408520597366, 0.013742249316793145, 0.01412952607270293, 0.012989151633967685, 0.013588957539608775, 0.0012767160365808698, 0.0019405543759390463, 0.002491158204285346, 0.0028612521239186972, 0.00311467216927763, 0.0031347865027278993, 0.0030709522758109857, 0.004407945816803521, 0.00563757336099153, 0.006468021143090092, 0.007002547765731732, 0.00712324541687532, 0.006920764982753192, 0.008417687511423287, 0.009659093905088621, 0.010374174861367193, 0.010643556587996535, 0.010329986162663595, 0.012231632234805515, 0.013043344367718205, 0.013355796817583695, 0.013122972140637738, 0.015001621367521647, 0.015079481874026317, 0.014866713683062907, 0.01598361749458508, 0.01505953378126753, 0.01464094265211411, 0.014877951223886063])
    cub_degree = 32
    tol = 5e-14
  elseif q <=34 # 243 nodes
    SymCubatures.setparams!(cub, T[0.008266884663220645, 0.04310539347820909, 0.10320070352110743, 0.18316244621832906, 0.2759477414810517, 0.37455992886544587, 0.4709371267081208, 0.5585175627225917, 0.6336401082579783, 1.9484254482320544, 0.008230015101628716, 1.8863226804336506, 0.008112260550386396, 1.7987260632509292, 0.007899433264885989, 1.6880120740940099, 0.007521999123551921, 1.5569794305334872, 0.007239642962484003, 1.4095276757503772, 0.006797977573397284, 1.2497136942971385, 0.0063754096297478, 1.082248054683094, 0.005887132541464266, 1.8527061161550051, 0.042473222935736396, 1.7669605266659114, 0.041334409598243575, 1.6591137799552735, 0.039425624254973604, 1.5311058424966202, 0.03785913617188616, 1.3874205092901606, 0.035599168339345184, 1.2315766628609721, 0.033347277606141146, 1.0683094893167695, 0.030759669470693193, 1.7111528725530942, 0.10035089186752755, 1.6081325285838692, 0.09594991114796082, 1.4858617618424308, 0.09184368863164219, 1.3489982945087098, 0.08651993918848004, 1.2003123801697135, 0.08094012454965716, 1.044104866646568, 0.07453969882726715, 1.5370647010807508, 0.17556152257723015, 1.4226797235466642, 0.16749694900232995, 1.294730962757294, 0.1579618291446958, 1.1559821070745186, 0.14763374660094614, 1.010013953862113, 0.1358670638559294, 1.3433962152098216, 0.2626016093572206, 1.2261981071292385, 0.24761730275720864, 1.0996700502142038, 0.23129429571396143, 0.9666508982978322, 0.21321886363083725, 1.1469105035995908, 0.35296961102085095, 1.0344380220617546, 0.3295915726628115, 0.9154805216708898, 0.3049588983607187, 0.9615552963469377, 0.4397795647130259, 0.8577902891058287, 0.4087382951178693, 0.7954746195992292, 0.5207403217458688])
    SymCubatures.setweights!(cub, T[0.0004490736894731839, 0.002353509414391198, 0.0053787667986239, 0.008785836333567505, 0.011897714520556721, 0.013812518333920144, 0.014270980942292032, 0.013290015986616084, 0.010784825125106253, 0.0010301575353532935, 0.0015660458771092451, 0.002021551116424048, 0.0023471940245202542, 0.0026002459746966335, 0.002695180814366277, 0.0026908890085105535, 0.00256246796622874, 0.003565040684912445, 0.004585095832694069, 0.005331502220736804, 0.005873087181601218, 0.0060970840890844915, 0.006068291050840314, 0.005764721454663133, 0.006884024629509667, 0.008002738532050092, 0.008734054898568485, 0.009075055798232373, 0.009023658053522583, 0.008573049607592685, 0.010222077902239274, 0.011071883306450722, 0.011456197115944143, 0.011355395668320846, 0.010825950504412407, 0.01282483002224189, 0.013150321061159462, 0.012940406599444496, 0.012413417633923843, 0.014096157865032958, 0.013812048517516821, 0.013344940407911135, 0.013922788428771158, 0.013516273633489348, 0.012777887742983824])
    cub_degree = 34
    tol = 5e-14
  elseif q <=36 # 271 nodes
    SymCubatures.setparams!(cub, T[0.0074329151494914964, 0.03894664900031569, 0.09332392115437489, 0.16678794374319456, 0.25386538708753836, 0.34580356805602774, 0.43951219325244845, 0.5271871484794897, 0.6052181329251465, 0.9974077140871066, 0.9865912650771974, 0.9678161657873201, 0.9415804776317647, 0.9081713102606398, 0.8678600591505573, 0.8221981241807722, 0.7760597907266883, 0.7285501609355293, 1.9535492241550318, 0.00743782907636874, 1.8975291487758088, 0.007300458481062214, 1.8181120517037932, 0.007216810248954028, 1.7175488221962911, 0.006837289939779979, 1.5979610915480882, 0.006626911491703685, 1.4627327082599755, 0.006211716461495279, 1.3150133776495185, 0.005873412467632962, 1.1585664785876384, 0.005499997086206304, 1.8669985176274586, 0.03829493898931899, 1.788977239127937, 0.03772361003122493, 1.6910986318544003, 0.0359409581881325, 1.574057834750127, 0.03461555402269072, 1.4417007349392037, 0.03262820106472629, 1.2969263313134105, 0.030694269267997634, 1.1438358472002628, 0.028778369483048814, 1.7375787219622623, 0.09144888031068091, 1.6430623229221404, 0.08778839549644009, 1.5310274671668727, 0.08388061354991899, 1.4045762192379398, 0.07957181230102374, 1.2667017089508745, 0.07460200297711231, 1.1197597365673417, 0.06978205519029841, 1.5769173039216158, 0.16123040118649054, 1.472206922336207, 0.15319746882300397, 1.3523488802089663, 0.1456625449736199, 1.2222733696653088, 0.13662119962176675, 1.0843661078196707, 0.12711126224297117, 1.395827407848706, 0.24098691273575235, 1.2856409622058869, 0.22875642056947904, 1.1660960617629215, 0.2152651442733251, 1.0397391650590795, 0.19996025485318075, 1.2088908340166638, 0.32731380432680135, 1.100699702591106, 0.30802013947117984, 0.9865427796588228, 0.28757318335252, 1.026199288176677, 0.41198687001161294, 0.9249896401159321, 0.38615425989029345, 0.8608941308475684, 0.49128891585743056])
    SymCubatures.setweights!(cub, T[0.0003632293546291326, 0.001920598286523193, 0.004438515016655851, 0.00729146593956933, 0.010179353567066972, 0.012238304432605396, 0.013264283946144078, 0.012979104138570171, 0.011697409676598292, 0.00214710426989189, 0.004738597433514434, 0.00703844347403328, 0.008940696862508962, 0.010477850514022958, 0.011593275382610952, 0.011418630053435507, 0.009931961680669066, 0.010457662264935622, 0.0008389305313777749, 0.0012744270339628684, 0.0016744727420152886, 0.001943661106209231, 0.0021758301480206666, 0.002268130426598389, 0.0023010553352944096, 0.002251457206718828, 0.0029118249492411613, 0.0037865433775298826, 0.004439554990415523, 0.004919140203123368, 0.00518440035318034, 0.005207638952568111, 0.005077013553518061, 0.005716357349435303, 0.006745134245889448, 0.007332703189925042, 0.007744135147722406, 0.007781902975338868, 0.007566896589458804, 0.008657664123990235, 0.009423351233248793, 0.009904164586398133, 0.009974515905086571, 0.009566937770684243, 0.011069466991700022, 0.011416548669194298, 0.011507673866546875, 0.01118509576753513, 0.012601060851594188, 0.012386663438057825, 0.012208031449467449, 0.012873480470148035, 0.0121611877005836, 0.011735945994848137, 0.011772545758280963])
    cub_degree = 36
    tol = 5e-14
  elseif q <=38 # 300 nodes
    SymCubatures.setparams!(cub, T[0.006731596857070331, 0.03531201365487037, 0.08496307922888283, 0.15188252628294518, 0.23181120888976062, 0.3193112819047896, 0.40840491367881343, 0.4935793715808355, 0.5710675103383652, 0.6371874283907188, 1.9579371561710661, 0.0067341783443620315, 1.9071434789967765, 0.00665960576527176, 1.8351787714579995, 0.00649544898617028, 1.7435097857347497, 0.006286244884718527, 1.6341183123418332, 0.006055336259857226, 1.50950905064432, 0.005763012635280278, 1.372516981671184, 0.005455080013071444, 1.2263439908426206, 0.005168128296851167, 1.074596286436765, 0.004738035727340244, 1.8790577923131033, 0.03489979858027476, 1.8081617660635154, 0.034040631080864386, 1.718078794420698, 0.03294545412979968, 1.6107388197456183, 0.031728502139476096, 1.4887581233908465, 0.030189970094257734, 1.3547006005723192, 0.028606636069766486, 1.2116749884336402, 0.027008989932633794, 1.0633340308798958, 0.02486771857403886, 1.7614652812173424, 0.08288850153453421, 1.6745042464218005, 0.08023476971376778, 1.5710311348383437, 0.07724608502526638, 1.4537609687597626, 0.07349516343453891, 1.3246750578599311, 0.06971469297100577, 1.1869218958780976, 0.06552596265120975, 1.04370909260954, 0.060631789121714234, 1.613599728154317, 0.14705516238987598, 1.5153636250117455, 0.14148771738620325, 1.4043823543477592, 0.13463401366250016, 1.2821240959622966, 0.12775718024816254, 1.152012024755516, 0.11968323527098738, 1.0160235369034314, 0.11112091186786444, 1.4447118060269004, 0.22282551528409303, 1.3412386442870696, 0.21209794018570036, 1.2273728499931094, 0.20110522994299102, 1.1068502332501473, 0.1882276618253417, 0.980563576844986, 0.17501474400756936, 1.2671779082284715, 0.3041352193258943, 1.1636793420324896, 0.28796430037826265, 1.053840412035123, 0.2698723780913897, 0.9386128519713935, 0.2510099086513202, 1.0914148227023015, 0.38613504813636956, 0.9934988592530242, 0.3626260290827232, 0.8906800865999976, 0.3378416923181421, 0.9278406064378751, 0.464125759565422, 0.8378557801847831, 0.4337889148301778, 0.7814467304564812, 0.5356289878925659])
    SymCubatures.setweights!(cub, T[0.00029789443395435676, 0.0015840139492238952, 0.0036716636125487872, 0.00615083894043689, 0.008595111233411919, 0.010473493203882949, 0.011549066886249224, 0.0115139665910034, 0.010548944695035058, 0.00861913896841509, 0.0006880663707766115, 0.0010532989571108748, 0.0013700555208241544, 0.001629195831669243, 0.0018264563573648456, 0.0019431815845965881, 0.001991048856168613, 0.0019828383041757304, 0.0018662563327425737, 0.0024160111141604116, 0.0031342830978667246, 0.0037177382389705707, 0.004156721185237865, 0.0044130160676743475, 0.004523511884312036, 0.004472254039613938, 0.0042425969657468214, 0.004752861895980843, 0.0056231224223146515, 0.006265851077925274, 0.00664217594011798, 0.00681266495830398, 0.006695480402441433, 0.006401524054792458, 0.007270941545387679, 0.008069994012335975, 0.008533812515834288, 0.00872080449046012, 0.008553271630841606, 0.00820625496229179, 0.009508631099704257, 0.010033402409041143, 0.010163181106056312, 0.009965101834893951, 0.00953624658884174, 0.01102826041802053, 0.01110652051347293, 0.010928693253593618, 0.010459182464347665, 0.011549793299009174, 0.011329593889370609, 0.010935271018037929, 0.011212145724188303, 0.01091096158385937, 0.010188991255774923])
    cub_degree = 38
    tol = 5e-14
  elseif q <=40 # 331 nodes
    SymCubatures.setparams!(cub, T[0.006135685729451402, 0.0321136480587115, 0.07745811098569379, 0.13970823552924658, 0.21384724580476475, 0.2958520277238817, 0.38089192294119606, 0.4646055871910022, 0.5430479697691739, 0.6107463845144343, 0.9979060004936974, 0.9890474223972371, 0.9734091826228977, 0.9513546413637778, 0.9233365499569893, 0.8899752462138527, 0.8516154271825221, 0.8089095605922344, 0.7663197453663518, 0.7224351021813991, 1.961663705655487, 0.006122932949158813, 1.915302837559188, 0.006052690819855725, 1.8494452149961076, 0.005987751723608712, 1.7655598189742396, 0.005748368128754208, 1.665136847970966, 0.0055505565477453, 1.55029939262898, 0.00531466282663926, 1.4233670758631343, 0.005074957347440149, 1.2869834363675954, 0.004829354426662089, 1.1441205256549771, 0.004519690663398513, 1.8899470331398969, 0.031759871276173454, 1.824901494386268, 0.03136457261289625, 1.7426232837057691, 0.030185666828891764, 1.6439092307553365, 0.02908881053842229, 1.5309752483633552, 0.02789054096428968, 1.4062600382846708, 0.026569741267999516, 1.2723703697215054, 0.025263803856722757, 1.1323842948760259, 0.023636427866187977, 1.7815352591320992, 0.07630226065053693, 1.701822684107487, 0.07368818532243279, 1.6065176780239712, 0.07085839938591848, 1.497466822903815, 0.06801344052462478, 1.3772553850609806, 0.06466912139777221, 1.2478737062682284, 0.0613123369333541, 1.11237334908983, 0.05739380734620685, 1.6444053005557195, 0.1353645066705022, 1.5538125525660758, 0.13003494781068886, 1.4499507393571303, 0.1247180335610021, 1.3357929204665453, 0.1185644943632337, 1.2132216535035003, 0.11192452867763082, 1.0843119809423432, 0.10483331689554794, 1.4871189284438062, 0.2055324629613967, 1.3896696058437472, 0.19665190638116592, 1.2823605105409879, 0.18704127512757582, 1.1679446256029582, 0.17608322945670077, 1.0477253300385754, 0.1647981713870123, 1.318707334231755, 0.2825159941691993, 1.2198240461675756, 0.26856633509129096, 1.1146295495487772, 0.25312161041553954, 1.0043168449516537, 0.23689353534610974, 1.1491454457164336, 0.3615350123977895, 1.0538428347647906, 0.3413310642313453, 0.9540017089004155, 0.32075971632424527, 0.9865805463572325, 0.43830302637276514, 0.8980941610277411, 0.41287486565273285, 0.8396672163508632, 0.5079639685895548])
    SymCubatures.setweights!(cub, T[0.00024753277763597467, 0.0013108708187609655, 0.0030689586665066815, 0.005202982105657326, 0.007395490366515022, 0.00925943395294559, 0.010442142069853161, 0.010881403296342743, 0.010534050067291596, 0.009426354976374892, 0.0015769067411551417, 0.003566782275800542, 0.005366727984172805, 0.006935016354046063, 0.008113434046584593, 0.008939404790850586, 0.009547279198362362, 0.00941832280602743, 0.008011032875408147, 0.008584213126873115, 0.0005707224435069417, 0.0008745653906498418, 0.001154997740914783, 0.001366117912669378, 0.0015391544967378612, 0.0016557468439244924, 0.001720579618351336, 0.0017363671060577955, 0.0016819708732891308, 0.0020057563928707433, 0.002635919585775433, 0.003130611144822596, 0.003513852199151729, 0.0037854237396669757, 0.003911327288851869, 0.003930251569261701, 0.0038031534379505505, 0.004008563770916177, 0.004770142012860521, 0.005319560924415143, 0.0057168481925545765, 0.005898135288676811, 0.00589065620750086, 0.005714854334003427, 0.006202167591560661, 0.006911715958672927, 0.0073716571297188374, 0.007598118811217097, 0.007528666328265829, 0.00731324290975237, 0.008256979458340242, 0.008738275033167242, 0.008972675896139495, 0.008872732849733802, 0.00854721208311355, 0.009782897010128163, 0.009934773577176636, 0.009873688442640888, 0.009569152782255695, 0.010580669388274157, 0.01039347053755265, 0.010175449593180372, 0.01057851996794738, 0.0099093319858995, 0.009344753225336438, 0.009766407655778103])
    cub_degree = 40
    tol = 5e-14
  elseif q <=42 # 363 nodes
    SymCubatures.setparams!(cub, T[0.005603025541436755, 0.029389618240164832, 0.07106134565331212, 0.12799066581760082, 0.19704393480980834, 0.2741069049530455, 0.3549916003201649, 0.43583420723331046, 0.512173279085049, 0.5804199662098773, 0.6395635804205347, 1.964974482889962, 0.005600191093930497, 1.922561390632371, 0.005555104251262358, 1.8622696010809545, 0.005446424580575802, 1.7851093396172117, 0.0053076348774683015, 1.6924825929672205, 0.005133359884367125, 1.586116661048414, 0.0049393911560631455, 1.4680994969541794, 0.004686713819492624, 1.3406804217314319, 0.004466183632542521, 1.2064620952909484, 0.004168470454249375, 1.0680289481679393, 0.003917134470843066, 1.8991457660559183, 0.02914072820939629, 1.8396411342666588, 0.028574554240826576, 1.7636792510387056, 0.027840453441638345, 1.672680768399397, 0.026934206553699287, 1.5682920883106803, 0.02589964658674142, 1.4525423819503764, 0.024623924813853955, 1.3273703902748555, 0.023408420063568434, 1.1954186671809723, 0.021935651729776223, 1.0590423739299306, 0.02055601157171618, 1.7998180108079236, 0.06969237121026628, 1.7258606199031834, 0.06788251179729023, 1.6374099711656935, 0.06568482947028417, 1.5361645176535803, 0.06310329631201007, 1.4240227478098035, 0.06014715332391158, 1.3028938973625848, 0.05701895825440189, 1.1750496873704874, 0.05367535476365192, 1.0428063130448144, 0.050177511751978196, 1.673059208294875, 0.12465535390359804, 1.5883132860184146, 0.12060911914126199, 1.491545921017872, 0.11580119009842516, 1.3842986044953307, 0.11060403414268394, 1.2687800751220821, 0.10468669443220191, 1.1464345123116026, 0.09880757717009243, 1.0198352292346329, 0.09229363203729346, 1.5258259897221649, 0.1905434162380074, 1.4343977341175984, 0.1829161652055589, 1.3330220866129376, 0.17484992577233346, 1.2243761577726635, 0.1655382500845974, 1.1093581886043022, 0.1562724034867959, 0.9902689662631681, 0.14610932910522298, 1.3671397603120312, 0.26322128159376096, 1.2731388852027603, 0.25161716597281414, 1.1722554054278147, 0.2385513165632887, 1.0656562827380345, 0.2249018046211617, 0.9549953774321521, 0.21048884062121745, 1.2044894213822535, 0.3392157372070796, 1.1126339635841989, 0.32204272558677977, 1.0159199731305575, 0.3034394673184409, 0.9148623072585282, 0.28410745490410505, 1.0463801333061937, 0.4138778035563614, 0.9603315309796133, 0.39017488624934477, 0.8701307110760793, 0.3656987187749486, 0.9008430843957381, 0.4832229076148667, 0.8217035688212547, 0.45402121052273486, 0.7709167681780463, 0.5469346461305029])
    SymCubatures.setweights!(cub, T[0.00020644515243425376, 0.0010994368636233001, 0.002583810129810963, 0.0044165470107815005, 0.006325475413626751, 0.007965442764429701, 0.009150480270921838, 0.0096583626178396, 0.009413292620532712, 0.008576049576692077, 0.00722492559515265, 0.0004770800425069569, 0.0007344791978758126, 0.0009641752742964974, 0.0011607162493833386, 0.0013158611819623725, 0.001427549311989774, 0.0014821996881457656, 0.0015044958410986202, 0.0014647915183074662, 0.001403412508566716, 0.0016875223386419193, 0.0022101426390681396, 0.0026531684252442934, 0.0030037838618271616, 0.0032516714920995697, 0.003388541859055468, 0.003426574845083184, 0.0033629875199789492, 0.003210192365838931, 0.0033792767718264175, 0.004047513850047306, 0.004572926931207644, 0.004933223042901251, 0.005150543358095162, 0.0051837726144143, 0.005115665421352807, 0.004876470929149357, 0.005283891199190715, 0.00595251824971902, 0.006406306836156863, 0.006685140499138853, 0.006724748437312275, 0.006632854723697899, 0.006341799630771279, 0.007110822296076363, 0.007644005030392193, 0.007935379983530331, 0.007979790674009259, 0.007812354730994895, 0.007501513417446152, 0.008554419820311592, 0.008856089807945635, 0.008923056532343751, 0.008690122448394693, 0.008332265568718954, 0.009420145098533344, 0.009447387011502592, 0.00921640242351258, 0.008850403148903339, 0.009597897408117766, 0.00934677345485714, 0.009027841151536917, 0.009120302845744443, 0.008883527487789063, 0.008356700328797694])
    cub_degree = 42
    tol = 5e-14
  elseif q <=44 # 379 nodes
    SymCubatures.setparams!(cub, T[0.005118934656555084, 0.026981577062377395, 0.06548283214890119, 0.11823727193391301, 0.1824075989280304, 0.2548974240245368, 0.33220493927346273, 0.4096337964941431, 0.4855847177857312, 0.5556637935037934, 0.6160322326569201, 0.9982579480603438, 0.990837025597508, 0.9776376056036035, 0.9590370835594006, 0.935359541908797, 0.9068658642353296, 0.874094585095366, 0.8374069587334187, 0.797631499389794, 0.7580696178660901, 0.7168363574020191, 1.9679695790759575, 0.005138954062278084, 1.9291554676274911, 0.0051144975432414995, 1.8739127485602627, 0.00502824943350591, 1.803087468166818, 0.004896559138601282, 1.7178060614874098, 0.00473859967020734, 1.6194789401693759, 0.0045648933466524285, 1.5097675211597466, 0.004390780020787801, 1.3906038204117717, 0.004194882985322652, 1.264097397723317, 0.003966558612496696, 1.1324787243972265, 0.0037729002296847136, 1.9072715375009492, 0.02683771280654856, 1.8524019119078299, 0.02637820887887859, 1.7822711105733218, 0.025688919764709946, 1.6979805860245951, 0.02486569093112259, 1.6009570158098365, 0.02393805030624331, 1.492784600842366, 0.02304304729721822, 1.3755869256970883, 0.021951418813144764, 1.2513377867792475, 0.02083543494298027, 1.1222819106661506, 0.019693988332604863, 1.815313650881646, 0.06434838140168932, 1.746917672672282, 0.06268716810593204, 1.6649055651060558, 0.06070773154653456, 1.5707497715597099, 0.058405616898994095, 1.4656919768122538, 0.05625842141462546, 1.3521968680744934, 0.05343858504329744, 1.2314300774564386, 0.050878513470915446, 1.1059719295250592, 0.04780241067958718, 1.6977157358979973, 0.11524245969249437, 1.6188599953644993, 0.11165549956071004, 1.5285023433588911, 0.1073918914285912, 1.4274871985887674, 0.10337093691235905, 1.3187743619225063, 0.09808895525800919, 1.2027791011766544, 0.09333179532093408, 1.0824685491842534, 0.08758503371268016, 1.5602715961189082, 0.17678518257797182, 1.4745241331729366, 0.17010535906498225, 1.3787963430255832, 0.16337006570631585, 1.27596005734609, 0.15524993527828734, 1.1665088235773877, 0.14720776986968065, 1.0525550220255169, 0.13841773830912824, 1.4099186835351278, 0.24552188050422485, 1.3207560050979945, 0.23512497568779112, 1.2243160049938016, 0.22392308316710163, 1.1225268064651757, 0.21173356544960473, 1.0162902234058409, 0.19934886181464198, 1.253917102420724, 0.3175143693273625, 1.164967538598104, 0.30269265303571885, 1.0716370271368996, 0.28627687499657783, 0.9741916912967002, 0.26976455555671913, 1.1003431710848728, 0.39047628527072664, 1.0159327060299672, 0.3696817381347457, 0.9275275353286734, 0.3491877659604576, 0.954076727889512, 0.45962729412115527, 0.8762460655055284, 0.4343783657879079, 0.8229906359123972, 0.5218739021044413])
    SymCubatures.setweights!(cub, T[0.00017233843722329505, 0.0009277140281514884, 0.002198478261671846, 0.003780142514881547, 0.005461280516183394, 0.007003589460657137, 0.008246718724097048, 0.008894700555231533, 0.0090601498451396, 0.00866639368500492, 0.00773176198064724, 0.0012056132304046016, 0.0027470888901324265, 0.004169816276136446, 0.005383293612775605, 0.006427517507335172, 0.007221790972925222, 0.007701779547647685, 0.008036223793204138, 0.007775808629456111, 0.006758231791593035, 0.007116878553078321, 0.00040028618975540656, 0.000618982608483182, 0.0008161400654015579, 0.0009843522189193275, 0.0011204820279659678, 0.0012226641861208332, 0.0012940273913425574, 0.0013259921179911616, 0.0013186499084018345, 0.0012893842316007364, 0.00143057215758393, 0.0018816625433858913, 0.002265254356695761, 0.002574479981015604, 0.002800245948076269, 0.0029612098499710257, 0.0030125064866929545, 0.0030078818397677684, 0.0029035313595224183, 0.0028840799562047208, 0.0034647660708310023, 0.00392968950299316, 0.004260619720346221, 0.004497658541421614, 0.004560125389477888, 0.004563938927615977, 0.004392466362603489, 0.004541155469294184, 0.005146442270710955, 0.00557473562136653, 0.005853216665794111, 0.005942612730933487, 0.005894208234095051, 0.005706284447367059, 0.0061870972125522864, 0.006705152337932016, 0.006974218536967979, 0.007100717416581014, 0.006975312231342241, 0.006789550689149731, 0.0076022060484322055, 0.007868916365371417, 0.008004236041145396, 0.007846706690170575, 0.007594109936942699, 0.008526189316011528, 0.008594880364905124, 0.0084570504960371, 0.008186406746955732, 0.008948499213288406, 0.008743565676644005, 0.008499500200647622, 0.00876356470391502, 0.008150527699685913, 0.007720889216894442, 0.007860452431334273])
    cub_degree = 44
    tol = 5e-14
  elseif q <=46 # 432 nodes
    SymCubatures.setparams!(cub, T[0.0047276569507481055, 0.024858144759086823, 0.06028855697064196, 0.10920938151502049, 0.1691358169079838, 0.2372081734008447, 0.310207184628501, 0.3851640582912604, 0.4585523704555884, 0.5272726294647907, 0.5887392032686983, 0.6420845549893645, 1.9704245757736985, 0.0047337173603027435, 1.9345495913132902, 0.004703836922943672, 1.883414703910566, 0.004631223364360326, 1.8177458063977256, 0.004527446155069051, 1.7385194147470124, 0.00440762680107076, 1.6469971224705593, 0.004252873439018215, 1.544636014393279, 0.0040865413372277275, 1.4331269187793867, 0.0038990339718496985, 1.3143060452507818, 0.0037001693929981363, 1.1901295319426808, 0.0035103330676026787, 1.062696470687294, 0.0032809135855909016, 1.9145826042700704, 0.02469235031893351, 1.8639157958012909, 0.02431237063384404, 1.7990395291424315, 0.02376559825969557, 1.7208945272998428, 0.02313482370046384, 1.630791301350361, 0.022323274090107502, 1.5300210155985614, 0.021460382170875425, 1.4202448182201934, 0.020475901858519576, 1.3031860267400872, 0.01944368776218315, 1.1807487024258547, 0.018429656610093857, 1.0550887749470481, 0.01725911391083946, 1.8297358258538605, 0.05936644288471424, 1.766293856309427, 0.058027245599735086, 1.6899649104382353, 0.05647274335296799, 1.6021898374534722, 0.05449783224946333, 1.5040897389094325, 0.05242216499304815, 1.3973877634409497, 0.05002770175851989, 1.283532675101696, 0.04754589167025436, 1.1642637172302364, 0.045019577018484636, 1.041593483873777, 0.042260147732181355, 1.7202491741153647, 0.10674969117630115, 1.6465076497065187, 0.10384104324756444, 1.5618509727553387, 0.10023511818747552, 1.4672309457122354, 0.09645633734000687, 1.3645345013867347, 0.09208699975598503, 1.254994284666779, 0.08757638883695038, 1.1403088920205977, 0.08286623973619373, 1.0221557302347755, 0.07793291445357844, 1.5914383664303293, 0.16442692042350482, 1.510721128386192, 0.15878069020226762, 1.4204818608904837, 0.1528207563947109, 1.32274283308765, 0.1459922368921729, 1.218583521521646, 0.13888736132182555, 1.1096122558622048, 0.13139596044480267, 0.9971965201789403, 0.12365992775938175, 1.4498506378060976, 0.22918341794125194, 1.365147480674066, 0.22059640882813314, 1.2733972323179814, 0.21092365756224227, 1.1756802417779202, 0.20066724098019864, 1.073311327618612, 0.18989217868835595, 0.9674676559775824, 0.1786184956570096, 1.3010466691183593, 0.29849827982485233, 1.215897505094087, 0.28561868754557707, 1.1255214337331036, 0.2716913371345464, 1.030900839672786, 0.2571983067377617, 0.9330712882989015, 0.24184034171920196, 1.1521177775877345, 0.36873988809821695, 1.0699609129128091, 0.3508208799888231, 0.9838541344890434, 0.33204004955724326, 0.8946241473814437, 0.3123800607019776, 1.0096723942589452, 0.43664017547194817, 0.9332535560402271, 0.41320013054027077, 0.8529762644477507, 0.38929318611314456, 0.8787810203490878, 0.4992609208501592, 0.8082964051143625, 0.47137755736401116, 0.761462581999484, 0.5570240879581907])
    SymCubatures.setweights!(cub, T[0.00014701499050156892, 0.0007877794452320271, 0.00186778956223647, 0.0032420816816595734, 0.004726828091032836, 0.006103737018443741, 0.0072485765064839724, 0.007948128447926169, 0.008131995177044438, 0.007845132829365697, 0.007155544724002234, 0.0059937564425873845, 0.000340730168570429, 0.0005265637910118002, 0.0006962705839636064, 0.0008444495795432031, 0.0009688712761667022, 0.001061560085290703, 0.0011252592270516832, 0.0011561541323693947, 0.0011575314926823327, 0.001136704875397184, 0.0010815129049989354, 0.0012142478706855865, 0.0016020532203012547, 0.0019383693390481675, 0.002219546731092341, 0.002429504520319999, 0.002576524901001423, 0.0026476764917516883, 0.0026542068708373633, 0.002602200938126511, 0.002484845735723617, 0.002461279864178378, 0.0029729379721938564, 0.0033955204618935095, 0.0037104319470256278, 0.003931049894003802, 0.004037554049667483, 0.0040525864608907, 0.003971575036750879, 0.003810562815087516, 0.003913984996016317, 0.004460032985161071, 0.004869296171086892, 0.005150446619154343, 0.00528491580792503, 0.0052974874682032965, 0.005185045819103452, 0.004983873053729913, 0.005377775766257213, 0.00586888923498486, 0.006193148447066716, 0.006347975690015123, 0.006346426127691523, 0.006207787375075502, 0.00595435184536097, 0.006655889586941313, 0.007010959103050481, 0.007185974598089423, 0.0071689631218275615, 0.007017896790967016, 0.006722054117686235, 0.007609938902965881, 0.0077821111947968084, 0.007731772399919157, 0.007537202009485338, 0.007233494436474263, 0.00809083674129125, 0.00804072794863499, 0.007801396766394829, 0.0075095921323745886, 0.008045416721684484, 0.007815398410524371, 0.007593576543280526, 0.0075870390565284055, 0.0074111184072079, 0.006933071240492395])
    cub_degree = 46
    tol = 5e-14
  elseif q <=48 # 469 nodes
    SymCubatures.setparams!(cub, T[0.004374403267739703, 0.022953071777193152, 0.0557316670528537, 0.10135012284992394, 0.15747172272636498, 0.22163772526216516, 0.2907341880124402, 0.36210791796945235, 0.43372210375817727, 0.5023814781001197, 0.5657218585564955, 0.6214518599648744, 0.998531948124674, 0.9922773097540275, 0.9810772208385705, 0.9650828676697426, 0.9446091032558428, 0.9200370099439127, 0.8917692264854484, 0.8602550479954347, 0.8254519513789932, 0.7875645456071046, 0.7508993242301097, 0.7141632878035474, 1.9726433296991923, 0.004369504117273375, 1.9394344172309992, 0.004344342270942659, 1.8920650458897201, 0.004288824439387421, 1.8311564290504827, 0.004227793304139139, 1.757679387112075, 0.004047346864599653, 1.6724866415003747, 0.00398275573614045, 1.577034597434429, 0.0038320146817607094, 1.472674794701072, 0.0036712366804263876, 1.3609458616689176, 0.0035026614290687547, 1.2435116722278512, 0.0032950117160055205, 1.1220773754499702, 0.003101071844859245, 1.9210844902978492, 0.02281276073229102, 1.8741281385609, 0.0225279113432623, 1.8138287582826333, 0.022178289875991405, 1.7414896593506544, 0.021300112840920957, 1.6572153830274599, 0.020874829816517913, 1.562929063425372, 0.020138008718889253, 1.459852863720225, 0.019258131164613692, 1.3494274528732588, 0.018403055773944046, 1.2335734823722433, 0.017332291930225428, 1.113909036741237, 0.016345419279426236, 1.842344686799491, 0.05505139284006577, 1.7831773337207246, 0.05409699398719893, 1.7125606593397371, 0.05219621223900051, 1.6301489703688865, 0.050863135686446036, 1.5380280797931842, 0.049223795294764804, 1.437770395854461, 0.047010369439309585, 1.3301989856670946, 0.044982162023378035, 1.2173156747805316, 0.04249483144681671, 1.1003301106664782, 0.04013018380394932, 1.7399027251826182, 0.0994253946855884, 1.6715271575317256, 0.09638062007691398, 1.5921416406054052, 0.09343397780688636, 1.503019417586394, 0.09056959688034898, 1.4064787898011442, 0.08654940899876169, 1.3026925522013237, 0.08276525354170834, 1.1936521520246186, 0.07851072037772376, 1.0805841204412747, 0.07418796869479911, 1.6193969178344307, 0.1531950489891539, 1.543733514180707, 0.14805563088116708, 1.4581917131088822, 0.14332372974776675, 1.3656607666581826, 0.13730369839078765, 1.2666692271289983, 0.13106828998295666, 1.1625679726401248, 0.12470072069431574, 1.0548093653373705, 0.11793915139392415, 1.4852963801979078, 0.21409263836064696, 1.404765915122779, 0.20658905170623432, 1.3170688789960248, 0.19845945554509428, 1.2238306588047256, 0.18930759621195728, 1.1256307344238297, 0.18009423884214207, 1.0239666630251159, 0.1704468645308705, 1.3437035442712164, 0.27970220430947373, 1.2613368557028113, 0.2688249960325205, 1.1742049180525391, 0.25669862560537643, 1.0831638991424783, 0.24381378107007384, 0.9887536192268952, 0.23055421972856216, 1.200177662489211, 0.3475236099532179, 1.1192650515387692, 0.3320935860455671, 1.0353322131552405, 0.3153065859857094, 0.9490326060507615, 0.297850530927641, 1.0595598777836852, 0.4141853959152407, 0.983155422949612, 0.3934744194779852, 0.9046152707275072, 0.37264773092698006, 0.9278462055675119, 0.47695350666219766, 0.857572437952181, 0.4527600028360214, 0.809181875303464, 0.5341339127052507])
    SymCubatures.setweights!(cub, T[0.000125869349252052, 0.0006721911819214603, 0.0015994229574246585, 0.0028016760184970163, 0.0041001497002233715, 0.005417322304210139, 0.006532472836351598, 0.0072832612603159355, 0.0076712855473958855, 0.007658350557729647, 0.007322783978684914, 0.006601853862341765, 0.0009339687196601844, 0.002140286959911097, 0.0032966888431143417, 0.0043307231452388495, 0.005179749540836551, 0.005844200456612081, 0.0063127126245901026, 0.0065949031293567354, 0.006793591352762243, 0.00658828685541234, 0.005475953847764807, 0.00587259566632568, 0.0002910863856046578, 0.00045027076482358046, 0.0005974692916055362, 0.0007312764390748318, 0.0008270968018952332, 0.0009252949807960514, 0.0009857008482064532, 0.0010206126278872682, 0.0010331734563762344, 0.001013296468478407, 0.0009785565690051535, 0.0010378095346639754, 0.0013757261409999058, 0.0016773811864634568, 0.0019073755072449463, 0.0021171469213993274, 0.0022668568341936554, 0.002340334445373825, 0.0023734005700355187, 0.0023298220928557527, 0.0022517886293027894, 0.0021185494456768816, 0.0025723779017676977, 0.002943329470380113, 0.0032302911180665094, 0.0034643641188970946, 0.003569021136541826, 0.0036191652357721355, 0.0035813390086957162, 0.0034720649862441267, 0.003389250667362178, 0.003891321147209804, 0.004238198124643864, 0.0045297906769193345, 0.00468585199719133, 0.004731205538090608, 0.004710246886454709, 0.0045703548045643696, 0.004721384909679583, 0.0051514023334920675, 0.005455495917522983, 0.005663296707058535, 0.005680371395478073, 0.005633202323689877, 0.005474040282214396, 0.00592215551862148, 0.006206413431880633, 0.0064341306392687305, 0.006467866406142498, 0.006347718608960759, 0.006147179848524022, 0.006851695001210488, 0.0070168449990386265, 0.0070494972337832845, 0.0068761505597732244, 0.0066067141030988095, 0.00745571169628502, 0.00745051099816164, 0.007258750675111826, 0.006966356062579333, 0.007643593151542233, 0.007411240553812757, 0.007165547303985224, 0.0073718550392302675, 0.006891100968264401, 0.0064842147353374, 0.006633256925121181])
    cub_degree = 48
    tol = 5e-14
  elseif q <=50 # 507 nodes
    SymCubatures.setparams!(cub, T[0.004044718961736702, 0.021277438243104217, 0.05179418541321082, 0.09423298423786153, 0.1464829057658257, 0.20681803252806388, 0.27265532017516936, 0.3413678195787203, 0.4104837860592545, 0.47739669610878727, 0.5397233146740809, 0.5954365007943827, 0.6439890566340963, 1.9746892384174533, 0.004049773435059005, 1.943940106111981, 0.004034208542846271, 1.9000342095840939, 0.003987320093424047, 1.8435220672413613, 0.0039026511715255963, 1.7751055058002516, 0.0038064042019867617, 1.6956799603508605, 0.0037130437960210014, 1.6063651605616827, 0.00357789595075233, 1.5083510168490863, 0.0034413935134932493, 1.4030082724049986, 0.0032900199192522545, 1.291797183389569, 0.0031371102676219574, 1.176286905430423, 0.002988582921680949, 1.0581532669479543, 0.0028178799599012767, 1.9267830885280672, 0.021190815592442656, 1.8831694817989015, 0.020940572711769325, 1.8272140082249226, 0.020497428133272638, 1.7595284572195538, 0.019994395768968214, 1.6809563736717255, 0.019494862216213008, 1.5927471018227772, 0.018792391459090868, 1.4959319633903185, 0.018078550535997964, 1.3919265085486356, 0.017280735137043295, 1.282132315241171, 0.016490025842407677, 1.168131174106432, 0.015684267685981208, 1.0515700343161256, 0.0148011568163126, 1.8534207289685878, 0.051168661307779806, 1.7984850151328233, 0.050094181260520156, 1.7321032505895457, 0.048872366944745116, 1.6551193359916498, 0.04761627251777919, 1.5689608411757814, 0.04593205411492247, 1.4744480304627752, 0.04420064405080811, 1.372981925641073, 0.04225235825383797, 1.26571250811385, 0.04035058459591831, 1.1542365467457913, 0.03831001434898953, 1.0400612852772424, 0.0361770567222787, 1.7582015239596407, 0.09228287474030279, 1.6937528123742993, 0.09004581090807642, 1.619030307885766, 0.08765184197102374, 1.5355453667606331, 0.08461935246923237, 1.4440276712427818, 0.08144644213231267, 1.3459219119219088, 0.07788223217946309, 1.242143986429625, 0.07439563748463518, 1.1343159957666031, 0.07056864771647861, 1.0236300041358444, 0.0666510186873882, 1.6451458815534192, 0.14295228111884134, 1.5733897178380805, 0.13902376552022985, 1.493182982458744, 0.13431376763603725, 1.4052971398328422, 0.1292920838669229, 1.3111897029286683, 0.12370808132004205, 1.211713321217047, 0.11811433780639219, 1.1084487396721292, 0.11205914632615457, 1.002386178122539, 0.10582656656987029, 1.5183110648654363, 0.20097886960357714, 1.4421537504357131, 0.19428037180045538, 1.3587892918808757, 0.18705177589871633, 1.2695802627476036, 0.17912790231811904, 1.175416782264303, 0.17089886348310548, 1.0774865049811662, 0.16226095108569985, 0.9767982273446288, 0.1532213953963404, 1.383309778723518, 0.26364530980341905, 1.3051295489837034, 0.253877364709623, 1.2214061802645826, 0.2433146879398133, 1.1333181853672736, 0.2319963889121938, 1.0415706352210012, 0.22039596361951866, 0.9471301967294229, 0.20817942277297136, 1.2449286526006953, 0.32874192367952476, 1.1673850616908115, 0.31524334448050784, 1.085992410806533, 0.3005586016806002, 1.0011782502939752, 0.28550299116707917, 0.9137573124973838, 0.26982553630431216, 1.1083762755889337, 0.3937141913625807, 1.0343052524256207, 0.37560222628099105, 0.95721754670282, 0.3566436797649436, 0.8773240368896642, 0.33721290625263967, 0.9788820161819283, 0.4558378507204741, 0.9101700578008799, 0.4328741917414395, 0.8383231072686156, 0.4095337974167157, 0.8601995296542159, 0.5129451370159395, 0.796871210572556, 0.4860314212540999, 0.7537704058509973, 0.5654356776216062])
    SymCubatures.setweights!(cub, T[0.0001076268438018305, 0.0005779091452547007, 0.0013835558112076048, 0.0024271064098401888, 0.0035836783236265564, 0.0047377292926316615, 0.005735032465255828, 0.006499373245997843, 0.006920287906882842, 0.006938043163940708, 0.006618480181691711, 0.006003552662172691, 0.00509582693880631, 0.00024962730878997605, 0.0003872809762236535, 0.0005151361514687561, 0.0006273996823026902, 0.0007240930022477023, 0.0008061327314526002, 0.0008624969144831796, 0.0009005373175571087, 0.000916691578494102, 0.0009152154806900204, 0.000898062987175781, 0.0008597651563497489, 0.000894909029858437, 0.0011880104327640312, 0.0014450002980487314, 0.001666118212314397, 0.0018516360759718966, 0.001981476818458488, 0.002068428126156071, 0.0021037372021230185, 0.002101405978692136, 0.002054725443304689, 0.0019694053161145864, 0.0018335295340564453, 0.0022272452370833196, 0.002563377135027804, 0.0028386766160081876, 0.003035976863876513, 0.00316672270324564, 0.00322087765392896, 0.003221023373894537, 0.003145131287989664, 0.0030183849027073154, 0.0029478166077040384, 0.003389853921604954, 0.00374343209873324, 0.004002792532803811, 0.004168939579356466, 0.00423847874396161, 0.004231699800365234, 0.004134858474492832, 0.003970578682354509, 0.004121729486265713, 0.0045442132480069276, 0.004858866799362313, 0.005054000924535614, 0.005136092999862892, 0.0051067822294558396, 0.004993455024731368, 0.0047925204959277275, 0.0052144039633039395, 0.005567232440660688, 0.005781221551332217, 0.005872455161081112, 0.0058228544332916505, 0.005699147651981842, 0.00547557603567639, 0.006115456881675805, 0.006344694682869211, 0.006436161637037278, 0.006366819835006978, 0.0062138302348848495, 0.00598473869301454, 0.006723859766757872, 0.006804646138116097, 0.006731735829722992, 0.006533619220849005, 0.006292991678671555, 0.00696295061394431, 0.0068892045840643745, 0.0066740822453501564, 0.006428539545048986, 0.006829935714412209, 0.006631617795085216, 0.006426592408950292, 0.006399044774873687, 0.006248856928256274, 0.005852612513473457])
    cub_degree = 50
    tol = 5e-14
  elseif q <=52 # 547 nodes
    SymCubatures.setparams!(cub, T[0.003758951652370829, 0.01976390747738702, 0.04816712260884112, 0.0879117814003429, 0.13700111503706267, 0.19366537341646545, 0.2562395628923429, 0.32179698640144994, 0.3885720185341922, 0.4540902094393828, 0.5165717297361165, 0.5744527435225399, 0.6247106136600531, 0.998716044980445, 0.9932939332297255, 0.9836800523171771, 0.9699186175722885, 0.9521385655073502, 0.9307431750989825, 0.9060868034464943, 0.8783179111122528, 0.8477770040047248, 0.8146372757796494, 0.7792369699977398, 0.744970727635657, 0.710080444933196, 1.9764787053606385, 0.0037607067468720333, 1.9478875897329504, 0.0037491262677770997, 1.9070343012623479, 0.0037150731431686537, 1.854401143950142, 0.0036498755959552694, 1.7906236572436969, 0.0035461196571300794, 1.7164383530473624, 0.0034477772164228718, 1.6327473446782728, 0.003376289573938611, 1.5406788160906626, 0.00324650956322026, 1.4413526937356103, 0.0031084138164049905, 1.3360587851134573, 0.002955888720344187, 1.2261757486526688, 0.002821004889715533, 1.1132438802991516, 0.0026401939299320505, 1.931954468206645, 0.01969819181402893, 1.8913183367909774, 0.01951658838094576, 1.8391134256934112, 0.019168953383907762, 1.7759994246472905, 0.018635936080135224, 1.7025637864829317, 0.018118435766500105, 1.6196426099198693, 0.017720968261489493, 1.5286574387810299, 0.017054314322793498, 1.4305899475005317, 0.016331317043127846, 1.3267029355458086, 0.01555681767292029, 1.2182518243634297, 0.014821879985285275, 1.1067076833761598, 0.013970609253012887, 1.8635037915312835, 0.04771089538169327, 1.8121308144014212, 0.04684524137974735, 1.750201502302064, 0.04558876201316976, 1.6781607045150002, 0.044323449774482215, 1.596761346331105, 0.04325913754820014, 1.507663595599397, 0.041695164723951406, 1.411694513428896, 0.03993439551621881, 1.3099550669698348, 0.038124483365686024, 1.203783766498716, 0.03627147196340335, 1.0945721104670225, 0.03442012904211177, 1.7741296922854994, 0.0862950612509703, 1.713911873556373, 0.08408195253851267, 1.6438957786644532, 0.08176100394502882, 1.5648346658243215, 0.07957709857862617, 1.4783913299699707, 0.07684269036487414, 1.3855569100005696, 0.07363765843327622, 1.2870660311588022, 0.07042384218264075, 1.1841771380209702, 0.06703844657972895, 1.0779324790868672, 0.06374250662461178, 1.668058549512383, 0.13365751689896516, 1.6006974806403709, 0.1300103769855037, 1.524762276376498, 0.1261780291781443, 1.441420139605095, 0.12199933084051066, 1.3520921110463044, 0.11703625747285387, 1.2573446995274848, 0.11195826454658478, 1.1583068974067587, 0.10680399032325427, 1.0562363862051598, 0.10136597186501375, 1.5482706292476942, 0.188434242944157, 1.475912803175389, 0.18247158062992327, 1.3961451914495482, 0.1764123722567676, 1.3107915828027503, 0.1695011131274835, 1.220704769478125, 0.16206649408986715, 1.1264357292767713, 0.15483708904335514, 1.0296657484560736, 0.14674026374889432, 1.4195716477764322, 0.2479107023877979, 1.344516639928882, 0.2394166480985748, 1.2638753536732947, 0.23033666348363263, 1.179105219688018, 0.22025933794419206, 1.0905006073029209, 0.210174583539504, 0.9993294453083451, 0.19934295066759378, 1.286608876909394, 0.3103496485165898, 1.2111614639648325, 0.2984871472166665, 1.1317725151242848, 0.28571868265511197, 1.0496739441664726, 0.27196218061385685, 0.9649752108556312, 0.25817527643899374, 1.1537372945704247, 0.37322446065424, 1.080282426202634, 0.35745072229957725, 1.0047436829577696, 0.34013519125794883, 0.9270627111742141, 0.3228091874305521, 1.0254847401456428, 0.43453453273823667, 0.9565766121174077, 0.41406749933079134, 0.8858769855376684, 0.39367109482276397, 0.9053700140319799, 0.49222672608068613, 0.8420647810574794, 0.468652373343398, 0.7971327379098591, 0.5441243224860327])
    SymCubatures.setweights!(cub, T[9.295951697904649e-5, 0.0004989374886136669, 0.0011984066809653002, 0.0021177136005076165, 0.0031404704636416314, 0.004196118692130455, 0.005154151674612369, 0.005923189991227376, 0.006427402536074561, 0.006620600681396203, 0.006535111825744637, 0.006214076452830316, 0.00555098627065631, 0.000754833948633958, 0.0017180226740439856, 0.0026219318002614566, 0.003488822457729071, 0.004233127978072337, 0.004789886020422311, 0.0052497224245022876, 0.0055773856454628444, 0.005720475867444771, 0.0058155657219709864, 0.005583304858718017, 0.004694858695648014, 0.00502553494425424, 0.00021548909402453194, 0.00033471138459313123, 0.0004466905946798245, 0.0005466591804992865, 0.0006295710990616532, 0.0007001933655914863, 0.0007631482560530689, 0.0007994195607076886, 0.0008183322684455302, 0.0008186161481881789, 0.0008085587639153914, 0.0007738666679707453, 0.0007738740146904678, 0.001031322578014772, 0.0012601255764575373, 0.0014516895673237357, 0.0016130956341056373, 0.0017521579592062945, 0.0018366173573793182, 0.001878338078663883, 0.0018834826491620903, 0.0018551891372374728, 0.0017990728332896392, 0.0015944544164585256, 0.0019444795062929273, 0.0022405222759333365, 0.0024857306949785537, 0.002684899301191892, 0.002820190417165637, 0.0028844708000039785, 0.0029005563085025944, 0.0028521158296563624, 0.00277899184478408, 0.0025790111719392805, 0.002973193880376287, 0.0032944980535624903, 0.0035339347608366613, 0.003713819216333963, 0.003798377544795061, 0.0038211471972431793, 0.0037796359175229933, 0.0036706697940628668, 0.003625730651980261, 0.004018253846298645, 0.004292416689866055, 0.004506903393623427, 0.0046201510913751005, 0.0046295667201474155, 0.004592031366371901, 0.004427732921957754, 0.004653860646272044, 0.004966314711557351, 0.005182371660107702, 0.005309669595200183, 0.00529663866527238, 0.00522042292619335, 0.0050448141860131015, 0.0055078292817087385, 0.005715857477662046, 0.005837797133444532, 0.005832536550717394, 0.005677212763868368, 0.005520535993578779, 0.006151493298238659, 0.006232434673463063, 0.006224093268902536, 0.006016486303719636, 0.005798770207884925, 0.006500832548196402, 0.006445586303007323, 0.006277757217160346, 0.006012193933105497, 0.006540003067409229, 0.0063487297690334125, 0.006132721023388242, 0.006291762129550324, 0.005845314194743745, 0.005482174330493251, 0.0056496758078770325])
    cub_degree = 52
    tol = 5e-14
  elseif q <=54 # 588 nodes
    SymCubatures.setparams!(cub, T[0.00350244032024077, 0.018421339176714558, 0.04493022524744672, 0.08202950178154739, 0.1280806026069353, 0.18169859863993915, 0.24096303061449909, 0.30380610240956146, 0.3682130122048638, 0.4320755646719692, 0.49331293949517396, 0.5503769375708557, 0.60098350698178, 0.6456505842098564, 1.9780801892249344, 0.003504648554066982, 1.951420964751685, 0.0034953282046491043, 1.9133029630675147, 0.0034630111366873873, 1.8641443624683214, 0.0034003632270480045, 1.804471116134648, 0.00332679682312236, 1.7349541647551523, 0.0032496026312659933, 1.656411908513708, 0.0031548902608898536, 1.5697598916720024, 0.0030510928467482983, 1.4760349667128858, 0.002933664780348348, 1.3763533842542321, 0.002820690011065584, 1.2719450209233472, 0.0026981233960854593, 1.1640897262050973, 0.002578586099768896, 1.0541440792707841, 0.0024416761871348443, 1.9365544470765632, 0.018368339786844004, 1.8986321833905457, 0.01819509756064234, 1.8498583594974716, 0.017867411728754903, 1.7906982841071686, 0.017481484052258467, 1.7218039229994924, 0.017072886642752003, 1.6440456231943337, 0.016575509672505897, 1.5583044637864107, 0.016032157726677113, 1.4656288112705274, 0.015416601703411211, 1.3670603943931328, 0.014820970544343297, 1.2638509575665113, 0.014175339851945554, 1.1572155916346902, 0.01353559356281375, 1.048507391590455, 0.012816539647706857, 1.872648741325593, 0.04449422542012841, 1.8246156593959915, 0.04369990737084967, 1.7664263873260595, 0.04275945034358956, 1.6987216922652868, 0.041748109584657335, 1.622435704382958, 0.04053552375202703, 1.5383830960031306, 0.03921398185170152, 1.4476016080584164, 0.037719320836221926, 1.3509990945928982, 0.036249495540433295, 1.249830784636372, 0.034667659036453975, 1.1452419621551724, 0.033057968748311274, 1.0385522894855068, 0.03130432554534254, 1.7890808327675427, 0.08058538409126627, 1.7323442102912938, 0.07886201533279256, 1.666361028161139, 0.07696665847935903, 1.5921086082424358, 0.07474173424410076, 1.510381014400283, 0.07231448966608021, 1.4222162121569852, 0.06959440403665212, 1.3284045598891068, 0.06684515880552347, 1.230104681838236, 0.06392913141200683, 1.128389645801391, 0.060880228790495344, 1.0244384420450323, 0.057666420903245384, 1.6890450399414763, 0.12536742672140042, 1.6253151296931667, 0.12230078113861666, 1.553565440242029, 0.11877978130889413, 1.474608144666234, 0.11492263194377508, 1.3894861514708639, 0.11066475167998159, 1.2990530395294277, 0.10622097770929918, 1.2043501940528956, 0.10158226944171644, 1.1064198170521535, 0.09667577328596748, 1.0061615279209968, 0.09162238601331595, 1.5756633546762984, 0.17718588624314524, 1.5070085916783018, 0.1720954567801407, 1.431486453974145, 0.1665114832282442, 1.3500751101877173, 0.16043968448882207, 1.263733072528653, 0.15393068826508324, 1.1733065081370972, 0.14717277287253802, 1.0797887452052077, 0.14009132329381088, 0.9839486550010595, 0.13284390153264083, 1.4531214089247064, 0.23403123984643018, 1.3816089543982912, 0.2264510465694623, 1.3044697388118527, 0.21829067728124435, 1.222823317536553, 0.20943079379779414, 1.1374050194184802, 0.2001719886342554, 1.0489335161705138, 0.19067773595991727, 0.9581364043215416, 0.18088449928628567, 1.3254631540604278, 0.29398230171735434, 1.2531942238362281, 0.28345861941417827, 1.1767103359829048, 0.27203760385121806, 1.0968038593362115, 0.25996212892646053, 1.0139699437807796, 0.2477637227474629, 0.9289573323381753, 0.2350955624558308, 1.19674934209386, 0.3550480169141854, 1.12601563034467, 0.34091647937119646, 1.0521709548720626, 0.32588267868641885, 0.9755707250513976, 0.31055322308958905, 0.8968664371613108, 0.2946801364718864, 1.0713529243922566, 0.41505508932846746, 1.0040740480395205, 0.39709090769464955, 0.9344523223433941, 0.3782298668485425, 0.8624192119171765, 0.35887985353618296, 0.9528217354477, 0.4723783756868767, 0.8906847310186503, 0.4498634464359883, 0.8258069552906485, 0.42709548106303963, 0.8445284335678638, 0.5243763697306262, 0.7870583565740741, 0.4987431277732318, 0.7472258772069508, 0.5727578504589478])
    SymCubatures.setweights!(cub, T[8.07115369024744e-5, 0.0004336379196902444, 0.0010439069799027334, 0.0018484372566179676, 0.0027616581765375404, 0.003705721289089597, 0.004576735709506493, 0.005299820986587521, 0.005807356650199201, 0.006032656838634098, 0.005971982324238649, 0.005647654450089076, 0.00511976850034328, 0.004375576450413571, 0.00018718941480647068, 0.00029105185556078635, 0.0003886940302216334, 0.0004760296304160933, 0.0005528240377621099, 0.000618841307270086, 0.0006703588987439086, 0.0007079114134470083, 0.0007298022349210083, 0.0007403814128394303, 0.0007364974654232674, 0.0007220362827315001, 0.0006924818627258663, 0.0006733441674098321, 0.0008979763711588099, 0.001098786512291473, 0.0012749245763717938, 0.0014254469586226636, 0.001542965033600552, 0.0016285623251999324, 0.001678176143151243, 0.0017009737911709356, 0.0016910585884116263, 0.0016547279546282967, 0.0015871888676452256, 0.0013898756206274926, 0.001698765346333203, 0.0019681801678768543, 0.0021957147906725635, 0.002373929946064694, 0.00250369075235151, 0.0025800980556875985, 0.002612095977680143, 0.0025965577010322737, 0.0025351478443660642, 0.002433365583355748, 0.002258601836947142, 0.0026149543687504135, 0.0029113659550678763, 0.003143746452260111, 0.0033106279453445535, 0.0034112268110032496, 0.0034474645011802457, 0.0034274677705350594, 0.0033440249206455118, 0.003215323366604009, 0.003198734177941799, 0.003558289831846611, 0.0038407828590893266, 0.004040215705032703, 0.004160508374714369, 0.004192780443738319, 0.004159493566706491, 0.004058984882325698, 0.003909911894756441, 0.004118870733734338, 0.004440587095020341, 0.004665454623041069, 0.004800260058769973, 0.004831630993465763, 0.004781725245787788, 0.0046714225793223095, 0.0045017404834723765, 0.004929838493774656, 0.005176272872744049, 0.0053175885770557286, 0.005347439452070921, 0.005278928945603864, 0.005160663811965364, 0.004976239534262198, 0.005556832893421344, 0.005700823094224183, 0.005733190681883684, 0.005650878781522794, 0.005502402705220109, 0.005302759310107801, 0.005935927960324273, 0.005961955247273043, 0.005883406484364541, 0.0057000689049504566, 0.005482186460026362, 0.006028519818900386, 0.005949336753433454, 0.005757023456628438, 0.005551255237299216, 0.005856557547361985, 0.0056681269865943105, 0.0055152972290721716, 0.005442783954590249, 0.005332909232938009, 0.0050094559716599265])
    cub_degree = 54
    tol = 5e-14
  elseif q <=56 # 631 nodes
    SymCubatures.setparams!(cub, T[0.0032850508915074927, 0.017149044330110365, 0.04195314617597781, 0.07701223792385879, 0.12017351090308624, 0.17063858455862427, 0.22710367573277718, 0.28696915557124647, 0.348719315051913, 0.41095483936677185, 0.4711342547613978, 0.528750292759843, 0.581447142031664, 0.6283618486438389, 0.998876014325462, 0.9941176439117984, 0.9857115907620959, 0.9738172925446147, 0.9584437206515887, 0.9395936461551411, 0.9176255492107568, 0.8930435884251994, 0.8660686136480654, 0.8368618262061288, 0.8053525269909534, 0.7717855634752475, 0.7395360051450671, 0.7076760944115794, 1.9794668762225849, 0.003261697708542561, 1.954459750595459, 0.0032610502448526093, 1.9187039977468936, 0.0032491438496285202, 1.8726206229027749, 0.003190185979512375, 1.8166972903756071, 0.003106802631890957, 1.7515286406431416, 0.0030337995970748373, 1.6778166963786685, 0.002995711206470109, 1.5964709590847823, 0.0028708595290183212, 1.5082405304525082, 0.002781458823316524, 1.4140972960199087, 0.0026576717927226834, 1.3149624277578753, 0.0025487258858050616, 1.2118952541997439, 0.002452376668571891, 1.10610193004731, 0.002352199914312819, 1.9408743831312123, 0.017142327844736363, 1.9054061465489203, 0.017076058064609302, 1.8598326157518454, 0.016766872173500407, 1.8044892659202423, 0.01633789756430171, 1.739780217399211, 0.015959128074729174, 1.6663170326292862, 0.015715956901213583, 1.585366080459248, 0.015106845444449143, 1.4974371968657167, 0.014599764251437724, 1.4037835934934464, 0.01398066262940157, 1.3054378095863561, 0.013389489840010023, 1.2035644374265193, 0.012888481374443192, 1.0994469690111923, 0.012305444404913447, 1.8808158486543385, 0.041765344273659515, 1.8358624850922693, 0.04100004593078519, 1.7814178481130984, 0.039991070215110974, 1.7177701928691653, 0.03908567197198691, 1.6456058882924218, 0.03834529508921606, 1.5664169381100725, 0.037022397566114175, 1.4805188159321103, 0.03569959895264029, 1.3888988777041928, 0.03428317013503656, 1.2925178620624456, 0.032811223595792735, 1.1922941856175568, 0.031536762136490856, 1.0896988557979124, 0.029990556582121887, 1.802004698742429, 0.07559045485776231, 1.748784968108282, 0.07382763202006438, 1.6865504169787282, 0.07217770288679301, 1.6161780914120716, 0.07052212387406893, 1.5389436054014607, 0.06836016158737111, 1.4555361229210562, 0.06584668654111639, 1.3665469882728671, 0.06332762125379716, 1.2729624256468093, 0.060675605922870156, 1.175536272634142, 0.05812309477435778, 1.0755291921156231, 0.05524311658865132, 1.7082404048420243, 0.11754764702007899, 1.6480265348353567, 0.1149162706753731, 1.5800588273741132, 0.11191241027643839, 1.5049534983298718, 0.10871541856396404, 1.4240260815347676, 0.10479384243772194, 1.3377004888782629, 0.10075156756170192, 1.2469216241037633, 0.09668358282891103, 1.1528909697121257, 0.09231923932900193, 1.0564240437888488, 0.08783625367429619, 1.6007714393712023, 0.16676134212391155, 1.5355977938093985, 0.16208162332962983, 1.4633300928100392, 0.15743723241171698, 1.3855572508675493, 0.1520792296663921, 1.303033721666484, 0.14617463581610848, 1.216164323541379, 0.1403389522216804, 1.1262384601706967, 0.13395278202271524, 1.033762538071122, 0.1274494353666316, 1.4840966774953002, 0.22062319789064125, 1.415612968464832, 0.2139907512017479, 1.3413828860216381, 0.20705857768363228, 1.262880925774777, 0.1991743549287692, 1.1805504119057715, 0.1909521506998954, 1.0951083547576952, 0.18253862109344265, 1.007346907893327, 0.17369850230380293, 1.361802783907403, 0.27788938182111017, 1.2915709706911154, 0.2688699372663649, 1.2171788898088836, 0.25898594249197593, 1.1398119388481713, 0.24797795823336885, 1.0594649566823697, 0.23711015977577954, 0.9770394167346841, 0.2259761969353087, 1.237417126826401, 0.3370053557620933, 1.1677741826083765, 0.32470819188207617, 1.0953082403360888, 0.31115083217423, 1.020581216072317, 0.2970622958488072, 0.943853065652581, 0.2832050769300032, 1.1144293245462644, 0.3954034204062702, 1.0471192132735505, 0.3794888023157476, 0.978270875246148, 0.36229353986635715, 0.9081477533892129, 0.34481710028378676, 0.9966845462727, 0.45214958701591357, 0.9334534681579171, 0.4321377048747554, 0.8695328473138613, 0.41160555386553144, 0.8863075932352416, 0.5048612983754374, 0.82837710534372, 0.4820399176270046, 0.7876608312311865, 0.5530331911861013])
    SymCubatures.setweights!(cub, T[7.100169588939963e-5, 0.000376046932861409, 0.0009117526188051608, 0.0016311489721039236, 0.002434991693155489, 0.003301205605260779, 0.004109391622090789, 0.004822284866748426, 0.0053403372873247635, 0.005658336442045272, 0.005730820537355288, 0.005649311085433417, 0.0053194827966055545, 0.004833476076309558, 0.0006194320382113281, 0.0013998702899521013, 0.002134056900763214, 0.002822268071222333, 0.0034574492732770434, 0.004021614066423628, 0.004434284989024414, 0.0046982390182267846, 0.004872878856618303, 0.004959323715001552, 0.005029796371965175, 0.0047807816313075215, 0.0039924808882005535, 0.00429286485147744, 0.00016339146053141726, 0.00025462658545917543, 0.0003419030584814374, 0.0004186488850838554, 0.00048399785809808196, 0.0005418311678069962, 0.0005970475930865998, 0.0006263234329435663, 0.0006520992106133052, 0.0006607721228979685, 0.0006626096873404492, 0.0006587099946071785, 0.0006434994849222338, 0.0005859550305232484, 0.0007868372683943916, 0.0009644415930968323, 0.0011180634784487195, 0.0012548684322267096, 0.0013779814511527162, 0.0014551384999716927, 0.0015077079612075442, 0.0015301771282930819, 0.0015248626629226843, 0.0015086801168033262, 0.0014574407981351484, 0.0012211950865183073, 0.001493250790406064, 0.0017302280097191499, 0.0019379948097911813, 0.0021083221796370268, 0.002237632882767586, 0.0023109247874665132, 0.002355717303797377, 0.0023541007696129415, 0.0023233177874984203, 0.002242134807703568, 0.0019918639448520234, 0.0023092769485922463, 0.002580870233958095, 0.0027863809536801897, 0.002963129408134213, 0.0030585776966189403, 0.003111233002560246, 0.0031216400033389556, 0.003063450873182867, 0.0029761301203929187, 0.002830950998746558, 0.00316418977056122, 0.003411520987405255, 0.0036246503811181627, 0.00376045757381402, 0.0038099337504378574, 0.0038168912591330004, 0.0037295524445898525, 0.0036174226366275895, 0.003686754570937482, 0.003973622915587414, 0.0041869777170352685, 0.004349395259094731, 0.004402418013449365, 0.0043850457921484854, 0.004317843133856076, 0.004179452141214162, 0.004440487041624129, 0.004656750604623788, 0.004827523608725018, 0.004903979076040785, 0.004844485843584506, 0.004775439352292027, 0.004642349012712461, 0.005059717244975912, 0.005205736768025572, 0.005283818110771739, 0.005221526814686268, 0.005077438870123536, 0.004942425486982535, 0.005487394287088306, 0.005530085810205154, 0.0055085157254242635, 0.005313407684456508, 0.005081378743724835, 0.005674263968739643, 0.005638840116339741, 0.0054735413084453025, 0.0051867638589929805, 0.005676045283695648, 0.005493629561260174, 0.005271560606417316, 0.005393656855657771, 0.005031618449197114, 0.004730521765396585, 0.004835401003558482])
    cub_degree = 56
    tol = 5e-14
  elseif q <=58 # 675 nodes
    SymCubatures.setparams!(cub, T[0.0030639025446898512, 0.0161020263848769, 0.039320733306980056, 0.07203116773967717, 0.11288381969374145, 0.16070788111695364, 0.21409254585455156, 0.2715722935605522, 0.331244631264744, 0.39145750224987563, 0.45061279087490347, 0.5069433318043616, 0.5591732563595755, 0.6060125517483922, 0.6470891344142562, 1.9808245751105804, 0.0030624197231182543, 1.9574801488737783, 0.00305583436243931, 1.9240602898937247, 0.0030349205645119905, 1.880892553764992, 0.0029890893885460164, 1.8283845641775933, 0.0029296087543714665, 1.767047297374496, 0.002865454456594171, 1.6974974290074887, 0.0027988317510518855, 1.6204618902423829, 0.00271576400555429, 1.5367364955246579, 0.0026243968696535815, 1.4471859417198287, 0.0025476179690313355, 1.3528000907505866, 0.002434635490791717, 1.2545358851010897, 0.0023407199118293635, 1.1534604549166914, 0.0022370052513559983, 1.0506662839492324, 0.002115910055869704, 1.9445108484343832, 0.016064689191803278, 1.9112451005309923, 0.015952227352576853, 1.868382146090487, 0.01571202020248426, 1.816296169392256, 0.015400082900094772, 1.75546881617413, 0.01506351450211631, 1.6865160151152072, 0.014709140497254108, 1.6102013509772244, 0.014276592546950958, 1.527276869755151, 0.013801536816473146, 1.4385299440449002, 0.013375191959113835, 1.3450479509622855, 0.012809852796920849, 1.2476875443574615, 0.012293385145222947, 1.1475623484036606, 0.0117542111836926, 1.0457835253568446, 0.011130608066089174, 1.8883859710138202, 0.03903554242071041, 1.8460655568964377, 0.038449665912163894, 1.7947229926609365, 0.03769076363435625, 1.7347961147571078, 0.03686736581103955, 1.6669083012008195, 0.03598558656343298, 1.5918615784189418, 0.03494174285526688, 1.5103788644953444, 0.03379685419948777, 1.4232037735198046, 0.03267984802040564, 1.3314060993996932, 0.03138189593980421, 1.2357949798593668, 0.030064436772055416, 1.1373440732538185, 0.028743826011855978, 1.0371686898000247, 0.02727066926384465, 1.8144639041854296, 0.07095497243484161, 1.7642125300032265, 0.06957060583664448, 1.7055896390711442, 0.06804478738888574, 1.6392371301903053, 0.06639138139500997, 1.5659770915035123, 0.06449373133058205, 1.486515333207938, 0.062409250962934884, 1.4015724740833921, 0.0602331851099063, 1.312018996747753, 0.057938047271527494, 1.218780168469913, 0.05547970581417395, 1.1227031232607763, 0.05298877723827064, 1.0248083532958387, 0.05037967810873415, 1.7255000953243755, 0.1107225740569898, 1.6686464416886104, 0.1082767244376233, 1.6042507042357044, 0.10561255899450277, 1.5331252894468157, 0.10262575711688447, 1.4560431183488896, 0.09932603468757631, 1.3737769787995038, 0.09577437755506915, 1.2870437791392755, 0.09213150333446857, 1.196762871621383, 0.08829342666262167, 1.103782754756272, 0.08423775522107571, 1.0088369859193111, 0.08018979865205669, 1.6238887966471196, 0.15712975875709095, 1.5618879903563638, 0.15322914400299378, 1.493389879726893, 0.1489216369640182, 1.4192429483323263, 0.1441299103005912, 1.3401506376392378, 0.13899120459496292, 1.2568843652669903, 0.13358724314006654, 1.1701107166428364, 0.12812811964339885, 1.0808263258541992, 0.12225375802914731, 0.9894906939226734, 0.11635786909586944, 1.5127119864251921, 0.20875300657220552, 1.4473512805223525, 0.20288530918346193, 1.376627134192269, 0.19634416322182102, 1.3010959853326935, 0.18944658785808854, 1.2217666365429574, 0.1819692885370762, 1.1390319854743225, 0.174481770763503, 1.0538582573453714, 0.16668585682622117, 0.9668634526346251, 0.1585375794716483, 1.3950911970527275, 0.2638752118580303, 1.3282064954028392, 0.25538858986117613, 1.2567117297663144, 0.2465080742224714, 1.1817763911882488, 0.23686357722001727, 1.1038492616130018, 0.22690645287348765, 1.0233582276855353, 0.21698015046464392, 0.9412040753792152, 0.2063769356596977, 1.2750299226738913, 0.3207037263174321, 1.2082178364833762, 0.3095870822787559, 1.1379754946487213, 0.29772535779817777, 1.0651434711328172, 0.28508212544485656, 0.9898544455716907, 0.27250196928199477, 0.9128458718798084, 0.259366164064485, 1.1552757594747536, 0.37787091032029374, 1.0902787263026297, 0.3635958780573183, 1.0229078575681516, 0.34837057998302884, 0.9534256259136753, 0.3326846767843517, 0.8820477969062542, 0.3168004047203417, 1.0394452799885727, 0.4336028306420557, 0.9779728614093465, 0.4158184648607399, 0.9147353246738655, 0.3970190230008276, 0.8492917549593114, 0.37805769538701184, 0.9307519897469217, 0.48637234684544983, 0.8740735074628183, 0.4646029649649481, 0.8148901461141091, 0.4426051291599224, 0.8312757929333625, 0.5344701227144323, 0.7787835952581795, 0.5097695178788924, 0.7414664803366856, 0.5788421788986758])
    SymCubatures.setweights!(cub, T[6.177158797170404e-5, 0.0003315809797526677, 0.0008012778158299902, 0.0014312829150622776, 0.0021580327660120466, 0.002929662147099171, 0.003672040558795897, 0.004336232336069512, 0.004840453868789263, 0.0051852291044034065, 0.0053005400510222, 0.005173415186053726, 0.00487403056163657, 0.004438507618387641, 0.00380836510381301, 0.00014317148896608053, 0.0002229174239064678, 0.00029882058730817723, 0.0003677561652666971, 0.00042889298084416007, 0.000482269185315206, 0.0005275264272004082, 0.0005614366083679786, 0.0005847707883186932, 0.0006023203816074067, 0.0006032228074322268, 0.0005997169042476812, 0.0005862861631281083, 0.0005610015175160643, 0.000515782864513203, 0.000690701362373195, 0.0008494800761628717, 0.0009900864638738387, 0.001112601864713123, 0.0012156099237892626, 0.001293859942908588, 0.0013482705160292802, 0.0013842076125868653, 0.0013915634806151375, 0.0013786341835621957, 0.0013478569933127677, 0.0012922782693960432, 0.0010715017015689974, 0.0013163946999377426, 0.0015328168510222225, 0.001720286286959417, 0.0018762131016221634, 0.001996299786293177, 0.0020798800544340615, 0.002126151053716339, 0.002144146132289082, 0.002121123673629057, 0.0020726929533755136, 0.0019958140919906836, 0.0017573536169275166, 0.0020453355942424193, 0.002291702535328548, 0.0024943349889277, 0.0026516836625071934, 0.0027608758622729183, 0.0028163581033951576, 0.0028395266934563425, 0.0028137643963950317, 0.002742322510985319, 0.00264887860746843, 0.0025133205613288943, 0.0028147816410980945, 0.003062167480805074, 0.003253887315628533, 0.0033827573871029163, 0.0034488930113187534, 0.0034615366001255485, 0.0034367493211422935, 0.003348317836053518, 0.0032316871808101655, 0.003279999508372013, 0.003565416990402673, 0.0037840344039180383, 0.00392744901199445, 0.004007589690788752, 0.004007152383834072, 0.003969821069906625, 0.003882099352343266, 0.00373721391350801, 0.003991089593592098, 0.004231563761706433, 0.004389562903615381, 0.004479027806306262, 0.004479696216662471, 0.0044108048945873035, 0.0043203839743461555, 0.004156496347498159, 0.004589242516057382, 0.004760889691811913, 0.004845067748924773, 0.004851572762652878, 0.004758325476972264, 0.004642611981319742, 0.004483717280201757, 0.00501848697358558, 0.005102121041628876, 0.005114704432633533, 0.005025836024489198, 0.004862564347077076, 0.004705597366268557, 0.005251677980762575, 0.005240256147444477, 0.005167524926987906, 0.00498825915477286, 0.00480906503930155, 0.005250742053008144, 0.005156146480601669, 0.005006274981778688, 0.004830731989377443, 0.0050414552472499265, 0.004894601194077141, 0.004776519238729736, 0.004706440229323193, 0.004600451643188559, 0.004311237981953578])
    cub_degree = 58
    tol = 5e-14
  elseif q <=60 # 721 nodes
    SymCubatures.setparams!(cub, T[0.0028745808779343796, 0.015110666869097948, 0.03689272996126802, 0.06765313845532146, 0.10622881767423487, 0.15166452536159966, 0.2025199452941509, 0.2570466017146392, 0.31439943976188345, 0.372357181078354, 0.4301515326796739, 0.4860097127398859, 0.5389149108123699, 0.5878165588771554, 0.6311823736384986, 0.9990369014151242, 0.9949206344001694, 0.9875166170217733, 0.9769357962311246, 0.9633539033444174, 0.9468335791664807, 0.9274268261048211, 0.9054639074691362, 0.8813004673837582, 0.8550279314116602, 0.8269147834826351, 0.7971471697312235, 0.765362058158803, 0.734956948496075, 0.7049618397924876, 1.9820082849459126, 0.0028733017523614743, 1.9600993855221205, 0.0028669086246095638, 1.928721449827852, 0.0028451661261097356, 1.8881460405285633, 0.002813617359323861, 1.8387498065763597, 0.002761593845312197, 1.7809957616878367, 0.0026856628074009136, 1.7153737056298073, 0.0026605397992374553, 1.6426571631679157, 0.002566481512731806, 1.5635238633958826, 0.002475348367547734, 1.478746266550299, 0.0024113235046818323, 1.3891672945970706, 0.002325096565205157, 1.2955542695575646, 0.0022323252066821972, 1.1987010998142753, 0.0021402474205636218, 1.0994959038665375, 0.0020599400256310885, 1.9479194370115827, 0.015071583367710447, 1.916670669755647, 0.014960642183202593, 1.876319481996453, 0.014788786785481968, 1.8272950433297142, 0.014515451058743793, 1.7700690572370872, 0.014137413791571777, 1.7048890703474244, 0.013960811608161101, 1.6327942053595437, 0.013506412049304053, 1.5541991244364133, 0.01301889941825777, 1.4696949246962565, 0.012672884186239798, 1.3803743402288766, 0.012208453065084123, 1.2871604839441542, 0.011727980421271112, 1.1911471808302627, 0.01124192939708486, 1.0933906069935997, 0.010782060522364795, 1.89525486843758, 0.03663361334252306, 1.8553895220415253, 0.03618852422049795, 1.8070288093656128, 0.03552718236635016, 1.7506082313099531, 0.03467063673744378, 1.6862447752224703, 0.03407995331773212, 1.615205449592109, 0.03310038785786708, 1.538056035026281, 0.0319069633664815, 1.455167749830453, 0.031013054184164505, 1.36778884972873, 0.02987731027576984, 1.2764952246289398, 0.028714674828887274, 1.1820528566459234, 0.02753468370066252, 1.085392991246467, 0.02626432551832737, 1.8255210949663134, 0.06677386340521381, 1.7780501007297755, 0.06558758671816925, 1.7227796249781435, 0.06412000106276984, 1.6599586073366734, 0.06274723883779322, 1.5905266172699941, 0.061149160817850795, 1.5153611148462196, 0.05903009460589017, 1.4344200510543461, 0.0572195040947774, 1.3489188293127798, 0.055152029301897165, 1.2595653688533792, 0.05298001246293531, 1.1671680374074815, 0.050815464598621504, 1.0727981354925, 0.04832684179340293, 1.7414182158951612, 0.10443293685172313, 1.6877026412478193, 0.10221521332170425, 1.626804377953394, 0.0997025268915262, 1.558979522409239, 0.09729425488186862, 1.4856737609012876, 0.09416102027427507, 1.4070394577342624, 0.09101356962829248, 1.3239670547852316, 0.0877537356250526, 1.2373030434482715, 0.08429324209587254, 1.1477321481315934, 0.08077577489880401, 1.0562567756732713, 0.07692045933787679, 1.6446780934174088, 0.1485047108795956, 1.5860398707834873, 0.14463568802355517, 1.5205136723248114, 0.14107092429021145, 1.4496528491801584, 0.13689057513002176, 1.3741081710171155, 0.13216769401855455, 1.2942661075830484, 0.12735111726439094, 1.2108948859659694, 0.1224382512326314, 1.124861021392325, 0.1171549331540784, 1.0365407218555465, 0.11188250934485526, 1.5390384373703285, 0.19724497122237744, 1.4765293195757816, 0.19211446365203397, 1.4085522299523352, 0.18665235412661169, 1.3361986246990512, 0.18039972691341336, 1.2599774020208423, 0.1736029309095279, 1.180084653272059, 0.16698125240453843, 1.097807787728888, 0.1597653670276651, 1.0133071044092647, 0.15267982862477064, 1.4266258013462152, 0.25003301402044903, 1.3619235306280517, 0.24271211842011006, 1.2926862883835653, 0.23501215519961743, 1.220358674356329, 0.2261347897179471, 1.1447674004906054, 0.21728799827570375, 1.066729885390322, 0.20826451223492917, 0.986849543175085, 0.19884729089166467, 1.3107353910836768, 0.30467484701390424, 1.2452133195671244, 0.29519809350169374, 1.1767274510464407, 0.2845207193831154, 1.105838227198624, 0.27300833040146094, 1.0324735894524222, 0.26181356240089204, 0.9575975657160386, 0.25006096467280586, 1.1945814642075094, 0.360343501308235, 1.130082423649742, 0.34780369246963566, 1.0635073492877594, 0.33391146809271743, 0.9952775506673277, 0.3195563865136498, 0.9256023698525188, 0.30538116835319074, 1.0806452392222416, 0.4149573704701197, 1.0185795632585921, 0.3990949860255431, 0.9555752576892778, 0.38167036445423935, 0.8917229046801584, 0.3642358779725977, 0.9717086226552848, 0.4673308788141996, 0.9135367303968588, 0.44754742500591327, 0.8550138076220749, 0.4275403961577622, 0.870352815857875, 0.515974191902912, 0.8169578749372197, 0.49409922283877544, 0.7787608996189476, 0.5603004318201514])
    SymCubatures.setweights!(cub, T[5.437518271745563e-5, 0.00029214649416493425, 0.0007056323974798909, 0.0012667060467892821, 0.0019121752612240642, 0.0026216689558157674, 0.003306355118082787, 0.003943293006545505, 0.004455225309451305, 0.0047978089742161745, 0.004972956809119636, 0.004998505929905516, 0.004886892477058277, 0.004679779767877706, 0.004229003297105809, 0.000498318609780815, 0.001140841045087574, 0.0017820206728828606, 0.002356186202655275, 0.0028621657310992556, 0.003345692100144343, 0.0037498454606397346, 0.004019982129536121, 0.004218206861687554, 0.0043328577452535506, 0.004331494145259286, 0.004319102290889414, 0.004167892291616759, 0.003431148351745416, 0.0037047076006364497, 0.00012605627535304873, 0.00019631206144388175, 0.0002631472979346438, 0.0003254353914588002, 0.0003805172658224127, 0.0004261298698099951, 0.0004729820251950502, 0.0005013575837221518, 0.0005216810936163518, 0.0005403964019995291, 0.0005473284906962487, 0.0005464878156712364, 0.0005393114062217612, 0.0005282920297944369, 0.0004544884227126598, 0.0006090787504982672, 0.000751926130278317, 0.0008783219507720991, 0.0009855467255001362, 0.0010872900432907053, 0.0011601069308185507, 0.0012085666573704737, 0.0012521020903499419, 0.0012664414435305223, 0.0012615322811495113, 0.0012374530189158373, 0.0011960397667422462, 0.0009452405586725311, 0.0011646467326509428, 0.0013608892717797193, 0.0015313638693284825, 0.0016751416302280585, 0.0017945426512001044, 0.0018663312935833984, 0.0019204764053035006, 0.0019422328505641513, 0.0019383502317506704, 0.001912973928463594, 0.0018444552888541893, 0.0015575846862185914, 0.001819343353904524, 0.0020457707985623405, 0.0022200772151432426, 0.002381282403429984, 0.0024889549179837674, 0.002551106021898365, 0.002585748820277953, 0.002575655800393419, 0.0025350563177876273, 0.002447939140150797, 0.0022374138570235265, 0.0025192026577190297, 0.0027345191438506522, 0.002927320795266417, 0.0030729418600355836, 0.0031347890629591624, 0.0031655449816359346, 0.0031540711567415686, 0.0030879943608021674, 0.0030074323410596224, 0.002947329304093401, 0.003201225725336982, 0.0034024543295984367, 0.0035710441387726086, 0.0036524010370957726, 0.0036679431249079346, 0.003659956015200426, 0.0035813220170760308, 0.0035020435741051643, 0.0036010934747488952, 0.003813794264721051, 0.003983427242115544, 0.00410299641864944, 0.004114862992024931, 0.004088039099437421, 0.004020097945231531, 0.0038934431920211684, 0.0041817552282982185, 0.0043309190931078104, 0.004458677563569905, 0.004487563668641186, 0.004411486793178473, 0.004348604283489666, 0.0042028023961104026, 0.004608039349742035, 0.004706265235718013, 0.004764035041323805, 0.004679268496995411, 0.004540186972018532, 0.00441516070885065, 0.00488285385078283, 0.0049277554426080845, 0.0048939471110029335, 0.004679073705447188, 0.004478319300482252, 0.004984692369992116, 0.004959942994000132, 0.004797347976035749, 0.0045270042389582614, 0.00493494641033973, 0.00479323899854941, 0.0045847320977177065, 0.00469784077179623, 0.004361476180049018, 0.004098699097084524, 0.004146043819077433])
    cub_degree = 60
    tol = 5e-14
  elseif q <=62 # 768 nodes
    SymCubatures.setparams!(cub, T[0.002697396635832736, 0.01420151968888748, 0.03472447082964377, 0.06369961937801394, 0.10019228118036184, 0.14310136819777314, 0.19137109251396517, 0.2437763877935655, 0.29892543421591355, 0.35534614854685076, 0.4117057406396788, 0.4666430983622796, 0.5187655311806209, 0.5672323812810857, 0.610742435642228, 0.6485288281962954, 1.9831108081706552, 0.002700260582601698, 1.962536572138339, 0.0026967024520938754, 1.9330522732559172, 0.002679217625404378, 1.894901584895951, 0.0026474990996807985, 1.8484077478511174, 0.002600788048167338, 1.7939661441768002, 0.002545667576367977, 1.7320456035485323, 0.0024944036305862356, 1.6632156769374618, 0.002434062108783147, 1.5881079103818887, 0.0023600694942550408, 1.507399873571235, 0.002287221524446419, 1.4218401008874646, 0.0022140549290321086, 1.3322308708533244, 0.002131305961227443, 1.2394043089580367, 0.002044383202523729, 1.1442379959975217, 0.0019508712866855371, 1.0476365201535907, 0.0018616124093478372, 1.9510326220028011, 0.014180219798576248, 1.92163228538315, 0.014087089323026474, 1.8836583697591853, 0.01391916018497956, 1.8374496817914274, 0.013674251285838957, 1.783383872040172, 0.013386398422779687, 1.7218896618888542, 0.013113028595460144, 1.653577585000304, 0.012795292684696698, 1.5790852478457071, 0.01241097326669597, 1.4990356105527447, 0.012024467633574832, 1.4141722527206007, 0.011637869010375442, 1.325331249862324, 0.011203485504388358, 1.2333302870681329, 0.010748108786563675, 1.1390423300852384, 0.010265389761661867, 1.0433282707968883, 0.009795614088254743, 1.9013685666751974, 0.03449192661579387, 1.863772782181998, 0.034076926282668965, 1.8181136787017682, 0.03348143384435835, 1.7647396185511288, 0.032783921098955754, 1.7040375686534706, 0.03209942591912035, 1.636672249396783, 0.031320128400460286, 1.56330508932746, 0.03039689370049244, 1.4845083074532708, 0.029440843352357377, 1.4009780402919307, 0.02848499633833049, 1.3135549960262198, 0.027425935261750533, 1.2229882501331886, 0.0263183564091342, 1.1301063216118346, 0.02516342171907458, 1.0357425535338716, 0.02401110552924274, 1.8356440999915686, 0.06292686647229191, 1.79084089529809, 0.0618408744954517, 1.7385072078862955, 0.06056852164767347, 1.6790025066116798, 0.0592713784463205, 1.613009588645091, 0.05782894445963163, 1.541212251712954, 0.05615891506096333, 1.464154304293288, 0.0543825981988139, 1.3824586406854344, 0.05258734369289331, 1.2969557720980058, 0.050640566173336, 1.2083429567882158, 0.048614016701646466, 1.117394783783008, 0.046519133030969405, 1.024902877763642, 0.04439175044608038, 1.7559885459104965, 0.09849391470743296, 1.7050411527211042, 0.09649409705258619, 1.647108959414263, 0.0943771860195045, 1.5828352982516893, 0.09206818817374887, 1.51294179789244, 0.08945217148770178, 1.4379971791445212, 0.08662774971584322, 1.3585923006351515, 0.08370925282862188, 1.275490485355251, 0.08061599882906821, 1.1893442194211807, 0.07742554301170611, 1.1008962541813903, 0.07411686640925975, 1.0108625062263468, 0.07073703163040224, 1.6645288284030908, 0.14023197602803592, 1.6085071195201062, 0.13709664376628633, 1.546292190129501, 0.13371137310091788, 1.47863152618356, 0.1299446880229945, 1.4061315278104407, 0.12587486882415944, 1.3294443752057996, 0.12156240766384532, 1.2492309332434324, 0.11705868706112017, 1.1660888771728262, 0.1124719463211995, 1.0807446942213152, 0.10767495786988868, 0.9937601064835203, 0.10277226069679703, 1.563538871011885, 0.1870497987305218, 1.5038417107317916, 0.18238132936575155, 1.4388863280445847, 0.17725745145588817, 1.369262385201003, 0.17176656840904686, 1.2957106335666724, 0.16584656452830498, 1.218791844446638, 0.1596670028856184, 1.1390457658026305, 0.15342019149089042, 1.0572315553973057, 0.14688464611275742, 0.9737865110436557, 0.14018914499070625, 1.4556409207905248, 0.23762025426987401, 1.3937418328760809, 0.2309234275944948, 1.327321479298815, 0.22381169111606927, 1.257225880817797, 0.2161344037104936, 1.1840352532708198, 0.2080647573717594, 1.1082285582969837, 0.19986337991414477, 1.0304527823464285, 0.19138912952301101, 0.9510749057194783, 0.18267841532447313, 1.3438015375366872, 0.2904721067705455, 1.281095695890735, 0.28152736876467616, 1.2148200590468694, 0.27197009038937675, 1.1456518082648497, 0.2618700247173205, 1.074151288494937, 0.25142818920151044, 1.000767969089003, 0.24079493209484276, 0.9258328500884605, 0.2299073568673555, 1.2307158383994399, 0.34434578572218116, 1.168660959595166, 0.3327410170207518, 1.103858297509025, 0.32051367992305313, 1.0370046340904782, 0.3076875650902965, 0.9684210596406083, 0.2946457166145207, 0.8982992434510529, 0.28142842116605643, 1.1191656403873567, 0.39783244904581205, 1.0591065741095598, 0.3833050115813323, 0.997134718744069, 0.36807509335900174, 0.9336739957480109, 0.35242718073048773, 0.8687558994555266, 0.3366582858414375, 1.0118342965972267, 0.44957042986732637, 0.9550780676518974, 0.43189531775675577, 0.8970314992141397, 0.4135513355998484, 0.83756719698514, 0.39501528872835345, 0.911192583002531, 0.49852255314033816, 0.8588697338306434, 0.4773953676735839, 0.8050485873794785, 0.4560382001551164, 0.8188783961852303, 0.5433577104343509, 0.7711142385491463, 0.5193311201279591, 0.7360881382673587, 0.5844088223351406])
    SymCubatures.setweights!(cub, T[4.788245861448093e-5, 0.0002580929907518513, 0.000625919191980983, 0.0011233752430201643, 0.001708548312401227, 0.002340311369490357, 0.002972913207685822, 0.0035579464015742302, 0.00404918525482742, 0.004419880040391597, 0.0046368540647748506, 0.0046879639776064115, 0.004561443190039606, 0.004301990543276285, 0.003876459682722278, 0.003298490057232257, 0.0001112054763491203, 0.0001734442272310742, 0.00023289931044447891, 0.0002880806070797087, 0.0003375038743237744, 0.0003808263895951607, 0.00041924703024859674, 0.0004503927740048761, 0.0004727795077770118, 0.0004888422785958621, 0.0004985070745814413, 0.0004997864359201552, 0.000494013549951542, 0.0004809701432332615, 0.0004635397798145613, 0.0004021636375470252, 0.0005395300265705898, 0.000666700445732782, 0.0007805098719943196, 0.0008802630089656843, 0.0009678838958109841, 0.0010391847655192138, 0.0010912223199361052, 0.001127361875582075, 0.001148818746327791, 0.0011513606993358574, 0.0011377690038691204, 0.0011089743329953637, 0.0010685677251002512, 0.0008388502418564221, 0.0010354295498634595, 0.0012113121060638062, 0.00136522980711275, 0.0014982067463125278, 0.0016067015954538803, 0.0016871426112666555, 0.0017408198459055667, 0.0017716333281500088, 0.0017757861000428003, 0.0017560828533142296, 0.0017149192085184719, 0.0016531180599979125, 0.001385617973381561, 0.0016204493456863601, 0.00182528448451648, 0.001999022509928268, 0.0021410550732368742, 0.0022479994260898696, 0.0023180424450896736, 0.002355063389532459, 0.0023603887423969527, 0.002335935767788881, 0.0022831811162037415, 0.0022028340367808708, 0.0019987806246195695, 0.0022516441716978843, 0.0024632522477444495, 0.002635521403581341, 0.0027655882714750423, 0.002850710978781757, 0.0028902181812608735, 0.0028938310549564, 0.002864595139806666, 0.002798563770988948, 0.0027016844754617236, 0.002637597159478765, 0.002884507103914168, 0.003083146372702512, 0.0032319214338700974, 0.00332978651988988, 0.0033697563210275294, 0.003366133398326751, 0.0033282709406546106, 0.0032496034869395687, 0.003138067046946008, 0.003251009611377, 0.0034709352325644712, 0.0036341565684547693, 0.003742519900335215, 0.0037876402763670687, 0.0037772816033502625, 0.003724401635835898, 0.0036339413575731646, 0.0035089146936048386, 0.0037973504442969375, 0.003972746400417764, 0.004085009601063672, 0.004133084058364022, 0.004116073290999455, 0.004042639232601965, 0.003940750014151773, 0.0038090445813911992, 0.004232978235778475, 0.004347227157961506, 0.00439856190186626, 0.004380698343528404, 0.004290814814535367, 0.004170145815089036, 0.004034134563184176, 0.004530926906672339, 0.004576344073131938, 0.004554855805991882, 0.004460509427061219, 0.004318428691380095, 0.004173364626459425, 0.004671010453481649, 0.004636411914705439, 0.004541417704861058, 0.004384146546565608, 0.004221491976106044, 0.004633310710475852, 0.0045280102183208985, 0.004365728698824588, 0.004196701005982506, 0.004444003674779544, 0.004274091463867159, 0.004120625456250823, 0.004124665175728683, 0.0039694426217756085, 0.0037265235699705797])
    cub_degree = 62
    tol = 5e-14
  elseif q <=64 # 817 nodes
    SymCubatures.setparams!(cub, T[0.0025486348235609618, 0.013359596531524583, 0.03265340206989669, 0.060109169046802594, 0.09474494459989659, 0.13530964930943504, 0.18115436328767565, 0.2312379266723597, 0.28417052497833734, 0.3388749677024106, 0.3936467649338445, 0.4473731852595884, 0.4991597850612063, 0.5481614920724025, 0.5938764397392041, 0.6331380305154178, 0.9991506967117468, 0.9955363362686358, 0.9890166720444393, 0.97954301484334, 0.9673216040402502, 0.9527431489602736, 0.9357548150730145, 0.9160084049424106, 0.8939120528991422, 0.8703135430188786, 0.8451801420698655, 0.8182301302077843, 0.7897386334523115, 0.7599730714735069, 0.7346716270621299, 0.7063747071892973, 1.9840541177883555, 0.0025398493710011086, 1.9646188345236588, 0.002534275284089833, 1.936756322580118, 0.002527004287908815, 1.9007109269223015, 0.0025022636059210397, 1.8567873947531255, 0.0024538107008139254, 1.8053367763171957, 0.0023989962059846118, 1.7467787236210224, 0.0023537287518937, 1.6816148052723376, 0.002312304491356434, 1.6104135107945665, 0.002254626040948167, 1.5337835985701698, 0.002168497452195527, 1.4523398931917295, 0.0020866610841242497, 1.3667319427083593, 0.0020647236572672823, 1.2778245628674896, 0.0019556984567288436, 1.1863161624304097, 0.0019099282484256737, 1.0931578238755641, 0.0017902704717231888, 1.9539383439933882, 0.013329271501723512, 1.9262070214470923, 0.013288531924756361, 1.890391288436227, 0.013157893745505375, 1.8468146597397832, 0.012905775525333273, 1.7957501819744324, 0.01262112872280354, 1.7375552057957597, 0.012377425144074456, 1.6727406887973642, 0.012156952469809412, 1.601963284903815, 0.011851646356004643, 1.5258863067697357, 0.011408080734074098, 1.4450427400104429, 0.010992591846762203, 1.3600033387031818, 0.010811193256353484, 1.2719424695507149, 0.010315718129570944, 1.1812485573065017, 0.009973433481187767, 1.088877136147486, 0.009466522271506534, 1.9071128057582019, 0.03254199813718169, 1.8716139934789373, 0.03221696896605575, 1.8285680542559135, 0.031609500538401884, 1.7781410791787358, 0.030926596451419194, 1.7206360655727333, 0.030310443065400976, 1.6565671930199966, 0.02975763787153538, 1.586738448654188, 0.029003181551225877, 1.5118093150621235, 0.027960593849879624, 1.4321065531143993, 0.026984376836842842, 1.3481409774013273, 0.026327357937857697, 1.261183741534412, 0.025318822159179615, 1.1717974428764684, 0.024247446023686836, 1.0807566773746329, 0.023239612643876768, 1.8447684441992689, 0.059496310168548915, 1.802500577128123, 0.05840195486834098, 1.7530323862098602, 0.057171952622632784, 1.6966400105937511, 0.05599584803203181, 1.6337799573378937, 0.05493687780802382, 1.5654213115776787, 0.05353184333403286, 1.4921777049356826, 0.05172069945740167, 1.4142817601155786, 0.04996237453005197, 1.3322856474936784, 0.048434546886496135, 1.2471454701563434, 0.04675257713486458, 1.1596513202336527, 0.044689796900745045, 1.070140138813396, 0.04280294199150009, 1.7692875569981326, 0.09306205907439989, 1.7210415719769057, 0.09115185275131588, 1.666042961559328, 0.0892210597161203, 1.6046654498108521, 0.08744531070271701, 1.5380516018676385, 0.0851828740887211, 1.4667008759861724, 0.08249645413026589, 1.390909813867318, 0.07970573485464164, 1.311123005707742, 0.0770120764708597, 1.2281660705093222, 0.07422451972330518, 1.1428216169920866, 0.07123354153966523, 1.0558499988189525, 0.06791197505774263, 1.6826383409952406, 0.1326021995281085, 1.6293122224283039, 0.1297313049362935, 1.569657267249802, 0.12698517881717591, 1.504881375848768, 0.12364826174005923, 1.4353888757923545, 0.11997722117226163, 1.3619064442464988, 0.11591708042648706, 1.2847911978625162, 0.11187932549871156, 1.2048746152064662, 0.10757084398435252, 1.1223246068608206, 0.10352489105195425, 1.0383121806249138, 0.09877185381226236, 1.586445438428102, 0.17720275434517596, 1.529142922910701, 0.17322670808506674, 1.4668809998012409, 0.1686413082007314, 1.3999206298625764, 0.16381660778923654, 1.3291402254476816, 0.1584029916534715, 1.254898481569751, 0.15275678447073346, 1.1778937217233947, 0.14688322975982337, 1.0985872811509751, 0.1410961128617252, 1.0175713464107325, 0.13529744304952937, 1.4828863276545443, 0.2257506511159763, 1.4232150919017537, 0.21975531886409944, 1.3590475094277392, 0.21341137902651292, 1.2911581673925687, 0.20669093811241557, 1.220525678281409, 0.19917370215678998, 1.1469493905413186, 0.191873581116324, 1.0715579056516304, 0.18393241407565908, 0.9940834582015973, 0.1765953410647343, 1.375306418125975, 0.2767141893608741, 1.314750779629407, 0.26848318522747894, 1.249919167009498, 0.2603210528402829, 1.1826746965694142, 0.25090948984634837, 1.112771339045399, 0.24174196602718903, 1.0413645193691436, 0.23204131488143098, 0.9684729746970473, 0.22201501206456625, 1.2655915477409359, 0.3286057492603276, 1.204840545263626, 0.31844467051356123, 1.1412460447414965, 0.30750882762119636, 1.0755199999765932, 0.29576163171062053, 1.0079395153425759, 0.2842465260771629, 0.9399142270064105, 0.2717264763437451, 1.1563462710759422, 0.38066525935777246, 1.0965632712769664, 0.3679675699159608, 1.0354346289175227, 0.3538886183556785, 0.9728309363896182, 0.3393909963718818, 0.9090661125031959, 0.32521275631731233, 1.0503323523358912, 0.43168994182849457, 0.9930393285784841, 0.415697931474996, 0.9356750537295394, 0.3982717969146247, 0.8772316249477786, 0.381523618170783, 0.9494214199379338, 0.4801980878364396, 0.8953036241858717, 0.46106293894677647, 0.842130808633103, 0.4417010347333659, 0.8547147829247783, 0.5260304139238664, 0.8050548724091747, 0.5047609618467293, 0.7709771770891091, 0.5654558949292985])
    SymCubatures.setweights!(cub, T[4.2746031773129734e-5, 0.00022846923355252427, 0.0005541953136029664, 0.0010012487193992456, 0.0015297251178330328, 0.00210193523836387, 0.0026851980766191184, 0.003239754654605982, 0.0036908849836180865, 0.00408377104303437, 0.004373891287527147, 0.004487833743821573, 0.004467955076227746, 0.0043845275850856735, 0.004142021610247399, 0.0038988191018807807, 0.00041014282260510497, 0.0009485525936438797, 0.0014810265754261179, 0.0020032000027930537, 0.0024243139622239195, 0.0027685665154956727, 0.003148184491350159, 0.003497654548620421, 0.0036792293997566076, 0.003724184665918306, 0.003785102886570224, 0.003822526313768326, 0.003786721613933286, 0.0033816937395766238, 0.0024493117842982796, 0.0036188283313534583, 9.882712975287154e-5, 0.00015399390727938752, 0.0002075164738374713, 0.0002572284301567449, 0.00030090379393246517, 0.0003392890670048434, 0.00037426320325368846, 0.0004052315579705256, 0.00042841139740123504, 0.00044070262643834587, 0.00044816677575226945, 0.0004625742485575382, 0.0004539410757150592, 0.00045243482042161794, 0.00043146110645205784, 0.00035593897709998776, 0.00047947064001884604, 0.0005943026973975635, 0.0006955784333136074, 0.0007849095535398054, 0.0008652800029745916, 0.0009365612635720001, 0.0009896686198039297, 0.0010191873727491542, 0.0010379033651750594, 0.0010578722524522753, 0.0010507528698766437, 0.0010295138552908597, 0.0010034153453214137, 0.0007456492100170128, 0.000923341856530544, 0.0010805946252933766, 0.0012194255855426256, 0.0013416083787372271, 0.0014492580461559932, 0.0015295011790198967, 0.001579442719257713, 0.001611664617268094, 0.0016233050636772282, 0.0016265467452826973, 0.001578659312062378, 0.0015420399497631091, 0.001237880893685043, 0.0014479944998441487, 0.0016328489746885214, 0.0017919081130175483, 0.0019301471290953246, 0.002034070091751288, 0.00210695756939881, 0.0021479918513891616, 0.002151631746726913, 0.002149327521863987, 0.0021136940844336425, 0.0020381615558115822, 0.001791079028128671, 0.002020363595004782, 0.0022139428774091167, 0.0023769698267616032, 0.0025004114974735373, 0.0025966674774463732, 0.002644725059975149, 0.0026539547924452934, 0.0026234838929252834, 0.002604590275713439, 0.002502053976858177, 0.0023749235502583117, 0.002604342246591618, 0.002791219513681267, 0.0029343614098884985, 0.0030458120385339306, 0.0030954556041956453, 0.0030994042680063277, 0.0030590195755991916, 0.003013699163561822, 0.0029551836654889323, 0.002945505383951965, 0.003147068431259058, 0.0033060241853187524, 0.003424076174631585, 0.003498437219933154, 0.0034961758386128257, 0.0034758558772307966, 0.003369346161995462, 0.003324706664557767, 0.0034603223967232693, 0.003639102494378333, 0.003743289874522367, 0.00382454942358002, 0.0038137683088956967, 0.0037971480508349176, 0.003704410500922186, 0.0035694349970547727, 0.003889135985522046, 0.004004301551672431, 0.004085617741772835, 0.004098844162086131, 0.004014115120925449, 0.003955073557936538, 0.003774033384589026, 0.004205150319720624, 0.00424050914755026, 0.004295894760715096, 0.004196823637190979, 0.004056977874624313, 0.003933054877538728, 0.004390254918979633, 0.004370383815344808, 0.004318967064572085, 0.004107407381779082, 0.003961641375366591, 0.004437233437102001, 0.004356982798952508, 0.004200005810748641, 0.003966105889351012, 0.004383667926864817, 0.004247554258168967, 0.004002101028667051, 0.004192537995457968, 0.003748901202273806, 0.003465437733495755, 0.0041020352335533255])
    cub_degree = 64
    tol = 5e-14
  elseif q <=66 # 867 nodes
    SymCubatures.setparams!(cub, T[0.0023950506702950956, 0.012615316872951494, 0.03087491288494519, 0.0567353985405659, 0.08944260974038083, 0.12812835677260426, 0.1719114475134776, 0.21981939645829923, 0.2706752282205768, 0.3233809437336214, 0.3766974727321983, 0.42948951365412863, 0.4806156067971628, 0.5289427637501485, 0.5739193300041471, 0.6143528876699889, 0.649649207941718, 1.9850010001816947, 0.002398171027517095, 1.966718178282592, 0.002396126403313788, 1.9404948139057883, 0.002383286322921836, 1.9065243309574575, 0.002358440487914191, 1.8650608160922622, 0.002322169638952118, 1.816416991673326, 0.0022790092616262157, 1.7609668065095367, 0.0022343499878052673, 1.6991493851458315, 0.0021887240317492117, 1.6314702595712245, 0.002132873379482222, 1.558476604897976, 0.002068732257162426, 1.4807508503694065, 0.0020124998096735375, 1.3989478411263758, 0.0019469823934233116, 1.313742418442944, 0.0018734999839298315, 1.2258378460084702, 0.0018005268512749377, 1.1359756473437197, 0.001726960365406886, 1.0449237449792688, 0.0016490634854668647, 1.956484100925683, 0.012602755264150303, 1.9303129966629013, 0.012534081031657698, 1.8964672660520536, 0.012402565522713554, 1.855208598860221, 0.012212194779674418, 1.8068394157369851, 0.011986019016789218, 1.7517179089296662, 0.011750254571996209, 1.6902821360877789, 0.011508038818812167, 1.6230671663892844, 0.011216210529724956, 1.5506105736657108, 0.01088082861847598, 1.473437129914622, 0.010580391836635552, 1.392254702192738, 0.010235561392802517, 1.3077298268377502, 0.009854338726576452, 1.220542372965106, 0.009470829157172662, 1.131422920134406, 0.009084202623773104, 1.0411223026580976, 0.008681256820513487, 1.9122337665593256, 0.030702400200946424, 1.8786698395452763, 0.030377526714972063, 1.8378239946116706, 0.029913715061082945, 1.789982780203278, 0.029363145583171248, 1.7354862991164466, 0.028782004105296247, 1.6747740445289774, 0.028179948825451707, 1.608419982550153, 0.027472770858188344, 1.536960374001866, 0.02666028017264365, 1.4608372990407912, 0.025905678110730117, 1.3807950985213608, 0.025060912556012085, 1.2974617500955006, 0.024146402770782133, 1.211480811297775, 0.023207742168739946, 1.123547607889453, 0.02226269654958592, 1.0344000654725856, 0.021292994991706293, 1.853458145599384, 0.056130582330081134, 1.8132743057550593, 0.055282164194071934, 1.7662444451766373, 0.05427347809694873, 1.7126919519202441, 0.053191428506676304, 1.6530539716439625, 0.05205910492663197, 1.5879279820556995, 0.050767485752773296, 1.5178635462581362, 0.04929024796273889, 1.4432328781352357, 0.04785340794246141, 1.3647683847456515, 0.046290994630975864, 1.2830572204046835, 0.0446377352757302, 1.1987296005076862, 0.04290741969072112, 1.1124098471087243, 0.041167272951985405, 1.0248098243858395, 0.0393932833024645, 1.781858922856735, 0.08811090300134701, 1.7359337642061579, 0.08652041709674144, 1.683636989283988, 0.0847834847547884, 1.6253870999052396, 0.08294434694146342, 1.561777933214297, 0.08090169428795572, 1.4933948506015957, 0.07859335677814512, 1.4206146697135835, 0.07623678553365469, 1.3440971493790193, 0.07373908389300407, 1.2644001930599642, 0.07114657824515604, 1.1821666701120659, 0.06840746481194801, 1.0979379956605122, 0.06564222201437868, 1.0123789771008844, 0.06282012393881023, 1.6991603066665966, 0.12584144442455245, 1.6483811209398582, 0.12330158881126498, 1.5918002650604302, 0.12057627975976207, 1.529977134324533, 0.11760858433418989, 1.4635272074870451, 0.11431769083641961, 1.3929267643663235, 0.11081861155986321, 1.318722307839605, 0.10717044450831417, 1.2414683650194525, 0.10342199849224966, 1.161795950620723, 0.0994869377508514, 1.0801849141537594, 0.09546679834892921, 0.9971969343542496, 0.09135597056429128, 1.607226835192126, 0.16843835938369373, 1.5526877838215176, 0.16465685126511576, 1.4930444242612655, 0.1605884876833616, 1.4288963653657576, 0.15616789312821394, 1.3608480803949454, 0.15134916859689318, 1.289310213567111, 0.14634833349723822, 1.214852323351767, 0.14119732896437914, 1.138048745878758, 0.1358848123554066, 1.059432820940922, 0.13038093976373863, 0.9794286900943654, 0.12475571404402412, 1.5081615871681464, 0.2148223656043538, 1.4509986404287372, 0.20947377184778002, 1.3894412217922547, 0.2037487038798301, 1.32421220025546, 0.1974785534975663, 1.2556809030238407, 0.19095363548698718, 1.1844548135157742, 0.18416396736240337, 1.1109678890779655, 0.177266153540192, 1.0357851752736302, 0.17009621458155097, 0.9591989883584328, 0.1627626067214885, 1.4043709036157233, 0.2638921105573205, 1.345831395334205, 0.2566728066691938, 1.2837287814472849, 0.24884634829937363, 1.2184706984268117, 0.2406560749684159, 1.1507278008587296, 0.23204499013158836, 1.0808844995068723, 0.22328646221136872, 1.0094383741495245, 0.21430797127946907, 0.9366445379592623, 0.20510672372844138, 1.2980721540170497, 0.3144699272519598, 1.2394325288231678, 0.3049648573524441, 1.177797245110718, 0.2949907805044457, 1.1138360702783496, 0.28448260664447983, 1.0480212368706292, 0.27361450627406747, 0.9806280787585318, 0.262668684170185, 0.9119502613670365, 0.251468998807903, 1.191777319442324, 0.3653532233664612, 1.1340909894493483, 0.35344674260491776, 1.0741290963312058, 0.34098431289308057, 1.0125596630064109, 0.32790492097032026, 0.9495421958249194, 0.3147579249092928, 0.885310451868797, 0.30142588713765694, 1.0876340774589441, 0.415485190714095, 1.0320781737377471, 0.4009322803683362, 0.9749742709764182, 0.3857221757869989, 0.9166231522756756, 0.3701556750732641, 0.8570561614945801, 0.35447655550995294, 0.9879291850966427, 0.46370533801422037, 0.9355014429746554, 0.446371213778456, 0.8821058959315604, 0.4283329091768226, 0.82745529129712, 0.41013738159662033, 0.8945675480378497, 0.5092965762154917, 0.8462505662457899, 0.488752506982747, 0.7966687889088951, 0.4680531176943008, 0.8088928129825161, 0.5509623026138608, 0.7646697935187546, 0.5279392124468004, 0.7317597038846532, 0.5893376372715906])
    SymCubatures.setweights!(cub, T[3.775316666362345e-5, 0.00020376198689250893, 0.0004955567077084097, 0.0008936028353755123, 0.0013678181166645797, 0.001888236904296481, 0.0024213795267460643, 0.0029328817002281737, 0.0033836723544466395, 0.003752984680384274, 0.004009905600115169, 0.00414826602349061, 0.004150033292589912, 0.004012890696678458, 0.003773159446754633, 0.0034064115183782048, 0.0029065778959600024, 8.773429527157661e-5, 0.0001369936900725632, 0.00018434433606265586, 0.0002286533420326912, 0.00026895687043007214, 0.00030492846135663523, 0.00033675958747919464, 0.0003642707891564383, 0.00038564876278698683, 0.0004007563098954073, 0.00041258218960326775, 0.00041787288413987954, 0.00041683761490315463, 0.0004113600020034517, 0.0004015237651508663, 0.00038686660722047814, 0.00031791725914616446, 0.00042745450004231893, 0.000529755370450601, 0.0006227319726334878, 0.0007056382374845851, 0.0007787108561038914, 0.0008415613611723939, 0.0008908448106017925, 0.0009257249127359936, 0.000951845769695471, 0.0009635994457330456, 0.000961741788264201, 0.0009487605507910909, 0.0009258932436572275, 0.0008932784640160183, 0.0006656932457126375, 0.0008242379181086459, 0.000968224504641952, 0.0010963419689263528, 0.0012084327985146968, 0.001303990963063199, 0.001379814393510266, 0.0014337462389446551, 0.0014713416749738552, 0.00148890351276539, 0.0014878128679254633, 0.001467961209820525, 0.0014331870992350562, 0.0013844308359119203, 0.001105734889072003, 0.0012984749518322894, 0.0014695124986206654, 0.0016178028670081246, 0.001742823063320065, 0.0018430874741711537, 0.001915463977096975, 0.0019617191125699237, 0.0019841215181034596, 0.0019844819997453756, 0.00195903240722407, 0.001913931019820999, 0.0018498396151819454, 0.0016065532214621862, 0.0018182409450078283, 0.00200061889509987, 0.0021527129143625225, 0.002274994267152182, 0.0023648449669325636, 0.0024175386677451013, 0.002442560309835886, 0.0024419747594548887, 0.0024119459983961825, 0.0023564132893321944, 0.002277362715126188, 0.0021380091333799336, 0.002352390180578431, 0.002529174960867808, 0.002670297781350702, 0.0027751846864333828, 0.00283321958019968, 0.0028575929404952008, 0.0028504199299704256, 0.002815691839721855, 0.002749215676674752, 0.002657459199431375, 0.002664290954335047, 0.0028624495361838074, 0.003019028126184896, 0.003135791207972484, 0.0032017068970504746, 0.0032263775350853566, 0.0032100321811923483, 0.0031674412204241594, 0.0030898523846973175, 0.0029873320192419952, 0.003150638754035521, 0.0033207326405665016, 0.0034444045795987608, 0.0035169935317599186, 0.003540225992090539, 0.003513906243782805, 0.0034569500600985454, 0.0033723141082691605, 0.0032631503025101817, 0.0035648626685567315, 0.0036933512577703416, 0.0037732269162692672, 0.0037977146098942074, 0.003766927235159724, 0.0036891970556591257, 0.0035950814984714286, 0.0034810198105193653, 0.0038821043562384874, 0.0039624772376791245, 0.003985692278240181, 0.003955864153724828, 0.0038643751208988467, 0.0037542904441563297, 0.0036357974597625974, 0.004085976671362502, 0.004105236144403418, 0.004073552376567322, 0.003981394282439093, 0.003847599894561807, 0.003718278946168524, 0.004154989770425244, 0.004108901834302496, 0.004024469678119942, 0.00388028468510121, 0.0037361718882328856, 0.00408463055664422, 0.003987603770709919, 0.003843629305850737, 0.003701342776372107, 0.0039040643174404603, 0.0037545124780060312, 0.0036245685659381352, 0.0036223069622885857, 0.0034904286058196577, 0.0032792336828855024])
    cub_degree = 66
    tol = 5e-14
  elseif q <=68 # 919 nodes
    SymCubatures.setparams!(cub, T[0.0022668166125661996, 0.01191991772299086, 0.029153187437044162, 0.0536385337917292, 0.0847857201887196, 0.12159575955654112, 0.16325171058433083, 0.20902126548893282, 0.25788766013949604, 0.3085762987496843, 0.36020627368864083, 0.4117379928673805, 0.4623301078397741, 0.5102494641867084, 0.5560967210825067, 0.5986259639001822, 0.6349767797173258, 0.9992123843999092, 0.9958989619083326, 0.9900925220456637, 0.9818502238199476, 0.9710441151677618, 0.9577268037829213, 0.9422792563177775, 0.9248944146661496, 0.9052525105203783, 0.8835216175767255, 0.860443292418528, 0.8360696632177981, 0.8102008201695918, 0.7832432258802142, 0.7550190585044889, 0.7309549461726995, 0.7039195067985586, 1.9858074619295378, 0.002265829527323669, 1.9685028063336136, 0.002261501285097067, 1.9436713209150558, 0.002252295915464058, 1.911494871240952, 0.002232717638347705, 1.8722045832755578, 0.002203727587351616, 1.8260936505724528, 0.002158299219555545, 1.7734905898541535, 0.002107663273945808, 1.7147484485639057, 0.00209410386259779, 1.6503994803802893, 0.002024658825250713, 1.5808741366340187, 0.0019683394822871663, 1.506725170920122, 0.0019116187130989136, 1.4285111423441694, 0.0018778829357261285, 1.3468915341506948, 0.001802881669290822, 1.2624660125946718, 0.0017363842508850648, 1.175906518471671, 0.0016813008199433952, 1.0879554032995233, 0.0015626398750608094, 1.9588902347318762, 0.011896703055331572, 1.9341321246940215, 0.011846025873119033, 1.9020905067991842, 0.011744319395356027, 1.863003498323119, 0.011589454827507398, 1.8171893124211194, 0.011353398917845816, 1.7649372094824725, 0.011093483575398997, 1.7064672033462647, 0.011000432371670527, 1.6426116835567863, 0.010659831356631445, 1.5735957239747413, 0.010349353801482974, 1.4999522137271666, 0.010067642289908953, 1.4222335654214018, 0.009845434786400839, 1.3412002259990947, 0.009496347291355846, 1.2574132611747315, 0.009111778370089663, 1.1714325513105348, 0.00883920390493401, 1.0841783191367025, 0.008275196730273736, 1.917078947772799, 0.029020666371239172, 1.8852794523924592, 0.028775258092989513, 1.846573009486462, 0.028387498889152592, 1.8012930175643989, 0.027821260862190507, 1.7496577083865725, 0.027207961737930345, 1.691739129230029, 0.026896553518032433, 1.6286680977319283, 0.02615051621056134, 1.560571226100227, 0.025351662761523095, 1.4878110257182477, 0.024699931912148823, 1.4111047855023608, 0.024023531470482487, 1.3310720975865578, 0.023277568193396145, 1.2484365222623819, 0.02229074639825456, 1.1634744842644371, 0.021596333806619208, 1.077313892509414, 0.020468933136196375, 1.8613081961374907, 0.05318966761370522, 1.8231804468470598, 0.052458577772934556, 1.7786478837207618, 0.05144124802174646, 1.727861759900845, 0.050353367081903286, 1.67085823898166, 0.049601504008168186, 1.6088245927249556, 0.04838887547543319, 1.5420606255484748, 0.04688493671996933, 1.4705584150296405, 0.04568797644334287, 1.3953143825091145, 0.04427106698852979, 1.3166609614170384, 0.04294737618311383, 1.2355562016370139, 0.041228368639873655, 1.1521697202905168, 0.039733322357438976, 1.0673769118170298, 0.03811727479077638, 1.79308553368788, 0.08360646545768141, 1.7495173121883776, 0.08203867537369595, 1.6998186651902099, 0.08036539924840587, 1.6441158366044777, 0.07890828378132544, 1.5833372070700542, 0.07716755591422783, 1.5182281565169709, 0.07483568343955334, 1.4485115766149255, 0.0728263149576561, 1.3752500573676953, 0.07049136413829946, 1.298625749582679, 0.06824325424410622, 1.2193489845681136, 0.06582095403344015, 1.138025733875467, 0.06314751902234166, 1.0548582946199558, 0.06084252236708486, 1.7144695004547756, 0.11940298184648915, 1.6661432180305267, 0.11702536873232434, 1.612088738203573, 0.11462811234908235, 1.552767545779146, 0.11219900873693212, 1.4892764745904434, 0.10902057950304295, 1.4214462571985402, 0.10587577653875077, 1.3500670062966198, 0.10248755410832308, 1.2756115212143975, 0.09898370032790536, 1.1984515457771951, 0.09568213558577984, 1.1196489718494993, 0.09183539186337933, 1.0392003864352775, 0.08816879219824579, 1.6267363267820736, 0.1600270023348699, 1.5745122036212742, 0.15654434541589468, 1.5170197697958843, 0.15313338967854503, 1.4553997579886624, 0.1491108819434456, 1.390021778247703, 0.14466709218889734, 1.3212431978500232, 0.14002284032335086, 1.249537450002039, 0.13523105441596847, 1.1752042559020965, 0.13054636410812406, 1.0989160057119003, 0.12566138447730815, 1.02133329075452, 0.12024372533041004, 1.5319725305158807, 0.20442578647407403, 1.476977909976728, 0.199706478456275, 1.4175723167893004, 0.19470127031050574, 1.3546041029237836, 0.188992153995485, 1.2884403702991327, 0.18279746845340655, 1.2192245787942528, 0.17672369185414324, 1.1479290978475771, 0.17029790025344407, 1.0747495469411525, 0.16398664369913982, 1.0003392581236101, 0.1572958763431773, 1.4319066456991028, 0.25158806942868706, 1.3750821903990913, 0.24523816107210825, 1.3147465377387102, 0.23834415015953445, 1.2516898443435074, 0.23059954062063615, 1.1856949336792344, 0.22300228964082805, 1.1177333524825235, 0.21497294337051537, 1.048176627475603, 0.2066134908655834, 0.9769565819636076, 0.1986937650609893, 1.3293145098191994, 0.3005037221965657, 1.2719020717760812, 0.29221084898581606, 1.2116612255940085, 0.2831503953142519, 1.1490573933827157, 0.27358016662234425, 1.0845081246035433, 0.2638891518654642, 1.0188094526197202, 0.2536637657397752, 0.951724217111131, 0.24350535656648176, 1.2259314335753875, 0.35006767451005905, 1.168874300472285, 0.33974542585474266, 1.1098406286816516, 0.3282168141295506, 1.0490130480148494, 0.3162307138516238, 0.9868117066790869, 0.30439311974290983, 0.9242840295794844, 0.2917477571428014, 1.1234403610856294, 0.39952113060766964, 1.06807260649596, 0.38619265977430356, 1.0111761759684625, 0.37188240722355026, 0.9531468241095357, 0.35745320774585837, 0.8946387092331433, 0.3431588658207121, 1.0247030268677315, 0.4467325389142843, 0.9716007164124381, 0.4308208403938761, 0.9181871675395848, 0.41357497497034573, 0.8642985733046012, 0.3969595863819924, 0.9307906101251064, 0.49209707380918466, 0.880450290121062, 0.4731378009688773, 0.8309093585689126, 0.4542492821628245, 0.8422731024720952, 0.5345692330795209, 0.7960269508863682, 0.5141342033084638, 0.7641645911163861, 0.5711949090647838])
    SymCubatures.setweights!(cub, T[3.381943049050127e-5, 0.0001819479370076078, 0.00044221533520813147, 0.0007995387417695878, 0.0012315709025748404, 0.0017030391444416548, 0.002197696532114307, 0.002669555457504199, 0.003113293188855489, 0.003461510003813475, 0.0037201574537423857, 0.0038868079683849247, 0.003930554273168023, 0.003918487733868056, 0.003867287306121565, 0.003669072974083625, 0.0034374183005243942, 0.0003581622951928171, 0.0008110631911314635, 0.001230368903018792, 0.0016463327663960668, 0.0020717675388857207, 0.0024179331842699107, 0.0026825100227657085, 0.002936144757595323, 0.003199488708684168, 0.0033376557759707457, 0.0033600067940383754, 0.003377365469057866, 0.0033594887666466974, 0.0032820023741594877, 0.0029547089659298456, 0.002146836940106807, 0.003193961781349927, 7.846361154746041e-5, 0.00012241045754871993, 0.00016496670352559887, 0.00020505283537110058, 0.00024188634164504778, 0.0002738738347525823, 0.00030160366827910643, 0.00033114154384423876, 0.00034849368552871307, 0.0003634088211972808, 0.00037438529958092715, 0.0003853575770165873, 0.0003848978450595733, 0.0003812336513653601, 0.00037708056964667896, 0.0003553591651560254, 0.00028372730867824963, 0.00038208112091421645, 0.00047484592827108267, 0.0005597181784995092, 0.000633788052975396, 0.0006984018744803433, 0.0007634846525922919, 0.0008066764949867415, 0.0008390280811613128, 0.0008667023155715418, 0.0008846731437456087, 0.0008905013791921671, 0.0008767278012520108, 0.0008688034033130843, 0.0008310660637081254, 0.0005949137270164897, 0.0007389542953856064, 0.0008699141575357505, 0.0009851161812807954, 0.0010864925843599316, 0.0011804909101473622, 0.0012535917157260195, 0.0013009067401365635, 0.001344671910210316, 0.0013617489203409197, 0.0013764276184520057, 0.0013573142137603001, 0.0013352179405892518, 0.0013056421385656658, 0.0009926335099581737, 0.0011675740887462807, 0.001322622444669008, 0.0014592288271852993, 0.0015752178205413606, 0.0016783455020423347, 0.0017429981149472496, 0.0017963295322282022, 0.001816060075427626, 0.0018258221324085127, 0.0018186895694055545, 0.001768649051313106, 0.0017483749231574663, 0.0014478401171510562, 0.0016408032390730382, 0.001809995209861752, 0.0019449905292595391, 0.0020712971584815236, 0.002157358841216227, 0.0022098354869715476, 0.0022392348214667706, 0.0022351962411013776, 0.0022424758526442155, 0.0021863454566999293, 0.0021324758632605633, 0.0019326772539230484, 0.0021320316321422877, 0.0022897337754746567, 0.002431843358962115, 0.002545831110329685, 0.0026011747102776582, 0.00263515630382152, 0.0026259506545087065, 0.00260726851121544, 0.0025706486307784434, 0.002461678349057888, 0.0024243688215846763, 0.002607861979748161, 0.002753134908255028, 0.0028811280065258924, 0.002946623404684293, 0.002972091817838245, 0.0029802826240331723, 0.0029323065884280827, 0.002901840029380229, 0.0028111853278055985, 0.002877696846377732, 0.0030295256811881974, 0.003165077449178015, 0.0032596031773260274, 0.003286317234648916, 0.003293809434337844, 0.003238853380118991, 0.003151580839477258, 0.003108102403324153, 0.0032763778071919526, 0.0033981773248818936, 0.003497944459749652, 0.0035429516476305424, 0.0035256171062251415, 0.0034927549129013455, 0.0033870689303943796, 0.003294126027123625, 0.003587249705731255, 0.0036804126493895114, 0.0037542062092970653, 0.0037195177206319945, 0.0036485206274069395, 0.0035764399915457195, 0.003414107396707108, 0.0038020910540348635, 0.003871920030407605, 0.0038695801408132583, 0.0037627401897465927, 0.003647973736496454, 0.0035157765897702977, 0.003924652473401832, 0.003929718654176294, 0.003865946037063217, 0.003658127673031082, 0.0035141720615233336, 0.003914130282304548, 0.0038909311024247955, 0.0037353876176674447, 0.0035002892911415373, 0.0038878282979561743, 0.003788201013427329, 0.0035460740998105394, 0.0037146604609608197, 0.0033004733028740153, 0.0030447425498352494, 0.0036632134820301435])
    cub_degree = 68
    tol = 5e-14
  elseif q <=70 # 972 nodes
    SymCubatures.setparams!(cub, T[0.0021433273125948025, 0.011276187024038467, 0.027616418811839117, 0.05085347510348085, 0.08032824368718612, 0.11531748579016954, 0.1551610678285746, 0.19901126692240684, 0.2459946152554396, 0.2950807333406103, 0.3452759220356024, 0.395640311559258, 0.44514392833129524, 0.49287774745621415, 0.537893989721353, 0.5797547762649649, 0.6175961208524305, 0.6506101149360312, 1.9865789194896577, 0.0021432170059397897, 1.9702086658382976, 0.002141997080408097, 1.9467121920877606, 0.002133918674124924, 1.9162503094056065, 0.0021143698346870014, 1.8790277044519248, 0.002084966170939325, 1.835296069281524, 0.0020501730629971624, 1.7853573251953414, 0.002013242281928909, 1.7295619922362684, 0.001976444962803042, 1.6683145229945078, 0.0019334335764654741, 1.6020601207892848, 0.0018804574696830007, 1.531266090212204, 0.0018348334485619206, 1.4564700651155251, 0.0017765051708565286, 1.3782104750853306, 0.0017211025003536385, 1.2970716557741921, 0.0016619024503168962, 1.2136554563215998, 0.0015981064416240144, 1.1285762406887727, 0.001543114512798299, 1.0424934876425476, 0.0014750041293900896, 1.9610907624905856, 0.011268517563945485, 1.937647779714647, 0.01122505950059859, 1.9073049608908434, 0.011121732960675843, 1.8702646231922504, 0.010967432429539467, 1.8267637300956512, 0.010784996736275889, 1.7770941969364138, 0.010590680936278811, 1.721601473224793, 0.010394747866421482, 1.6607118301444612, 0.010169503945624965, 1.594884555416564, 0.009893757495486545, 1.524527575910886, 0.009646456151668186, 1.4502252876897972, 0.009347945377629842, 1.3724792127942842, 0.009051962225244946, 1.2918827510803188, 0.008741782424430241, 1.2090457324050528, 0.008411391628004148, 1.1245610597492746, 0.008110198694555615, 1.0391092948347576, 0.007761124750957645, 1.9214263236610531, 0.027505757534089095, 1.89129916169119, 0.027250454558776706, 1.8545851175000958, 0.02687395928656615, 1.811501859789648, 0.02642955307462069, 1.762330143713, 0.02595284682253201, 1.7074094711141794, 0.025463948383918335, 1.6471837252143307, 0.024915441018689113, 1.582138947019349, 0.024251301780902314, 1.5126175476995631, 0.02362004564719945, 1.439223429880611, 0.022915516668319303, 1.3624594113202217, 0.0221793809591561, 1.282880483229433, 0.02142112552026628, 1.2010735208088426, 0.020629817272613343, 1.1176036285282893, 0.019854467291250115, 1.033127021851222, 0.019025155579411753, 1.8685315506764137, 0.050377539446748235, 1.83232508984511, 0.04968664480453786, 1.7898643113773784, 0.04887211092120322, 1.7414255316981773, 0.04798805216579378, 1.6873521379701715, 0.04706579671637349, 1.6280988917759351, 0.04605693197406805, 1.5641954030451248, 0.044855235713224216, 1.4959400036372572, 0.043643935373052126, 1.4238543167152935, 0.04238120161375986, 1.348482972886396, 0.041015671333319355, 1.2703065314962532, 0.039608846077152045, 1.1898723874101147, 0.038175415906806245, 1.1077762302725067, 0.03669071162017883, 1.0245938577345195, 0.035186988591922666, 1.8038847633943886, 0.07924040574241646, 1.7623071760650166, 0.07795713858812661, 1.7148769292880082, 0.07654164142540758, 1.6619190249335172, 0.07504292077696402, 1.6038672913705496, 0.07343353781376673, 1.5412967533202782, 0.07155655406821401, 1.4745269351314219, 0.06957571109262303, 1.403975833243781, 0.06758032798154724, 1.3302666642880618, 0.06542263145388327, 1.2538571906592153, 0.06316732908952859, 1.1752088037382935, 0.060902656925265634, 1.094932507885164, 0.05851480515577861, 1.0135025951241667, 0.056126218885359236, 1.7288793893184626, 0.11347534098793931, 1.6826901332069486, 0.1114073468816263, 1.6310920971557705, 0.10919157985989829, 1.5744944670805456, 0.1068326786774084, 1.5135045798848292, 0.10414754728339601, 1.4484975348518896, 0.10124457047057928, 1.379817671978941, 0.09830876920743035, 1.3080774556352122, 0.09520981980931159, 1.2337821196123038, 0.09193241674736902, 1.1573161793514104, 0.0886214944474742, 1.0792253736049453, 0.08518489803327041, 0.9999537760250189, 0.0816944537184755, 1.6449710126594617, 0.1523299781325661, 1.595016312934746, 0.14926623230198177, 1.540193879272693, 0.14600333316680936, 1.4810930260237034, 0.14237479374237225, 1.4181282105966744, 0.13843924914633685, 1.3516616056553021, 0.134353557619636, 1.2821896935764274, 0.13013262476287255, 1.2102734173558642, 0.12569477243478744, 1.1363273494820916, 0.12112550000531869, 1.0607621993591336, 0.11648443154643791, 0.98403499143187, 0.11171034152796486, 1.5539231091435504, 0.19497990366971593, 1.5011570141577817, 0.19065623245110994, 1.4442084891802909, 0.18593743521782363, 1.3835102515583473, 0.18087114520107705, 1.3195460246314223, 0.17547870430266385, 1.2527068252064197, 0.1699155166379011, 1.183498286978701, 0.16417474715543648, 1.112406918972861, 0.15820483127866486, 1.0397357338823554, 0.15214820726707554, 0.9658787693328852, 0.14594270434159898, 1.4576028906296776, 0.24046651169407182, 1.4031193566044444, 0.23451072932282704, 1.3449651567439287, 0.22819474479127014, 1.2837309871745717, 0.22140755367260498, 1.2198178148474963, 0.21430224477059318, 1.1536266308663925, 0.20704778434795978, 1.085665401861952, 0.19958641184810103, 1.0162723586206346, 0.19191946987443884, 0.9456371583803885, 0.18411602636229207, 1.358144107995621, 0.2877537296611511, 1.302872218092006, 0.2800356956963996, 1.244604039620284, 0.27178903429722734, 1.1838256127879943, 0.26304005545535053, 1.1208973824823512, 0.25403509176207156, 1.056262445600276, 0.24490406136904252, 0.9903918315241468, 0.2355258873505743, 0.9233866698856641, 0.22595571386666233, 1.2573572135728464, 0.33598210521942284, 1.2023552212359183, 0.32616815462578347, 1.1449720878302743, 0.3157655869382598, 1.0856198578551373, 0.30491803096516557, 1.024673698757298, 0.29384985656854906, 0.9625024503826912, 0.28265226687196493, 0.8993200463150478, 0.2712135586478737, 1.157301236372667, 0.3840915123368327, 1.1035628855648851, 0.37197736448280366, 1.0479700788438102, 0.35932489059744893, 0.9910152103779608, 0.3461833612877917, 0.9328842353759786, 0.33292802231276597, 0.8736733633673931, 0.31956220215552456, 1.0599046857948493, 0.43115338031503975, 1.0082902516090548, 0.4166421631307594, 0.9554023734414802, 0.401527326356992, 0.9015628554518673, 0.3860145907410214, 0.8466163323497237, 0.3705443613833775, 0.9669203653421299, 0.4762529721855151, 0.9182605617648875, 0.4592153262193044, 0.868823242095962, 0.4415239100925702, 0.8183567072954239, 0.4237515257292552, 0.8800267241938006, 0.5187424180043569, 0.8351282564990445, 0.4989076136464215, 0.7891806272960713, 0.4788738817947392, 0.8000777299977473, 0.5577096118936079, 0.7589969777595176, 0.535551275667805, 0.7279301098240469, 0.5935358445717752])
    SymCubatures.setweights!(cub, T[3.023612091195227e-5, 0.0001628717081091782, 0.00039696985328674064, 0.0007195419318884649, 0.0011074033428338285, 0.001537779397397616, 0.0019881821589297186, 0.0024299476971052594, 0.0028384466761680995, 0.003185672787446136, 0.003461037644040299, 0.0036402361932349834, 0.003716393552251078, 0.0036886492704305204, 0.003549065339185308, 0.003342814838967559, 0.003026711781754969, 0.002591594167254663, 7.01889609116668e-5, 0.0001096804945252113, 0.00014793074671520513, 0.00018390174055440012, 0.00021691674803894133, 0.0002467996035463271, 0.0002735493655562577, 0.00029725098800494174, 0.0003167942918529314, 0.0003311957990158825, 0.00034322735240407276, 0.0003494378541066804, 0.0003525036367444727, 0.0003514018235891847, 0.0003460399731098008, 0.0003392832191885153, 0.0003269722066445213, 0.0002543752281621511, 0.0003428963919331676, 0.0004260636667823326, 0.0005023859829855994, 0.0005714412961342484, 0.0006331354326247507, 0.000687439001715061, 0.0007325058343629131, 0.0007660078770282342, 0.000792505321226204, 0.0008079817548808726, 0.0008141967155238817, 0.0008115767862204825, 0.0007998025965572012, 0.000781775636157717, 0.0007548227998393887, 0.0005346538630077154, 0.0006637682382555577, 0.0007821403118539513, 0.0008890822477953044, 0.000984213071670155, 0.001067194197091, 0.0011366457080635164, 0.0011889740892996768, 0.0012272673583368365, 0.0012527049753846928, 0.001261205065256619, 0.0012567577727234697, 0.0012403362009608576, 0.0012094767127721331, 0.0011701601567838667, 0.0008927621089426728, 0.0010515611156171295, 0.001194746639668632, 0.0013211273253983065, 0.0014300766162368026, 0.0015214803175785329, 0.0015914469286168102, 0.0016396464393674357, 0.0016739824000883026, 0.0016863669419363411, 0.0016799874023108326, 0.0016596238225715057, 0.0016179299282978476, 0.001566147016438327, 0.0013045590619521793, 0.0014823491239710587, 0.0016384469465377178, 0.001772160989377781, 0.0018840899250393444, 0.0019710645413464883, 0.0020292522158027277, 0.0020677556132473116, 0.0020836552872089114, 0.002074326514137343, 0.0020471686851745025, 0.0019992500377512567, 0.0019342081135371766, 0.0017482875216835601, 0.001932460446422174, 0.002089339167294333, 0.002219096392226603, 0.002320692938575604, 0.0023896133458681025, 0.0024284016778561074, 0.0024451142540035807, 0.0024341675591152663, 0.002397654805352886, 0.002345586515741198, 0.002269387549241992, 0.0021976424432120093, 0.00237513052029215, 0.002519899636770823, 0.002633591017995324, 0.002713857423195656, 0.0027539589343553543, 0.002766818194027064, 0.002754350891940241, 0.0027097173540253316, 0.0026482794269412817, 0.0025650536053046236, 0.0026263849121366993, 0.002784599965851446, 0.002907909059148797, 0.002996540070590101, 0.0030402468929051183, 0.003045376510912508, 0.003026076763210736, 0.0029787766914930545, 0.0029048798593427037, 0.0028155758029830306, 0.0030078408128006496, 0.003138879223089889, 0.0032323071256755275, 0.003282393231783815, 0.0032836907087489997, 0.0032500452850123553, 0.003196198951232896, 0.003114706584349293, 0.0030184062140424747, 0.0033224402945481206, 0.003417680780372259, 0.0034732991204815835, 0.003479832304474775, 0.003436490971193264, 0.0033630616323988977, 0.003273516033925578, 0.0031701438082577605, 0.0035532923113520526, 0.0036057428588534366, 0.003616446631337838, 0.003576156242004029, 0.003485767598267702, 0.003381530128597778, 0.0032749359386694973, 0.003684300992897163, 0.0036905134515635416, 0.0036544930174716296, 0.0035638425021092606, 0.0034398001913265555, 0.003330093019864756, 0.0037100507400315294, 0.003664858913592633, 0.003584529017001945, 0.00345367295976861, 0.003328972315015195, 0.0036237481101984584, 0.0035344673862233084, 0.0034125289954227735, 0.0032832234105619625, 0.00344918357138088, 0.0033197864315203944, 0.003207783426998228, 0.0032073536977279718, 0.0030817118453168335, 0.002900630011203609])
    cub_degree = 70
    tol = 5e-14
  elseif q <=72 # 1027 nodes
    SymCubatures.setparams!(cub, T[0.0020343201981817926, 0.01069799291994009, 0.02616337722097255, 0.0481459471991805, 0.07625963259220457, 0.10979163310468501, 0.14767657213897203, 0.18963027933635038, 0.2346060242939451, 0.2817865471718109, 0.3306397252408776, 0.3797767305699981, 0.4282477508900253, 0.4753386051737315, 0.5202333587808573, 0.5628294782677561, 0.6028451900848582, 0.6369275565721573, 0.9992813592103068, 0.9962939973405244, 0.9911175927578466, 0.9837350177870724, 0.9740333968103462, 0.9620521092616117, 0.9480322996809314, 0.9322583944057199, 0.9147376283817538, 0.8951582223480218, 0.8737127564574257, 0.8513457716524793, 0.8279606470727137, 0.803042726966684, 0.7771966953273195, 0.750435702156241, 0.7281443242866507, 0.7020333563272789, 1.9872614795575825, 0.0020333113167277025, 1.9717223959001142, 0.0020285364755874523, 1.9494084777683942, 0.0020196454923130614, 1.9204664844782322, 0.0020026994014855793, 1.8850704629876112, 0.001992505702482857, 1.8435060540932904, 0.001934002276763243, 1.7959619120096353, 0.0019150124536699803, 1.7427982468088172, 0.0018917443157672274, 1.6844030073468113, 0.0018389376076754105, 1.6211615665689083, 0.0017926103926421272, 1.553551946510177, 0.001744269402697717, 1.4820638254056835, 0.001707183192473356, 1.4072168050930243, 0.0016546304946564651, 1.3294755499667303, 0.0015957223516895771, 1.2492759875880775, 0.001546158217566695, 1.1671223213780553, 0.0014926828884080487, 1.0836303662434639, 0.00136942411713117, 1.9631023632285631, 0.010673192163523302, 1.9408633887566145, 0.010624082081507865, 1.9120450411379901, 0.010538741844691965, 1.8767801834880053, 0.010475829094485588, 1.8355508795242539, 0.010183607108995754, 1.7882423560452663, 0.010074934457787272, 1.7353891450000647, 0.009941515562060667, 1.6774412738323166, 0.009685694151506804, 1.6146657713669366, 0.009426610532667752, 1.5474350781786712, 0.009185321954589603, 1.4761712206460134, 0.008966970972796438, 1.4014122195903105, 0.008714015546555194, 1.3237673343810086, 0.008385878376943963, 1.2437779668573483, 0.00813384743354839, 1.1621905752793933, 0.007830549437821159, 1.0796812319699172, 0.007323663889582721, 1.925570852302602, 0.02603538095023418, 1.8969554063268834, 0.02583865589096441, 1.861957048201815, 0.025649403496337556, 1.8212766073360096, 0.02499429257123793, 1.774390750417704, 0.02469259134385301, 1.7220105016401954, 0.02432433450756655, 1.6645877982136688, 0.023770673665693035, 1.6024609725706944, 0.023093322964239642, 1.5358886163555205, 0.022529771462600764, 1.4655113756435536, 0.02194005948428518, 1.391807418783168, 0.02137158197342164, 1.3154311076205605, 0.020562533265060214, 1.236537923462079, 0.019919107208925098, 1.155809655600255, 0.019182716613643257, 1.0738887503572685, 0.018249526607243584, 1.8753958959603678, 0.047800571857037295, 1.8408085783648174, 0.047378683349105416, 1.8007472894703875, 0.04630133282470769, 1.7544545240283433, 0.045661580323540674, 1.7028303285784883, 0.0448993364371739, 1.6462742834197128, 0.044007665403815927, 1.585412814735559, 0.04272207030034999, 1.5200698960042358, 0.04168391751259823, 1.4509437216338799, 0.04055316830740465, 1.3783495844244316, 0.039492786661006814, 1.3031165548845611, 0.038090843802798806, 1.22549014721878, 0.036759967083540436, 1.1459973069231641, 0.03550884066707385, 1.0652809029090842, 0.034057556483060344, 1.813514965551407, 0.07549067050083696, 1.7742564707981197, 0.07398330254888609, 1.7289455363847368, 0.07283591598105119, 1.6784216291885425, 0.07150746110824657, 1.622815099394243, 0.07022784596335452, 1.563205152279674, 0.06821482019551725, 1.4991261501194841, 0.06648151176163, 1.431351621395142, 0.06469403754134198, 1.3603270810425585, 0.06287385242936866, 1.2866303326366493, 0.060826150272591244, 1.2106736345575397, 0.05856849746294543, 1.1326525305677053, 0.056635902266580336, 1.0535775530543423, 0.05446244407360917, 1.7420674120662472, 0.10785458254508465, 1.6979274518873162, 0.10603237884055339, 1.6486652576656518, 0.10396961110753201, 1.5941968820326442, 0.10217978730321566, 1.5360446228845634, 0.09941038308911482, 1.473619835207568, 0.09677406387723368, 1.4074542717590743, 0.09418469429810185, 1.3382686035911355, 0.09137610619372596, 1.2665002080541388, 0.08847512947562479, 1.1927429148636295, 0.08533558528367867, 1.117098575457867, 0.08231973196047075, 1.0400343379501154, 0.07919856891859815, 1.6618289622468843, 0.14504951000314018, 1.614033890418892, 0.14211066695143673, 1.5609984337334235, 0.13958512839636442, 1.5044040855458156, 0.13606683732136857, 1.4439396647094689, 0.13242711483668532, 1.3799607398198241, 0.12880254906250252, 1.313082362198592, 0.12491880600428794, 1.243616900210557, 0.12080220774530899, 1.1717018740257619, 0.11683738366892295, 1.0983196280096843, 0.112416283784707, 1.0236931306388406, 0.1080119315784288, 1.5748104567635766, 0.18574297371058146, 1.5237730022758618, 0.18217429071418484, 1.4689189837337322, 0.1778882894380211, 1.4103279392224723, 0.17322739978192803, 1.3483674597600073, 0.16830700598124212, 1.2835727371546557, 0.16317981608682755, 1.2166590435977482, 0.15769686162179325, 1.1475582994917133, 0.15257623532357387, 1.0771094384656243, 0.14692930029260773, 1.0053232593316992, 0.14100354327695916, 1.4820279780470564, 0.22963493521464373, 1.4290608434683527, 0.22446405402473116, 1.3726880617797048, 0.21876826099424518, 1.3133801155378533, 0.21255706847945924, 1.251445976687348, 0.20593524997615026, 1.1870675138231193, 0.19918705436599246, 1.120780700673696, 0.1923208949237423, 1.053077837747025, 0.18543747965939272, 0.9843528968355305, 0.17830479551060982, 1.386143042922063, 0.2755572949442648, 1.3323035037973392, 0.26869046339684977, 1.2753285112888835, 0.26123841414007004, 1.216063425190884, 0.25295191841693665, 1.1544414156238134, 0.24478307654164758, 1.09130862344127, 0.23621810030913806, 1.026797307982125, 0.2275497627216459, 0.9609962032182763, 0.21934168505351934, 1.2880365372496363, 0.32244551806286287, 1.23428244802781, 0.3137244196580037, 1.1780153312924038, 0.30410843428131346, 1.1195203992020384, 0.2940418905118154, 1.0594127613037962, 0.2839633118056177, 0.9986115920844257, 0.2734694717153918, 0.9366857439436637, 0.263190312954512, 1.190152168252434, 0.3693916462967947, 1.137181156515845, 0.3586777111868373, 1.08235555372719, 0.3468115394345797, 1.025521669039761, 0.33469231111464903, 0.9676081604056803, 0.3225977726856906, 0.9099665072281421, 0.30978205803163855, 1.0940388760315785, 0.41571402927763734, 1.0425105735599318, 0.4024580635558766, 0.9897549705046039, 0.38817040137671105, 0.9358163860150276, 0.37378763489059336, 0.8816538106209403, 0.35939161428046407, 1.001793110952995, 0.4602375728686977, 0.9523521141521669, 0.4443716950968946, 0.9027336926403373, 0.4273297019245618, 0.8527971797280974, 0.41101363416378267, 0.9141895803892408, 0.5024314541881847, 0.8672988889346868, 0.48391341512857444, 0.8212754482315994, 0.46556774544269974, 0.8315708999256756, 0.5422033930652448, 0.7884203490975132, 0.5224426512698988, 0.7578882737100152, 0.5763010152231002])
    SymCubatures.setweights!(cub, T[2.723993169547436e-5, 0.0001466012721047343, 0.0003565402150669151, 0.0006457245983689979, 0.0010005519358923995, 0.0013911706242107282, 0.0018159670827638352, 0.0022195414490605828, 0.002614205333724754, 0.0029159702211011924, 0.003188466550258188, 0.003386302558977106, 0.0035023447458884553, 0.0035052091449863554, 0.0034659614762080526, 0.0034254759253869375, 0.003301781447631705, 0.0030704433579455763, 0.0003096487594476438, 0.0006846244869113544, 0.0010398736078015027, 0.0014085609083232867, 0.0017606252205422273, 0.002086427362003224, 0.0023308354448864053, 0.002530124530396224, 0.0027428776941054002, 0.0029398966431727117, 0.003037055015098869, 0.002985531183167572, 0.002999131704644486, 0.0029988670100804973, 0.002909835065939314, 0.002591288833982568, 0.0018640404260170228, 0.0028292811472652706, 6.321486710953929e-5, 9.863586025526453e-5, 0.0001329900764203222, 0.00016553637422543168, 0.00019705078944743672, 0.00022153465673698542, 0.00024775654946727633, 0.00027111423728418783, 0.00028756405285787967, 0.0003013977449903387, 0.0003118344891082601, 0.0003208160925951768, 0.0003244837316638371, 0.00032365626820999866, 0.00032255307551000434, 0.0003175076461482068, 0.000296770730773018, 0.00022865037103621908, 0.0003080327105658478, 0.00038354877963597887, 0.00045564408253787644, 0.0005135502571197725, 0.0005729745722104711, 0.0006251908245559789, 0.0006656744321474193, 0.0006961585593722433, 0.0007232557235362792, 0.0007425299084007926, 0.000755453123806077, 0.0007507497376441038, 0.0007462000124307923, 0.0007288180972365191, 0.0007006219726788194, 0.0004799048406638583, 0.000597542801240989, 0.0007075524420939919, 0.000800583279629014, 0.0008908539521033464, 0.0009696573274070392, 0.0010379011189092582, 0.001082582617186552, 0.0011234498592742265, 0.0011471658804376621, 0.0011661315081739675, 0.001162508944914019, 0.0011518838094186268, 0.0011346456531735427, 0.0011106541631851496, 0.0008043111328324907, 0.0009505392708564485, 0.0010803475233341542, 0.0011973235805077054, 0.0012971871886040861, 0.0013899968631319636, 0.0014492661806179918, 0.0015009437109958612, 0.00153696630814932, 0.00155735045951061, 0.0015632794122214433, 0.001536782887137493, 0.0015207259644222926, 0.0014899921619949434, 0.0011795750938837903, 0.001343637499194605, 0.0014846394287310887, 0.0016060723303364703, 0.0017230028803008309, 0.001803690081625027, 0.0018606058596937693, 0.0019049884135311902, 0.0019189469007944084, 0.0019308324254202177, 0.0019076183539769673, 0.001873412737289186, 0.0018248788380729522, 0.0015890242472079768, 0.001755337043321866, 0.001897093048782707, 0.002028136886870807, 0.0021303530808308704, 0.002196471865124422, 0.002244106106387955, 0.002258205168403053, 0.0022506967828788155, 0.002243069380286745, 0.0021838381251314733, 0.0021281374855990813, 0.0020058923513077563, 0.00216699915412893, 0.002304124577146998, 0.002423318375636395, 0.0024994589621372687, 0.0025399648354268623, 0.002563383833574691, 0.0025502828711313277, 0.00254657356649308, 0.002483215644008329, 0.0023939037084407503, 0.002402997626077465, 0.0025439514519916763, 0.0026823523726886393, 0.002781516861654648, 0.002827909276860076, 0.0028376124424574255, 0.0028200895109886928, 0.0027704461151115546, 0.00273663413245083, 0.002662191001422362, 0.0027607797043394744, 0.002901648953432091, 0.0030007319705971786, 0.0030589940116094406, 0.003065079790193343, 0.003068103795728568, 0.0029971661053741637, 0.0029275540177747324, 0.0028870383861900684, 0.0030664132573248767, 0.003183989149320091, 0.0032639704234809475, 0.0032776398307566377, 0.0032469943839342134, 0.0032004157098957957, 0.003098110821691953, 0.003012082215712564, 0.003293199958695482, 0.003372528000691446, 0.0034391728753380667, 0.0033972489525227172, 0.0033255096773994446, 0.0032357851915973124, 0.003080826420921872, 0.0034361907773846757, 0.0034909600017741025, 0.003506226963270196, 0.0034200075067578246, 0.0033027484065695062, 0.0031429995182992577, 0.003511645041552288, 0.0035255419391118734, 0.003464379534041329, 0.003288389467477963, 0.00315111429733226, 0.003496558492765939, 0.00347398043582166, 0.0033249713879488904, 0.003112696057548909, 0.0034450853218008647, 0.0033674088790711485, 0.0031469442228341704, 0.0033155217913512867, 0.002926264669605203, 0.0027181994097547303, 0.0032508849872969177])
    cub_degree = 72
    tol = 5e-14
  else
    error("polynomial degree must be <= 40 (presently)\n")
  end

  mask = Int[]
  append!(mask, 1 : cub.numparams+cub.numweights)
  vtx = T[-1 -1; 1 -1; -1 1]
  Cubature.solvecubature!(cub, cub_degree, mask, tol=tol, hist=false)
  return cub, vtx
end


function getTriCubatureSparse(q::Int, T=Float64;
                              tol=10*eps(typeof(real(one(T)))))
  @assert( mod(q,4) == 0)
  function getGuess(ξ,η)
    return (ξ - ξ*η/3 -1), (η - ξ*η/3 - 1)
  end  
  cub_degree = q
  mask = zeros(Int64, (0))
  # get the underlying LGL nodes
  xlgl, wlgl = OrthoPoly.lglnodes(div(q+2,2), T)
  xlgl[:] += 1.0
  xlgl[:] *= 0.5
  # figure out how many edge, S21, and S111 nodes
  numedge = div(q,4)
  numS21 = numedge
  numS111 = div((numedge-1)*numedge,2)
  params = zeros(T, (numedge+numS21+2*numS111))
  weights = zeros(T, (1+numedge+numS21+numS111))
  ptr = 1
  wptr = 1
  weights[wptr] = wlgl[1]*wlgl[1]
  wptr += 1
  # S21 parameter guesses
  for i = 1:numS21
    x, y = getGuess(2*xlgl[i+1],2*xlgl[i+1])
    params[ptr] = (x + 1)
    weights[wptr] = wlgl[i+1]*wlgl[i+1]
    ptr += 1
    wptr += 1
  end
  # edge parameters
  for i = 1:numedge
    params[ptr] = xlgl[i+1]
    weights[wptr] = wlgl[i+1]*wlgl[1]
    ptr += 1
    wptr += 1
  end
  # S111 parameters
  for i = 2:numedge
    for j = i:numedge
      x, y = getGuess(2*xlgl[i+1],2*xlgl[j+1])
      params[ptr] = (x + 1)
      params[ptr+1] = (y + 1)
      weights[wptr] = wlgl[i+1]*wlgl[j+1]
      ptr += 2
      wptr += 1
    end
  end
  cub = SymCubatures.TriSymCub{T}(vertices=true, numedge=numedge,
                                  numS21=numS21, numS111=numS111)
  #weights[:] = 2./cub.numnodes
  SymCubatures.setweights!(cub, weights)

  wts = SymCubatures.calcweights(cub)
  weights[:] *= 2.0/sum(wts)

  weights[:] = 0.1/cub.numnodes
  
  SymCubatures.setweights!(cub, weights)
  println("sum wts = ",sum(SymCubatures.calcweights(cub)))

  SymCubatures.setparams!(cub, params)
  mask = zeros(Int64, (0))
  mask = SymCubatures.getInternalParamMask(cub)
  param_mask = zeros(Int64, (size(mask)))
  param_mask[:] = mask[:]
  append!(mask, (cub.numparams+1):(cub.numparams+cub.numweights))
  vtx = T[-1 -1; 1 -1; -1 1]
  success = false
  count = 0
  while success == false
    count += 1
    if count > 20
      break
    end
    SymCubatures.setweights!(cub, 0.05*rand(cub.numweights)/cub.numnodes)
    rand_params = params
    rand_params[param_mask] .*= (1.0 + 0.1*randn(size(param_mask)))
    SymCubatures.setparams!(cub, rand_params)
    try Cubature.solvecubature!(cub, cub_degree, mask, tol=tol, hist=true)
    catch
      success = false
    end
    if success
      if minimum(cub.weights) < 0.0
        success = false
        continue
      end
      break
    end
  end
    
  #Cubature.solvecubatureweights!(cub, cub_degree, tol=tol, hist=true)
  return cub, vtx
end

function getTriCubatureStaggered(q::Int, T=Float64; tol=10*eps(typeof(real(one(T)))))
  @assert(q>=2)

  cub_omega_degree = q-1
  cub_diage_degree = q
  mask_omega = zeros(Int64, (0))
  mask_diage = zeros(Int64, (0))

  if q<=2 
    #3 nodes 
    cub_omega = SymCubatures.TriSymCub{T}(vertices=false, numS21=1)      
    SymCubatures.setweights!(cub_omega, T[2/3])
    SymCubatures.setparams!(cub_omega, T[4/5])
    cub_omega_degree = 1

    #9 nodes
    cub_diage = SymCubatures.TriSymCub{T}(vertices = true,
                                    midedges = true,
                                    centroid = false,
                                    numedge = 0,
                                    numS21 = 1,
                                    numS111 = 0)                
    SymCubatures.setparams!(cub_diage, T[0.8])
    SymCubatures.setweights!(cub_diage, T[0.02906846240179568, 0.5337822671156006, 0.10381593714927037])
    cub_diage_degree = 2
    tol = tol

  elseif q<=4 
    #6 nodes
    cub_omega = SymCubatures.TriSymCub{T}(vertices=false, numS21=2, numS111=0)
    SymCubatures.setparams!(cub_omega, T[0.32232557541453105, 0.9443569106417902])
    SymCubatures.setweights!(cub_omega, T[0.5467137619817344, 0.11995290468493017])
    cub_omega_degree = 3

    #16 nodes
    cub_diage = SymCubatures.TriSymCub{T}(vertices = true,
                                    midedges = false,
                                    centroid = true,
                                    numedge = 1,
                                    numS21 = 2,
                                    numS111 = 0)                
    SymCubatures.setparams!(cub_diage, T[0.32232557541453044, 0.9443569106417901, 0.7236067977499789])
    SymCubatures.setweights!(cub_diage, T[0.02530741638066958, 0.2886037014101086, 0.1253889340745886, 0.05326140946739244, 0.36253138759954523])
    cub_diage_degree = 4
    tol = tol
  elseif q<=6
    #10 nodes
    cub_omega = SymCubatures.TriSymCub{T}(vertices=false, centroid=true,
                                    numS21=1, numS111=1)
    SymCubatures.setparams!(cub_omega, T[0.13862330627662678;0.14215944055500324;0.6226442585632832])
    SymCubatures.setweights!(cub_omega, T[0.11550472674301035;0.20924480696331949;0.39801697799105223])
    cub_omega_degree = 3

    #25 nodes
    cub_diage = SymCubatures.TriSymCub{T}(vertices = true,
                                    midedges = true,
                                    centroid = true,
                                    numedge = 1,
                                    numS21 = 2,
                                    numS111 = 1)                
    SymCubatures.setparams!(cub_diage, T[0.13862330627662678, 0.4765740799486239, 0.8273268353539885, 0.14215944055500324, 0.6226442585632832])
    SymCubatures.setweights!(cub_diage, T[0.0019263750791820421, 0.024548784813778073, 0.08983957252485644, 0.17471428465551592, 0.01736767414771074, 0.1514579304646647, 0.11395932110574997])
    cub_diage_degree = 6
    tol = tol
  elseif q<=8
    #15 nodes
    cub_omega = SymCubatures.TriSymCub{T}(vertices=false, numS21=3, numS111=1)
    SymCubatures.setparams!(cub_omega, T[0.08433122881886425;0.9485893782350207;0.4841719475189572;0.4231241172761849;0.0959626827429292])
    SymCubatures.setweights!(cub_omega, T[0.045386157905236965;0.1458284149509071;0.2543369199180239;0.11055758694624956])
    cub_omega_degree = 7

    #36 nodes
    cub_diage = SymCubatures.TriSymCub{T}(vertices = true,
                                    midedges = false,
                                    centroid = false,
                                    numedge = 2,
                                    numS21 = 5,
                                    numS111 = 1)                
    SymCubatures.setparams!(cub_diage, T[0.08433122881886425, 0.9485893782350207, 0.4841719475189572, 0.1988964028330802, 0.9207133051101312, 0.8825276619647324, 0.6426157582403226, 0.4231241172761849, 0.0959626827429292])
    SymCubatures.setweights!(cub_diage, T[0.0025791062643030858, 0.03102947540295087, 0.0835900339018613, 0.24334829933607555, 0.02777498089803061, 0.06557422568594258, 0.0031972306734735154, 0.009031313309203166, 0.09415672860607459])
    cub_diage_degree = 8
    tol = tol
  elseif q<=10
    #21 nodes
    cub_omega = SymCubatures.TriSymCub{T}(vertices=false, numS21=3, numS111=2, centroid=false)
    SymCubatures.setparams!(cub_omega, T[0.09037801956875349, 0.9630396695666223, 0.8072079596358801, 0.43658014194276157, 0.27398240252980893, 0.06084872345763977, 0.4441263310746356])
    SymCubatures.setweights!(cub_omega, T[0.05198714206463918, 0.10323440513804151, 0.18816014691671232, 0.09093907609523856, 0.07070341017839829])
    cub_omega_degree = 9

    #48 nodes
    cub_diage = SymCubatures.TriSymCub{T}(vertices = true,
                                    midedges = true,
                                    centroid = false,
                                    numedge = 2,
                                    numS21 = 6,
                                    numS111 = 2)                
    SymCubatures.setparams!(cub_diage, T[0.09199464041197784, 0.9547402313677645, 0.353676230440559, 0.9077053997995523, 0.19338251790806957, 0.5261554525794916, 0.9151119481392835, 0.7344243967353571, 0.4612850345245523, 0.06105167151116454, 0.1497962945573922, 0.5569626271469988])
    SymCubatures.setweights!(cub_diage, T[0.0013736390774556591, 0.025830304422168313, 0.01760407612922746, 2.7641301651682635e-5, 0.05424690442194781, 0.11921245370140136, 0.06034214056503969, 0.15865142740276447, 0.010769947391019101, 0.0094654627709123, 0.029664228393072672, 0.06478940126750116])
    cub_diage_degree = 10
    tol = 5e-15
  end

  vtx = T[-1 -1; 1 -1; -1 1]
  mask_omega = 1:(cub_omega.numparams+cub_omega.numweights)
  Cubature.solvecubature!(cub_omega, cub_omega_degree, mask_omega, tol=tol)

  mask_diage = SymCubatures.getInternalParamMask(cub_diage)
  append!(mask_diage, (cub_diage.numparams+1):(cub_diage.numparams+cub_diage.numweights))
  Cubature.solvecubature!(cub_diage, cub_diage_degree, mask_diage, tol=tol)

  return cub_omega, cub_diage, vtx
end

"""
### Cubature.tetcubature{T}

This high-level function computes and returns a symmetric cubature of requested
accuracy on the right tetrahedron.

**Inputs**

* `q`: maximum degree of polynomial for which the cubature is exact
* `T`: the data type used to represent the cubature
* `internal`: if true, all nodes are strictly internal (default false)
* `facequad`: if true, the cubatures' face nodes coincide with a quadrature
* `tol`: tolerance with which to solve the cubature

**Outputs**

* `cub`: a symmetric cubature for the right tetrahedron
* `vtx`: vertices for the right tetrahedron

"""
function tetcubature(q::Int, T=Float64; internal::Bool=false,
                     facequad::Bool=false,
                     tol=10*eps(typeof(real(one(T)))))

end

"""
### Cubature.getTetCubatureGamma{T}

Returns a cubature rule and vertices for the SBP Gamma operators on tetrahedra;
these are operators with (p+1)(p+2)/2 nodes on each face, where, typically, p =
(`q`+1)/2.

**Inputs**

* `q`: maximum degree of polynomial for which the cubature is exact
* `T`: the data type used to represent the cubature
* `tol`: tolerance with which to solve the cubature

**Outputs**

* `cub`: a symmetric cubature for the right tetrahedron
* `vtx`: vertices for the right tetrahedron

"""
function getTetCubatureGamma(q::Int, T=Float64;
                             tol=10*eps(typeof(real(one(T)))))
  cub_degree = -1
  if q <= 1
    # P1 (vertices only); 2nd order cubature
    cub = SymCubatures.TetSymCub{T}()
    SymCubatures.setweights!(cub, T[1/3])
    cub_degree = 1
  elseif q <= 3
    # P2 + 1 bubble node; 4th order cubature
    cub = SymCubatures.TetSymCub{T}(midedges=true, centroid=true)
    SymCubatures.setweights!(cub, T[1/45 4/45 32/45])
    cub_degree = 3
  elseif q <= 5
    # P3 + 4 bubble nodes; 6th order cubature
    cub = SymCubatures.TetSymCub{T}(facecentroid=true, numedge=1, numS31=1)
    SymCubatures.setweights!(cub, T[0.004421633248304776 0.20653163611605146
                                    0.06935370366814568 0.0176754534336105])
    SymCubatures.setparams!(cub, T[0.45720884759834435 0.30480589839889616])
    cub_degree = 5
  elseif q <= 7
    # P3 + 11 bubble nodes; 8th order cubature
    cub = SymCubatures.TetSymCub{T}(midedges=true, centroid=true, numedge=1,
                                    numfaceS21=1, numS31=1, numS22=1)
    # SymCubatures.setweights!(cub, T[0.0015106273303336273,0.060490542374353584,
    #                                 0.004038881996228382, 0.10344930834722398,
    #                                 0.005696088152131421, 0.02424296133613638,
    #                                 0.08113091859465722])
    # SymCubatures.setparams!(cub, T[0.28418700275470193,0.21742832019555544,
    #                                0.25737274681480826,0.45008848310824695])
    SymCubatures.setweights!(cub, T[0.0015106273303336273,0.060490542374353584,
                                    0.004038881996228382, 0.10344930834722398,
                                    0.02424296133613638,0.005696088152131421,
                                    0.08113091859465722])
    SymCubatures.setparams!(cub, T[0.28418700275470193,0.21742832019555544,
                                   0.45008848310824695,0.25737274681480826])
    cub_degree = 7
  else
    error("polynomial degree must be 1, 3, 5, or 7 (presently)\n")
  end
  mask = 1:(cub.numparams+cub.numweights)
  vtx = T[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
  Cubature.solvecubature!(cub, cub_degree, mask, tol=tol)
  return cub, vtx
end

"""
### Cubature.getTetCubatureOmega{T}

Returns a cubature rule and vertices for the SBP Omega operators on tetrahedra;
these are cubatures that are analogous to Gauss-Legendre in 1D, and they are
strictly internal to the tet. 

**Inputs**

* `q`: maximum degree of polynomial for which the cubature is exact
* `T`: the data type used to represent the cubature
* `tol`: tolerance with which to solve the cubature

**Outputs**

* `cub`: a symmetric cubature for the right tetrahedron
* `vtx`: vertices for the right tetrahedron

"""
function getTetCubatureOmega(q::Int, T=Float64;
                             tol=10*eps(typeof(real(one(T)))))
  cub_degree = -1
  if q<=1 
    cub = SymCubatures.TetSymCub{T}(vertices=false, numS31=1)
    SymCubatures.setweights!(cub, T[1/3])
    SymCubatures.setparams!(cub, T[3/5])
    cub_degree = 1
  elseif q <= 2
    # P1; 3rd order cubature
    cub = SymCubatures.TetSymCub{T}(vertices=false, numS31=1)
    SymCubatures.setweights!(cub, T[1/3])
    SymCubatures.setparams!(cub, T[(1 - sqrt(5)/5)*3/4])
    cub_degree = 2
  elseif q <= 3 
    #10 nodes 
    cub = SymCubatures.TetSymCub{T}(vertices=false, numS31=1, numS22=1)
    SymCubatures.setweights!(cub, T[0.15107940134333817, 0.1215026213266635])
    SymCubatures.setparams!(cub, T[0.3540971084487333, 0.18203145362159834])

    # 10 nodes, P2; 4th order cubature
    # cub = SymCubatures.TetSymCub{T}(vertices=false, numS31=1, numS22=1)
    # SymCubatures.setweights!(cub, T[0.06483158243276162;
    #                                 0.17900116726703835])
    # SymCubatures.setparams!(cub, T[0.22511815489558668;
    #                                0.18771315212883505])

    #SymCubatures.setweights!(cub, T[0.1302091416313459;
    #                                0.13541612780132486])
    #SymCubatures.setparams!(cub, T[0.33398409622579817;
    #                               0.18658191164952043])

    #8 nodes (too few nodes to support degree 2 operators)
    # cub = SymCubatures.TetSymCub{T}(vertices=false, numS31=2)
    # SymCubatures.setweights!(cub, T[0.20211619905201966, 0.13121713428131376])
    # SymCubatures.setparams!(cub, T[0.36740882067683256, 0.9956016068076595])

    cub_degree = 3
  elseif q <= 5 
    #20 nodes
    cub = SymCubatures.TetSymCub{T}(vertices=false, numS31=2, numS211=1)
    SymCubatures.setparams!(cub, T[0.27308286881050636, 0.931780068315618, 0.09856058998260914, 0.7856072614781071])
    SymCubatures.setweights!(cub, T[0.0905678982461762, 0.14714726910635742, 0.03187272199359993])

    # # 20 nodes
    # #P3; 6th order cubature
    # cub = SymCubatures.TetSymCub{T}(vertices=false, numS31=2, numS211=1)
    # #SymCubatures.setweights!(cub, T[0.061630217648090097;
    # #                                0.13793513058238085;
    # #                                0.04458932836762084])
    # #SymCubatures.setparams!(cub, T[0.24722530396402584;
    # #                               0.9298082909679131;
    # #                               0.11664936229736803;
    # #                               0.6505900754758551])
    # SymCubatures.setweights!(cub, T[0.1330522438256425;
    #                                 0.027128275352375174;
    #                                 0.05771760471843853])
    # SymCubatures.setparams!(cub, T[0.9291909288071043;
    #                                0.18740453102806648;
    #                                0.12393623056047982;
    #                                0.5651871191159269])
    
    #14 nodes (too few nodes to support degree 2 operators)
    # cub = SymCubatures.TetSymCub{T}(vertices=false, numS31=2, numS22=1)
    # SymCubatures.setweights!(cub, T[0.1502505676240212, 0.09799072415514933, 0.05672802770277528])
    # SymCubatures.setparams!(cub, T[0.9326577577899019, 0.27820575093267375, 0.9089925917487007])
    cub_degree = 5
  elseif q <= 7
    #35 nodes
    # cub = SymCubatures.TetSymCub{T}(vertices=false, centroid=true, numS31=1, numS22=1, numS211=2)
    # SymCubatures.setweights!(cub, T[0.056439441613289405, 0.0425292371104768, 0.04960950763777952, 0.01081436110653779, 0.12731371928550783])
    # SymCubatures.setparams!(cub, T[0.9471034493346083, 0.10097964519679277, 0.3776676620520021, 0.09432140072199578, 0.04253094508296646, 0.29327762763696985])

    # 38 nodes
    # P4; 8th order cubature
    cub = SymCubatures.TetSymCub{T}(vertices=false, numS31=2, numS22=1,
                                    numS211=2)
    # minimum node distance for this one is 0.17819
    # SymCubatures.setweights!(cub, T[0.02832965568227839;
    #                                 0.048583147669377845;
    #                                 0.039177188602071006;
    #                                 0.03021820473246459;
    #                                 0.03566671096039236])
    # SymCubatures.setparams!(cub, T[0.18167711419341304;
    #                                0.5398647398205032;
    #                                0.7170540544966304;
    #                                0.0881323679975843;
    #                                0.5992257377201948;
    #                                0.4688384403943167;
    #                                1.0098301020743294])
    # minimum node distance for this one is 0.17937
    # SymCubatures.setweights!(cub, T[0.05138889021100641;
    #                                 0.028389644845730116;
    #                                 0.03630699901869689;
    #                                 0.03614773635849856;
    #                                 0.030217030224352005])
    # SymCubatures.setparams!(cub, T[0.5465722016176291;
    #                                0.18183572000191664;
    #                                0.7209450878473552;
    #                                0.46859987367504574;
    #                                0.0537103011531615;
    #                                0.08816163364683494;
    #                                0.5991524190716274])
    # minimum node distance for this one is 0.2000058
    SymCubatures.setweights!(cub, T[0.0680746674820208;
                                    0.028423617986180167;
                                    0.018310980956037618;
                                    0.025461800630426842;
                                    0.044327724846598505])
    SymCubatures.setparams!(cub, T[0.5821093910011628;
                                   0.18262906602002502;
                                   0.17262216834048155;
                                   0.08099388793344552;
                                   0.5908415591286749;
                                   0.4612405506278837;
                                   0.07057889678019784])
    cub_degree = 7
  elseif q <= 9 
    #65 nodes 
    cub = SymCubatures.TetSymCub{T}(vertices=false, centroid=true, numS31=1, numS22=0, numS211=5)
    SymCubatures.setweights!(cub, T[0.06845875962635535, 0.03996133054703136, 0.004908706029927262, 0.012435872193517524, 0.018315521918100697, 0.011842283605048317, 0.009933723304410471])
    SymCubatures.setparams!(cub, T[0.5659668188394116, 0.7807167565199028, 0.0846814812674472, 0.04046668289033655, 0.20685562107102665, 0.061511434675817794, 0.6927211311947881, 0.14342099002644818, 0.3786528536651681, 0.3596513879533776, 0.02013276646659172])

    #59 nodes (Vandermonde matrix only has 54 (instead of 56) nonlinearly independent columns for degree 5 SBP operator)
    # P4; 10th order cubature
    # cub = SymCubatures.TetSymCub{T}(vertices=false, centroid=true, numS31=4, numS22=1, numS211=3)
    # SymCubatures.setweights!(cub, T[0.007355949263994053, 0.0318527786973356, 0.0400902106173932, 0.003960361288674464, 0.04886596434433109, 0.010769970544765473, 0.013078308790859482, 0.029122923275671307, 0.07144591646220326])
    # SymCubatures.setparams!(cub, T[0.11923697815179675, 0.5157517849125112, 0.9668538043763143, 0.17887605340644788, 0.21738047101393035, 0.9186717278941127, 0.16011659785744337, 0.0656000901160755, 0.43100467872680415, 0.3627623183625254, 1.2011673771365483])

    #65 nodes (Vandermonde matrix only has 54 (instead of 56) nonlinearly independent columns for degree 5 SBP operator)
    # cub = SymCubatures.TetSymCub{T}(vertices=false, centroid=true, numS31=4, numS22=2, numS211=3)
    # SymCubatures.setweights!(cub, T[0.010503533716038056, 0.03669835311481454, 0.031833839676230956, 0.006710436752072768, 0.002944793014934955, 0.04356959384164538, 0.01294188126646907, 0.028765065998142422, 0.013737475454947027, 0.04592928652252455])
    # SymCubatures.setparams!(cub, T[0.12996543337694574, 0.9673676540009212, 0.5162909518831329, 0.8593193682825483, 0.9704688160101619, 0.7644670309084176, 0.06536028636815383, 1.4354151453150563, 0.35663275933640914, 0.07397620486118074, 0.8880310794454519, 0.02401759075563747])

    cub_degree = 9
    tol=1e-14
  end
  mask = 1:(cub.numparams+cub.numweights)
  vtx = T[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
  Cubature.solvecubature!(cub, cub_degree, mask, tol=tol)
  return cub, vtx
end

function getTetCubatureOmegaLG(q::Int, T=Float64;
                             tol=10*eps(typeof(real(one(T)))))
  cub_degree = -1
  mask = zeros(Int64, (0))
  p=convert(Int, floor(q/2))
  ne = p+1
  numedge = convert(Int,(ne-mod(ne,2))/2)
  tol=5e-14
  cub = SymCubatures.TetSymCub{T}(vertices=false,
                                  midedges=false,
                                  numfaceS21=0,
                                  numedge=0,
                                  numfaceS111=0,
                                  facecentroid =false,
                                  numS31 = (1+mod(ne,2))*numedge,
                                  numS22 = mod(ne,2)*numedge,
                                  numS211 = convert(Int, (1+2*mod(ne,2))/(1+mod(ne,2)) * (numedge*(numedge-1))),
                                  numS1111 = convert(Int, 1/6 *((numedge-1)^3 - (numedge-1))), 
                                  centroid = convert(Bool, mod(ne,2)))
  
  if q<=2 
    SymCubatures.setparams!(cub, T[0.41458980337503154])
    SymCubatures.setweights!(cub, T[0.33333333333333337])
    cub_degree = 2
    tol = 1e-14
  elseif q <= 4
    SymCubatures.setparams!(cub, T[0.258411222979547, 0.9479685828265005, 0.8548881205292567])
    SymCubatures.setweights!(cub, T[0.0839659833549715, 0.07732176889453111, 0.10006762726446286, 0.08777656074854054])
    cub_degree = 4
  elseif q <= 6
    SymCubatures.setparams!(cub, T[0.1736418582760552, 0.6090285425552456, 0.09821576246902468, 1.2137233208634202, 0.47084657041102673, 0.08683024714379607])
    SymCubatures.setweights!(cub, T[0.025779302085158017, 0.05398405242331654, 0.03629419843699386, 0.04822912783795905])
    cub_degree = 6
  elseif q <= 7
    # cub = SymCubatures.TetSymCub{T}(vertices=false,
    #                               midedges=false,
    #                               numfaceS21=0,
    #                               numedge=0,
    #                               numfaceS111=0,
    #                               facecentroid =false,
    #                               numS31 = 1,
    #                               numS22 = 0,
    #                               numS211 = 3,
    #                               numS1111 = 0, 
    #                               centroid = true)
    # SymCubatures.setparams!(cub, T[0.1819592217914153, 0.05924209898012428, 1.3059663722514356, 0.43643295260855886, 0.10211288185342225, 0.843729781027974, 0.07575686050558157])
    # SymCubatures.setweights!(cub, T[0.02747072575893043, 0.016463810523532317, 0.0579151962914137, 0.01672715877835158, 0.1301764431780407])
    cub = SymCubatures.TetSymCub{T}(vertices=false,
                                  midedges=false,
                                  numfaceS21=0,
                                  numedge=0,
                                  numfaceS111=0,
                                  facecentroid =false,
                                  numS31 = 2,
                                  numS22 = 1,
                                  numS211 = 2,
                                  numS1111 = 0, 
                                  centroid = true)
    SymCubatures.setparams!(cub, T[0.18291677528623593, 0.5692831253357231, 0.8155856559865274, 0.07948262421311397, 1.2472768800723302, 0.4604765824853171, 0.07044535143802749])
    SymCubatures.setweights!(cub, T[0.02857209016010688, 0.05773460666891144, 0.02061191857272869, 0.02504823647158195, 0.044380319217177305, 0.031292366315777455])
    cub_degree = 7
  elseif q <= 8 
    SymCubatures.setparams!(cub, T[0.12570564553075184, 0.9831219906422596, 0.465608785401795, 0.8945656159690933, 0.9426862477840718, 0.7424672320740971, 0.07255261807032559, 1.4330029713468562, 0.36888830237894754, 0.05962095898671999, 0.8350664216465402, 0.05048623399323985])
    SymCubatures.setweights!(cub, T[0.009612842761491094, 0.01456910080856551, 0.03782620833766524, 0.03179143512660305, 0.011257965118796539, 0.03316921331538733, 0.015140889804851157, 0.02267722431734249, 0.01786049866374207, 0.023428561159700174])
    cub_degree = 8
  elseif q <= 9
    cub = SymCubatures.TetSymCub{T}(vertices=false,
                                  midedges=false,
                                  numfaceS21=0,
                                  numedge=0,
                                  numfaceS111=0,
                                  facecentroid =false,
                                  numS31 = 3,
                                  numS22 = 2,
                                  numS211 = 3,
                                  numS1111 = 0, 
                                  centroid = true)
    SymCubatures.setparams!(cub, T[0.12678586138792855, 0.4812443504978587, 0.9285845347634994, 0.95979440432313, 0.8170957452233956, 0.07347200166352832, 1.4284984182739606, 0.36632845951209003, 0.05663433748603191, 0.835565009120679, 0.0007538063415509122])
    SymCubatures.setweights!(cub, T[0.009845077724761996, 0.046498539894039057, 0.0509861563361617, 0.007695573686468984, 0.04058411645430292, 0.015347138948068498, 0.02150469555501786, 0.009586679079486007, 0.05707393367798336])
    cub_degree = 9
    tol = 1e-14
  elseif q <= 10
    SymCubatures.setparams!(cub, T[0.08872609101537908, 0.3737939913588266, 0.6645643197060291, 0.057862724047801045, 1.577677719099153, 0.049858612480221616, 1.1857406042552787, 0.29013021437417563, 0.04782437766934095, 0.5474332392465812, 0.04091749119678137, 0.2313766578185174, 0.9598627989847222, 0.5105686795460095, 0.2007678166765484, 0.655346348282158, 1.0511509509721237, 0.2476053990773419])
    SymCubatures.setweights!(cub, T[0.00347893406868075, 0.019481222382805825, 0.014365093530277609, 0.007238136940769622, 0.007314901031380741, 0.01145345263812562, 0.010359851131662852, 0.022248784528790727, 0.016322103350822105, 0.01186606574781905])
    # cub = SymCubatures.TetSymCub{T}(vertices=false,
    #                               midedges=false,
    #                               numfaceS21=0,
    #                               numedge=0,
    #                               numfaceS111=0,
    #                               facecentroid =false,
    #                               numS31 = 3,
    #                               numS22 = 0,
    #                               numS211 = 5,
    #                               numS1111 = 1, 
    #                               centroid = false)
    # SymCubatures.setparams!(cub, T[0.09757422691102068, 0.3922726784233176, 0.639555450266003, 0.047185196130793604, 1.5812978340470165, 0.06235285427674928, 1.1654153128534939, 0.27622151474996287, 0.0592084839193545, 0.2150533517407038, 0.9874854295286944, 0.5397057087640437, 0.10782171505432232, 0.6553613833589721, 1.0209141724998967, 0.3096939308709218])
    # SymCubatures.setweights!(cub, T[0.004419163160480125, 0.010380654813474326, 0.0311671556547637, 0.005379620447785844, 0.010820119029501769, 0.015097910061858235, 0.029743498060344084, 0.01910651814557053, 0.00782056041157232])       
    cub_degree = 10
    tol = 5e-14
  elseif q <= 12
    SymCubatures.setparams!(cub, T[0.0692215234998783, 0.986193517827988, 0.30079042382029686, 0.9430488613787327, 0.5633201087288852, 0.85711593356254, 0.9691933850172209, 0.8534715937413939, 0.6742758291498514, 0.04503185158907004, 1.6721058867408332, 0.04042481010467566, 1.360145161974983, 0.23049239642304453, 0.03802737005104062, 0.9053957133736139, 0.02913949506966734, 0.46533112493990203, 0.02897411392562341, 0.7978073960625174, 0.026785089980937524, 0.17803879977009437, 1.165165263380566, 0.4264328481689769, 0.1541393802798876, 0.7553402876826483, 0.14304867521129053, 0.5332544823607697, 1.2298610867796809, 0.20363269165278214])
    SymCubatures.setweights!(cub, T[0.001636952898635648, 0.004862224070312301, 0.010752389756258229, 0.008782534636178774, 0.0159574533303034, 0.01361242034065208, 0.0026404038669245924, 0.009356805249132913, 0.011315662125607197, 0.0034440012226790054, 0.003907890870696643, 0.0060707200657213285, 0.004793462890843898, 0.00618023704558517, 0.005447890843921627, 0.012561679707330227, 0.013541714484928656, 0.01192648244597215, 0.006253542677190052, 0.0064662265692781926])
    cub_degree = 12
    tol = 1e-14
  elseif q <= 13
    cub = SymCubatures.TetSymCub{T}(vertices=false,
                                  midedges=false,
                                  numfaceS21=0,
                                  numedge=0,
                                  numfaceS111=0,
                                  facecentroid =false,
                                  numS31 = 2,
                                  numS22 = 0,
                                  numS211 = 5,
                                  numS1111 = 5, 
                                  centroid = true)
    SymCubatures.setparams!(cub, T[0.9688668328238884, 0.3049040881460405, 0.1503459999947402, 0.03655673805222385, 0.41548749872704677, 0.19741924888398893, 0.01378741117939847, 1.8585705727655366, 0.043758254989289085, 1.523476086804378, 0.7006586963427017, 0.2157450074988152, 0.3742893627319491, 1.376963564094081, 0.22951358334780664, 0.7499955475809141, 1.1549825951300592, 0.08537696258483667, 0.6303784510062743, 0.9740912009544433, 0.3676604466945102, 0.7582576123673513, 0.959080562171818, 0.19131110641116408, 0.47023158026323814, 1.228722571595374, 0.09198475992361727])
    SymCubatures.setweights!(cub, T[0.013995590796848917, 0.012538295363476733, 0.003665810918921296, 0.018984459877616267, 0.0005855242258414606, 0.004020070603317004, 0.017499238077053066, 0.0028439048700980205, 0.0029218284990632755, 0.007091075578684134, 0.0061666729289612315, 0.008156268610558765, 0.037822532562271534])
    
    # cub = SymCubatures.TetSymCub{T}(vertices=false,
    #                               midedges=false,
    #                               numfaceS21=0,
    #                               numedge=0,
    #                               numfaceS111=0,
    #                               facecentroid =false,
    #                               numS31 = 2,
    #                               numS22 = 2,
    #                               numS211 = 4,
    #                               numS1111 = 5, 
    #                               centroid = true)
    # SymCubatures.setparams!(cub, T[0.9338078891037311, 0.34112559281098653, 0.8314894335624173, 0.6946191494984467, 0.18067968097605358, 0.04258098400970393, 0.024200719939166476, 1.8335330387438205, 0.042059951266458076, 1.5281394110044122, 0.8597237093780349, 0.053903729546270866, 0.4451388766644147, 1.2940377259324143, 0.24220701965317093, 0.7504895284284939, 1.1544729133879537, 0.07901159084889105, 0.623388242360332, 0.9131571130730572, 0.4335911615129601, 0.486369134162889, 0.9289853377876371, 0.3886762049556831, 0.4941550373361054, 1.203844537713373, 0.2125428835865013])
    # SymCubatures.setweights!(cub, T[0.01734049859822038, 0.01484969783668714, 0.009896828804573805, 0.01720458204249585, 0.005413717432789543, 0.0009710352622061783, 0.0038092989709655307, 0.008317441223056525, 0.0033275817972999075, 0.0029212140606067263, 0.0060359983698216, 0.011267472107712348, 0.009127978797295835, 0.03550028465739826])
    
    cub_degree = 13
    tol = 3e-14
  elseif q <= 14
    SymCubatures.setparams!(cub, T[0.051271884101886864, 0.25735535662514486, 0.4900287215057178, 0.690316089629788, 0.0362869645445464, 1.7463811261860935, 0.03394529009462525, 1.4940609162683463, 0.028252804546461275, 1.1610141982814228, 0.18637495921923083, 0.0336856891163145, 0.3932109043019217, 0.025514658444716336, 0.5807847323328273, 0.023162740211350286, 0.15440052349969263, 1.2901033439571483, 0.13980511919065022, 1.0202125662689325, 0.36746943069477306, 0.13390059597564266, 0.5544255086161157, 0.1139994624737987, 0.3033159336869329, 0.8201685858116599, 0.5131919675806913, 0.26236757472779787, 0.43626571490526983, 1.3640647820287295, 0.17058978751152823, 0.7499988472134831, 1.0767755556363754, 0.14549309644233047, 0.6760665838114465, 0.9576119374338002, 0.341513302651534, 0.6370823860831846, 0.9107175833731902, 0.3227293675312617])
    SymCubatures.setweights!(cub, T[0.0006811579314556321, 0.006346613615371775, 0.01083026631702925, 0.0061197722968127, 0.0017536820402800886, 0.0022662381841157507, 0.0019397210002682227, 0.0035631975226471254, 0.004122041433944076, 0.003307842644291322, 0.007769067912818856, 0.007125461647705706, 0.008697793802586737, 0.006328218515286998, 0.009483026976832044, 0.0071923988134652405, 0.003904379964256275, 0.0036554203363115298, 0.004016920995583651, 0.008208187318838125])
    # cub = SymCubatures.TetSymCub{T}(vertices=false,
    #                               midedges=false,
    #                               numfaceS21=0,
    #                               numedge=0,
    #                               numfaceS111=0,
    #                               facecentroid =false,
    #                               numS31 = 2,
    #                               numS22 = 0,
    #                               numS211 = 10,
    #                               numS1111 = 4, 
    #                               centroid = false)
    # SymCubatures.setparams!(cub, T[0.045975948972661405, 0.25711822174514515, 0.03713459008520997, 1.752653176319518, 0.036859441217342095, 1.4954321820653618, 0.016764917366921834, 1.1768679911239268, 0.19292345386501383, 0.0341546997588657, 0.3852957669601553, 0.05028087452341492, 0.1523020102411455, 1.2855425238100968, 0.35732266036891663, 0.23409849627611704, 0.5886258079198655, 0.057012540652161936, 0.36050497169381884, 0.6802074836303502, 0.4973399099037253, 0.1836707909892999, 0.44649857682216615, 1.3601587173405376, 0.1732757689129945, 0.7478349555906457, 1.0823826096512206, 0.1334657424720997, 0.6702153480890423, 0.9445683748505124, 0.3733428439831655, 0.662412807962613, 0.936852688803749, 0.27508661319008554])
    # SymCubatures.setweights!(cub, T[0.0005246868271175226, 0.006216322411427566, 0.0018350712284524377, 0.002618366624379609, 0.0010807947821359365, 0.0038141869561627594, 0.008016791500675456, 0.009933486302484066, 0.008258638122729442, 0.005607376727217563, 0.01040903152317131, 0.012443750173156425, 0.002945720760238543, 0.005195586184370721, 0.0034271065468949318, 0.010854893554011358])
    cub_degree = 14
    tol = 1e-14
  elseif q <= 15
    cub = SymCubatures.TetSymCub{T}(vertices=false,
                                  midedges=false,
                                  numfaceS21=0,
                                  numedge=0,
                                  numfaceS111=0,
                                  facecentroid =false,
                                  numS31 = 2,
                                  numS22 = 0,
                                  numS211 = 12,
                                  numS1111 = 4, 
                                  centroid = false)
    SymCubatures.setparams!(cub, T[0.05867939297014181, 0.23881256695190567, 0.03354387198694401, 1.7331713434004834, 0.03604795781246912, 1.470028274624915, 0.022564579625897484, 1.1628630223698426, 0.1797765824831259, 0.03168956409050956, 0.3784754441080742, 0.048267348508157885, 0.5751473620737454, 0.031131808550884852, 0.14641437390822062, 1.3280128449769013, 0.13320993217627147, 1.0629706186570664, 0.36077194714699723, 0.18872742799145478, 0.5349409955572286, 0.14379477777422503, 0.2948966800208439, 0.8382823416326981, 0.4946517686426776, 0.327120837863124, 0.43591921578308823, 1.3665687579268502, 0.17631230170149217, 0.7624562317884868, 1.0766637101880505, 0.13239507652063823, 0.6638097993139133, 0.9796885604515414, 0.33995356270116395, 0.6669819752556102, 0.9218737734450502, 0.2978044725381843])
    SymCubatures.setweights!(cub, T[0.0009869030694832491, 0.00487076088902247, 0.0016146044262710568, 0.002442459871655366, 0.0013238366739470255, 0.0033226982261078774, 0.006546372183214952, 0.004567430694494952, 0.008392943954079782, 0.006250784237481711, 0.010036752819221235, 0.009475593456813484, 0.010676482993558048, 0.007589146491939629, 0.0032342910542169675, 0.003624870659337799, 0.0034717343512185063, 0.008128829149972157])
    cub_degree = 15
    tol = 1e-14
  elseif q <= 16
    # SymCubatures.setparams!(cub, T[0.040676809856772994, 0.9933594234293571, 0.214162059138513, 0.9595853797364323, 0.42512357224984765, 0.9103497997015512, 0.6133207466981552, 0.8379581962430264, 0.9794265714385896, 0.8997478722200393, 0.7714200686237902, 0.6331335258127481, 0.029746395996696637, 1.7947515022409164, 0.028061162821753926, 1.587126369682274, 0.024969577104493403, 1.3063824099999732, 0.1530711125876708, 0.027842118482124308, 0.9382794793201881, 0.01784707605116398, 0.33589439361816126, 0.022773216056137805, 0.8648785302352959, 0.018699236429031683, 0.5152530116647441, 0.020180307247112214, 0.7733903836489562, 0.017202724210816787, 0.12943964703419505, 1.4068114624273973, 0.12308198110051574, 1.1707095114529582, 0.31392794867774254, 0.11699998523821613, 0.8335770203486704, 0.09718370628408288, 0.4952537158755062, 0.09501363588766701, 0.7471717264846178, 0.08838366828682932, 0.25501120017276085, 0.9802122022111784, 0.4525176003746238, 0.22832242073492398, 0.7029408600731651, 0.20950601390261728, 0.36081316740942265, 1.4715721649805555, 0.14345819074134508, 0.6303775706255182, 1.2178466399964454, 0.12699217052217215, 0.5819820200548265, 1.100970954734693, 0.2959404361367043, 0.5563115214073501, 1.0515340951725114, 0.2820637067583443])
    # SymCubatures.setweights!(cub, T[0.0003449237231291464, 0.0017922220630630535, 0.0037317826419158077, 0.003943841028181184, 0.007431190363376696, 0.005393235597926561, 0.006704061398506255, 0.006130167454876543, 0.0009418626977721714, 0.004037645532983101, 0.005146004121836477, 0.005630309218313758, 0.0009595155230276578, 0.0012925466683243478, 0.0012758269653130555, 0.0020061261659943684, 0.0016783671450054372, 0.00280961843705534, 0.002197758447125428, 0.002570529986144292, 0.002060531738117513, 0.004725440482813177, 0.004679816419484507, 0.00599096798129004, 0.004471850645497315, 0.004871861095710251, 0.004380517250377522, 0.006260283814445084, 0.007132305247809295, 0.0056264781998181836, 0.0023624243767770836, 0.0024795926493089617, 0.002692242766692862, 0.005524990841468916, 0.00360658504181843])
    cub = SymCubatures.TetSymCub{T}(vertices=false,
                                  midedges=false,
                                  numfaceS21=0,
                                  numedge=0,
                                  numfaceS111=0,
                                  facecentroid =false,
                                  numS31 = 6,
                                  numS22 = 1,
                                  numS211 = 14,
                                  numS1111 = 4, 
                                  centroid = true)
    SymCubatures.setparams!(cub, T[0.045281545587205964, 0.2537904464163216, 0.9959481351143199, 0.9274532548505848, 0.6257762999594808, 0.8472186415707206, 0.9671453508185685, 0.026847241295899802, 1.7952845229842747, 0.029966113528883363, 1.5814078989333515, 0.023099401698838076, 1.303424836901987, 0.15004701952493738, 0.03390380524152307, 0.9342924509134978, 0.0017127023801264114, 0.8368259161034622, 0.007908331396807586, 0.10937882631351452, 1.4278227566200636, 0.31649671032706744, 0.05161543591193, 0.874721734267193, 0.07750207881482551, 0.4970036876492346, 0.05065481551208896, 0.7617508686532033, 0.07244962897970288, 0.22854692234632518, 1.1837750712865132, 0.444144228948042, 0.21823469992700653, 0.7333600931178887, 0.2358243330423619, 0.37095466822671086, 1.474581578477629, 0.14973521768858936, 0.6222435728859834, 1.2060370210345144, 0.13483895952285774, 0.5733642915557914, 1.0969017543221986, 0.32495705356012616, 0.558174612336641, 1.0662336598452646, 0.258627113009147])
    SymCubatures.setweights!(cub, T[0.00044135960686352563, 0.004462509356016514, 0.003464700122049525, 0.01034828978884195, 0.00939351900131461, 0.010842952623993695, 0.0020651082763411144, 0.0008181222837940155, 0.0014173509373744358, 0.0012291161984204816, 0.002584303072054426, 0.001008583409162352, 0.0017137710909476485, 0.005360744874736632, 0.006179889089930036, 0.005983453514277345, 0.005839228067957987, 0.00755208314911875, 0.006252135873092713, 0.012045025366791176, 0.008479067133012027, 0.0012654158195992348, 0.0039380175083981155, 0.0018265404066583146, 0.008080132621195509, 0.0049323204104993815])
    cub_degree = 16
    tol = 1e-14
  elseif q <= 18
    SymCubatures.setparams!(cub, T[0.03218632427874295, 0.18277834759129588, 0.37570305564040046, 0.5647728421058482, 0.7085766054616001, 0.024271614562100016, 1.8347398818402152, 0.0239470153664658, 1.662653798378602, 0.021215630179879867, 1.4264443143152872, 0.018486254825594545, 1.1383753665232237, 0.12666590483919624, 0.02382112219616801, 0.28621908896357395, 0.01986233466244978, 0.45570953414436866, 0.016528632833258484, 0.6000669570705549, 0.01544990394988216, 0.11420677083219842, 1.4854948026644648, 0.10493575994171536, 1.2854006124210504, 0.09395959070054728, 1.0437744416032544, 0.27288796018255734, 0.10387945391954902, 0.43892252103870893, 0.08648511345118752, 0.5819340432591967, 0.07566956433301017, 0.227898131588059, 1.0975429410061033, 0.2089689954468006, 0.9069693288543188, 0.4102995883623609, 0.20863652575468242, 0.5537863695579436, 0.17199377807867644, 0.34375127715192777, 0.7437405116958161, 0.5139763057507128, 0.308170027129841, 0.302497448948452, 1.5537342523975055, 0.12192692767984621, 0.5357384882233825, 1.3354897690270395, 0.10860290008924896, 0.8054112178493144, 1.0806637094306275, 0.09559315226024476, 0.5019946184882687, 1.2215929898239508, 0.2589149746221707, 0.7517448872368152, 1.0034989567679484, 0.22804479409837386, 0.6834688869043647, 0.8977251075936078, 0.40304674185298905, 0.4840102355621522, 1.174571468378619, 0.24903182274152397, 0.7263019929058256, 0.9657310710529147, 0.22072945422717807, 0.6595833058586322, 0.8697859714788556, 0.38927673288027886, 0.6215875856865059, 0.8153103124500438, 0.36944247516421147])
    SymCubatures.setweights!(cub, T[0.000172372260784541, 0.0023522720089646817, 0.004969208471540756, 0.005479651117466608, 0.0027346288970878756, 0.0005204089036787224, 0.000779542297556491, 0.0007940159278560783, 0.0006922016921975886, 0.0012076577097583255, 0.0018535991996488212, 0.0018308260400811472, 0.001403886151979098, 0.0030895580263752844, 0.0032194354869903956, 0.002879813302738201, 0.0039828459409418145, 0.003854536042095059, 0.0026975840998188097, 0.0052914381918962965, 0.004680445582917884, 0.005378948318348577, 0.0034436837218239404, 0.004565703901976264, 0.0037860053128784977, 0.0015434168343348654, 0.001565613640083025, 0.0013598106358658946, 0.001811218337204813, 0.0016738002797087588, 0.001645988527766185, 0.003929497124508638, 0.003466183275027603, 0.0034730916339655, 0.0044928452153375036])
    cub_degree = 18
    tol = 5e-14
  elseif q <= 20
    SymCubatures.setparams!(cub, T[0.026235073523158316, 0.9941783523159263, 0.15410772850631407, 0.9732668635573675, 0.32896043457408386, 0.9351927775671707, 0.5020399198509936, 0.8855842600231202, 0.6435482283389237, 0.8215038205845979, 0.9849527787033872, 0.9278721265699311, 0.8386195967994342, 0.7264382267655505, 0.6070311844223174, 0.02060721875308032, 1.861703086501912, 0.02021658597707565, 1.7164665220417223, 0.018224644848237374, 1.513101955333914, 0.01687027182514579, 1.262089625647263, 0.1079903205207026, 0.020139138471380706, 0.9548872974122032, 0.014304805524436939, 0.24422357252269133, 0.018132869410138287, 0.9055612001283477, 0.014521994174411312, 0.4001464304258677, 0.014901534536737664, 0.8392617882507785, 0.012522818974668599, 0.5475108675680717, 0.01197629753925832, 0.758216396926336, 0.01267710460507457, 0.0970986142214088, 1.5640232297116998, 0.09368148951092686, 1.3841862129011795, 0.08335083485229011, 1.1699586099263646, 0.23576460011057415, 0.09159789061834618, 0.8790743821750115, 0.0711826858052932, 0.3889221120587519, 0.07739911636975251, 0.8166608464322267, 0.06757465607312098, 0.5250418916977998, 0.0657695712097856, 0.7389070743002298, 0.06108186787626463, 0.20342475036660304, 1.197347447921784, 0.1885570880030539, 1.0236816270935085, 0.36355681695073144, 0.18839206282941356, 0.7776700948409768, 0.15685902886418251, 0.5064909377356093, 0.15110290512276883, 0.7036517263385236, 0.14220789491089447, 0.3065202187585157, 0.8755935769901346, 0.465307556531469, 0.2760667452222141, 0.6626738639361359, 0.2558569844260604, 0.2585267973643046, 1.6205209133196463, 0.10284311567343286, 0.45853872137622803, 1.4283834219164546, 0.09451962293108782, 0.6970070891224232, 1.2020985537687785, 0.085671088052477, 0.4351552662211405, 1.3232379659817521, 0.2256419186411567, 0.6580047861006134, 1.1223030702980397, 0.20336745379967733, 0.6088483929928525, 1.016404005445834, 0.36118664387398014, 0.42181382497981473, 1.2742503233695972, 0.2190284326977373, 0.6376169798609427, 1.0823740263265629, 0.19913121509428755, 0.5908362723094897, 0.9888136832430534, 0.3495690403883576, 0.5639322705827271, 0.9332395486569621, 0.3315508507274537])
    SymCubatures.setweights!(cub, T[9.531177198288232e-5, 0.0009223829890891882, 0.001399775799708288, 0.0018480256297753403, 0.003543650580037429, 0.0023881875046595034, 0.003548965369882417, 0.002776352420425414, 0.0034915173537336933, 0.0030221592683990035, 0.00041816795356553675, 0.0016495098271490938, 0.0026081308925320758, 0.0032041866743648074, 0.003150677274619523, 0.0003159768974764574, 0.0004723509795519229, 0.0005099952582302842, 0.0004957696553083157, 0.000744694841246093, 0.0007500051589230677, 0.0012845738427005717, 0.0010237603186453565, 0.0013690617767236673, 0.0010525007734989193, 0.0011130076071856174, 0.0009792295940196619, 0.001994265332602829, 0.0021473792471465715, 0.0022118119629363363, 0.0026855294286059826, 0.0019748259441236503, 0.0028924548746001267, 0.002229298411446503, 0.0023000927816112446, 0.002249938320518623, 0.0035689462818847433, 0.003266034819611738, 0.004055511641167838, 0.002836659708458956, 0.003049010766328443, 0.0027604883365003646, 0.003547314319981309, 0.0037166390603892343, 0.003541024854114663, 0.0009375702702154895, 0.0010872825591694924, 0.0009259506312575793, 0.0013094563850556253, 0.0013064617417338584, 0.0012318968252004216, 0.002829289670083863, 0.002555242220228537, 0.0025355964269501595, 0.0035947259280581048, 0.0018228055118310695])
    # cub = SymCubatures.TetSymCub{T}(vertices=false,
    #                               midedges=false,
    #                               numfaceS21=0,
    #                               numedge=0,
    #                               numfaceS111=0,
    #                               facecentroid =false,
    #                               numS31 = 3,
    #                               numS22 = 1,
    #                               numS211 = 21,
    #                               numS1111 = 10, 
    #                               centroid = true)
    # SymCubatures.setparams!(cub, T[0.0351029783124773, 0.19461341959235667, 0.32774762283388953, 0.014286075297803923, 0.012430173120742122, 1.8519938211366438, 0.019663475690596822, 1.685396536198994, 0.02108074368432553, 1.4780332697700753, 0.017354361868649065, 1.2211917929930691, 0.0957808182405573, 0.02084228060880907, 0.21140746664823737, 0.032862437239031574, 0.945533353396032, 0.02368552468444976, 0.5824912719160749, 0.020947661542596197, 0.7966930092799975, 0.027796667368427588, 0.06038600764620161, 1.662980832757256, 0.11436507870580273, 1.4532226323563582, 0.071363351076409, 1.1978332738133997, 0.44083149754337636, 0.03861996540223646, 0.5510846560134239, 0.11233473831300875, 0.16761737447335312, 1.1888291541313027, 0.30392300421827567, 0.07821522385155352, 0.3634161663158641, 0.19793504095352574, 0.7754311504920931, 0.1058409199315675, 0.21219509937551403, 0.9015319283480187, 0.41761590141259985, 0.35607655787623155, 0.6437052550117989, 0.271907380706364, 0.26790481916926856, 1.6286613248039081, 0.1004960909727349, 0.43061597011458763, 1.428002717679103, 0.11257177069831549, 0.6667810437908496, 1.2372129642465035, 0.09481023383450514, 0.40209460844668776, 1.331672541742261, 0.2630822943714447, 0.750283198295939, 1.0243348008310926, 0.2090998435372668, 0.6016333753611804, 1.0343493676488507, 0.36218578102485893, 0.5227526301706532, 1.219593174714512, 0.21306912904242833, 0.7630977638207153, 0.9779786045479847, 0.17339875979166186, 0.5737238495398432, 1.0123678238081923, 0.3279887529722974, 0.5660352953447593, 0.8699375986837853, 0.35332587145145483])
    # SymCubatures.setweights!(cub, T[0.0002116488732658184, 0.0024042166151886546, 0.005402648948389856, 0.0003314476842644244, 0.00018531468433973463, 0.00040853932211106186, 0.0007044791204848382, 0.0005151915652043191, 0.0007975140418094187, 0.002034585141909045, 0.0016642974836611783, 0.0027283010529853256, 0.0030467634495825965, 0.001174709048270329, 0.0029878347961463, 0.0027772974418439045, 0.003086099486886805, 0.006542443008389504, 0.005701155688953274, 0.004250884375574712, 0.007301836009197132, 0.004058412617176143, 0.004750250810954529, 0.0037390821029571132, 0.008012797507361722, 0.0004250642001805095, 0.001685563685392287, 0.0005583169978516564, 0.0007766771665145538, 0.0016063409808587312, 0.0008123952645587789, 0.0025970949094773826, 0.0030182374825687692, 0.003969249277439666, 0.004921736743392272, 0.012760883413153432])
    cub_degree = 20
    tol = 5e-14
  elseif q <= 21
    cub = SymCubatures.TetSymCub{T}(vertices=false,
                                  midedges=false,
                                  numfaceS21=0,
                                  numedge=0,
                                  numfaceS111=0,
                                  facecentroid =false,
                                  numS31 = 5,
                                  numS22 = 0,
                                  numS211 = 19,
                                  numS1111 = 20, 
                                  centroid = false)
    SymCubatures.setparams!(cub, T[0.1399724595324041, 0.32439084609205104, 0.4199631796638862, 0.5680758608308248, 0.7182730618031284, 0.014083116803629956, 1.9322195013973986, 0.013614171063652425, 1.829250170381165, 0.015096318561196951, 1.6441211216384866, 0.015851002566230295, 1.4010175519239694, 0.07977355385395295, 0.015891115203808148, 0.20660364239887133, 6.33364749425542e-5, 0.35039899550363734, 0.020190132533226417, 0.5043459582950012, 0.018992770812368606, 0.6243132371610652, 0.021059289878076363, 0.11071488061555305, 1.5420474083519684, 0.08740765001646521, 1.4151743827334533, 0.06293889377207833, 1.2319799304714176, 0.05383075373976252, 1.031687944168834, 0.20770206536374444, 0.04996454869895909, 0.32231013143599685, 0.8795693137398066, 0.45310418944183956, 0.16087698541555334, 0.5707360087017441, 0.12868785468102212, 0.3653239245556906, 0.57780510396999, 0.5193219192861724, 0.2704759137362186, 0.21289604035967777, 1.6817775068546361, 0.08394740712695047, 0.40141928014535677, 1.5018051153879017, 0.07981309054724227, 0.6075018460243857, 1.3074407761136686, 0.07830307719428116, 0.8396914954214879, 1.1224457231191671, 0.035486058649577795, 0.37261210857994165, 1.413709686433335, 0.19425414746117903, 0.5800810655208827, 1.2200609746889752, 0.17699080606951903, 0.8113603482863389, 1.030459672047222, 0.14162544310072245, 0.5391687995855652, 1.128839658356885, 0.3234672841496651, 0.7550998186419303, 0.947142820871834, 0.28633530696129433, 0.6924884806648465, 0.8546780991481556, 0.4524251072416156, 0.36033100502970844, 1.3077115280053104, 0.2290183049274993, 0.553067240874932, 1.162347280420692, 0.18212437959243608, 0.7714717075277714, 0.97023524039153, 0.1676846345787812, 0.49515303783612397, 1.0890823402260932, 0.3432723549283889, 0.7127992111891461, 0.9335535049485753, 0.2978272951347123, 0.6502582143292099, 0.834513673927766, 0.45941539891385524, 0.4276781937160051, 1.087179506067291, 0.2828607906218557, 0.6630990351951757, 0.8981216462526161, 0.2684300067946759, 0.6579129993629108, 0.7953180915517027, 0.3929103338143285, 0.5960352827357231, 0.7855158816438622, 0.36314903909607166])
    SymCubatures.setweights!(cub, T[0.0013239124836094134, 0.004111182249972409, 0.00157958021238787, 0.004517096541270943, 0.0030966965505332114, 6.647883321260461e-5, 0.00020241367341127653, 0.0003308499891298548, 0.00041490212599776163, 0.0004134906033302473, 0.0003076674301486554, 0.0014605237306898567, 0.0014784740975514851, 0.0010004518940844717, 0.0020023543649371985, 0.0020288673886851983, 0.0017893816488504888, 0.0012090090253824627, 0.0015569173354053588, 0.0022340966791976444, 0.004574738367089111, 0.0032497794132193554, 0.00288012346027303, 0.004539043674528827, 0.0008808902283165901, 0.0007907747821736168, 0.0003962212833009055, 0.00032573614175860656, 0.0012741822015780759, 0.0015490681271176701, 0.0012034757629581796, 0.0009140633242184705, 0.0008678650185551705, 0.0004488200303248261, 0.0031648503675318564, 0.002662618755315703, 0.0024936691393412512, 0.003081923098748991, 0.002298357351407556, 0.0026424204031094915, 0.0036900324632204934, 0.0041015471215820765, 0.0022951453620384805, 0.0021660347190994946])
    cub_degree = 21
    tol = 5e-14
  elseif q <= 22
    # SymCubatures.setparams!(cub, T[0.022568778262303177, 0.1329941580672685, 0.290725004356311, 0.45828631530939534, 0.6068101201890869, 0.7144170346780221, 0.017273929008860357, 1.8831494195861254, 0.01740995342659044, 1.7595444493769559, 0.016099118228653793, 1.5863434449768774, 0.01454462664257378, 1.368668761795122, 0.013260576208279041, 1.1187804496902032, 0.09087198126656501, 0.01720122677141029, 0.21325592820632408, 0.015321908852180056, 0.3549746870987866, 0.01303279649770051, 0.4957460866543715, 0.01214030712046733, 0.6128651969322795, 0.010013918954395601, 0.08544591699058052, 1.6188114974042678, 0.08052647230758149, 1.46383821122218, 0.07453886562841912, 1.272341830846867, 0.06636652366447772, 1.0535464331332345, 0.20549425964624504, 0.08027131645123692, 0.34544661871718746, 0.06929321962786232, 0.48215725733221115, 0.06035527533094078, 0.5992506050740232, 0.05266415834738498, 0.18088519231207656, 1.2903323484305202, 0.16915344392671916, 1.1298479228339358, 0.15182105402729132, 0.9529153782141526, 0.32773524221849404, 0.16910296280268536, 0.46284791411395443, 0.14221426858890698, 0.5770064857433621, 0.12378120673741005, 0.2814270532522427, 0.9677708214049192, 0.2598984883739099, 0.8284602795687365, 0.4350727956244424, 0.26034065642059084, 0.5499530183446922, 0.21761339396998522, 0.37197872140883376, 0.6965370168915882, 0.5143644529447843, 0.34058321626438076, 0.21966776873023316, 1.6743528409432498, 0.08962459541017036, 0.39596800038775015, 1.5060444123302696, 0.08261771100884524, 0.6079036170246478, 1.3025659522687318, 0.0751978559297125, 0.8413211852468553, 1.0778756231098787, 0.06815966928197197, 0.3796947362709243, 1.4078046657787235, 0.19808272047749162, 0.5789549118351445, 1.2266591511237326, 0.1808625400717894, 0.8013195662979654, 1.0240181523731484, 0.16245532857667125, 0.5422782272894071, 1.1213543464315374, 0.324415160012613, 0.7487203364742582, 0.9490658378624506, 0.2908743583776048, 0.6853549538302612, 0.8580231055091357, 0.4454728674439535, 0.36922766819563935, 1.3631257282428337, 0.19291426563352637, 0.563874744119591, 1.1884456775826222, 0.1772529038430123, 0.7807403231747513, 0.9970467048977982, 0.15914210226445225, 0.5285072624361372, 1.091575993056567, 0.3169819748121223, 0.7307313431113375, 0.9251196590286022, 0.2843414127124351, 0.6678630311488047, 0.8401616647253202, 0.43490761391773103, 0.5070991027458067, 1.0370992770701735, 0.30200113517185306, 0.7006303196534314, 0.8814215878896651, 0.2731715504127282, 0.6418646384481467, 0.8041648223343756, 0.41955748163261514, 0.6062441088986841, 0.7552651649958647, 0.39764038432618193])
    # SymCubatures.setweights!(cub, T[5.964904267739489e-5, 0.0009351444556801282, 0.0024180270488575197, 0.0032056159889337, 0.00287593725545261, 0.0014879673290900464, 0.00018735230700927527, 0.0002978892339787836, 0.00033806696167951633, 0.00033112322608528874, 0.00029966947025242827, 0.0004576032901517503, 0.0008372878838454209, 0.0009888739680401833, 0.0009327645602329524, 0.0006397905715206124, 0.0013082275077485053, 0.0014848974704321889, 0.0014886585999709647, 0.0013156212995660548, 0.0018331162330358201, 0.0021738327802753358, 0.0018249989699756946, 0.0014126669367032687, 0.002690299606303707, 0.0025872345971021437, 0.002342357824538052, 0.0030872871502557893, 0.002629349841781244, 0.0018166862503205664, 0.0031827011317704717, 0.002797049850364766, 0.003190936875118394, 0.002072644974036301, 0.0025549227722429066, 0.002086945752572134, 0.0006396369041351533, 0.0007071312900273752, 0.0006855042543400853, 0.0005761766300230089, 0.0009214948143410591, 0.000904588517795911, 0.0007761318288887891, 0.0009264530317322682, 0.0008337900290214272, 0.0008088106777448954, 0.0019873682975585323, 0.0019261158087567727, 0.0016601423028331551, 0.0020192665615295004, 0.0018179712944170405, 0.0017179633207217895, 0.002858715880612928, 0.00247559928667182, 0.002280214326520459, 0.0026066613626464855])
    cub = SymCubatures.TetSymCub{T}(vertices=false,
                                  midedges=false,
                                  numfaceS21=0,
                                  numedge=0,
                                  numfaceS111=0,
                                  facecentroid =false,
                                  numS31 = 4,
                                  numS22 = 0,
                                  numS211 = 24,
                                  numS1111 = 20, 
                                  centroid = false)
    SymCubatures.setparams!(cub, T[0.15118069211949475, 0.27961753074391343, 0.44321041759828844, 0.5807192446700915, 0.009080101689626977, 1.9342275312715924, 0.01641020074600556, 1.8052091603737495, 0.01667582825339645, 1.63737743848467, 0.014950517362582556, 1.413560626204907, 0.013999739273899555, 1.1361734024597363, 0.07338283325780708, 0.017323435576961352, 0.22323554658557043, 0.006996932850474797, 0.3549832123304931, 0.023440427623243297, 0.4947003736654832, 0.025708456626247166, 0.08991651820899127, 1.5623155381395448, 0.07487584446438124, 1.3916493246457755, 0.06563060590422383, 1.1900155687365068, 0.06338053607564827, 0.987817109120371, 0.2118438085315738, 0.06192767862593782, 0.17787917659763433, 1.2877475942806542, 0.14824527158728124, 1.1440612299377, 0.14165408316293623, 0.9489208723896108, 0.3346481751974456, 0.12110187815352576, 0.47837416728394827, 0.11785404666185946, 0.637360576160364, 0.03430803205225767, 0.438358628689922, 0.26382899765921747, 0.5491434576856672, 0.198427212269647, 0.3624280409928045, 0.7227878962223815, 0.5143757537162006, 0.36006120891588966, 0.20134308856456662, 1.6883591372422064, 0.09098411900012839, 0.3855450114415276, 1.516784066583109, 0.0824939183570844, 0.6204895711391556, 1.2932121648827484, 0.07477192878372936, 0.8588665339903636, 1.0590925591432816, 0.07008067856320671, 0.38804521639576384, 1.4089150057247146, 0.1940960006869073, 0.5798433888329263, 1.2316857140083526, 0.17638696231289772, 0.8056971037829067, 1.029680478735717, 0.16181130880400765, 0.543894224531589, 1.126537994327492, 0.3218752636000436, 0.745683460269364, 0.9490094223446411, 0.2923365955975051, 0.6744747469965282, 0.83204610816805, 0.48945695648216575, 0.3657436150741299, 1.3779814548672136, 0.1924227404667691, 0.5535352343600697, 1.2144940452263497, 0.1805119785878897, 0.778110947537035, 1.0152269246441876, 0.15990396621956995, 0.5273913799564095, 1.088084004602805, 0.3093606438222266, 0.7350998042156626, 0.916571874518789, 0.28161017394645316, 0.674039752246415, 0.8640967857875383, 0.4204387923213435, 0.48871209462237075, 1.0090020561432271, 0.3016587235264824, 0.6847592206404705, 0.8826140053819652, 0.2773390889521058, 0.6548112742621096, 0.7691806090263563, 0.47499810307828166, 0.6404091225704139, 0.7633870556326567, 0.3630098494693315])
    SymCubatures.setweights!(cub, T[0.0015985152471982927, 0.0031273852588997782, 0.005276560578390946, 0.0024612145869206655, 6.080608503381889e-5, 0.00025474680983506444, 0.000354689843769119, 0.000371655481039303, 0.00037464135754604637, 0.0003985335678239596, 0.0005506311939182579, 0.0018197399910844427, 0.002016713990838128, 0.001600444587104679, 0.0016159531386577847, 0.0013835091067092575, 0.0007378013268880064, 0.0016684189230159025, 0.0034019934223756703, 0.00335857228251502, 0.0025012032673452783, 0.003594354110881271, 0.0036613265976918018, 0.0011191953812811014, 0.004081709656693302, 0.0035344585640820865, 0.0029852150075828266, 0.0029318622146172362, 0.0007945382169429894, 0.0007434135246630654, 0.0006244120791516598, 0.00048483111532073725, 0.0006298004961098784, 0.0006693238714522125, 0.000372415367651363, 0.0008331835692158337, 0.0010135448299164042, 0.0005883922979078986, 0.0019316357072847899, 0.0015211843944129967, 0.0018687381056953242, 0.003035432465727825, 0.001763259738549441, 0.0015538409331515206, 0.0038708147706666895, 0.0034550828661530308, 0.0024925371718531554, 0.003042806800995812])
    cub_degree = 22
    tol = 5e-14
  elseif q <= 24
    SymCubatures.setparams!(cub, T[0.019208060242372746, 0.9962986723475374, 0.11436849611497896, 0.9789286434305156, 0.25564223248268775, 0.9518758413201647, 0.412187339986219, 0.9135342784786636, 0.5490926542651874, 0.8673307977123852, 0.6595001373052675, 0.8096779407126868, 0.9890439093827207, 0.945543755761708, 0.8770206356505909, 0.7905304190948762, 0.690914660559547, 0.5903943930371192, 0.014895121535287163, 1.8998990829136486, 0.015086585204704794, 1.7935598730322875, 0.014067903617779081, 1.6429967915078947, 0.012959540243718315, 1.4510139552620744, 0.01185605531975773, 1.2286425114705266, 0.07871596360660654, 0.014763097723873504, 0.9664975986699033, 0.010900384436104723, 0.18595996431031786, 0.012985538283105719, 0.9292939553614482, 0.010433922107649212, 0.3142845764437629, 0.011627702884637757, 0.879700211260951, 0.009549411416389698, 0.4484130758900593, 0.010933950906263994, 0.8173840874992038, 0.009378128466492898, 0.5676578100860077, 0.00905934571088567, 0.7418488609617189, 0.008174069181269915, 0.07484492899895327, 1.6686338955951598, 0.07092219725731164, 1.5314781828533364, 0.06672958675468067, 1.3596323272268063, 0.06028359788522793, 1.1615308516884375, 0.18011296312935124, 0.06938791579002591, 0.9101678723365413, 0.05309806678714871, 0.3080857279989616, 0.06268222483160156, 0.8620293346114263, 0.051170873219174996, 0.4382551887729326, 0.055393935789738113, 0.8014906808397163, 0.046729937437590394, 0.5553557118671406, 0.046525574131138094, 0.7293611131556732, 0.04337467780765945, 0.1618326820847265, 1.3686243374690965, 0.15409661629359203, 1.2199879257016735, 0.13947014243997138, 1.0584201630365304, 0.29333091514437853, 0.15290888417998974, 0.8293736234147765, 0.12074025923123033, 0.42225000673458846, 0.13048781739269402, 0.771459981348387, 0.1133819974777295, 0.5372814487102507, 0.11282325559724518, 0.7082187104602004, 0.10560401704182287, 0.2586609505541962, 1.0547338502398587, 0.23853244111801827, 0.9237036933762415, 0.3968684880161852, 0.237380931783098, 0.7349766304944383, 0.20404039617720704, 0.508666730059637, 0.19955360001325423, 0.672588299049731, 0.1873186462983291, 0.3422907473389629, 0.8013627269917379, 0.4732142200946035, 0.3099889162910724, 0.6351100800084143, 0.29280464729613725, 0.19104499336831635, 1.7168170965421934, 0.07778759048563232, 0.3457127041653812, 1.5685563652560175, 0.07233437744688058, 0.533025924700097, 1.3868975532757568, 0.06712576358138754, 0.743785885773334, 1.183843532929811, 0.06112622120790487, 0.33313301577106047, 1.4791535013294879, 0.17419127304667575, 0.5120370049979497, 1.3146076600512433, 0.16146075844903082, 0.7131615231810626, 1.1280115248957483, 0.14723037114901094, 0.48385298474466326, 1.2147650821729983, 0.2908050370051031, 0.6733367815822678, 1.0512284858656051, 0.26511889117399023, 0.6248945416277928, 0.9575922696226423, 0.40816527871498925, 0.3242409784128155, 1.4362231137939827, 0.1712778230620884, 0.4996766446293356, 1.2781008124808753, 0.1586035047049433, 0.6972115092180895, 1.0993164385874046, 0.14528887355967637, 0.4738062246996265, 1.184211065621208, 0.2858880741092558, 0.6580147206590077, 1.0269932100764294, 0.2601516470002545, 0.6119504734586226, 0.9395215376427296, 0.39942746309171334, 0.45561925163808686, 1.1300467372675402, 0.273643328999373, 0.6350926819928324, 0.9811861809248973, 0.2510641057951001, 0.5912712611945683, 0.905822342395188, 0.3852134540144929, 0.5617243365386921, 0.8554905430880942, 0.3653022388343082])
    SymCubatures.setweights!(cub, T[3.686493120359383e-5, 0.00039080348277488095, 0.0006006364170991906, 0.0009693639369705167, 0.0017411871168184032, 0.001322194491162504, 0.00247283878445424, 0.0014934515340426597, 0.0020610039773475354, 0.0016879801129693432, 0.0020599136627266084, 0.0018503934151035443, 0.0001907816947365532, 0.0007806003521172828, 0.0014166562117559852, 0.0017180756738293932, 0.0018840631612427215, 0.0017746977712151922, 0.00011950530461360321, 0.00019311889718666, 0.00022573373938866438, 0.000232963342545602, 0.00021511906053222458, 0.00029697997886553356, 0.0003747233672065986, 0.0005528086539883565, 0.0004952358668589647, 0.0007221033413848855, 0.0005413978691796412, 0.0007475294583688245, 0.0005654998535525354, 0.0005393468893146356, 0.0005162801682933058, 0.0008745359771180758, 0.0010407900333649087, 0.0010564160031263246, 0.001012595482721358, 0.0012500304750372886, 0.0009993914984682724, 0.0016200412048760363, 0.0012308010742079277, 0.0014964306554653335, 0.0012247723728819773, 0.0012261748770318831, 0.0011181923549337556, 0.0018944801724333875, 0.001928786800853871, 0.0018480044572891708, 0.0022754051309582923, 0.001609028438374714, 0.0021573291991888693, 0.001687370935147055, 0.0017116303846957935, 0.0015432353866595403, 0.0025399740995886524, 0.002196415621821548, 0.0024068276055515926, 0.0019319043251833433, 0.0019235400006542953, 0.0017872282278852479, 0.0022087587064853445, 0.0019685577814808953, 0.0022436470315930768, 0.00042710022865596344, 0.00047844347234132274, 0.0004924853851174957, 0.0004288452594639622, 0.0006806458964820225, 0.0006525042818170535, 0.000616375515338743, 0.0006891917677860919, 0.000656238112344369, 0.0006184522486700219, 0.0014182407706819686, 0.0014656048339801388, 0.0012741699545787695, 0.0015536463332504282, 0.0014401607276436261, 0.00133236676641173, 0.0022748146454793683, 0.001986836598082793, 0.0018402915753033633, 0.002279051491943338, 0.0012586706480216968])
    cub_degree = 24
    tol = 5e-14
  elseif q <= 26
    SymCubatures.setparams!(cub, T[0.01671092202861797, 0.10047410419345693, 0.22644248177177745, 0.37416812794865184, 0.5125852881231817, 0.6339179701677572, 0.7187619642647416, 0.012866873649813645, 1.9133465992441507, 0.013166240060283696, 1.8209947650077598, 0.012518620582586492, 1.689741575198336, 0.0114729685308201, 1.5218347051100773, 0.010592334251076117, 1.3227611826551007, 0.00978415373758277, 1.1037942792261088, 0.06811466281692735, 0.012932898195581915, 0.16370547466337323, 0.011740808452166917, 0.27985383321870994, 0.010565902972636394, 0.4056691857784461, 0.010380671763885834, 0.5246403354536154, 0.008450606525889809, 0.6221706426501302, 0.007784702444720195, 0.06582673033849962, 1.7082883636485429, 0.06297304473578201, 1.5860887438160256, 0.059180093622638964, 1.4335062175926996, 0.054996371496045124, 1.25326396553493, 0.049444127540877916, 1.0557687663609168, 0.15900145658083645, 0.06192173404322375, 0.27525804866098075, 0.05589770825781443, 0.3976859311287926, 0.05197072274519868, 0.5135168241188887, 0.04391730328922433, 0.6113119252357265, 0.03905872782819193, 0.14499480695794406, 1.436766424609656, 0.1372332439028884, 1.3025317746796794, 0.12775820465049972, 1.148922922093968, 0.11516956844769013, 0.9808106449603907, 0.26279544238247926, 0.1370822948431567, 0.3844673367067596, 0.12116308160281623, 0.49909901856716343, 0.10497779334387283, 0.5955065759059033, 0.09196846892049261, 0.2353380293109235, 1.142402991835322, 0.21887061593467594, 1.0182497916140427, 0.2012515523588726, 0.880781413768773, 0.3667109070965957, 0.21871011634876783, 0.47690321220199233, 0.18986804456563625, 0.5728736321918158, 0.1655605932577224, 0.32008186613279443, 0.8792596468987337, 0.29705121050247985, 0.7729128700087675, 0.4517506291010222, 0.2977679709134176, 0.5446811389072059, 0.2559218771903714, 0.3917138114640003, 0.6618112168401423, 0.5125695078266984, 0.36312261470655616, 0.16611893826975557, 1.753025785455154, 0.06825228758842365, 0.30280877987320753, 1.6208570088318377, 0.06438207903344634, 0.47095423476324577, 1.4583542496630748, 0.0593819373836258, 0.6613104591667803, 1.2730791892545739, 0.05505446044379276, 0.8661016327460088, 1.0739837733194524, 0.05052730940113493, 0.29432318130061663, 1.53907299297304, 0.15491628847630043, 0.45469059689486097, 1.390449314604267, 0.14399412570995246, 0.6378199020252409, 1.2189395661197366, 0.1331885640304643, 0.8358814189574507, 1.033608941064993, 0.12118510567738491, 0.4328738373185493, 1.2959602715457454, 0.2615292472491743, 0.6065549647170392, 1.1429257249684237, 0.24092094248566334, 0.7945252896796872, 0.977968815938125, 0.21898953647421504, 0.5688023520250333, 1.0489357607548169, 0.3734042346863617, 0.7435422650915204, 0.908679205057138, 0.33938732094240276, 0.6859676917361343, 0.8293885622167496, 0.47679808979591093, 0.2877518521007122, 1.4997615072041413, 0.15217526347195123, 0.44565680861899537, 1.3555652773956361, 0.14206618339826363, 0.6246724277745362, 1.1913238626495974, 0.1313587186809469, 0.8193338942280245, 1.0132423823102228, 0.11945184585646786, 0.4253362240370706, 1.2657762820869183, 0.25726241469460653, 0.5943736853921348, 1.118868796111545, 0.2371162388321773, 0.7791951535545645, 0.9598855141882484, 0.2153804500058896, 0.5583530112573662, 1.0283397104863075, 0.36686962918722377, 0.7306178056214777, 0.8927968880954735, 0.33340555631143937, 0.6730575248467542, 0.8161183124849726, 0.46943298699410196, 0.40961709322524803, 1.214388720940895, 0.24789587195684898, 0.5769190360330314, 1.073956269912603, 0.22945146574605135, 0.7549844945813005, 0.9251256895527286, 0.2085433318411698, 0.5408749366682704, 0.9922239971614697, 0.3553028310456484, 0.7084673132618383, 0.8624641145545486, 0.32477690827684924, 0.6533425104243651, 0.7910640569801235, 0.45642363740772307, 0.5170197825111772, 0.9410320271263966, 0.33959688993741294, 0.6762966116637747, 0.8217590897178921, 0.31114047605230805, 0.6268222039488626, 0.7578889533585297, 0.43916320056209657, 0.594154063123596, 0.7149296384470413, 0.4168400968498244])
    SymCubatures.setweights!(cub, T[2.416289785648004e-5, 0.00041110135701568246, 0.0012054310612554723, 0.0019238308668742636, 0.0020413494050353677, 0.0017158993754515272, 0.0009385471101715419, 7.717446930059795e-5, 0.00012773972603697504, 0.00015569423972777515, 0.00016125775312450676, 0.00015646543393606056, 0.00014105736754374765, 0.0001960525982134805, 0.0003901166844339673, 0.0005253203414260054, 0.0006067723619930958, 0.0004863712336805875, 0.0003497745010891599, 0.0006016654798085059, 0.0007220960527437869, 0.0007624172701727613, 0.0007394395050899027, 0.0006492142484260169, 0.0008722871792629323, 0.0011793701234469284, 0.0012203462675881376, 0.0010336351101838232, 0.0007361331500356934, 0.0013924048742764697, 0.0014648469121127355, 0.0014158115581440409, 0.0012436934862493784, 0.0016898981089483198, 0.0017029018399487186, 0.0014542535972767235, 0.0010356168965676266, 0.0020218535137176657, 0.001900654944839502, 0.0016349856207261734, 0.0020657993842711133, 0.001742814438934748, 0.0012169841400879956, 0.002086941069106874, 0.0017561921814808201, 0.001874728599178645, 0.0013071568394607523, 0.0015427674650414257, 0.0012471425915453137, 0.0002900734659105286, 0.0003379378164909399, 0.0003472172036156539, 0.0003287282055981346, 0.0002801362712866217, 0.0004702398966258614, 0.00048393152299534375, 0.00045491370414710796, 0.00039713702895503253, 0.000528529894820051, 0.000519287127610448, 0.0004430655589139014, 0.0005197032838854126, 0.00045791874553392293, 0.00041905142722732766, 0.0010271619596915168, 0.0010570592643979618, 0.000991378780130222, 0.0008584454372828779, 0.0011966650590230891, 0.0011156719000714717, 0.0010039312500114406, 0.0011290811630527058, 0.0009941089573914437, 0.0009449473673050499, 0.001765091789334425, 0.0015941499526470262, 0.0014092743620644282, 0.001591008636576582, 0.0014263190119434682, 0.0012745991579039514, 0.00189613742175059, 0.0016783251661021009, 0.0014925388529629198, 0.0016111439870962233])
    cub_degree = 26
    tol = 5e-14
  elseif q <= 28
    SymCubatures.setparams!(cub, T[0.014817468908585226, 0.9969367953606998, 0.08810517946667024, 0.9841016892736195, 0.20172258126978543, 0.962163058087366, 0.33900527924424995, 0.9335159147844759, 0.4729921502783736, 0.89436435998092, 0.5840021148344644, 0.8545850175108215, 0.6762893104165691, 0.8022063034958896, 0.9918100619140696, 0.9580050697965401, 0.9030808176265941, 0.8313829070898759, 0.7527906397489831, 0.6643768235585055, 0.5789886719285755, 0.011296973350796973, 1.9237075458784858, 0.011592381048835637, 1.8424660809776197, 0.01112721700183509, 1.7264660089819512, 0.010291805424325215, 1.577182380415096, 0.009566072116946593, 1.398800434037678, 0.009032540837164607, 1.2001295593805272, 0.05987556030986066, 0.01133003114067537, 0.9746567246190725, 0.008240472750088467, 0.14484142156806984, 0.010480257363261919, 0.9455438578655893, 0.007611409504168889, 0.25045893533251884, 0.009676967948594278, 0.9061407592527448, 0.006952007617279838, 0.3667644592530695, 0.009115988136492309, 0.8558619455340061, 0.007258747371804916, 0.4825391813543612, 0.008098285893421904, 0.7983809497262963, 0.006531105515512001, 0.5844665851138239, 0.007102664579709546, 0.7344922913990485, 0.006027846598906525, 0.05811171875608232, 1.7432072977916915, 0.05604683598001038, 1.634323353713125, 0.05326195289939965, 1.4974107375108738, 0.049650126565255, 1.3344617371968426, 0.0446254113205983, 1.1533672223337723, 0.14079681948556969, 0.05523287414134355, 0.9304553616918161, 0.04026104139281574, 0.24675983687637096, 0.05093098473492836, 0.8920245718975137, 0.03869914265434077, 0.3617329760796588, 0.0465173515542218, 0.8439528374930384, 0.03695661068532598, 0.47407370911394203, 0.04144629625716769, 0.7863101504818444, 0.035918260501776086, 0.5746253489844142, 0.03512175806185804, 0.7257474226225519, 0.033755218956322317, 0.13065571257203953, 1.4948726387446707, 0.12525231704204676, 1.3705471650314185, 0.11621271526044662, 1.2297677649031469, 0.10603793012889909, 1.074274278515551, 0.23737531533406384, 0.12404859783515501, 0.8662708811949239, 0.09520067021363314, 0.34924297397627563, 0.11087142597660041, 0.8198728986290014, 0.09107485234998533, 0.4602527006805732, 0.0978306280925192, 0.7671028340935678, 0.08577657431677063, 0.5582054293731415, 0.08424479227035253, 0.7070973654819859, 0.08093814586599765, 0.21411461537045257, 1.2193149739748943, 0.20206840214691663, 1.098950639861483, 0.18602909304752327, 0.9710977614429737, 0.33500299027714064, 0.20202006420089938, 0.7880317987080738, 0.16446256328674075, 0.442911751819494, 0.17561680592291534, 0.7367646963361237, 0.15378192077005706, 0.5372976417486796, 0.15409067491808873, 0.6824596353478888, 0.14445214646307244, 0.29776436708218046, 0.9585197302660784, 0.27599329346491075, 0.8588099246960124, 0.418263232748361, 0.2744756407922223, 0.7026950643671407, 0.24130430448846155, 0.5106753431588098, 0.2334025156629649, 0.6490635145863097, 0.22548266673441797, 0.3667547915598209, 0.7545003224185305, 0.4800701721189029, 0.33641108494371663, 0.6160723895253571, 0.32011264260181627, 0.14645265581152075, 1.7821885131876432, 0.06024778033064757, 0.26794440507190975, 1.6641146674877594, 0.05733742123425308, 0.41839476203073855, 1.518067105361466, 0.05334428811747204, 0.5906022684805868, 1.350230999299905, 0.04974296196948773, 0.7784123547427365, 1.166910322109157, 0.04640436380892695, 0.26127666671963634, 1.590000494667383, 0.13826439924156156, 0.40567886603021586, 1.4548269085170975, 0.12939350097515978, 0.5726199386033359, 1.297902252768007, 0.12068911882479592, 0.7544131303693881, 1.1255868653401973, 0.11093200901994318, 0.3889666949150796, 1.3670298197960442, 0.23525481713427376, 0.5477874616349236, 1.2242870079079362, 0.21924919235021362, 0.721727528724115, 1.068534814544433, 0.2017154649364227, 0.5177590403176787, 1.132954145489987, 0.3415335068485304, 0.6811542041731449, 0.9970323670837431, 0.3140348316060283, 0.6348762794517111, 0.9149213482406752, 0.44285899735604284, 0.2560372715008203, 1.553909113268017, 0.13591137307016468, 0.39855937699385197, 1.4212191323718153, 0.12785376897111167, 0.5623505493815578, 1.2713753826994592, 0.11935800827821992, 0.7408484645762936, 1.1046457935863898, 0.10941064387546272, 0.38321174309582173, 1.3374593275048723, 0.23252217393525648, 0.5384699321089593, 1.1998973332149607, 0.21660771436761675, 0.7081790531238072, 1.0501291944263071, 0.19901837031858743, 0.5099146341139608, 1.1118696796209715, 0.3364731375878068, 0.6696978654755924, 0.980495220178429, 0.3098406408482719, 0.6250765743244087, 0.9019611008212357, 0.43602714353697614, 0.3706706429989081, 1.2883294039206379, 0.2246606501758862, 0.523233349267284, 1.156805925682232, 0.21050436000080303, 0.6894816785641609, 1.013828041040478, 0.193506441751141, 0.49509218624863977, 1.0751608662786833, 0.32641881840200765, 0.6522718938171024, 0.9492300287449045, 0.30236294703643346, 0.6096678503765135, 0.8774805821501719, 0.4248974257059669, 0.4747345860479952, 1.0221858051261457, 0.3140130786420872, 0.6264707669242554, 0.906109068051824, 0.29107698461583814, 0.5870527942484262, 0.8436594047290688, 0.40865553547905226, 0.5575076948552236, 0.7975675254904018, 0.3886009972240311])
    SymCubatures.setweights!(cub, T[1.6736682980309203e-5, 0.0002443276777154614, 0.0002781501313528614, 0.0005985974923222581, 0.0008554533996575519, 0.0008134324834639237, 0.001475944274316531, 0.0009901154634081681, 0.0016695162312579974, 0.0010360418320048712, 0.0012611867427931049, 0.001047338009551903, 0.0011514388073584668, 0.0011713345092123416, 9.26249195021732e-5, 0.00042230943657434075, 0.000815165078779053, 0.0010503499001742874, 0.001190564654291237, 0.0012237631309001084, 0.0011960453392767176, 5.230595122834084e-5, 8.721537732636338e-5, 0.00010886969318507361, 0.00011577714177584349, 0.00011489506128718076, 0.00010936495700943008, 0.00013327615102280807, 0.00019444160687300514, 0.0002749154381233567, 0.0002559692986035435, 0.0003913811524606396, 0.000289792764319763, 0.0004568130844600298, 0.0003257620561466447, 0.00042157364934765514, 0.00030462388772144904, 0.0003225618909447832, 0.0002698699149556797, 0.0004148717395190405, 0.000510949642485861, 0.0005529597573397664, 0.0005572885057476426, 0.000509552250751695, 0.0006178529946110513, 0.0005685226657218233, 0.0008779998267877899, 0.0007004272969549056, 0.000960578233566603, 0.0007209676239300428, 0.0008890129782521174, 0.0007154362299859992, 0.0006812856012695617, 0.0006341984489938116, 0.0010326234142334126, 0.001083272076625385, 0.0011264889321579096, 0.0009686935351746869, 0.0012575844006613732, 0.000944280655420437, 0.0013784282365164058, 0.001069320650188437, 0.0012247050076120332, 0.0009320062382740746, 0.0010006510780814183, 0.0008875732463812095, 0.0015795022447274802, 0.0015615985386725772, 0.0013652871457991865, 0.0016853851164383446, 0.001157057313020604, 0.001479846768415818, 0.0011999414588612183, 0.0011766174494856941, 0.0010991958832486158, 0.0016703415311904546, 0.0014843552532419885, 0.0015074855314406413, 0.0013346693134464087, 0.0012063841742120452, 0.0012234388189317495, 0.0013587165071586801, 0.0013354103111458682, 0.0013582702830932972, 0.00020043986233557465, 0.000238741578195507, 0.0002529321447128219, 0.00024292369100340053, 0.00021009813747887616, 0.0003384702293254671, 0.0003668907993897757, 0.00033506011205455295, 0.0003271213591444401, 0.00039607933850229767, 0.000401072080604601, 0.00036674322039389595, 0.0004031505619425702, 0.00037792901257010755, 0.00035809606027620793, 0.0007494026208846461, 0.0007968714737203279, 0.0007632056520061786, 0.0006788567242781722, 0.000912957072434304, 0.0008710025598807034, 0.0008180774278940014, 0.0009137822910964895, 0.0008193938632469258, 0.0007331555785044965, 0.0013275324814024626, 0.0012796254289475894, 0.0011642782869463785, 0.0013425207861102661, 0.0011666387735778268, 0.0010767629176518965, 0.0015798277191276665, 0.0014489109166426122, 0.0012814583306617548, 0.0014682838462584534, 0.0007855852470863119])
    cub_degree = 28
    tol = 5e-14
  elseif q <= 30
    SymCubatures.setparams!(cub, T[0.012785080268792548, 0.07875303676990886, 0.18073113412651826, 0.3083524121177318, 0.4362778749293397, 0.55166483889298, 0.6533839356094001, 0.7240249406661444, 0.009948667100624187, 1.9333015298239309, 0.010270685435186808, 1.8616856069751826, 0.009932520602234292, 1.7588747536646043, 0.009272325843874432, 1.625948427072661, 0.008632999104922615, 1.465871198511623, 0.00817936693300928, 1.284775612471486, 0.007476363250213737, 1.091762958733151, 0.052885624978502126, 0.01012328642156545, 0.12899137922801918, 0.009344974809312248, 0.22478330704720617, 0.008616517160200824, 0.3328572523427233, 0.008189058464354398, 0.44312296222951975, 0.007641674532132713, 0.5451756061512869, 0.006507928904938145, 0.6288699250196296, 0.005928444607023675, 0.05173767297189735, 1.7707579880961668, 0.05023222557388265, 1.6726145932750065, 0.04794112498047123, 1.5496216613572023, 0.04508298384532212, 1.4025784445590086, 0.04202410945921873, 1.2353557325292692, 0.03857428952688367, 1.0556695763911952, 0.12563727787140241, 0.049274563506447114, 0.22204832144043304, 0.04569463963454748, 0.32867527313100603, 0.04251816672711641, 0.4374461812119311, 0.03882593338730367, 0.5369246877918415, 0.03349712751433863, 0.6204984353274967, 0.03041827175300439, 0.11821011482532891, 1.5438711725730543, 0.11285962572319197, 1.430748536709384, 0.10565533743901238, 1.3008581468630132, 0.09874953387067932, 1.1543110032366197, 0.09093712492470155, 0.995807004092092, 0.2142199356871994, 0.11243366710089957, 0.3192015973910344, 0.10220468591325528, 0.42553551604487233, 0.09187597579891375, 0.5246683866645944, 0.0805637045823593, 0.607737914542357, 0.07108588761611045, 0.19507056139954004, 1.2890140422262089, 0.18503406833168087, 1.1763239647185473, 0.1730458225679025, 1.0518531100215718, 0.15986783879530767, 0.9174998968730014, 0.3071820786976028, 0.1853280290449244, 0.41186193065858967, 0.16525002527289379, 0.5073213862333015, 0.14544920419651786, 0.5903886670282279, 0.12844793325236123, 0.27529961923768453, 1.036300464562905, 0.2572960282991571, 0.9380821362586924, 0.2397216572975956, 0.8287859139234801, 0.3920408060328657, 0.2576178216850809, 0.4867926024550492, 0.22696855672322272, 0.5664511860878081, 0.20247236725895426, 0.34741394875877635, 0.8165722603645041, 0.3234103185995005, 0.7335110646190753, 0.46138694527457014, 0.3252627124111853, 0.541291046927406, 0.2847841570218865, 0.40679169347098887, 0.6386063582221316, 0.5121137975999835, 0.37986812026579264, 0.12964811683072414, 1.80699569878314, 0.05348042052443019, 0.23791680292039055, 1.7012425224696992, 0.05130523057510023, 0.37337329440164546, 1.569454497462164, 0.048033971019859636, 0.5301733656274059, 1.4162649218786814, 0.044962248243647585, 0.7021786889776722, 1.2474571325679735, 0.042374357215349566, 0.8845583665670307, 1.0692492683293358, 0.03883803759382034, 0.23362067964947666, 1.63302539933248, 0.12384368141441807, 0.3640085578995245, 1.5104116356279418, 0.11672910653599358, 0.5156287139724045, 1.3667081807037207, 0.10945898787361881, 0.6830105036137758, 1.2069102516569588, 0.10222714184214113, 0.8607583663656444, 1.0380303864734344, 0.09379249921481218, 0.3503598021966636, 1.4285826840146636, 0.21288752279234463, 0.496058151898632, 1.2961193718664175, 0.19971388497062267, 0.6567956607479225, 1.1500930349417011, 0.18564929599132574, 0.827605841383433, 0.9952065680525938, 0.17045644605918353, 0.47187298630545793, 1.2084448445373053, 0.31228766207468206, 0.6244006943170108, 1.0786020232963571, 0.2897301379537195, 0.7862111782209102, 0.9410901058568092, 0.2660668838069627, 0.587006637250353, 0.9954654136372354, 0.41096306227724216, 0.7384052377118175, 0.8779752921105417, 0.37716044923950637, 0.6855112775747141, 0.8076435518426058, 0.5006489550016754, 0.22889781615564528, 1.5999816431559553, 0.12190325929553836, 0.35802508627667745, 1.4800014374999024, 0.11554861089086517, 0.5076054577720772, 1.3409570179812706, 0.10829717602841311, 0.6718699758371256, 1.1859959745219513, 0.10130350084533601, 0.8470640806634796, 1.0216874507104536, 0.09311754550501936, 0.346135128402168, 1.4003981768552614, 0.21052480894589234, 0.4885983584065371, 1.2722994887971233, 0.19695478623108198, 0.6457370693892552, 1.1315007480365415, 0.18360982508929452, 0.8148554915327442, 0.9809050175499356, 0.16840769483533713, 0.46533431698712163, 1.1869303176093398, 0.30854814057040164, 0.6147779310725139, 1.0619483469298248, 0.285746585407585, 0.7750103048730517, 0.9284631107600824, 0.2618869650487465, 0.5789994715569515, 0.9807719074313429, 0.4054970775538639, 0.7277880774503748, 0.8664398692986822, 0.3722901799672962, 0.6754240114351981, 0.7973628771780312, 0.4954635570261268, 0.33504660773676154, 1.3545913223925758, 0.2048696734815762, 0.4757468956588062, 1.2303240624073541, 0.19245395750929017, 0.6308621974234098, 1.0953832249394317, 0.17885230440373512, 0.7948641880452498, 0.9531034543405564, 0.16415554922784736, 0.4529750546417083, 1.1505202333509903, 0.30041570327940653, 0.6005197367921002, 1.0300352732115539, 0.27931524905572247, 0.7562740560485484, 0.9031199601673844, 0.2565118476054074, 0.5650212978566675, 0.9548689711455588, 0.3960936614641029, 0.7113806393921878, 0.8438386590493807, 0.36459025705662335, 0.6606470944279766, 0.7792046451539488, 0.4845804037965987, 0.4367858260768075, 1.0981098772997062, 0.2896430456284555, 0.5792868266449179, 0.9872554739872449, 0.2697142852216555, 0.729654299584774, 0.8681104192916653, 0.2486296801560747, 0.5462522048287126, 0.9177795751177091, 0.3830748892778252, 0.6868584851089603, 0.8146822863205236, 0.35347135473469155, 0.6396272562293617, 0.7544494787592504, 0.47035545601587603, 0.5218159994369087, 0.8704213878762441, 0.3660516128152735, 0.6571431487797371, 0.7760531389287072, 0.33957184491801273, 0.6145411151885816, 0.7221691245140242, 0.45282600177647925, 0.5835308856381565, 0.6830647247189565, 0.4298579767630695])
    SymCubatures.setweights!(cub, T[1.0854339806771791e-5, 0.00019901175418043587, 0.0006324938443672249, 0.001123026160546735, 0.0014251892401601794, 0.0013874552271409214, 0.0011063141767099188, 0.0006086867726608548, 3.56051444280322e-5, 6.044819816327548e-5, 7.70486737187541e-5, 8.382739909533573e-5, 8.469200121988763e-5, 8.287411597475585e-5, 7.228450097486774e-5, 9.322924074802029e-5, 0.0001967892765298743, 0.00028503843572983816, 0.0003479885000623679, 0.00035312456922736524, 0.00028461111001466637, 0.00020103560940413478, 0.0002959577398715932, 0.0003665997246697164, 0.00040457935982334686, 0.00041438353836056467, 0.00039918982954627454, 0.0003523959546897458, 0.00044239537860325264, 0.0006554986065046095, 0.0007583267121961636, 0.0007382745168654321, 0.0006067161104266619, 0.00042980175990421947, 0.0007621288651097889, 0.0008427849611754982, 0.0008512395287540369, 0.0007893063387198798, 0.0006890740872442787, 0.0009403873075154679, 0.0010769191576192, 0.0010462197684255224, 0.0008788030826329615, 0.0005996350458070864, 0.001247727050963307, 0.0012129978546696744, 0.001133027766188652, 0.0009997154218658521, 0.001342382203224506, 0.0012663698730199663, 0.0010440556553141726, 0.000766838562838257, 0.0014034140024301513, 0.0012962856703428553, 0.0011622202393704753, 0.0013953192494538362, 0.0011683529282775279, 0.000829153114472294, 0.0013320794880038905, 0.001134663055713877, 0.0012032024316774702, 0.0008202669720679585, 0.0009545580772300045, 0.0008064925331662638, 0.00014088683685228303, 0.00017210503946348533, 0.0001855273685653018, 0.00018365049923653478, 0.0001714263326734769, 0.00015074474033293265, 0.00024921845906557364, 0.00026478276553287724, 0.00026026201493612356, 0.0002471307110556189, 0.0002213837767280252, 0.0003066683987258996, 0.0003195949493149044, 0.0002917989228053216, 0.0002514924839828565, 0.00033021866918324643, 0.00031550268856900406, 0.00027313701640610733, 0.00029940561897524875, 0.00027206823631505936, 0.0002490852784754977, 0.000550864312927374, 0.0005932790322808154, 0.0005793980054601604, 0.0005438158869603324, 0.00047898472116059403, 0.0006888413063685563, 0.0006942885041668063, 0.0006441028191666998, 0.0005719103029977951, 0.0007487838307874767, 0.0006814035627427715, 0.0006084450458229345, 0.0006649606214638517, 0.000606905118832717, 0.0005406469408755948, 0.0010196039031304567, 0.0010122658673700276, 0.0009336164669749449, 0.0008218713542548457, 0.0010855779958526266, 0.0009856371559317163, 0.0008823471127852474, 0.0009565465893505764, 0.0008488467414885513, 0.0007658699303332844, 0.0013179405443867393, 0.0011962786206808813, 0.0010513469206907369, 0.001162109512831111, 0.0010382829094835932, 0.0009191259969764408, 0.0013294272682142093, 0.0011314910993298324, 0.001024067371306568, 0.0010792413056873742])
    cub_degree = 30
    tol = 5e-14
  elseif q <= 32
    SymCubatures.setparams!(cub, T[0.011171979148375082, 0.9975001179165643, 0.07037881861920747, 0.9874661181652318, 0.16251014155330526, 0.9705710547703618, 0.28140232907785795, 0.9457498865499001, 0.40403329518395276, 0.916355309010942, 0.5167658736835482, 0.8804353893315977, 0.6092417284763598, 0.8437925229193075, 0.6880251171488341, 0.7959842484206325, 0.9933657317199155, 0.9668323398131057, 0.9226861905884634, 0.862829265317449, 0.7939124598226615, 0.7212234767320032, 0.6447041194061423, 0.5694654570144637, 0.008841916263185989, 1.9410786865652438, 0.009182960171947136, 1.877398155826167, 0.008913887442837642, 1.7856399055837786, 0.008386961123512105, 1.6661627797326768, 0.007858603067709654, 1.5212542874861366, 0.0074631989875815945, 1.3560644798180337, 0.006949490785066169, 1.1782047933428075, 0.0471148314249694, 0.009031490772227622, 0.9798383304527233, 0.006320281210997118, 0.1157074087158498, 0.00840913435607211, 0.9567383798965058, 0.006124144008928092, 0.20272998417094587, 0.008066603045147024, 0.9252738314845996, 0.0057919325764713234, 0.30257991507430443, 0.007447105425353105, 0.8848467431667761, 0.006017397584415245, 0.40742527491109826, 0.006866367572053721, 0.8364188739412001, 0.005639179834341451, 0.5076611575665806, 0.0061045580620398155, 0.7834465606188031, 0.005714403148129404, 0.5975876403412514, 0.005102557017248732, 0.730452784619378, 0.005605707450991512, 0.046298600883570476, 1.7948278238397102, 0.04517401899279016, 1.7064084698217643, 0.043412169390546415, 1.595440666473308, 0.040990085628022525, 1.4622321613908535, 0.03818079354268722, 1.3095780788405242, 0.03548748836747063, 1.14295392918848, 0.11258491144652304, 0.044318382707624504, 0.9444808527406833, 0.03182248855781652, 0.200706501950974, 0.04232520981971399, 0.9134849612449651, 0.03090470290928404, 0.2997738947882057, 0.03858164905114515, 0.874355999791093, 0.0301808142475201, 0.4030026439034952, 0.03534984162996726, 0.8274820569972196, 0.029483558827623452, 0.5013021846931901, 0.03150437488990831, 0.776394971609795, 0.02802198748500909, 0.5875103656051429, 0.02765450472409815, 0.7207236381278224, 0.026458809204198377, 0.10684235149053598, 1.5884170618014188, 0.1028835215479602, 1.4838222375546126, 0.09685990456323693, 1.3629108529418132, 0.09042850985394527, 1.2270293553939535, 0.08438488086372985, 1.0787902929661104, 0.19438182846037222, 0.10313638461935329, 0.8922436505512208, 0.07555462257352026, 0.2912720944132182, 0.09340523693852422, 0.8546145534852153, 0.0728038620509571, 0.3929028972188648, 0.08465332057115975, 0.8098984469488562, 0.06997079945331564, 0.49031025859645255, 0.07603303037825378, 0.7593016831114754, 0.06651148497328642, 0.5763963748720224, 0.0661991274821913, 0.7063730434089878, 0.0632941180780893, 0.1783832326662177, 1.3509543049229653, 0.17116243316999058, 1.2441202994487057, 0.16013316287968726, 1.1274797809539143, 0.14925452757655097, 0.998772742234182, 0.281462349823799, 0.1710184712604719, 0.8273151045738752, 0.1337318743046608, 0.3808695366300461, 0.15315874578950486, 0.7853880461891495, 0.12785741571643086, 0.4743180969480624, 0.136870783289815, 0.737883264540954, 0.12123282979348649, 0.5580822164489296, 0.11959705858936358, 0.6859585392973552, 0.11483529732785323, 0.25528903643565687, 1.105692787128379, 0.24130149758521194, 1.0098204445464882, 0.22439898176924414, 0.908074416913776, 0.36408288416309603, 0.24101603755425585, 0.7546488446031333, 0.2004918563688193, 0.45674939302452583, 0.2115860520467363, 0.7089281036543539, 0.18915900534927557, 0.5356641140566203, 0.1887003258631458, 0.6613418401923223, 0.17963706811173683, 0.32658578384792575, 0.8883642026135896, 0.30587232517650337, 0.808301402131929, 0.4322807725024418, 0.30313831631212496, 0.6776675212301408, 0.2716180302661259, 0.510911195897636, 0.2635299539528834, 0.6315492448036338, 0.2553858910188466, 0.3852088703170925, 0.7173896452697384, 0.483824640449393, 0.35624054169919167, 0.6014946657958228, 0.34048157776448007, 0.1157316013362656, 1.8275403328655337, 0.047892146025045865, 0.21296389676064384, 1.7323813951994769, 0.04608873844146323, 0.3353639407313408, 1.6129047754428858, 0.04346957225412439, 0.4779678919051018, 1.4732885316208588, 0.04094010880286005, 0.6358322687451124, 1.3182534617592654, 0.038719981040708014, 0.8048767950900017, 1.1523501053006937, 0.0360237476312784, 0.20994450058928205, 1.6701440258751996, 0.11140263202068028, 0.327858519631918, 1.5582721237382269, 0.10570206324887349, 0.4662752630162992, 1.4265034575201347, 0.09974911687603145, 0.6205447101760324, 1.2786001694929416, 0.0937221820948923, 0.7855008829402438, 1.1207007754649627, 0.08690744144706579, 0.31711845585634574, 1.4824794300105997, 0.1929989440259221, 0.4505869856031221, 1.3596503040600358, 0.1823592493591714, 0.5991866701069867, 1.2230146430377977, 0.17082980565632333, 0.7588336849567789, 1.0768046870625751, 0.1583314071751655, 0.43083944594681955, 1.2764256390013726, 0.2860713715734354, 0.573038908827816, 1.152664925710153, 0.26757777797596816, 0.7247902844790477, 1.0208846444030792, 0.24812880683649832, 0.5421915942903393, 1.070642159463173, 0.38079262529087715, 0.6853449183233726, 0.9555940047669238, 0.3530707831815261, 0.6420216058416361, 0.882734744153185, 0.46993842491769333, 0.20577781290880848, 1.6400265685279671, 0.10995387052277876, 0.3230829038262045, 1.5295661117218786, 0.10479561687385236, 0.4598848967833175, 1.401958239818815, 0.09867450588884513, 0.6115402758028056, 1.2586572832092642, 0.0926493834338205, 0.7745545368026204, 1.1040849267037613, 0.08599782492973138, 0.3139222290575911, 1.4555777135680323, 0.19133501073387163, 0.4451354701073036, 1.3362532734015764, 0.18003563650831528, 0.5903942927434582, 1.204208785005297, 0.16916010136828707, 0.7473223391052987, 1.063037363823344, 0.15672482577212957, 0.4259805924671556, 1.255411821296135, 0.2829237014005598, 0.5651013103856433, 1.1355263033608856, 0.26456244070483437, 0.7147278641628081, 1.0074918486887188, 0.24541427070061408, 0.5356532024914514, 1.0551988421878373, 0.37647554021439256, 0.6760815507338023, 0.9432422057573747, 0.3494225275016532, 0.6340618683467946, 0.8729188472696695, 0.4642470168395634, 0.30517436232119155, 1.4119268993511793, 0.18643458942385624, 0.4340355880590535, 1.2958668754759113, 0.1765784542696571, 0.5778101888165887, 1.1697734076930486, 0.1650938488748741, 0.731163454188469, 1.0337236998822315, 0.15327054520493316, 0.4151418749366107, 1.2201862746212677, 0.27637348422233066, 0.5530381974036154, 1.103921485284279, 0.2591909743042126, 0.6996030589521259, 0.9812550029503094, 0.24059728116857818, 0.5238122470787389, 1.0292428933346174, 0.36882763833622806, 0.6629188509025681, 0.9201358035463248, 0.34223040445343206, 0.6214690034902431, 0.85439058267092, 0.4548276846464696, 0.4020516965163155, 1.1679945845405388, 0.2675362007660664, 0.5351683334629871, 1.060867012793136, 0.251234208194368, 0.6774481195238984, 0.9457045091048271, 0.2334632764289877, 0.5080830708876876, 0.990653237946924, 0.35774446116474007, 0.6427603559202858, 0.889846184565145, 0.3335661060388491, 0.6037530220787382, 0.8296956528421378, 0.44212935003875736, 0.48621246261944023, 0.9428660510637442, 0.3432590828350156, 0.6175487725970539, 0.8500595143808525, 0.3205073461914499, 0.5812603715247506, 0.7967012239951841, 0.4266255722622919, 0.5538132039669326, 0.7555094808187145, 0.40540645384384677])
    SymCubatures.setweights!(cub, T[7.301544772195012e-6, 0.00016756237046129741, 0.0001428066556718224, 0.000348585211835348, 0.00046162232598915035, 0.0005188496395923017, 0.000839390427727145, 0.0006294661424594369, 0.0011761109112851363, 0.0007272944385193228, 0.0011660528363873377, 0.0007247911653444256, 0.0008345233703874511, 0.0007149052354400774, 0.0007541514666577747, 0.0007437180793509171, 5.328153196871179e-5, 0.0002479905706053417, 0.0004956882253005511, 0.0007088197525961302, 0.000803936245336496, 0.0008372049405252288, 0.0007949050455869735, 0.0008367699238196625, 2.4910914856114476e-5, 4.300632769159766e-5, 5.553270499350718e-5, 6.18045877570645e-5, 6.376247051277095e-5, 6.326346282432454e-5, 5.777539513084398e-5, 6.622000745954439e-5, 0.00010604165427229918, 0.00014321784478815272, 0.0001456091292010134, 0.00021958576412277664, 0.00017270974203069197, 0.0002672647148543085, 0.0002052925834363985, 0.0002798927140872489, 0.00020861767001583844, 0.00025041273966400923, 0.00018890085849596107, 0.00019399250129704807, 0.00016717501304950203, 0.00021277220204266568, 0.0002668921568427411, 0.00029978693747257317, 0.0003125252730920517, 0.00030689278515700675, 0.00028051219392290997, 0.00032160127326612384, 0.00032466920716668586, 0.000497640864232404, 0.0004014977887389896, 0.0005870049975460098, 0.00043274426284660683, 0.0006053645769794869, 0.00045182228009678366, 0.000537769891715963, 0.0004014515256698585, 0.00042789168150005584, 0.0003688544142864881, 0.000573324643342978, 0.0006370825432736306, 0.0006659878543269514, 0.0006290922853702161, 0.0005590829377334518, 0.0007168857840285991, 0.0005828730627061555, 0.0008529304476480269, 0.0006450012915399674, 0.0008659002930855217, 0.0006360338864380253, 0.000780989370916166, 0.0006087069447767129, 0.0005946083463528612, 0.0005505186503637064, 0.0009584253799353781, 0.0009683780461251142, 0.0009510661594219184, 0.0008372238722461425, 0.0010868964471464752, 0.0007567788353957863, 0.0010629403722911001, 0.0008202375718380654, 0.0008914116187162245, 0.000747111622390203, 0.0007671472180674464, 0.0006849626026824631, 0.0011422555823176594, 0.001104599513329289, 0.0010011124467466175, 0.001194508558326213, 0.000854097509857298, 0.000997701765769668, 0.0008628841015712929, 0.0007976800376393809, 0.0008007484340951913, 0.0010820052302534626, 0.0010300522967275523, 0.0010204421893642626, 0.0009027470359415924, 0.0008229036916541655, 0.0008343791754441236, 0.0009072907203260115, 0.0008702789500688721, 0.0009264092135400851, 0.00010124399700938557, 0.00012528629658881908, 0.0001377979496045256, 0.00013871773983874735, 0.00013091313979183915, 0.00011966391405012892, 0.00018184377188353997, 0.00020107543695895637, 0.00019933558086096745, 0.00019196575512427796, 0.0001788520892966476, 0.00023205608954754352, 0.00024746495318969317, 0.00023588924649738409, 0.00020105903357749217, 0.0002580447963317067, 0.000257665973873819, 0.00022942425632259728, 0.0002581488175082885, 0.00023159645938308502, 0.00020606085819471333, 0.00040872032263425314, 0.00044898234313088433, 0.00044833394527227206, 0.00042810577968111445, 0.00038867497044226706, 0.0005303164241715435, 0.0005453642690005615, 0.0005141172526007905, 0.00047659684470315896, 0.0005952398835804661, 0.0005607784309567437, 0.0005133548164590458, 0.0005518147784874352, 0.0005107302175363627, 0.00046932063675491037, 0.0007818621517329088, 0.0008010951664065421, 0.0007562448569389358, 0.0006955664650388605, 0.0008808275783092766, 0.0008236099992586066, 0.0007446566286759701, 0.0008097422395584106, 0.000714418621310439, 0.0006510928149017122, 0.001074755611463994, 0.0010159126073422758, 0.0009126532421695471, 0.0010202257030080416, 0.0008708633823200213, 0.0007942131117617588, 0.001133884190972759, 0.0010515934916763353, 0.0009035404213751143, 0.0009855138120358153, 0.0005127967699135106])
    cub_degree = 32
    tol = 5e-14
  elseif q <= 34
    SymCubatures.setparams!(cub, T[0.010630914833370465, 0.06279282162636195, 0.14726174216515345, 0.25439692945379117, 0.37253400616459464, 0.4830704394456631, 0.5807537053641586, 0.6669803507898387, 0.7256240414734728, 0.007946352664858568, 1.945917498317762, 0.008186202422950591, 1.8883631891733639, 0.008032499483850267, 1.8056111775582049, 0.007622634269473475, 1.6983889820744216, 0.007144400190977871, 1.5681566467460917, 0.0068414678724712805, 1.4180089467141186, 0.00645813810308232, 1.2541112276591115, 0.005789314682918124, 1.0822619578024153, 0.04215444162922966, 0.008056396729347303, 0.10389362536644475, 0.007621146108473452, 0.18362048812600248, 0.007210467065768172, 0.2763894102940846, 0.006954888598832471, 0.3749553017186231, 0.00644631090595126, 0.4720255111704739, 0.005697476887750893, 0.5605285645208578, 0.004935343146387225, 0.6327581594093227, 0.004725186499043654, 0.04173315338324623, 1.815903655366778, 0.040945281694246674, 1.7358862467272522, 0.03939572303935757, 1.6354228091619487, 0.03743381627332282, 1.5142360010640348, 0.035541816923823466, 1.3730707581281003, 0.03303740469549388, 1.2185117510768393, 0.030564285106489272, 1.0547419643379758, 0.10157749927651845, 0.040114592665272855, 0.18205384221917287, 0.037761133015647284, 0.27375403494492095, 0.03574858314195864, 0.37128078546469884, 0.03285883224557227, 0.4671186469273159, 0.02959247356745385, 0.5544263523644043, 0.026299004395406347, 0.6274818261796976, 0.02392709650624502, 0.09718902789686595, 1.6257538268201916, 0.0936577822981032, 1.529997669303864, 0.08860023249413546, 1.4190680393975197, 0.08430511472452304, 1.2918447748992974, 0.07829842142848126, 1.1545662520887388, 0.07362662975861026, 1.0057366736803184, 0.17696528950274373, 0.09227294359280093, 0.26756454178100636, 0.0859916257716117, 0.36333352937495444, 0.07875399960121125, 0.4578338845920754, 0.07132602555928182, 0.5450033199089986, 0.06338836508225688, 0.6168596928815062, 0.05701557678368321, 0.1633986711121827, 1.4067583294521016, 0.15782503924665597, 1.3045103799776039, 0.14872773061883823, 1.1932765018362506, 0.13957420708740884, 1.0711574011810214, 0.12943900248136803, 0.9423497094679552, 0.25814479549663405, 0.15759521836895984, 0.3526171155182269, 0.14369599177683975, 0.44460687898204915, 0.12963473453630556, 0.5295782605444481, 0.11510582299947503, 0.6016906283990117, 0.10301516978087685, 0.2363910343916465, 1.1728333567863038, 0.22541984095994008, 1.0781224092076012, 0.21159271196319415, 0.9778968366743351, 0.19683444805427722, 0.8686381124417739, 0.3392398852238519, 0.22490186489146619, 0.4297458783805426, 0.2022962052233736, 0.5130758738665798, 0.18035334585155724, 0.5828608134725651, 0.16142903776005468, 0.30580882293603856, 0.9563225458108123, 0.28786431340838337, 0.8746817918451977, 0.27050326505195016, 0.7863731133422143, 0.40988908640661303, 0.2878036584344298, 0.4921644683981602, 0.25740265613981755, 0.5616308589043102, 0.23169693591431412, 0.36768015288089223, 0.7706168332217536, 0.3450540570098977, 0.7019358209284271, 0.46888368111276585, 0.34730175957984355, 0.5379662258621606, 0.30763051440347905, 0.41788841383761904, 0.6202002326495251, 0.5115833477465637, 0.39425897620715655, 0.1035667984423451, 1.8456481830659057, 0.042827482497160035, 0.1909727206107468, 1.7595391396825815, 0.041702105879357694, 0.3019068542086988, 1.6511051056757529, 0.03952363946386999, 0.4324525492803771, 1.5230740234003195, 0.037321913148078756, 0.5779774016702697, 1.3797521578552616, 0.035540669940478775, 0.7346493989416144, 1.2256797556122838, 0.033438773459834786, 0.8990219213907884, 1.0646703238199096, 0.030461008381224133, 0.18936121091465433, 1.7020794947873485, 0.1008549539848424, 0.2967888296448669, 1.5996697796728447, 0.09614036768473207, 0.4232949377013703, 1.4787779247459498, 0.09113083976785408, 0.5652003270135956, 1.341893552952335, 0.08619760443111646, 0.7186447663708895, 1.1944292655818514, 0.08074298678212732, 0.8794316510180981, 1.0401013246163917, 0.07448100934153695, 0.2877495040783061, 1.5295753321077947, 0.17597883647740187, 0.4103450411010091, 1.415689217948483, 0.16686611953901992, 0.5480024258650634, 1.2882750789304795, 0.1574694844771642, 0.6966575947238877, 1.150281328889438, 0.14715442613099533, 0.8524406918110612, 1.0057171097975628, 0.13630498596956814, 0.3945136237680319, 1.3368600520104363, 0.26244311093923073, 0.5265456573017891, 1.2201191531399411, 0.24712950493101973, 0.6688198259120561, 1.0943062853966055, 0.23111030040178632, 0.8181696063793181, 0.9625401912289383, 0.21381690167101516, 0.5010706450850274, 1.1401416212456155, 0.3529612000331241, 0.6364463696291819, 1.0284451032274586, 0.3295359502884257, 0.7781192129350819, 0.9114596426785965, 0.30506390042821263, 0.6002186033964527, 0.95415580540738, 0.44047230868974907, 0.7331924453330246, 0.8538834304626005, 0.40764831882686, 0.6849292356922357, 0.7911263148974008, 0.5191303498099614, 0.18581154642401948, 1.6744641613792988, 0.09953156367301078, 0.29234576896820386, 1.573651958161726, 0.09540804009452969, 0.41779865266554955, 1.4559989357606102, 0.09015007230524526, 0.5575362972436687, 1.3219262909858187, 0.08574120205570811, 0.7090608754169798, 1.1783232036986406, 0.08037347082735348, 0.8683142661485717, 1.0268364010066409, 0.07396581692893245, 0.2854984029170972, 1.5044023704568221, 0.17438272100399949, 0.40600266956821857, 1.3927318125828323, 0.16485514559697775, 0.5409774156766906, 1.269882418751564, 0.15593169028353493, 0.6867703130785249, 1.1354643414711831, 0.14661919960582834, 0.8416367399361911, 0.9941871197994088, 0.13500116185515365, 0.3905075748150885, 1.3164248775391263, 0.2598027265048251, 0.5202203326455179, 1.2031490605646487, 0.24443987428014924, 0.6601608756890159, 1.0810328004517364, 0.22871056425810066, 0.8083672189000539, 0.9517178235619304, 0.21152583118510515, 0.4956509511057422, 1.1243586195215824, 0.34934078840366894, 0.6285564393362572, 1.0163505493376797, 0.32587410690169244, 0.7689676587045128, 0.9018297670884005, 0.30138294840338764, 0.5937403586688487, 0.9432738469713708, 0.43594061443907367, 0.724973221964899, 0.8449454454948297, 0.4035304354927568, 0.6772150018322796, 0.7826388151491237, 0.5150335340958019, 0.2777754368567904, 1.463100408733662, 0.1708182754189528, 0.39652584053980383, 1.3547347857398275, 0.16189944236368767, 0.5301853288830303, 1.235351797130856, 0.15325608522910034, 0.6738971523427418, 1.1069975338929106, 0.14288200920645364, 0.824801002344124, 0.971724862197686, 0.13246170509283234, 0.381609294433819, 1.2818729246132134, 0.2544982448084608, 0.5096602144434702, 1.172197866271541, 0.2398382802451839, 0.6477811112535676, 1.0544274675338754, 0.22479413799343967, 0.7927374910151018, 0.9312996515363308, 0.20730582576269274, 0.4854215579418377, 1.098111418366368, 0.34225446751283384, 0.6167896022836449, 0.9929644156935815, 0.319756161979387, 0.7545056309844849, 0.8825562349795851, 0.29666904709156217, 0.5824364912308989, 0.9237107213758706, 0.42762445306003993, 0.7117665867379609, 0.8283493920937526, 0.39725339386710373, 0.6653808428376791, 0.7685935313083311, 0.5060216995304072, 0.3701255117648392, 1.2328359138692702, 0.247104285799722, 0.49455980426449225, 1.1287075066935661, 0.23307471333466903, 0.629429694305174, 1.018241523360308, 0.21833713773353974, 0.7700940088104853, 0.902001684901299, 0.202503657824852, 0.47215802704390225, 1.0595009966309727, 0.3331362306585573, 0.6000654961367738, 0.9606328325812777, 0.31190133312211565, 0.7341399434338831, 0.8568520441261314, 0.2897802604149922, 0.5668721795092739, 0.8961236152961551, 0.4167327620471011, 0.69248726149588, 0.8054613264308983, 0.38789574889935874, 0.6497731771579284, 0.7487321383102836, 0.49374417006332505, 0.45448659564925303, 1.0112748774287412, 0.32088234383717257, 0.578559539003431, 0.9206228619736098, 0.3007410598594452, 0.7074733005832956, 0.8239837853904699, 0.28046133825701686, 0.5476686556171984, 0.8611361622116516, 0.4032881942049661, 0.6687396176263709, 0.776099438149572, 0.3761726344041806, 0.6284471337658212, 0.7238342890256956, 0.4790077247929883, 0.5236737730739855, 0.817951954220035, 0.38617999076451653, 0.6414802274980228, 0.7412682775390101, 0.36096286121532645, 0.6034027420599962, 0.6942223460673612, 0.46125302833914633, 0.5751399438850172, 0.6593247588812284, 0.4395283614839984])
    SymCubatures.setweights!(cub, T[6.127199109054335e-6, 0.0001016765784919534, 0.0003474644980466799, 0.0006609428891262158, 0.0009329842649459088, 0.0010443765474662944, 0.0009819536324209353, 0.0007418471646045045, 0.00041876512888266656, 1.8329621655703924e-5, 3.086971016703621e-5, 4.052223877462063e-5, 4.573129632094533e-5, 4.7635438695949956e-5, 4.86567492495152e-5, 4.622465164642692e-5, 3.8757082179924014e-5, 4.723303561110897e-5, 0.00010602245706121389, 0.00016172486741706678, 0.00021053581783753712, 0.00022827272838096585, 0.00021151925819953157, 0.00017019461363205045, 0.00012630173974559821, 0.0001555751433839619, 0.00019704110729313153, 0.0002252787913530077, 0.000237064502847571, 0.00024333429315308227, 0.0002284845716149582, 0.0001997860684587871, 0.0002375165996929673, 0.00036982729471147325, 0.00045957543970257447, 0.0004912890973698017, 0.0004655069473217036, 0.00038666035680217465, 0.0002631855714702198, 0.0004344137792102604, 0.0004906663962165678, 0.0005124984034821098, 0.0005004395753219363, 0.0004655641707414013, 0.00041011940541568265, 0.0005371350549494959, 0.0006676474098597133, 0.0007174827045270913, 0.0006779158762896742, 0.0005511854316686129, 0.00038742620167706937, 0.0007550887576145525, 0.0007712942690556769, 0.0007620130716339133, 0.0007042177576110891, 0.0006225807160321661, 0.0008621580283102892, 0.0008920547335408057, 0.0008253942997131022, 0.0007011647258173893, 0.0004885691072939845, 0.0009655206818021566, 0.0009005236745203429, 0.0008596655883195234, 0.0007814069012637062, 0.0010017066538647683, 0.0009160877262094111, 0.0007572033075578126, 0.0005747729270239625, 0.0009906528370978577, 0.0009178910557678994, 0.0008337307186114322, 0.0009599673117797507, 0.0008062430830523589, 0.0005717975526104478, 0.0008999526499993531, 0.0007722657848832455, 0.0008222241146797923, 0.00056468260132914, 0.0006282224214673961, 0.0005448200022249712, 7.333535802143412e-5, 9.296803205066186e-5, 0.00010326662695690406, 0.00010679755848453821, 0.0001039852346891078, 9.579589469036835e-5, 8.535472546556689e-5, 0.00013589199633940192, 0.00015173012707193906, 0.00015256050694965404, 0.00015430110362838587, 0.00014047049267000035, 0.00013048558126607565, 0.00017668895541662986, 0.00020085738140119983, 0.0001840367315325803, 0.00017095418760206603, 0.0001523721259722305, 0.0002060037812932944, 0.00020807388332371663, 0.00019025012017327516, 0.000170875095875121, 0.00021285984071504676, 0.0001960428890160109, 0.00017800017245593882, 0.0001848773844338834, 0.0001739019341816745, 0.00015133437535840447, 0.0003053652857088765, 0.0003403430119462343, 0.000347244445288992, 0.0003396700024529003, 0.0003178049463185055, 0.00028347711744935696, 0.0004085424914247499, 0.00043181505397080834, 0.0004143810593606186, 0.00038495008841115, 0.00034683295733370557, 0.000479768703666609, 0.000458522831877361, 0.00041825232178708887, 0.00037823506456873744, 0.000475025532655496, 0.0004380182680678698, 0.00039394140767281014, 0.00041466459392708517, 0.0003704497097247391, 0.0003392317744954032, 0.0006102357768170724, 0.0006291087503997503, 0.0006132897214921347, 0.0005735645199022745, 0.0004937427918330546, 0.000701929117640445, 0.0006828139770870416, 0.0006218085358518492, 0.000555646034420752, 0.0006903429388580195, 0.0006254242518362369, 0.0005569723553654271, 0.000619075384255034, 0.0005349645257777683, 0.0004761793539046686, 0.0008636132336795996, 0.0008482649354398481, 0.0007792765100484852, 0.0006804459624129158, 0.000860676516224544, 0.0007727005737475672, 0.0006917089278516673, 0.0007296716130683827, 0.0006732127597832311, 0.0005874961873413039, 0.0009728093365564699, 0.000878213645635733, 0.0007967994657069574, 0.0008550779188205088, 0.0007495703141667252, 0.0006596708615631706, 0.0009386614080713931, 0.0007830047426098054, 0.0007155318119846747, 0.0007182400988418958])
    cub_degree = 34
    tol = 5e-14
  elseif q <= 36
    SymCubatures.setparams!(cub, T[0.009070169900198951, 0.9981273091263114, 0.056698262497157356, 0.9899702957328278, 0.13361477351352852, 0.9756770836404604, 0.23469318005021833, 0.9564267770526282, 0.34425240059920187, 0.9310492640346771, 0.45307447231788395, 0.90223365274455, 0.5489326761857459, 0.8687634330983043, 0.6284957382693334, 0.8337258102906036, 0.6958061673804251, 0.7914648377155905, 0.9946296963010024, 0.9730513429169267, 0.9372398731999511, 0.8879117355315589, 0.8278717093872527, 0.7623825717429314, 0.6958538176340734, 0.6292662580403661, 0.5614425343260985, 0.007132780940094448, 1.9524148896610842, 0.007421635292314713, 1.9010486963593325, 0.007262756783328213, 1.8267093736872773, 0.006925669997472389, 1.7292389301921318, 0.0065638923908188365, 1.6099620550791252, 0.006263472131637265, 1.471879855037128, 0.00596295397653971, 1.3200459641726074, 0.005519860656415913, 1.159738530531422, 0.0380452684528002, 0.007265707798212239, 0.9836580431205008, 0.005144902147177505, 0.09419985520229782, 0.006863934344990947, 0.9647437252420898, 0.005164271679979169, 0.16693617901803828, 0.006688262304201151, 0.9391286261057958, 0.004924565444646091, 0.25274136544337344, 0.006198821252672051, 0.906508103899654, 0.004775431594881229, 0.3455257210155083, 0.005820724243281995, 0.8660675998256187, 0.004471059492645211, 0.4389076055093046, 0.00541402153494503, 0.8194928040379961, 0.004549868916166495, 0.527036859796575, 0.004764491161080851, 0.7731899508384353, 0.004397864592222994, 0.6039299673612517, 0.004473299433052728, 0.7247529872757299, 0.0038253474142360623, 0.03783228091010142, 1.8333580328796435, 0.03716104805979778, 1.7607779915890198, 0.03597886659108545, 1.6692199929397187, 0.034347971441492375, 1.5585153481205944, 0.03249723169943833, 1.4298053339849859, 0.030400705009246097, 1.287013954511315, 0.028452196507243146, 1.1337169465502053, 0.0921015564800844, 0.03626880397797753, 0.954832946183074, 0.02613286220334582, 0.1657803314715076, 0.03512447261683057, 0.9292591089720249, 0.02575193747826697, 0.2509207291607992, 0.03243320932079416, 0.8967748014839876, 0.024867581702911007, 0.34268920340232595, 0.030223357219035137, 0.8578345773132035, 0.024161300559829724, 0.43505605029501887, 0.027724181557673717, 0.813920312488662, 0.02327897057553928, 0.5217032760833856, 0.02498756589118742, 0.7657654209397985, 0.02275584778810012, 0.5995638609581896, 0.02193508485560876, 0.7175575063643942, 0.021511859856701083, 0.08848574819043105, 1.6593052351894833, 0.08579326473423791, 1.5701491270167278, 0.08119952844775959, 1.4673804505995114, 0.07751062346682966, 1.3490383924899987, 0.07320276981643482, 1.219899756899908, 0.06813459484357709, 1.0818793616135423, 0.16143561647828555, 0.08567383013140582, 0.9117961919818564, 0.06178786395185505, 0.2458230107257011, 0.07889416363059952, 0.8800079497261765, 0.06016393998540656, 0.3360649644112602, 0.07295959511581207, 0.8427177081838045, 0.058715310680344794, 0.42705512716759253, 0.06635109328665457, 0.7999880355044794, 0.056074635893784856, 0.513046457224444, 0.059856406927788194, 0.7538175834199682, 0.05389109233501056, 0.5882171000483203, 0.05335819629933035, 0.7050399142710629, 0.051227179912290824, 0.15039905888991967, 1.4538696734994379, 0.1456386881655483, 1.3590040586684653, 0.13780843790044608, 1.2545581659579392, 0.13012206513196145, 1.1380329927288746, 0.12195829598663373, 1.0147467051978767, 0.23788081740930808, 0.14504553878921214, 0.8573508028046234, 0.10969157834292236, 0.32669757398652216, 0.13360908531586915, 0.8210806161507919, 0.1058028968965767, 0.415917358339372, 0.12106885658608313, 0.7801636952226589, 0.10198907719365413, 0.4989945597425267, 0.10874642197696571, 0.7357321573204914, 0.09741906005752578, 0.5731357482714089, 0.09718038786280006, 0.6887682460704626, 0.09321562329844223, 0.21966347782233261, 1.2323320834620068, 0.21095300131313974, 1.1408779027289724, 0.19802227539683415, 1.0456113562456923, 0.18567205150812433, 0.940376772269086, 0.3156518850911789, 0.21075744039962344, 0.7942346032023816, 0.16711781335645656, 0.4023887950847897, 0.18978125785156663, 0.7559911906200968, 0.1601228693739842, 0.4839935500047855, 0.17083213017155324, 0.7141655803098167, 0.15401122043980844, 0.554418236892913, 0.15103365702601498, 0.6678085015649915, 0.14633356381697735, 0.28726634563008013, 1.020502319039716, 0.2735097611315729, 0.9404356743916187, 0.25569466328109136, 0.8579276415338198, 0.3855789558505656, 0.2724153495897438, 0.7271570468461066, 0.23115505669763975, 0.46552626762068533, 0.242033325798672, 0.6867979366058545, 0.2189901003997417, 0.533977755094153, 0.21753294526491293, 0.6448652995860875, 0.208095093779033, 0.3486478355288108, 0.8346161439235915, 0.3287499758154186, 0.768537213864539, 0.442120576844606, 0.3252733123838498, 0.6580462081484206, 0.29590915600990947, 0.511175701230086, 0.2889122972260679, 0.6173649259516026, 0.2786515622577497, 0.39968330254054507, 0.6876826299569125, 0.4858094898578419, 0.3709143985293763, 0.5899847794679947, 0.35727168588519087, 0.09368028676223318, 1.8602756986995643, 0.0388248609070485, 0.17304640948087235, 1.7822144587800066, 0.03768767786120602, 0.27402426807264163, 1.6831506730647476, 0.03598638808211897, 0.3933288597489773, 1.5659578654251045, 0.034188793658879696, 0.5271059207263166, 1.4340906134319489, 0.03265023103032498, 0.6722265540722745, 1.2911586173579006, 0.030900634409058172, 0.8258485879033133, 1.1399279470839616, 0.028812835391235206, 0.17181366974493453, 1.7297806204610058, 0.09137851117772082, 0.2697330223235605, 1.6359183031641384, 0.08762557650494214, 0.38565833321450627, 1.5245849180540718, 0.08346319333359097, 0.5166400311789335, 1.3979295089899175, 0.07938110609556737, 0.6591782787839016, 1.2600976100437484, 0.074806134106358, 0.809756801206215, 1.1150996885104616, 0.06982560547661217, 0.2623796701921216, 1.5708358192368022, 0.1604967341571048, 0.3752244324962438, 1.465284609858568, 0.1530921384747177, 0.5025340962788282, 1.3463680063151549, 0.14525306423609724, 0.6409897263372523, 1.2168401459849267, 0.13668339029399665, 0.7874103723836383, 1.079974056383627, 0.12749505734135902, 0.3618796316611277, 1.3911787724825364, 0.2411785340705476, 0.4847117962752077, 1.2809469144088759, 0.22868848327242705, 0.6181458888980201, 1.1614064562423507, 0.21503478578020624, 0.7588671612322397, 1.0352835770164106, 0.2004919834845693, 0.4634132300992605, 1.2037627195884422, 0.327330312043801, 0.5911019855076631, 1.0959020099163044, 0.30776910918809053, 0.7253680827455532, 0.9824948755302055, 0.287080959555006, 0.5604116675860532, 1.0222169054773989, 0.4125289799270431, 0.6877821665725686, 0.9226225168577576, 0.38482411556957585, 0.6468377924748652, 0.8575762396042523, 0.49111524053492966, 0.16867356781231307, 1.704217993028845, 0.09052623970164278, 0.26588513881693304, 1.611693485938788, 0.08715521087335658, 0.38097421404124493, 1.503128881194068, 0.08279987109716949, 0.5104404366088778, 1.3790310807613646, 0.07873360844845037, 0.6515517144690327, 1.2441202044567983, 0.07389416777072563, 0.8012380729871468, 1.1015432360705824, 0.06924902826799906, 0.26064858519105033, 1.54673385680864, 0.15952952418261804, 0.37197333778505054, 1.4434027068849686, 0.15130048534781532, 0.49708870116664855, 1.3281421607489514, 0.14414877290462305, 0.6332021665708021, 1.2022076848431418, 0.13552326286170432, 0.7779141105392154, 1.068243132421188, 0.1265277358066581, 0.3587820470720731, 1.3717238929473212, 0.23893309125626847, 0.4796804370435458, 1.2641561443765648, 0.226475860893049, 0.6108195943132722, 1.1476608554696925, 0.21324819121912392, 0.7500795993534386, 1.0243053934567528, 0.19855347292346065, 0.4591732942539566, 1.1879586711674677, 0.32416138231857766, 0.5846586701973384, 1.083045375772759, 0.30500894483861657, 0.7170445276859939, 0.9722959511099108, 0.2843985590913445, 0.555118512045183, 1.0105552426390543, 0.4085577946615939, 0.6803802854536773, 0.913395684735532, 0.3815607566803597, 0.6406418695545615, 0.8499454939969285, 0.48658531200904215, 0.25409812824268274, 1.5080279922233413, 0.15657062083267725, 0.3636432179158387, 1.4072302855395673, 0.14897926500464856, 0.4879024594477808, 1.2956776941236796, 0.14155540402607358, 0.6219963036823458, 1.1740080906025618, 0.13303041621679565, 0.763445433053049, 1.0450013215629048, 0.12443664061299003, 0.3514938446561123, 1.3383225656315214, 0.2347762620482738, 0.47044953067885076, 1.2343370098005837, 0.2227082115437738, 0.6001275153655609, 1.121733454252774, 0.20960954157804393, 0.7370264946784432, 1.0031623361354483, 0.19543085991642617, 0.4504674718766605, 1.161630244222667, 0.31841653783741, 0.5742377069921129, 1.0597938794289368, 0.3000035813463357, 0.705140287165713, 0.9526480594913433, 0.27973732780730853, 0.5453027540508022, 0.9903259272433981, 0.40138710179816833, 0.669047066220114, 0.8959558876978255, 0.3756712724632499, 0.630240239955964, 0.8358435491943178, 0.47892919823095753, 0.3416683120744848, 1.2909626996978807, 0.22832617255562937, 0.45802821751239386, 1.1917002366596212, 0.21666798983855518, 0.5847637517858112, 1.0859328706809028, 0.20436651260154376, 0.7182263949785196, 0.9728728055198833, 0.19084875021440892, 0.4389158411814783, 1.1239335759170734, 0.3107185482866461, 0.560311781484533, 1.0269003135524442, 0.292745704445207, 0.6874417307499, 0.9260661314717822, 0.2738977299759346, 0.531842734993245, 0.9618383967990414, 0.3923367734308234, 0.6532857760920603, 0.8722609566529156, 0.3674024902676665, 0.6161137150719637, 0.8154568291652676, 0.46816037961522894, 0.42428150280721155, 1.0757625001625322, 0.30031800347318416, 0.5412444757675662, 0.9865469227724772, 0.28349814788483757, 0.6649617721551265, 0.8918597184852052, 0.2654681328276862, 0.5150608709739427, 0.9260659837582275, 0.3805194997106757, 0.6332780349663476, 0.8426130367303568, 0.35757983062987403, 0.5974888424325534, 0.7903875305702373, 0.4549088718737326, 0.4945727292839349, 0.8831507478866232, 0.36574743455401176, 0.6088682938366365, 0.8060068448740247, 0.34385257705086425, 0.5754733190607846, 0.7601458218264251, 0.4380932991482773, 0.5503715804217278, 0.72368678113467, 0.41852640827770066])
    SymCubatures.setweights!(cub, T[3.886083071053066e-6, 0.00010810278165385951, 7.515874264107571e-5, 0.00023887201563494283, 0.00026300826510193076, 0.00033138248623035786, 0.0005073045305534883, 0.0004411043246172677, 0.0007655055407634197, 0.0005013071566228895, 0.0008700404022069206, 0.0005249847285772345, 0.0008057766548891601, 0.0005093260945073036, 0.0005631745158145163, 0.0004898748636151713, 0.0005203447971631899, 0.0005395093396160784, 3.106847168110692e-5, 0.00014768012079305508, 0.0003006404792339889, 0.00045342693485733276, 0.000587102973716865, 0.0006098463456504874, 0.0005658440489076591, 0.0005431536778811469, 0.0005689448825707474, 1.3071223224596582e-5, 2.269141733786588e-5, 2.993502188840146e-5, 3.44841555845125e-5, 3.6809148187233764e-5, 3.773323861467014e-5, 3.6609006931835376e-5, 3.2776675883518294e-5, 3.485505902057857e-5, 6.321682046901766e-5, 7.866592457772319e-5, 9.038194533383355e-5, 0.00012588377181536883, 0.00010847311723000414, 0.00016068312026081572, 0.00012561866695088298, 0.00018110786380795722, 0.00013498861757758776, 0.00018007394876498728, 0.00013186262748526333, 0.00015569960371174702, 0.00011222699290108221, 0.00012380737537877237, 0.00010568253313428579, 0.00011563141472831415, 0.00014779236823064077, 0.00017030620486530382, 0.00018323057920469117, 0.00018730196355086192, 0.00017990947256951676, 0.00016503251056342895, 0.0001776005950955874, 0.00019204170832354573, 0.0002879965453281891, 0.0002484103193397145, 0.0003598244196755857, 0.00028244661041396514, 0.0004019627047425168, 0.00030345097450827836, 0.00039050064239018754, 0.00028999037697894796, 0.00034655789257200914, 0.0002673002107211964, 0.00027215933447187685, 0.0002369510817765285, 0.0003323148595642319, 0.0003816609887622506, 0.0004021387856309632, 0.00039688204678999184, 0.0003725653180533579, 0.000348331553664088, 0.0004176388138151359, 0.0003648126139162236, 0.0005304119250258779, 0.00040796085556000356, 0.0005874204126162893, 0.00042323559884629785, 0.0005732978982940994, 0.00042963220390149687, 0.0004950450966872411, 0.00037721057044542386, 0.00040161063138563155, 0.0003524092464126996, 0.0005910753095961284, 0.0006079963634256181, 0.0006251791207555191, 0.0005913490224365378, 0.0005162652981951803, 0.0006834174372976349, 0.0005108735501057162, 0.0007413691859588538, 0.0005369649240528991, 0.0007058040284125523, 0.000520795404240657, 0.0006144480434454624, 0.0004926504303770046, 0.00047708920404133844, 0.00046109302806421494, 0.0007878503706300684, 0.0007628604966232015, 0.0007540535578819224, 0.0006681676067377076, 0.0008505086827940896, 0.0006098727411433364, 0.0007983620278565114, 0.0006046247246080649, 0.000660573034224233, 0.0005590520588530501, 0.0005443295924022432, 0.0005435852780643234, 0.0008359134364381656, 0.0008163994683558118, 0.0007525019209724701, 0.0008620843389625517, 0.0006375581486982812, 0.0006899286629969413, 0.0006180858329571642, 0.0005802790333780978, 0.00059143195334495, 0.0007785832596742558, 0.0007368537864456463, 0.0007086195510571507, 0.0006311976099042361, 0.0005851674783392712, 0.0006015398337296288, 0.0006591448730891366, 0.0005747556632990825, 0.0006700096091495865, 5.4687032988916415e-5, 6.933761997300757e-5, 7.859164691888854e-5, 8.183649703591253e-5, 8.080714977173538e-5, 7.583002445975926e-5, 6.984755815629731e-5, 0.00010256292736009522, 0.00011520650493577006, 0.00011891467298094935, 0.00011967230747557665, 0.00011650073532805943, 0.00010226199857312642, 0.00013890077733510805, 0.00015489517847744422, 0.00014768513959001727, 0.00013939045526926503, 0.00012598610123434745, 0.00016397590091615704, 0.00016683358039683058, 0.0001588741491725702, 0.0001497583692956833, 0.0001768822158818094, 0.00016630355996354786, 0.00015422335164840288, 0.00016066165810473454, 0.00014754887261606294, 0.00013376847054960838, 0.0002305705777582581, 0.0002617635319276331, 0.00026846225410250067, 0.00027125670792412405, 0.0002565217069306666, 0.00023295908812253035, 0.0003171327326526032, 0.00034128853043440905, 0.00033018877510582217, 0.0003175079517579136, 0.00029110420075579495, 0.0003783440437551694, 0.0003728799271985556, 0.00035368330583645884, 0.0003181305835570283, 0.00039495347717054, 0.00036816274367155387, 0.00034025405740088126, 0.0003647953359659455, 0.00032484103585851925, 0.0002895885584196321, 0.0004746729363069271, 0.0005012773215169735, 0.000497658847454874, 0.00046547068831814115, 0.00042736325299394024, 0.0005642803851249634, 0.0005597180338693875, 0.0005231029638839835, 0.000467445430779183, 0.0005739622587053986, 0.0005395297934040871, 0.0004835522920044635, 0.000532886198969197, 0.0004731250029255267, 0.00042733132868492984, 0.0007065366590712515, 0.0007088834937801795, 0.0006505076751502686, 0.0005946600835450333, 0.0007250910161847849, 0.0006743966711751422, 0.0006063286351736413, 0.0006452539653172381, 0.0005705693675470717, 0.0005323467402837545, 0.0008199847024479732, 0.0007733000896749448, 0.0007054677126018175, 0.0007505434398535203, 0.000660753207552903, 0.0005976136680977513, 0.0007970345543031147, 0.0007518248028137795, 0.0006370683162375709, 0.0006752362436943432, 0.00036759723438345573])
    cub_degree = 36
    tol = 5e-14
  elseif q <= 38
    SymCubatures.setparams!(cub, T[0.00862737903799628, 0.05133986510315037, 0.12218723640382904, 0.21383631635518835, 0.31845813100843756, 0.4224362945738618, 0.5179025101668452, 0.602002944251124, 0.677661440330468, 0.7274297644003183, 0.006468387694878754, 1.956081898705295, 0.00670121885379102, 1.9092704474722233, 0.006624113834473502, 1.8417182365254188, 0.00635193975756349, 1.7535624952148752, 0.006010767088271724, 1.6454978131364963, 0.005743925385723461, 1.5194657649587375, 0.005533590246721263, 1.3793467839058635, 0.005155522836670435, 1.2300227757118036, 0.004745368269407758, 1.074504207717183, 0.03438524992921262, 0.0065782506754793494, 0.08547582754032569, 0.00632934911019405, 0.15255378458247998, 0.006049996843794977, 0.23221419421563888, 0.005855875239609631, 0.3193031787164848, 0.005393629080591188, 0.40866401301675037, 0.004920063807187139, 0.4944954348527162, 0.004572687474281646, 0.5722460539624292, 0.004059834047343143, 0.6360449567222778, 0.0037491315587739014, 0.03431192128505527, 1.8489039320977052, 0.033837662528190274, 1.7825954652075893, 0.03284419248020252, 1.6988484558125447, 0.03152047457578716, 1.5973083041780765, 0.0301133720368696, 1.4779425254625829, 0.028378124538781437, 1.3455025422378537, 0.026253270704827584, 1.2031653227194006, 0.024715533830355612, 1.0529013827750133, 0.08366976410440627, 0.033285480724169994, 0.15141013898035277, 0.0316317281415504, 0.23041662356785683, 0.030087245420087443, 0.31673717119829387, 0.027893415013004536, 0.40510876927697037, 0.02584651436677703, 0.49062606777585244, 0.023652594096032774, 0.5677712868199858, 0.021077792149042308, 0.6319808616754125, 0.019471484992507684, 0.0809461971156054, 1.6884944627388951, 0.07864653308957759, 1.6069195494171966, 0.07489900916590854, 1.5122451731111333, 0.07200946171618461, 1.4023871011366214, 0.06827243240940514, 1.2817568808718545, 0.06373173961106525, 1.1520283204440764, 0.059508207420990895, 1.013587825596617, 0.14783462765073094, 0.07727540370222646, 0.2263652288511168, 0.0726285653715558, 0.3114556859422303, 0.0678526847104048, 0.39830231886928646, 0.06260923902168083, 0.48256692257060313, 0.05658174832932556, 0.558904688307751, 0.05124357638938666, 0.6226593971870155, 0.04622028279057291, 0.13874873962118572, 1.4974858890122842, 0.13453503374838416, 1.4083625074245436, 0.12840248206939245, 1.308856918569283, 0.12169105129246612, 1.1989070210004311, 0.11454268754730021, 1.0819722962908789, 0.10662681020386883, 0.9600500894208003, 0.21935839972061397, 0.1340324883658162, 0.30324571648314735, 0.12498783776655269, 0.3886291484989987, 0.11378705889608255, 0.4718709160838844, 0.10342379081599856, 0.5471632764807275, 0.09324297216816868, 0.6103098905063252, 0.08406779821694103, 0.20386239442757756, 1.2876845633412797, 0.19619545633521618, 1.1992188149009089, 0.1864988569186089, 1.1046724117450633, 0.17526316200471342, 1.0050120017203579, 0.1637155777105711, 0.897788049290156, 0.2941689914772817, 0.19624385352125262, 0.37770623894613875, 0.17894409725577742, 0.45789458267484084, 0.16304635512109855, 0.5323050153029614, 0.14643173637808204, 0.5947399642818744, 0.1324130522333495, 0.26942445780320234, 1.0813485263901723, 0.25733537548880897, 1.0008067642246519, 0.24260972077104578, 0.9187466244218042, 0.22781411434270477, 0.8281612987506057, 0.36263308532789856, 0.25706971111558835, 0.44247551735040586, 0.23179740171666963, 0.5153394974568388, 0.21080757921478846, 0.5771476277439198, 0.18946071445299023, 0.3300393423335487, 0.8964865070734083, 0.31256766022362437, 0.8253465684223706, 0.29490949953143636, 0.7532249124458454, 0.42267872003748436, 0.3123678882569139, 0.49547525915041507, 0.2818516243341372, 0.5562245128364717, 0.25559700713733946, 0.3835310429851673, 0.7348308641802014, 0.3621697564117793, 0.6758834795473111, 0.47339106095502703, 0.36378247490422555, 0.535184470169523, 0.32596551549904085, 0.42671702424312064, 0.6057702942320087, 0.5107239251918496, 0.4056064495035273, 0.08471579992449887, 1.873606846619408, 0.035141123011881326, 0.15679561441039208, 1.8023058563520344, 0.0344657757980861, 0.2490258041545216, 1.7117752364292331, 0.03298281966148709, 0.35865108857088945, 1.6039623543558834, 0.031366328036317954, 0.4822591866744156, 1.4820898004151883, 0.02996697591660192, 0.6170563186624951, 1.348823671593745, 0.028743901713825754, 0.7609284259266794, 1.2073413098911105, 0.026815952391730267, 0.9102124356792451, 1.060259252520273, 0.024830813680327603, 0.1563154526885006, 1.7537720208094394, 0.08354226486037414, 0.24594009245203938, 1.6674870583574044, 0.08034245986399799, 0.35252008524807876, 1.5651630621549981, 0.07663643741450044, 0.47363735137203217, 1.4474607684359782, 0.07314486750552672, 0.6063666678716375, 1.3188154529202252, 0.06953420588971079, 0.7473111444795689, 1.1824677534240648, 0.0651481285019266, 0.8940000514273814, 1.040710646401072, 0.06056978120649507, 0.24003743917414422, 1.6068671192495922, 0.14738206335256004, 0.3440898425186433, 1.509268260341783, 0.14065081057455953, 0.461988742237371, 1.3982705291789272, 0.13435181279554545, 0.5910539176895068, 1.2766389218611676, 0.12698043345931542, 0.7285982558558151, 1.1472448059939047, 0.11931380685805983, 0.8714633284148069, 1.012740398663051, 0.11123054573979753, 0.3328841411562329, 1.4395628458573244, 0.22221815731795527, 0.44702222613219916, 1.3358844172781816, 0.2117987255798131, 0.572036999720454, 1.2226825225647184, 0.2002904018703028, 0.7047894850185961, 1.1023574557750324, 0.18813198014663254, 0.842938525184701, 0.9771524960627875, 0.17535233213450868, 0.4292953095158836, 1.2618643484405303, 0.3038096438526283, 0.5493213797456709, 1.1583425775739944, 0.287392032363698, 0.6766533623962196, 1.0488132871199756, 0.26994538567388454, 0.8089505546779432, 0.9352439272448607, 0.25145658677156596, 0.5233919039306754, 1.0858116094849877, 0.3860572216027223, 0.6448366137793823, 0.9882908139082008, 0.36243614515786493, 0.7706640784673834, 0.8870889106071793, 0.3379241105010368, 0.6100037082256162, 0.9217695212446964, 0.4642114882701555, 0.7288065996169639, 0.8345079381499741, 0.4325490086045207, 0.6841889398085256, 0.7781133812763548, 0.5337228011149512, 0.15354411051019945, 1.730634667717488, 0.08246873362088575, 0.24260714182084891, 1.6451678458743424, 0.07979549696471874, 0.34870769256291245, 1.5449515561220424, 0.07603316773951914, 0.46820338363637415, 1.429086935753297, 0.0729354209035634, 0.5994627378737557, 1.303744395976856, 0.06887404092210758, 0.7398114103282439, 1.1695820597736428, 0.06435394652731521, 0.8853686382870901, 1.029892044239677, 0.06019071506619334, 0.23851263125812533, 1.5848139619824104, 0.14617988378862976, 0.34115138302227155, 1.488698367678292, 0.1394118511034496, 0.4573887106843634, 1.3806916842731973, 0.13323752254212642, 0.5844017653924225, 1.2616536677003238, 0.1264343751048263, 0.720292469812986, 1.1354627881539348, 0.1185514718228809, 0.8626719966449489, 1.002947488849055, 0.11031165497098509, 0.3304849773860963, 1.4205200500860773, 0.22040157164356794, 0.4430746166071567, 1.3195188045929631, 0.20983885371347294, 0.565805644114653, 1.2089307605285224, 0.1988405560935943, 0.6970305874651623, 1.0912852411828529, 0.1869539297363601, 0.8343610494896703, 0.9681670973933999, 0.174019480590529, 0.4259687254061599, 1.24631721304131, 0.3011965331788109, 0.5441965871028467, 1.1455501743138379, 0.2846010490856416, 0.6693575220748199, 1.038799940715609, 0.26775788150616264, 0.801041806028025, 0.9265455620601643, 0.24951022845459642, 0.518945482992081, 1.0738305035702764, 0.38271753834557404, 0.6385898210073568, 0.9787531486590139, 0.3594131454779532, 0.7632525531840106, 0.8797710178873233, 0.3348302851063268, 0.604985846685228, 0.9132060290539655, 0.46021073431554693, 0.7219935468746355, 0.8275470428955126, 0.4290090060804406, 0.6781069310325007, 0.7716049678898887, 0.5299404276727466, 0.23276172782259796, 1.548145842852985, 0.1438239030559104, 0.33421720945088207, 1.4546995439241515, 0.13747658317267858, 0.44949322540108116, 1.3488871602880608, 0.13137427057756884, 0.575236096154936, 1.2342006815793378, 0.12416372070847666, 0.7085914131978512, 1.1121326797614692, 0.11653407615469184, 0.8478425919008141, 0.9848545448408523, 0.10878970412708605, 0.3240992929466443, 1.3885939005380326, 0.21698326174170499, 0.4351475874282241, 1.2906205987073627, 0.20675657989117283, 0.5567512189722339, 1.1831974814772492, 0.19552833126581054, 0.6859338995789646, 1.0696903056061284, 0.18393078270723212, 0.8209928046652968, 0.9509669369512241, 0.17139595486540563, 0.4184106361572162, 1.2206892504715503, 0.2962045778023287, 0.5350767079106563, 1.1226814860538359, 0.2804305165025022, 0.6590611110581294, 1.0187573732005826, 0.2635471312483895, 0.7882096681048938, 0.9108143766947513, 0.24583775778127173, 0.5101653359461459, 1.0537824339645108, 0.3766924800840004, 0.6289198554505382, 0.9611191530385711, 0.3539028213180966, 0.7518655005534012, 0.8644372744115842, 0.33049931903015844, 0.5954617843498471, 0.8977557769495408, 0.4537094339577678, 0.7116280738655714, 0.8140495890240758, 0.4231937366917767, 0.6682896083125343, 0.760387262756227, 0.5225142665486524, 0.3156875717208477, 1.3439350474162264, 0.21171749075661084, 0.424271802154873, 1.2493913326186934, 0.20168093017341698, 0.5433397198677125, 1.147495870317398, 0.19108789168248774, 0.6705380357533643, 1.0386873974110205, 0.17969838805282312, 0.8019956051550025, 0.9263966449964013, 0.16773729134780307, 0.408439716066018, 1.1834306243289032, 0.289308951974559, 0.5231404100542159, 1.0899454518583784, 0.2739714831554656, 0.6438502709624608, 0.9913277334664365, 0.2578744215890205, 0.7703655132951613, 0.8883996628148367, 0.2411377713487128, 0.49855162708286155, 1.0246091470710235, 0.36862867669727256, 0.6152562626804997, 0.9361493540641274, 0.3466004184007058, 0.7349960405613973, 0.8444699216915751, 0.32418951884115116, 0.5827808213640172, 0.8765530818922169, 0.4445204024547904, 0.6961056292489566, 0.7959407534751203, 0.4153047259156953, 0.6547846466540085, 0.7447644014447942, 0.5122560892974708, 0.39584958823478855, 1.1357643351743387, 0.2803830313356274, 0.5063753542970286, 1.0483064697718116, 0.2665655121931819, 0.624877254779684, 0.9566038566418138, 0.2504626828870918, 0.7470978097024676, 0.8598626330124445, 0.23470664544079248, 0.4842686852646821, 0.9883360326513716, 0.3575351868914664, 0.5971041778948444, 0.9052044728232842, 0.33762529510950234, 0.7140115108509981, 0.8186738377340885, 0.31651555529711284, 0.5665446842354587, 0.8491892482691605, 0.4330357688512156, 0.6769642575772176, 0.772481712365671, 0.4052282757154331, 0.6382377045684434, 0.7246652194179285, 0.49946836593752403, 0.46660154336200965, 0.9436225431040435, 0.3447485759914316, 0.5764345130041014, 0.868061803113585, 0.3257846486532807, 0.6883537933632832, 0.7876198987050782, 0.306072490468803, 0.5470228076188501, 0.8159053365269854, 0.41864015499424484, 0.6538770692507893, 0.7455618541093683, 0.39308111892227104, 0.6172005301015533, 0.7002367842757887, 0.4844425763376754, 0.5245894587011785, 0.7776548990224793, 0.40098517999860067, 0.6283849006878046, 0.7131320233062907, 0.3780314778879857, 0.5943460090820747, 0.6715380416384631, 0.46765608348855675, 0.5690585415361065, 0.6411351166602715, 0.44694492071875236])
    SymCubatures.setweights!(cub, T[3.2719271313082825e-6, 5.586675355325429e-5, 0.00019977872260599483, 0.0004066441185812722, 0.0006113110159954992, 0.0007442309830536526, 0.0007237250395841257, 0.0006737319696482052, 0.0005078161737416776, 0.0002901500923638866, 9.858900681482128e-6, 1.683170921575801e-5, 2.2546564290396962e-5, 2.619982634847146e-5, 2.8097888319417403e-5, 2.9102695654966605e-5, 2.9272974797727198e-5, 2.6743271583184833e-5, 2.3362581732885822e-5, 2.5766794959187633e-5, 6.0115093242109214e-5, 9.508780551004748e-5, 0.0001279071488274041, 0.00014575173400921968, 0.00014770282821944581, 0.00013735126169145745, 0.00011011400718164452, 8.26297542462283e-5, 8.675930735719455e-5, 0.00011146131288938815, 0.00013063372462057256, 0.00014114695932409765, 0.00014929490708509108, 0.00014419268473661403, 0.00013370784086220744, 0.00012045509646288228, 0.00013464086714866917, 0.0002182083094211129, 0.0002826859277684064, 0.00032488537663238177, 0.0003332739811116, 0.0003024186212419062, 0.0002482225254693872, 0.00017698533422394146, 0.00025494171817024296, 0.0002929409539279364, 0.0003159963629422833, 0.0003155590157666518, 0.0003091456157017859, 0.0002891129795843755, 0.0002607759715618227, 0.0003195751421259658, 0.00041881076486075636, 0.00048532528087595183, 0.00048503586495467606, 0.0004388102829633639, 0.0003617778873019976, 0.00026058964373557047, 0.0004648546699475561, 0.0004979904352913158, 0.0005003962983966075, 0.0004876773659527273, 0.0004446212206681516, 0.0003973946364384667, 0.0005465101591946391, 0.0006087848689074305, 0.0006030831526383381, 0.000560957698977965, 0.00044513524434202094, 0.00034462238785831734, 0.0006633109488916284, 0.0006356007568712715, 0.0006248837652009318, 0.0005743817623039469, 0.0005247183477315806, 0.0006975432832036552, 0.0007038875075657503, 0.0006285172994280415, 0.0005225141776840562, 0.000393542440845961, 0.0007221669318856896, 0.0006815385715689482, 0.0006566603712534479, 0.0006022749542189338, 0.0007484911532720743, 0.0006708752981851903, 0.0005466122448310578, 0.00042737674558552445, 0.0007002640656650868, 0.0006831268291791804, 0.0005997764630924765, 0.0006947164528182116, 0.0005675794733554913, 0.0004205735202870612, 0.0006415038015816846, 0.0005241768513160336, 0.0006011405432713408, 0.00040829090013110963, 0.0004471694417343044, 0.000389030488333512, 4.0703868310632936e-5, 5.265133535947842e-5, 6.007106759007905e-5, 6.38505676221702e-5, 6.375650126300805e-5, 6.177326800050951e-5, 5.5813294699842595e-5, 5.053001660071662e-5, 7.781046960340765e-5, 8.974574673953857e-5, 9.140714303088638e-5, 9.76531555683959e-5, 9.102701893948147e-5, 8.563689045826225e-5, 7.665036671326828e-5, 0.00010746288033130129, 0.00012313731201646014, 0.00011852743909622753, 0.00011738239451381644, 0.00010598825874613238, 9.541246521421958e-5, 0.00013153677928767087, 0.00013552045349643694, 0.0001306170234437544, 0.00012022273152343918, 0.00010950020661033098, 0.00014339997732133104, 0.00014087538510713116, 0.0001272831058391115, 0.00011520730960082064, 0.00014105341954898456, 0.00012791661065458972, 0.00011651261043115836, 0.00011731969631692553, 0.00011038811428388161, 9.925545816279351e-5, 0.00017634418177857385, 0.0002018853642873283, 0.00021179624683215028, 0.0002143745289124401, 0.00021020363184045532, 0.00019240119986681218, 0.0001703280676550032, 0.0002496108367687665, 0.0002673679706889319, 0.00027014307880151927, 0.00025860361161906186, 0.0002444788965775184, 0.00021744289699096022, 0.0003065127557220846, 0.0003022679740443356, 0.0002941750958075251, 0.00027080830473067314, 0.0002433185313054651, 0.0003256734913883249, 0.00030841202140863446, 0.0002866221693121094, 0.0002600197691947609, 0.00031377200943543684, 0.0002841754820915072, 0.0002562253784921966, 0.00027130785793565355, 0.00024575178220833966, 0.00021804010954919086, 0.0003743756182011479, 0.00039757848668563834, 0.00040308943474452996, 0.000382138840980024, 0.0003572045770843399, 0.00031773837594862854, 0.00045110888058021614, 0.0004593562261618135, 0.00043871433065972105, 0.00040580590909382856, 0.0003598583210042021, 0.0004803091359516338, 0.0004574183072337277, 0.00042303629134197816, 0.000375866697811508, 0.0004636721659871044, 0.00042026780803838104, 0.0003748772008630468, 0.000400742690261997, 0.00035088671190139223, 0.0003206879462142952, 0.0005663908544559531, 0.0005797802659901327, 0.0005517989192362435, 0.0005079956575450666, 0.00045097868097547435, 0.0006143883383077764, 0.0005821277869595179, 0.0005270529037913316, 0.00047737597800814896, 0.0005818066773588316, 0.0005249612659593677, 0.00046191933527146384, 0.0004917096880624012, 0.00044619973577431446, 0.0003911928454379417, 0.0006931340710903759, 0.0006662208032631113, 0.0006080881158863594, 0.000552530035463541, 0.0006587857832132346, 0.0006029606281887187, 0.000542152657980133, 0.000578781064248337, 0.0004987077411406144, 0.0004545313508682264, 0.00072433802429692, 0.0006529039878179267, 0.0005919505270330603, 0.0006322655749792029, 0.0005448141586666867, 0.0004934623160494074, 0.0006695082384519306, 0.000563853536167134, 0.0004926886027184743, 0.0004970040791906683])
    cub_degree = 38
    tol = 5e-14
  elseif q <= 40
    SymCubatures.setparams!(cub, T[0.00780330525525622, 0.998493752865077, 0.04645725589611356, 0.9917644702538027, 0.11174125604356362, 0.980224191438278, 0.19711768630471135, 0.9636440633642875, 0.2953832122443872, 0.9432574598344768, 0.3962210336631515, 0.9177684131869567, 0.49152901827090406, 0.8903708290323409, 0.5742677642603903, 0.8593034295516675, 0.6443184359345666, 0.8268957794482568, 0.701852790353926, 0.7888554625061005, 0.9958440067364807, 0.9784124162903821, 0.9482482433707597, 0.9072874428195205, 0.8563095234372619, 0.7981838828618165, 0.7369218364173937, 0.6776892051204217, 0.6190459664782587, 0.5574412454896516, 0.0058828858132688455, 1.9601728001741805, 0.006109485459602238, 1.9176775033568327, 0.006036395021978035, 1.8563064794379012, 0.005801038863192481, 1.7759005238841472, 0.005537941187594949, 1.6768088325859962, 0.005319640090786907, 1.5607999731790514, 0.005091482086901414, 1.4310646471562276, 0.004797544221997284, 1.291529754406188, 0.004489272746612871, 1.145435323859664, 0.031297717944662184, 0.005944642503171424, 0.9872664129398502, 0.003976139640472778, 0.0780735583462331, 0.005765505479523608, 0.9715911489197082, 0.0042694502091038585, 0.13969820302080832, 0.0055962357945114225, 0.9493612293971047, 0.004026567493826906, 0.2136654662015896, 0.005346987697957068, 0.9214490571531605, 0.003932906682549396, 0.29532516517481355, 0.005058675611490092, 0.8890464655990614, 0.0037249938090554305, 0.3805813130074721, 0.004684598763412118, 0.8507682326892949, 0.0036826378721035143, 0.46431604379352087, 0.004156079153129802, 0.807587934592507, 0.003510554194803897, 0.5422440547398109, 0.003852431509136228, 0.7628426535224752, 0.00345213371319157, 0.6121591818954955, 0.0030953684521964823, 0.7204216109896704, 0.003704234898837557, 0.031379299626650765, 1.8624187851139853, 0.031061462059898473, 1.8017227237023616, 0.030185382412241932, 1.7251252775195411, 0.029066537930402297, 1.631903553838148, 0.027875493668759216, 1.5219341374313855, 0.026421222548616917, 1.3990748298482503, 0.024731152756800965, 1.2657182198975099, 0.023254737806550917, 1.124519873402946, 0.07653300084047618, 0.030375176086618508, 0.9627342385714175, 0.021587622700813328, 0.13891465107510692, 0.0292568681818942, 0.9412548465748816, 0.0213810989423556, 0.21221891867943393, 0.027751842881520713, 0.9142988336238946, 0.02054528460933853, 0.2932735385616038, 0.026099416187304384, 0.8813476818825909, 0.019973557194804956, 0.37765412541437426, 0.024163722501851965, 0.8441207397661412, 0.019194450921857457, 0.4610069602614224, 0.021847100755721313, 0.8024595755991555, 0.01876181172143078, 0.5385574816241506, 0.019829740179515054, 0.7592279249353137, 0.01817794485654486, 0.6058289998480808, 0.017684643266719622, 0.7125967583581571, 0.01684757644680894, 0.0742865329746897, 1.71433225005506, 0.0723183968466602, 1.6387876607166394, 0.0690980960191792, 1.550665554683651, 0.06641526344689719, 1.448369029728259, 0.06350041068336028, 1.3350920780853586, 0.05993354672228551, 1.2128268341663957, 0.055989022313176394, 1.083188399842525, 0.13597408376823608, 0.07137647810502702, 0.9269516285016883, 0.05122784981979309, 0.20900331882079984, 0.06747663764896816, 0.8999081929941254, 0.050235347489872614, 0.2890037948976987, 0.06337012983673694, 0.8679316891223846, 0.04883184678534311, 0.37192688878095226, 0.05846960205142487, 0.8315057473107447, 0.047186468900143046, 0.4540137345442682, 0.05315489482271383, 0.7918281922415308, 0.045750833384987644, 0.5299503168242145, 0.04833865446527458, 0.7480035448526876, 0.0436297024650081, 0.5967788559286459, 0.04294298032449289, 0.7023162717220909, 0.04201661210619631, 0.12805907349996906, 1.5363106289191182, 0.12462455674367147, 1.4524494189651618, 0.11907435858868169, 1.3589883363373856, 0.11368869755571498, 1.2545999482603964, 0.107370875511136, 1.1439028221718874, 0.1005594119586508, 1.0277628906966114, 0.20275194280242587, 0.12451739960191574, 0.8805506775280213, 0.09154588224049484, 0.2813952432215588, 0.11674601237323147, 0.8487824846520727, 0.08904064186090631, 0.36407396064559544, 0.10667607086958295, 0.8141753537851664, 0.0870650937077143, 0.4446771382228984, 0.09707943691358042, 0.7762222087349298, 0.08354888587983, 0.5194998134150152, 0.08810216065430934, 0.7345648080065215, 0.08006488445084142, 0.583157271043959, 0.0791427397764501, 0.6888076585653634, 0.07649730315713078, 0.1897919421901424, 1.3378326119600155, 0.1833345347706287, 1.253638617966338, 0.17439370137012333, 1.1630885553198387, 0.1648741942717928, 1.0654351677583465, 0.15480667397350278, 0.9623140087174925, 0.27412543578858084, 0.18379612396740838, 0.8260164402274832, 0.14044793596510366, 0.35404448083469986, 0.16798954040106148, 0.7918494173419627, 0.13621034973413215, 0.4322676502294918, 0.15361221468670047, 0.7559749588497939, 0.1320355383391761, 0.5054950225556303, 0.13856578215741094, 0.7157278179303439, 0.12646103973318007, 0.5689541731043642, 0.1252639251455935, 0.6726159491725323, 0.12004871850279192, 0.25321867883481763, 1.1373738956201747, 0.24263615994312812, 1.0601746839286135, 0.22934415535112604, 0.980417802019005, 0.21641185489472292, 0.8919613044938237, 0.34183421851693996, 0.2427829312580798, 0.7664331738039935, 0.19651133604548954, 0.4186904751163051, 0.2197127334376235, 0.7321279583811412, 0.18871833577150646, 0.4892915771060803, 0.1996533107528508, 0.6949697459809466, 0.1818096161756553, 0.5508237711476772, 0.1773574752647182, 0.6527465129833759, 0.17278695397976734, 0.31256690205842413, 0.9532495295585618, 0.29851914997744317, 0.8859592314571256, 0.2808368161536759, 0.8177156478837411, 0.4011256925246551, 0.29760967675221534, 0.7048098421275736, 0.2560150841810342, 0.4707999251333685, 0.2657189074041359, 0.6696855817467668, 0.24379066078430178, 0.5313512910605486, 0.24066340214691576, 0.6314268775704904, 0.23199294921356722, 0.3656745973420324, 0.7928281498445406, 0.34691924743236274, 0.7385772313884305, 0.44664981499916395, 0.34330128098238544, 0.6433281817608787, 0.3154366121268928, 0.5093517235486231, 0.30527978571932957, 0.6061268022119838, 0.29710694392160486, 0.41129839743695984, 0.6618898172184555, 0.4845531919245438, 0.37873563124070153, 0.5812587303749379, 0.3685703012808796, 0.07715568682702717, 1.8848007938758473, 0.03206143216772718, 0.14294945115745722, 1.8197124330629229, 0.03143226632855319, 0.22744222659264407, 1.7366684496698153, 0.030171864815979565, 0.32839674735712443, 1.6371360108470996, 0.028924595510418344, 0.44286407764286184, 1.5241718692888868, 0.027707961276903385, 0.568121652976517, 1.4003052663228184, 0.02655175326114717, 0.7023342543371837, 1.2679946774331972, 0.025016033953573116, 0.843041734477596, 1.1291280930540017, 0.023424282767248856, 0.14298044187557063, 1.7748529462429838, 0.07632910415927463, 0.22531985041661529, 1.695409992774357, 0.07359101779513941, 0.32344288801690096, 1.6005905229930208, 0.07067325936015707, 0.4354664218082627, 1.491681134909094, 0.06763067586646822, 0.5590164707768701, 1.3714502100557413, 0.06456224807686642, 0.691207464382688, 1.2431563159092196, 0.060941915182895004, 0.829396503935381, 1.1091785780675236, 0.057007628404800746, 0.2202530684972627, 1.6391631999155778, 0.13525412051351177, 0.31634061509980466, 1.5483589294489595, 0.1297829369259209, 0.425795260190744, 1.4448490134425356, 0.12438828882298822, 0.5462185696520654, 1.3306018763502807, 0.1182317435480304, 0.6753541390265299, 1.2084240853516597, 0.11174426172102411, 0.8103888296891036, 1.0804054736355493, 0.10473842020085145, 0.3070655126719125, 1.4828321811087262, 0.20517580034222194, 0.41322601335352943, 1.385368115692493, 0.19657081930142126, 0.5301847933774515, 1.278453586698146, 0.18671735617107552, 0.6553877797627854, 1.1638280446956877, 0.17631414353305694, 0.7863678324659987, 1.043990963310704, 0.1653934228151322, 0.3981357119745141, 1.3147131640533323, 0.28249427784999886, 0.5110610572091873, 1.215859592351918, 0.2684442110604749, 0.63162915203647, 1.1106155509546374, 0.2534571956863742, 0.7574713689189585, 1.0007645105681242, 0.23769244877318585, 0.4888477318986891, 1.1452523831189336, 0.3616144116509917, 0.6045368199039223, 1.0499853546552187, 0.3413658820306683, 0.7248994296487496, 0.9509581547607739, 0.32018968177005896, 0.574456703822437, 0.9834230940414245, 0.43821772698297573, 0.6890371579727865, 0.8962415517196841, 0.4110956615488823, 0.6505750756450451, 0.8376901532021493, 0.5084120580902062, 0.14052125148718372, 1.7532365425054162, 0.0756748445779972, 0.2223049318252422, 1.6746762155016006, 0.07332668101535529, 0.3201008661553161, 1.5816408129240458, 0.07018600741984427, 0.4310835381919044, 1.4741128632976308, 0.06759284665873853, 0.5532599702487863, 1.3566629948003184, 0.06395849913849418, 0.684634656394984, 1.2305989553757706, 0.06024819975985494, 0.8222232454061247, 1.0980057945889996, 0.056736613310884414, 0.21909795684288358, 1.6182074109710212, 0.13443631899721356, 0.3141248119152104, 1.5287273739511884, 0.12870118133335276, 0.42220686576134636, 1.427710109650065, 0.12365897118873875, 0.5409563380736065, 1.316002820774294, 0.11745523255851714, 0.6684096667395474, 1.1968148331169435, 0.11089432220000923, 0.8021475117221667, 1.0707594216164762, 0.10402276261965863, 0.3049973937900846, 1.465002527609978, 0.2036690293118738, 0.4101674435100587, 1.369637752771689, 0.19491765973281122, 0.5253691468567119, 1.2646932221422307, 0.18532909415785648, 0.6487696204454886, 1.1527289644987584, 0.17521284076402396, 0.7786074593699304, 1.035191060265319, 0.16420069732333886, 0.3954945932852071, 1.299705921229436, 0.28008887045611247, 0.5070630596314076, 1.2029040939661089, 0.26616265307450415, 0.6254990875433489, 1.1001604435884378, 0.2516487673942523, 0.750097522315124, 0.9925338782146623, 0.23591094945938781, 0.48523526593714233, 1.13340256951507, 0.3587041141092554, 0.5995441503485262, 1.0401689333365836, 0.3387929483159857, 0.7183778088546487, 0.9432759884733457, 0.3178692728588705, 0.5702092137375855, 0.9745153805435947, 0.43487211786629315, 0.6831610560010081, 0.8890479174929353, 0.40854379052142775, 0.6457078780223664, 0.8316960643341577, 0.5046516292609412, 0.21450821859228816, 1.5834646890407893, 0.13257601119307955, 0.30832597960038616, 1.496433210901876, 0.1269256345529843, 0.4153951837799907, 1.3978950980802627, 0.1218823389731849, 0.5328529335183829, 1.2895796921208182, 0.11560497776943837, 0.6584564580216028, 1.173554908071805, 0.10947879544441946, 0.7901346802490147, 1.0518059540267792, 0.10268954943206571, 0.2999086441593928, 1.4343520893363197, 0.20092575104519245, 0.40334969607694937, 1.3419353205035542, 0.1923563084154565, 0.5173383655345657, 1.2401185613371073, 0.1826481677647537, 0.6393425314087763, 1.1312457667784064, 0.1724528075689049, 0.7674034943520455, 1.0174354730300665, 0.1620357280215891, 0.3892041481509033, 1.274558561522094, 0.276067765459665, 0.49888992480668354, 1.1806599340725201, 0.2626936379462337, 0.6165661528746503, 1.080423873871388, 0.24793475419010758, 0.739689517297324, 0.9757666024597649, 0.23254300326288432, 0.47786996667402537, 1.113141761390223, 0.35371476424842396, 0.5908776055701332, 1.0225503014091128, 0.3342800296835896, 0.7085684955051411, 0.9278142134768731, 0.313921035691027, 0.5620420874446619, 0.9589295643376913, 0.4292173755071477, 0.6739912788773836, 0.8754078109234771, 0.403374982498494, 0.6373868078570746, 0.8204941618025243, 0.49797933161687363, 0.29240022671251376, 1.3922966905541316, 0.19596775129885013, 0.3939341519551235, 1.3021677058823573, 0.18799217893748746, 0.5061163647926148, 1.2051692105850984, 0.1789517450093806, 0.6260093807520594, 1.1005273773784479, 0.16891903608459768, 0.7512198501158079, 0.9921195348707843, 0.15872122458295246, 0.3802927602139976, 1.2387987868895565, 0.27053980426675917, 0.4888314952296137, 1.148238773164083, 0.2571336734694442, 0.603851493355202, 1.0526853029343575, 0.24288421496993332, 0.7241222028126749, 0.9527038939610596, 0.22873200611520675, 0.4680073338902228, 1.0842944854139172, 0.34644883210380045, 0.5789734135009265, 0.997344903237675, 0.3279612131895608, 0.6938488128573709, 0.9070261009970451, 0.3081091182991775, 0.5508379493779388, 0.9366800138170718, 0.4207484349929641, 0.6607847123810955, 0.8565315737559065, 0.39640804163669063, 0.6251143634626606, 0.8048932323864286, 0.4885239427326082, 0.3698588006434931, 1.191926412156377, 0.26272319786218734, 0.4744502702437166, 1.1073985756996452, 0.25058033157171067, 0.5871044073478657, 1.0178790378637845, 0.2366702957781978, 0.7042185403915016, 0.9234015875024685, 0.2224420842955376, 0.455350648848083, 1.0479892225305172, 0.3373802947901438, 0.5631884774564808, 0.9657043267046809, 0.31945190220049063, 0.6754698377522926, 0.8804989952614686, 0.3013829848784477, 0.5359950534480435, 0.9089742152974956, 0.4108143262788168, 0.6446989369543636, 0.8328761012481055, 0.38720968535028, 0.61017058322892, 0.7844153891406758, 0.4770994013289756, 0.44041229029225437, 1.0034367032395912, 0.3262427095471687, 0.5442223803900844, 0.9274326494269797, 0.30958268832135005, 0.653321726390288, 0.8488115553225367, 0.2916290209696817, 0.519767669496071, 0.8760075939067802, 0.3976760165113464, 0.6238417499220121, 0.8057144499785194, 0.3765915284647388, 0.5914469906332879, 0.7613454762163736, 0.46281848845344775, 0.4995895236022162, 0.8382558966255147, 0.38277808943628616, 0.6011921634381396, 0.7721065612620002, 0.36166134882565476, 0.5712286846135588, 0.7328432591462286, 0.447655298008404, 0.5476629836999174, 0.7009877484704158, 0.4282591354605115])
    SymCubatures.setweights!(cub, T[2.423642464851422e-6, 7.195986953917698e-5, 4.164606084048165e-5, 0.0001535242371805211, 0.00015410537948731454, 0.00023150069103827686, 0.0003171883138609388, 0.0002929553419545396, 0.00048751000443365757, 0.00034827087026887327, 0.0006241784162421916, 0.000368875737833155, 0.0006511899764665889, 0.00038470812137795705, 0.0005580993302171582, 0.0003613775908410493, 0.0004006673557999964, 0.0003687678338341045, 0.0003917309541186645, 0.0004253522266026853, 1.6913597045635203e-5, 9.174658946158981e-5, 0.00018979119939338323, 0.00029924299316253937, 0.0003992096777182415, 0.00046923177406325806, 0.00046280501504571234, 0.0003951476207845053, 0.000394922282792041, 0.0004307772354034397, 7.393444672093994e-6, 1.2702102577644239e-5, 1.7042404275977357e-5, 1.9982190265088145e-5, 2.189669835266259e-5, 2.3015560617523906e-5, 2.3114864499728757e-5, 2.1722408849895718e-5, 1.9662821065757053e-5, 1.9313137566254865e-5, 3.574694087929181e-5, 4.5877927935487416e-5, 5.797030993117181e-5, 7.4434754922545e-5, 7.2087404463712e-5, 0.00010068264850209669, 7.894227974794404e-5, 0.00011895561563034872, 8.386309938674132e-5, 0.00012391553539483954, 8.96429963590888e-5, 0.0001147887480623884, 8.601201572943318e-5, 0.00010148725256459309, 7.469927868987284e-5, 7.917244521486891e-5, 6.742437930783743e-5, 6.60259190959837e-5, 8.573169196788499e-5, 0.00010086147503849125, 0.00011047280782122583, 0.00011828935740418825, 0.0001157917632565184, 0.0001098929914276632, 0.00010128219853628055, 0.0001030939305398462, 0.0001230385998877804, 0.00017025289027108574, 0.00015900568651430614, 0.000225402384897008, 0.00018388526660766786, 0.000265034371642331, 0.00020069134440412484, 0.000275744693771443, 0.0001984369216447328, 0.00026190981319165665, 0.00019688825692827658, 0.00022414933730287152, 0.00017245624125085697, 0.00017774763472718416, 0.0001604578190264321, 0.00019886442796249705, 0.00023148596879694167, 0.0002489291854845464, 0.00025572916453104136, 0.0002483757836182154, 0.00024116045671771787, 0.00021989878913274594, 0.0002515934442787987, 0.00023189090498593328, 0.0003369948759868514, 0.00027463098157738723, 0.0003967003232051097, 0.00029250663373756383, 0.00040918735645039474, 0.00029971032973394096, 0.0003830968168621212, 0.00028191519783746905, 0.00033215325720208696, 0.0002708114794909905, 0.00026132098201593515, 0.00024891410930796604, 0.00036995355521233717, 0.0003960686297174715, 0.0004127360163013731, 0.00040271339691074143, 0.00037699663471796904, 0.00033837928882089887, 0.0004367615755903434, 0.000348170487057004, 0.0005043460063358334, 0.00036820438795250764, 0.0005097453232149307, 0.0003601547036712761, 0.0004874367587615911, 0.0003682347427421079, 0.00040433714148208636, 0.00034127723121920466, 0.000335373381400618, 0.00032534432141597514, 0.000543430858028091, 0.0005371995761673537, 0.0005431703117717094, 0.0004989188544559907, 0.0004488909823981631, 0.0005819287640264218, 0.0004368071069286406, 0.0005986464212728102, 0.0004223873299188917, 0.0005533945427910939, 0.0004005433260293946, 0.00046806862078827265, 0.0004070382198588219, 0.00037117380036740037, 0.0003868962448410896, 0.0006253680817975281, 0.000602750909579928, 0.0005928560214556486, 0.0005287549926147999, 0.000661272339548215, 0.00047083722745431465, 0.0005939602694897986, 0.00044822806009663194, 0.0004812273642407062, 0.00043046401463153993, 0.00041359201307623634, 0.0004419676553291329, 0.00063056335375048, 0.0005990247730651117, 0.0005668790443570373, 0.0006227385608865493, 0.000467881868663277, 0.00048544366954072674, 0.0004542755225634419, 0.00042047873447552457, 0.0004665417475680599, 0.0005665825376099381, 0.0005358518175430689, 0.0004963756703822888, 0.00045063115665927047, 0.00041559782150091485, 0.0004691063991433209, 0.00047600575939381687, 0.0004376935741818605, 0.000513449627125142, 3.10072082009303e-5, 4.0355631752195054e-5, 4.64613966058361e-5, 5.01085748604645e-5, 5.058207879229078e-5, 4.981205684771526e-5, 4.631049524969775e-5, 4.2321601896264604e-5, 5.9966195059552035e-5, 6.920294116864251e-5, 7.239633859100445e-5, 7.6236207672248e-5, 7.478009220186292e-5, 7.02106739121592e-5, 6.400331809863122e-5, 8.508260153241094e-5, 9.724490119940859e-5, 9.472558526194885e-5, 9.554935071170533e-5, 8.711778662707237e-5, 8.336471556461486e-5, 0.00010457442195903162, 0.00010810864277210731, 0.00010742906868962859, 0.00010158047368292405, 9.292622821671022e-5, 0.00011725256005181206, 0.00011769500001244147, 0.00010867392749133607, 9.875208856503368e-5, 0.00011661725061328852, 0.000109170830350418, 9.948641665007307e-5, 0.00010502341246877764, 9.26041208167884e-5, 8.217186100399946e-5, 0.00013601080051418895, 0.00015715811275648724, 0.00016688361114123765, 0.00017106824308341724, 0.000169504858225154, 0.00015929822796308918, 0.00014401244203226393, 0.00019598973579698368, 0.0002133914055398549, 0.00021624468587197124, 0.00021138183453791218, 0.00020174180970877523, 0.00018380580138422324, 0.00024386593120814183, 0.00024409226509558088, 0.00024301988410061707, 0.00022801697212913622, 0.00020663931286647716, 0.00026877870874524593, 0.0002577222062437904, 0.00024457712448014993, 0.00022444003891337977, 0.0002664840596534675, 0.00024422084757282753, 0.00022304697648653258, 0.00023619677710331373, 0.00021198578026951345, 0.00019142198245010884, 0.00029513186056798625, 0.00031752500044335594, 0.00032631929901213346, 0.00031559859145547324, 0.0002972253234602571, 0.00026933401714650853, 0.0003639087258114058, 0.00037626258462802214, 0.0003622427127669694, 0.0003435866922522155, 0.00030658037687659816, 0.00039297411579118885, 0.00038732101169235215, 0.0003641899771513878, 0.0003275740150030641, 0.0003945961121397952, 0.00036671602420752455, 0.00033525801882466644, 0.00035580296238417475, 0.00031701425171807066, 0.00028433296694728265, 0.00046160793058577424, 0.00048374954834342764, 0.00046349578377012966, 0.0004299459501076951, 0.00039825914175206363, 0.0005135921021882168, 0.0004939892742221552, 0.000458753276373684, 0.0004144023131505741, 0.0004982777622755685, 0.0004656100084474529, 0.0004173674568657814, 0.0004491003781844341, 0.0003931924155577909, 0.00035898516599408196, 0.0005902187034681309, 0.0005761414005489118, 0.0005296672165671958, 0.0004951579914755016, 0.000568867191662482, 0.0005461139207048167, 0.00048322199809639984, 0.0005104360077936755, 0.0004487771617324709, 0.00040821672648732904, 0.0006321722511301668, 0.0005855819868628697, 0.000545045947260563, 0.000546921385550289, 0.0004975215626594298, 0.0004393034732055227, 0.0005810424805707117, 0.0005533102267042188, 0.0004587557042410595, 0.0005095621825637839, 0.0003310066675168652])
    cub_degree = 40
    tol = 5e-14
  end
  mask = 1:(cub.numparams+cub.numweights)
  vtx = T[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
  Cubature.solvecubature!(cub, cub_degree, mask, tol=tol)
  return cub, vtx
end

"""
### Cubature.getTetCubatureDiagE{T}

Returns a cubature rule and vertices for the SBP DiagE operators on tetrahedra;
these are cubatures that have nodes on the boundary that are analogous to LG or
LGL quadrature rules, which then leads to diagonal E.

**Inputs**

* `q`: maximum degree of polynomial for which the cubature is exact
* `T`: the data type used to represent the cubature
* `faceopertype`: the operator type on the facets of the tetrahedron
* `tol`: tolerance with which to solve the cubature

**Outputs**

* `cub`: a symmetric cubature for the right tetrahedron
* `vtx`: vertices for the right tetrahedron

"""
function getTetCubatureDiagE(q::Int, T=Float64; faceopertype::Symbol = :DiagE, tol=10*eps(typeof(real(one(T)))))
  cub_degree = q
  mask = zeros(Int64, (0))

  if faceopertype == :Omega
    if q <= 2
      # 13 nodes
      cub = SymCubatures.TetSymCub{T}(vertices=false, numfaceS21=1,
                                      centroid=true)
      SymCubatures.setparams!(cub, T[1/3])
      SymCubatures.setweights!(cub, T[2/30, 8/15])
      cub_degree = 2
    elseif q <= 4
      # 36 nodes
      # the boundary nodes are strictly internal to the face
      cub = SymCubatures.TetSymCub{T}(vertices=false, numfaceS21=2, numS211=1)
      SymCubatures.setparams!(cub, T[0.8918969818319298;
                                     0.18315242701954149;
                                     0.40398819659496876;
                                     0.18710499145686446])
      SymCubatures.setweights!(cub, T[0.02922867858673424;
                                      0.012914918864852366;
                                      0.06896751365952453])
      cub_degree = 4
      tol = 1e-14
    elseif q <= 6
      # 69 nodes
      cub = SymCubatures.TetSymCub{Float64}(vertices=false, centroid=true,
                                            numS31=2, numS22=0, 
                                            numfaceS21=2, numS211=1,
                                            numfaceS111=1,
                                            numS1111=0)
      SymCubatures.setparams!(cub, T[0.21759254267395756; 0.9034356836878749;
                                     0.1261780289830045; 0.49857349034182097;
                                     0.16668051173596254; 0.5518039588837356;
                                     0.10629009968963399; 0.6207049020675687])
      SymCubatures.setweights!(cub, T[0.020370182220389683; 0.08061671634468645;
                                      0.0027518135350883413; 0.01252025357308984;
                                      0.04988851015482234; 0.005364555718326347;
                                      0.01870947467719025])
      cub_degree = 6
      tol = 1e-13
    elseif q <= 8
      # 99 nodes
      numS31 = 1       # 4  symmetry
      numS22 = 1       # 6  symmetry
      numfaceS21 = 3   # 12 symmetry
      numS211 = 2      # 12 symmetry
      numfaceS111 = 1  # 24 symmetry
      numS1111 = 0     # 24 symmetry
      cub = SymCubatures.TetSymCub{Float64}(vertices=false, centroid=true,
                                            facecentroid=true,
                                            numS31=numS31, numS22=numS22, 
                                            numfaceS21=numfaceS21,
                                            numS211=numS211,
                                            numfaceS111=numfaceS111,
                                            numS1111=numS1111)
      SymCubatures.setparams!(cub, T[0.1836596600107025,0.8275896865477325,0.3411386155035204,0.9185851765854463,0.10109445663406191,0.10888330398707555,0.4087055731374049,0.4521142310502033,0.1419370312058845,0.5262256592692762,0.0167895548199152])
      SymCubatures.setweights!(cub, T[0.011702201915075788,0.007003762166809546,0.02826881609586945,0.005737494407876313,0.009215129850744767,0.001683701003412158,0.019901864184903844,0.0455776334197571,0.0008896293665706435,0.08215560123254996])
      
      cub_degree = 8
      tol = 1e-14      
    else
      error("polynomial degree must be <= 4 (presently)\n")
    end
  else
    if q <=1 # 6 nodes
      cub = SymCubatures.TetSymCub{T}(vertices = false,
                                      midedges = true,
                                      facecentroid = false,
                                      numedge = 0,
                                      numfaceS21 = 0,
                                      numfaceS111 = 0,
                                      centroid = false,
                                      numS31 = 0,
                                      numS22 = 0,
                                      numS211= 0,
                                      numS1111 =0)
      SymCubatures.setweights!(cub, T[0.22222222222222235])
      cub_degree = 1
      tol = 1e-14
    elseif q <= 2 #7 nodes
      cub = SymCubatures.TetSymCub{T}(vertices = false,
                                      midedges = true,
                                      facecentroid = false,
                                      numedge = 0,
                                      numfaceS21 = 0,
                                      numfaceS111 = 0,
                                      centroid = true,
                                      numS31 = 0,
                                      numS22 = 0,
                                      numS211= 0,
                                      numS1111 =0)
      SymCubatures.setweights!(cub, T[0.1333333333333337, 0.5333333333333319])
      cub_degree = 2
      tol = 1e-14

    # elseif q <= 3
      # #23 nodes based on 7 node facet operator
      # cub = SymCubatures.TetSymCub{T}(vertices = false,
      #                                 midedges = true,
      #                                 facecentroid = true,
      #                                 numedge = 0,
      #                                 numfaceS21 = 1,
      #                                 numfaceS111 = 0,
      #                                 centroid = true,
      #                                 numS31 = 0,
      #                                 numS22 = 0,
      #                                 numS211= 0,
      #                                 numS1111 =0)
      # SymCubatures.setparams!(cub, T[0.22222222222222215])
      # SymCubatures.setweights!(cub, T[0.01225481148524747, 0.02799913080121597, 0.16219856069976812, 0.27502065200818465])
      # cub_degree = 3
      # tol = 1e-14

    elseif q <=4 #23 nodes (for q=3 and q=4)
      cub = SymCubatures.TetSymCub{T}(vertices = true,
                                      midedges = true,
                                      facecentroid = false,
                                      numedge = 0,
                                      numfaceS21 = 1,
                                      numfaceS111 = 0,
                                      centroid = true,
                                      numS31 = 0,
                                      numS22 = 0,
                                      numS211= 0,
                                      numS1111 =0)
      SymCubatures.setparams!(cub, T[0.3771609693928902])
      SymCubatures.setweights!(cub, T[0.0017328144234984297, 0.010017745385041531, 0.07166219974832357, 0.40634920634920585])
      cub_degree = 4
      tol = 1e-14

      # #26 nodes based on 7 node facet operators
      # cub = SymCubatures.TetSymCub{T}(vertices = false,
      #                                 midedges = true,
      #                                 facecentroid = true,
      #                                 numedge = 0,
      #                                 numfaceS21 = 1,
      #                                 numfaceS111 = 0,
      #                                 centroid = false,
      #                                 numS31 = 1,
      #                                 numS22 = 0,
      #                                 numS211= 0,
      #                                 numS1111 =0)
      # SymCubatures.setparams!(cub, T[0.4842953457301756, 0.22222222222222215])
      # SymCubatures.setweights!(cub, T[0.18290740424677793, 0.022508956569022617, 0.01644731543334552, 0.06732054793298498])
      # cub_degree = 4
      # tol = 1e-14
    elseif q <= 5 #44 nodes
      cub = SymCubatures.TetSymCub{T}(vertices = true,
                midedges = false,
                facecentroid = false,
                numedge = 1,
                numfaceS21 = 2,
                numfaceS111 = 0,
                centroid = false,
                numS31 = 1,
                numS22 = 0,
                numS211= 0,
                numS1111 =0)
      SymCubatures.setparams!(cub, T[0.5008941915142769, 0.8506802519794945, 0.23722737279318576, 0.3077459416259917])
      SymCubatures.setweights!(cub, T[0.0015673886232196292, 0.17081759879508043, 0.033441261076507856, 0.01477813407660693, 0.00543005348522964])
      cub_degree = 5
      tol = 1e-14
    elseif q <= 6 #51 nodes
      cub = SymCubatures.TetSymCub{T}(vertices = true,
                midedges = false,
                facecentroid = false,
                numedge = 1,
                numfaceS21 = 2,
                numfaceS111 = 0,
                centroid = true,
                numS31 = 1,
                numS22 = 1,
                numS211= 0,
                numS1111 =0)
      SymCubatures.setparams!(cub, T[0.393148466877178, 0.25786391856319707, 0.8506802519794945, 0.23722737279318576, 0.3077459416259917])
      SymCubatures.setweights!(cub, T[0.0008311992138364736, 0.067074516283967, 0.08700624630786795, 0.025281032406421995, 0.01448700395673796, 0.0046508273850214945, 0.006646628516734983])
      cub_degree = 6
      tol = 1e-14

    elseif q<=7
      #76 nodes
      cub = SymCubatures.TetSymCub{T}(vertices = true,
                                      midedges = true,
                                      facecentroid = true,
                                      numedge = 1,
                                      numfaceS21 = 1,
                                      numfaceS111 = 1,
                                      centroid = false,
                                      numS31 = 2,
                                      numS22 = 1,
                                      numS211= 0,
                                      numS1111 =0)
      # SymCubatures.setparams!(cub, T[0.9465552688408668, 0.326390241974141, 0.7461816901961349, 0.160991838340075, 0.8003678928805428, 0.6058255660767269, 0.215183643569735])
      # SymCubatures.setweights!(cub, T[0.0008255000130940284, 0.01874509335958128, 0.06004170145017678, 0.0019891243235202693, 0.09315724230205255, 0.0030561455241973827, 0.0021889210067103236, 0.013884472504428686, 0.011959453952826883])
      SymCubatures.setparams!(cub, T[0.3649069878051798, 0.9069464284900626, 0.19709340992078445, 0.160991838340075, 0.8003678928805428, 0.6058255660767269, 0.215183643569735])
      SymCubatures.setweights!(cub, T[0.0003131195741790005, 0.06872580049477045, 0.0770182765903434, 0.001297874862616682, 0.05122298827420443, 0.004967708556549984, 0.001966184574009691, 0.013198465506081449, 0.00850236954064094])
      cub_degree = 7
      tol=1e-14
    elseif q<=8
      #89 nodes
      cub = SymCubatures.TetSymCub{T}(vertices = true,
                                      midedges = true,
                                      facecentroid = true,
                                      numedge = 1,
                                      numfaceS21 = 1,
                                      numfaceS111 = 1,
                                      centroid = true,
                                      numS31 = 2,
                                      numS22 = 1,
                                      numS211= 1,
                                      numS1111 =0)
      SymCubatures.setparams!(cub, T[0.26967567337250364, 0.9186565269236643, 0.7878184728172669, 0.16099183834007277, 0.800367892880542, 0.19236824251877926, 1.0726528991802897, 0.6058255660767288, 0.2151836435697346])
      SymCubatures.setweights!(cub, T[0.00017788261247842277, 0.022946869343355818, 0.05842047515941764, 0.0025968334566281045, 0.012039117612706244, 0.00471464641426865, 0.0014485783932065817, 0.041098991084533104, 0.010486947759748596, 0.01215056364118055, 0.051901126953532224])
      cub_degree = 8
      tol=1e-14

    elseif q<=9
      #121 nodes
      cub = SymCubatures.TetSymCub{T}(vertices = true,
                midedges = true,
                facecentroid = true,
                numedge = 1,
                numfaceS21 = 3,
                numfaceS111 = 1,
                centroid = true,
                numS31 = 1,
                numS22 = 1,
                numS211= 2,
                numS1111 =0)
      SymCubatures.setparams!(cub, T[0.2595307420691634, 0.8441512685303726, 0.40849548761177584, 0.1693510205976773, 0.895541492424858, 0.9192858150590572, 0.1362955833435948, 0.5043846222657248, 0.4701004391173425, 0.1604772007159291, 0.5199457984883094, 0.07662224108118755])
      SymCubatures.setweights!(cub, T[0.00015279459961010358, 0.01975401788383916, 0.001314202393676677, 0.021390096706851512, 0.006798995776920046, 0.0039159324894427845, 0.006885934192465157, 0.0005069027474549871, 0.020477522536811944, 0.04175001086337507, 0.0023605549505845433, 0.009102594817477005, 0.060393007434791986])
      cub_degree = 9
      tol=1e-14

    elseif q<=10
      #145 nodes
      cub = SymCubatures.TetSymCub{T}(vertices = true,
                midedges = true,
                facecentroid = true,
                numedge = 1,
                numfaceS21 = 3,
                numfaceS111 = 1,
                centroid = true,
                numS31 = 1,
                numS22 = 1,
                numS211= 4,
                numS1111 =0)
      SymCubatures.setparams!(cub, T[0.5700558294241873, 0.8819914623822424, 0.40849548761177584, 0.1693510205976773, 0.895541492424858, 0.9192858150590572, 0.3674432698076532, 0.12171080589025479, 0.7415577501465243, 0.1346989794996326, 0.10308419269846046, 0.5052883884978763, 0.11223871825337459, 0.203319497342745, 0.5199457984883094, 0.07662224108118755])
      SymCubatures.setweights!(cub, T[8.342292637650742e-5, 0.02995257400472967, 0.001178939330700515, 0.014024455736117564, 0.005665640945854503, 0.0025319848901773656, 0.005490413383261791, 0.0006215155149072649, 0.021063871513121817, 0.02981187131916029, 0.013879210846470981, 0.005258641145889365, 0.0016466522072965662, 0.0074534718610153455, 0.04075764008270213])
      cub_degree = 10
      tol=1e-14
      # #145 nodes
      # cub = SymCubatures.TetSymCub{T}(vertices = true,
      #           midedges = true,
      #           facecentroid = true,
      #           numedge = 1,
      #           numfaceS21 = 3,
      #           numfaceS111 = 1,
      #           centroid = true,
      #           numS31 = 1,
      #           numS22 = 1,
      #           numS211= 4,
      #           numS1111 =0)
      # SymCubatures.setparams!(cub, T[0.5831510926071354, 0.855295628963342, 0.40849548761177584, 0.1693510205976773, 0.895541492424858, 0.9192858150590572, 0.34277718068052015, 0.11217384327295098, 0.7459376222015264, 0.1271800044965199, 0.09162790449408253, 0.5994079435205514, 0.07787143633861693, 0.24259229663883067, 0.5199457984883094, 0.07662224108118755])
      # SymCubatures.setweights!(cub, T[0.00021816263981043935, 0.03921759527783465, 0.0011043100191260003, 0.009765154777214433, 0.005209712411026045, 0.0015442679048136543, 0.00532449820711722, 0.0003501468689026116, 0.022900629369277666, 0.0305175435641846, 0.013142498713694687, 0.006879720547610928, 0.0012142200424160243, 0.006991964941085668, 0.022836161062856575])
      # cub_degree = 10
      # tol=1e-14

      # #139 nodes
      # cub = SymCubatures.TetSymCub{T}(vertices = true,
      #           midedges = true,
      #           facecentroid = true,
      #           numedge = 1,
      #           numfaceS21 = 3,
      #           numfaceS111 = 1,
      #           centroid = true,
      #           numS31 = 1,
      #           numS22 = 0,
      #           numS211= 4,
      #           numS1111 =0)
      # SymCubatures.setparams!(cub, T[0.5803333201819036, 0.40849548761177584, 0.1693510205976773, 0.895541492424858, 0.9192858150590572, 0.330429102391519, 0.10730337252239175, 0.7539038145106844, 0.12433214114387407, 0.09670411438362254, 0.638953493007011, 0.06568839621878922, 0.23886759837771818, 0.5199457984883094, 0.07662224108118755])
      # SymCubatures.setweights!(cub, T[0.0003277864005158543, 0.04215578248651099, 0.000992857927615538, 0.00513620782204975, 0.001120039985381197, 0.004759201791331538, 4.93013126328509e-5, 0.022445330296939237, 0.03376662009972908, 0.015646220106639838, 0.006871141551676767, 0.0013019041677660326, 0.006840093117228722, 0.021307082127670487])
      # cub_degree = 10
      # tol=1e-14

    else
      error("polynomial degree must be <= 10 (presently)\n")
    end
  end

  mask = SymCubatures.getInternalParamMask(cub)
  append!(mask, (cub.numparams+1):(cub.numparams+cub.numweights))
  vtx = T[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
  Cubature.solvecubature!(cub, cub_degree, mask, tol=tol, hist=false)
  return cub, vtx
end

function getTetCubatureDiagELGL(q::Int, T=Float64; tol=10*eps(typeof(real(one(T)))))
  @assert(mod(q,2)==0, "Quadrature degree must be even.")

  cub_degree = q
  mask = zeros(Int64, (0))
  p=convert(Int, q/2)
  numedge = convert(Int,(p-mod(p,2))/2)

  # total number of nodes can be computed as:
  # n = p % 2 + sum(4 * (3 * ((p - 2 * i - p % 2) ÷ 2)^2 + 3 * (p % 2) * ((p - 2 * i - p % 2) ÷ 2) + (p % 2)) + 6 * (2 * ((p - 2 * i - p % 2) ÷ 2) + (p % 2)) + 4 for i in 0:div(p, 2))

  cub = SymCubatures.TetSymCub{T}(vertices = true,
                                  midedges = convert(Bool,mod(p,2)),
                                  facecentroid = convert(Bool, mod(p,2)),
                                  numedge = numedge,
                                  numfaceS21 = convert(Int,(1+mod(p,2))*numedge),
                                  numfaceS111 = convert(Int,1/2*(numedge^2 - numedge)),
                                  centroid = convert(Bool, mod(p,2)),
                                  numS31 = (1+mod(p,2))*numedge,
                                  numS22 = mod(p,2)*numedge,
                                  numS211= convert(Int, ((1+2*mod(p,2))/(1+mod(p,2)))*(numedge^2 - numedge)),
                                  numS1111 = convert(Int, (1/6)*((numedge-1)^3 - (numedge-1))))

  if q <=2 # 15 nodes
    SymCubatures.setweights!(cub, T[0.04848134481000349, 0.00011963663749342253, 0.16312953184124754, 0.4861720069033677])
    cub_degree = 1
    tol = 7e-15
  elseif q <=4 # 32 nodes
    SymCubatures.setparams!(cub, T[0.5199944491972114, 0.4257087142236166, 0.7236067977499789])
    SymCubatures.setweights!(cub, T[0.0031901936438669945, 0.15875259779608566, 0.043764724432287364, 0.013365456198839486])
    cub_degree = 1
    tol = 7e-15
  elseif q <=6 # 65 nodes
    SymCubatures.setparams!(cub, T[0.3778034871005324, 0.9186994092128735, 0.7952309920327962, 0.2891640824232583, 0.8750401312391274, 0.8273268353539885])
    SymCubatures.setweights!(cub, T[0.0009037335687971389, 0.06721830571168358, 0.06443360329716863, 0.0034238614529817873, 0.0550656183216355, 0.01538252157366906, 0.016103636737984227, 0.003588165446107915, 0.0011732397812532385, 0.02658904015688455])
    cub_degree = 1
    tol = 7e-15
  elseif q <=8 # 108 nodes
    SymCubatures.setparams!(cub, T[0.27533015932115545, 0.6219130396871575, 0.5306627609684195, 0.20735501628561037, 0.8825276619647324, 0.6426157582403226, 0.16950355726408317, 1.0748718697318453, 0.5077931286678652, 0.11291803951430215, 0.17654792120316234, 0.6492809445301031])
    SymCubatures.setweights!(cub, T[0.00028951549745859543, 0.027071133314835596, 0.04010512193137021, 0.004343361688919254, 0.005948822363053062, 0.001221435876252704, 0.0014021737097675549, 0.03428298138514228, 0.027369603517278764, 0.00702707116140467])
    cub_degree = 1
    tol = 7e-15
  elseif q <=10 # 175 nodes
    SymCubatures.setparams!(cub, T[0.2128785129992084, 0.9673676540009212, 0.5162909518831329, 0.8593193682825483, 0.8923932601020271, 0.733783076383442, 0.4127652030196255, 0.1577404405372036, 0.9383113442093055, 0.828770971730687, 0.9151119481392835, 0.7344243967353571, 0.12845378081353867, 1.2838828878296134, 0.4076479469911398, 0.07616274419943524, 0.7628782111929711, 0.1180398395722123, 0.5020748121582398, 1.3787314834940385])
    SymCubatures.setweights!(cub, T[0.00012587350281596975, 0.012244256930795668, 0.010856385828906114, 0.03437297074454886, 0.019933209306358283, 0.0004056776588533049, 0.012929763269911168, 0.012660585309612265, 0.0027277570409687943, 0.0028197770596327886, 0.003107921681073908, 0.003520838412260152, 0.0004675679444793231, 0.0005974016026886628, 0.01777124054365211, 0.01563787368302707, 0.018637092172174203, 0.0025937227508652174, 0.0014600888678078594, 0.015687047477891198])
    cub_degree = 1
    tol = 7e-15
  elseif q <=12 # 256 nodes
    SymCubatures.setparams!(cub, T[0.16714176919719123, 0.4376749036080209, 0.6655414830682306, 0.11946629335997022, 0.3416879952650081, 0.5747680456926114, 0.9358700742548033, 0.7958500907165711, 0.6046496089512394, 0.10143350406260081, 1.4348953166280882, 0.09291922438923644, 1.0995011713761764, 0.3267808232930995, 0.09059747390998994, 0.5504827645379977, 0.08338273047787416, 0.256637743345451, 0.8997479550832372, 0.5072568404058005, 0.23483688434081124, 1.5087470795599724, 0.10666597054658086, 1.15236027848101, 0.09573484698587671, 1.025654579962154, 0.30250431498538716, 0.6472795692230837, 0.985759608785055, 0.2932666748247045])
    SymCubatures.setweights!(cub, T[4.453736565450326e-5, 0.006294866838787169, 0.014468535490294646, 0.011369469267879198, 0.0011859157119268702, 0.0020397021188218708, 0.0019171330256122347, 0.00023188226221746013, 0.00029472347955574976, 0.00029494805484376105, 0.00854920439718218, 0.00869706490993752, 0.011186871763453005, 0.009079459537671438, 0.015771112345600732, 0.010976716961156688, 0.0016048978372864354, 0.0016107330379964728, 0.001670406019104794, 0.010194249883408871])
    cub_degree = 1
    tol = 7e-15
  elseif q <=14 # 369 nodes
    SymCubatures.setparams!(cub, T[0.1296892655169506, 0.9626988404922239, 0.35065610516353674, 0.929590842715833, 0.5930804702560417, 0.8453389071505838, 0.936282537613451, 0.8149234672156173, 0.6532780076474521, 0.9540069927604912, 0.09149922166977983, 0.5272375124699329, 0.8540637650586674, 0.29746375857253127, 0.7370347878007151, 0.94987899770573, 0.8385931397553689, 0.6815587319130891, 0.08700635067935361, 1.5425010805191477, 0.07963269048821799, 1.2672198623262216, 0.2860626424668289, 0.06714697397403416, 0.8566734539712058, 0.06596916988795214, 0.48453428504516755, 0.061562129559997854, 0.76316082556966, 0.05296246295844539, 0.21279873206383815, 1.0805122303883352, 0.43806340001431066, 0.18873206704283635, 0.7219064904355826, 0.17830187626501884, 0.08118879249452095, 0.6045976879703241, 0.2993309303393121, 0.09452178262313014, 0.5659292908830563, 0.2642822989570226, 0.5641483544960307, 1.135834147615453, 0.24861800080758142])
    SymCubatures.setweights!(cub, T[1.7069193136043886e-5, 0.003066150463265946, 0.0037466402760866716, 0.008740577821885484, 0.005814668907525145, 0.01112415029509867, 0.008279429053787412, 0.000173736117290934, 0.004207343465927521, 0.0067427753285156395, 0.009090625742179932, 0.0010594018449768324, 0.0005191448026041, 0.0010344762640983227, 0.0010988214803226801, 0.0009404827894709548, 0.0006283740888676858, 0.00011632365022825327, 0.00017028080303852386, 0.00021050723333464365, 0.004996724693852042, 0.005693312918251379, 0.00692111680933792, 0.005993728906829222, 0.006732184335287888, 0.004534141435385863, 0.011711051191065452, 0.009288944631538554, 0.00868579692008093, 0.0009157584201311722, 0.0010206163582135396, 0.0009570117877959086, 0.000535897946374851, 0.005412490226951414, 0.003389292964146634])
    cub_degree = 1
    tol = 7e-15
  elseif q <=16 # 500 nodes
    SymCubatures.setparams!(cub, T[0.10798541049256605, 0.3075779309172736, 0.5403701524697062, 0.7274856731293831, 0.07695246604559251, 0.23155195440690438, 0.4214410348442252, 0.5984805733893661, 0.9597669540832294, 0.8693869325527526, 0.7389624749052223, 0.5826394788331934, 0.07023829900189839, 1.623019813572776, 0.06518659285145366, 1.3889467233721564, 0.06028617750998111, 1.0986051159966757, 0.2287896125698148, 0.061760560505221024, 0.41618647382669427, 0.05495801808090896, 0.587558735478257, 0.03362800343153835, 0.1896550417731606, 1.2058246615376147, 0.18281129518289194, 0.9578607006416764, 0.3861348935917819, 0.18390888782000817, 0.562329439277371, 0.15403175956136744, 0.3107369671387058, 0.7952396519958825, 0.502794123403988, 0.2916032136531869, 1.6780802584948735, 0.070868673270381, 1.4310000152356053, 0.06475029714620971, 1.1307663966420336, 0.059871502734567955, 1.3223892173002267, 0.21047585794328214, 1.0527058565418306, 0.19351851169866596, 0.9390386879493979, 0.38284849575318963, 0.45896972131739533, 1.2736783661960833, 0.2105739512057541, 0.7363423349757549, 1.0203073240345906, 0.19025299751359556, 0.6679986539263091, 0.9194392430302105, 0.36953255942219404, 0.6335917110194244, 0.862667724645658, 0.37498707404759024])
    SymCubatures.setweights!(cub, T[1.1061003295760145e-5, 0.0017187261469557658, 0.005592196115798773, 0.005659205332586884, 0.0035484341356814184, 0.0003237766614955429, 0.0006883458400272007, 0.0008092190741809945, 0.00013348366900240352, 5.922514086723272e-5, 8.621535693661223e-5, 8.984934277903591e-5, 8.858358247904911e-5, 0.0028059651527839505, 0.003169275580873389, 0.003026039857225944, 0.003955386054627854, 0.004622744274278587, 0.0034041787508385527, 0.006472418887075142, 0.007358451103190142, 0.007758326308930371, 0.00437441329806314, 0.00819557030279324, 0.005814277934647648, 0.0004917744627490272, 0.0005552236822797203, 0.0005372843915624547, 0.000729543320470393, 0.0006702813998471054, 0.0007841504954145658, 0.004338561616787332, 0.0041244213273086, 0.0029875445308260784, 0.005963960119375802])
    cub_degree = 1
    tol = 7e-15
  elseif q <=18 # 671 nodes
    SymCubatures.setparams!(cub, T[0.08775494621928843, 0.9831537799293543, 0.27022817610911204, 0.949298614773821, 0.4655702533694297, 0.899666920361583, 0.6286877458943603, 0.8272600112737918, 0.9610118137945878, 0.8711002981987827, 0.7537566144330214, 0.6203459168442255, 0.06262408792230445, 0.19777933032042158, 0.37882457759149224, 0.5522678470067186, 0.7417773969871918, 0.8384879881534788, 0.9192489794614319, 0.9751369967510601, 0.9670007152040296, 0.8922417368315723, 0.7826176634981026, 0.6478790677934697, 0.05954192818223979, 1.6871080734634474, 0.05595787134445857, 1.4889093571399676, 0.05112887316168394, 1.237709023301744, 0.1946647014919956, 0.05586959787456791, 0.9124153813065555, 0.03974009731573857, 0.3637421060174282, 0.0479542092350331, 0.8417360583463076, 0.039429835537427696, 0.5260430787857253, 0.038600028073434045, 0.7621121572976974, 0.035378913129723244, 0.16482614506204957, 1.3083104021670124, 0.1553451132224468, 1.1024503201562428, 0.3420554048220757, 0.15185834747599677, 0.8078131718721219, 0.12578551544800565, 0.5015468572912704, 0.12410005364733051, 0.7259165335807718, 0.11433984116458397, 0.2831538590232135, 0.9162701881828489, 0.46192485193511656, 0.25121135935548583, 0.6792951941644226, 0.23330862099563554, 0.3969247802147992, 0.19029720152483445, 0.17280786654060148, 0.6444759719518841, 0.34728044663382485, 0.6010071331698775, 0.4182247647898829, 0.05843657476537568, 0.20577249022713864, 0.06088553891723073, 0.6831211550516773, 0.05282478749715142, 0.39305803844554077, 1.3778608367517136, 0.18177643769204604, 0.6355156526924747, 1.1520845758239402, 0.1637045336845858, 0.5880104153315449, 1.0446749055376783, 0.3243153051338426, 0.5608052216089431, 0.9900181655959259, 0.30831789590006403])
    SymCubatures.setweights!(cub, T[5.911486313202136e-6, 0.0009331876101695682, 0.002196020067304392, 0.003565036791969426, 0.0028746748348889104, 0.005763258943373063, 0.0036627237496656075, 0.005253109220909293, 0.004367633949792308, 5.1450408320104254e-5, 0.001368991669966781, 0.003280485945789359, 0.0038515707290594885, 0.004040021954782587, 0.00017502547320340754, 0.0004539501240271268, 0.0005874360894126949, 0.0004350633745699377, 0.0002675297617943632, 0.0004671417560834312, 0.0004278998138208395, 0.0002475591873914189, 3.2184458230159365e-5, 5.284562164558444e-5, 6.4504573967051e-5, 6.000991673531189e-5, 0.0016847030234403235, 0.002022281506178841, 0.0019785579604394226, 0.0026066772789108887, 0.0021945438027684157, 0.003176210006371219, 0.0023266577198767153, 0.0026884439957394865, 0.002249609409511851, 0.004528534358692779, 0.004195042346283999, 0.005067296691385673, 0.0037601636420916116, 0.00404972414040079, 0.0035861404532634995, 0.005125902066550957, 0.004708594768079708, 0.0047747905167315535, 0.0004602783945228078, 0.000530447363192729, 0.0005024487666867299, 0.00035874771931221957, 0.00029688737779064167, 0.00035087653443877846, 0.00017729947139099876, 0.0028059360471188394, 0.0027897843074485796, 0.002926700633652137, 0.004493367594528651, 0.0022711048029175596])
    cub_degree = 1
    tol = 7e-15
  elseif q <=20 # 864 nodes
    SymCubatures.setparams!(cub, T[0.07320987508233036, 0.23708956166533762, 0.42052720637628144, 0.5872167397348284, 0.7163530406203162, 0.05310181454886258, 0.16669909638072203, 0.31890625932107397, 0.4769062635958667, 0.6122291577074322, 0.9724496361114411, 0.9096396608220033, 0.8164380765159303, 0.6997654704826745, 0.5682764664274638, 0.05005680424615204, 1.7375892606156815, 0.04920633339168344, 1.56549404379135, 0.04518808801458637, 1.347489906719589, 0.03909183282301191, 1.0965672957046337, 0.16492578415912287, 0.049228666995134354, 0.32076753699349486, 0.03947471875412321, 0.4721330872067248, 0.036662099371676096, 0.6029728591269837, 0.023945558351164957, 0.14784094061655323, 1.3884934390862578, 0.13319844207875914, 1.2159832169078397, 0.12044787173888412, 1.002551152684079, 0.30335594825215423, 0.13403571321529564, 0.4505713028393479, 0.1173371559209213, 0.5829090958667188, 0.08990538092441049, 0.25111300133189446, 1.0343829219748817, 0.229204510800901, 0.8720545738724137, 0.427030452629397, 0.23238146219811348, 0.5559778001082458, 0.194781224633483, 0.35443065745948305, 0.7200232432430729, 0.517279834571761, 0.32704520834797957, 1.7744542331519155, 0.05080126431977686, 1.5955445904517396, 0.04801899242996552, 1.3704216313092366, 0.043836148777747434, 1.114531215089879, 0.03979079119920485, 1.5034294955202192, 0.15742672682529332, 1.2985603401515144, 0.14404642861985076, 1.062482496970096, 0.13060240593213884, 1.186451789763432, 0.29251018581533667, 0.9830892364157693, 0.2651997011457629, 0.8833988726389188, 0.43207240512538564, 0.3360625495213574, 1.4604614612527553, 0.15967305287498548, 0.5501397238283253, 1.2641179092130337, 0.1466873978275155, 0.7942010272112414, 1.0422064863366503, 0.12757015400694582, 0.5206730008120326, 1.1515717100349108, 0.2939140715957976, 0.7454109215898914, 0.9641157809556234, 0.25835904252523134, 0.6780158233738014, 0.8704829512360621, 0.4218001058096557, 0.5009293738147287, 1.108861076485701, 0.2747957406587572, 0.7167894205502473, 0.9264512410559093, 0.24786498672484442, 0.6543131262442629, 0.83854639439804, 0.4077546990918673, 0.6184255909495217, 0.7882072389858968, 0.3851804327426033])
    SymCubatures.setweights!(cub, T[3.457666136620096e-6, 0.0005446371070159081, 0.002402748443943679, 0.004300824211474075, 0.004704952064185493, 0.0018789616690896257, 0.0001049134370571845, 0.00029060596786208425, 0.0003581772371044784, 0.0003685228732702409, 0.00016460479355809155, 1.9268177091255335e-5, 3.081853438974878e-5, 3.663184583211621e-5, 3.6028789119741226e-5, 3.243108646460085e-5, 0.0010092271669264845, 0.0013360338061897863, 0.0013766942983641367, 0.0011345458698144453, 0.0016953375680144336, 0.0021710739093550417, 0.002068401925506505, 0.00130378578708043, 0.0030567197187023995, 0.003064631726539437, 0.0027897372081293844, 0.0036969857030823977, 0.0032815887116856815, 0.0024139998386612277, 0.004574944249237393, 0.0038629224859970655, 0.00425037201840505, 0.0027178949576739868, 0.003808137753399647, 0.002698317977458079, 0.0001778085831669669, 0.0002253029214597641, 0.00022378611159659504, 0.00019036533592946581, 0.00032144267197103976, 0.00031376869383474695, 0.0002842905074447705, 0.0003354411174404457, 0.000304188179093097, 0.0002870020109996699, 0.001964269438925013, 0.0018846731654573587, 0.0016516988146342186, 0.0019231865982985924, 0.0017677905379361383, 0.001670584389441466, 0.0034332195089700978, 0.002990374311124493, 0.0028024572865725855, 0.003621297466631965])
    cub_degree = 1
    tol = 7e-15
  elseif q <=22 # 1105 nodes
    SymCubatures.setparams!(cub, T[0.06267501653640925, 0.9886515835360495, 0.20389212670713358, 0.9651370128723687, 0.373157064885905, 0.9233972171118542, 0.5275534389398668, 0.8772046108153062, 0.657657726346992, 0.8154649233744519, 0.9697702902403434, 0.9040002171582837, 0.8148228213833496, 0.7069444236844128, 0.5994576122809905, 0.0452671188791389, 0.14251445923625328, 0.28073020926730075, 0.42127731315299516, 0.558591685525777, 0.9808784475256531, 0.941755226640668, 0.8924521179816877, 0.8344191924303693, 0.7628077735811749, 0.976654923321082, 0.9231737823259362, 0.8430942345408787, 0.741454910545668, 0.62464346505312, 0.04301634090076417, 1.7748620474802164, 0.04207730929415001, 1.626939041570611, 0.03877355397537654, 1.4373997170005115, 0.03466603521311862, 1.2136462483062873, 0.1420641499208082, 0.04178674422407022, 0.9363690885607485, 0.028755582519160505, 0.27884435962643006, 0.03405441050335307, 0.8854842409514766, 0.03018317172149535, 0.42126563334329903, 0.03189659881578633, 0.8228445230876009, 0.026593458903117703, 0.5544696138450975, 0.025864153268350308, 0.7495574761125045, 0.024255809905436375, 0.12877282765659626, 1.468755179463406, 0.1208720595693909, 1.3096316049910526, 0.11184975928907294, 1.1177123488419374, 0.2651425282457156, 0.11814893662330131, 0.8560683258053479, 0.09363931840918585, 0.4065931091494727, 0.10492788738867206, 0.79520708353927, 0.08732558482078812, 0.53257434528926, 0.08681819932276882, 0.7217404783687638, 0.08105523837760362, 0.2287208450485482, 1.127254975724832, 0.21273455669404578, 0.9738153864019438, 0.3831561340944387, 0.211108902571711, 0.7571338450637239, 0.1796857932385684, 0.511792567653549, 0.1748361449773906, 0.6882612490524181, 0.1643024395117324, 0.3219390745821825, 0.8365931906727629, 0.47510343303640834, 0.2948734812200478, 0.6484196606510297, 0.2757477477689316, 1.80711191161222, 0.04340190963100611, 1.652396972843757, 0.04108975716835016, 1.453269828612064, 0.04044688560715337, 1.2280021473835083, 0.033797064982277754, 1.5726788533048874, 0.13666382887880732, 1.390247092801606, 0.12988285384097822, 1.179582780981497, 0.11621906093358203, 1.2785482904192165, 0.2569247204326685, 1.0899451542969714, 0.24050646929898697, 0.9970083540636144, 0.3906008774317397, 0.29103737924286505, 1.532832670659594, 0.1374497154982262, 0.48056743352445297, 1.3571634594368764, 0.1263500987344382, 0.698926104984377, 1.1525424756859348, 0.114053480813203, 0.4557705190971788, 1.252069863445409, 0.25644685415315205, 0.6598266167093648, 1.0760648167043299, 0.23291864788864505, 0.6103836820869334, 0.9787563718318739, 0.38198780186195114, 0.43932324996167105, 1.2009475359460595, 0.2472496377426039, 0.6385792569992292, 1.0317625391751015, 0.22505241165885617, 0.5912869469615453, 0.9446905549979533, 0.3687632431542732, 0.5633693802375813, 0.8916647337534536, 0.3501548161000538])
    SymCubatures.setweights!(cub, T[2.0756550332744058e-6, 0.00034245763374506636, 0.0008873992068828519, 0.0015983007529660378, 0.0017305159207780033, 0.0032701495778036235, 0.0017267001132495543, 0.0032922880554504706, 0.0020340896796033742, 0.002691791941916502, 0.002257969684322491, 2.7359854094459054e-5, 0.000630280080315763, 0.0016760755697486022, 0.0021911358296262656, 0.002460672514284436, 0.0022456074876082634, 6.53287403238021e-5, 0.0001790923912460344, 0.00024390622395091544, 0.00026562140793698726, 0.00019090073993537762, 0.00012869166307239802, 0.00016629290361885592, 0.00020959848469178407, 0.00017826149468073884, 0.00018779930734935695, 1.1907154584507953e-5, 1.910948876134341e-5, 2.322957284858701e-5, 2.7438947366605114e-5, 2.061310981561832e-5, 0.0006399357460202354, 0.0008531177755975936, 0.0008882793781463526, 0.0008195338153678536, 0.0010863135849309188, 0.0009685486760516867, 0.0015004638303425487, 0.001191690027295624, 0.0016370066425377442, 0.001201292854553737, 0.0013012618751383528, 0.0011297049237715565, 0.0020191320916072607, 0.0021815710645412467, 0.0020784192826637723, 0.002617470335953165, 0.0018102843537313156, 0.0024969942460498994, 0.002010058644972833, 0.001971095670819797, 0.0018638477808491914, 0.0030810589806465646, 0.0027437369539764886, 0.0032413162770384032, 0.00243486462221746, 0.002425210176074058, 0.0023828811987207958, 0.0029367842634549267, 0.002991599353993106, 0.0025085098970065024, 0.00011254559226623966, 0.00014259439723228603, 0.00015428246355692017, 0.00012871153761765221, 0.00022042541943505798, 0.00021839722235679251, 0.00023905993999581945, 0.0002877788152764785, 0.00023469026506236078, 0.00025081762828592407, 0.0002637602106177805, 0.0013110556925019238, 0.0013614880886752693, 0.0012624295420573908, 0.0015933950242007866, 0.001484830485264371, 0.0013928276815738853, 0.0025430777885440637, 0.002319924493250731, 0.002306225979672744, 0.0028095731475255028, 0.001420091234291673])
    cub_degree = 1
    tol = 7e-15
  elseif q <=24 # 1372 nodes
    SymCubatures.setparams!(cub, T[0.054272317591233306, 0.17740394572978535, 0.3334815027774708, 0.4853460048942257, 0.6195204919782251, 0.717589084457966, 0.03880020746539761, 0.12469808500950672, 0.24546909601117012, 0.3824152837478771, 0.5133193869937424, 0.6213749471902272, 0.9799675226336304, 0.9339005269151737, 0.8644342995456631, 0.7753197014643236, 0.6713620066713564, 0.5581659344418519, 0.03709924918904916, 1.8052943958466718, 0.036688280399653966, 1.676094935128791, 0.03443744876776544, 1.5091403197289726, 0.030662377302941447, 1.3112100017550559, 0.027177985844068567, 1.089769080943189, 0.12298650555014924, 0.036250705046658026, 0.24508442187238352, 0.031627704098978734, 0.3791338328662843, 0.030243694992127687, 0.5084152081320412, 0.022611825649661004, 0.6144362889395394, 0.019387501043631912, 0.11440071411670867, 1.5318132361076024, 0.10782679179937568, 1.389197046242167, 0.10018749872379176, 1.216248110650283, 0.09472292841307627, 1.0149572143191454, 0.2352744945977336, 0.10659916575335687, 0.36599028259888844, 0.09544301626772941, 0.49029836256067627, 0.08080756548428675, 0.5995169329427684, 0.07199978244624564, 0.20629800677453525, 1.2160148796722325, 0.1939077515134735, 1.0716096563102089, 0.17714851949389207, 0.9155012689101941, 0.3484208923805293, 0.19162423671395706, 0.4705169598086495, 0.16704967599644213, 0.5726373358646704, 0.14324434319449741, 0.29927490546008206, 0.9193976033163173, 0.2766105977925864, 0.7973096124593809, 0.4448168108647823, 0.27872814187815864, 0.5510251423292437, 0.23721362434187473, 0.38090147940345787, 0.6744214959927355, 0.5118864047508731, 0.35316337433645584, 1.8338443095729768, 0.03780922263195278, 1.6994790270415137, 0.03640951307210373, 1.526293098464035, 0.03429184874986117, 1.3234587981554498, 0.031190632422757072, 1.1015360584910252, 0.027579295831113855, 1.6242840516199397, 0.12000396109006645, 1.4624546933179234, 0.11317147001885593, 1.2735879053298942, 0.10322501854397217, 1.0653951368868222, 0.09170711961572527, 1.3626513277391037, 0.2316819648099439, 1.194470962768185, 0.21217173127961658, 1.0086952027493412, 0.189774133359083, 1.0911070815110906, 0.35183079053218175, 0.9340684516928085, 0.3169984087390115, 0.8466357861808691, 0.46557399075669426, 0.25332462299295466, 1.5919735075581651, 0.12029167814204232, 0.4211248924074315, 1.4337297909782272, 0.1128381141793399, 0.6177161585439129, 1.2512783919175985, 0.10079470198977941, 0.8315073682904743, 1.047722551931521, 0.09119640368770532, 0.40381577914012157, 1.3363639046098073, 0.23060455188133236, 0.5892741108166348, 1.1737381847031887, 0.20756797938433696, 0.7930676594498164, 0.9940937583938299, 0.1893443724662453, 0.5517472009458536, 1.0761851178847928, 0.34702308863736164, 0.740086867878258, 0.9222056790681263, 0.31503746824850437, 0.6838781582231325, 0.8384062818908468, 0.46062635760620146, 0.39114675695221873, 1.289572924988927, 0.22113236499944103, 0.5704136529225154, 1.1301338363664983, 0.20396700547447255, 0.7694620901406508, 0.9620444180048763, 0.18676951014079043, 0.537032127420324, 1.0413552518855198, 0.3361330842633645, 0.7212335606195697, 0.8945261896459107, 0.30558504917002477, 0.6650452550569614, 0.8222076635647253, 0.45194186479777626, 0.5137744346351704, 0.9861685415478656, 0.3212461320206679, 0.6897534832584541, 0.8517380286265132, 0.2956991306255154, 0.6393759960089184, 0.7854532770633915, 0.4316679702213831, 0.6013774455924963, 0.7357834647744356, 0.410231435737946])
    SymCubatures.setweights!(cub, T[1.3203999664634947e-6, 0.00022233387074639324, 0.001075692369970965, 0.0023142997198701155, 0.0029486582542071174, 0.002489363837513305, 0.0011896708763397003, 4.16429006650473e-5, 0.00012059170622744767, 0.00017776276026942734, 0.00022214350504183713, 0.00015288158506796521, 8.716307122577691e-5, 7.4928816549237526e-6, 1.2524291868869741e-5, 1.5919966823808693e-5, 1.720044231522618e-5, 1.5925240008802963e-5, 1.2965751611201807e-5, 0.0004143354673086455, 0.0005651489947306139, 0.0006227821000846054, 0.0005683725379559052, 0.00048673274531808793, 0.0007144258198725049, 0.001078387313396924, 0.0012385485583481547, 0.0010321833310374588, 0.0007549015673624986, 0.001413161660232122, 0.0015507648228727765, 0.0015500311975714433, 0.0013764858621524625, 0.0018396076621712326, 0.0019339247819809716, 0.0017881434243964972, 0.0012283356568920725, 0.0024578339653224183, 0.002169907299296085, 0.002186009536427966, 0.00257614474707917, 0.0021258705491747822, 0.001538678782871474, 0.002608674408600378, 0.0023119770997932277, 0.0025493673938483835, 0.0016010099903913002, 0.002001551793076393, 0.0016073084697206565, 7.273896653016784e-5, 9.572754652904861e-5, 0.00010307461973454017, 9.43711901735695e-5, 7.85107447974832e-5, 0.0001500446009004079, 0.00016232185095297325, 0.0001572224644829598, 0.00014842732222625487, 0.00018685076662515319, 0.00019768847374829925, 0.00014682530469774238, 0.00017877916597577374, 0.00015124123860439214, 0.00013184102563748408, 0.0009075234077163679, 0.0009731765731832298, 0.000915565987434123, 0.0008589938879708946, 0.0011250539987449114, 0.001140734943646346, 0.0009143999024120266, 0.001113206285348649, 0.0009758902762408646, 0.0006472017187688928, 0.0019100838603704742, 0.0017933415846694601, 0.0015572566787899267, 0.0018607579201142912, 0.001627939181920845, 0.0015060242849763651, 0.0023823993234685257, 0.0019955259085156943, 0.002031392894923965, 0.0021691209455873446])
    cub_degree = 1
    tol = 5e-14
  elseif q <=26 # 1695 nodes
    SymCubatures.setparams!(cub, T[0.047542773969328475, 0.9922601533745666, 0.15574706872767344, 0.9720917065115615, 0.2971546461010158, 0.9436462541604224, 0.4428004024896039, 0.9041908947442963, 0.5660235517754111, 0.8609189428813461, 0.6681546230088244, 0.8061079220759985, 0.976950779086389, 0.9274466365692201, 0.8552558318675292, 0.7691513302549177, 0.6771929414949461, 0.5836636716082532, 0.03393923611308314, 0.10873649956316811, 0.21482988698982655, 0.3398694855014765, 0.4633755619477792, 0.5765169904802263, 0.9866066659382519, 0.957310865314053, 0.9167523374350691, 0.8686175130656929, 0.8131598018284998, 0.7482527988848714, 0.9826229632519192, 0.9425410221114882, 0.8817598449759076, 0.8031266027349229, 0.7103190273568363, 0.6076769776818971, 0.032325826672350774, 1.829800672303596, 0.03210254394837036, 1.7162604525908451, 0.030678559720125093, 1.5682070736376512, 0.028029380252420723, 1.3912129443503642, 0.02523007201098723, 1.1911626571193605, 0.10745389926403225, 0.03181176625283458, 0.9514380735216003, 0.021966583795962356, 0.21629171387046428, 0.02825219628568353, 0.9130640243494226, 0.02263877510547311, 0.34059597659624297, 0.025987478969196457, 0.8642837773789623, 0.02016936694715854, 0.4629617017859494, 0.023166376596852035, 0.8030411319436598, 0.01969722157292543, 0.5724501036270055, 0.019054340519405657, 0.734878819479639, 0.01826528916540262, 0.10137730511202311, 1.5866385371971037, 0.09647146259471681, 1.4583934102700549, 0.08991125649045013, 1.3028443348373524, 0.08419359217955509, 1.1212553797866338, 0.20924653207946806, 0.09466386593995839, 0.8899104727407678, 0.07189781287424638, 0.33027903814863563, 0.08452060215614, 0.841498923678987, 0.06861887978582411, 0.4489906898409905, 0.07491431549847598, 0.7848544571460777, 0.06347262067400022, 0.5572044641158977, 0.06467103187740449, 0.7184785814062606, 0.05969828084341939, 0.1860686826680774, 1.2949276287725129, 0.17719078313514186, 1.1577534861553191, 0.16301664558380513, 1.0137795850239864, 0.3140615016199303, 0.17463047157989287, 0.8094167657225045, 0.14073490417343504, 0.4342954106057527, 0.15308839520584575, 0.7541699239019988, 0.13269231595137798, 0.5404666821910308, 0.1328297600086847, 0.6948453150052952, 0.12454781706928476, 0.27586897903073115, 1.0070190661301286, 0.25719245121724194, 0.887886809222656, 0.40848399861481766, 0.25654067543778275, 0.7184374284712464, 0.22308473248749078, 0.5138124028081285, 0.21603567873654417, 0.6611162254858067, 0.20817885051653537, 0.35288552950582064, 0.773799260205942, 0.478757393158449, 0.32522903861140173, 0.6257287543920163, 0.30563060026650857, 1.8547247311426134, 0.03286746714508377, 1.7366818435509916, 0.031520503014385394, 1.5832211302981727, 0.029980899742798908, 1.4010761443451814, 0.028752652405317998, 1.2005374709609065, 0.025115089079989942, 1.6716587080221865, 0.10434769233046208, 1.5270795860587258, 0.09920116448130761, 1.354249854084208, 0.09440983064095608, 1.1644201761813233, 0.08468894753424339, 1.4380112429828442, 0.20422577449383864, 1.280888624803641, 0.19215819989412208, 1.1065564041141398, 0.17677134397464672, 1.1837996453343822, 0.3161530200288055, 1.0291162067290274, 0.29455381556320365, 0.9443063768751075, 0.4286992022860749, 0.22232303409636162, 1.641622783643567, 0.10563768544366371, 0.37173298208639455, 1.498627165007368, 0.10069858894588063, 0.5487981310096703, 1.3323500599939369, 0.09209132153675954, 0.7444448191076446, 1.146228179322952, 0.08355850025745597, 0.35835392632172075, 1.408944934815291, 0.20609157342436213, 0.5260297108326948, 1.2574993544288187, 0.18918951518821522, 0.7125891911832346, 1.0916195385639416, 0.17317253929116108, 0.4980987181177975, 1.1631801436254336, 0.31493351051413226, 0.6725463091395113, 1.0170323969157848, 0.2889202237475548, 0.6254486973343394, 0.9305791441930017, 0.4234579237891051, 0.34894792218250775, 1.3636912086739377, 0.19820695655286646, 0.5116949684573664, 1.2162173535307261, 0.18442210169921644, 0.694890039381511, 1.0573831907627456, 0.1703173202304459, 0.4859152704528159, 1.1269073985837414, 0.30620111047525533, 0.6564256069813624, 0.989327757104708, 0.2818135965959616, 0.6096846934899125, 0.9089521303541026, 0.41386705779012667, 0.4654193640312123, 1.0737809067093729, 0.29361664698719536, 0.6341050760722612, 0.9424532970191494, 0.27073347469145537, 0.5893937243657038, 0.8732414158656747, 0.3986827745345298, 0.5610610450009167, 0.8255933271748032, 0.37756252746516017])
    SymCubatures.setweights!(cub, T[8.627919873415628e-7, 0.00014941246102938215, 0.0005305650674186445, 0.0007323573086667097, 0.000881744360890124, 0.00169639887721994, 0.0011424166567645835, 0.0022795711815730858, 0.0012124366312364205, 0.0018587173068908183, 0.0012423267875998735, 0.0014909394474546351, 0.0014754550496454983, 1.1886877195664484e-5, 0.0003237703187209379, 0.0008595304708178885, 0.001327209844885543, 0.0014754140105383215, 0.0014359038113490555, 0.0014489046012772682, 2.7925942282198708e-5, 8.11939552194382e-5, 0.00012569451166091585, 0.00015785534814169574, 0.00015243916729601963, 0.00011381517054105804, 5.922369496445363e-5, 9.147978055208415e-5, 0.00011101847515646787, 9.757894526254777e-5, 0.0001138440385729111, 0.00010552504814237452, 4.97471368453526e-6, 8.242997170469378e-6, 1.0423675157149577e-5, 1.1676759327350718e-5, 1.2243131690800659e-5, 9.817303803472048e-6, 0.00027625673950643347, 0.00038086401724069945, 0.00043763082165022984, 0.00042468167358269103, 0.00037817526255198825, 0.0004826123780644111, 0.0004794711179465042, 0.0007632615386268363, 0.0006294756637878089, 0.0009014162422327122, 0.0006931165959535975, 0.0008687513374059021, 0.0006718089774743637, 0.0006626288601368487, 0.0006028126809699694, 0.000990906763214186, 0.0011110634425075792, 0.00115201290069364, 0.0010696231391070132, 0.0013119936933846449, 0.001003192410971056, 0.0015102645342577484, 0.0011636014309895513, 0.0013881597610060183, 0.0011155198899406148, 0.0011344778378863459, 0.0010397453088579255, 0.0018156095296418871, 0.0016708656878969253, 0.001637596263131074, 0.0020137845779386668, 0.001393456844779394, 0.001830116655170204, 0.0014786274062824445, 0.001396251041457776, 0.0012933543169207334, 0.0020007928290862688, 0.0018095578938466998, 0.0020398587187918796, 0.0015376226396823095, 0.0015489172346277865, 0.001462182527480776, 0.0018078577402601474, 0.0017063891937952897, 0.001658089199013831, 4.831468315445846e-5, 6.405824681593872e-5, 7.133630220154571e-5, 7.118386326923511e-5, 6.114079702205854e-5, 0.00010223787219128464, 0.00011487530280338325, 0.00011297134620580652, 0.00011067815239644849, 0.0001366149865370557, 0.00015122463793903476, 0.00012372460068037786, 0.0001416371257435199, 0.00013111248713929165, 0.00012375819481195223, 8.616975456074076e-5, 0.0006289823849729984, 0.0006984274403166871, 0.000678310737625329, 0.0006310737243428655, 0.0008343034991250964, 0.0008735252310700412, 0.0007473638706176485, 0.0008961768464056792, 0.000779832938868989, 0.000743329360489549, 0.001412061717597786, 0.001387476607657943, 0.0012704340619455231, 0.00150920985192356, 0.0013245562672888754, 0.001233522897691453, 0.0019275204740899674, 0.0017336053491906014, 0.0015896133041380057, 0.0018471517362639981, 0.0009746576395282034])
    cub_degree = 1
    tol = 5e-14
  elseif q <=28 # 2048 nodes
    SymCubatures.setparams!(cub, T[0.04158870045717802, 0.1392091338564074, 0.26927612518805905, 0.4032069171280793, 0.5320355246602935, 0.6437728335384955, 0.7236171210643413, 0.0297684342510763, 0.0955990655973306, 0.19291754638860226, 0.3071372298707168, 0.42368986691738036, 0.5320960490073281, 0.6246702357765312, 0.984784023135109, 0.949600266546736, 0.8960041459309076, 0.8261943514412465, 0.7430297109435688, 0.6499152344503816, 0.5506631367609747, 0.02845677921001763, 1.850609865803159, 0.028408895820192084, 1.7503585980993346, 0.02729416042431138, 1.6190444910495994, 0.025340355608688303, 1.4603001077920041, 0.02306231693450154, 1.2789969037649007, 0.021094120455130688, 1.0815717808766856, 0.0947048293055917, 0.028237956621989115, 0.19277537269335251, 0.02568982882415432, 0.3067896862019649, 0.023760561474237062, 0.42203018682091364, 0.021993056221571466, 0.5322226256876401, 0.017646810855223825, 0.6223977817809468, 0.015438174102042117, 0.0900886905557347, 1.6305796695880115, 0.0861917952942873, 1.514236084753114, 0.08127949490828822, 1.373940927999722, 0.07604255227333451, 1.2094473085670387, 0.07042825460421237, 1.0273131730771814, 0.185882692697119, 0.0863865455068635, 0.29712863273732687, 0.07765314466730477, 0.41361801617749105, 0.07035543748422751, 0.5201204155766331, 0.060444428464199264, 0.6104313262441384, 0.05194212702875735, 0.16802427588164717, 1.363308240539159, 0.16108298159293888, 1.236468451370892, 0.1492394375634891, 1.100323834651408, 0.13735268338522422, 0.9489227179266444, 0.2853732753080144, 0.1590977400437381, 0.3994180652084759, 0.14204973975457427, 0.5043491135122099, 0.12517009179387314, 0.5936113922637425, 0.11029853517292955, 0.2540938525576177, 1.0895868739164292, 0.2375766637517794, 0.9747525551656407, 0.21920658446418073, 0.854527524933342, 0.3777442181194952, 0.23759246472723106, 0.48205478135817836, 0.20800968413940946, 0.5690384205913904, 0.184280018002248, 0.3329976139098257, 0.8463903655718918, 0.30945655443685555, 0.75094905217955, 0.4564689761658077, 0.31169412010503617, 0.5440408053821607, 0.27041016880950486, 0.3995088082620552, 0.6495180742744258, 0.5123759428918147, 0.3707687599302176, 1.872342032026377, 0.028871097803943917, 1.7675811074108279, 0.028158691401833948, 1.6306934015947248, 0.02719784703705251, 1.4677125664505128, 0.025321800402231868, 1.2846879400321705, 0.02347760065222769, 1.0891749012037502, 0.022023861231761973, 1.7093628889667856, 0.09344851689576034, 1.5801489061325036, 0.08990447329553126, 1.4262344614567835, 0.08409466147453422, 1.2515984975183934, 0.07779747113489333, 1.0624982098653115, 0.07221835901759747, 1.4931838100842858, 0.18443446862241503, 1.3509572181334975, 0.1729159926865194, 1.191532901609541, 0.1600935095819077, 1.0192254255841438, 0.14693939098728045, 1.2587688670743364, 0.2877527514030317, 1.1161287832133815, 0.2667227445499398, 0.963965368480241, 0.2427509005204612, 1.0307078228125615, 0.39285142735909817, 0.8993858298065905, 0.35757238978691713, 0.826436434580063, 0.4877366783995886, 0.19648204320532384, 1.6825035800088408, 0.09381786912282739, 0.3300452513440478, 1.5544408095879079, 0.08977504639167833, 0.48991838351656175, 1.402346887980505, 0.08336938179388355, 0.6688918612256559, 1.2322787258540262, 0.07603016653538328, 0.8587487554842956, 1.0495493452597453, 0.0700629165552182, 0.3206249505356876, 1.4708887755807216, 0.184512834623024, 0.4722644371279332, 1.3324926901400342, 0.171653148469879, 0.6432411231148438, 1.17739866575729, 0.1574568049185672, 0.8269060721272822, 1.0085560822164283, 0.14494972521757973, 0.45132576952173115, 1.2406125803752732, 0.285812462121139, 0.6116150559097279, 1.103762598848345, 0.2643889553148952, 0.7851432120970435, 0.9547648159474906, 0.2418482282514706, 0.5741893430692758, 1.016032564933992, 0.39179875571208206, 0.7367136092323692, 0.8885951534228153, 0.3569938684846282, 0.681959168144022, 0.8148406762661212, 0.48691931821462864, 0.31183103145851104, 1.4299876352341703, 0.17806224589573175, 0.46128607525869364, 1.2933141136839224, 0.16734514788623475, 0.6281280915154388, 1.1439574574846467, 0.15576669669987867, 0.8074038363365205, 0.9845851217397601, 0.14240406903404035, 0.43909769573697016, 1.2060541359109587, 0.2806314125638258, 0.5999445578283571, 1.0740415217494326, 0.25890874409994663, 0.7701262906868256, 0.9312895625226492, 0.23640913764166338, 0.5625481488344894, 0.9928958450085585, 0.382784929544118, 0.7212219417010991, 0.8698652223819244, 0.3503922952109442, 0.6682878056913185, 0.7988782558166245, 0.4785741057459089, 0.4245371338531299, 1.15379125849042, 0.26891574854266626, 0.5801311859496757, 1.0304297414337051, 0.24928291689438542, 0.7435686249225335, 0.8967619615434923, 0.22999665648093337, 0.5449137398734734, 0.9542855299376315, 0.36995161064665283, 0.6983914409156355, 0.8387079189278283, 0.34002486311119956, 0.6482441402024439, 0.7738338273737353, 0.4647991830006785, 0.5201961939374949, 0.9040376111739913, 0.35296674261690797, 0.6686292710988843, 0.7978088836178812, 0.326238808320408, 0.6208061099653465, 0.7390163625403874, 0.4466285975908461, 0.5891333581237516, 0.6974941828691256, 0.4236723603264375])
    SymCubatures.setweights!(cub, T[5.767984612051227e-7, 0.00010021963130054694, 0.0005366731844440196, 0.0012331518598919522, 0.0017603659254977146, 0.001796097372428086, 0.001373273231291659, 0.0007554744665503432, 1.8819337304816677e-5, 5.562744376651123e-5, 9.26735117512961e-5, 0.0001181737209268278, 0.0001239718552667664, 9.063969937820853e-5, 7.238096067833714e-5, 3.3478524917876287e-6, 5.585810991421275e-6, 7.365982152174759e-6, 8.560272906984345e-6, 8.506220628458909e-6, 8.004129379537835e-6, 7.347205869821417e-6, 0.00018805384713014306, 0.00026423863066709697, 0.00030809623066520536, 0.0003127140587693477, 0.00028754791289668487, 0.00025492960144164316, 0.00033558506102010055, 0.000561487976401711, 0.000698527719489261, 0.0007234399272098342, 0.0005734667387989546, 0.00039503420255281247, 0.0007019290104789134, 0.0008038656185034544, 0.0008404335955696766, 0.0008299183377941924, 0.0007213668478555634, 0.00095478801569021, 0.001141203868953517, 0.0011268372697789918, 0.0009989158788810413, 0.0007083004806839604, 0.0013719603572986806, 0.001343324419315638, 0.00127542327641519, 0.0011817406941196835, 0.0015503466723833619, 0.001494519611283823, 0.001249088733777519, 0.0009138375840632124, 0.0016675638073368015, 0.0016027633058242183, 0.001430532061838893, 0.0017540923290007664, 0.001426222046829823, 0.0010016371429113617, 0.0016611906177421486, 0.001425794886254029, 0.0015825764532400872, 0.0010214951020202122, 0.0012328934636477006, 0.0010014923837542355, 3.296144398097881e-5, 4.480020140522377e-5, 5.1432944398440226e-5, 5.129912209639039e-5, 4.726748811399397e-5, 4.223633562689029e-5, 7.296091977444859e-5, 8.176972366720507e-5, 8.480467933179649e-5, 8.051205630508531e-5, 7.407367007646002e-5, 0.00010153500491617255, 0.00010721827815071764, 0.0001008748905523194, 8.447093511303674e-5, 0.00011613747274447413, 0.00010660273198330812, 8.752576310108459e-5, 9.619576922901025e-5, 9.475080770030532e-5, 8.421741934332235e-5, 0.00044475697588659303, 0.0005027854607238196, 0.0005077824126276755, 0.0004720852733924957, 0.000425673710258127, 0.0006110015387785212, 0.0006394899566120506, 0.0006005629943681483, 0.0005157518652257133, 0.0006981323646314305, 0.0006369197523112427, 0.0005597806990750513, 0.000611771660653284, 0.0005491099235064752, 0.0004933049427631179, 0.001049336669606262, 0.0010836369051484109, 0.0009835299023527471, 0.0008653983515923574, 0.0012067478573646185, 0.0010751409099962281, 0.0009475329609903365, 0.0010633919298325895, 0.0009500833695059442, 0.0008336568397899625, 0.0015242733204674354, 0.00142434018385195, 0.001218781228500237, 0.0014211350044032857, 0.0012360769109664388, 0.0011278236828327728, 0.001646129271275087, 0.001370986328802364, 0.0012904095091306217, 0.0012951919271860738])
    cub_degree = 1
    tol = 7e-15
  elseif q <=30 # 2465 nodes
    SymCubatures.setparams!(cub, T[0.037048223259546195, 0.9933426144021754, 0.12333194474938634, 0.9786504301490093, 0.24055279981822036, 0.9549836289235876, 0.3699205805860287, 0.9249158266942779, 0.493358076569877, 0.8882677208324776, 0.595858739063228, 0.8481857806501849, 0.6824847618563065, 0.7998204931958816, 0.9821303868006912, 0.9428466550949708, 0.8846725372979491, 0.8124965015612639, 0.7359957118960878, 0.6546423904806773, 0.5730186315677996, 0.026263712933993107, 0.08517773895617918, 0.17279615006913418, 0.27589945905594404, 0.38841435055824713, 0.4974390314519614, 0.5848434660310842, 0.9888569686145411, 0.9661320178091347, 0.939352201678907, 0.9019197972335102, 0.8491616753371906, 0.7917391940597742, 0.73363931848329, 0.9865660883157092, 0.9554399979577868, 0.9078481256108851, 0.8455144903138423, 0.7706926996650507, 0.6860872167827385, 0.5947559867591588, 0.025217944803116732, 1.8671914605144386, 0.025224197566544437, 1.7777746812989759, 0.02441405182375261, 1.6601755408213765, 0.022907249573853525, 1.5173034076153527, 0.021084917555640688, 1.3530254630108385, 0.019429287422530677, 1.1722173026395653, 0.08404028534849839, 0.024989105772438378, 0.9625886815266075, 0.01647647827677876, 0.17207262729953338, 0.022792071179089306, 0.9317324174277491, 0.017573523364220935, 0.2772670400323661, 0.021189187909605955, 0.8922514239776869, 0.016310156875218997, 0.38552734300526936, 0.019817523822734908, 0.8436394700699704, 0.014906937393683374, 0.4928672114002477, 0.0172014696841062, 0.7880049522205852, 0.01430150904416084, 0.5873062368419248, 0.014858409810738316, 0.7278117776617278, 0.013706470991133318, 0.08087904850447089, 1.670232985172078, 0.0778162921490587, 1.5650977086650268, 0.07379996494148194, 1.4374484268463736, 0.06960075485683571, 1.287155445108531, 0.06435123016736963, 1.1204536816843598, 0.16739226713339367, 0.07674561460936621, 0.9125236472613715, 0.056669417421042985, 0.2688379383836915, 0.07079190521616503, 0.8743946967117932, 0.054253782742413104, 0.3799655491384078, 0.06510053719638796, 0.8272584470222711, 0.051463870094896426, 0.48307639788109036, 0.056225426886057105, 0.7725748632813912, 0.04903288356206592, 0.5760146167786909, 0.04869154528561676, 0.7167854406599811, 0.047096421240653295, 0.15409860979087833, 1.4219043774258937, 0.14693926271857558, 1.3055352453987517, 0.13759266957440613, 1.1764032456071416, 0.12710378666822408, 1.0365307301010855, 0.2600016290445755, 0.1458794912398464, 0.8476366507630012, 0.11292232443984694, 0.36592375700385404, 0.13150654627504382, 0.8024876325428267, 0.10807674591408868, 0.4690137913028312, 0.1160393075881066, 0.7526005368683352, 0.10283059059372435, 0.5589942538068315, 0.1016318305303808, 0.6957990750662147, 0.09752836479962287, 0.2342693301043778, 1.1617871643391358, 0.22220654228924183, 1.0504598383974204, 0.20541159998457592, 0.9384515310972853, 0.34876509589350335, 0.22034891413098917, 0.7712508459279767, 0.18215023020299706, 0.4495084984823967, 0.19371850774905303, 0.7230024740349676, 0.17128972769314463, 0.538012229799591, 0.1704739693446314, 0.6717786559483014, 0.1615205152875771, 0.3118010040824791, 0.921523170846175, 0.2904558785954642, 0.8270927596650787, 0.42589198138596057, 0.28956455251309954, 0.6891390172691241, 0.2568228858236239, 0.5113230703638945, 0.2503950152269634, 0.6403843241103415, 0.23967105737194722, 0.37568761801989803, 0.7338291154990482, 0.4814704449476996, 0.34548913080574134, 0.6086319296888677, 0.330336433301088, 1.8870221434136396, 0.02570366713733081, 1.793983492536147, 0.025119958538073715, 1.671668422695357, 0.02432310058038001, 1.5249412598428131, 0.02255867065673869, 1.3578052693232254, 0.021693571868118818, 1.1783277514824677, 0.019013639177892776, 1.7405684419522667, 0.0834737046514834, 1.6245714394071344, 0.08045190568031999, 1.48581264880407, 0.07531756667200452, 1.3259706249257068, 0.07132176062782467, 1.152384095100876, 0.06529596804242158, 1.544801702685509, 0.1652266243762263, 1.4137371720619094, 0.15610982592370423, 1.2659378667488168, 0.14566531636274646, 1.1041469369299484, 0.13811994187867957, 1.329023725898816, 0.26259117857381364, 1.1969481364030816, 0.24335464892304964, 1.0461028692407564, 0.23226408985909297, 1.108279879880839, 0.3623302386259286, 0.9830644472062909, 0.3357220260332416, 0.9068645543524504, 0.4456950580906604, 0.17477161422093132, 1.7173718796757236, 0.08348019932656423, 0.2945427133369693, 1.6016690232789963, 0.0805382712304214, 0.43907225255919485, 1.4633971790110152, 0.07542471778990997, 0.602642589891286, 1.3071005827903217, 0.06963773179350687, 0.778727644412937, 1.1369832928338124, 0.06443803835575168, 0.2872937582920823, 1.5241964751774404, 0.16611987185410662, 0.4255643309850218, 1.3975031804118185, 0.15526505710448132, 0.5821638852543752, 1.2530280185978537, 0.14464071050782532, 0.7523259133121397, 1.0962211127225419, 0.13338120636293693, 0.4096767031899597, 1.3110139303390507, 0.2593715516568011, 0.5577620397145326, 1.1813825845625938, 0.2429357619792322, 0.7190504848792741, 1.0398470949639693, 0.22420333741441306, 0.5274461916928402, 1.0946075256204277, 0.36004999871114024, 0.6792293492958136, 0.9707280344476192, 0.33296061290746093, 0.6356554949281281, 0.8945763025635896, 0.4545431397663953, 0.28105130826113206, 1.4845211369947724, 0.16041201901001972, 0.4168136707953892, 1.3587493161818829, 0.15268130184537013, 0.5705183607793907, 1.2205667496381183, 0.14296491291400393, 0.7370070605001497, 1.0706792017706455, 0.13129299394714952, 0.39916879927675747, 1.277183378011966, 0.25730846793633533, 0.5471907246515181, 1.1522704124365952, 0.23849434282666607, 0.7055680290086259, 1.0166019340031802, 0.22024908073322805, 0.5182570784803734, 1.0702543407402798, 0.3530497916498024, 0.6666488020392265, 0.951173109962792, 0.3273425586904497, 0.6237603092936128, 0.8780290589236487, 0.4464701742746096, 0.3869003997116288, 1.2284432648917127, 0.24515358516566438, 0.5303769398335952, 1.1077869695484166, 0.23099523985633405, 0.6861479117635724, 0.9781509505570543, 0.21422896801871738, 0.5033572662323721, 1.031831150792892, 0.34299092309020746, 0.6483655059334752, 0.9193861056822852, 0.31858547473647914, 0.6064215985800504, 0.8531486913063959, 0.43470563481266905, 0.48108661308708595, 0.9828593471069001, 0.32875130294337274, 0.6238341686214706, 0.8762126416856025, 0.3060910031304059, 0.5845283520166847, 0.8192913658125607, 0.41760048926311183, 0.5566992616489097, 0.7774267509451381, 0.39736925106238014])
    SymCubatures.setweights!(cub, T[3.965461741034693e-7, 7.087906212211371e-5, 0.0002988022327901447, 0.00037538750814287823, 0.0005672343402050146, 0.0009077088990181246, 0.0007341375058349935, 0.0013859299684558372, 0.0007925368823683418, 0.0014719368900558642, 0.000861996381210137, 0.0010070789275717408, 0.00084412825163972, 0.001013463720431327, 0.0009681874548699814, 7.336720987205299e-6, 0.00017107247871097893, 0.0005057113974688011, 0.0007759215289986043, 0.0009633351232703756, 0.0009366770749093939, 0.0009700065825003335, 0.0009353289624214487, 1.3048392049330872e-5, 3.915505968413717e-5, 6.653718408812856e-5, 8.513708528873408e-5, 9.389157955924884e-5, 8.937418803679832e-5, 7.36518771027195e-5, 3.284797553186828e-5, 4.120924740569072e-5, 5.277328187219344e-5, 6.84151840059153e-5, 7.081900947925282e-5, 5.7272586454322424e-5, 5.990369588826064e-5, 2.3045984712680876e-6, 3.916761879517422e-6, 5.216260960951863e-6, 6.117164323124529e-6, 6.0730288767240855e-6, 6.285658356716969e-6, 4.810000064658842e-6, 0.00013165913342551614, 0.0001860782646533145, 0.00022089655320039597, 0.00023041417921993668, 0.00021852460314393255, 0.0002001855994919186, 0.00023505421840595885, 0.00026098804289542335, 0.0004030119111010987, 0.00035490183979957013, 0.0005268021041981447, 0.0004018037171620031, 0.000585287410505873, 0.0004134844734119369, 0.0005064013731564314, 0.0003817537596937393, 0.00038861428683632424, 0.00034228528898161284, 0.0005088780215100273, 0.0005895934049477864, 0.0006342917480091036, 0.0006443216159056995, 0.000567286063250928, 0.0006969163765858905, 0.0005841554822155197, 0.0008889389473557324, 0.0006795318171380688, 0.0009180997442508423, 0.0007324115543499704, 0.0008451010796384262, 0.0006767414840192419, 0.000672647599326924, 0.0005938737109549258, 0.0010434533812490557, 0.0010205285045306748, 0.0010396881362962605, 0.0009447888286714825, 0.001168277535401594, 0.0008689675135635068, 0.0012341774797674987, 0.0009276719777834636, 0.0010736703503734805, 0.0008572011893275781, 0.0008543826157149001, 0.0008118253711252402, 0.001349046504154141, 0.0013111149207844929, 0.0011924884031559732, 0.001424448918134951, 0.0010094285055937955, 0.0012503753559140304, 0.0010224681881604773, 0.0010002137888524257, 0.000916655465737192, 0.00144254381746779, 0.001247867858802824, 0.001260737646049255, 0.0010618384362525244, 0.0010139280867645796, 0.0010148769966689094, 0.0011951212656369312, 0.0010372752639633788, 0.0011594030719273398, 2.3023697925639468e-5, 3.1576817241222676e-5, 3.694287434917333e-5, 3.7483750987849144e-5, 3.6451929459996667e-5, 3.1857145773822446e-5, 5.2360249030470636e-5, 5.934129216679983e-5, 6.303277388459008e-5, 5.976172437841223e-5, 6.0641820332334563e-5, 7.805500363715648e-5, 8.239905487586215e-5, 7.872957073458994e-5, 7.005838014017106e-5, 9.030796223979015e-5, 7.875928846404414e-5, 7.428023984246964e-5, 9.018684926387437e-5, 7.571011415031929e-5, 6.454974943671703e-5, 5.987211731503044e-5, 0.00031823806097146885, 0.00036872608687070637, 0.00037885948564400336, 0.0003655605095373658, 0.00033208855373809435, 0.0004670652954950174, 0.0004910747608720107, 0.0004670936587521292, 0.0004184670457148038, 0.0005278080279048188, 0.0005028637676577212, 0.0004592419701721722, 0.0005186839069584399, 0.00046868121560022087, 0.00043857605458185074, 0.0007864525702217564, 0.0008368900398998376, 0.000776967162751651, 0.0007096150086018871, 0.0009372504200026877, 0.0008870048088225463, 0.000801742977301409, 0.0008601676408603608, 0.0007743077298483614, 0.0007179252557614198, 0.0012247871359885737, 0.0011884077028258256, 0.0010657137593217826, 0.0011790431950650325, 0.0010355181944595874, 0.0009425723244607048, 0.001374811741069551, 0.0012277144644810035, 0.0010875224910095833, 0.0011808438166945362, 0.0006012290159227857])
    cub_degree = 1
    tol = 5e-14
  elseif q <=32 # 2916 nodes
    SymCubatures.setparams!(cub, T[0.03292392932219886, 0.11096229399878556, 0.2187588190751789, 0.33793941097761626, 0.45858100716788414, 0.5668821592027113, 0.660117447230065, 0.7239844359090146, 0.02321871305376844, 0.07661618717846833, 0.1559791421368959, 0.2521738623055426, 0.3545692287061555, 0.45767816195407424, 0.549463769469837, 0.6306162242106939, 0.9880527787060993, 0.9603245926737669, 0.917796767609045, 0.8618396646416213, 0.7942524171593308, 0.717207518456062, 0.6331813264391405, 0.544874546742326, 0.022523593891196188, 1.881656949304226, 0.022589112180450763, 1.8016991536041418, 0.021937303878080378, 1.6961683654632131, 0.020717571282793754, 1.567239370370338, 0.019304270921206002, 1.4177804121419817, 0.01796687552792962, 1.2517374294009684, 0.016942911264399327, 1.07434433735354, 0.07513296810311446, 0.02242205911214415, 0.15488082837553172, 0.020838643046910903, 0.2504024584295928, 0.019124458464056495, 0.35343964686596324, 0.01825351502107181, 0.4562875795031461, 0.015742925168652153, 0.5510346131572275, 0.013835927667872514, 0.6290500221834444, 0.011855329175573038, 0.07264535599529108, 1.7029302051371862, 0.0703660911668078, 1.6063404921791082, 0.06734465861400131, 1.4895312304477362, 0.06370630775977852, 1.3531758681624746, 0.059428726072949835, 1.19922510042264, 0.055179130794721226, 1.032688819898488, 0.15043053333759399, 0.06989951706346584, 0.24479244308176534, 0.06382311578759144, 0.348052722849642, 0.06015713444753263, 0.4481502125965576, 0.05344863855046989, 0.5421402676582915, 0.046190479498127, 0.6199656774743899, 0.04135228477686586, 0.13935894422499082, 1.475406121458416, 0.1340941315659495, 1.3671454458749506, 0.12736049382919806, 1.2466081317388025, 0.11870649677079262, 1.1148719411204164, 0.10951256916103519, 0.9707327695514582, 0.2356654136047586, 0.13305834805525835, 0.3369375240580414, 0.1229096892047817, 0.43693038064266776, 0.11008370345560838, 0.5291750089934982, 0.09701694776174834, 0.6060176053553311, 0.0868926262238164, 0.21611348900889948, 1.2305400916704858, 0.2053362598373555, 1.124520615544315, 0.19156790179883398, 1.0124745811166238, 0.1773525641703557, 0.8926839343718739, 0.3231762055131826, 0.20455781874381024, 0.4210923801729993, 0.1830050501185761, 0.5109779260124654, 0.16308849304227413, 0.58736143154645, 0.14574261558482202, 0.29025982614802437, 0.9968505077844313, 0.2725053860078486, 0.9038608737142517, 0.2543354976310317, 0.8059491650836064, 0.40107831566409224, 0.2731370972742169, 0.489291470234981, 0.24204908722710436, 0.5637637874625003, 0.2165554413797762, 0.3575604954426512, 0.7921584570840499, 0.3348273279091471, 0.7163154386251188, 0.46533275476524766, 0.3369707950533015, 0.539346161066513, 0.2962912168626489, 0.41257200212428247, 0.6274111603872519, 0.5114288339252974, 0.3870938301572597, 1.8996496108572756, 0.023092166887234995, 1.81676414176096, 0.022715838377297376, 1.7074073590386896, 0.02190016787386545, 1.5748173470551243, 0.02064368665731146, 1.422688000649648, 0.019585055825585337, 1.2561048718941046, 0.01861624970654637, 1.0808477405119028, 0.0172820077297201, 1.766051495211983, 0.07538386482064902, 1.6609851555541308, 0.07267709679349796, 1.5346215483282108, 0.06852877906576661, 1.388707690087819, 0.06506031014508493, 1.2284430082860829, 0.061407475271634704, 1.0591012767277843, 0.05706415698149338, 1.5879766214823252, 0.15025750903992952, 1.4691418539032406, 0.14168359490021026, 1.3318603782356464, 0.13452633663066355, 1.1827648341624186, 0.12591884275049853, 1.0247834740282877, 0.117066670587145, 1.3868287877740464, 0.23790955375732145, 1.2603790725236477, 0.225354316669128, 1.1248510950616952, 0.209547690202714, 0.9807770173315405, 0.1944139815353357, 1.1785158683308932, 0.33441715076647405, 1.057898055579627, 0.3099719835521217, 0.9289685153380944, 0.2868439182948843, 0.9825174624912336, 0.42421743735884065, 0.870007736902752, 0.3930708327461422, 0.8057659998170382, 0.5102884711042318, 0.15647155512162672, 1.746747890284088, 0.07490303826648698, 0.26428824386837535, 1.6421843756312777, 0.07244737445756406, 0.39526410365890635, 1.516140814627538, 0.06844393104018116, 0.5450210882490384, 1.3720674199971898, 0.06378509694012593, 0.7080615354494714, 1.2146698747508287, 0.05943942301407873, 0.8784401739601426, 1.0491116662858635, 0.055745126017547976, 0.25963181661460627, 1.5706635680377845, 0.14951440810169653, 0.3856692877814089, 1.4533693506276686, 0.14165893986153164, 0.529121798245735, 1.319984064867113, 0.13226742426958046, 0.6861287515785289, 1.173445526108176, 0.12327739928181954, 0.8522742265472546, 1.017330111047808, 0.11450309359243557, 0.37160313362563613, 1.372185157096738, 0.23778902740525049, 0.5093053214768704, 1.2504834814964105, 0.22249567146605947, 0.659407917202795, 1.117053537323575, 0.20735349202158831, 0.8189091495239873, 0.9750015127200697, 0.1913262029688561, 0.4847863225360965, 1.1666852286358806, 0.33174628819895435, 0.6275553763128542, 1.0482850174159946, 0.30888739789474096, 0.7787525768313562, 0.9223610118546859, 0.2842789562773595, 0.5909844507615741, 0.9698085129790861, 0.42487483809084453, 0.7324335045069025, 0.8625339680951052, 0.39153539526298375, 0.6822243866525769, 0.7962038700053276, 0.5091892913731942, 0.25309169021962497, 1.534171606133728, 0.14570383302099207, 0.37737211722827124, 1.4189490343145277, 0.13928640558058775, 0.5190341300724176, 1.2887533532565736, 0.13068000040739872, 0.6732542381904079, 1.1476256571361745, 0.12209089126399521, 0.836370941029582, 0.997885323167505, 0.11307162845118932, 0.3638628255929324, 1.340371121671835, 0.23442674523942283, 0.4998383410546837, 1.2220657950938592, 0.2191350289453249, 0.6486596194925739, 1.093114450719928, 0.204103722935886, 0.8054287174802123, 0.956248342992197, 0.18838704172327422, 0.47725792536734896, 1.141519240808622, 0.32630220691313605, 0.6167393411599261, 1.028335756909739, 0.3034526191862394, 0.765423129346693, 0.9067905817297237, 0.27983113528678266, 0.5808850221624084, 0.9531344867819499, 0.4183188645610809, 0.7203451365172568, 0.8483410311963062, 0.38578511896497286, 0.6716984139032423, 0.7836048953820015, 0.5027220220680679, 0.35270005185894354, 1.2934139270245917, 0.22604973402225695, 0.48689762988405855, 1.1781445895278997, 0.2129796545695412, 0.6323447601341917, 1.0558237028298494, 0.1987012130670041, 0.7835319527986442, 0.9277698679437383, 0.18409487897549745, 0.46382227602280646, 1.1042931895779105, 0.317623281413211, 0.6018834153436915, 0.9953573768539736, 0.2959062232577004, 0.7463815355981862, 0.879624436988916, 0.2740168014458375, 0.5669712806093654, 0.9256948227543378, 0.40742546276256725, 0.702214345610701, 0.8259319375045897, 0.3770074001845192, 0.656260697863104, 0.7642501475874974, 0.4905079378419401, 0.44588509185840486, 1.0542039968153485, 0.30627112041978033, 0.5809081028246161, 0.9532277551593504, 0.28523077479210174, 0.7194907845342883, 0.8454356261756723, 0.2652476537561868, 0.5468255994363604, 0.8891673714636801, 0.39383090316409197, 0.678559004717918, 0.7956374976535279, 0.36562344140257175, 0.63438093942691, 0.7388095646495751, 0.47491470405507086, 0.5232016956768855, 0.8434189105232128, 0.3767343786099358, 0.6497813504567407, 0.7580749028954953, 0.3508721424260339, 0.6085612889048596, 0.7078261985337029, 0.4573266873274583, 0.5798430283028467, 0.6712977627572584, 0.43559813204580106])
    SymCubatures.setweights!(cub, T[2.780278763279196e-7, 4.980586832984121e-5, 0.0002770356963340783, 0.0006871926636965323, 0.0010552691414203702, 0.0012313982813138633, 0.001164717112636442, 0.0008854814750448663, 0.0005134889280729157, 9.074710485781593e-6, 2.851842926026058e-5, 4.972304244155218e-5, 6.599457922521704e-5, 7.642704263269801e-5, 6.648918720076546e-5, 5.7741032598145085e-5, 3.927030704522834e-5, 1.6012139856530567e-6, 2.8167880829180796e-6, 3.81921083663993e-6, 4.444303340819056e-6, 4.624564741201614e-6, 4.634913062223644e-6, 4.504063187306187e-6, 4.015911058805839e-6, 9.370173386343186e-5, 0.00013363182167225193, 0.0001604649545042098, 0.00017067336582065153, 0.0001675708543745398, 0.00015746206792308808, 0.00014492064078336752, 0.00016914322698732075, 0.0003015284366758477, 0.00039523981977983107, 0.0004568385229723104, 0.0004303064303563598, 0.0003435423609328787, 0.00023819478288418947, 0.0003741030892081357, 0.00043959651858296086, 0.0004783343347746024, 0.00048428117385094507, 0.0004651809928354501, 0.0004141521039076269, 0.000513593099065234, 0.0006778920510774924, 0.0007593343258683295, 0.0007308678322826667, 0.0005954740818762267, 0.00042202373587704424, 0.0007766934704263713, 0.0008212239771553452, 0.0008009724167172452, 0.0007683867045571581, 0.0006986243646255568, 0.0009117406590151537, 0.0009855700333553184, 0.0009395949684291361, 0.0008139808817032287, 0.000548015037429383, 0.0011164307779225943, 0.00108835018061284, 0.0009992693937360237, 0.0009087754590117037, 0.0011491281143650978, 0.0010667186530007513, 0.0008925041205639027, 0.0006421120845794251, 0.0012120727702282991, 0.001149503135604183, 0.0009791378354075913, 0.0011979363349543522, 0.0009623951745041691, 0.0006915176448279995, 0.0011179371002483858, 0.0009484696954818394, 0.00101522596098049, 0.0006976145131317656, 0.0007975604369522487, 0.0006572187691559951, 1.6393079800206103e-5, 2.278914601664292e-5, 2.6877264112372235e-5, 2.8062297930360863e-5, 2.7739460540036496e-5, 2.6271473147873257e-5, 2.3951149560166187e-5, 3.824490589360905e-5, 4.4427402001856526e-5, 4.707668050435677e-5, 4.754174470244523e-5, 4.451477231551199e-5, 4.0281224839431635e-5, 5.8255666955285214e-5, 6.0830213028607983e-5, 6.187964769003188e-5, 5.581375463728759e-5, 5.002661385042155e-5, 6.990877525393183e-5, 6.903281747974676e-5, 6.15001532580495e-5, 5.262603291482136e-5, 7.301953890220114e-5, 6.265841151518408e-5, 5.852960122718696e-5, 6.14060770093894e-5, 5.444706019268683e-5, 4.9396304961392855e-5, 0.00023113845025563067, 0.00027154937837191, 0.0002882941220165046, 0.00028328248878675406, 0.0002616075456471628, 0.0002331380950904389, 0.00034800143470983537, 0.0003691752242959021, 0.00036785850196283625, 0.000339403791693146, 0.00029728520557971167, 0.00041972171752273065, 0.0004190574670929112, 0.0003790058167212635, 0.00033143080226746455, 0.00043092004353100476, 0.0003933356312673984, 0.0003484521979248466, 0.00037107251951271753, 0.00033628793207202503, 0.00028987600268796693, 0.0005928122346252597, 0.000639853573014277, 0.0006210778297832203, 0.0005729141826303015, 0.0005045267213925907, 0.000724728588808197, 0.0007144348121369651, 0.0006491384602382493, 0.0005683002350843994, 0.0007186795652005503, 0.0006698318157569413, 0.0005937773192630831, 0.0006435359407799692, 0.0005873581895028261, 0.0005142375162676856, 0.0009760604006192424, 0.0009461367050732857, 0.0008645476691559972, 0.0007582699015596673, 0.0009866067559614498, 0.0008855031668914643, 0.0007897678380489191, 0.0008628852608726156, 0.0007645928091745051, 0.000700856059834038, 0.0011547677192141183, 0.0010531374480286728, 0.0009190843332236636, 0.0010256658767697987, 0.0008935801958782, 0.0008090701665339042, 0.0011136535576182153, 0.0009624224263285319, 0.0008405565151347767, 0.0008641274301413914])
    cub_degree = 1
    tol = 5e-14
  elseif q <=34 # 3439 nodes
    SymCubatures.setparams!(cub, T[0.029342246371486344, 0.994749260044068, 0.09950972820461466, 0.9828155897815947, 0.198712573002591, 0.9644315680237172, 0.3109510729112185, 0.9384741340499597, 0.428089234429479, 0.910188536214542, 0.5346471887476704, 0.873409919173587, 0.6190880058292592, 0.8387680811994246, 0.6912899963460113, 0.7936876971583103, 0.9852741031992712, 0.9541709185013736, 0.90598719030508, 0.8460410508110564, 0.7779551063857626, 0.7088896923221368, 0.6365384685935566, 0.5660307977643578, 0.02087996297483963, 0.0687657959397549, 0.1400942075337306, 0.22904399887821697, 0.3244246918508021, 0.4225035034208078, 0.5164179478612008, 0.5968962009656194, 0.9909754186855135, 0.9714893060655033, 0.9467140834714922, 0.9187810654688159, 0.8798110837026452, 0.8314878536112061, 0.7787816679235553, 0.726123093532526, 0.98930588311104, 0.9644507640762932, 0.926230288898323, 0.8757471012763065, 0.8144540686326103, 0.7441146428403568, 0.6667524239122493, 0.5845930117046407, 0.02026577413529052, 1.894043578088919, 0.020414610160103177, 1.8221714105647329, 0.01983731317718484, 1.7270663347402413, 0.018733626034514405, 1.6103079663050517, 0.017531377682412662, 1.4738213113011114, 0.01656641895798358, 1.321071083525595, 0.01593748487113416, 1.1560409539754544, 0.06769513750632393, 0.020062970685288402, 0.969384063146181, 0.013526905470925719, 0.14015522506002556, 0.01899141128858455, 0.9449756080307137, 0.01390913131860271, 0.2277171089736141, 0.01770004306584177, 0.9129673831858282, 0.013576346272357768, 0.3234134002019973, 0.016729478843354493, 0.8727164102382805, 0.012856929464227358, 0.42196777228754334, 0.014927266605974258, 0.825241222854672, 0.012199424904844125, 0.5156064012559421, 0.013040909869195915, 0.776080542312801, 0.011602368487350898, 0.5980171810431978, 0.011607932969039607, 0.7252587996694394, 0.011417099557440174, 0.0655556055465828, 1.7324449623200122, 0.06383751412812075, 1.6444961812772305, 0.06152124752162862, 1.5379469538565624, 0.05849978270700428, 1.4130119953002946, 0.054795674026639714, 1.270808143775504, 0.05085017715333722, 1.1156518075448254, 0.136011218077533, 0.06357541168125834, 0.9294416769415584, 0.04538723635311136, 0.22339741610526076, 0.058930555375619316, 0.8976028105353963, 0.044071420644933273, 0.3195389368657632, 0.054818323206478864, 0.859327576701754, 0.0426041393723541, 0.41592083048351536, 0.04958268972527626, 0.8138620342793665, 0.04042980406259035, 0.5078867437043131, 0.043820889686742186, 0.7639649569226232, 0.03877529071201577, 0.5883153971517883, 0.039206440778559455, 0.7156436348957722, 0.03686379766634618, 0.12801764140496866, 1.5206316893535177, 0.12302370925551248, 1.4210588840238254, 0.11659843899418615, 1.3077777561899193, 0.10989599116050877, 1.1823258095713922, 0.10166762272089001, 1.0488149577300765, 0.21629307903497216, 0.12265828901455393, 0.875235952828312, 0.09138328448992394, 0.30888677335426995, 0.11357711534086162, 0.8386612977915361, 0.08808773181920723, 0.4054472261613168, 0.10236911071513828, 0.7954299812192378, 0.08467273908968995, 0.495930201524675, 0.09141675951242778, 0.7468473906535745, 0.08142744397277041, 0.5756753453558913, 0.08078322716985253, 0.6983617489252516, 0.07795537963073143, 0.20044147272007154, 1.28949184917805, 0.19090517435757262, 1.1906688752409957, 0.17845464330159472, 1.084906582574412, 0.1669819536623483, 0.969136215411845, 0.2995235627398126, 0.1911818344271906, 0.8108717322584742, 0.15034932214445437, 0.39204692466882696, 0.17079090645531866, 0.7708128372122467, 0.14378803221571787, 0.47917746658989135, 0.1536615448446252, 0.7256256761152, 0.13769076769576516, 0.5568061921715343, 0.13494016681223284, 0.6768784710720931, 0.13042053334647621, 0.27155310686670414, 1.0630368402730268, 0.25671509262729725, 0.9737458980147804, 0.2395724429864425, 0.8797961825921544, 0.3758165881122988, 0.25656350019772545, 0.7402590137874027, 0.2162477750535999, 0.4616663503724882, 0.22691974981885404, 0.6980494030863261, 0.20420069452847006, 0.5354864874779006, 0.20253470298801923, 0.6527805161014754, 0.1941725242974572, 0.3379448802267046, 0.8600114073981319, 0.3167862774468422, 0.788588486413622, 0.43818115631422966, 0.31500804163086854, 0.667684991267596, 0.28394869732363986, 0.5095063025771248, 0.2751772636130997, 0.6243673089055545, 0.2660774183436894, 0.39284692807565186, 0.7018475116515689, 0.48405757530688426, 0.36333532290155873, 0.5955403855136036, 0.3493217861771322, 1.9097809206961698, 0.020701947099993952, 1.8350983922783348, 0.020317844321620204, 1.7360431634277984, 0.019794122468429938, 1.6155468191667306, 0.018796019451068163, 1.4767619143526631, 0.01759796701910264, 1.3233048593378425, 0.01693455896928361, 1.1611521286328597, 0.014501061604756487, 1.7899081894668152, 0.06752946481275546, 1.6942579383447123, 0.06579642432776667, 1.5792257729511674, 0.062356861707573014, 1.4461675962696445, 0.058864788836341886, 1.2981244140858863, 0.055880025454212334, 1.1411760531336008, 0.050078261953227406, 1.6277875980656276, 0.13631496767508472, 1.5187331981980707, 0.12891768995190417, 1.3918465629393817, 0.12246242971781919, 1.2525337540363704, 0.11516622190496596, 1.105052441785466, 0.10633720288126963, 1.441741186638359, 0.21687345876510444, 1.3245514845186277, 0.20616666989104557, 1.1961373214124558, 0.1940723697547783, 1.059271392227268, 0.18054317027492742, 1.2453601669765948, 0.30691184233155966, 1.1274913457938105, 0.2904222097435548, 1.0053609064657212, 0.2690501205016952, 1.0512083267736851, 0.39935883483120954, 0.9434781388674804, 0.3693737003918547, 0.8785179611793953, 0.4785347431708478, 0.1411584590523768, 1.7713992826063358, 0.06772110565289186, 0.23882386417647503, 1.6765351257856043, 0.06558288423330239, 0.3580136941387545, 1.5617633176430668, 0.0618196892758698, 0.49512317584997945, 1.428897752502194, 0.058394454984533085, 0.6459763610037661, 1.2829881119672553, 0.05467039229770881, 0.8051428224603242, 1.1274196112233699, 0.05222971932737169, 0.23536722526415865, 1.610472116339389, 0.1356179339909656, 0.35059954161931844, 1.5029775627659951, 0.12814371195079613, 0.4826794329537563, 1.379465243621183, 0.12177465168655474, 0.6277817425278849, 1.2420631480412352, 0.1136451187314416, 0.7830003746033334, 1.096335692591769, 0.1066018335041359, 0.3391777034517775, 1.4279551771616952, 0.21614471228499207, 0.4661636720761086, 1.3125302001792531, 0.20512966113657083, 0.6060327456332066, 1.1865705349008244, 0.19223079343528474, 0.7557564902353316, 1.0522150614196788, 0.17823891968102665, 0.446168322312181, 1.232161590712875, 0.3059265600927851, 0.5800250159349687, 1.1184711791118063, 0.2872284752904891, 0.7217166514248049, 0.9975298745948978, 0.26685593877111624, 0.5493625337713138, 1.0412033976266326, 0.3956080163048625, 0.6836219336315624, 0.9350304732549491, 0.368695979960234, 0.6423763966555887, 0.866695355909381, 0.47888231905654177, 0.22959582685615218, 1.5761603345785296, 0.13264121985841612, 0.3437517708745209, 1.4692164114584554, 0.12713685937121502, 0.47490463602961536, 1.349877777136444, 0.12014208477791588, 0.6171897991060361, 1.2168504521318895, 0.11253618181798127, 0.7697408617089171, 1.0767631252541108, 0.10519877637883768, 0.33281941020094746, 1.396402708927549, 0.2142929198149679, 0.4584304003688542, 1.2858235270001206, 0.2019862693004477, 0.5966384615558337, 1.1635902854396638, 0.18923466104822712, 0.7435347195580245, 1.0334221955135006, 0.17601436373038296, 0.4395610820170381, 1.2074966854561247, 0.301304170436784, 0.5707094020035489, 1.0981610543304936, 0.28306041881462907, 0.7109191464032858, 0.9808776539609045, 0.2628048526317255, 0.5414585637688932, 1.0228969000269363, 0.39002260965024804, 0.673071767015262, 0.9200119615711366, 0.3643776620851791, 0.6329084961316888, 0.8544875921890124, 0.4720507585815483, 0.3238959886173777, 1.3513004806494233, 0.20624737753461092, 0.4470180261095249, 1.2434456052271086, 0.197041762897604, 0.5829034232594739, 1.1272318365063607, 0.18515052930968914, 0.7265031078269136, 1.0027659205977086, 0.17174789230888202, 0.42812157635598436, 1.1716413748276369, 0.2943686606505854, 0.5577634850291239, 1.0647440580173275, 0.27676393919965037, 0.6942611250471539, 0.9536367259276864, 0.2579518194134113, 0.528952515877634, 0.995111109532388, 0.3813525599246105, 0.6589317543221005, 0.8964474427300055, 0.3564510534545794, 0.6192427762947775, 0.8349765440533423, 0.4620340874309504, 0.41385959071435585, 1.12238799974317, 0.2844227357344164, 0.5390405766181463, 1.0225810320043578, 0.26812267878476287, 0.6721938862290309, 0.9181260894623057, 0.250147056495663, 0.5123653614013709, 0.9582539490628205, 0.36945099869034265, 0.6384635044316825, 0.8654861594161238, 0.34640891565714205, 0.6011706863612518, 0.8101339118959419, 0.44850071867553454, 0.49111479883561315, 0.9127895616597855, 0.35526101976157737, 0.6138676946553207, 0.8265306350057091, 0.3327860738244814, 0.5792009917476904, 0.7788469323288272, 0.43225131771694003, 0.5520645123768301, 0.7401062110759693, 0.41242427119762726])
    SymCubatures.setweights!(cub, T[1.9912998467485588e-7, 3.5295434491601796e-5, 0.00021711247025701523, 0.00020159200809844602, 0.0003732581861209702, 0.0005152174187986969, 0.000519171426446207, 0.0008138516310988468, 0.0005527173585893715, 0.0009832867187760568, 0.0006083003372038188, 0.000988498262455316, 0.0005936962944825147, 0.000650435042190643, 0.0006266888259244219, 0.0006109633167344667, 0.00065689592510655, 4.225078919333619e-6, 0.00010347532329599887, 0.00028950512444034005, 0.0004990900551615237, 0.0006531111051047228, 0.0006801692310193614, 0.0006704671343953672, 0.000689844624669115, 0.0006825242923744503, 6.550438205965679e-6, 2.060287043936785e-5, 3.7018731323065966e-5, 5.091751599697767e-5, 6.101445347011549e-5, 5.853203914757648e-5, 4.774516545591513e-5, 4.093040680862794e-5, 1.9984401781725996e-5, 2.846341871524523e-5, 3.165221410546044e-5, 3.7315504353881796e-5, 4.506624319803882e-5, 4.5656839270894655e-5, 4.063073416478815e-5, 3.9177186267196336e-5, 1.1563937792622417e-6, 2.0300396127487383e-6, 2.7469810719446626e-6, 3.2825405522231424e-6, 3.4932670342531e-6, 3.405459757915415e-6, 3.4699659760358672e-6, 2.428263572273966e-6, 6.792361065771395e-5, 9.810486077088883e-5, 0.00011870840303280345, 0.000126659194152073, 0.00012765295629563677, 0.0001229545878966626, 0.00011963825594257303, 0.00012331069624177033, 0.0001515022522287307, 0.000225899554338137, 0.000204143473590089, 0.0003065185116441621, 0.00024472822576695135, 0.00036334848533108347, 0.0002725492506371013, 0.0003570905458452723, 0.0002606799328161308, 0.0003120772369368194, 0.00022154475617933004, 0.00024912728736002247, 0.0002093376341813553, 0.00027591717797411833, 0.0003304016096709192, 0.0003631116001069157, 0.0003748590154528255, 0.0003677912361040409, 0.00033956577320780915, 0.00038524645623597083, 0.0003622116070950602, 0.0005266997351692214, 0.0004093647590027147, 0.0006034303985082181, 0.00044132988860249965, 0.0006117203195950386, 0.00046316825903675864, 0.0005354689127145436, 0.00039541342367147014, 0.0004353781691150377, 0.0003602485149279577, 0.0005998120659057249, 0.0006246119343225412, 0.0006675981111317621, 0.0006149506552385072, 0.0005699613681744138, 0.0007088742168947227, 0.0005495531212736911, 0.0008255911623293346, 0.0005907380900934331, 0.0007900233120222531, 0.0006107172803410964, 0.0007066547509241312, 0.0005598769562633435, 0.0005527764285708994, 0.0005009234148396633, 0.0008886353629422977, 0.0008906331815075586, 0.0008693709844405632, 0.0007613417475799683, 0.0009453563769736016, 0.0006937823009472512, 0.0009271424679177871, 0.0007142529443382871, 0.0007894618663031648, 0.0006482904380484864, 0.0006516900968992556, 0.000609513655107206, 0.0010405654082792412, 0.0009836711294783996, 0.0008857812248279053, 0.001030084411652146, 0.0007276980102894953, 0.0008278320555507074, 0.0007232536717939517, 0.0006979351480432321, 0.0006895516673924907, 0.0009607869145559685, 0.0008784680809038548, 0.0008252080541018061, 0.0007620831752744662, 0.0006700891432810024, 0.0007468366206174332, 0.0007632300737582749, 0.0007311668451679651, 0.0007819974607049982, 1.1915365386058294e-5, 1.663819867657342e-5, 1.987743251611421e-5, 2.1022067946318046e-5, 2.070172420163359e-5, 2.0408394890314477e-5, 1.8157803452401706e-5, 2.7929142365181533e-5, 3.2950557980765116e-5, 3.5559268610743426e-5, 3.669421172150535e-5, 3.451817393102084e-5, 3.1811564293269295e-5, 4.426951782296347e-5, 4.877069868593164e-5, 4.4452852008238735e-5, 4.773768525006437e-5, 3.8573678914802054e-5, 5.27979163019708e-5, 5.455184633099069e-5, 5.1620323553457556e-5, 4.5189953658639475e-5, 5.9049791041127184e-5, 5.389684139042792e-5, 5.081319312110603e-5, 5.3841131793634344e-5, 4.6244530863796525e-5, 4.324890813915638e-5, 3.4721897608414376e-5, 0.00017055969929780397, 0.00020257957649755625, 0.00021747632848968077, 0.00022109834927587377, 0.0002063534845517534, 0.0001862513407297873, 0.00026569642577323877, 0.0002895028815085534, 0.00027844185746087237, 0.00027816712268756937, 0.00023519448557994912, 0.00032750339317772214, 0.00032767833015840347, 0.00031439310141310556, 0.0002774388518921736, 0.0003566297956530771, 0.000323804248012919, 0.0003029940837777027, 0.00032559435465889417, 0.00028724231745351015, 0.00026494672988505475, 0.00045769852128528776, 0.0004925173751769604, 0.000493320648568659, 0.0004528935717418126, 0.0004221208873069081, 0.0005740534697965577, 0.0005669885826468717, 0.0005320935532872687, 0.00048616560101890635, 0.0005877262272852191, 0.0005696420692325319, 0.0005115631960263997, 0.0005461480247238033, 0.0004963088783019366, 0.0004466046638990977, 0.0007640856703409359, 0.0007745882647213911, 0.0007198613524571971, 0.0006536917028468471, 0.0007971619057332943, 0.0007590410308397126, 0.0006889975704782214, 0.0007398632098348963, 0.0006522272309020055, 0.0005839392534013979, 0.0009518084229540507, 0.0009032537488887851, 0.0008126659800157187, 0.0008653020158596169, 0.0007820172377430854, 0.0006931164395918919, 0.0009730001646247824, 0.0008811641948465469, 0.000746968500107917, 0.0008253769347446382, 0.0004912272140030157])
    cub_degree = 1
    tol = 5e-14
  elseif q <=36 # 4000 nodes
    SymCubatures.setparams!(cub, T[0.026754465690901147, 0.0908041479524711, 0.18079093481364383, 0.2847098870512747, 0.3969510899198101, 0.5017489096718918, 0.5921798080058434, 0.6732806594101794, 0.7265348867224719, 0.01894176483035853, 0.061816073159419206, 0.12643600829156484, 0.20850063414635683, 0.30321120631787235, 0.39991144153174507, 0.49159217052238263, 0.5809793286938126, 0.6447965478902975, 0.9903718524469571, 0.9679672494063327, 0.933438989044975, 0.8876841304760279, 0.8318882011451556, 0.7674964320159432, 0.6961765918569547, 0.6197758529614933, 0.5402729686194109, 0.018197545533218864, 1.9040429994683161, 0.01835179155187852, 1.8388023604468993, 0.018037138564491417, 1.7522442092830042, 0.01724856500733124, 1.645915568029415, 0.016238743150110996, 1.5215920444938622, 0.015251686512809282, 1.381588597247137, 0.014305929299716852, 1.2292757395059435, 0.0134776781758311, 1.0686938494451916, 0.060792769208544095, 0.018283390546154285, 0.12682133380238433, 0.017276629050546975, 0.20791715006863426, 0.01595197855329089, 0.2977718154665299, 0.015279177479443849, 0.3912778783262341, 0.013664449705002822, 0.48266197947519623, 0.011861728113778105, 0.5661873341694587, 0.010558650126566312, 0.6340024295988647, 0.009616427139578733, 0.05978462916986172, 1.7556026296801002, 0.05815987941039526, 1.6747130710418998, 0.055900184696949344, 1.5768487776205715, 0.053538486449972916, 1.461984071780823, 0.050986675797660176, 1.3309512447604428, 0.047850359328596534, 1.187098805288655, 0.04440518269051769, 1.0349628382575464, 0.12405620358572268, 0.0577729980604054, 0.203773138656222, 0.05333981973463903, 0.2936634085965918, 0.05075822812377151, 0.38605884223423326, 0.04600686584348682, 0.4752992105022292, 0.041248721489810174, 0.5573726358399269, 0.036839728470106047, 0.6250438305129388, 0.033125536036669415, 0.11580680380203057, 1.5647076026937188, 0.11288321878382031, 1.4707446521325624, 0.10864361436231305, 1.3652855713755516, 0.10315170439235931, 1.2490482032096804, 0.09587754945908375, 1.1225687960166923, 0.08928637676845337, 0.9853482904262426, 0.19659762538572612, 0.11167316359569848, 0.2860835207102877, 0.10514701028314032, 0.37740231951275244, 0.09563458101432724, 0.4657151082301722, 0.08680577267412738, 0.5469241892111689, 0.07705681072584718, 0.614398286494469, 0.06988862687326358, 0.18445257438828028, 1.3457225463067792, 0.17685046082213823, 1.2503145316373903, 0.16785970451931376, 1.1458141170419962, 0.15559521784375585, 1.0379722677871304, 0.14628674366928454, 0.9198406830432562, 0.27617818269554234, 0.17693448140221854, 0.36543752601376384, 0.16107998530439727, 0.45186554208231944, 0.1462171889180436, 0.5322660128363557, 0.1313577603553357, 0.599518779496862, 0.11758630177466994, 0.2530400018353012, 1.127986629326051, 0.2415844321374558, 1.0385548272983731, 0.22653561081236784, 0.9456826590480261, 0.2118826533507639, 0.8470748331093213, 0.35167357407865, 0.24085086225441818, 0.43629176153085186, 0.21702311967394086, 0.5145235479101767, 0.19570317891220704, 0.5820890747834244, 0.1763106943762475, 0.31778374076389554, 0.9252532272074785, 0.30011982734457493, 0.8498053903274939, 0.28239226938603296, 0.7696931826798016, 0.41719013829762547, 0.3006364696643953, 0.4947578443415419, 0.2689549933872242, 0.5592344038311882, 0.24451803697426436, 0.3758313558969692, 0.7513513765300593, 0.353116119545926, 0.6897413821684255, 0.4708421395277051, 0.3557916668149538, 0.5365032753619372, 0.31668222253909706, 0.4228639217544909, 0.6122578697908995, 0.5108600714691002, 0.3995118656282075, 1.9183549033456793, 0.018571681116591108, 1.8504983815398175, 0.018308977800245974, 1.7604201504525863, 0.017942212872515937, 1.65040574936511, 0.017427672026863218, 1.5233685705884894, 0.016588533230919168, 1.3821944273495466, 0.015932908007251292, 1.2307658587225636, 0.015168906392015537, 1.0734132411859731, 0.013578807143490055, 1.81057839628979, 0.060899702444417045, 1.7232545175922793, 0.059663272186598844, 1.616578872729244, 0.05803280037805141, 1.4937559918624357, 0.05525690315823585, 1.3566074432759891, 0.052838950513637024, 1.2090095362007585, 0.05043856759487502, 1.0561123392325074, 0.045723136769648634, 1.6631059751082655, 0.12380169064276847, 1.5608943750162068, 0.12056324959425499, 1.444653013184646, 0.1149084436767297, 1.3148827077957512, 0.10925376682930046, 1.174047537156922, 0.10452589518367211, 1.0283930220217068, 0.09626773427103735, 1.485580142216857, 0.20313929586660695, 1.3766726349187943, 0.1938618040283219, 1.2569224171974605, 0.18323485028526507, 1.1259205826542678, 0.17533447001195032, 0.9903943984236928, 0.16416980204736373, 1.2922312602203085, 0.2895955678281579, 1.183090376406563, 0.27254445495898655, 1.0640427143587934, 0.2598230255879276, 0.9421013632971927, 0.24690690455392922, 1.102522816613522, 0.3761951569474201, 0.9968975134323241, 0.3563094643728103, 0.8881042879362973, 0.3402668446718541, 0.9239301286091602, 0.46542992626003427, 0.833587851232978, 0.4386134304285687, 0.7691290963988814, 0.5396788226945799, 0.12699090670463933, 1.7939522273082669, 0.061063137751673595, 0.21546874794141258, 1.7073514849255569, 0.059745435760419525, 0.3241886847033874, 1.6020298433681988, 0.0570956769095629, 0.4502901581874317, 1.4798677138315834, 0.05378549870955195, 0.59003067032535, 1.3440872129578405, 0.05058240762710454, 0.7391560515681299, 1.1991129529176638, 0.04740259385889485, 0.893591228650802, 1.0483248357320687, 0.044672699018182334, 0.21375760085246742, 1.6459878302231359, 0.12365049036439854, 0.31941093493283534, 1.545963002602668, 0.11835743542015603, 0.4410424534198637, 1.431378758832877, 0.111784597259868, 0.575715871405377, 1.3039843192010412, 0.10533274984964597, 0.7207266015806241, 1.1668504349529363, 0.09873584690636009, 0.8720946094633313, 1.0224804367606701, 0.09241688918786004, 0.31022463644955595, 1.4746559953012872, 0.19934282902853334, 0.4279970830516797, 1.3679035460729216, 0.18868394313083595, 0.5581980673082674, 1.25004862969609, 0.17806966648698194, 0.6974092481039887, 1.12244583635643, 0.1670143246787732, 0.8441284180908682, 0.9884789912396992, 0.15556613355048252, 0.41092960559385083, 1.2918167880278661, 0.28268918667498016, 0.5361035556383577, 1.1839421501526872, 0.26665731946150883, 0.6700998857996368, 1.0675686594933425, 0.2501600855943411, 0.8105800754812915, 0.9457214071305635, 0.23191655970471273, 0.5108314283616625, 1.1073197213285317, 0.3688923548334204, 0.6383409420937394, 1.0043827213851473, 0.3452496556262815, 0.7716444429644497, 0.8962595991959739, 0.32052961573593125, 0.6027958932400064, 0.9342930005765624, 0.45128440674649006, 0.7281680084326677, 0.8415843722317907, 0.41975916902394594, 0.6826963123415778, 0.7816180656909841, 0.5251598226163093, 0.20930848464841853, 1.6144785176704215, 0.1207910839821355, 0.31327608286892417, 1.516368354944333, 0.11617474160210392, 0.43350661814718466, 1.4035520435584488, 0.11057124738188184, 0.5665028613870032, 1.2788915794231834, 0.10491325156096634, 0.7092133573685045, 1.1467154481864668, 0.09847647278876571, 0.8584559486952343, 1.0068380962899262, 0.0916163175215255, 0.304494138550265, 1.4466679479862243, 0.19633762067818894, 0.42056108633451916, 1.3417876575589192, 0.18653817859312619, 0.5498655719110555, 1.2275589831564844, 0.17616341926858906, 0.6883078863357005, 1.1031037970406115, 0.16508722511023444, 0.8328046466265313, 0.9730301026468345, 0.15349629728234668, 0.4057003877815498, 1.267043110800546, 0.2793759589260874, 0.5289628230678202, 1.1631776461692092, 0.2628094024106854, 0.6607677217582877, 1.0512220110467962, 0.24642898483067407, 0.7994451435935094, 0.9322594799245747, 0.22912045352126256, 0.5040013973360322, 1.0894890808021043, 0.36351909698833623, 0.6291363211873016, 0.9897259444058982, 0.3410458159836395, 0.7607519874503399, 0.884237741744259, 0.3168363226879933, 0.5953558381493838, 0.9206404274517177, 0.44641840495622176, 0.718943988760101, 0.8303314374362466, 0.41454362309231174, 0.673917972454777, 0.7718552594269688, 0.519968872300118, 0.29558280353885963, 1.4043339499963896, 0.19087409969344274, 0.4109410313215303, 1.3013431659307182, 0.1822891602670474, 0.5389209189462285, 1.1904318701624965, 0.17233427783839453, 0.6746791001783963, 1.0731777800786568, 0.16161102987402493, 0.8140884642379208, 0.9490452108778555, 0.15056770970981706, 0.3960603677410883, 1.2319834325974517, 0.2727206484953099, 0.5174240478025348, 1.1307834689499012, 0.25718343606038835, 0.647297891314082, 1.0227413412079889, 0.24165597619481227, 0.7826619225720058, 0.9099860680623019, 0.22519162484229927, 0.49345119597484904, 1.0609176233719828, 0.3557246067173962, 0.6165421265480938, 0.9655087722848994, 0.33435456540157854, 0.7451714445236081, 0.8645520044628282, 0.3109887873730867, 0.583669437399411, 0.8999092871784652, 0.43732792345143817, 0.7048030117205609, 0.8124440612824935, 0.4069476785039418, 0.6606880298910431, 0.7572884590857182, 0.5104576664334265, 0.38290957658495184, 1.1850193749121916, 0.2640970440880439, 0.5011310828781481, 1.0873557292968745, 0.2507531293123569, 0.6286297616607178, 0.9862482336669807, 0.23512054579249728, 0.7596189042213037, 0.8804066813481598, 0.2192544513543838, 0.4787614635029351, 1.0232772866390847, 0.34579069952999386, 0.5995990262507038, 0.9320233007527382, 0.3253022456557349, 0.7244704502061456, 0.8374133703032717, 0.3036863515974658, 0.5671617699807964, 0.8729367244394474, 0.42480308730631783, 0.6852568118154008, 0.7897101502639652, 0.3965816970790194, 0.6440655198244278, 0.73651475068209, 0.4973145045944067, 0.46109826773860685, 0.9769199175322039, 0.33331321525704005, 0.5779782706621683, 0.8926619927641714, 0.31390358387110834, 0.6987641515364923, 0.8043789586879218, 0.29350650574702486, 0.547286743287771, 0.8385287391994461, 0.4105626714453101, 0.6611545163633816, 0.7607831509461249, 0.384260637587503, 0.62346689253479, 0.7112704363499482, 0.4817814145490059, 0.524619134777727, 0.7977389941975636, 0.39383921281487555, 0.6337056226962524, 0.7265203753758646, 0.36971029001763556, 0.5992617932955134, 0.6833417015108338, 0.46501129073864605, 0.5720945339261151, 0.6500702040185572, 0.4442815923782751])
    SymCubatures.setweights!(cub, T[1.4562826837329852e-7, 2.6738528640848785e-5, 0.00015429406362154518, 0.0003909485890768612, 0.0006566593384827444, 0.0008126228921371345, 0.000900957355741184, 0.0007947971987179538, 0.0005950972067549413, 0.0003514316886399991, 4.905712957912532e-6, 1.5249461620547559e-5, 2.7546255296119016e-5, 3.852220483178303e-5, 4.8112212518180824e-5, 4.793086414572952e-5, 4.0752253975146393e-5, 3.247318719404759e-5, 1.6372686959953974e-5, 8.59142701063408e-7, 1.4710991406843036e-6, 2.019343027589663e-6, 2.451735866430391e-6, 2.734924997189703e-6, 2.794286885638834e-6, 2.8135990723816778e-6, 2.6660552753937685e-6, 2.2005220734834503e-6, 4.976794080378921e-5, 7.201945223845209e-5, 8.897423652499935e-5, 9.775636638131961e-5, 9.93130558505453e-5, 9.681130280193087e-5, 9.061765913216366e-5, 8.33844246680532e-5, 9.090498276343645e-5, 0.00017031918504320916, 0.00023334150412605067, 0.00028344216449170033, 0.00029299623800134743, 0.0002654239603554628, 0.00022132358652232517, 0.00015129826796479742, 0.00021118021705143506, 0.00025131709330691203, 0.00027753062041246945, 0.0002893239596692865, 0.0002920939109759685, 0.00027709295762541977, 0.0002457899524306705, 0.0002907653027761461, 0.00040489930115670934, 0.0004890176427530207, 0.0005019549750343047, 0.00047591655811552694, 0.000379046475835544, 0.000279260264540438, 0.000456660662207228, 0.000499718588726745, 0.0005075643701499922, 0.0005003462096721712, 0.00048618771444540395, 0.00043255893554470994, 0.0005511427466754103, 0.0006538927320208735, 0.0006750997618663875, 0.0006348479607304005, 0.0005099367410950926, 0.0003663230381924519, 0.0007117096090098148, 0.0007141138396056875, 0.0007164768080134648, 0.0006792532696394775, 0.0005809959985031408, 0.0007635121920721387, 0.0007964266512588645, 0.000728417770678046, 0.0006325802003164003, 0.0004267433930598031, 0.0008636392405112176, 0.0008188830147312088, 0.0007837191494652136, 0.0006808282043050265, 0.0008799199336842759, 0.0007851424248047646, 0.0006606058587938165, 0.0004719505302375921, 0.0008503940071710778, 0.0007869192816075546, 0.0007001028780277797, 0.0008400273464174014, 0.0006865753037996108, 0.0004787983424778078, 0.0007898740924562781, 0.000624997576395438, 0.0006911111215422757, 0.000471800425144867, 0.000541586132656699, 0.00045619234654971087, 8.712516300697183e-6, 1.2228345569722455e-5, 1.4907526410049551e-5, 1.6376464290421583e-5, 1.659035755916955e-5, 1.628438347408185e-5, 1.548031710548281e-5, 1.3444383117429498e-5, 2.0834147170596318e-5, 2.5015698041864088e-5, 2.7624247228546468e-5, 2.878279882864869e-5, 2.8506919477540664e-5, 2.6853554410483968e-5, 2.4086789989275378e-5, 3.269704228321e-5, 3.684540991378742e-5, 3.865748357653321e-5, 3.740471731518909e-5, 3.48263202334025e-5, 3.2081043365254265e-5, 4.365053345842367e-5, 4.5678324971668314e-5, 4.038289948331195e-5, 3.908301648102746e-5, 3.457923062887195e-5, 4.910013958051332e-5, 4.368476258078592e-5, 3.8377742292279335e-5, 3.78892567816412e-5, 4.6551019877616556e-5, 4.012033841195317e-5, 3.5465828993995446e-5, 4.1367592675039314e-5, 2.8336713432039006e-5, 2.6806992572489475e-5, 0.00012725471576460566, 0.00015334387563738247, 0.0001673482219333436, 0.00017145226281659712, 0.0001666593651014413, 0.0001535986740070922, 0.0001369818728700118, 0.00019962011544480532, 0.0002207361015885259, 0.00022771426707771957, 0.00022196233132768206, 0.00020101762802674048, 0.0001818191736630698, 0.0002600389667229122, 0.00027059125075615487, 0.0002474840973733381, 0.0002324633724013991, 0.0002036048840789814, 0.00028676445894549067, 0.0002711365463277524, 0.00024257781830488864, 0.00021937633228431264, 0.0002757613969390301, 0.00024838808909818504, 0.00022639382537474305, 0.00023521682753563445, 0.00021127928279365577, 0.00018661446870221393, 0.0003432987848981795, 0.0003821441057489099, 0.00039140152393363984, 0.00037757226568276285, 0.000344908740689134, 0.0003094713866170045, 0.00044959983751956527, 0.00046318365533727803, 0.00043904978154475855, 0.0003924629302592646, 0.00036158319013788673, 0.00048684139446535633, 0.0004707715441295045, 0.00043419143232035834, 0.0003816906535852336, 0.00047412362742419193, 0.00042862891133240756, 0.00038840754547577323, 0.00040453923123304225, 0.00037484629832599376, 0.0003213917972978349, 0.0006113481987918792, 0.0006186477076652098, 0.0006038211186992972, 0.0005391456464767247, 0.0004892159922902839, 0.0006590957943376966, 0.0006414046268412501, 0.0005900707374328927, 0.0005201546415542233, 0.0006458717181127154, 0.0005859350646128007, 0.0005215204848964281, 0.0005652910051480652, 0.0004967673982045334, 0.0004491364542923841, 0.0007873203458003095, 0.0007646689225320042, 0.0007068292252825738, 0.000623896094102932, 0.0007606673669120502, 0.0007010705603346936, 0.0006207255887850219, 0.0006692847269896878, 0.0005929963938000916, 0.0005394720103930157, 0.0008361649203891661, 0.0007609730906429257, 0.0006815226424945602, 0.0007280193082411766, 0.0006566117065221034, 0.0005755083974424068, 0.0007780454481719327, 0.0006678742324827166, 0.0005831672085895834, 0.0006199038864598796])
    cub_degree = 1
    tol = 7e-15
  elseif q <=38 # 4641 nodes
    SymCubatures.setparams!(cub, T[0.023990094862258646, 0.9959218306481689, 0.08231260142667242, 0.9858778173688648, 0.16549844681690498, 0.9698990679797042, 0.2616111238127977, 0.9502597878876895, 0.37043422040278545, 0.9242695708167444, 0.4724467432891474, 0.8959728594368369, 0.5630901176595424, 0.862876647901953, 0.6369581756276775, 0.8304789993398476, 0.6989071366102237, 0.7900935604427468, 0.9883671116433795, 0.9635294190099523, 0.9235129974471623, 0.8733644974745259, 0.8133847258531196, 0.7498231260788922, 0.6866034739306096, 0.6218034848949681, 0.5585204812492363, 0.017130461689212957, 0.056255950615688645, 0.11525556122196583, 0.19049292290426412, 0.27500927170437267, 0.3653216631347806, 0.4503942845194703, 0.5321603379829839, 0.6086339024072717, 0.9922198977344079, 0.9749681231366449, 0.9530078521863322, 0.9322140444198274, 0.9012549058521676, 0.8601544588818287, 0.8151184876160785, 0.7693430258227226, 0.7242809798928879, 0.991286148302274, 0.9709881484798728, 0.9396473776617953, 0.8980009630388561, 0.8470255130311116, 0.7879159801309153, 0.7220578916395011, 0.6509949282543824, 0.5763927579010928, 0.01652184788851134, 1.9134897361074064, 0.016720724787375112, 1.8545627676663377, 0.01644225528149931, 1.7761088530887665, 0.015737672453974852, 1.6791420447493648, 0.014854767341251199, 1.564981496178023, 0.014034309454596827, 1.4355078778969645, 0.013433438052078342, 1.293595720805757, 0.01315529899571407, 1.142734985473923, 0.05531643850985449, 0.016576104136548905, 0.9753203598480703, 0.010324943620815935, 0.11557512401727912, 0.015748718624983515, 0.9543964728256177, 0.011815458825957723, 0.19021975188689333, 0.014855140591827808, 0.9276169297001904, 0.010798186643139822, 0.27407965447091853, 0.014411076606964288, 0.8948735406164595, 0.010361763975426232, 0.3630062516608466, 0.01309455281816576, 0.8561291291195542, 0.009881036112957655, 0.45093672090742704, 0.011770218952038043, 0.8120770627548138, 0.009808688058020149, 0.5334608883902732, 0.009727721088504145, 0.7658195838087276, 0.009243962234649994, 0.6055842957215843, 0.009536905755858102, 0.7191839730540613, 0.008795733768299011, 0.05453811940079126, 1.7777702841511105, 0.05327193340351346, 1.7039165904316687, 0.0515040199482008, 1.6142169207338257, 0.049453330392802274, 1.5085043355226644, 0.04689887947474692, 1.3875682187358407, 0.04387162730725795, 1.2537656168669291, 0.04134594936525573, 1.1098273080155567, 0.11330552525257881, 0.05275066875552214, 0.9417672394459713, 0.037528574897610655, 0.18747751701453538, 0.049442963590494445, 0.9148928465566767, 0.036786786190526805, 0.2706224530687066, 0.047332494787658566, 0.8831917695720277, 0.03530052048426714, 0.35770537367008853, 0.043450940151668116, 0.8452572243035225, 0.03443603216297384, 0.44498780738124033, 0.038917922179089494, 0.8028477049960361, 0.03262795153661771, 0.526884346101555, 0.03471850784897761, 0.7574915784785169, 0.03165495453054062, 0.5979422555078616, 0.03131552457172249, 0.7093007926536641, 0.03064815079096278, 0.10685628760142381, 1.5992990703731453, 0.10434444082100984, 1.5115246931130326, 0.09997523078743718, 1.4135224677683496, 0.09550387262800626, 1.3038234646414943, 0.09014588591743326, 1.1846637010477405, 0.08330559084461076, 1.0574685541053552, 0.18130263188815468, 0.10298640842774381, 0.8964073367505376, 0.07585149708256561, 0.2640935331032218, 0.09741413052452665, 0.8645383748213028, 0.0739592453633065, 0.35078091707329184, 0.08944683537035358, 0.8291616549238084, 0.07146397238721319, 0.43646657302853187, 0.0807524627532933, 0.7880001369842718, 0.06940170868373087, 0.516780405079056, 0.073339315583006, 0.7437254347729815, 0.06612378429757954, 0.5849142095319416, 0.06601179318252726, 0.6972643084712135, 0.06278270709137046, 0.17052931998149767, 1.3970748205019319, 0.1650440331773624, 1.305441595923735, 0.1556230853245676, 1.2067538888317273, 0.14674946151884555, 1.0980688126494729, 0.13762867784364868, 0.9857056944557046, 0.25649338688003814, 0.16433258043227078, 0.8413114930622198, 0.12446990217665074, 0.34064367878845386, 0.1506178685529867, 0.8063757759532179, 0.12078007921737055, 0.4248365580901501, 0.13752495292466904, 0.7686239753789254, 0.11708606837338287, 0.5030495056467492, 0.12450182890463556, 0.7248442705684585, 0.11212549302081504, 0.5723146198024157, 0.11054959122406811, 0.6804768148944281, 0.10672375895227193, 0.23650087050109792, 1.1850552194580442, 0.2266896941769336, 1.1007141736521502, 0.21388536592055638, 1.0111790292723546, 0.2011254068093387, 0.9143828581263734, 0.32968527382585583, 0.22606635317901605, 0.7797799519566648, 0.1825503315074003, 0.4107718363510602, 0.205479045119416, 0.7438514218438445, 0.17514865532558768, 0.4868778776803427, 0.18516074118728912, 0.7037042497895932, 0.1670443843441454, 0.5534472743624664, 0.16479331091688382, 0.6603915085793434, 0.15905654809955372, 0.30026192839909627, 0.9841560118051685, 0.285292742117286, 0.9118053754645411, 0.2688275796299732, 0.8373603726880354, 0.3941004907692921, 0.28447755923684515, 0.715318840847165, 0.24385118449213666, 0.46807136991423237, 0.25412083047024486, 0.677814119651809, 0.2316397053201795, 0.5327045857794355, 0.2288314228867426, 0.6377217886294108, 0.22080740238619792, 0.35750513234507164, 0.8139592299695135, 0.33766471835310247, 0.7543545826978385, 0.4466199433354183, 0.3348493589233822, 0.6500326045931178, 0.30602631814440157, 0.5095421416844387, 0.2977881295295873, 0.6114021780735976, 0.28874108906245205, 0.4054479032739489, 0.6775793802012414, 0.4860877210854043, 0.3774125256305269, 0.5856581662438417, 0.36426308497535914, 1.9260189539740704, 0.016896739142096075, 1.8644642012883639, 0.01664716288365133, 1.7824850475909926, 0.016319288532510072, 1.6820682824023627, 0.015732724086128032, 1.565528115609385, 0.014831963088829664, 1.4349615951543966, 0.014442479192257869, 1.2942265142051905, 0.013410364916528264, 1.1467622601757383, 0.01130838277257212, 1.8275301869038898, 0.05542781403266054, 1.7477965296845848, 0.05429789300619424, 1.6504309064587088, 0.052423969153442926, 1.5377802922450499, 0.04946146802257783, 1.4103479520687638, 0.04787249153033048, 1.2735832834487555, 0.044574685113981015, 1.1301650703594774, 0.03987713700916072, 1.6925471728489834, 0.11283137996688535, 1.5993905725086583, 0.108990480925924, 1.4926436496439215, 0.10316314830002635, 1.3716640739680208, 0.09875046761353254, 1.241060626805168, 0.09291914012525138, 1.1031663098324538, 0.08599390045395436, 1.5306779043076437, 0.18376066231745639, 1.429427928705733, 0.1748138133304228, 1.31680852752447, 0.16477160937559418, 1.192992036221837, 0.1576918914627673, 1.065209376259242, 0.14683841726772687, 1.3561890451272893, 0.26334303820584115, 1.2555865233160035, 0.24539602517020204, 1.1379877320912346, 0.23726684337456486, 1.0206205763182443, 0.21950500341316162, 1.1805578360991846, 0.3402977014644074, 1.075202109917414, 0.32746881419190976, 0.9715641532305442, 0.30356432332239647, 1.0074227794717605, 0.4242583994648633, 0.9116694974531468, 0.39767029102624096, 0.8542668216872868, 0.5004385346933509, 0.11577020778297441, 1.8121985781366472, 0.0556297583969074, 0.1967564447505714, 1.7328567542349285, 0.05446690319739, 0.29644815206184044, 1.6360555361267664, 0.05209764623037249, 0.4124078052560299, 1.5235219178288348, 0.049264422764643184, 0.5416199187784911, 1.3976442419358925, 0.04670702206450719, 0.6806852260767164, 1.261806537373057, 0.04445591987509049, 0.8260129689128608, 1.1195740292865108, 0.04239564751688782, 0.19502852790868472, 1.6766714674784207, 0.11290798415372588, 0.29211282926919413, 1.5844263382148838, 0.10824193161697221, 0.40462688042476175, 1.478300207096887, 0.10274672334484948, 0.5298137451560757, 1.3588916046805963, 0.09756759412052525, 0.6648997890968523, 1.2294972398064838, 0.09219074088547466, 0.8073364981004524, 1.0941912112972711, 0.08661596595655456, 0.284582620167658, 1.518487525744311, 0.18287182581473485, 0.39366354176963714, 1.4182188786227803, 0.17409206141487904, 0.5148693571469414, 1.3062645054176067, 0.16523012379228044, 0.6461200593687602, 1.1863771189352696, 0.15575884381986177, 0.7840920522384689, 1.0587181158399137, 0.1455906605110941, 0.3797379172405233, 1.3454758716091157, 0.2616762654167215, 0.4968503119337686, 1.242423240047397, 0.24788484993472, 0.6228105339024734, 1.1318892392556446, 0.23367046970782282, 0.7551780368652787, 1.0148847401088843, 0.2183586615053323, 0.4755206358571073, 1.168571847677879, 0.3437774381145543, 0.5959567434119198, 1.0687530995630519, 0.32413644492946736, 0.7223982739646335, 0.963532522546387, 0.30319631832456817, 0.5659756799104844, 0.9985030555424234, 0.4246932756493861, 0.6859897536257918, 0.9064916193858257, 0.39757992610546156, 0.6467182156014208, 0.8449911405016292, 0.4982450849263137, 0.1914617270291166, 1.646382044865666, 0.11091112510008974, 0.28730297264575144, 1.555480963108397, 0.10683476973517614, 0.398876900689469, 1.4512795445793014, 0.1020004086357402, 0.5225356328464124, 1.334959300776224, 0.09669922132884555, 0.6553275663720357, 1.2098766827030425, 0.09096372924847299, 0.7963178148036103, 1.0783478570732996, 0.08602528091109128, 0.28072473707788026, 1.4911947875311964, 0.180653549129523, 0.3879689271751353, 1.3930960876847436, 0.17196738837859482, 0.5075597380944057, 1.2845879919959193, 0.16339881475681467, 0.6373647092450386, 1.167609278365358, 0.15395285398373026, 0.7737831658514265, 1.0434718011315838, 0.1441264944668377, 0.37507946200649217, 1.3232058104162117, 0.2578419403547792, 0.4900231581832107, 1.2221820873637672, 0.2451677425407883, 0.6145848417537013, 1.1150073255311472, 0.2310163888189108, 0.7462418421143729, 1.0008692086341746, 0.2159449132121912, 0.4694470616858065, 1.150415456228881, 0.33972605169667563, 0.5887897471013616, 1.0528311919674616, 0.3203872330634397, 0.7129833228522775, 0.950645009055851, 0.30021602479453324, 0.5590509604661639, 0.9849066350505524, 0.4198882141724312, 0.6779438366082098, 0.8946310518894789, 0.39344175273732235, 0.6395431566487751, 0.835407682314156, 0.49262593952280664, 0.2728938075773403, 1.4510089535987478, 0.17608666836485798, 0.37948878231238686, 1.3539966481206829, 0.1687339926473916, 0.49849579181687914, 1.2492686414930163, 0.1606650229015726, 0.6251853454759514, 1.1368844322972804, 0.15092547134345588, 0.7579810062776021, 1.0196750939969352, 0.14173422267092678, 0.36684216236274053, 1.2880151463371303, 0.2530876858893762, 0.4807066012016706, 1.1906238340831556, 0.24034816712218007, 0.6029385790213954, 1.0870258869016458, 0.22653375359608277, 0.7314747272494017, 0.9787586229036458, 0.21260186018888227, 0.460129392008179, 1.1226261375540576, 0.33287838139893283, 0.5774815716120444, 1.0285927793114626, 0.3144049595761841, 0.7002565210829096, 0.9297425244755235, 0.2948906310106038, 0.5485764905234549, 0.9630994923478129, 0.4122923481403998, 0.6656791940914748, 0.8759343981722598, 0.38680919299510885, 0.6280242453845138, 0.8200172750108352, 0.48415104829725425, 0.3562928131220296, 1.2409544074444327, 0.24624368744327552, 0.46676685442987137, 1.1489682249947906, 0.23426037676370126, 0.5866843332635648, 1.0511177977199175, 0.22022962224033302, 0.7121927209156812, 0.9474503408712301, 0.20725861388396113, 0.44751333470083093, 1.0853006975794686, 0.3249404070461038, 0.5626751820800306, 0.9951765425385538, 0.3069491652250456, 0.6820171793366065, 0.902859084435128, 0.2875533059048039, 0.5346090515124726, 0.9351616160829186, 0.4020245112266916, 0.6489567147936385, 0.8519793391025159, 0.3782693397992473, 0.6134206890222683, 0.7995515843830872, 0.47238754419181084, 0.4331604338503452, 1.0388870416107128, 0.3135598510123516, 0.5425085980885869, 0.9563300279893205, 0.29718054297448504, 0.6599460987240283, 0.8687888026609281, 0.278550873744392, 0.5181643838033964, 0.8999223408732047, 0.38946358210963405, 0.6285306655731685, 0.8231814473397543, 0.36763759896517406, 0.5948263159029695, 0.7753677864312063, 0.4585296043892216, 0.49701266142179235, 0.8590731864264474, 0.3749259125642577, 0.6050970514350454, 0.7886864078023129, 0.35286937625557824, 0.5737561452770511, 0.7450749021640463, 0.44385391758674503, 0.5485398003355776, 0.7104412133425148, 0.4241517315415163])
    SymCubatures.setweights!(cub, T[1.0768878301484183e-7, 1.9279487447672886e-5, 0.00013382919519483566, 0.00011488987283086573, 0.00024169948328079988, 0.00030481109379418503, 0.0003313907296471001, 0.0004999199160118321, 0.00038844743353571367, 0.0006864136830454952, 0.00043200900734167123, 0.0007905892228546043, 0.0004618404399635471, 0.0007135298624265211, 0.0004528310805706878, 0.00047491848371295403, 0.00045160604733119177, 0.0004395342166900476, 0.0004726917927342642, 2.702982059731306e-6, 6.149209878464388e-5, 0.0001827397978437967, 0.00031247272194871625, 0.00042444716542802485, 0.0005339750489484777, 0.000538127110174706, 0.0005073898441983236, 0.0004861229700181642, 0.0005032919165458012, 3.60603721375237e-6, 1.1463532161524064e-5, 2.0959265450980906e-5, 3.0408246495974413e-5, 3.8146306059289136e-5, 4.02816739981076e-5, 3.998878299297334e-5, 2.6613104267536673e-5, 2.763283697321029e-5, 1.1850772429189469e-5, 1.7770780770420428e-5, 2.1574178162853035e-5, 1.7343946911257426e-5, 2.877079010508515e-5, 2.7159296563306917e-5, 2.9684284756858866e-5, 2.287857176892461e-5, 2.0537936451699036e-5, 6.34356310823869e-7, 1.1035567152166906e-6, 1.5157286619947306e-6, 1.847132394872013e-6, 2.0382497120596605e-6, 2.0532144429539833e-6, 2.147049297112308e-6, 1.9552835527288293e-6, 1.195676750748235e-6, 3.69549459504855e-5, 5.401842331975917e-5, 6.723374220493863e-5, 7.452171618330929e-5, 7.664688510931791e-5, 7.62656155102047e-5, 7.444673217356701e-5, 7.220399543625622e-5, 6.853462516964323e-5, 8.94314770742328e-5, 0.0001293914597520049, 0.00013212138577659274, 0.00018351084776428592, 0.00015694886295578402, 0.00023079196283892504, 0.00016507211267224848, 0.0002458106619785616, 0.00018066350810361334, 0.000232914127480728, 0.00016622059488405133, 0.00019781735978802804, 0.0001522417084708289, 0.00015939177689948186, 0.00014143289758144755, 0.00015974685380618997, 0.00019313045490414052, 0.00021477230857824137, 0.00022783547860839826, 0.00023081904978186117, 0.00022147248231878285, 0.00021161565002753658, 0.0002233123331097015, 0.00022626138297605826, 0.000318741559532457, 0.00026658326975517884, 0.00038867693162751305, 0.00029595506867869113, 0.0004192147915857368, 0.0003117418412086016, 0.0003907971089873783, 0.0003016989721351929, 0.00034908485255303657, 0.0002699648651653721, 0.00027283969981098995, 0.0002489158396716293, 0.0003592011390267991, 0.0003929510159377007, 0.0004119443121842644, 0.00040913115149673933, 0.0003806615491160814, 0.00037423903475206554, 0.0004295621225276752, 0.00036821116728429537, 0.0005248994935348457, 0.0003858100094051777, 0.0005549011144530277, 0.0003966488361082117, 0.0005417703335183296, 0.0004100533899966452, 0.00046217626389282037, 0.0003711315480540256, 0.00036610013148675026, 0.0003478140827962338, 0.0005800347181637052, 0.0005962525640006816, 0.0006097996662447274, 0.0005416215498241498, 0.0005083087411947196, 0.0006360378828797298, 0.00048441789710666517, 0.0006634942662177376, 0.00046993364803653536, 0.0006395214898278357, 0.0004725347068629482, 0.0005509993632448457, 0.00045358171189927433, 0.0004365284305954249, 0.00040882652675875927, 0.0007188378725077856, 0.0006992659632594282, 0.0006703783284188955, 0.0006012721890183396, 0.0007599228545374989, 0.0005179383621237922, 0.0006932530266331406, 0.000531807852104478, 0.0005661618212377434, 0.0004923217638074252, 0.0004919673648796802, 0.00047224275700392975, 0.0007676705433416923, 0.0006737809236596144, 0.0006406048999695179, 0.0007142855688951093, 0.000540632973064918, 0.0005785483106366463, 0.0005338747057862019, 0.0004897772955712028, 0.0005209507036053521, 0.000665664359829141, 0.000636350653007317, 0.0005768038994786533, 0.0005304845290919822, 0.0004848818763190974, 0.0005243644439841112, 0.0005242638502617119, 0.000497745905391426, 0.0005611474939621888, 6.521135070504003e-6, 9.207860964403337e-6, 1.1264732454629012e-5, 1.2344546963893246e-5, 1.2490570644116533e-5, 1.255892085140708e-5, 1.202143664931723e-5, 1.0978505077314173e-5, 1.578838729273589e-5, 1.8996386096290223e-5, 2.119637675525467e-5, 2.1963377194777013e-5, 2.1979942611017074e-5, 2.0295407744251782e-5, 1.8080272399531747e-5, 2.5439352448881615e-5, 2.8841562263798238e-5, 2.927251554096104e-5, 2.8610055226766142e-5, 2.933052729668441e-5, 2.55574696397596e-5, 3.2510382057166576e-5, 3.5362799046393164e-5, 3.492171198349051e-5, 2.93355729607732e-5, 2.959410390772986e-5, 3.781345190201151e-5, 3.6913580028723115e-5, 3.465766875030565e-5, 3.3003565945027633e-5, 3.8689808766704285e-5, 3.18661398146085e-5, 3.3225311334361227e-5, 3.5099397527290095e-5, 2.8505108809145053e-5, 3.1892187752834275e-5, 2.3042183148775606e-5, 9.664185953782125e-5, 0.00011732044852319015, 0.0001295273555825144, 0.00013406759100622416, 0.00013170659806499563, 0.00012292290573716718, 0.00011021154535944736, 0.00015517457994805383, 0.0001744656520789207, 0.0001780754224654782, 0.0001774878184956212, 0.00016990498516153914, 0.00014511453513631865, 0.00019963102220319653, 0.0002147183722006714, 0.00021118702809528705, 0.0001907563354613396, 0.00017451346302960454, 0.00022841376086937633, 0.00022861940940286884, 0.00020597487070169116, 0.00019306498318697412, 0.0002311821700017202, 0.0002154294442786894, 0.00019592309038229738, 0.00020519456318951958, 0.00018655711722664862, 0.00016815896197908508, 0.0002687762248814082, 0.0003000216563468613, 0.0003094673295152433, 0.00030424915175656476, 0.00028577309397654786, 0.00025557289557795785, 0.00035264315490815247, 0.00036830069148332465, 0.0003553580620324632, 0.00034655402753741804, 0.0003010098855618292, 0.0003964927354693064, 0.0003890197897529844, 0.0003716801691582825, 0.00032040960843397764, 0.0004041513651157397, 0.0003738524605671413, 0.0003420836297611925, 0.00036907961409571896, 0.0003231823001456337, 0.00029436765339503537, 0.0004867058965646404, 0.0005056706960204164, 0.0004837450173630336, 0.0004566670869003791, 0.0004156027377037701, 0.0005455556428205308, 0.0005279746076505549, 0.000505196933971419, 0.0004582201656412018, 0.0005423048149618452, 0.0005121577245124286, 0.0004578492592026912, 0.0004966323973529231, 0.000444728313497226, 0.00039528112671854244, 0.0006450947240199076, 0.0006500694890247152, 0.0006060551567271032, 0.0005491283943380436, 0.0006496700235194724, 0.0006204876424445784, 0.0005589400376082861, 0.000577366041770339, 0.0005146355771982385, 0.0004659050656028497, 0.000718973681994016, 0.0006855709048884845, 0.0006218821050164441, 0.0006446183738621824, 0.0005822157803769427, 0.0004984357019365748, 0.0006968332132453695, 0.0006388397600920674, 0.0005514316899476489, 0.000597573162698055, 0.0003050873572930493])
    cub_degree = 1
    tol = 7e-15
  elseif q <=40 # 5324 nodes
    SymCubatures.setparams!(cub, T[0.022199827847323, 0.07478252550423985, 0.15106939339241093, 0.24225356859529124, 0.344556382408867, 0.44305591027101804, 0.5331739439552663, 0.6117145068888113, 0.6824715154448738, 0.7290458345784716, 0.015520195916893246, 0.051361835022348444, 0.10586161437431742, 0.17545300117935672, 0.25619669042302295, 0.34284594418975295, 0.4293872002621678, 0.5124420709366497, 0.5900107185979331, 0.6481200613115444, 0.992076219228823, 0.9736021419996144, 0.9450311450954523, 0.9069744638059606, 0.8602436199806011, 0.8058347191421296, 0.7449074375949511, 0.6787603550694598, 0.6088032925796425, 0.5365272700054492, 0.015058736719031067, 1.9204399763780426, 0.01519188343478328, 1.866210910348032, 0.015014877545343759, 1.7939622158837594, 0.014476814762625327, 1.7048348022868012, 0.013739496463927005, 1.6000000125214235, 0.013059219475206904, 1.4808491913432686, 0.012506699971488466, 1.3495476929541064, 0.011892776824890374, 1.209268784805901, 0.011185607631157739, 1.063231044757345, 0.050354266120189284, 0.015037919765395211, 0.10541388810451569, 0.014454289015742128, 0.17479169490513802, 0.013425042637498124, 0.2531026435274504, 0.01302646603880126, 0.33693169344275253, 0.012063564668226045, 0.4217438702338358, 0.010928348182659166, 0.5028651340405039, 0.009855326550255238, 0.5748239090635932, 0.008879566684262187, 0.6354744078983405, 0.00768699303377101, 0.04979302146926438, 1.7973592362358215, 0.04895203574962178, 1.7290723950626996, 0.047411915613733736, 1.6465205218454921, 0.04550899471400995, 1.5492143058897037, 0.04348606422075897, 1.437301375548764, 0.041240125103594455, 1.312336636956345, 0.03841283671629551, 1.1777972907276517, 0.035461512425585505, 1.0374382927044243, 0.10365423373129168, 0.048239843544192695, 0.17225774920617767, 0.04488924520193362, 0.24977824997810708, 0.043384870356906934, 0.33311799968983613, 0.04031448112305941, 0.4169818406148851, 0.036828939165659706, 0.49758015316466897, 0.03333757886257213, 0.5698550158812379, 0.029527471074885306, 0.6304334352267259, 0.026944551448796462, 0.09890312239888267, 1.6317190096213567, 0.09601646204495269, 1.552183728043746, 0.09240726733929444, 1.460566785796808, 0.08907928241707506, 1.3566751108305009, 0.08471715703519032, 1.2439059641636288, 0.07834265627815076, 1.1246085848630047, 0.07279291614288, 0.9962969082235592, 0.16767537096396717, 0.09434202774148401, 0.2444965951575429, 0.09067081997213261, 0.32704151963903005, 0.084059803039886, 0.40959464819377206, 0.07667623252850034, 0.4891626019914355, 0.06960922619952042, 0.5603199739807515, 0.06301509454016181, 0.620470717039661, 0.05732252569675701, 0.15830375998204688, 1.4393067864084046, 0.1525813897206073, 1.353837496055499, 0.14524251699413643, 1.259528518987923, 0.13820731126641234, 1.1567820632002621, 0.13001160922525964, 1.0507089405403305, 0.12138412066028134, 0.9399226179104404, 0.23765512101564198, 0.15366911563117455, 0.3177293904083171, 0.1418858332532202, 0.3995280382335048, 0.12968477889098137, 0.4772363072415848, 0.11822388446711149, 0.5487850025221668, 0.10703904145211898, 0.6091048705716396, 0.09711836244324593, 0.2208179313197539, 1.2407758097924115, 0.21238896870203514, 1.1582171587548935, 0.2021844708105463, 1.0685397275980957, 0.18966308756996053, 0.9758057575586496, 0.1784369872383618, 0.8774034063879292, 0.3084530536106153, 0.21187479937657572, 0.38783934541251774, 0.1944346834326873, 0.46335622663090287, 0.17712295713812085, 0.5333653080190562, 0.16089723853271193, 0.5932416896568782, 0.14596051749268643, 0.28264680417221233, 1.044911871806447, 0.2703984223561586, 0.9699334869071983, 0.25556719214749163, 0.8931764574786264, 0.24064396999702922, 0.8104455177314845, 0.3728625706213381, 0.2694609481746303, 0.4474184622542271, 0.24555406332412624, 0.5164150247583206, 0.2230761854449754, 0.5752580250561184, 0.20238290426892896, 0.33954331163685103, 0.8708146509368562, 0.32224690425060376, 0.8067438914073636, 0.3055765941994229, 0.7398220552657658, 0.4288580606826542, 0.32260591048171205, 0.4970612024308597, 0.2917180045997581, 0.5541556721046997, 0.2665658077444518, 0.38983378662191576, 0.7206819921701644, 0.36871454617675103, 0.6671935798589851, 0.47577636051242855, 0.371482490224035, 0.533559031276394, 0.3336680473729438, 0.430599282804034, 0.6001254672177957, 0.5105238032298987, 0.4101381668396862, 1.9328122124888099, 0.015422948202503575, 1.87676324921528, 0.015270747393298109, 1.8019422285286835, 0.014990275739266935, 1.7098170420035164, 0.014588484060144732, 1.602255837288058, 0.014090639605248229, 1.4815519124855931, 0.013457665634837926, 1.3501924770536882, 0.012837971455718829, 1.2110827019324186, 0.012142854715495081, 1.0673633108293215, 0.011137468970963366, 1.8422260290816894, 0.05086522748851043, 1.7691905317471874, 0.04991204833554682, 1.679553136565619, 0.048595097268011346, 1.5750337031979533, 0.04695383891511884, 1.4578292464741107, 0.04480631014918908, 1.3298091807397527, 0.042749167202004806, 1.1939138360347499, 0.04046860208607936, 1.0533422877744056, 0.03740904514137683, 1.7171768822425573, 0.10380633717448347, 1.6302624641466859, 0.10107928806113173, 1.529519425842615, 0.09767473292024843, 1.4174823173481104, 0.09311246760488218, 1.2951201915526183, 0.08879797208386385, 1.1651674058219412, 0.0842628238297553, 1.0306068277786176, 0.07859062752831351, 1.5659628061625543, 0.17083588671670435, 1.469889065097346, 0.16507458674346523, 1.3639545366717825, 0.15726725255304355, 1.2486795194846063, 0.14963869803640817, 1.1261406821085884, 0.1425712911966415, 0.9996827145834297, 0.13399948780146256, 1.3973616702595404, 0.2474868638716685, 1.2985419024092064, 0.23589110420097514, 1.1918136001921487, 0.22362199224160173, 1.077658513621494, 0.2138944472861949, 0.9611926851714845, 0.2023080066841032, 1.223718742192, 0.32739888034830167, 1.1273516141147384, 0.3095991276917263, 1.022908057505157, 0.29617571170406565, 0.9167704298899507, 0.28157897118524744, 1.0546955259769366, 0.4065760213185983, 0.9626075774216891, 0.38721344016843695, 0.8683483205798233, 0.36893446854103373, 0.8932110316253533, 0.48541232469064305, 0.8147171532779959, 0.4603246894147091, 0.7573104809598555, 0.5543869544469994, 0.10535630351764214, 1.8290663497449329, 0.05062608490653875, 0.17922542009611936, 1.7562420525640632, 0.04986784130762788, 0.27062459418544754, 1.6671689035824795, 0.0480295016652198, 0.37764442351882926, 1.5631770965491423, 0.04560352519492038, 0.49771222087465083, 1.4458974968962859, 0.043393071202662024, 0.6277399362060965, 1.3184369283086443, 0.04150507593045156, 0.7647286436184453, 1.1847736247145377, 0.03917890468998616, 0.9057105290704494, 1.0468297271392732, 0.03698489904129193, 0.17828273533082012, 1.7039784656399655, 0.1035485125672938, 0.26764475812535393, 1.6186775230994943, 0.09984821281813563, 0.37167555756017434, 1.520031484468149, 0.09503634430562215, 0.48810180159143624, 1.4084999785334058, 0.090507871010419, 0.6145542872504612, 1.286904789467667, 0.08631554794894525, 0.7489220625483608, 1.1587489760934113, 0.08105902818046039, 0.8878890370224347, 1.0257043143149194, 0.0763075724280058, 0.26195394959896123, 1.5557983287096953, 0.1686966114213096, 0.3632500510400304, 1.4629060603786623, 0.1610078556173752, 0.47613237084632176, 1.3579150342829454, 0.15336105008540588, 0.5986172093858025, 1.2441935427789772, 0.14587218206720298, 0.7285358595772634, 1.1235857533209401, 0.1371056715151547, 0.8640672406473606, 0.997609641222369, 0.12816130587975383, 0.35114671692604776, 1.3940476433823343, 0.2422654115211699, 0.4606243104200303, 1.2967151259869443, 0.2307305328318151, 0.579330456189484, 1.19101157717531, 0.21869811647961995, 0.704821368303868, 1.0787568964268865, 0.2060278275197527, 0.8356479780548125, 0.9622865865002376, 0.1922101407398799, 0.44273476959104063, 1.225175795403263, 0.3208536917980045, 0.556787703623731, 1.1287145681763862, 0.3035714728597984, 0.677051752839557, 1.0268504339096391, 0.2857705509894095, 0.8022471607158146, 0.9207811914897205, 0.2674107823599544, 0.5309423253605858, 1.0595246303511765, 0.3991584103501594, 0.6458156443057542, 0.968675170258828, 0.37584635843993763, 0.7648269972258516, 0.8748216533781237, 0.35124136465831324, 0.6123736351480571, 0.9054110405059101, 0.4729302069848738, 0.7250111311038234, 0.824103162778729, 0.44237309467146996, 0.6817483538433291, 0.7714291213775442, 0.5386018482731619, 0.1754555165814396, 1.6753495027485896, 0.10180755552236068, 0.2638897011770189, 1.5917961493519988, 0.09830469454363241, 0.36696622322485867, 1.4944025976790263, 0.09437801032684583, 0.4820172227273116, 1.385365456300101, 0.08981194928906533, 0.6065378420192457, 1.2675534149325638, 0.08522836548515447, 0.7389494520857235, 1.1434720177510271, 0.08018478095867412, 0.877214110948095, 1.0138586103077392, 0.0746440619373313, 0.2581985754434773, 1.5303883055929413, 0.1662072778767396, 0.35812853669336087, 1.4391630879246033, 0.15983216982823117, 0.4697879587473987, 1.3366388592384881, 0.151896980285559, 0.5914530873837477, 1.22611324952043, 0.1440828928006998, 0.7201942907404371, 1.107963654880301, 0.1356427310154304, 0.8539289818028355, 0.9848831816504896, 0.1270976384341323, 0.3469438574666096, 1.3715724136137124, 0.23997705841956665, 0.45518756464428917, 1.2764927123989755, 0.22836344137279657, 0.5724464829279465, 1.1741958294703958, 0.2163009155096108, 0.6964921395375564, 1.0647631926116452, 0.203603911056667, 0.8258076508680118, 0.9505735894374236, 0.1909134005729043, 0.4376815194668834, 1.207180106102637, 0.31751517637129034, 0.5504646142044854, 1.1129723265146672, 0.30037775457925087, 0.6690523045052623, 1.013748466689661, 0.28338634924286993, 0.7930154757602726, 0.9103006566207881, 0.26457313303474583, 0.5251600103355804, 1.045436794141143, 0.3949647257399785, 0.6385830250258058, 0.9568272263823749, 0.372061395657747, 0.7563982144628548, 0.8650174074041627, 0.34786646587938325, 0.6062356872647028, 0.8950079741869693, 0.46811105382759816, 0.7169530503798512, 0.8153925151053393, 0.43810710028145416, 0.674622971264058, 0.7636570736231301, 0.5337491819018174, 0.2519635669494599, 1.4915856713342799, 0.16258785060684708, 0.3508741821256997, 1.4030517605835133, 0.15622270890041692, 0.46142697955544637, 1.3027680616293857, 0.1497768668342465, 0.5813765222119195, 1.1955281752912421, 0.14212021058689717, 0.707994284496233, 1.0828787355973242, 0.13290457669484573, 0.8387056037740831, 0.9655014802094626, 0.12452601949282513, 0.3400295920942309, 1.3388281746695603, 0.23529157284947272, 0.4466672657370754, 1.2455982407538457, 0.22410065330823517, 0.5623901284549055, 1.1463793830889013, 0.21289255154911615, 0.6841533711238237, 1.041494361480246, 0.20076068501840735, 0.8109847231622332, 0.9319799308069359, 0.18816140493471226, 0.42976168284598343, 1.1797857807484204, 0.31136418722231207, 0.540715461872331, 1.0890359672295955, 0.29533537118547243, 0.6578709127219846, 0.9927169453547844, 0.278749700368811, 0.7790772039151991, 0.8929630460755866, 0.2612225537755143, 0.5159563082093906, 1.0236647403548895, 0.3881070015887864, 0.6281075100745023, 0.9382906629161287, 0.3659837582710967, 0.7440066216642981, 0.8490860151580677, 0.3420848740604499, 0.5961647018400976, 0.8790602685206745, 0.4606993869631361, 0.7050738885043755, 0.8013095528614327, 0.43112099884934557, 0.6641891433753024, 0.7505539230149688, 0.5261739827119455, 0.3310901106497247, 1.294449231114047, 0.22904278121463964, 0.43479059804283754, 1.2049148527780387, 0.2186844618725084, 0.548144150726483, 1.110111989305712, 0.20794051004119254, 0.6682123491532278, 1.0103469282019755, 0.19612407507166685, 0.7912228471217572, 0.9073094411507153, 0.18354518726381702, 0.4185981916330194, 1.1423600031878698, 0.30427316754576555, 0.5279175888303262, 1.0556366292768808, 0.28871305284032184, 0.6424501825262324, 0.9646786822581822, 0.2726273039422319, 0.7608384629319399, 0.8697037511126379, 0.25582433698896395, 0.5041460813056181, 0.9945349062789647, 0.3791112715698882, 0.6137678329651999, 0.9134901079028793, 0.3574084519121725, 0.7261340831458876, 0.8285127794381959, 0.33547697870113086, 0.5822413368688307, 0.8570408659590153, 0.45078730294173325, 0.6895996318795011, 0.7830641097668527, 0.42206512769026155, 0.6494296823253827, 0.7340640609093809, 0.5155113536495927, 0.40548210014384245, 1.0968148290828885, 0.2946295453803807, 0.5105597483243456, 1.015084891361534, 0.28043974781432296, 0.6230876467161637, 0.9307296363071579, 0.26410617379596424, 0.737461215866017, 0.841305653935221, 0.24874733036772986, 0.4890053874870176, 0.9590784599469842, 0.36796875787246197, 0.5954216701015541, 0.8819206216288656, 0.34810759534490837, 0.7052586023479718, 0.8027943218913044, 0.3271835266675848, 0.5653462027895884, 0.8300302953448745, 0.438485154732223, 0.6702040029848944, 0.7597447449467837, 0.4113871477163099, 0.6327026591764323, 0.7141603435179287, 0.5022820207779932, 0.4715103780225529, 0.9163340400571376, 0.3546411264642131, 0.5743616069761267, 0.8461307249867823, 0.335996006471769, 0.6804183585354865, 0.7722864553964527, 0.316292646835758, 0.5468630545802307, 0.7983363079540797, 0.42433659207999985, 0.6467094153233678, 0.7324713645516923, 0.39977907748730124, 0.6133294361598335, 0.6901286749775495, 0.48702197062855906, 0.5249399097977397, 0.7610709831646563, 0.40756864330261644, 0.6219983548183591, 0.7023025085336244, 0.3845972075533262, 0.590449328783493, 0.663122653743548, 0.47078424816620484, 0.5658154539955681, 0.6331450699831942, 0.4504943661976875])
    SymCubatures.setweights!(cub, T[8.099357962650956e-8, 1.5285618162618987e-5, 8.676452833338081e-5, 0.00023066554449620224, 0.00042205256462239427, 0.0005281786791954733, 0.0006741079601189232, 0.0006673381843461289, 0.0005899962602208206, 0.0004312059236587771, 0.0002465943363077258, 2.7348772222432073e-6, 8.673872436270056e-6, 1.6317969396959142e-5, 2.3417150710049823e-5, 3.0256551407673075e-5, 3.2765176215497654e-5, 3.176799560272849e-5, 2.8890123706277612e-5, 2.3351564156960867e-5, 1.0056624760007575e-5, 4.751380822847034e-7, 8.374684685107456e-7, 1.1629419720658885e-6, 1.425711308293452e-6, 1.6096842811342544e-6, 1.7102087910686022e-6, 1.7194506594021565e-6, 1.6800630191121523e-6, 1.5724875644607408e-6, 1.3332343543776101e-6, 2.833068561123775e-5, 4.10880877895661e-5, 5.1529820202868786e-5, 5.783701030944102e-5, 6.018725154267533e-5, 6.066926114135691e-5, 6.0213920339878534e-5, 5.690749559103444e-5, 5.1764084445950796e-5, 5.146465146678093e-5, 9.939585513678845e-5, 0.00014161211197280424, 0.00018028379605741852, 0.00019999827647687992, 0.00019646314646224779, 0.0001768067243052076, 0.0001384188716158169, 0.00010572180018108306, 0.00012236718308723088, 0.00014912986870061815, 0.00016716032335248863, 0.00017968699832080534, 0.0001840719738386483, 0.00018150759309337343, 0.00017114130524235733, 0.00015062948480308255, 0.0001715027746065851, 0.00024802310766317345, 0.0003155078534731524, 0.0003456947361059327, 0.00034149289488753456, 0.00030970710359693387, 0.0002488046574104587, 0.00018991627671913658, 0.0002812674092721895, 0.00031409072357766447, 0.0003283172131675627, 0.0003327774524017609, 0.00032380216064203487, 0.0003039806181719344, 0.0002834583373082634, 0.0003417324432454158, 0.0004307463012516493, 0.000469683742440004, 0.0004637109693858714, 0.0004202182370291312, 0.0003512754476446315, 0.0002549841496966802, 0.00046061123595421253, 0.0004925309363309547, 0.0004878628501337113, 0.00046825764140353686, 0.00043035563020407756, 0.0003885226264223902, 0.0005129389399223852, 0.0005563764903908291, 0.0005551962448068527, 0.0005040749346020444, 0.0004073923126167009, 0.00029697645616305466, 0.0006110334840223968, 0.000589068912197084, 0.000582041369923114, 0.0005343191829884704, 0.0004802112274163804, 0.0006264990316445798, 0.0006122340771696744, 0.0005500562987430913, 0.0004666440254610062, 0.00033164260081326366, 0.0006498510533042224, 0.0005996892437437464, 0.0005849134480589704, 0.0005183629367825912, 0.0006549424139512169, 0.0005890550825164644, 0.00047692080035629646, 0.0003616690439069085, 0.000618627073194044, 0.000562559769860783, 0.0005136672458642031, 0.0005853078114082747, 0.00048389408689938635, 0.00035643516762681586, 0.0005537865563280545, 0.0004493874001723517, 0.0005025905265219626, 0.0003450619847378998, 0.0003740878477308911, 0.0003300270668284514, 4.923271654814432e-6, 6.988799253072274e-6, 8.646537207616682e-6, 9.685140338901805e-6, 1.0140135029400883e-5, 1.0168354331686587e-5, 9.96922444560236e-6, 9.423478294716586e-6, 8.368803003366873e-6, 1.2068139799288586e-5, 1.4766322280119653e-5, 1.6602138974539604e-5, 1.756505957580225e-5, 1.7744624839100828e-5, 1.7269423889723988e-5, 1.5619116573400127e-5, 1.379028688883717e-5, 1.9907314187928127e-5, 2.2527645721542486e-5, 2.3609567758604013e-5, 2.4152555391085747e-5, 2.326609403021199e-5, 2.1541219352104883e-5, 1.849310586189742e-5, 2.7231880876901997e-5, 2.8296508162566623e-5, 2.8941312773347616e-5, 2.576480091721282e-5, 2.4701711998353663e-5, 2.2813867681142024e-5, 3.1734841545934864e-5, 3.119152319795262e-5, 2.838930530502212e-5, 2.648542630172424e-5, 2.4998144303271345e-5, 3.1879579701022764e-5, 3.072395568296874e-5, 2.8659967841361275e-5, 2.4719122592219584e-5, 3.076446874203842e-5, 2.5342116166538322e-5, 2.341361463966437e-5, 2.6296765986848972e-5, 1.8568195717164886e-5, 1.540861271360577e-5, 7.32395935870186e-5, 9.052010932252854e-5, 0.00010123781917429168, 0.00010549912820745537, 0.00010575781035060808, 0.00010173791306656416, 9.156721333824598e-5, 8.156696572175628e-5, 0.00012121823569188276, 0.0001353955270030282, 0.00014155517351223097, 0.00014376933895810263, 0.0001385779891718405, 0.0001244381528665675, 0.00010733500415327773, 0.00016403811206773148, 0.00016889382086128778, 0.00017215712991607756, 0.0001575938658080287, 0.00014747034560927926, 0.00013228182213895461, 0.00018932942444563095, 0.00018840791325993838, 0.0001737243625005562, 0.0001622189396138637, 0.00014512272068123693, 0.00019283889359662546, 0.00018641211175308438, 0.0001680816734817016, 0.00015298510213763178, 0.00018273108091474493, 0.00016779336406325823, 0.00014756795175772575, 0.00015398374273906291, 0.00013737589809222157, 0.0001240318510134958, 0.00021055056191566127, 0.00023470724886823997, 0.00024766454212632574, 0.0002469860164978328, 0.00023733706189302007, 0.00021670492066703614, 0.0001970003030411766, 0.00028342637799126403, 0.0002951515399969489, 0.0002947020753366659, 0.00028144522229109917, 0.00025856945431278184, 0.00023435682835808142, 0.0003233504952080224, 0.00032889642557506885, 0.0003094426356563159, 0.00028596445945594344, 0.00025204769564516103, 0.0003380252257490662, 0.00031727319875610806, 0.00029118414419577164, 0.00026446642864937813, 0.00032053910399326886, 0.0002853080695433358, 0.00026453633785297596, 0.00026753959452402963, 0.00025050651226947, 0.0002206900888664173, 0.00038405444188315154, 0.0004088301005170719, 0.00039695696966044424, 0.0003878580681188495, 0.0003606295185298699, 0.0003201921642807433, 0.0004414945094317493, 0.00044701276207353594, 0.00042592037762958785, 0.000390324270346898, 0.00035549161814783624, 0.0004654668694247773, 0.00043439498040514146, 0.0003999745562842932, 0.0003548835235215982, 0.000438430396800801, 0.00039361473864428195, 0.00035661445748634416, 0.00037103271628289767, 0.00033836909845246286, 0.0003075020632495842, 0.0005351815013044161, 0.0005388759031398061, 0.0005203425894605272, 0.0004753639418292483, 0.0004292512609236976, 0.0005660602172498112, 0.0005393994688376891, 0.0004826618645231336, 0.0004406752812582228, 0.0005197202783794025, 0.0004785605665399891, 0.00041987707747019663, 0.00046678314615397485, 0.0004019618371100607, 0.00036995683198477606, 0.0006200014980027648, 0.0005999257196798463, 0.0005567559810970381, 0.0004927562637228745, 0.0005854247152743407, 0.0005398796931747147, 0.0004741282748919356, 0.0004991462417402437, 0.0004553253426283307, 0.0004044490669296852, 0.0006235052404755572, 0.0005629225726913256, 0.0005193008274506842, 0.0005345844501465073, 0.00047484719476894357, 0.0004222338495940242, 0.0005787025030494676, 0.0004794377266221715, 0.00042939476302579694, 0.00043374936938922915])
    cub_degree = 1
    tol = 7e-15

  else
    error("polynomial degree must be <= 40 (presently)\n")
  end

  mask = SymCubatures.getInternalParamMask(cub)
  append!(mask, (cub.numparams+1):(cub.numparams+cub.numweights))
  vtx = T[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
  Cubature.solvecubature!(cub, cub_degree, mask, tol=tol, hist=false)
  return cub, vtx
end

function getTetCubatureStaggered(q::Int, T=Float64; tol=10*eps(typeof(real(one(T)))))
  @assert(q>=2)

  cub_omega_degree = q-1
  cub_diage_degree = q
  mask_omega = zeros(Int64, (0))
  mask_diage = zeros(Int64, (0))

  if q<=2
    #3 nodes 
    cub_omega = SymCubatures.TetSymCub{T}(vertices=false, numS31=1)
    SymCubatures.setweights!(cub_omega, T[1/3])
    SymCubatures.setparams!(cub_omega, T[3/5])
    cub_omega_degree = 1

    #10 nodes
    cub_diage = SymCubatures.TetSymCub{T}(vertices = false,
                                          midedges = true,
                                          centroid = false,
                                          facecentroid = false,
                                          numedge = 0,
                                          numfaceS21 = 0,
                                          numfaceS111 = 0,
                                          numS31 = 1,
                                          numS22 = 0,
                                          numS211= 0,
                                          numS1111 =0)
    SymCubatures.setparams!(cub_diage, T[0.6])
    SymCubatures.setweights!(cub_diage, T[0.15151515151515163, 0.12121212121212123])
    cub_diage_degree = 2
  elseif q<=4 
    #10 nodes
    cub_omega = SymCubatures.TetSymCub{T}(vertices=false, numS31=1, numS22=1)
    SymCubatures.setparams!(cub_omega, T[0.3540971084487333, 0.18203145362159834])
    SymCubatures.setweights!(cub_omega, T[0.15107940134333817, 0.1215026213266635])
    cub_omega_degree = 3

    #33 nodes
    cub_diage = SymCubatures.TetSymCub{T}(vertices = true,
                                          midedges = true,
                                          centroid = true,
                                          facecentroid = false,
                                          numedge = 0,
                                          numfaceS21 = 1,
                                          numfaceS111 = 0,
                                          numS31 = 1,
                                          numS22 = 1,
                                          numS211= 0,
                                          numS1111 =0)
    SymCubatures.setparams!(cub_diage, T[0.3540971084487333, 0.18203145362159834, 0.3771609693928902])
    SymCubatures.setweights!(cub_diage, T[0.004928516904817574, 0.041493887927438966, 0.004844353491327435, 0.07430723352373057, 0.03913863272373876, 0.20307059922909518])
    cub_diage_degree = 4
  elseif q<=6
    #20 nodes
    cub_omega = SymCubatures.TetSymCub{T}(vertices=false, numS31=2, numS211=1)
    SymCubatures.setparams!(cub_omega, T[0.27308286881050636, 0.931780068315618, 0.09856058998260914, 0.7856072614781071])
    SymCubatures.setweights!(cub_omega, T[0.0905678982461762, 0.14714726910635742, 0.03187272199359993])
    cub_omega_degree = 5
    
    #65 nodes
    cub_diage = SymCubatures.TetSymCub{T}(vertices = true,
                                          midedges = false,
                                          centroid = true,
                                          facecentroid = false,
                                          numedge = 1,
                                          numfaceS21 = 2,
                                          numfaceS111 = 0,
                                          numS31 = 3,
                                          numS22 = 0,
                                          numS211= 1,
                                          numS1111 =0)
    SymCubatures.setparams!(cub_diage, T[0.9292686931229018, 0.12323052758795768, 0.44753553811036784, 0.8506802519794945, 0.23722737279318576, 0.3077459416259917, 0.12199561057672483, 0.5329607392794149])
    SymCubatures.setweights!(cub_diage, T[0.00031819609592614824, 0.09650100779728245, 0.010804402645722555, 0.02413806959723523, 0.00938437497734488, 0.0022614523416125557, 7.752351481171119e-5, 0.050604466827047646, 0.05835281685886672])
    cub_diage_degree = 6
  elseif q<=8
    #38 nodes
    cub_omega = SymCubatures.TetSymCub{T}(vertices=false, numS31=2, numS22=1, numS211=2)
    SymCubatures.setparams!(cub_omega, T[0.5821093910011628;0.18262906602002502;0.17262216834048155;0.08099388793344552;0.5908415591286749;0.4612405506278837;0.07057889678019784])
    SymCubatures.setweights!(cub_omega, T[0.0680746674820208;0.028423617986180167;0.018310980956037618;0.025461800630426842;0.044327724846598505])
    cub_omega_degree = 7

    #115 nodes
    cub_diage = SymCubatures.TetSymCub{T}(vertices = true,
                                          midedges = true,
                                          centroid = true,
                                          facecentroid = true,
                                          numedge = 1,
                                          numfaceS21 = 1,
                                          numfaceS111 = 1,
                                          numS31 = 4,
                                          numS22 = 2,
                                          numS211= 2,
                                          numS1111 =0)
    SymCubatures.setparams!(cub_diage, T[0.5821093910011628, 0.18262906602002502, 0.3523929920095433, 0.9444339788155788, 0.17262216834048155, 0.2961158854913302, 0.16099183834007516, 0.800367892880542, 0.08099388793344552, 0.5908415591286749, 0.4612405506278837, 0.07057889678019784, 0.6058255660767269, 0.21518364356973504])
    SymCubatures.setweights!(cub_diage, T[9.290013171024362e-5, 0.013992807836720857, 0.0015333081173549005, 0.04511998898013164, 0.023503817647965487, 0.0011541612891429204, 0.027095118471877604, 0.03139349633821964, 0.0057388548217078075, 0.000666580936171083, 0.011601819303286062, 0.019128527258902996, 0.006356360583123295, 0.004628272425307567, 0.021810254345355308])
    cub_diage_degree = 8
  elseif q<=10
    #59 nodes
    cub_omega = SymCubatures.TetSymCub{T}(vertices=false, centroid=true, numS31=1, numS22=0, numS211=5)
    SymCubatures.setweights!(cub_omega, T[0.06845875962635535, 0.03996133054703136, 0.004908706029927262, 0.012435872193517524, 0.018315521918100697, 0.011842283605048317, 0.009933723304410471])
    SymCubatures.setparams!(cub_omega, T[0.5659668188394116, 0.7807167565199028, 0.0846814812674472, 0.04046668289033655, 0.20685562107102665, 0.061511434675817794, 0.6927211311947881, 0.14342099002644818, 0.3786528536651681, 0.3596513879533776, 0.02013276646659172])
    cub_omega_degree = 9

    #189 nodes
    cub_diage = SymCubatures.TetSymCub{T}(vertices = true,
                                          midedges = true,
                                          centroid = true,
                                          facecentroid = true,
                                          numedge = 1,
                                          numfaceS21 = 3,
                                          numfaceS111 = 1,
                                          numS31 = 3,
                                          numS22 = 1,
                                          numS211= 7,
                                          numS1111 =0)
    SymCubatures.setparams!(cub_diage, T[0.5659668188394116, 0.9032995869665994, 0.21639321563468925, 0.17463991159645564, 0.40849548761177584, 0.1693510205976773, 0.895541492424858, 0.9192858150590572, 0.7807167565199028, 0.0846814812674472, 0.04046668289033655, 0.20685562107102665, 0.061511434675817794, 0.6927211311947881, 0.14342099002644818, 0.3786528536651681, 0.3596513879533776, 0.02013276646659172, 0.10927759244854979, 0.4640791974658525, 0.4272982512809475, 0.16349869355999422, 0.5199457984883094, 0.07662224108118755])
    SymCubatures.setweights!(cub_diage, T[0.00013994493674627766, 0.009574509298507828, 0.021235812009821727, 0.007520872044858897, 0.0006245093964225804, 0.02060712701121673, 0.0030696098820162847, 0.0016902128933026987, 0.005162612925040695, 0.0002184046985262755, 0.012346273603626279, 0.0021689232147176474, 0.005390084892914799, 0.008896922246407089, 0.005698467517511197, 0.007585856549682867, 0.02652315981352278, 0.0011464801433576754, 0.007633016738826831, 0.049005032484643295])
    cub_diage_degree = 10
    tol=5e-15
  end

  vtx = T[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
  mask_omega = 1:(cub_omega.numparams+cub_omega.numweights)
  Cubature.solvecubature!(cub_omega, cub_omega_degree, mask_omega, tol=tol)

  mask_diage = SymCubatures.getInternalParamMask(cub_diage)
  append!(mask_diage, (cub_diage.numparams+1):(cub_diage.numparams+cub_diage.numweights))
  Cubature.solvecubature!(cub_diage, cub_diage_degree, mask_diage, tol=tol)

  return cub_omega, cub_diage, vtx
end

"""
### Cubature.equivalenceconstant{T}

Computes the equivalence constant for a given cubature; that is, it finds the
maximum eigenvalue for the matrix pk^T H pm, where H = diag(weights) and pk
denotes the orthogonal polynomial evaluated at the cubature points.

**Inputs**

* `cub`: symmetric cubature rule
* `vtx`: vertices of the right simplex
* `q`: maximum degree of polynomial for which the cubature is to be tested

**Outputs**

* `λmax`: maximum eigenvalue, which is the equivalence constant

"""
function equivalenceconstant(cub::TriSymCub{T}, vtx::Array{T,2}, q::Int) where {T}
  N = convert(Int, (q+1)*(q+2)/2) 
  P = zeros(T, (cub.numnodes,N) )
  x = SymCubatures.calcnodes(cub, vtx)
  ptr = 1
  for r = 0:q
    for j = 0:r
      i = r-j
      P[:,ptr] = OrthoPoly.proriolpoly(vec(x[1,:]), vec(x[2,:]), i, j)
      ptr += 1
    end
  end
  H = diagm(SymCubatures.calcweights(cub))
  A = P'*H*P
  return eigmax(0.5.*(A + A'))
end  

function equivalenceconstant(cub::TetSymCub{T}, vtx::Array{T,2}, q::Int) where {T}
  N = convert(Int, (q+1)*(q+2)*(q+3)/6 )
  P = zeros(T, (cub.numnodes,N) )
  x = SymCubatures.calcnodes(cub, vtx)
  ptr = 1
  for r = 0:q
    for k = 0:r
      for j = 0:r-k
        i = r-j-k
        P[:,ptr] = OrthoPoly.proriolpoly(vec(x[1,:]), vec(x[2,:]), vec(x[3,:]),
                                         i, j, k)
        ptr += 1
      end
    end
  end
  H = diagm(SymCubatures.calcweights(cub))
  A = P'*H*P
  return eigmax(0.5.*(A + A'))
end

end
