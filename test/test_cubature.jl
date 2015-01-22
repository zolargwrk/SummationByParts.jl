# Tests for Cubature module
using SummationByParts.Cubature

facts("Testing Cubature Module...") do

  context("Testing Cubature.cubatureresidual (TriSymCub method)") do
    # Check that a P1 vertices-only rule produces zero residual
    vtx = Float64[-1 -1; 1 -1; -1 1]
    cub = SymCubatures.TriSymCub{Float64}()
    x,y = SymCubatures.calcnodes(cub, vtx)
    SymCubatures.setweights!(cub, [2/3])
    F, dF = Cubature.cubatureresidual(cub, 1)
    @fact F => roughly(zeros(3), atol=1e-15)
  end

  context("Testing Cubature.cubatureresidual (TetSymCub method)") do
    # Check that a P1 vertices-only rule produces zero residual
    vtx = Float64[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
    cub = SymCubatures.TetSymCub{Float64}()
    x,y,z = SymCubatures.calcnodes(cub, vtx)
    SymCubatures.setweights!(cub, [1/3])
    F, dF = Cubature.cubatureresidual(cub, 1)
    @fact F => roughly(zeros(4), atol=1e-15)
  end

  context("Testing Cubature.solvecubature (TriSymCub method)") do
    # recover cubature for centroid only rule
    cub = SymCubatures.TriSymCub{Float64}(vertices=false, centroid=true)
    SymCubatures.setweights!(cub, [0.5])
    Cubature.solvecubature!(cub, 1, tol=1e-14)
    w = SymCubatures.calcweights(cub)
    @fact w => roughly([2], atol=1e-15)

    # recover P1 vertices-only rule
    cub = SymCubatures.TriSymCub{Float64}()
    SymCubatures.setweights!(cub, [0.5])
    Cubature.solvecubature!(cub, 1, tol=1e-14)
    w = SymCubatures.calcweights(cub)
    @fact w => roughly([2/3, 2/3, 2/3], atol=1e-15)

    # recover cubature for centroid+S21 (exact to degree 3)
    cub = SymCubatures.TriSymCub{Float64}(vertices=false, centroid=true, numS21=1)
    SymCubatures.setweights!(cub, [0.5, 0.5])
    SymCubatures.setparams!(cub, [0.25])
    Cubature.solvecubature!(cub, 3, tol=1e-14)
    w  = SymCubatures.calcweights(cub)
    @fact w => roughly([-18/16, 50/48, 50/48, 50/48], atol=1e-15)

    # create P2 element cubature with 1 bubble node at the centroid
    cub = SymCubatures.TriSymCub{Float64}(midedges=true, centroid=true)
    SymCubatures.setweights!(cub, [0.5, 0.5, 0.5])
    Cubature.solvecubature!(cub, 3, tol=1e-14)
    w  = SymCubatures.calcweights(cub)
    @fact w => roughly([1/10, 1/10, 1/10, 4/15, 4/15, 4/15, 18/20], atol=1e-15)

    # create P3 element cubature with 3 bubble nodes
    cub = SymCubatures.TriSymCub{Float64}(numedge=1, numS21=1)
    SymCubatures.setweights!(cub, [0.5, 0.5, 0.5])
    SymCubatures.setparams!(cub, [0.25, 0.25])
    Cubature.solvecubature!(cub, 5, tol=1e-14)
    w  = SymCubatures.calcweights(cub)
    @fact w => roughly([0.029745826049641155,0.029745826049641155,0.029745826049641155,0.097683362468102,0.097683362468102,0.097683362468102,0.097683362468102,0.097683362468102,0.097683362468102,0.44155411568082154,0.44155411568082154,0.44155411568082154], atol=1e-15)
  end

  context("Testing Cubature.solvecubature (TetSymCub method)") do
    # recover P1 vertices-only rule
    cub = SymCubatures.TetSymCub{Float64}()
    SymCubatures.setweights!(cub, [0.5])
    Cubature.solvecubature!(cub, 1, tol=1e-14)
    w = SymCubatures.calcweights(cub)
    @fact w => roughly([1/3, 1/3, 1/3, 1/3], atol=1e-15)

    # recover cubature for vertex+centroid (exact to degree 2)
    cub = SymCubatures.TetSymCub{Float64}(vertices=true, centroid=true)
    SymCubatures.setweights!(cub, [0.1 0.1])
    Cubature.solvecubature!(cub, 2, tol=1e-14)
    w = SymCubatures.calcweights(cub)
    @fact w => roughly([4/60, 4/60, 4/60, 4/60, 16/15], atol=1e-14)

    # recover cubature for S31+centroid (exact to degree 3)
    cub = SymCubatures.TetSymCub{Float64}(vertices=false, centroid=true,
                                          numS31=1)
    SymCubatures.setweights!(cub, [0.1 0.1])
    SymCubatures.setparams!(cub, [0.1])
    Cubature.solvecubature!(cub, 3, tol=1e-14)
    w = SymCubatures.calcweights(cub)
    #print("w = ",w,"\n")
    #print("params = ",cub.params,"\n")
    @fact w => roughly([-16/15, 3/5, 3/5, 3/5, 3/5], atol=1e-14)

    # create a P2 element with 1 bubble node
    cub = SymCubatures.TetSymCub{Float64}(midedges=true, centroid=true) #numS31=1)
    wuni = 0.1 #(4.0/3.0)/orbits.numnodes
    SymCubatures.setweights!(cub, [wuni wuni wuni]) #[0.02 0.05 0.23])
    Cubature.solvecubature!(cub, 3, tol=1e-14)
    w = SymCubatures.calcweights(cub)
    #print("w = ",w,"\n")
    #print("cub.params = ",cub.params,"\n")
    @fact w => roughly([2/90.*ones(4); 8/90.*ones(6); 32/45], atol=1e-14)

    # create a P3 element with 4 bubble nodes
    cub = SymCubatures.TetSymCub{Float64}(facecentroid=true,
                                          numedge=2, numS31=1)
    SymCubatures.setweights!(cub, [wuni wuni wuni wuni wuni])
    SymCubatures.setparams!(cub, [1/5 4/5 5/6])
    Cubature.solvecubature!(cub, 5, tol=1e-14)
    w = SymCubatures.calcweights(cub)
    #print("w = ",w,"\n")
    #print("orbits.params = ",cub.params,"\n")
    @fact w =>
    roughly([0.004421633248304779,0.004421633248304779,0.004421633248304779,0.004421633248304779,0.06935370366814568,0.06935370366814568,0.06935370366814568,0.06935370366814568,0.008837726716805293,0.008837726716805293,0.008837726716805293,0.008837726716805293,0.008837726716805293,0.008837726716805293,0.008837726716805293,0.008837726716805293,0.008837726716805293,0.008837726716805293,0.008837726716805293,0.008837726716805293,0.008837726716805204,0.008837726716805204,0.008837726716805204,0.008837726716805204,0.008837726716805204,0.008837726716805204,0.008837726716805204,0.008837726716805204,0.008837726716805204,0.008837726716805204,0.008837726716805204,0.008837726716805204,0.2065316361160515,0.2065316361160515,0.2065316361160515,0.2065316361160515], atol=1e-14)
  end

  context("Testing Cubature.tricubature") do
    # test using Float32, because this has not been done above
    w, x = tricubature(1, Float32)
    @fact w => roughly(Float32[2/3, 2/3, 2/3], atol=1e-7)
    w, x = tricubature(3, Float32)
    @fact w => roughly(Float32[1/10, 1/10, 1/10, 4/15, 4/15, 4/15, 18/20],
                       atol=1e-7)
    w, x = tricubature(5, Float32)
    @fact w =>
    roughly(Float32[0.029745826049641155,0.029745826049641155,0.029745826049641155,0.097683362468102,0.097683362468102,0.097683362468102,0.097683362468102,0.097683362468102,0.097683362468102,0.44155411568082154,0.44155411568082154,0.44155411568082154],
    atol=1e-7)
  end
  
    context("Testing Cubature.tetcubature") do
      # test using Float32, because this has not been done above
      w, x = tetcubature(1, Float32)
      @fact w => roughly(Float32[1/3, 1/3, 1/3, 1/3], atol=1e-7)   
      w, x = tetcubature(3, Float32)
      @fact w => roughly([Float32(2/90).*ones(Float32, 4);
                          Float32(8/90).*ones(Float32, 6);
                          Float32[32/45]], atol=1e-7)
      w, x = tetcubature(5, Float64)
      @fact w =>
      roughly([0.004421633248304779,0.004421633248304779,0.004421633248304779,0.004421633248304779,0.06935370366814568,0.06935370366814568,0.06935370366814568,0.06935370366814568,0.008837726716805293,0.008837726716805293,0.008837726716805293,0.008837726716805293,0.008837726716805293,0.008837726716805293,0.008837726716805293,0.008837726716805293,0.008837726716805293,0.008837726716805293,0.008837726716805293,0.008837726716805293,0.008837726716805204,0.008837726716805204,0.008837726716805204,0.008837726716805204,0.008837726716805204,0.008837726716805204,0.008837726716805204,0.008837726716805204,0.008837726716805204,0.008837726716805204,0.008837726716805204,0.008837726716805204,0.2065316361160515,0.2065316361160515,0.2065316361160515,0.2065316361160515], atol=1e-14)
    end

end