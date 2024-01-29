# if !("/home/z/zingg/workuzel/SummationByPartsPrivate/src/" in LOAD_PATH)
#     push!(LOAD_PATH, "/home/z/zingg/workuzel/SummationByPartsPrivate/src/") 
# end
# ENV["JULIA_NUM_THREADS"] = "40"
using SummationByParts
using SummationByParts.OrthoPoly
using SummationByParts.Cubature
using SummationByParts.SymCubatures
using SummationByParts.CubatureB
using SummationByParts.AsymCubatures
using SummationByParts.Optimizer
using LinearAlgebra

using Random
using PlotlyJS
using DelimitedFiles
using Latexify

# Define a custom show method for matrices in python format
function pyshow(m::Matrix)
    n1, n2 = size(m)
    print("[")
    for i in 1:n1
        print("[")
        for j in 1:n2
            print(m[i, j])
            if j < n2
                print(", ")
            end
        end
        print("]")
        if i < n1
            print(",")
        end
    end
    print("]")
end
function cppshow(m::Union{Matrix{T},Vector{T}}; lhs="x_cart = ") where {T}
    type_vec=false
    if typeof(m)==Vector{T}
        type_vec=true
        m = reshape(m,1, length(m))
    end
    n1, n2 = size(m)
    print("\n")
    print(lhs, "{")
    for i in 1:n1
        if !type_vec
            print("{")
        end
        for j in 1:n2
            print(m[i, j])
            if j < n2
                print(", ")
            end
        end
        if !type_vec
            print("}")
        end
        if i < n1
            print(",")
        end
    end
    print("};")
    print("\n")
end

function find_intersection(p1, p2, p3, p4)
    # Extract coordinates
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3
    x4, y4 = p4
    
    # Calculate slopes (m1 and m2) and intercepts (b1 and b2)
    m1 = (y2 - y1) / (x2 - x1)
    b1 = y1 - m1 * x1
    
    m2 = (y4 - y3) / (x4 - x3)
    b2 = y3 - m2 * x3
    
    # Calculate x-coordinate of the intersection
    x_intersect = (b2 - b1) / (m1 - m2)
    
    # Calculate y-coordinate using one of the line equations
    y_intersect = m1 * x_intersect + b1
    
    return (x_intersect, y_intersect)
end

# This file gathers functions that are used to search for minimum number of nodes in 
#SBP operators

function getOper(p)
    vertices = false
    faceoper =:Omega
    quad_degree=2*p-1
    oper = SummationByParts.getTetSBPDiagE(degree=p, faceopertype=faceoper, cubdegree=quad_degree)
    # oper = SummationByParts.getTetSBPOmega(degree=p)
    # oper = SummationByParts.getTetSBPGamma(degree=p)

    # oper = SummationByParts.getTriSBPDiagE(degree=p,vertices=vertices,quad_degree=quad_degree)

    # oper = SummationByParts.getTriSBPOmega(degree=p)
    # oper = SummationByParts.getTetSBPOmega(degree=p)
    # oper = SummationByParts.getTriSBPGamma(degree=p)

    x = SymCubatures.calcnodes(oper.cub, oper.vtx)
    #w = SymCubature.calcweights(oper.cub)
    err = truncErr(p, x, oper.w, oper.Q)
    quad_err = quadTruncErr(oper.cub, quad_degree)
    mindist = calcminnodedistance(oper, oper.vtx)
    println("min distance = ", round(mindist,digits=4))
    println("quadrature truncation err = ", quad_err)
    println("truncation err = ", err)

    #-------
    # Vfull = OrthoPoly.vandermonde_full(p,x)
    # S = svdvals(Vfull)
    # # println(prod(S))
    # println("|det(Vfull)| = ",abs(det(Vfull)))
    # println("min(S) = ", minimum(S), "     ","max(S) = ", maximum(S), "    ","cond(Vfull) = ", cond(Vfull) )
    #-------

    println("-------------------------------------------")
    # println(x)
    for i = 1:size(x)[1]
        println(join(x[i,:], ','))
    end
    println("-------------------------------------------")
    println(join(oper.w',','))
    println("-------------------------------------------")

    cub_facet, vtx_facet = SummationByParts.getTriCubatureForTetFaceDiagE(2*p, faceopertype=faceoper)
    numnodes_facet = cub_facet.numnodes
    println("numnodes_facet = ", numnodes_facet)
    xf = SymCubatures.calcnodes(cub_facet, vtx_facet)
    wf = SymCubatures.calcweights(cub_facet)
    println()
    println("===========================================")
    for i = 1:size(xf)[1]
        println(join(xf[i,:], ','))
    end
    println("-------------------------------------------")
    println(join(wf',','))
    println("-------------------------------------------")

    # Check if operator has volume nodes on the boundar or outside the domain
    checkInteriorNodeLocaton(oper.cub)

    return oper, x, oper.w'
end

function minTetCubatureDiagE(;q::Int=1, T=Float64)  

    tol=1e-14

    vtx = T[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
    
    # cub, vtx = getTetCubatureDiagE(q, vertices=false)
    # x = SymCubatures.calcnodes(cub, vtx)
    # w = SymCubatures.calcweights(cub)
    # return cub, vtx, x, w

    # vert_node = false
    # numS31_node = 0 # add 4
    # midedges_node = true
    # numS22_node = 0 # add 6
    # numfaceS21_node = 0
    # numedge_node = 0
    # numS211_node= 0  # add 12
    # numS21111_node =0 # add 24
    # facetcentroid_node = false
    # numfaceS111_node = 0
    # centroid_node = false
    

    vert_node = true
    numS31_node = 0# add 4
    midedges_node = true
    numS22_node = 0# add 6
    numfaceS21_node = 0
    numedge_node = 0
    numS211_node=0# add 12
    numfaceS111_node = 0
    facetcentroid_node = true
    numS21111_node = 0 # add 24
    centroid_node = true
    # vert_node = false
    # numS31_node = 1# add 4
    # midedges_node = true
    # numS22_node = 0 # add 6
    # numfaceS21_node = 0
    # numedge_node = 0
    # numS211_node=0 # add 12
    # numS21111_node =0 # add 24
    # facetcentroid_node = false
    # numfaceS111_node = 0
    # centroid_node = false

    cub = SymCubatures.TetSymCub{T}(vertices=vert_node, midedges=midedges_node, 
    facecentroid=facetcentroid_node, numedge=numedge_node, numfaceS21=numfaceS21_node, 
    numfaceS111=numfaceS111_node, centroid=centroid_node, numS31=numS31_node, 
    numS22=numS22_node, numS211=numS211_node, numS1111=numS21111_node)

    # q=4
    # SymCubatures.setparams!(cub, T[0.1, 0.22222222222222215])
    # SymCubatures.setparams!(cub, T[0.49999998838515475])
    # q=7
    # SymCubatures.setparams!(cub, T[0.1, 0.1, 0.1, 0.160991838340075, 0.8003678928805428, 0.6058255660767269, 0.2151836435697350])
    # q=8
    # SymCubatures.setparams!(cub, T[0.1, 0.1, 0.1, 0.160991838340075, 0.8003678928805428, 0.1, 0.1, 0.6058255660767269, 0.2151836435697350])
    # SymCubatures.setparams!(cub, T[0.1, 0.1, 0.22632295665280852, 0.34176809382044804, 0.08963891727960145, 0.2179650666798872, 0.644977982946902])
    # SymCubatures.setparams!(cub, T[0.1, 0.1, 0.1, 0.22632295665280555, 0.3417680938204483, 0.9103610827203965, 0.6449779829468978, 0.21796506667988766])
    # SymCubatures.setparams!(cub, T[0.1,0.1,0.1,0.915673676158322, 0.513718214523918, 0.11505536882281817, 0.6367019258463157, 0.867735418367292, 0.1, 0.1, 0.4420024375197799, 0.15638516725103452])
    # SymCubatures.setparams!(cub, T[0.22222222222222215])
    # SymCubatures.setparams!(cub, T[0.9465552688408668, 0.326390241974141, 0.7461816901961349, 0.160991838340075, 0.8003678928805428, 0.6058255660767269, 0.215183643569735])
    #q=9
    # SymCubatures.setparams!(cub,T[0.1, 0.1, 0.1, 0.915673676158322, 0.513718214523918, 0.11505536882281817, 0.6367019258463157, 0.867735418367292, 0.1, 0.1, 0.4420024375197799, 0.15638516725103452])
    # SymCubatures.setparams!(cub,T[0.1, 0.1, 0.1, 0.1, 0.40849548761177584, 0.1693510205976773, 0.895541492424858, 0.9192858150590572, 0.1, 0.1, 0.5199457984883094, 0.07662224108118755])
    # SymCubatures.setparams!(cub,T[0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1386668611165089, 0.9502911551207535, 0.32650013189680627, 0.9490995606077927, 0.6179556878200296, 0.1, 0.1, 0.1, 0.1,    0.6523549849744118, 0.15746672629638275, 0.6628851010581924, 0.41914124420319115, 0.05419324757802851, 0.3864445540346162])
    # SymCubatures.setparams!(cub, T[0.3771609693928902])
    # q=11
    # SymCubatures.setparams!(cub, T[0.1, 0.1, 0.1,
    #                                0.1, 0.1, 
    #                                0.9740380115784705, 0.09108541574493219, 0.8862300877546929, 
    #                                0.1431149128800494, 
    #                                0.1, 0.1, 0.1, 0.1, 
    #                                0.14954832600320392, 0.3552819365975969, 0.5734258260777191, 0.3234186242329474, 0.04978702493839072, 0.647438880058641])

    mask = SymCubatures.getInternalParamMask(cub)
    append!(mask, (cub.numparams+1):(cub.numparams+cub.numweights))
    # mask = 1:(cub.numparams+cub.numweights)    
    Cubature.solvecubature!(cub, q, mask, tol=tol, hist=true)
    
    return cub, vtx
end

function minTriCubatureDiagE(;q::Int=1, T=Float64)

    tol=1e-14

    vtx = T[-1 -1; 1 -1; -1 1]
    # vtx = T[0 0; 1 0; 0 1]

    # vert_node = true
    # midedges_node = false
    # numS21_node = 5
    # numedge_node = 5
    # numS111_node =8
    # centroid_node = false

    vert_node = false
    midedges_node = false
    numS21_node = 1
    numedge_node = 3
    numS111_node = 2
    centroid_node = false

    cub = SymCubatures.TriSymCub{T}(numedge=numedge_node, numS21=numS21_node, 
            numS111=numS111_node, vertices=vert_node, midedges=midedges_node, 
            centroid = centroid_node)
    
    # q=3
    # SymCubatures.setparams!(cub, T[0.7236067977499789])
    # q=13
    # SymCubatures.setparams!(cub, T[0.1, 0.1, 0.1, 0.1, 0.94987899770573, 0.8385931397553689, 0.6815587319130891, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
    # q=15
    # SymCubatures.setparams!(cub, T[0.1, 0.1, 0.1, 0.1, 0.1, 0.9597669540832294, 0.8693869325527526, 0.7389624749052223, 0.5826394788331934, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
    # q=5
    # SymCubatures.setparams!(cub, T[0.1, 0.8273268353539885])
    # q=7
    # SymCubatures.setparams!(cub, T[0.1, 0.8825276619647324, 0.642615758240322, 0.1, 0.1])
    # q=9
    # SymCubatures.setparams!(cub, T[0.1, 0.1, 0.9151119481392835, 0.7344243967353571, 0.1, 0.1])
    # q=10
    # SymCubatures.setparams!(cub, T[0.1, 0.1, 0.9151119481392835, 0.7344243967353571, 0.1, 0.1, 0.1, 0.1])
    # q = 11
    # SymCubatures.setparams!(cub, T[0.1, 0.1, 0.1, 0.1, 0.1, 0.9358700742548033, 0.7958500907165711, 0.6046496089512394, 0.1, 0.1])
    # q=12
    # SymCubatures.setparams!(cub, T[0.1, 0.1, 0.9358700742548033, 0.7958500907165711, 0.6046496089512394, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
    # q=17
    # SymCubatures.setparams!(cub, T[0.1, 0.1, 0.1, 0.1, 0.1, 0.9670007152040296, 0.8922417368315723, 0.7826176634981026, 0.6478790677934697, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
    # q = 19 
    # SymCubatures.setparams!(cub, T[0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.9724496361114411, 0.9096396608220033, 0.8164380765159303, 0.6997654704826745, 0.5682764664274638, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
    # q =20
    # SymCubatures.setparams!(cub, T[0.1, 0.1, 0.1, 0.1, 0.1, 0.9724496361114411, 0.9096396608220033, 0.8164380765159303, 0.6997654704826745, 0.5682764664274638, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])

    # q=6
    # SymCubatures.setparams!(cub, T[0.1, 0.1, 0.033765242898423975, 0.16939530676686776, 0.38069040695840156, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,])
    # SymCubatures.setparams!(cub, T[0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,])

    # q=20
    # SymCubatures.setparams!(cub, T[0.1, 0.1, 0.1, 0.010885670926971514, 0.05646870011595234, 0.13492399721297532, 0.2404519353965941, 0.36522842202382755, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
    # q=6
    # SymCubatures.setparams!(cub, T[0.1, 0.1, 0.06943184420297371, 0.33000947820757187, 0.1, 0.1])
    # q=9
    SymCubatures.setparams!(cub, T[0.4, 0.033765242898423975, 0.16939530676686776, 0.38069040695840156, 0.1, 0.1, 0.1, 0.1])
    # SymCubatures.setparams!(cub,T[0.1, 0.1, 0.1, 0.1, 0.019855071751231856, 0.10166676129318664, 0.2372337950418355, 0.4082826787521751, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,])
    
    mask = SymCubatures.getInternalParamMask(cub)
    append!(mask, (cub.numparams+1):(cub.numparams+cub.numweights))  
    # mask = 1:(cub.numparams+cub.numweights)  
    Cubature.solvecubature!(cub, q, mask, tol=tol, hist=true)

    return cub, vtx
end

function mingetLineSegSBPLobbato(;q::Int=1)

    cub, vtx = SummationByParts.Cubature.quadrature(q, internal=false)
   
    return cub, vtx
end

function mingetLineSegSBPLegendre(;q::Int=1)

    cub, vtx = SummationByParts.Cubature.quadrature(q, internal=true)
   
    return cub, vtx
end

function minTriCubatureGamma(;q::Int=1, T=Float64)

    tol=1e-14

    vtx = T[-1 -1; 1 -1; -1 1]

    vert_node = false
    midedges_node = false
    centroid_node = false
    numedge_node = 1
    numS21_node = 1
    numS111_node = 0

    cub = SymCubatures.TriSymCub{T}(numedge=numedge_node, numS21=numS21_node, 
            numS111=numS111_node, vertices=vert_node, midedges=midedges_node, 
            centroid = centroid_node)
            
    mask = SymCubatures.getInternalParamMask(cub)
    append!(mask, (cub.numparams):(cub.numparams+cub.numweights))
    # mask = 1:(cub.numparams+cub.numweights)
    Cubature.solvecubature!(cub, q, mask, tol=tol, hist=true)

    return cub, vtx
end

function minTriCubatureOmega(;q::Int=1, T=Float64)

    tol=1e-14

    vtx = T[-1 -1; 1 -1; -1 1]
    # vtx = T[0 0; 1 0; 0 1]

    vert_node = true 
    midedges_node = false
    centroid_node = true
    numedge_node = 0
    numS21_node = 5
    numS111_node = 2

    cub = SymCubatures.TriSymCub{T}(numedge=numedge_node, numS21=numS21_node, 
            numS111=numS111_node, vertices=vert_node, midedges=midedges_node, 
            centroid = centroid_node)
            
    mask = 1:(cub.numparams+cub.numweights)
    Cubature.solvecubature!(cub, q, mask, tol=tol, hist=true)

    return cub, vtx
end


# cub, vtx = minTriCubatureDiagE(q=9)
# cubtet, vtx= minTetCubatureDiagE(q=9)
# cub, vtx = minTriCubatureOmega(q=2)
# oper, x, w= getOper(2)
# cub, vert = mingetLineSegSBPLobbato(q=11)
# cub, vert = mingetLineSegSBPLegendre(q=7)
# x = SymCubatures.calcnodes(oper.cub, oper.vtx)

# p=5
# q=2p-0
# vtx = [1 -sqrt(3)/3 -sqrt(6)/6;
#        0 2*sqrt(3)/3 -sqrt(6)/6;
#        -1 -sqrt(3)/3 -sqrt(6)/6;
#        0 0 sqrt(6)/2]
# # xinit_param=[0.5659668188394116, 0.7807167565199028, 0.0846814812674472, 0.04046668289033655, 0.20685562107102665, 0.061511434675817794, 0.6927211311947881, 0.14342099002644818, 0.3786528536651681, 0.3596513879533776, 0.02013276646659172]
# # xinit_weight = [0.06845875962635535, 0.03996133054703136, 0.004908706029927262, 0.012435872193517524, 0.018315521918100697, 0.011842283605048317, 0.009933723304410471]
# # xinit = collect(Iterators.flatten([xinit_param,xinit_weight]))
# xinit= []
# cub, _ = SummationByParts.deriveTetCubatureOmega(q=q,
#                                                    vertices=false,
#                                                    numS31=5,
#                                                    midedges=false, 
#                                                    numS22=2,
#                                                    numfaceS21=0, 
#                                                    numedge=0, 
#                                                    numS211=3,
#                                                    numfaceS111=0, 
#                                                    facecentroid=false,
#                                                    numS1111=0,
#                                                    centroid=false,
#                                                    delta1=1e-3,
#                                                    delta2=1e-3,
#                                                    xinit=xinit,
#                                                    tol=5e-15,
#                                                    verbose=true)
# xyz = SymCubatures.calcnodes(cub,vtx)
# SummationByParts.plotly_tet_nodes(q=1, x=xyz, vtx=vtx)
# println("\n", cub.params,"\n")
# println(cub.weights,"\n")
# println(cub)
# V,_, _,_ = OrthoPoly.vandermonde(p, xyz[1,:], xyz[2,:], xyz[3,:])
# println("rank(V) = ", rank(V))
# checkInteriorNodeLocaton(cub)


# cub, vtx = SummationByParts.deriveTetCubatureGamma(q=5,
#                                                    vertices=true,
#                                                    numS31=1,
#                                                    midedges=false, 
#                                                    numS22=0,
#                                                    numfaceS21=0, 
#                                                    numedge=1, 
#                                                    numS211=0,
#                                                    numfaceS111=0, 
#                                                    facecentroid=true,
#                                                    numS1111=0,
#                                                    centroid=false,
#                                                    delta1=1e-3,
#                                                    delta2=1e-1,
#                                                    verbose=true)

# p = 2
# xinit = [0.8410998004251519, 0.3771609693928902, 0.0410802706918665, 0.00011786868691553593, 0.09988884100967022, 0.07775499410394296, 2.054710830645874e-5]
# cub, vtx = Cubature.getTriCubatureForTetFaceDiagE(2*p)
# xedge = collect(Iterators.flatten(collect(Iterators.flatten(vcat(cub.params, cub.weights)))))
# cub, vtx = SummationByParts.deriveTetCubatureDiagE(q=4,
#                                                    vertices=true,
#                                                    numS31=1,
#                                                    midedges=false, 
#                                                    numS22=0,
#                                                    numfaceS21=1, 
#                                                    numedge=1, 
#                                                    numS211=0,
#                                                    numfaceS111=0, 
#                                                    facecentroid=false,
#                                                    numS1111=0,
#                                                    centroid=false,
#                                                    xinit=[],
#                                                    xedge=xedge,
#                                                    delta1=1e-3,
#                                                    delta2=1e-2,
#                                                    verbose=true)

# qf = 8;
# cub_facet, vtx_facet = SummationByParts.deriveTriCubatureGamma(q=qf,
#                                                        vertices=true,
#                                                        midedges=true,
#                                                        numS21=2,
#                                                        numedge=1,
#                                                        numS111=0,
#                                                        centroid=true,
#                                                        xinit=[],
#                                                        delta1=1e-3,
#                                                        delta2=1e-1,
#                                                        verbose=true)

# xedge = cub_facet.params
# cub_tet, vtx_tet = SummationByParts.deriveTetCubatureDiagE(q=6,
#                                                    vertices=true,
#                                                    numS31=1,
#                                                    midedges=false, 
#                                                    numS22=1,
#                                                    numfaceS21=2, 
#                                                    numedge=1, 
#                                                    numS211=0,
#                                                    numfaceS111=0, 
#                                                    facecentroid=false,
#                                                    numS1111=0,
#                                                    centroid=true,
#                                                    xinit=[],
#                                                    xedge=xedge,
#                                                    delta1=1e-3,
#                                                    delta2=1e-1,
#                                                    verbose=false)
                                                   
# w, Q, E = SummationByParts.buildoperators(cub_tet, vtx_tet, 3, faceopertype=:DiagE, facecub=cub_facet, facevtx=vtx_facet)        


# sbp = TetSBP{Float64}(degree=2)


# # xinit = [0.2843222021131286, 0.06411074643388699, 
#          0.2962657715676411, 0.6436259905776708, 0.05923977897745962, 0.738293563655622, 0.3274034674743651, 0.05673533067987675,  

#          0.09191592720948945, 0.026705937626299157,
#          0.12780981279284817, 0.06836929632591884, 0.050595515414576735,
#          0.16348665829257197]
# xinit = []
# xinit = [0.3258262357481903, 0.057007000576775574, 0.29362301078786074, 1.0329852386556753, 
#          0.0673713973612208, 0.3066061103391226, 0.7267252339891408, 0.05861520900915901];
# cub, vtx = SummationByParts.deriveTriCubatureOmega(q=10,numS21=2,numS111=3,centroid=true,delta1=1e-3, delta2=1e-1, xinit = xinit, verbose=true,
#                                                     xinit_sym_group=[ "numS21","numS111","centroid"])

# xinit = [0.41469035132718185,0.29346955590904017, 0.02974582604964118,0.4415541156808217,0.09768336246810204]
# xinit = [0.014051428402751958, 0.3287022696410532, 0.117372049157595, 0.07570102262315924, 0.39953923368546157, 0.21669832279221568, 0.3411646771654866, 1.0108403093559752, 0.7489974554058978, 0.1069655164292194, 0.38625904523245186, 0.10265834844109815, 0.0006804547572204099, 0.0017378611770165666, 0.07750236796631678, 0.03817211660252702, 0.007120536305293825, 0.011767069907757862, 0.009369502376647165, 0.10172850185830666, 0.06858252249036582, 0.05540194870335351, 0.12190110864040862]
# xinit = [0.5749380303918797,0.11974220085793316,
#         0.33343014155598055,0.20437477140696825,
#         0.39432452613154606,0.06536359154212519,
#         0.38875869465672286,0.10195028749023816,
#         1.1491810826793598,0.09609353164480232,
#         0.6657786329998556,1.0308822535578346,
#         0.0013434826332230758,0.07754170837489897,
#         0.0392937103109862,0.08132263825474927,
#         0.009447719133224907,0.011139320053630379,
#         0.007018823640551441,0.055270427044833675,
#         0.06131670240670696,0.08938957126745721]
# xinit = []
# cub, vtx = SummationByParts.deriveTriCubatureGamma(q=13,vertices=true, midedges=false, numS21=3, numedge=3, numS111=3, centroid=false,
#                                                     delta1=1e-3, delta2=1e-4, 
#                                                     xinit = xinit, 
#                                                     verbose=true)

# xinit = []
# xedge = [0.8825276619647324, 0.6426157582403226]
# xedge = [0.9151119481392835, 0.7344243967353571]
# cub, vtx = SummationByParts.deriveTriCubatureDiagE(q=9,vertices=true, midedges=true, numS21=3, numedge=2, numS111=1, centroid=false,
#                                                     delta1=1e-2, delta2=1e-2, 
#                                                     xinit = xinit, 
#                                                     xedge = xedge,
#                                                     verbose=true)

# cub, vtx = SummationByParts.deriveTriCubatureOmega(q=12,centroid=false,numS21=5,numS111=3,delta1=1e-3, delta2=1e-1, verbose=true)
# cub, vtx = SummationByParts.deriveTriCubatureOmega(q=15,centroid=true,numS21=6,numS111=5,delta1=1e-3, delta2=1e-3, verbose=true)

##--------------------------
# use the following to check unisolvency; it should give value close to zero
# In 2D
# V, _,_ = OrthoPoly.vandermonde(oper.degree+10,x[1,:],x[2,:]); n=oper.numnodes; V = V[1:n,1:n]; norm(I(n)-V*inv(V))
# in 3D
# V, _,_ = OrthoPoly.vandermonde(oper.degree+10,x[1,:],x[2,:],x[3,:]); n=oper.numnodes; V = V[1:n,1:n]; norm(I(n)-V*inv(V))
##--------------------------

# println(cub)
# println(x[1,:])
# println(x[2,:])
# # println(x[3,:])
# println(w)

# oper = SummationByParts.getTetSBPDiagE(degree=3,vertices=false)
# xinit = 50*[0.45 0.3; 0.12 0.8;0.73 0.34;0.33 0.97;0.21 0.76]
# xinit=[]
# Optimizer.pso(Optimizer.rosenbrock, 100, xinit=xinit, np=20, xL=-5.0, xR=5.0,maxiter=10000000, tol=1e-15, save_iter=false,verbose=1)
# Optimizer.pso(Optimizer.rastrigin, 10, xinit=xinit, np=1000, xL=-5.0, xR=5.0,maxiter=10000000, tol=1e-15, save_iter=false,verbose=1)
# f =  Optimizer.rosenbrock([1.3,0.22])
# println(f)
# function testing(fun::Function, n::Int64)
#     print("working")
# end
# testing(Optimizer.rosenbrock,2)

# r = [0.030542254461691786815435989410616,
# 0.6554971842160468575855247763684,
# 0.036263912183347235540509245765861,
# 0.56234675429593172779618726053741,
# 0.7266761512591353167067609319929,
# 0.24015365562874879667987215725589,
#   0.47923880853381994882056460483,
# 0.59145130844516946577726912437356,
# 0.12159404179037625048920290282695,
# 0.61814494410729081685929031664273]
# s = [ 0.44321393070976644601444149884628, 
# 0.078046683005259298582245719444472, 
# 0.79974793251154241424671909044264,
#  0.88941338089262178545624237813172,
#  0.81700570546183981512911032041302,
# 0.026319253778551265909868561720941,
#  0.54824629826500270723954599816352,
#  0.42279009872944128822780385235092,
#  0.88426213454773316957613360500545,
#  0.66282072528512681053314281598432]
# # V, Vdx, Vdy, Vinteg = OrthoPoly.vandermonde_monomial(4, r, s)
# V, Vdx, Vdy = OrthoPoly.vandermonde(4, r, s)

# # x = [-0.8168   0.6337  -0.8168  -0.1081  -0.1081  -0.7838]
# # y = [-0.8168  -0.8168   0.6337  -0.7838  -0.1081  -0.1081]
# # V, Vdx, Vdy, Vinteg = OrthoPoly.vandermonde_monomial(2, x, y)

# # V =  [1 0.1 0.01; 1 0.5 0.25; 1 0.9 0.81]
# M = copy(Matrix(V'))
# M = Optimizer.preconditioner(M)

# x = Matrix(reshape(-1:0.5:1, :,1))
# V,Vdx,Hes = OrthoPoly.vandermonde_arnoldi(1,x,compute_grad=true)

# x = [-1.0;1.0;-1.0;-1.0]; y = [-1.0;-1.0;1.0;-1.0]; z = [-1.0;-1.0;-1.0;1.0]; 
# V,Vdx,Vdy,Vdz,Vinteg = OrthoPoly.vandermonde_monomial(1,x,y,z,compute_grad=true,compute_integ=true)

# x = [-1.0;1.0;-1.0]; y = [-1.0;-1.0;1.0];
# V,Vdx,Vdy,Hes = OrthoPoly.vandermonde_arnoldi(1,x,y,compute_grad=true)

# x = [-1.0;1.0;-1.0;-1.0]; y = [-1.0;-1.0;1.0;-1.0]; z = [-1.0;-1.0;-1.0;1.0];
# V,Vdx,Vdy,Vdz,Hes = OrthoPoly.vandermonde_arnoldi(1,x,y,z,compute_grad=true)

# x= [-1.0,1.0,-1.0,0.0,0.0,-1.0,-0.48890019595307077,-0.48890019595307077,-0.02219960809385846,-0.09107681185926941,-0.09107681185926941,-0.8178463762814612,-0.8000557112982302,-0.8000557112982302,0.6001114225964602,-0.8302238962785671,0.8302238962785671,0.8302238962785671,-0.8302238962785671,-1.0,-1.0,-0.4688487934707142,0.4688487934707142,0.4688487934707142,-0.4688487934707142,-1.0,-1.0,-0.42841695365174737,-0.8534906794713483,-0.8534906794713482,-0.42841695365174737,0.28190763312309564,0.2819076331230956];
# y = [-1.0,-1.0,1.0,-1.0,0.0,0.0,-0.02219960809385846,-0.48890019595307077,-0.48890019595307077,-0.8178463762814612,-0.09107681185926936,-0.09107681185926936,0.6001114225964603,-0.8000557112982302,-0.8000557112982302,-1.0,-1.0,-0.8302238962785671,0.8302238962785671,0.8302238962785671,-0.8302238962785671,-1.0,-1.0,-0.4688487934707142,0.4688487934707142,0.4688487934707142,-0.4688487934707142,0.28190763312309564,0.28190763312309564,-0.4284169536517474,-0.8534906794713483,-0.8534906794713483,-0.4284169536517474];
# # V,Vdx,Vdy,Hes = OrthoPoly.vandermonde_arnoldi(9,x,y,compute_grad=true)
# V,Vdx,Vdy,Vinteg = OrthoPoly.vandermonde_monomial(9,x,y,compute_grad=false,compute_integ=false)

# vtx = Float64[-1 -1; 1 -1; -1 1]
# condA = SummationByParts.computeConditionNumber(2,:DiagE,vtx)

# x=[-0.9369109855084978, -0.9369109855084978, 0.8738219710169957, -0.7507132548290896, -0.7507132548290896, 0.5014265096581791, -0.9468549501551831, -0.6896475489662155, -0.6896475489662155, -0.9468549501551831, 0.6365024991213986, 0.6365024991213986]
# y=[-0.9369109855084978, -0.9369109855084978, 0.8738219710169957, -0.7507132548290896, -0.7507132548290896, 0.5014265096581791, -0.9468549501551831, -0.6896475489662155, -0.6896475489662155, -0.9468549501551831, 0.6365024991213986, 0.6365024991213986]
# w=[0.10168981274041364, 0.10168981274041364, 0.10168981274041364, 0.23357255145275874, 0.23357255145275874, 0.23357255145275874, 0.16570215123674714, 0.16570215123674714, 0.16570215123674714, 0.16570215123674714, 0.16570215123674714, 0.16570215123674714]
# V,Vdx,Vdy,Vinteg = OrthoPoly.vandermonde_monomial(6,x,y,compute_grad=false,compute_integ=true)
# SummationByParts.getTetCubatureDiagE(4)
# SummationByParts.getTriCubatureForTetFaceDiagE(9)

# x = [-1.0, 1.0, -1.0, -0.44721359549995787, 0.44721359549995787, 0.44721359549995787, -0.44721359549995787, -1.0, -1.0, -0.3333333333333333]
# y = [-1.0, -1.0, 1.0, -1.0, -1.0, -0.44721359549995787, 0.44721359549995787, 0.44721359549995787, -0.44721359549995787, -0.3333333333333333]
# V, Vdx, Vdy = OrthoPoly.vandermonde(3,x,y,compute_grad=true)

#--------------------
# oper = SummationByParts.getTetSBPDiagE(degree=2, faceopertype=:DiagE, cubdegree=3)

tol=1e-14
mask = Int64[]
xinit = []
T = Float64
p=1;
q=2*p;
qf = q+mod(q,2)+1;
vtx_equilateral = T[0 0; 1 0; 1/2 sqrt(3/4)];
vtx_right = T[-1 -1; 1 -1; -1 1];
vtx_right_tet = T[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1];

# cub,vtx = SummationByParts.CubatureB.getTriCubatureDiagE(q, T, vertices=true);
# cub,vtx = SummationByParts.Cubature.getTriCubatureDiagE(q, T, vertices=true);
# cub,vtx = SummationByParts.CubatureB.getTetCubatureDiagE(q, T)
# cub,vtx = SummationByParts.Cubature.getTetCubatureDiagE(q, T)
# cub,vtx = SummationByParts.getTriCubatureGamma(q, T);

# cub_lgl,_ = SummationByParts.Cubature.quadrature(qf, internal=false); # get LGL nodes for the facet nodes
# xedge = cub_lgl.params

# xinit = [0.05314302102577261, 0.1666635950827463, 0.3145441089637701, 0.5094616454919219, 0.649816429368135, 0.9724496361114411, 0.9096396608220033, 0.8164380765159303, 0.6997654704826745, 0.5682764664274638, 0.05076655478490776, 0.1748658371401517, 0.15516865764729046, 0.3397740410570893, 0.047353516792609344, 0.35667773358319066, 0.31157712648071023, 0.5100287066627915, 0.15342665422017426, 0.5556000187641795, 0.04677187946920449, 0.5848582361562104, 0.510714684764098, 0.7074111605301531, 0.3158343317412534, 0.7279133295290301, 0.1562597309467402, 0.7965756807588819, 0.047762578570958385, 0.8421927194138071, 0.0002293568317083877, 0.00783875765951049, 0.021217180172630416, 0.03225655922925358, 0.04462993108992245, 0.003345027288644441, 0.0013377082718243282, 0.0022027174065093365, 0.0027659799779818895, 0.0032319726448620167, 0.0035682247898442686, 0.012954441477701537, 0.026586968785604882, 0.01633205973767796, 0.037746712865852614, 0.031089892534822434, 0.019144384397105894, 0.025645875054780776, 0.040947212218559376, 0.03389196423996988, 0.02112881279440114]
# xinit = [0.053139749801927004, 0.16668013760910033, 0.31456408264181895, 0.5092584781120019, 0.6575446980290367, 0.9724496361114411, 0.9096396608220033, 0.8164380765159303, 0.6997654704826745, 0.5682764664274638, 0.05077333145877637, 0.17485590647837276, 0.15517673991374759, 0.33979568808665384, 0.04735631662825568, 0.35666361674808805, 0.31149974556523374, 0.5100665247597754, 0.15338550495595896, 0.5556329842834971, 0.04675822774571211, 0.5848499558524128, 0.5105064158752478, 0.7065224622265588, 0.315690911435267, 0.7279773583684073, 0.15620091712376277, 0.7966103928710151, 0.04774670572255327, 0.8421945087727019, 0.00022935539418505812, 0.007837847151051467, 0.021219805667139453, 0.032260455291629936, 0.044527762572560574, 0.0032747908791620913, 0.0013376223002578604, 0.0022030273844742708, 0.002766117976534284, 0.0032309913489269974, 0.0035671030848210723, 0.012955494711603432, 0.026588559653793753, 0.016332897280745596, 0.03773686875624017, 0.03108243742956658, 0.019139493855889395, 0.025795173409900618, 0.04092458515061569, 0.03387561178369854, 0.021122340728400843]
# xinit = [0.05314083654939277, 0.16667506142371177, 0.3145561531413489, 0.5093061923989223, 0.6549792136742238, 0.9724496361114411, 0.9096396608220033, 0.8164380765159303, 0.6997654704826745, 0.5682764664274638, 0.05077121806372943, 0.1748591625081172, 0.1551726372074361, 0.33978941969431686, 0.04735486725061155, 0.3566683111752631, 0.31151775866380504, 0.510054294404, 0.1533956574428573, 0.555623755696761, 0.046761711552671954, 0.5848529031561664, 0.5105543189657, 0.706733020425517, 0.3157268703944357, 0.7279601416263844, 0.15621570450619568, 0.7966012113076019, 0.047750683865428564, 0.8421944177704404, 0.00022935586625155854, 0.00783814786086399, 0.021219026100221524, 0.03225936399260833, 0.04455220101558799, 0.003296332269910699, 0.0013376510523029352, 0.002202930009144862, 0.0027660390954023675, 0.0032312449049490756, 0.0035673833482377515, 0.01295517601165268, 0.026587921013393242, 0.016332453814833368, 0.03773921674052688, 0.03108421346808243, 0.019140697570692502, 0.025756853875834113, 0.040930581070070034, 0.033879844892502425, 0.021123912912986757]
# xinit = [0.053137030101403276, 0.1666915779213917, 0.31458753514287785, 0.5092021638875648, 0.6666668405345657, 0.9724496361114411, 0.9096396608220033, 0.8164380765159303, 0.6997654704826745, 0.5682764664274638, 0.05077820761400746, 0.17484789114373706, 0.15519113891527953, 0.3398086369535158, 0.04736145840713413, 0.3566518693716474, 0.3114787893327528, 0.5100954321200096, 0.15337109149682662, 0.5556508213036703, 0.04675279836594615, 0.5848419543641816, 0.5104427075211018, 0.7062588739799792, 0.3156389071224948, 0.7280067057908648, 0.15617981793355662, 0.7966256507246169, 0.04774113829849305, 0.8421930095751233, 0.00022935420935502203, 0.00783709998406271, 0.021221476161012897, 0.032262492434599534, 0.04449731844585891, 0.0032370623409345464, 0.0013375498099657509, 0.0022032543002081977, 0.002766411829736784, 0.0032305838474907924, 0.003566716200274423, 0.01295619922055755, 0.026590475456830104, 0.016334489192857883, 0.03773393588501787, 0.03108018187020854, 0.01913779655433308, 0.02584977633418476, 0.04091427525828663, 0.033868927348109375, 0.021120358437359668]
# xinit = [0.05294490986531675, 0.16773289449842865, 0.31622279360326333, 0.4879667756898061, 0.6513493837496125, 0.9724496361114411, 0.9096396608220033, 0.8164380765159303, 0.6997654704826745, 0.5682764664274638, 0.051195998723422086, 0.17426772855417563, 0.15559719406273465, 0.3410916700392985, 0.047462163492437585, 0.3558632464690134, 0.3015188695888322, 0.5149248126378921, 0.1485561447226546, 0.5580518026197324, 0.045282999147861325, 0.5845980081189597, 0.47766585164678715, 0.6759879691102543, 0.2984203296477869, 0.7360991066041834, 0.14866968526150143, 0.8002558630925941, 0.04562318536898432, 0.8428283926597013, 0.0002292509138976338, 0.007783575266946965, 0.021386288997735533, 0.03282340926195889, 0.03737509863710259, 0.01232514825527602, 0.0013325681607189903, 0.0022223561249509146, 0.002768972623204779, 0.0031274574625577323, 0.0034137825308797613, 0.013023905490374535, 0.02669382193873107, 0.016372911768879273, 0.036527379608411736, 0.0301333208225669, 0.01857476823754868, 0.032598392548470634, 0.03828529775483018, 0.032083926453061136, 0.02021308614168822]

# xinit = [0.5527864045000421, 0.7236067977499789]
# xinit = [0.5306627609684195, 0.20735501628561037, 
#          0.8825276619647324, 0.6426157582403226, 
#          0.17654792120316234, 0.6492809445301031]
# xinit = [0.7147684835193548, 0.23494467607053515, 
#          0.8825276619647324, 0.6426157582403226, 
#          1.05028684041011, 0.23494467607053515]
# cub = SymCubatures.TriSymCub{T}(vertices=true, 
#                                 midedges=false, 
#                                 numS21=1, 
#                                 numedge=1, 
#                                 numS111=0, 
#                                 centroid=false) 
# SymCubatures.setparams!(cub, xinit[1:cub.numparams])

# SymCubatures.setweights!(cub, [0.0007073734397208365, 0.06644745454101596, 0.1550028531199623, 0.028539132896818628, 0.033740426796698954, 0.15997493308946634])
# param_target=Int64[6,7,8]
# weight_target=collect(29:(29+23-1))
# mask=vcat(param_target,weight_target) 
# mask = Int64[]
# append!(mask, (cub.numparams+1):(cub.numparams+cub.numweights))

# cub, vtx = SummationByParts.deriveTriCubatureDiagE(q=q,
#                                                     vertices=true, 
#                                                     midedges=false, 
#                                                     numS21=1, 
#                                                     numedge=1, 
#                                                     numS111=0, 
#                                                     centroid=false,
#                                                     delta1=1e-2, 
#                                                     delta2=1e-2, 
#                                                     xinit = xinit, 
#                                                     xedge = xedge, 
#                                                     mask = mask,
#                                                     verbose=true)

# oper_lgl = getLineSegSBPLobbato(degree=p+1)
# Q1 = oper_lgl.Q[:,:,1]
# H1 = diagm(oper_lgl.w)
# D1 = inv(H1)*Q1[:,:,1]
# x1 = vec(SymCubatures.calcnodes(oper_lgl.cub, oper_lgl.vtx))
# x1_sort = sort(x1)
# perm = []
# for i=1:length(x1)
#     push!(perm, findall(a->x1[i] in a, x1_sort)[1])
# end

# x1 = SymCubatures.calcnodes(oper_lgl.cub, oper_lgl.vtx)
# y1 = 0 .*x1
# V1, Vdx1, Vdy1 = OrthoPoly.vandermonde(p, reshape(x1,(length(x1),1)), reshape(y1,(length(x1),1)))


# xy = SymCubatures.calcnodes(cub, vtx_right)
# x = (xy')[:,1]
# y = (xy')[:,2]
# z = (xy')[:,3]
# N = ones(cub.numnodes, 1)

# checkInteriorNodeLocaton(cub)
# SummationByParts.plot_tri_nodes(x=xy, vtx=vtx_right, q=q, n=cub.numnodes,write_title=true)

# save_dir="/Users/zelalemaregaworku/OneDriveUofT/UTIAS/RPI/efficient_sbp/figures/"
# SummationByParts.plot_tri_nodes(x=xy, vtx=vtx_equilateral,q=q,n=cub.numnodes,save_dir=save_dir,write_title=true)


# oper = SummationByParts.getTriSBPDiagE(degree=p,vertices=true,quad_degree=q)
# S = oper.Q - 0.5.*oper.E
# Dx = inv(diagm(oper.w))*oper.Q[:,:,1]
# # println("\n","Sx = ")
# # writedlm(stdout, round.(S[:,:,1],digits=4))
# nz_s= count(i->(abs(i)<1e-4),S[:,:,1]) #length(findall(abs.(S[:,:,1]) .< 1e-13))
# println("\n","nzero(Sx) = ", nz_s, "  : ", "size(Sx)=", size(S[:,:,1]), " = ", size(S[:,:,1])[1]^2, "  : r = ", nz_s/size(S[:,:,1])[1]^2)       

# ines = SummationByParts.sparse_stencil(cub, p);
# writedlm(stdout, lines)

# node_lines = SummationByParts.node_sparse_stencil(cub, p);
# node_lines[1,:,:]


# metric = SummationByParts.sparse_metric_terms(cub, vtx_right, p);

# D, Q, H, E, S, D2, Q2, H2, ddxi, ddeta, dxid, detad, jac, D1= SummationByParts.sparse_sbp_operator(cub, vtx_right, p);
# D, Q, H, E, S, D2= SummationByParts.sparse_sbp_operator(cub, vtx, p);


# errs_x, errs_y, rate_x, rate_y, h_vec, errs_val= SummationByParts.test_accuracy(cub, p) 
# println(round.(errs_x',sigdigits=4))
# println(round.(rate_x',sigdigits=4))
# println(round.(errs_y',sigdigits=4))
# println(round.(rate_y',sigdigits=4))


# ml, xi_x, eta_x = SummationByParts.accurate_metric_terms_newton(p)
# dxid, detad = SummationByParts.accurate_metric_terms_frobenius(p, D2)

# V, Vdx, Vdy = OrthoPoly.vandermonde_monomial(p, x, y)
# D, Q, H, E, S, D2, Q2, H2, ddxi, ddeta, dxid, detad, jac, D1= SummationByParts.sparse_sbp_operator(cub, vtx, p)
# xi_x = diagm(dxid[1,:])
# eta_x = diagm(detad[1,:])
# D_xi = D2[:,:,1]
# D_eta = D2[:,:,2]

# err = ((xi_x*D_xi + eta_x*D_eta)*V - Vdx)
# writedlm(stdout,round.(err,sigdigits=6))

# V, Vdx, Vdy = OrthoPoly.vandermonde_monomial(p, x, y)
# err1 = ((diagm(dxid[1,:])*D2[:,:,1] + diagm(detad[1,:])*D2[:,:,2])*V - Vdx)
# err2 = ((diagm(dxid[2,:])*D2[:,:,1] + diagm(detad[2,:])*D2[:,:,2])*V - Vdy)
# # err1 = norm(inv(diagm(oper.w))*oper.Q[:,:,1]*V - Vdx)
# # err2 = norm(inv(diagm(oper.w))*oper.Q[:,:,2]*V - Vdy)
# writedlm(stdout, round.(err1, sigdigits=6))
# println("-------")
# writedlm(stdout, round.(err2, sigdigits=4))

# err_x, err_y = SummationByParts.test_accuracy_D2(cub, p)

# oper = SummationByParts.getTriSBPDiagE(degree=p,quad_degree=q)
# w = oper.w
# Q = oper.Q
#-----------
# # xy = SymCubatures.calcnodes(cub, vtx_right)
# xy = SymCubatures.calcnodes(cub, vtx_right_tet)
# # xy = AsymCubatures.calcnodes(cub)
# x = (xy')[:,1]
# y = (xy')[:,2]
# z = (xy')[:,3]
# w, Q = SummationByParts.buildsparseoperators(cub, vtx, p)
# H = diagm(w)
# # V, Vdx, _,_ = OrthoPoly.vandermonde_monomial(p, x, y)
# V, Vdx, _,_ = OrthoPoly.vandermonde_monomial(p, x, y, z)
# Dx = inv(H)*Q[:,:,1]
# Sx = Q[:,:,1] - 0.5.*(Q[:,:,1]+Q[:,:,1]')
# println(norm(Dx*V - Vdx))
# # writedlm(stdout, round.(Dx*V - Vdx,sigdigits=4))
# nz_s= count(i->(abs(i)<1e-10),Sx[:,:,1]) 
# println("\n","p = ", p, "  :  n = ", cub.numnodes, "  :  nzero(Sx) = ", nz_s, "  : ", "size(Sx)=", size(Sx), " = ", size(Sx)[1]^2, "  : r = ", round(nz_s/size(Sx)[1]^2,digits=4))    

#-----------
# check if OMP can solve for S for a given sparsity
# face = TriFace{T}(p, cub, vtx, vertices=true)
# Q = zeros(T, (cub.numnodes,cub.numnodes,2) )
# SummationByParts.boundaryoperator!(face, 1, view(Q,:,:,1))
# SummationByParts.boundaryoperator!(face, 2, view(Q,:,:,2))
# Q .*= 0.5 #scale!(Q, 0.5)
# A, bx, by = SummationByParts.accuracyconstraints(cub, vtx, p, Q)
# idx = SummationByParts.independent_rows(A)
# A = A[idx,:]
# M = SummationByParts.compute_mutual_coherence(A)
# s = (cub.numnodes^2-10)
# println(M)
# println(1/(2*s-1))
#-------------

#-------------
# N = 2*(q-1)+ 0
# p=1;
# q=2*p;
# N = 22#q+mod(q,2)+1 +0
# cub_1d, vtx_1d= SummationByParts.Cubature.quadratureUniform(2*p, N, T); 
# # cub_1d, vtx_1d = SummationByParts.Cubature.quadrature(q+1, internal=false);
# xedge = SymCubatures.calcnodes(cub_1d, vtx_1d)

# # cub = AsymCubatures.TriAsymCub{T}(numedgenodes=N)
# # cub = AsymCubatures.TriAsymCub{T}(numedgenodes=length(xedge))
# cub = AsymCubatures.TetAsymCub{T}(numedgenodes=length(xedge))
# # # # cub = AsymCubatures.TetAsymCub{T}(numedgenodes=length(xedge))
# AsymCubatures.setparams!(cub, xedge)

# # # xy = AsymCubatures.calcnodes(cub)
# # # SummationByParts.plot_tri_nodes(x=xy, vtx=vtx_right, q=q, n=cub.numnodes,write_title=true)
# # # SummationByParts.plot_tet_nodes(x=xy, q=q, n=cub.numnodes)
# # cub_g,vtx_g = SummationByParts.Cubature.getTriCubatureGregory(q-1,N)
# cub_g,vtx_g = SummationByParts.Cubature.getTetCubatureGregory(q-3,N)
# xinit = zeros((cub.numparams+cub.numweights,1))
# xinit[1:cub.numparams] = cub.params
# xinit[(cub.numparams+1):(cub.numparams+cub.numweights)] = cub_g.weights

# append!(mask, (cub.numparams+1):(cub.numparams+cub.numweights))
# # # append!(mask, 1:(cub.numparams+cub.numweights))
# CubatureB.solvecubature!(cub, q, mask, tol=tol, delta1=1e-2, delta2=1e-2, xinit=xinit, verbose=true, hist=true)
# println("min(w) = ", minimum(cub.weights),"\n")
# println(cub.edgeparams,"\n")
# println(cub.weights,"\n")

# xy = AsymCubatures.calcnodes(cub)
# SummationByParts.plot_tri_nodes(x=xy, vtx=vtx_right, q=q, n=cub.numnodes,save_dir=save_dir,write_title=true,label_nodes=false)
# SummationByParts.plot_tet_nodes(x=xy, q=q, n=cub.numnodes)

# cub_oper = CubatureB.getTriCubatureDiagE(q, T)

# cub_1d, vtx_1d= SummationByParts.Cubature.quadratureUniform(q, q+1, T); 
# x = sort(vec(SymCubatures.calcnodes(cub_1d, vtx_1d)))
# x2 = zeros(T,size(x))
# x2[1]=-1.0
# x2[end] = 2*sqrt(2)-1
# for i = 2:length(x)-1
#     x2[i] = x2[1]+ (x2[end]+1)-sqrt((x[i]-1)^2+(x[length(x)-(i-1)]-(-1))^2)
# end
# x2 = (x2.+1)./(x2[end]-x2[1]).*2 .- 1
# println(norm(x2 .- x))

# # oper = SummationByParts.getTriSBPDiagE(degree=p,quad_degree=q,asym_cub=true)
# # w = oper.w
# # Q = oper.Q
# # w, Q = SummationByParts.buildoperators(cub, vtx_right, p)
# w, Q = SummationByParts.buildsparseoperators(cub, vtx, p)
# H = diagm(w)
# # xy = AsymCubatures.calcnodes(oper.cub)
# # xy = SymCubatures.calcnodes(cub,vtx)
# xy = AsymCubatures.calcnodes(cub)
# x = (xy')[:,1]
# y = (xy')[:,2]
# V, Vdx, _,_ = OrthoPoly.vandermonde_monomial(p, x, y)
# # V, Vdx, _,_ = OrthoPoly.vandermonde_monomial(p, x, y, z)
# Qx = Q[:,:,1]
# Dx = inv(H)*Q[:,:,1]
# Sx = Q[:,:,1] - 0.5.*(Q[:,:,1]+Q[:,:,1]')
# println("accuracy: ", norm(Dx*V - Vdx))
# # writedlm(stdout, round.(Dx*V - Vdx,sigdigits=6))
# # println("sbp_property: ", norm(Q[:,:,1]+Q[:,:,1]'-oper.E[:,:,1]))
# nz_s= count(i->(abs(i)<1e-10),Qx[:,:,1]) 
# # println("\n","p = ", p, "  :  n = ", oper.cub.numnodes, "  :  nzero(Sx) = ", nz_s, "  : ", "size(Sx)=", size(Sx), " = ", size(Sx)[1]^2, "  : r = ", round(nz_s/size(Sx)[1]^2,digits=4))    
# println("\n","p = ", p, "  :  n = ", cub.numnodes, "  :  nzero(Qx) = ", nz_s, "  : ", "size(Qx)=", size(Qx), " = ", size(Qx)[1]^2, "  : r = ", round(nz_s/size(Qx)[1]^2,digits=4))    
# # errs_x, errs_y, rate_x, rate_y, h_vec, errs_val= SummationByParts.test_accuracy(oper, p) 
# # println(errs_x)
# # println(rate_x)

#-----------
# # Gregory quadrature
# q=9
# N = 2*(q-1)+10
# quad,vtx = SummationByParts.Cubature.getLineCubatureGregory(q,N)
# xgreg = SymCubatures.calcnodes(quad,vtx)
# perm = sortperm(vec(xgreg))
# wgreg = SymCubatures.calcweights(quad)
# # println(wgreg[perm]')
# # sum(wgreg)
# # V, _, Vinteg=OrthoPoly.vandermonde_monomial(q, (xgreg.+1.0)./2.0, compute_grad=false, compute_integ=true);
# # println((0.5 .*wgreg'*V)' - Vinteg)
# errs=[]
# rates=[]
# exact = exp(1.0)-1.0;
# # nn = [11,22,44,88,176]
# nn = [N,2*N,4*N,8*N]
# for i=1:length(nn)
#     n = nn[i] #10*10^i
#     quad_,vtx_ = SummationByParts.Cubature.getLineCubatureGregory(q,n)
#     xgreg_ = SymCubatures.calcnodes(quad_,vtx_)
#     wgreg_ = SymCubatures.calcweights(quad_)
#     V_, _, Vinteg_=OrthoPoly.vandermonde_monomial(q, (xgreg_.+1.0)./2.0, compute_grad=false, compute_integ=true);
#     integ_err = 0.5 .*(wgreg_'*V_)'- Vinteg_
#     # integ_err = 0.5 .*(wgreg_'*exp.((xgreg_.+1.0)./2.0)') .-exact
#     append!(errs, norm(integ_err))
#     # println(integ_err)
#     # V_, _=OrthoPoly.vandermonde(q, xgreg_, compute_grad=false);
#     # integ_err = (wgreg_'*V_)'
#     # integ_err[1] -= 2.0/sqrt(2.0)
#     # println(integ_err)
#     # append!(errs, norm(integ_err))
#     if length(errs)>1
#         append!(rates, log10(errs[i]/errs[i-1])/log10((1/nn[i])/(1/nn[i-1])))
#     end
# end
# println("errs: ", errs)
# println("rates: ", rates)

#-----------
# q=2
# N = 2*(q-1)+4
# quad,vtx = SummationByParts.Cubature.getTriCubatureGregory(q,N)
# xgreg = AsymCubatures.calcnodes(quad)
# wgreg = AsymCubatures.calcweights(quad)

# errs=[]
# rates=[]
# # exact = 0.5*(exp(2)+exp(-2)-2)
# # exact = 0.0
# # exact = 1.0
# nn = [N,2*N,4*N,8*N]
# for i=1:length(nn)
#     n = nn[i] #10*10^i
#     cub_,vtx_ = SummationByParts.Cubature.getTriCubatureGregory(q,n)
#     xgreg_ = AsymCubatures.calcnodes(cub_)
#     wgreg_ = AsymCubatures.calcweights(cub_)
#     # V_, _,_, Vinteg_=OrthoPoly.vandermonde_monomial(q, (xgreg_[1,:].+1.0)./2.0, (xgreg_[2,:].+1.0)./2.0, compute_grad=false, compute_integ=true);
#     # integ_err = 0.25 .*(wgreg_'*V_)'- Vinteg_
#     V_,_,_=OrthoPoly.vandermonde(q, xgreg_[1,:], xgreg_[2,:], compute_grad=false);
#     integ_err = (wgreg_'*V_)'
#     integ_err[1] -= 2.0/sqrt(2.0)
#     # integ_err = (wgreg_'*exp.(xgreg_[1,:].+xgreg_[2,:])) .-exact
#     # integ_err = (wgreg_'*sin.(2*pi*(xgreg_[1,:].+xgreg_[2,:]))) .-exact
#     # println(integ_err)
#     # xgreg_ = (xgreg_.+1)./2
#     # integ_err = (0.25 .*wgreg_'*exp.(xgreg_[1,:]+xgreg_[2,:])) .-exact
#     # integ_err = (0.25 .*wgreg_'*sin.(2*pi*(xgreg_[1,:].+xgreg_[2,:]))) .-exact
#     append!(errs, norm(integ_err))
#     # println(integ_err)
#     # V_, _=OrthoPoly.vandermonde(q, xgreg_, compute_grad=false);
#     # integ_err = (wgreg_'*V_)'
#     # integ_err[1] -= 2.0/sqrt(2.0)
#     # println(integ_err)
#     # append!(errs, norm(integ_err))
#     if length(errs)>1
#         append!(rates, log10(errs[i]/errs[i-1])/log10((1/nn[i])/(1/nn[i-1])))
#     end
# end
# println("errs: ", errs)
# println("rates: ", rates)
#-----------
# p=2
# q=2*p
# N = 2*(q-1)+1
# cub,vtx = SummationByParts.Cubature.getTriCubatureGregory(q,N)
# # cub,vtx = SummationByParts.CubatureB.getTriCubatureDiagE(q, T, vertices=true);
# # xgreg = AsymCubatures.calcnodes(cub)
# # wgreg = AsymCubatures.calcweights(cub)

# errs=[]
# rates=[]

# nn = [N,2*N,4*N,8*N]
# for i=1:length(nn)
#     n = nn[i] #10*10^i
#     cub_,vtx_ = SummationByParts.Cubature.getTetCubatureGregory(q,n)
#     xgreg_ = AsymCubatures.calcnodes(cub_)
#     wgreg_ = AsymCubatures.calcweights(cub_)
#     # V_, _,_, Vinteg_=OrthoPoly.vandermonde_monomial(q, (xgreg_[1,:].+1.0)./2.0, (xgreg_[2,:].+1.0)./2.0, compute_grad=false, compute_integ=true);
#     # integ_err = 0.25 .*(wgreg_'*V_)'- Vinteg_
#     V_,_,_,_=OrthoPoly.vandermonde(q, xgreg_[1,:], xgreg_[2,:], xgreg_[3,:], compute_grad=false)
#     integ_err = (wgreg_'*V_)'
#     integ_err[1] -= 2.0/sqrt(3.0)
#     # integ_err = (wgreg_'*exp.(xgreg_[1,:].+xgreg_[2,:])) .-exact
#     # integ_err = (wgreg_'*sin.(2*pi*(xgreg_[1,:].+xgreg_[2,:]))) .-exact
#     # println(integ_err)
#     # xgreg_ = (xgreg_.+1)./2
#     # integ_err = (0.25 .*wgreg_'*exp.(xgreg_[1,:]+xgreg_[2,:])) .-exact
#     # integ_err = (0.25 .*wgreg_'*sin.(2*pi*(xgreg_[1,:].+xgreg_[2,:]))) .-exact
#     append!(errs, norm(integ_err))
#     # println(integ_err)
#     # V_, _=OrthoPoly.vandermonde(q, xgreg_, compute_grad=false);
#     # integ_err = (wgreg_'*V_)'
#     # integ_err[1] -= 2.0/sqrt(2.0)
#     # println(integ_err)
#     # append!(errs, norm(integ_err))
#     if length(errs)>1
#         append!(rates, log10(errs[i]/errs[i-1])/log10((1/nn[i])/(1/nn[i-1])))
#     end
# end
# println("errs: ", errs)
# println("rates: ", rates)

# nn = [N,N+4,N+8,N+12]
# ndof = [(i+1)*i/2 for i in nn]
# for i=1:length(nn)
#     n = nn[i] #10*10^i
#     cub_,vtx_ = SummationByParts.Cubature.getTriCubatureGregory(q,n)
#     w_, Q_ = SummationByParts.buildoperators(cub_, vtx_, p)
#     H_ = diagm(w_)
#     xgreg_ = AsymCubatures.calcnodes(cub_)
#     x_ = (xgreg_')[:,1]
#     y_ = (xgreg_')[:,2]
#     V_, Vdx_, _,_ = OrthoPoly.vandermonde_monomial(p, x_, y_)
#     Qx_ = Q_[:,:,1]
#     Dx_ = inv(H_)*Q_[:,:,1]
#     err_ = Dx_*V_ - Vdx_
#     append!(errs,sum(err_'*H_*err_))
    
#     if length(errs)>1
#         append!(rates, log10(errs[i]/errs[i-1])/log10((1/ndof[i])/(1/ndof[i-1])))
#     end
#     println(errs)
#     println(rates)
# end
# println("errs: ", errs)
# println("rates: ", rates)

# # w, Q = SummationByParts.buildsparseoperators(cub, vtx, p)
# w, Q = SummationByParts.buildoperators(cub, vtx, p)
# H = diagm(w)
# # xy = AsymCubatures.calcnodes(oper.cub)
# xy = SymCubatures.calcnodes(cub,vtx_equilateral)
# xy = AsymCubatures.calcnodes(cub)
# x = (xy')[:,1]
# y = (xy')[:,2]
# V, Vdx, _,_ = OrthoPoly.vandermonde_monomial(p, x, y)
# # V, Vdx, _,_ = OrthoPoly.vandermonde_monomial(p, x, y, z)
# Qx = Q[:,:,1]
# Dx = inv(H)*Q[:,:,1]
# Sx = Q[:,:,1] - 0.5.*(Q[:,:,1]+Q[:,:,1]')
# println("accuracy: ", norm(Dx*V - Vdx))
# # writedlm(stdout, round.(Dx*V - Vdx,sigdigits=6))
# # println("sbp_property: ", norm(Q[:,:,1]+Q[:,:,1]'-oper.E[:,:,1]))
# nz_s= count(i->(abs(i)<1e-10),Qx[:,:,1]) 
# # println("\n","p = ", p, "  :  n = ", oper.cub.numnodes, "  :  nzero(Sx) = ", nz_s, "  : ", "size(Sx)=", size(Sx), " = ", size(Sx)[1]^2, "  : r = ", round(nz_s/size(Sx)[1]^2,digits=4))    
# println("\n","p = ", p, "  :  n = ", cub.numnodes, "  :  nzero(Qx) = ", nz_s, "  : ", "size(Qx)=", size(Qx), " = ", size(Qx)[1]^2, "  : r = ", round(nz_s/size(Qx)[1]^2,digits=4))    
# errs_x, rate_x, h_vec, errs_val= SummationByParts.test_accuracy(Qx, w, x, y, refine=6) 
# println(errs_x)
# println(rate_x)

# SummationByParts.plot_tri_nodes(x=xy, vtx=vtx_equilateral, q=q, n=cub.numnodes,write_title=true,label_nodes=false)
# SummationByParts.plot_tet_nodes(x=xy, q=q, n=cub.numnodes)

#--------------
p =2
q = 6*p
vtx = T[-1 -1/sqrt(3); 1 -1/sqrt(3); 0 2/sqrt(3)]
# cub = SummationByParts.init_tri_nodes(p, q)
# # # cub = SummationByParts.init_tri_uni_nodes(p, q)
# xy = SymCubatures.calcnodes(cub,vtx)
# SummationByParts.plot_tri_nodes(x=xy, vtx=vtx, q=q, n=cub.numnodes,write_title=true,label_nodes=false)

# coef = SummationByParts.get_cut_rect_to_equi_tri_mapping(p)
# cut_rect = SummationByParts.generate_cut_rect(p)
# xy_tri = SummationByParts.map_cut_rec_to_equi_tri(p)
# xlgl_2d = SummationByParts.get_lgl_nodes_2d(p, Array(T[-1,-1/sqrt(3)]), Array(T[1, -1/sqrt(3)]), half=false)
# parab = SummationByParts.get_parabola_tri(p)
# p1=coef[:,2];
# p2=coef[:,3];
# points = SummationByParts.get_interior_points_tri(p)

# cub, _= SummationByParts.Cubature.getTriCubatureDiagE(q, vertices=true)
# xy = SymCubatures.calcnodes(cub,vtx)
# SummationByParts.plot_tri_nodes(x=xy, vtx=vtx, q=q, n=cub.numnodes,write_title=true,label_nodes=false)
# mindist = calcminnodedistance(cub, vtx)
# println("\n","mindist = ", mindist)
# println("\n","minweight = ", minimum(SummationByParts.SymCubatures.calcweights(cub)))
# println("\n","midedges = ", cub.midedges,",")
# println("numS21 = ", cub.numS21,",")
# println("numedge = ", cub.numedge,",")
# println("numS111 = ", cub.numS111,",")
# println("centroid = ", cub.centroid)

# cub, _= SummationByParts.Cubature.getTriCubatureDiagE(q)
# xy = SymCubatures.calcnodes(cub,vtx)
# SummationByParts.plot_tri_nodes(x=xy, vtx=vtx, q=q, n=cub.numnodes,write_title=true,label_nodes=false)


# # check accuracy
# # vtx = T[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
# vtx_01 = T[0 0 ; 1 0; 0 1]
# xy = SymCubatures.calcnodes(cub,vtx_01)
# w = SymCubatures.calcweights(cub).*(1/4)
# V, _, _, Vinteg=OrthoPoly.vandermonde_monomial(q, xy[1,:], xy[2,:], compute_grad=false, compute_integ=true);
# integ_err = norm((w'*V)' - Vinteg)
# println("\n", "integ error = ", integ_err,"\n")
# checkInteriorNodeLocaton(cub)
# mindist = calcminnodedistance(cub, vtx)
# println("\n","mindist = ", mindist)
# println("\n","minweight = ", minimum(w))

# cub, _= SummationByParts.Cubature.getTriCubatureDiagE(q)
# # w, Q = SummationByParts.buildoperators(cub, vtx, p)
# w, Q = SummationByParts.buildsparseoperators(cub, vtx, p)
# H = diagm(w)
# # xy = AsymCubatures.calcnodes(oper.cub)
# xy = SymCubatures.calcnodes(cub,vtx_right)
# x = (xy')[:,1]
# y = (xy')[:,2]
# V, Vdx, _,_ = OrthoPoly.vandermonde_monomial(p, x, y)
# # V, Vdx, _,_ = OrthoPoly.vandermonde_monomial(p, x, y, z)
# Qx = Q[:,:,1]
# Dx = inv(H)*Q[:,:,1]
# Sx = Q[:,:,1] - 0.5.*(Q[:,:,1]+Q[:,:,1]')
# println("accuracy: ", norm(Dx*V - Vdx))
# # writedlm(stdout, round.(Dx*V - Vdx,sigdigits=6))
# # println("sbp_property: ", norm(Q[:,:,1]+Q[:,:,1]'-oper.E[:,:,1]))
# nz_s= count(i->(abs(i)<1e-10),Qx[:,:,1]) 
# # println("\n","p = ", p, "  :  n = ", oper.cub.numnodes, "  :  nzero(Sx) = ", nz_s, "  : ", "size(Sx)=", size(Sx), " = ", size(Sx)[1]^2, "  : r = ", round(nz_s/size(Sx)[1]^2,digits=4))    
# println("\n","p = ", p, "  :  n = ", cub.numnodes, "  :  nzero(Qx) = ", nz_s, "  : ", "size(Qx)=", size(Qx), " = ", size(Qx)[1]^2, "  : r = ", round(nz_s/size(Qx)[1]^2,digits=4))    
# errs_x, rate_x, h_vec, errs_val= SummationByParts.test_accuracy(Qx, w, x, y, refine=6) 
# println(errs_x)
# println(rate_x)

#-------------------------------
#3d quads
# p=4
# q=2*p
# vtx = [1 -sqrt(3)/3 -sqrt(6)/6;
#        0 2*sqrt(3)/3 -sqrt(6)/6;
#        -1 -sqrt(3)/3 -sqrt(6)/6;
#        0 0 sqrt(6)/2]

# cub, _= SummationByParts.Cubature.getTetCubatureDiagE(q)
# xyz = SymCubatures.calcnodes(cub,vtx)
# SummationByParts.plotly_tet_nodes(q=1, x=xyz, vtx=vtx)
# cub_lgl, vtx_lgl =  SummationByParts.Cubature.quadrature(2*(p+2)-3, internal=false)
# vtx_lgl[1]=0; vtx_lgl[2]=1;
# xlgl_1d =sort(vec(SymCubatures.calcnodes(cub_lgl,vtx_lgl)))

# xlgl_3d, wlgl_3d = SummationByParts.get_lgl_nodes_3d(p,[1,0,0],[2,0,0],half=true)
# nodes,weights = SummationByParts.get_interior_points_tet(p)
# cub = SummationByParts.init_tet_nodes(p,q)
# xyz = SymCubatures.calcnodes(cub,vtx)
# SummationByParts.plotly_tet_nodes(q=q, x=xyz, vtx=vtx)

# println(cub)

# @time begin
#     cub = SummationByParts.init_tet_nodes(p,q)
# end

# @profview SummationByParts.init_tet_nodes(p,q)


# cub, _= SummationByParts.Cubature.getTetCubatureDiagELGL(q)
# xyz = SymCubatures.calcnodes(cub,vtx)
# SummationByParts.plotly_tet_nodes(q=1, x=xyz, vtx=vtx)

# # check accuracy
# # vtx = T[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
# vtx_01 = T[0 0 0; 1 0 0; 0 1 0; 0 0 1]
# xyz = SymCubatures.calcnodes(cub,vtx_01)
# w = SymCubatures.calcweights(cub).*(1/8)
# V, _, _, _, Vinteg=OrthoPoly.vandermonde_monomial(q, xyz[1,:], xyz[2,:], xyz[3,:], compute_grad=false, compute_integ=true);
# integ_err = norm((w'*V)' - Vinteg)
# println("integ error = ", integ_err)
# checkInteriorNodeLocaton(cub)
# mindist = calcminnodedistance(cub, vtx)
# println("\n","mindist = ", mindist)
# println("\n","minweight = ", minimum(w))


# w, Q = SummationByParts.buildsparseoperators(cub, vtx_right_tet, p)
# H = diagm(w)
# # xy = AsymCubatures.calcnodes(oper.cub)
# xyz = SymCubatures.calcnodes(cub,vtx_right_tet)
# x = (xyz')[:,1]
# y = (xyz')[:,2]
# z = (xyz')[:,3]
# V, Vdx, _,_,_= OrthoPoly.vandermonde_monomial(p, x, y, z)
# # V, Vdx, _,_ = OrthoPoly.vandermonde_monomial(p, x, y, z)
# Qx = Q[:,:,1]
# Dx = inv(H)*Q[:,:,1]
# Sx = Q[:,:,1] - 0.5.*(Q[:,:,1]+Q[:,:,1]')
# println("accuracy: ", norm(Dx*V - Vdx))
# # writedlm(stdout, round.(Dx*V - Vdx,sigdigits=6))
# # println("sbp_property: ", norm(Q[:,:,1]+Q[:,:,1]'-oper.E[:,:,1]))
# nz_s= count(i->(abs(i)<1e-10),Qx[:,:,1]) 
# # println("\n","p = ", p, "  :  n = ", oper.cub.numnodes, "  :  nzero(Sx) = ", nz_s, "  : ", "size(Sx)=", size(Sx), " = ", size(Sx)[1]^2, "  : r = ", round(nz_s/size(Sx)[1]^2,digits=4))    
# println("\n","p = ", p, "  :  n = ", cub.numnodes, "  :  nzero(Qx) = ", nz_s, "  : ", "size(Qx)=", size(Qx), " = ", size(Qx)[1]^2, "  : r = ", round(nz_s/size(Qx)[1]^2,digits=4))    


# staggered SBP?
# p = 1
# q = 2*p
# vtx = T[-1 -1/sqrt(3); 1 -1/sqrt(3); 0 2/sqrt(3)]
# mask = Int[]
# xinit = []

# cub = SummationByParts.init_tri_staggered_nodes(p,q)
# xy = SymCubatures.calcnodes(cub,vtx)
# SummationByParts.plot_tri_nodes(x=xy, vtx=vtx, q=q, n=cub.numnodes,write_title=true,label_nodes=false)
# println(cub)

# cub_omega, _ = SummationByParts.getTriCubatureOmega(q-1, T)
# xy = SymCubatures.calcnodes(cub_omega,vtx)
# SummationByParts.plot_tri_nodes(x=xy, vtx=vtx, q=q, n=cub_omega.numnodes,write_title=true,label_nodes=false)

# xinit = [0.944532850510785, 0.061319596747426944, 0.3670399210288128,
#         0.9643317137065957, 0.6791559176167298, 0.2800879551160052, 0.1108945564420906, 0.11759149845575509, 0.5898408487390895,
#         0.1, 0.8,]
# append!(mask, (9+1):(9+2)) 
# append!(mask, (10+1):(10+7)) 
# xinit = [0.08010585884468961, 0.4897462888312121, 0.9464080771293774, 0.9381417384706295, 0.9078579671102778, 0.10037096010366602, 0.4791188903183325, 0.04274212618449988, 0.24373071220777096, 0.1495880562226898, 9.894989818027191e-6, 0.11529299103603494]
# xinit=[4/5,1/3]
# cub, _ = SummationByParts.deriveTriCubatureOmega(q=q-1,vertices=false,
#                                                      midedges=false,
#                                                      centroid=false,
#                                                      numedge=0,
#                                                      numS21=1,
#                                                      numS111=0, 
#                                                      xinit=xinit, 
#                                                      mask=mask, delta1=1e-2, delta2=1e-2, verbose=true)
# xy = SymCubatures.calcnodes(cub,vtx)
# SummationByParts.plot_tri_nodes(x=xy, vtx=vtx, q=q-1, n=cub.numnodes, write_title=true,label_nodes=false)


# p =7
# q = 2*p-1
# vertices =false
# midedges = false
# numedge = 0
# numS21 = 4
# numS111 = 4
# centroid = true
# mask = Int[]
# # xinit =[0.11946629335997022, 0.3416879952650081, 0.5747680456926114, 0.7,
# #         1.5087470795599724, 0.10666597054658086, 1.15236027848101, 0.09573484698587671, 1.025654579962154, 0.30250431498538716, 
# #         0.03881033382161821, 0.08154407642608524, 0.07780339023362927, 0.04,
# #         0.05674688135939928, 0.06187386430322001, 0.0869035589388141,
# #         0.05]

# # numparams = numedge+numS21+numS111*2
# # numweights = convert(Int, vertices)+convert(Int, midedges)+convert(Int, centroid)+ numedge+numS21+numS111
# # append!(mask, (numparams+1):(numparams+numweights)) 
# xinit = 2.0.*[0.489076946452539,0.22137228629183290,0.426941414259800,0.0215096811088,
#          0.16359740106785,0.087895483032197,0.02437018690109,0.110922042803463,
#          0.0680122435542, 0.3084417608921, 0.27251581777342, 0.0051263891023]
# cub, _ = SummationByParts.deriveTriCubatureGamma(q=q,vertices=vertices, 
#                                                      midedges=midedges,
#                                                      numedge=numedge,
#                                                      numS21=numS21,
#                                                      numS111=numS111,
#                                                      centroid=centroid,
#                                                      xinit=xinit, mask=mask, delta1=1e-3, delta2=1e-3, verbose=true)
# xy = SymCubatures.calcnodes(cub,vtx)
# SummationByParts.plot_tri_nodes(x=xy, vtx=vtx, q=q, n=cub.numnodes, write_title=true,label_nodes=false)

# println("\n", cub.params,"\n")
# println(cub.weights,"\n")
# println(cub)
#-------------------------------
#3d quads
# p=1
# q=2*p
# vtx = [1 -sqrt(3)/3 -sqrt(6)/6;
#        0 2*sqrt(3)/3 -sqrt(6)/6;
#        -1 -sqrt(3)/3 -sqrt(6)/6;
#        0 0 sqrt(6)/2]



# p = 2
# q = 2*p
# vtx = [1 -sqrt(3)/3 -sqrt(6)/6; 0 2*sqrt(3)/3 -sqrt(6)/6; -1 -sqrt(3)/3 -sqrt(6)/6; 0 0 sqrt(6)/2]
# mask = Int[]
# xinit =[]
# T=Float64
# cub = SummationByParts.init_tet_staggered_nodes(p,q)
# xyz = SymCubatures.calcnodes(cub,vtx)
# SummationByParts.plotly_tet_nodes(q=q, x=xyz, vtx=vtx)
# println(cub)
# checkInteriorNodeLocaton(cub)
# mindist = calcminnodedistance(cub, vtx)
# println("\n","mindist = ", mindist)
# println("\n","minweight = ", minimum(cub.weights))

# cub_omega, _ = SummationByParts.getTetCubatureOmega(q-1, T)
# xyz = SymCubatures.calcnodes(cub_omega,vtx)
# SummationByParts.plotly_tet_nodes(x=xyz, q=q-1, vtx=vtx)

# xinit = [3*6.198169944546508e-3, 3*0.16077453539526, 3*0.32227652182, 3*0.084510891834,
#          2*0.15112296546004, 
#          2*0.45887144875245, 2*0.0879702523262, 2*0.043377587068, 2*0.214097932187, 2*0.1836413698099, 2*0.5983013498019]
# xinit = [0.063016995050116347, 0.49566269198076035, 0.967171316315405, 0.13671661443396887, 
#          0.2362102644786683, 
#          0.89173236608943954, 0.16033488492260764, 0.06636582995628885, 0.42998261112028685, 0.36434367871752205, 1.19975284321548]


# cub_tri,_ = SummationByParts.Cubature.getTriCubatureDiagE(10, T, vertices=true)
# cub_tri,_ = SummationByParts.Cubature.getTriCubatureForTetFaceDiagE(2, T, faceopertype=:DiagE)
# xedge = cub_tri.params
# vertices = true
# numS31=6
# midedges=true
# numS22=1
# numfaceS21=3
# numedge=1
# numS211=6
# numfaceS111=1
# facecentroid=true
# numS1111=0
# centroid=true
# T=Float64

# cub_tri,_ = SummationByParts.Cubature.getTriCubatureDiagE(10, T, vertices=true)
# cub_tri,_ = SummationByParts.Cubature.getTriCubatureForTetFaceDiagE(10, T, faceopertype=:DiagE)
# xedge = cub_tri.params
# vertices = true
# numS31=3
# midedges=true
# numS22=1
# numfaceS21=3
# numedge=1
# numS211=7
# numfaceS111=1
# facecentroid=true
# numS1111=0
# centroid=true
# T=Float64

# cub_tri,_ = SummationByParts.Cubature.getTriCubatureForTetFaceDiagE(10, T, faceopertype=:Omega)
# xedge = cub_tri.params
# vertices =false
# numS31=6
# midedges=false
# numS22=3
# numfaceS21=2
# numedge=0
# numS211=6
# numfaceS111=3
# facecentroid=false
# numS1111=0
# centroid=true
# T=Float64


# xinit_param = T[0.11472743760821089, 0.5318204555864082, 0.965576980321931, 0.17916702722411257, 0.6134897791764771, 0.8553678267622892, 0.21738047101393035, 0.8312104497931875, 0.6589744120902207, 0.9540069927604912, 0.09149922166977983, 0.5272375124699329, 0.8540637650586674, 0.29746375857253127, 0.7370347878007151, 0.94987899770573, 0.8385931397553689, 0.6815587319130891, 0.07895898982199448, 1.5377270240632477, 0.051716990985179706, 1.163799449769477, 0.288060240995325, 0.08389843163197085, 0.0656000901160755, 0.43100467872680415, 0.3627623183625254, 1.2011673771365483, 0.7445836246919536, 0.056291285091435266, 0.9186717278941127, 0.16011659785744337, 0.37610345455698846, 0.18248445516545148, 0.7126333016167006, 0.16565383779695866, 0.08118879249452095, 0.6045976879703241, 0.2993309303393121, 0.09452178262313014, 0.5659292908830563, 0.2642822989570226, 0.5696340778390068, 1.1426176900435507, 0.23954262607339435]
# xinit_weight = T[3.288504964424688e-5, 9.753041350813653e-5, 0.00011520484609071461, 0.00010476000783738548, 0.0044234582385603224, 0.0098176277012563, 0.0072899393477029945, 7.101794143069975e-6, 0.00017792719974400249, 0.017203249761615223, 0.012908744946507568, 0.000290557629825316, 0.0007841621838219072, 0.0016300200790781239, 0.0019961857450261464, 0.000899155551780939, 0.0002463563155429928, 0.00010833661646919797, 9.622789213770869e-5, 0.00020678264241689854, 0.005286165400814374, 0.005879526838087039, 0.0075787185314118485, 0.0015884164248744362, 0.0006875075779134724, 0.009204582917336123, 0.0003678116694141033, 0.017519924317651724, 0.012598790850331037, 0.00011030144673854798, 0.0008571460223324479, 0.000151793480467417, 0.00010931872793260245, 0.009667232918942639, 0.0010821707403974733]
# filter!(x -> !(x in [1,2,3,4,7,31,32,25,26,27,28]), mask)


# xinit_param = [0.11472743760821089, 0.5318204555864082, 0.965576980321931, 0.17916702722411257, 0.8626563866259739, 0.3081367269710948, #S31
#                0.21738047101393035, 0.2972099577496514, 0.041421657753863646, #S22
               
#                0.06411074643388708, 0.28432220211312814, #faceS21
               
#                0.9186717278941127, 0.16011659785744337, 0.0656000901160755, 0.43100467872680415, 0.3627623183625254, 1.2011673771365483, #S211
#                0.10741401623659595, 1.1814586290264468, 0.40459594052401315, 0.13803326752119652, 0.7791763504517007, 0.054103041947293255, 
               
#                0.296265771567641, 1.0601082378546884, 0.32740346747436544, 0.05673533067987688, 0.05923977897745965, 0.7382935636556222, #faceS111
#                ]


# xinit_weight = [0.000542407471580597, 0.01516837459029912, 0.006718956722451832, 0.00898046603227968, 0.024070742821193813, 0.013218405584518676, #S31
#                 0.013305571208276851, 0.01728372530149623, 0.003838847057560028, #S22
#                 0.0004983121100224209, 0.004156144071866689, #faceS21
                
#                 0.001522510382645644, 0.0030038628813855574, 0.006652271690799391,  #S211
#                 0.015488962371333698, 0.016678687525364597, 0.015340774007776238, 
                
#                 0.0021592460098505115, 0.0012182427564368925, 0.0004028430527669785,  #faceS111
#                 0.0011407948884045273, #centroid
#                 ]               

# xinit_param = [0.11472743760821089, 0.5318204555864082, 0.965576980321931, 0.17916702722411257, 0.8646715310822238, #S31
#                0.21738047101393035, #S22 
#                0.40849548761177584, 0.1693510205976773, 0.895541492424858, #faceS21
#                0.9192858150590572,  #edge 
#                0.9186717278941127, 0.16011659785744337, 0.0656000901160755, 0.43100467872680415, 0.3627623183625254, 1.2011673771365483, #S211 
#                0.1449174650817117, 1.4729029925354282, 
#                0.09522522712004036, 1.201643313568532, 0.38845937766444344, 0.1121114795276667, 0.7930646372397318, 0.12122513715140107, 
#                0.5199457984883094, 0.07662224108118755 #faceS111
#                ]

# xinit_weight =[7.345525580861541e-5, #vert 
#                0.0011609597417641204, 0.030694762804852258, 0.01586824410304884, 0.002188078543807912, 0.03288433111408857, #S31 
#                0.0011579437744899314, #mid 
#                0.0037597723764178936, #S22 
#                0.004355822004338593, 0.0027306354173336797, 0.0037148584900241663,#faceS21
#                0.0004962636440179242, #edge 
               
#                0.0020604877261022666, 0.005417960467572003, 0.006076974426895572, #S211 
#                0.005331043239013086, 
#                0.01093021303948351, 0.015967278409681235, 0.021050526873275653, 
               
#                0.0006233214335578109, #faceS111 
#                0.004934184581650462, #facecent 
#                6.649258956432584e-5, #cent 
#                ]
   

# xinit_param = [0.9471034493346083, 9.384845267719673e-5, 0.33124138425491634, 0.20150491430760925, 0.7366101005688043, 0.5820307535676872, #S31 
#                0.10097964519679277, 0.24939922192711989, 0.27082379296062964, 0.7061452963466225, #S22 
#                0.16099183834007516, #faceS21 
#                0.800367892880542, #edge 
#                0.3776676620520021, 0.09432140072199578, 0.04253094508296646, 0.29327762763696985, #S211 
#                0.052749604178576735, 0.4544880142731505, 
#                0.6058255660767269, 0.21518364356973504, #faceS111
#                ]
               

# xinit_weight = [0.0004406289388671523, 
#                 0.03265938927514943, 0.00012432802251942217, 0.00025016465308837915, 0.017745056616470755, 0.0012245031557392306, 0.03589629191515014, #S31
#                 4.5951359225529654e-5, 
#                 0.025344219097630277, 0.02993435321381897, 3.697962011359693e-5, 1.5397118713940708e-5, #S22
#                 1.6127233871850176e-5, #faceS21 
#                 0.0003073861422526284, #edge 
#                 0.03171821859261543, 0.002465518593346626, #S211 
#                 0.005744819210837762, 
#                 0.004878421112495816, #faceS111 
#                 0.00892401816444221, #facecent
#                 0.009406581518871167]


# xinit_param = [0.9471034493346083, 0.2018443265439173, 0.5853975925976131, 
#                0.10097964519679277, 0.2493912276328112, 
#                0.16099183834007516, 
#                0.800367892880542, 
#                0.3776676620520021, 0.09432140072199578, 0.04253094508296646, 0.29327762763696985, 0.052682687854065775, 0.45279652549893573, 
#                0.6058255660767269, 0.21518364356973504]

# xinit_weight = [0.0005645066429057503, 
#                 0.03294947416904794, 0.018075007721278555, 0.037426628927518545, 
#                 1.913274212772993e-6, 
#                 0.02559081188395829, 0.029522271850699323, 
#                 0.00014736810493388766, 
#                 0.000406514921311113, 
#                 0.03188006451598216, 0.002083547848478684, 0.005835780135782155, 
#                 0.004827897191134791, 
#                 0.008875576350711288, 
#                 0.01096973712917172]


# xinit_param = [0.5659668188394116, 0.5157517849125112, 0.9668538043763143, #S31 
#                0.17887605340644788, 0.8706611357229256, 0.27211387974453455, 
#                0.21738047101393035, #S22 
#                0.40849548761177584, 0.1693510205976773, 0.895541492424858, #faceS21 
#                0.9192858150590572, #edge 
#                0.0846814812674472, 0.04046668289033655, 0.0656000901160755, 0.43100467872680415, #S211 
#                0.20685562107102665, 0.061511434675817794, 0.0961068596984431, 1.198231392770597, 
#                0.3910740962304225, 0.09470203572608848, 0.7884648818635481, 0.11573213640110705, 
#                0.5199457984883094, 0.07662224108118755, #faceS111 
#                ]
# xinit_weight = [0.00010815435808274956, #vert 
#                 0.00016332853559682236, 0.03475208580378614, 0.011910895143152505, #S31 
#                 0.005215776924427568, 0.03333673684748539, 0.01140461840114329, 
#                 0.0011388792236087205, #mid 
#                 0.007081614468570737, #S22 
#                 0.00334191363035361, 0.002766261334638904, 0.003147799147282123, #faceS21 
#                 0.0005168026388415012, #edge 
#                 0.002513522013095886, 0.006025086692341678, 0.009424116574390842, #S211 
#                 0.011173954759548715, 0.011287034765991303, 0.020928765404119958, 
#                 0.0005013903568232711, #faceS111 
#                 0.005834464921372579, #facecentroid 
#                 0.007569675349055752, #centroid 
#                 ]

# xinit_param = [0.5659668188394116, 0.5178097678183219, 0.9628210763605295, #S31 
#                0.21329830110742443, 0.8742283150073534, 0.259968614066233, 
#                0.21738047101393035, #S22
#                0.40849548761177584, 0.1693510205976773, 0.895541492424858, #faceS21 
#                0.9192858150590572, #edge 
#                0.7807167565199028, 0.0846814812674472, 0.0656000901160755, 0.43100467872680415, #S211 
#                0.04046668289033655, 0.20685562107102665, 0.0961068596984431, 1.198231392770597, 
#                0.3910740962304225, 0.09213243598209726, 0.7975795133937595, 0.11165705230083325, 
#                0.5199457984883094, 0.07662224108118755 #faceS111 
#                ]
   
# xinit_param =[0.5659668188394116, 0.9038385869189673, 0.23315089190324234, #S31 
#               0.17626980760211586, #S22 
#               0.40849548761177584, 0.1693510205976773, 0.895541492424858, #faceS21 
#               0.9192858150590572, #edge 
#               0.7807167565199028, 0.0846814812674472, 0.04046668289033655, 0.20685562107102665, #S211 
#               0.061511434675817794, 0.6927211311947881, 0.14342099002644818, 0.3786528536651681, 
#               0.3596513879533776, 0.02013276646659172, 0.11457584452278376, 0.4753447793598157, 
#               0.43092643401573644, 0.1642273614407805, 
              
#               0.5199457984883094, 0.07662224108118755, #faceS111 
#               ]

# xinit_weight = [0.00019729611741621086, #vert 
#                 0.01124293827557943, 0.019327859474738186, 0.007304105899451999, #S31 
#                 0.0005813450370534768, #mid 
#                 0.019703716011470255, #S22 
#                 0.0027220095911565296, 0.0015022574567931319, 0.0050685022079231925, #faceS21 
#                 8.2466680692705518e-4, #edge 
#                 0.012918914540752316, 0.002947425992323903, 0.005651611481472151,#S211
#                 0.009092220233367846, 0.006300144344667722, 
#                 0.007432190958512254, 0.02582065522717608, 
                
#                 0.0011421957175394295, 0.00753620599719168, 0.04827850233716233]

# xinit = collect(Iterators.flatten([xinit_param,xinit_weight]))
# cub = SymCubatures.TetSymCub{T}(vertices=vertices,
#                                 numS31=numS31,
#                                 midedges=midedges,
#                                 numS22=numS22,
#                                 numfaceS21=numfaceS21,
#                                 numedge=numedge,
#                                 numS211=numS211,
#                                 numfaceS111=numfaceS111,
#                                 facecentroid=facecentroid,
#                                 numS1111=numS1111,
#                                 centroid=centroid)
#                                 mask = SymCubatures.getInternalParamMask(cub)

# SymCubatures.setparams!(cub,xinit_param)
# SymCubatures.setweights!(cub,xinit_weight)
# mask = SymCubatures.getInternalParamMask(cub)
# append!(mask, (cub.numparams+1):(cub.numparams+cub.numweights))
# a = -2
# filter!(x -> !(x in [1,7-a,8-a,9-a,10-a,11-a,12-a,13-a,14-a,15-a,16-a]), mask)

# p=5
# q = 2*p
# cub, _ = SummationByParts.deriveTetCubatureDiagE(q=q,
#                                                  vertices=vertices,
#                                                  numS31=numS31,
#                                                  midedges=midedges,
#                                                  numS22=numS22,
#                                                  numfaceS21=numfaceS21,
#                                                  numedge=numedge,
#                                                  numS211=numS211,
#                                                  numfaceS111=numfaceS111,
#                                                  facecentroid=facecentroid,
#                                                  numS1111=numS1111,
#                                                  centroid=centroid,
#                                                  delta1=1e-4,
#                                                  delta2=1e-4,
#                                                  verbose=true,
#                                                  xinit=xinit,
#                                                  xedge=xedge,
#                                                  mask=mask,
#                                                  tol=3e-15)
# println("\n", cub.params,"\n")
# println("\n", cub.weights,"\n")
# println(cub)
# xyz = SymCubatures.calcnodes(cub,vtx)
# SummationByParts.plotly_tet_nodes(x=xyz, q=q, vtx=vtx)


#-----------------
# p = 3
# q = 6 
# T = Float64

# cub_omega, vtx = SummationByParts.Cubature.getTriCubatureOmega(q-1, T)

# vertices = true
# midedges = true
# numedge = 1
# numS21 = 2
# numS111 = 1
# centroid = true

# cub_diage = SymCubatures.TriSymCub{T}(vertices=vertices,
#                                 midedges=midedges,
#                                 numS21=numS21,
#                                 numedge=numedge,
#                                 numS111=numS111,
#                                 centroid = centroid)
# SymCubatures.setparams!(cub_diage, T[0.13862330627662678, 0.4765740799486239, 0.8273268353539885, 0.14215944055500324, 0.6226442585632832])
# SymCubatures.setweights!(cub_diage, T[0.0019263750791820421, 0.024548784813778073, 0.08983957252485644, 0.17471428465551592, 0.01736767414771074, 0.1514579304646647, 0.11395932110574997])

# xy_omega = SymCubatures.calcnodes(cub_omega,vtx)
# w_omega = SymCubatures.calcweights(cub_omega)
# H_omega = diagm(w_omega)
# xy_diage = SymCubatures.calcnodes(cub_diage,vtx)
# w_diage = SymCubatures.calcweights(cub_diage) 
# H_diage = diagm(w_diage)

# Vp_omega, _,_ = OrthoPoly.vandermonde(p, xy_omega[1,:],xy_omega[2,:])
# Vp1_omega, _,_ = OrthoPoly.vandermonde(p-1, xy_omega[1,:],xy_omega[2,:])
# Vp_diage, _,_ = OrthoPoly.vandermonde(p, xy_diage[1,:],xy_diage[2,:])
# Vp1_diage, _,_ = OrthoPoly.vandermonde(p-1, xy_diage[1,:],xy_diage[2,:])

# np = binomial(p+2,2)
# Vlarge, _,_ = OrthoPoly.vandermonde(p+5, xy_omega[1,:],xy_omega[2,:])
# W_omega = Vlarge[:,np+1:cub_omega.numnodes]

# # W_omega = I(cub_omega.numnodes)[:,np+1:cub_omega.numnodes]

# W_diage = (W_omega' * H_omega*Vp1_omega * inv(Vp1_diage' * Vp1_diage) *Vp1_diage' * inv(H_diage))'

# Is2f = hcat(Vp_diage,W_diage) * inv(hcat(Vp_omega,W_omega))
# If2s = inv(H_omega) * Is2f' * H_diage

# norm(If2s*Vp1_diage - Vp1_omega)
# norm(Is2f*Vp_omega - Vp_diage)

#--------------
# p = 1
# q = 2*p
# cub_omega,cub_diage,vtx = SummationByParts.Cubature.getTriCubatureStaggered(q)
# facecub, facevtx = SummationByParts.Cubature.quadrature(2*p, T, internal=false)
# R_omega, perm = SummationByParts.buildfacereconstruction(facecub, cub_omega, vtx, p,faceopertype=:Omega)
# R_diage, perm = SummationByParts.buildfacereconstruction(facecub, cub_diage, vtx, p,faceopertype=:DiagE)

# xy_omega = SymCubatures.calcnodes(cub_omega,vtx)
# w_omega = SymCubatures.calcweights(cub_omega)
# H_omega = diagm(w_omega)
# xy_diage = SymCubatures.calcnodes(cub_diage,vtx)
# w_diage = SymCubatures.calcweights(cub_diage) 
# H_diage = diagm(w_diage)

# Vp_omega, _,_ = OrthoPoly.vandermonde(p, xy_omega[1,:],xy_omega[2,:])
# Vp1_omega, _,_ = OrthoPoly.vandermonde(p-1, xy_omega[1,:],xy_omega[2,:])
# Vp_diage, _,_ = OrthoPoly.vandermonde(p, xy_diage[1,:],xy_diage[2,:])
# Vp1_diage, _,_ = OrthoPoly.vandermonde(p-1, xy_diage[1,:],xy_diage[2,:])

# np = binomial(p+2,2)
# Vlarge, _,_ = OrthoPoly.vandermonde(p+10, xy_omega[1,:],xy_omega[2,:])
# W_omega = Vlarge[:,np+1:cub_omega.numnodes]

# # W_omega = I(cub_omega.numnodes)[:,np+1:cub_omega.numnodes]

# W_diage = (W_omega' * H_omega*Vp1_omega * inv(Vp1_diage' * Vp1_diage) *Vp1_diage' * inv(H_diage))'

# Is2f = hcat(Vp_diage,W_diage) * inv(hcat(Vp_omega,W_omega))
# If2s = inv(H_omega) * Is2f' * H_diage
# norm(R_omega*Vp_omega - Is2f[perm[:,1],:]*Vp_omega)
#--------------
# p = 1
# q = 2*p
# cub_omega,cub_diage,vtx = SummationByParts.Cubature.getTetCubatureStaggered(q)
# facecub,facevtx = SummationByParts.Cubature.getTriCubatureForTetFaceDiagE(q, faceopertype=:DiagE)
# R_omega, perm = SummationByParts.buildfacereconstruction(facecub, cub_omega, vtx, p,faceopertype=:Omega)
# R_diage, perm = SummationByParts.buildfacereconstruction(facecub, cub_diage, vtx, p,faceopertype=:DiagE)

# xy_omega = SymCubatures.calcnodes(cub_omega,vtx)
# w_omega = SymCubatures.calcweights(cub_omega)
# H_omega = diagm(w_omega)
# xy_diage = SymCubatures.calcnodes(cub_diage,vtx)
# w_diage = SymCubatures.calcweights(cub_diage) 
# H_diage = diagm(w_diage)

# Vp_omega, _,_ = OrthoPoly.vandermonde(p, xy_omega[1,:],xy_omega[2,:],xy_omega[3,:])
# Vp1_omega, _,_ = OrthoPoly.vandermonde(p-1, xy_omega[1,:],xy_omega[2,:],xy_omega[3,:])
# Vp_diage, _,_ = OrthoPoly.vandermonde(p, xy_diage[1,:],xy_diage[2,:],xy_diage[3,:])
# Vp1_diage, _,_ = OrthoPoly.vandermonde(p-1, xy_diage[1,:],xy_diage[2,:],xy_diage[3,:])

# # np = binomial(p+3,3)
# # Vlarge, _,_ = OrthoPoly.vandermonde(p+5, xy_omega[1,:],xy_omega[2,:],xy_omega[3,:])
# # W_omega = Vlarge[:,np+1:cub_omega.numnodes]
# # W_omega = I(cub_omega.numnodes)[:,np+1:cub_omega.numnodes]
# # Q,R=qr(Vp_omega)
# # W_omega = Q[:,size(R,1)+1:size(Vp_omega,1)]

# W_omega = nullspace(Vp_omega')
# W_diage = (W_omega' * H_omega*Vp1_omega * inv(Vp1_diage' * Vp1_diage) *Vp1_diage' * inv(H_diage))'

# Is2f = hcat(Vp_diage,W_diage) * inv(hcat(Vp_omega,W_omega))
# If2s = inv(H_omega) * Is2f' * H_diage
# norm(R_omega*Vp_omega - Is2f[perm[:,1],:]*Vp_omega)

#-----------
# p = 5
# q = 2*p
# # cub_omega,cub_diage,vtx = SummationByParts.Cubature.getTriCubatureStaggered(q)
# cub_omega,cub_diage,vtx = SummationByParts.Cubature.getTetCubatureStaggered(q)
# xy_omega = SymCubatures.calcnodes(cub_omega,vtx)
# w_omega = SymCubatures.calcweights(cub_omega)
# xy_diage = SymCubatures.calcnodes(cub_diage,vtx)
# w_diage = SymCubatures.calcweights(cub_diage)

# println("\n")
# println("// ",cub_omega.numnodes, " SBP-omega nodes and ", cub_diage.numnodes, " SBP-DiagE nodes."  )
# println("p_cub_o = ", q-1,";")
# println("p_cub_e = ", q,";")
# cppshow(xy_omega,lhs="x_cart_o = ")
# cppshow(w_omega,lhs="cub_w_o = ")
# cppshow(xy_diage,lhs="x_cart_e = ")
# cppshow(w_diage,lhs="cub_w_e = ")


#--------------
# p = 12
# q = 2*p+1
# T = Float64
# vtx = T[-1 -1/sqrt(3); 1 -1/sqrt(3); 0 2/sqrt(3)]

# cub = SummationByParts.init_tri_nodes_omega(p, q)
# xy = SymCubatures.calcnodes(cub,vtx)
# SummationByParts.plot_tri_nodes(x=xy, vtx=vtx, q=q, n=cub.numnodes,write_title=true,label_nodes=false)
# println(cub)

# println("------------")
# mindist = calcminnodedistance(cub, vtx)
# w = SymCubatures.calcweights(cub)
# println("\n","mindist = ", mindist)
# println("\n","minweight = ", minimum(w))
# checkInteriorNodeLocaton(cub)

# #--------------
# p = 2
# q = 2*p+1
# vtx = [1 -sqrt(3)/3 -sqrt(6)/6; 0 2*sqrt(3)/3 -sqrt(6)/6; -1 -sqrt(3)/3 -sqrt(6)/6; 0 0 sqrt(6)/2]

# cub = SummationByParts.init_tet_nodes_omega(p, q)
# xyz = SymCubatures.calcnodes(cub,vtx)
# SummationByParts.plotly_tet_nodes(q=q, x=xyz, vtx=vtx)
# println(cub)

# println("------------")
# mindist = calcminnodedistance(cub, vtx)
# w = SymCubatures.calcweights(cub)
# println("\n","mindist = ", mindist)
# println("\n","minweight = ", minimum(w))
# checkInteriorNodeLocaton(cub)

# ----------------
# p = 19
# q = 2*p+0
# qq = 2*p+2
# T = Float64
# vtx = T[-1 -1/sqrt(3); 1 -1/sqrt(3); 0 2/sqrt(3)]

# cub, _ = SummationByParts.Cubature.getTriCubatureOmegaLG(qq)
# cub, res = SummationByParts.eliminate_nodes(cub, p, q)

# xy = SymCubatures.calcnodes(cub,vtx)
# SummationByParts.plot_tri_nodes(x=xy, vtx=vtx, q=q, n=cub.numnodes,write_title=true,label_nodes=false)
# println(cub)

# mindist = calcminnodedistance(cub, vtx)
# w = SymCubatures.calcweights(cub)
# println("\n","mindist = ", mindist)
# println("\n","minweight = ", minimum(w))
# checkInteriorNodeLocaton(cub)
# println("-------------------")
# ---------------
function get_tet_omega_lg()
    p = 5
    q = 2*p-0
    qq = 2*p+2*mod(p,2)
    T = Float64
    vtx = [1 -sqrt(3)/3 -sqrt(6)/6; 0 2*sqrt(3)/3 -sqrt(6)/6; -1 -sqrt(3)/3 -sqrt(6)/6; 0 0 sqrt(6)/2]

    cuba, _ = SummationByParts.Cubature.getTetCubatureOmegaLG(qq)
    # cub, res = SummationByParts.eliminate_nodes(cuba, p, q)
    cubb, _ = SummationByParts.Cubature.getTetCubatureOmegaLG(qq-2)
    cub = SummationByParts.combine_tet_cubs(cuba, cubb)
    cub,res = SummationByParts.eliminate_nodes(cub, p, q)

    xyz = SymCubatures.calcnodes(cub,vtx)
    SummationByParts.plotly_tet_nodes(q=q, x=xyz, vtx=vtx)
    # dir = "/project/z/zingg/workuzel/quadrature/"
    dir =  joinpath(pwd(), "src/")
    file = joinpath(dir,"tet_omega_lg.dat")

    open(file, "a") do io 
        redirect_stdout(io) do
            println("------------------------")
            println("res = ", res)
            println("p = ", p, ":  q = ", q)
            println(cub.params)
            println(cub.weights)
            println(cub)

            mindist = calcminnodedistance(cub, vtx)
            w = SymCubatures.calcweights(cub)
            println("\n","mindist = ", mindist)
            println("\n","minweight = ", minimum(w))
            checkInteriorNodeLocaton(cub)
            println("------------------------","\n")
        end
    end

    mindist = calcminnodedistance(cub, vtx)
    w = SymCubatures.calcweights(cub)
    println("\n","mindist = ", mindist)
    println("\n","minweight = ", minimum(w))
    checkInteriorNodeLocaton(cub)
end

get_tet_omega_lg()

#---------------
# function derive_tri_omega_lg()
#     for i = 1:2
#         p = i
#         q = 2*p+1
#         T = Float64
#         vtx = T[-1 -1/sqrt(3); 1 -1/sqrt(3); 0 2/sqrt(3)]

#         cub, res = SummationByParts.init_tri_nodes_omega(p, q)

#         # dir = "/project/z/zingg/workuzel/quadrature/"
#         dir =  joinpath(pwd(), "src/")
#         file = joinpath(dir,"tri_omega_lg_2p1.dat")

#         open(file, "a") do io 
#             redirect_stdout(io) do
#                 println("------------------------")
#                 println("res = ", res)
#                 println("p = ", p, ":  q = ", q)
#                 println(cub.params)
#                 println(cub.weights)
#                 println(cub)

#                 mindist = calcminnodedistance(cub, vtx)
#                 w = SymCubatures.calcweights(cub)
#                 println("\n","mindist = ", mindist)
#                 println("\n","minweight = ", minimum(w))
#                 checkInteriorNodeLocaton(cub)
#                 println("------------------------","\n")
#             end
#         end

#         mindist = calcminnodedistance(cub, vtx)
#         w = SymCubatures.calcweights(cub)
#         println("\n","mindist = ", mindist)
#         println("\n","minweight = ", minimum(w))
#         checkInteriorNodeLocaton(cub)
#     end
# end

# derive_tri_omega_lg()