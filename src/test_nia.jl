if !("src/" in LOAD_PATH)
    push!(LOAD_PATH, "src/") 
end
using SummationByParts
println("Running on Niagara...")