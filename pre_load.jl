# This is the file for precompiling self-created packages in niagara system

if !("/home/z/zingg/workuzel/SummationByPartsPrivate/src/" in LOAD_PATH)
    push!(LOAD_PATH, "/home/z/zingg/workuzel/SummationByPartsPrivate/src/") 
end
ENV["JULIA_NUM_THREADS"] = "40"
using SummationByParts