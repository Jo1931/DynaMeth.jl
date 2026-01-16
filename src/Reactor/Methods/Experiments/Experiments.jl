struct DynamicExp <: AbstractReactor
    number
    files
end
function DynamicExp()
    numb = (22, 23, 24, 25, 28, 29, 30, 31, 32, 33, 34, 35, 37, 39, 40, 41, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 73, 74, 76, 77, 79, 81)
    numb = (19, 20, 21, 23, 24, 25, 26, 27, 28, 29, 31, 33, 35, 36, 37, 39, 40, 41, 43, 45, 46, 47, 49, 52, 54, 55)#,73,74,76,77,81)
    numb = (23, 24, 25, 26, 28, 29, 31, 33, 34, 35, 37, 39, 40, 41, 43, 45, 47, 49, 52, 54, 55)
    files = []
    path = "C:/Users/Johan/Documents/20_OVGU/10_SPP2080/21_CSTR_fit/Messdaten"
    for i in numb
        file = matread(path * "/Messdaten/MeOH_0" * string(i) * ".mat")
        push!(files, file)
    end
    DynamicExp(
        numb,
        files)
end
function DynamicVollbrecht()
    return DataFrame(CSV.File("data/exp_pro/experiment.csv"))
end
mutable struct Experiments <: AbstractReactor
    xin_mess
    xout_mess
    xout_sim
    files
    idx
    objectiveSS
    objectiveDyn
    subset
end
