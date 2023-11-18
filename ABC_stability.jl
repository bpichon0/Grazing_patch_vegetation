
using Distributed

addprocs(60, exeflags="--project=$(Base.active_project())")

@everywhere begin
    using StatsBase, RCall, Random, LaTeXStrings
    using BenchmarkTools, Images, Tables, CSV, LinearAlgebra, Distributions, DataFrames
end

@everywhere include("./Drylands_ABC_functions.jl")


@everywhere function Distance_tipping(id)

    step_size = 1 / 200
    n_sample=100
    posteriors=readdlm("../Data/Inferrence/param_inferred.csv", ';')[2:end,2:end]
    post_p = posteriors[:, id]
    post_q = posteriors[:, id+504]
    Keeping_data = zeros(Int((1/step_size))*n_sample, 3)
    index=1
    param=zeros(2)
    for sample_id in 1:n_sample

        if rand(Distributions.Uniform(0, 1)) < .5
          println("")
          GC.gc()
          ccall(:malloc_trim, Cvoid, (Cint,), 0)
          GC.safepoint()
        end
        sample_row=rand(1:length(post_p))
        param[1] = post_p[sample_row]
        param[2] = post_q[sample_row]

        p_to_desert = push!(collect(0:step_size:param[1]),param[1])

        for pcrit_id in 1:length(p_to_desert)

            param[1]=p_to_desert[pcrit_id]

            Keeping_data[index, 1:2] .= param
            fraction_cover = [0.8, 0.2]
            
            size_landscape = 80
            ini_land = Get_initial_lattice_Eby(frac=fraction_cover, size_mat=size_landscape)

            d1, land1 = IBM_Eby_model(time_t=50, param=copy(param), landscape=copy(ini_land),
            keep_landscape=true, burning_phase=2000, intensity_feedback=6)

            mean_cover = mean([length(findall(land1[:, :, k] .== 1)) / (size(land1)[1] * size(land1)[2]) for k in 1:size(land1)[3]])

            Keeping_data[index, 3] = ifelse(any([length(findall(land1[:, :, k] .== 1)) == 0 for k in 1:size(land1)[3]]), 0, mean_cover)
            
            index += 1
            println(index)
            end
    end


    CSV.write("../Data/Inferrence/Prediction/Dist_tipping_" * repr(id) * ".csv", Tables.table(Keeping_data), writeheader=false)
end


sites=readdlm("../Data/Inferrence/Keeping_sites.csv",';',Int64)[:,1]
print(sites)
pmap(Distance_tipping, sites)



