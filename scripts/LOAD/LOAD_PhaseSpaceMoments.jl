using DrWatson
@quickactivate "QuantumFPUT"


filename = "HPC_PhSpMoments__N=3_dim=15_n_mcwf=200_α=0.0_β=0.2_kbT=0.0_γ=0.1_state=coherent.jld2"
PlotSites = [1,2,3]



#####################################
##########  Loading File   ##########
#####################################

using JLD2
file = load(datadir("StateCharacteristics", filename))
# println(keys(file))

## Load keys
parameters      = file["parameters"]
periods         = file["periods"]

totalenergies   = real.(file["totalenergies"])

array_Nmeans    = file["array_Nmeans"]
array_Nvars     = file["array_Nvars"]

array_meansPhSp = file["array_meansPhSp"]
array_covarPhSp = file["array_covarPhSp"]

array_dimerrors = file["array_dimerrors"]



## Minor post-processing
import LinearAlgebra: eigvals, norm, det
function squeezing(Σ::Matrix{Float64})
    σ1,σ2 = abs.(eigvals(Σ))
    return 1-σ1/σ2
end

t_start, t_end = first(periods), last(periods)
array_absPhSp = norm.(array_meansPhSp)
array_covdets = det.(array_covarPhSp)
array_covsqzs = squeezing.(array_covarPhSp)







##########################################
##########   Print Parameters   ##########
##########################################

N = parameters["N"]
dim = parameters["dim"]
n_mcwf = get(parameters,"n_mcwf.", nothing)

κ = parameters["κ"]
α = parameters["α"]
β = parameters["β"]

γ   = get(parameters,"γ", nothing)
kbT = get(parameters,"kbT", nothing)
BathSites = get(parameters, "BathSites", nothing)

abstol = get(parameters, "abstol", "1e-8")
reltol = get(parameters, "reltol", "1e-6")
alg = get(parameters, "alg", "OrdinaryDiffEq.DP5()")

println("\n\n\n")
println("##################### ################ #####################")
println("#####################    Parameters    #####################\n")

println("Chain length is $N oscillators")
println("Dim. cutoff is $dim")
println("Coupling Potential: V(r)= r²*$κ/2! + r³*$α/3! + r⁴*$β/4!")
println("Time interval from $t_start*2πω to $t_end*2πω")

if γ === nothing
    println("No Markovian baths attached\n")
else
    println("Baths temperature is kbT=$kbT")
    println("Baths sites are $BathSites")
    println("Baths coupling strength is γ=$γ")
    println("Averaged over $n_mcwf MC trajectories\n")
end

println("Solved with $alg at abstol=$abstol and reltol=$reltol")

println("##################### ################ #####################")
println("##################### ################ #####################")
println("\n")





###################################################
##################   Plot Data   ##################
###################################################

using PyPlot
pygui(true)

fig, axs = subplots()
fig.suptitle("Phase space Means")
for p in eachindex(PlotSites)
    axs.plot(periods, array_absPhSp[:,p], alpha=0.5, label="osc. $(PlotSites[p])")
end
axs.grid()
axs.legend()
axs.set_xlim([t_start-0.01, t_end+0.01])
axs.set_title("absolute of displ. |μ|")
show()


fig, axs = subplots(nrows=2,ncols=1)
fig.suptitle("Phase space Covariances")
for p in eachindex(PlotSites)
    axs[1].plot(periods, array_covdets[:,p], alpha=0.5, label="osc. $(PlotSites[p])")
    axs[2].plot(periods, array_covsqzs[:,p], alpha=0.5, label="osc. $(PlotSites[p])")
end
axs[1].grid()
axs[1].legend()
axs[1].set_xlim([t_start-0.01, t_end+0.01])
axs[1].set_title("determinant of cov. mat. det(Σ)")
axs[2].grid()
axs[2].legend()
axs[2].set_xlim([t_start-0.01, t_end+0.01])
axs[2].set_title("squeezing of cov. mat. s=1-σ1/σ2")
show()


# fig, axs = subplots(nrows=2,ncols=1)
# fig.suptitle("Hilbert space characterstics")
# for p in eachindex(PlotSites)
#     axs[1].plot(periods, array_Nmeans[:,p], alpha=0.5, label="osc. $(PlotSites[p])")
#     axs[2].plot(periods,  array_Nvars[:,p], alpha=0.5, label="osc. $(PlotSites[p])")
# end
# axs[1].grid()
# axs[1].legend()
# axs[1].set_xlim([t_start-0.01, t_end+0.01])
# axs[1].set_title("mean value of number op. <N>")
# axs[2].grid()
# axs[2].legend()
# axs[2].set_xlim([t_start-0.01, t_end+0.01])
# axs[2].set_title("variance of number op. <N²>-<N>²")
# show()


fig,axs = subplots()
fig.suptitle("Toal Energy in the Chain")
axs.plot(periods, totalenergies)
axs.grid()
axs.set_xlim([t_start-0.01, t_end+0.01])
axs.set_ylim([-0.05, maximum(totalenergies)+0.05])
show()


# fig,axs = subplots()
# fig.suptitle("Dimension Cutoff Error")
# for p in eachindex(PlotSites)
#     axs.plot(periods, array_dimerrors[:,p], alpha=0.5, label="osc. $(PlotSites[p])")
# end
# axs.grid()
# axs.legend()
# axs.set_xlim([t_start-0.01, t_end+0.01])
# show()