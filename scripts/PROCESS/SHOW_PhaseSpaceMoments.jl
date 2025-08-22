using DrWatson
@quickactivate "QuantumFPUT"


filename = "Dbl__HPC-SHORT__N=4_dim=15_α=0.0_β=0.0_state=coherent.jld2"
PlotSites = [1]

colorlist = ["lightcoral", "midnightblue", "mediumvioletred"]

fs=20
ts=15



#####################################
##########  Loading File   ##########
#####################################

using JLD2
import QuantumOptics: Ket, Operator, FockBasis
file = load(datadir("KetTrajectories", filename))
println(keys(file))

parameters     = file["parameters"]
periods        = file["periods"]
t_start, t_end = first(periods), 20.0 # last(periods)

# states = file["states1"]
states = file["states1"]




##########################################
##########   Print Parameters   ##########
##########################################
using JLD2

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




######################################
##########   Compute Data   ##########
######################################

using ReduceSystems
using StateCharacteristics

@time array_dms = reduce_to(states, N, PlotSites)

array_positions = Array{Float64,2}(undef, length(periods), length(PlotSites))
array_momentums = Array{Float64,2}(undef, length(periods), length(PlotSites))

array_VarsQ  = Array{Float64,2}(undef, length(periods), length(PlotSites))
array_VarsP  = Array{Float64,2}(undef, length(periods), length(PlotSites))
array_CovsQP = Array{Float64,2}(undef, length(periods), length(PlotSites))

@time for p in 1:length(PlotSites)
    dms = array_dms[:,p]

    array_positions[:,p] .= meanposition.(dms)
    array_momentums[:,p] .= meanmomentum.(dms)

    array_VarsQ[:,p]  .= covariance_qq.(dms)
    array_VarsP[:,p]  .= covariance_pp.(dms)
    array_CovsQP[:,p] .= covariance_qp.(dms)
end




###################################################
##################   Plot Data   ##################
###################################################

using PyPlot
pygui(true)

fig, axs = subplots()
fig.suptitle("Quantum phase-space meanvalues", fontsize=fs)
for p in eachindex(PlotSites)
    axs.plot(periods, array_positions[:,p], color=colorlist[1], lw=3, label=L"$\langle q \rangle_\rho$")
    axs.plot(periods, array_momentums[:,p], color=colorlist[2], lw=3, label=L"$\langle p \rangle_\rho$")
end
axs.set_xlim([t_start-0.01, t_end+0.01])
axs.set_xlabel(L"time $t \cdot \frac{2\pi}{\sqrt{\kappa}}$", fontsize=fs)
axs.tick_params(labelsize=ts)
axs.legend(fontsize="xx-large", loc="upper right")
axs.grid()
show()

fig, axs = subplots()
fig.suptitle("Quantum phase-space (co-) variances", fontsize=fs)
for p in eachindex(PlotSites)
    axs.plot(periods,  array_VarsQ[:,p], color=colorlist[1], lw=3, label=L"$\sigma_{qq}^\rho$")
    axs.plot(periods,  array_VarsP[:,p], color=colorlist[2], lw=3, label=L"$\sigma_{pp}^\rho$")
    axs.plot(periods, array_CovsQP[:,p], color=colorlist[3], lw=3, label=L"$\sigma_{qp}^\rho$")
end
axs.set_xlim([t_start-0.01, t_end+0.01])
axs.set_xlabel(L"time $t \cdot \frac{2\pi}{\sqrt{\kappa}}$", fontsize=fs)
axs.tick_params(labelsize=ts)
axs.legend(fontsize="xx-large", loc="upper right")
axs.grid()
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



# fig,axs = subplots()
# fig.suptitle("Dimension Cutoff Error")
# for p in eachindex(PlotSites)
#     axs.plot(periods, array_dimerrors[:,p], alpha=0.5, label="osc. $(PlotSites[p])")
# end
# axs.grid()
# axs.legend()
# axs.set_xlim([t_start-0.01, t_end+0.01])
# show()