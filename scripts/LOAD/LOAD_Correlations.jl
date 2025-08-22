using DrWatson
@quickactivate "QuantumFPUT"


filename = "HPC-LONG_Corr__N=4_dim=15_α=0.0_β=0.2_state=coherent.jld2"
CorrSites = [1,3]
colors = ["k","lightcoral","teal","tab:orange","tab:cyan","tab:brown","tab:pink"]
colorlist = ["midnightblue","orange","crimson"]

fs=25
ts=20


#####################################
##########  Loading File   ##########
#####################################

using JLD2
import QuantumOptics: Ket, Operator, FockBasis

parameters    = load(datadir("Correlations", filename), "parameters")
periods       = load(datadir("Correlations", filename), "periods")
initialstate  = load(datadir("Correlations", filename), "initialstate")
TotalCorrDict = load(datadir("Correlations", filename), "TotalCorrDict")
ClassCorrDict = load(datadir("Correlations", filename), "ClassCorrDict")
QuantCorrDict = load(datadir("Correlations", filename), "QuantCorrDict")

totalcorrelations = TotalCorrDict[CorrSites]
classcorrelations = ClassCorrDict[CorrSites]
quantcorrelations = QuantCorrDict[CorrSites]

t_start, t_end = first(periods), last(periods)


##########################################
##########   Print Parameters   ##########
##########################################
using JLD2

N = parameters["N"]
dim = parameters["dim"]
n_mcwf = get(parameters,"n_mcwf", 1)
state = parameters["state"]

κ = parameters["κ"]
α = parameters["α"]
β = parameters["β"]

γ   = get(parameters,"γ", nothing)
kbT = get(parameters,"kbT", nothing)
BathSites = get(parameters, "BathSites", nothing)

println("\n\n\n")
println("##################### ################ #####################")
println("#####################    Parameters    #####################\n")

println("Chain length is $N oscillators, initialized at site 1 as "*state*" states")
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

println("##################### ################ #####################")
println("##################### ################ #####################")
println("\n")





###################################################
##################   Plot Data   ##################
###################################################

using PyPlot, LaTeXStrings
pygui(true)


# import QuantumOptics: wigner, reduced
# rng = LinRange(-5.0, 5.0, 100)
# W = transpose(wigner(reduced(initialstate, 1), rng, rng))
# Wmax = maximum(W)

# fig, axs = subplots()
# axs.pcolormesh(rng, rng, W, vmin=-Wmax, vmax=Wmax, cmap="seismic")
# axs.set_aspect(1)
# axs.grid()
# show()


# fig, axs = subplots()
# fig.suptitle("Correlations in FPUT-Chain")

# axs.plot(periods, totalcorrelations, label="full corr. ρ$(CorrSites[1]) and ρ$(CorrSites[2])")
# # axs.plot(periods, classcorrelations, label="class. corr. ρ$(CorrSites[1]) and ρ$(CorrSites[2])")
# # axs.plot(periods, quantcorrelations, label="quant. corr. ρ$(CorrSites[1]) and ρ$(CorrSites[2])")
# axs.grid()
# axs.legend()
# axs.set_ylim([-0.01, 1.01])
# axs.set_xlim([t_start-0.01, t_end+0.01])
# show()


fig, axs = subplots()
# fig.suptitle("Correlations in FPUT-Chain")

axs.plot(periods, TotalCorrDict[[1,2]], color=colorlist[1], lw=2,  label=L"$\mathcal{H}^{1+2}$")
axs.plot(periods, TotalCorrDict[[1,3]], color=colorlist[2], lw=2,  label=L"$\mathcal{H}^{1+3}$", alpha=0.7)
axs.plot(periods, TotalCorrDict[[1,4]], color=colorlist[3], lw=2,  label=L"$\mathcal{H}^{1+4}$", alpha=0.7)
axs.grid()
axs.legend(fontsize="xx-large", loc="upper right")
axs.set_ylim([-0.01, 0.51])
axs.set_ylabel(L"$\mathcal{C}_\mathcal{H}(\rho)$", fontsize=fs)
axs.set_xlim([t_start-0.01, t_end+0.01])
axs.set_xlabel(L"time $t \cdot \frac{2\pi}{\sqrt{\kappa}}$", fontsize=fs)
axs.tick_params(labelsize=ts)
show()


