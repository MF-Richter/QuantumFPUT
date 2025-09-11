using DrWatson
@quickactivate "QuantumFPUT"
filename1 = "HPC__N=3_dim=15_n_mcwf=1000_α=0.0_β=0.2_kbT=0.0_γ=0.1_state=coh+.jld2"
filename2 = "HPC__N=3_dim=15_n_mcwf=1000_α=0.0_β=0.2_kbT=0.0_γ=0.1_state=coh-.jld2"

PlotSites = [1,2,3]

colorlist = ["k","lightcoral","tab:orange","teal","tab:cyan","tab:brown","tab:pink"]

fs=25
ts=20



########################################################
##########  Loading and Printing Parameters   ##########
########################################################

using JLD2
import QuantumOptics: Ket, Operator, FockBasis

parameters     = load(datadir("ReducedTrajectories", filename1), "parameters")
periods        = load(datadir("ReducedTrajectories", filename1), "periods")
t_start, t_end = first(periods), last(periods)


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

abstol = get(parameters, "abstol", "1e-8")
reltol = get(parameters, "reltol", "1e-6")
# alg = get(parameters, "alg", "OrdinaryDiffEq.DP5()")
alg = "DP5()"


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
    println("Bath temperature is kbT=$kbT")
    println("Bath sites are $BathSites")
    println("Bath coupling strength is γ=$γ")
    println("Averaged over $n_mcwf MC trajectories\n")
end

println("Solved with $alg at abstol=$abstol and reltol=$reltol")

println("##################### ################ #####################")
println("##################### ################ #####################")
println("\n")




######################################
##########   Compute Data   ##########
######################################

using ReduceSystems, HamiltonianFPUT, StateCharacteristics
import QuantumOptics: tracedistance, tr, entropy_vn, FockBasis,fockstate, expect

array_dms1 = load(datadir("ReducedTrajectories", filename1), "array_dms")
array_dms2 = load(datadir("ReducedTrajectories", filename2), "array_dms")

array_tracedistances = Array{Float64,2}(undef, length(periods), length(PlotSites))

array_dimerrors   = Array{Float64,3}(undef, length(periods), length(PlotSites), 2)
array_cutoffratio = Array{Float64,3}(undef, length(periods), length(PlotSites), 2)

@time for p in eachindex(PlotSites)
    dms1 = array_dms1[:,p]
    dms2 = array_dms2[:,p]

    array_tracedistances[:,p].= tracedistance.(dms1,dms2)
    array_dimerrors[:,p,1]   .= tracenorm_diviation.(dms1)
    array_dimerrors[:,p,2]   .= tracenorm_diviation.(dms2)
    array_cutoffratio[:,p,1] .= relative_dimcutoffvalue.(dms1)
    array_cutoffratio[:,p,2] .= relative_dimcutoffvalue.(dms2)
end

sglbasis = FockBasis(dim)
HamFPUT = HamOpFPUT(N, κ,α,β, sglbasis)
groundstate = localKets(fockstate(sglbasis, 0), N)
groundenergy = expect(HamFPUT, groundstate)




###################################################
##################   Plot Data   ##################
###################################################

using PyPlot, LaTeXStrings
pygui(true)

fig, axs = subplots()
# fig.suptitle("Information Flow through Quantum FPUT-Chain", fontsize=11)
for p in eachindex(PlotSites)
    axs.plot(periods, array_tracedistances[:,p], c=colorlist[p], label="site $(PlotSites[p])")
end
axs.plot(periods, array_tracedistances[:,1], c=colorlist[1])
axs.grid()
axs.legend(fontsize="xx-large", loc="upper right")
axs.set_xlim([t_start-0.01, t_end+0.01])
axs.set_xlabel(L"time $t \cdot \frac{2\pi}{\sqrt{\kappa}}$", fontsize=fs)
axs.set_ylim([-0.01, 1.01])
axs.set_ylabel(L"$d_\mathrm{tr}(\rho_1,\rho_2)$", fontsize=fs)
axs.tick_params(labelsize=ts)
show()





# fig, axs = subplots(nrows=2)
# fig.suptitle("Dimension Error from Cutoff")
# for p in eachindex(PlotSites)
#     axs[1].plot(periods, array_dimerrors[:,p,1], c=colorlist[p], label="osc. $(PlotSites[p])")
#     axs[1].plot(periods, array_dimerrors[:,p,2], c=colorlist[p], ls="dotted")

#     axs[2].plot(periods, array_cutoffratio[:,p,1], c=colorlist[p], label="osc. $(PlotSites[p])")
#     axs[2].plot(periods, array_cutoffratio[:,p,2], c=colorlist[p], ls="dotted")
# end
# axs[1].grid()
# axs[1].legend()
# axs[1].set_xlim([t_start-0.01, t_end+0.01])
# axs[1].set_title(L"Deviation from trace norm $\mathrm{tr}[\rho]-1$")
# axs[2].grid()
# axs[2].legend()
# axs[2].set_xlim([t_start-0.01, t_end+0.01])
# axs[2].set_xlabel(L"time $t \cdot \frac{2\pi}{\sqrt{\kappa}}$")
# axs[2].set_title("Ratio between larges and outest matrix entry")
# show()


# totalenergies1 = real.(load(datadir("ReducedTrajectories", filename1), "totalenergies"))
# totalenergies2 = real.(load(datadir("ReducedTrajectories", filename2), "totalenergies"))
# groundenergies = groundenergy*ones(length(periods))

# fig, axs = subplots()
# fig.suptitle("Total energy in the chain")
# axs.plot(periods, totalenergies1, label=L"1st. trajectory")
# axs.plot(periods, totalenergies2, label=L"2nd. trajectory")
# axs.plot(periods, groundenergies, label=L"ground state", color="k")
# axs.grid()
# axs.legend()
# axs.set_xlim([t_start-0.01, t_end+0.01])
# axs.set_xlabel(L"time $t \cdot \frac{2\pi}{\sqrt{\kappa}}$")
# axs.set_ylim([-0.01, maximum(totalenergies1)+0.01])
# axs.set_ylabel(L"$E_\mathrm{tot.}$")
# show()