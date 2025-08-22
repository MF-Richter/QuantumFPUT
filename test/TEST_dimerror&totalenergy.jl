using DrWatson
@quickactivate "QuantumFPUT"

filename = "Dbl__HPC-LONG__N=4_dim=15_α=0.0_β=0.2_state=coherent.jld2"
# folder = "EnsembleTrajectories"
folder = "KetTrajectories"

colorlist = ["k","lightcoral","teal","tab:orange","tab:cyan","tab:brown","tab:pink"]

fs=15
ts=10



#####################################
##########  Loading File   ##########
#####################################

using JLD2
import QuantumOptics: FockBasis, Ket, Operator, expect

parameters   = load(datadir(folder, filename), "parameters")
periods      = load(datadir(folder, filename), "periods")
states = load(datadir(folder, filename), "states1")




##########################################
##########   Print Parameters   ##########
##########################################


N = parameters["N"]
dim = parameters["dim"]
n_mcwf = get(parameters,"n_mcwf", nothing)

κ = parameters["κ"]
α = parameters["α"]
β = parameters["β"]

γ   = get(parameters,"γ", nothing)
kbT = get(parameters,"kbT", nothing)
BathSites = get(parameters, "BathSites", nothing)

abstol = get(parameters, "abstol", "1e-8")
reltol = get(parameters, "reltol", "1e-6")
alg = get(parameters, "alg", "OrdinaryDiffEq.DP5()")

t_start, t_end = first(periods), last(periods)


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

using HamiltonianFPUT
using StateCharacteristics
using ReduceSystems

Ψ0 = states[1]
single_basis = reduce_to(Ψ0,N,1).basis_l

Ham = HamOpFPUT(N, κ,α,β, single_basis)
energy(ψ::Ket) = expect(Ham, ψ)

cutoffvalue = Array{Float64,2}(undef, N,length(periods))
cutoffratio = Array{Float64,2}(undef, N,length(periods))
for i in 1:N
    @time cutoffvalue[i,:] .= tracenorm_diviation.(reduce_to(states, N, i))
    @time cutoffratio[i,:] .= relative_dimcutoffvalue.(reduce_to(states, N, i))
end

@time energies = energy.(states)
energydiviation = energies[1]*ones(length(energies)) - energies
relative_energydiviation = ( energies[1]*ones(length(energies)) - energies )/energies[1]

using PyPlot
pygui(true)

fig, axs = subplots(nrows=3,ncols=1)
# fig.suptitle(L"Dimension Errors $\Delta = | 1 - \|\rho\|_\mathrm{tr} |$ and Energy Diviation $\epsilon = \frac{\langle \hat{H}_\mathrm{FPUT} \rangle(t_0) - \langle \hat{H}_\mathrm{FPUT} \rangle(t)}{\langle \hat{H}_\mathrm{FPUT} \rangle(t_0)}$")

for i in 1:N
    axs[1].plot(periods, cutoffvalue[i,:], ".", color=colorlist[i], label="site $i")
    axs[2].plot(periods, cutoffratio[i,:], ".", color=colorlist[i], label="site $i")
end

# axs[1].set_xlabel(L"time $t \cdot \frac{2\pi}{\sqrt{\kappa}}$", fontsize=fs)
axs[1].set_xlim([t_start-0.01, t_end+0.01])
axs[1].set_ylabel(L"$\Delta$", fontsize=fs)
axs[1].legend(fontsize="xx-large", loc="upper right")
axs[1].tick_params(labelsize=ts)
axs[1].grid()

# axs[2].set_xlabel(L"time $t \cdot \frac{2\pi}{\sqrt{\kappa}}$", fontsize=fs)
axs[2].set_xlim([t_start-0.01, t_end+0.01])
axs[2].set_ylabel(L"$\delta$", fontsize=fs)
# axs[2].legend(fontsize="xx-large", loc="upper right")
axs[2].tick_params(labelsize=ts)
axs[2].grid()

axs[3].plot(periods, relative_energydiviation, ".", color="midnightblue")
axs[3].set_xlim([t_start-0.01, t_end+0.01])
axs[3].set_xlabel(L"time $t \cdot \frac{2\pi}{\sqrt{\kappa}}$", fontsize=fs)
axs[3].set_ylabel(L"$\epsilon$", fontsize=fs)
axs[3].tick_params(labelsize=ts)
axs[3].grid()

show()