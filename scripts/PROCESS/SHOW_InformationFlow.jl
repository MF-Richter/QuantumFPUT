"""
This script loads a double state trajectory from the data folder and computes, for given
sites of the chain, the trace distance and the Kolmogorov distance between the Wigner functions
of the two trajectories. It also computes the von-Neumann entropy for the first one. The distances
and the entropy are plotted and saved as processed data in the folder:
    data/InformationFlow
"""

using DrWatson
@quickactivate "QuantumFPUT"
   

filename = "Dbl___N=2_dim=15_α=0.0_β=0.0_state=coherent.jld2"
folder = "KetTrajectories"
PlotSites = [1,2]

Prefix = ""


#####################################
##########  Loading File   ##########
#####################################

using JLD2
import QuantumOptics: Ket, Operator, FockBasis

parameters     = load(datadir(folder, filename), "parameters")
periods        = load(datadir(folder, filename), "periods")
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

abstol = get(parameters, "abstol", "1e-8")
reltol = get(parameters, "reltol", "1e-6")
alg = get(parameters, "alg", "OrdinaryDiffEq.DP5()")


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

using ReduceSystems
import QuantumOptics: tracedistance, tr, entropy_vn

states = load(datadir(folder, filename), "states1")
ket_init_1 = states[1,1]
@time array_dms1 = reduce_to(states, N, PlotSites)

states = load(datadir(folder, filename), "states2")
ket_init_2 = states[1,1]
@time array_dms2 = reduce_to(states, N, PlotSites)

array_tracedistances  = Array{Float64,2}(undef, length(periods), length(PlotSites))
array_wignerdistances = Array{Float64,2}(undef, length(periods), length(PlotSites))
array_entropy   = Array{Float64,2}(undef, length(periods), length(PlotSites))

import QuantumOptics: wigner
"""
Kolmogorov distance between the Wigner functions of density operators 'ρ1' and 'ρ2', given over a phase-space grid
defined by position range 'qrange' and momentum range 'prange'.
"""
function wignerdistance(ρ1::Operator, ρ2::Operator; qrange=LinRange(-5.0,5.0,100),prange=LinRange(-5.0,5.0,100))
    W1 = wigner(ρ1, qrange,prange)
    W2 = wigner(ρ2, qrange,prange)
    ΔW = abs.(W1-W2)
    dx = (maximum(qrange)-minimum(qrange))/length(qrange) * (maximum(prange)-minimum(prange))/length(prange)
    return 0.5*sum(ΔW)*dx
end


@time for p in eachindex(PlotSites)
    dms1 = array_dms1[:,p]
    dms2 = array_dms2[:,p]
    array_tracedistances[:,p]  .=  tracedistance.(dms1,dms2)
    array_wignerdistances[:,p] .= wignerdistance.(dms1,dms2)
    array_entropy[:,p]   .= entropy_vn.(dms1)
end




###################################################
##################   Plot Data   ##################
###################################################

using PyPlot, LaTeXStrings
pygui(true)
colors = ["tab:blue","tab:orange","tab:green","tab:purple","tab:cyan","tab:brown","tab:pink"]

fig, axs = subplots()
fig.suptitle("Information Flow through FPUT-Chain", fontsize=11)
for p in eachindex(PlotSites)
    axs.plot(periods, array_tracedistances[:,p], color=colors[p], label="osc. $(PlotSites[p])")
end
axs.grid()
axs.legend()
axs.set_xlim([t_start-0.01, t_end+0.01])
axs.set_xlabel(L"time $t \cdot \frac{2\pi}{\sqrt{\kappa}}$")
axs.set_ylim([-0.01, 1.01])
axs.set_ylabel(L"$d_\mathrm{tr}(\rho_1,\rho_2)$")
axs.set_title("Coupling Potential: V(r)= r²*$κ/2! + r³*$α/3! + r⁴*$β/4!", fontsize=9)
show()


fig, axs = subplots()
fig.suptitle("Information Flow through FPUT-Chain", fontsize=11)
for p in eachindex(PlotSites)
    axs.plot(periods, array_wignerdistances[:,p], color=colors[p], label="osc. $(PlotSites[p])")
end
axs.grid()
axs.legend()
axs.set_xlim([t_start-0.01, t_end+0.01])
axs.set_xlabel(L"time $t \cdot \frac{2\pi}{\sqrt{\kappa}}$")
axs.set_ylim([-0.01, 1.01])
axs.set_ylabel(L"$d_\mathrm{Kol.}(W_1,W_2)$")
axs.set_title("Coupling Potential: V(r)= r²*$κ/2! + r³*$α/3! + r⁴*$β/4!", fontsize=9)
show()


fig, axs = subplots()
fig.suptitle("Entropy in reduced systems")
for p in eachindex(PlotSites)
    axs.plot(periods, array_entropy[:,p], color=colors[p], label="osc. $(PlotSites[p])")
end
axs.grid()
axs.legend()
axs.set_xlim([t_start-0.01, t_end+0.01])
axs.set_xlabel(L"time $t \cdot \frac{2\pi}{\sqrt{\kappa}}$")
axs.set_ylim([-0.01, maximum(array_entropy)+0.01])
axs.set_ylabel(L"$\mathrm{tr}(\rho \mathrm{ln}\rho)$")
show()



###################################
##########   Saving Data ########## 
###################################

initialstates = [ket_init_1, ket_init_2]
data = @strdict initialstates parameters periods array_tracedistances array_wignerdistances array_entropy

if γ===nothing
    filename = datadir("InformationFlow", savename("Distances(t)"*Prefix, parameters, "jld2", accesses=["N","dim", "α","β", "state"], sort=false))
else
    filename = datadir("InformationFlow", savename("Distances(t)"*Prefix, parameters, "jld2", accesses=["N","dim","n_mcwf", "α","β", "kbT","γ", "state"], sort=false))
end
safesave(filename, data)
println("data saved\n")