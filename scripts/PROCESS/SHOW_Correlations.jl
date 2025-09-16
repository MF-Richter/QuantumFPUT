using DrWatson
@quickactivate "QuantumFPUT"


filename = "Dbl___N=2_dim=15_α=0.0_β=1.0_state=coherent.jld2"
# folder = "EnsembleTrajectories"
folder = "KetTrajectories"

CorrSites = (1,2)



#####################################
##########  Loading File   ##########
#####################################

using JLD2
import QuantumOptics: Ket, Operator, FockBasis

parameters  = load(datadir(folder, filename), "parameters")
periods     = load(datadir(folder, filename), "periods")
states      = load(datadir(folder, filename), "states1")

t_start, t_end = first(periods), last(periods)




##########################################
##########   Print Parameters   ##########
##########################################
using JLD2

N = parameters["N"]
dim = parameters["dim"]
n_mcwf = get(parameters,"n_mcwf", 1)

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

using Correlations

totalcorrelations, classcorrelations, quantcorrelations = @time correlations(states, CorrSites)



###################################################
##################   Plot Data   ##################
###################################################

using PyPlot
pygui(true)

fig, axs = subplots()
axs.plot(periods, totalcorrelations, label="full corr. ρ$(CorrSites[1]) and ρ$(CorrSites[2])")
axs.plot(periods, classcorrelations, label="class. corr. ρ$(CorrSites[1]) and ρ$(CorrSites[2])")
axs.plot(periods, quantcorrelations, label="quant. corr. ρ$(CorrSites[1]) and ρ$(CorrSites[2])")
axs.grid()
axs.legend()
axs.set_ylim([-0.01, 1.01])
axs.set_xlim([t_start-0.01, t_end+0.01])
show()