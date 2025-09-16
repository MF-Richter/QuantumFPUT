"""
This script loads some state trajectories from the data folder and computes the mean values and
covariance matrices for the chain's normal mode quadratures and plots them.
"""

using DrWatson
@quickactivate "QuantumFPUT"

filename = "Dbl___N=2_dim=15_α=0.0_β=0.0_state=coherent.jld2"
folder   = "KetTrajectories"



#####################################
##########  Loading File   ##########
#####################################

using JLD2
import QuantumOptics: Ket, Operator, FockBasis
file = load(datadir(folder, filename))

parameters     = file["parameters"]
periods        = file["periods"]
t_start, t_end = first(periods), last(periods)

states = file["states1"]




##########################################
##########   Print Parameters   ##########
##########################################
using JLD2

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

import NormalModes: QuadraturesNormalMode, NumberoperatorNormalMode, means
function varμ(varQQ, varPP, varQP)
    return sqrt(varQQ*varPP - varQP^2)
end


Qs = Vector{Operator}(undef, N)
Ps = Vector{Operator}(undef, N)

QQs = Vector{Operator}(undef, N)
PPs = Vector{Operator}(undef, N)
QPs = Vector{Operator}(undef, N)

@time for mode in 1:N
    Position,Momentum = QuadraturesNormalMode(N,mode,FockBasis(dim))

    Qs[mode] = Position
    Ps[mode] = Momentum

    QQs[mode] = 2 * Position*Position
    PPs[mode] = 2 * Momentum*Momentum
    QPs[mode] = Position*Momentum + Momentum*Position
end

array_QmeansNM = means(Qs,states)
array_PmeansNM = means(Ps,states)

array_μsNM = sqrt.(array_QmeansNM.^2 + array_PmeansNM.^2)

array_QQvarsNM = means(QQs, states) - 2*array_QmeansNM.^2
array_PPvarsNM = means(PPs, states) - 2*array_PmeansNM.^2
array_QPvarsNM = means(QPs, states) - 2*array_QmeansNM.*array_PmeansNM

array_detCovsNM = varμ.(array_QQvarsNM, array_PPvarsNM, array_QPvarsNM)




###################################################
##################   Plot Data   ##################
###################################################

using PyPlot, LaTeXStrings
pygui(true)


fig, axs = subplots(nrows=2,ncols=1)
fig.suptitle("Normal Modes: Phase-Space Means and (Co)variances")

for mode in 1:N
    axs[1].plot(periods, array_QmeansNM[:,mode],  label="normal mode $mode")
    axs[2].plot(periods, array_detCovsNM[:,mode], label="normal mode $mode")
end

axs[1].grid()
axs[1].legend()
axs[1].set_xlim([t_start-0.01, t_end+0.01])
axs[1].set_ylabel(L"$|\vec{\mu}|^2$")

axs[2].grid()
axs[2].legend()
axs[2].set_xlim([t_start-0.01, t_end+0.01])
axs[2].set_ylim([0.0, maximum(array_detCovsNM)+0.01])
axs[2].set_xlabel(L"time $t \cdot \frac{2\pi}{\sqrt{\kappa}}$")
axs[2].set_ylabel(L"$\sqrt{\mathrm{det} \Sigma}$")

show()