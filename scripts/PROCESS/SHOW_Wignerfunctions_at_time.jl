"""
This script loads a state trajectory and plots, for a given time (if 'T' is a float) or time index (if 'T' is an integer),
the Wigner functions at this point in time over the local phase-spaces set by 'PlotSites'.
"""

using DrWatson
@quickactivate "QuantumFPUT"


filename = "Dbl___N=2_dim=15_α=0.0_β=1.0_state=coherent.jld2"
folder = "KetTrajectories"
PlotSites = [1,2]

qrange = LinRange(-5, 5, 250)
prange = LinRange(-5, 5, 250)

# Time or index of plotting
T = 5.0



#####################################
##########  Loading File   ##########
#####################################

using JLD2
import QuantumOptics: Ket, Operator, FockBasis

parameters  = load(datadir(folder, filename), "parameters")
periods     = load(datadir(folder, filename), "periods")
states      = load(datadir(folder, filename), "states1")

t_start, t_end = first(periods), last(periods)




##########################################################
##################   Print Parameters   ##################
##########################################################
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

println("##################### ################ #####################")
println("##################### ################ #####################")
println("\n")


#############################################################
##################   Compute Partial DMs   ##################
#############################################################
using ReduceSystems
import QuantumOptics: wigner

if typeof(T)==Int64
    τ = T
elseif typeof(T)==Float64
    τ = round(Int, T/t_end * length(periods))
else
    error("No valid data type for time or time index")
end
t = round(periods[τ], digits=2)

dms  = @time reduce_to(states[τ],N,PlotSites)

Ws = Vector{Matrix{Float64}}(undef, length(PlotSites))
vmax = 0.0
@time for p in eachindex(PlotSites)
    Ws[p] = wigner(dms[p], qrange,prange)
    if vmax < maximum(Ws[p])
        global vmax = maximum(Ws[p])
    end
end


###################################################
##################   Plot Data   ##################
###################################################
using PyPlot
pygui(true)


if typeof(T)==Int64
    τ = T
elseif typeof(T)==Float64
    τ = round(Int, T/t_end * length(periods))
else
    error("No valid data type for time or time index")
end
t = round(periods[τ], digits=2)

array_dms = @time reduce_to(states[τ],N,PlotSites)

fig, axs = subplots(ncols=length(PlotSites))

t = round(periods[τ], digits=2)
fig.suptitle("Quantum Wigner fct. at t = $t "*L"$2\pi/\sqrt{\kappa}$")

for p in eachindex(PlotSites)

    site = PlotSites[p]

    if length(PlotSites)==1
        ax = axs
    else
        ax = axs[p]
    end

    ax.pcolormesh(qrange, prange, transpose(Ws[p]), cmap="seismic", aa=true, vmin=-vmax, vmax=vmax)
    ax.set_xlabel(L"$q$", fontsize=12)
    if p==1
        ax.set_ylabel(L"$p$", fontsize=12)
    end
    ax.set_title("Site $site")
    ax.set_aspect(1)
    ax.tick_params(labelsize=10)
    ax.grid()
end

show()