using DrWatson
@quickactivate "QuantumFPUT"

filename1 = "HPC__N=3_dim=15_n_mcwf=1000_α=0.0_β=0.2_kbT=0.0_γ=0.1_state=coh+.jld2"
filename2 = "HPC__N=3_dim=15_n_mcwf=1000_α=0.0_β=0.2_kbT=0.0_γ=0.1_state=coh-.jld2"
PlotSites = [1,2,3]
qrange = LinRange(-5, 5, 100)
prange = LinRange(-5, 5, 100)

# Time or index of plotting
# T = 210
# T = 500
T = 100.0



#####################################
##########  Loading File   ##########
#####################################

using JLD2
import QuantumOptics: Ket, Operator, FockBasis

parameters = load(datadir("ReducedTrajectories", filename1), "parameters")
periods    = load(datadir("ReducedTrajectories", filename1), "periods")
array_dms = load(datadir("ReducedTrajectories", filename1), "array_dms")
# array_dms = load(datadir("ReducedTrajectories", filename2), "array_dms")

t_start, t_end = first(periods), last(periods)
n_time = length(periods)


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

if γ === nothing
    println("No Markovian baths attached\n")
else
    println("Baths temperature is kbT=$kbT")
    println("Baths sites are $BathSites")
    println("Baths coupling strength is γ=$γ")
    println("Averaged over $n_mcwf MC trajectories\n")
end

# println("Solved with $alg at abstol=$abstol and reltol=$reltol")

println("##################### ################ #####################")
println("##################### ################ #####################")
println("\n")


#############################################################
##################   Compute Partial DMs   ##################
#############################################################
import QuantumOptics: wigner

if typeof(T)==Int64
    τ = T
elseif typeof(T)==Float64
    τ = round(Int, T/t_end * n_time)
else
    error("No valid data type for time or time index")
end
t = round(periods[τ], digits=2)

Ws = Vector{Matrix{Float64}}(undef, length(PlotSites))
vmax = 0.0
@time for p in eachindex(PlotSites)
    Ws[p] = wigner(array_dms[τ,p], qrange,prange)
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
    τ = round(Int, T/t_end * n_time)
else
    error("No valid data type for time or time index")
end
t = round(periods[τ], digits=2)

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
    # ax.grid()
end

show()