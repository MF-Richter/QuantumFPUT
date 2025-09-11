using DrWatson
@quickactivate "QuantumFPUT"


filename = "Distances(t)-SHORT__N=4_dim=15_α=0.0_β=0.2_state=coherent.jld2"
PlotSites = [1,2,3,4]

# colorlist = ["k","lightcoral","teal","tab:orange","tab:cyan","tab:brown","tab:pink"]
colorlist = ["midnightblue", "lightskyblue", "crimson", "orange"]

fs=25
ts=20


#####################################
##########  Loading File   ##########
#####################################

using JLD2
import QuantumOptics: Ket, Operator, FockBasis
file = load(datadir("InformationFlow", filename))

parameters     = file["parameters"]
periods        = file["periods"]
t_start, t_end = first(periods), 20.0 # last(periods)

initialstates = file["initialstates"]

array_tracedistances = file["array_tracedistances"]
array_wignerdistances = file["array_wignerdistances"]
array_entropy = file["array_entropy"]




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
# W1 = transpose(wigner(reduced(initialstates[1],1), rng, rng))
# W2 = transpose(wigner(reduced(initialstates[2],1), rng, rng))
# Wmax = max(maximum(W1), maximum(W1))

# fig, axs = subplots(nrows=1, ncols=2)
# axs[1].pcolormesh(rng, rng, W1, vmin=-Wmax, vmax=Wmax, cmap="seismic")
# axs[1].set_aspect(1)
# axs[1].grid()
# axs[2].pcolormesh(rng, rng, W2, vmin=-Wmax, vmax=Wmax, cmap="seismic")
# axs[2].set_aspect(1)
# axs[2].grid()
# show()


fig, axs = subplots()
# fig.suptitle("Information Flow through FPUT-Chain")
for p in eachindex(PlotSites)
    site = PlotSites[p]

    if site == 1
        alpha = 1.0
    else
        alpha = 1.0
    end

    if site == 1 || site == N
        ls = "solid"
    else
        ls = "dashed"
    end

    axs.plot(periods, array_tracedistances[:,site], color=colorlist[site], alpha=alpha, ls=ls, label="site $site")
end
axs.plot(periods, array_tracedistances[:,1], color=colorlist[1])
axs.grid()
axs.legend(fontsize="xx-large", loc="upper right")
axs.set_xlim([t_start-0.01, t_end+0.01])
axs.set_xlabel(L"time $t / \tau$", fontsize=fs)
axs.set_ylim([-0.01, 1.01])
axs.set_ylabel(L"$d_\mathrm{tr}(\rho_1,\rho_2)$", fontsize=fs)
show()


# fig, axs = subplots()
# # fig.suptitle("Information Flow through FPUT-Chain")
# for p in eachindex(PlotSites)
#     site = PlotSites[p]
#     axs.plot(periods, array_wignerdistances[:,site], color=colorlist[site], label="site $site")
# end
# axs.grid()
# axs.legend(fontsize="xx-large", loc="upper right")
# axs.set_xlim([t_start-0.01, t_end+0.01])
# axs.set_xlabel(L"time $t \cdot \frac{2\pi}{\sqrt{\kappa}}$", fontsize=fs)
# axs.set_ylim([-0.01, 1.01])
# axs.set_ylabel(L"$d_\mathrm{Kol.}(W_1^{s=0},W_2^{s=0})$", fontsize=fs)
# show()

# fig, axs = subplots()
# for p in eachindex(PlotSites)
#     site = PlotSites[p]
#     axs.plot(periods, array_entropy[:,site], color=colorlist[site], label="site $site")
# end
# axs.grid()
# axs.legend(fontsize="xx-large", loc="upper right")
# axs.set_xlim([t_start-0.01, t_end+0.01])
# axs.set_xlabel(L"time $t \cdot \frac{2\pi}{\sqrt{\kappa}}$", fontsize=fs)
# axs.set_ylim([-0.01, maximum(array_entropy)+0.01])
# axs.set_ylabel(L"$\mathrm{tr}(\rho \ \mathrm{ln}\rho)$", fontsize=fs)
# show()





# fig, axs = subplots()
# plot2, = axs.plot(periods, array_tracedistances[:,2], color=colorlist[2], lw=2, alpha = 0.5, label="site 2")
# plot3, = axs.plot(periods, array_tracedistances[:,3], color=colorlist[3], lw=2, alpha = 0.6, label="site 3")
# plot4, = axs.plot(periods, array_tracedistances[:,4], color=colorlist[4], lw=3, alpha = 1.0, label="site 4")
# plot1, = axs.plot(periods, array_tracedistances[:,1], color=colorlist[1], lw=3, alpha = 1.0, label="site 1")
# axs.grid()
# axs.legend(fontsize="xx-large", loc="upper right", handles=[plot1, plot2, plot3, plot4])
# axs.set_xlim([t_start-0.01, t_end+0.01])
# axs.set_xlabel(L"time $t \cdot \frac{2\pi}{\sqrt{\kappa}}$", fontsize=fs)
# axs.set_ylim([-0.01, 1.01])
# axs.set_ylabel(L"$d_\mathrm{tr}(\rho_1,\rho_2)$", fontsize=fs)
# axs.tick_params(labelsize=ts)
# show()


# fig, axs = subplots()
# plot2, = axs.plot(periods, array_wignerdistances[:,2], color=colorlist[2], lw=2, alpha = 0.5, label=L"$\Gamma_2$")
# plot3, = axs.plot(periods, array_wignerdistances[:,3], color=colorlist[3], lw=2, alpha = 0.6, label=L"$\Gamma_3$")
# plot4, = axs.plot(periods, array_wignerdistances[:,4], color=colorlist[4], lw=3, alpha = 1.0, label=L"$\Gamma_4$")
# plot1, = axs.plot(periods, array_wignerdistances[:,1], color=colorlist[1], lw=3, alpha = 1.0, label=L"$\Gamma_1$")
# axs.grid()
# axs.legend(fontsize="xx-large", loc="upper right", handles=[plot1, plot2, plot3, plot4])
# axs.set_xlim([t_start-0.01, t_end+0.01])
# axs.set_xlabel(L"time $t \cdot \frac{2\pi}{\sqrt{\kappa}}$", fontsize=fs)
# axs.set_ylim([-0.01, 1.01])
# axs.set_ylabel(L"$d_\mathrm{Kol.}(W_1^{s=0},W_2^{s=0})$", fontsize=fs)
# axs.tick_params(labelsize=ts)
# show()

# fig, axs = subplots()
# plot2, = axs.plot(periods, array_entropy[:,2], color=colorlist[2], lw=1, label="site 2")
# plot3, = axs.plot(periods, array_entropy[:,3], color=colorlist[3], lw=1, label="site 3")
# plot4, = # axs.plot(periods, array_entropy[:,4], color=colorlist[4], lw=1, label="site 4")
# plot1, = axs.plot(periods, array_entropy[:,1], color=colorlist[1], lw=1, label="site 1", ls="dotted")
# axs.grid()
# # axs.legend(fontsize="xx-large", loc="upper right", handles=[plot1, plot2, plot3, plot4])
# axs.legend(fontsize="xx-large", loc="lower right", handles=[plot1, plot2, plot3])
# axs.set_xlim([t_start-0.01, t_end+0.01])
# axs.set_xlabel(L"time $t \cdot \frac{2\pi}{\sqrt{\kappa}}$", fontsize=fs)
# axs.set_ylim([-0.01, maximum(array_entropy)+0.01])
# axs.set_ylabel(L"$\mathrm{tr}(\rho \ \mathrm{ln}\rho)$", fontsize=fs)
# axs.tick_params(labelsize=ts)
# show()