"""
With this script, one can compute two ensembles of Monte-Carlo wave function trajectories,
initially starting with two oppositely displaced coherent states, for a FPUT chain with
a thermal bath of a given temperature coupled to a given site of the chain. The two ensembles
are saved in data/EnsembleTrajectories for further processing.
"""

using DrWatson
@quickactivate "QuantumFPUT"

# Complex phase of initial states
r = 1.5
θ = 0.5*pi
φ1 = r * (cos(θ) + sin(θ)*im)/sqrt(2)
φ2 = -φ1

## Coupling parameters
# κ... general coupling of interaction Hamiltonian
# α... strength of S³-term in FPUT interaction
# β... strength of S⁴-term in FPUT interaction
κ = 1.0
α = 0.0
β = 0.0


## Bath parameters
# kbT... temperature of the bath
# γ...   coupling strength of the bath
# BathSites... sites where to attache the bath to
kbT = 0.0
γ   = 0.25
BathSites = [2]

## Parameters of the chain and numerical approximation
# Number of chain elements and dim. cutoff of single Hilbert spaces
N = 2
dim = 15
n_mcwf = 100

## Parameters of time
# t_0... initial time measured in single oscillator periodes
# t_end... final time measured in single oscillator periodes
# n_time... number of time steps
t_0 = 0.0
t_end =  10.0
n_time = 200

# optinal prefix for the name of the data file
Prefix = "TEST"







######################################################
##################   Compute Data   ##################
######################################################

using HamiltonianFPUT
import QuantumOptics: Ket, FockBasis,coherentstate, destroy, ⊗
using EnsembleMCWF

println("\n\n\n")
println("##################### ###################### #####################")
println("#####################  Computations started  #####################")

basis = FockBasis(dim)
times = LinRange(t_0, t_end*2*pi, n_time)
periods = LinRange(t_0, t_end, n_time)

ket1 = coherentstate(basis, φ1)
ket2 = coherentstate(basis, φ2)

ketFPUT1 = localKets(ket1, N)
ketFPUT2 = localKets(ket2, N)
Hamiltonian = HamOpFPUT(N, κ,α,β, basis)
Jops = JumpOperators(destroy(basis), N, BathSites)
Jrates = JumpRates(kbT,γ,BathSites)
println("\ninitial ket state and Hamiltonian and normal mode energies built")

using OrdinaryDiffEq
alg = OrdinaryDiffEq.DP8()
abstol = 1e-8
reltol = 1e-6
timesout, states1 = @time ensembleevolution_mcwf(times,n_mcwf, ketFPUT1,Hamiltonian,Jops; rates=Jrates, abstol=abstol,reltol=reltol, alg=alg)
timesout, states2 = @time ensembleevolution_mcwf(times,n_mcwf, ketFPUT2,Hamiltonian,Jops; rates=Jrates, abstol=abstol,reltol=reltol, alg=alg)

println("time evolution of ket ensemble computed")



#############################################
##########   Saving Trajectories   ##########
#############################################

using JLD2

parameters = @strdict N dim n_mcwf κ α β kbT γ BathSites state="coherent" abstol reltol alg
data = @strdict parameters periods states1 states2
filename = datadir("EnsembleTrajectories", savename("Dbl__"*Prefix, parameters, "jld2", accesses=["N","dim","n_mcwf", "α","β", "kbT","γ", "state"], sort=false))

safesave(filename, data)

println("states saved\n")

println("##################### Computations completed #####################")
println("##################### ###################### #####################\n\n")