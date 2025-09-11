using DrWatson
@quickactivate "QuantumFPUT"


filename = "Distances(t)-LONG__N=4_dim=15_α=0.0_β=0.0_state=coherent.jld2"
folder = "InformationFlow"


#####################################
##########  Loading File   ##########
#####################################

using JLD2
import QuantumOptics: Ket, Operator, FockBasis
file = load(datadir(folder, filename))

println("The file has the following keys: \n", keys(file))

println(file["parameters"])