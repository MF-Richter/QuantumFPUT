"""
This module contains only one function to compute an ensemble of Monte-Carlo trajectories
to approximate the Markovian time evolution of a high-dimensional system coupled to a Markovian
bath. It takes a vector 'times' with the time steps to display the states, a number
'n_mcwf' of trajectories to be averaged over later, as well as the system's Hamiltonian 'Ham'
and the jump operators, and returns the time list and the ensemble of trajectories.

It is essentially a wrapper around timeevolution.mcwf() of the QuantumOptics package for better
use in this project.
"""
module EnsembleMCWF

    import QuantumOptics: Ket,Operator, timeevolution, dagger, expect, reduced
    using Base.Threads

    export ensembleevolution_mcwf

    function ensembleevolution_mcwf(
        times::Any,     # vector of time steps to display
        n_mcwf::Int64,  # number of trajectories
        ket0::Ket,      # initial ket state of the system
        Ham::Operator,  # system Hamiltonian
        J;              # jump operators with the bath
        kwargs...
        )

        n_times = length(times)
        ensemble_trajectories = Array{Ket, 2}(undef, n_times,n_mcwf)
        timesout = times

        lock = ReentrantLock()
        @threads for i in 1:n_mcwf  # loop for n_mcwf trajectories, stored in 'ensemble_times'
            timesout, trajectory = timeevolution.mcwf(times,ket0,Ham,J; kwargs...)
            @lock lock ensemble_trajectories[:,i] = trajectory
        end

        return timesout, ensemble_trajectories
    end

end
