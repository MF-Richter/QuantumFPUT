module EnsembleMCWF

    import QuantumOptics: Ket,Operator, timeevolution, dagger, expect, reduced
    using Base.Threads

    export ensembleevolution_mcwf


    """
    This function computes an ensemble of Monte-Carlo trajectories to approximate the Markovian time evolution
    of a high dimesnional system coupled to a Markovian bath. It takes a vector 'times' with the time steps to
    display the states, a number 'n_mcwf' of trajectories to be averaged over later as well as the systems Ha-
    miltonian 'Ham' and the jump operators and returns the time list and the ensemble of trajectories.
    """
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

        # for i in 1:n_mcwf  # loop for n_mcwf trajectories, stored in 'ensemble_times'
        #     timesout, trajectory = timeevolution.mcwf(times,ket0,Ham,J; kwargs...)
        #     ensemble_trajectories[:,i] = trajectory
        # end

        return timesout, ensemble_trajectories
    end


    # function reducedevolution_mcwf(
    #     times::Any,     # vector of time steps to display
    #     n_mcwf::Int64,  # number of trajectories
    #     ket0::Ket,      # initial ket state of the system
    #     Ham::Operator,  # system Hamiltonian
    #     J;              # jump operators with the bath
    #     kwargs...
    #     )

    #     # computing single subsystem dimension and number of subsystems from chain basis object
    #     basis_chain = ket0.basis
    #     dim_chain = length(basis_chain)
    #     dim_subsys = length(reduced(basis_chain,1))
    #     N = Int64(log(dim_subsys,dim_chain))

    #     n_times = length(times)
        
    #     # prepare Array for reduced states and total energy function
    #     reducedstates = Array{Operator, 2}(undef, N, n_times)
    #     totalenergy(ket::Ket) = expect(Ham, ket)

    #     timesout, trajectory = timeevolution.mcwf(times,ket0,Ham,J; kwargs...)
    #     totalenergies = totalenergy.(trajectory)
    #     for n in 1:N
            
    #         reducedstates[n,:] .= reduced.(trajectory,n)
    #     end

    #     totalenergies = Vector{Float64}(undef, n_times)
    #     reducedtrajectories = Array{Operator, 2}(undef, N, n_times)

    #     lock = ReentrantLock()
    #     @threads for i in 1:n_mcwf  # loop for n_mcwf trajectories, stored in 'ensemble_times'
    #         timesout, trajectory = timeevolution.mcwf(times,ket0,Ham,J; kwargs...)
    #         @lock lock ensemble_trajectories[:,i] = trajectory
    #     end

    #     # for i in 1:n_mcwf  # loop for n_mcwf trajectories, stored in 'ensemble_times'
    #     #     timesout, trajectory = timeevolution.mcwf(times,ket0,Ham,J; kwargs...)
    #     #     ensemble_trajectories[:,i] = trajectory
    #     # end

    #     return timesout, ensemble_trajectories
    # end

end
