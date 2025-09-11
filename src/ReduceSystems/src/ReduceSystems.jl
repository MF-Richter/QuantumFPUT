module ReduceSystems

    import QuantumOptics: Ket, Operator, dagger, tr, tracenorm, ptrace, reduced
    using Base.Threads
    export reduce_to, reduce_partially, correct_dm


    """
    This function reduces the ket vector of a multipartite system to the reduced states of given subsystems
    """  
    function reduce_to(
        ket::Ket,  # ket state of whole chain
        N::Int64,  # number of chain elements
        subsystem::Int64  # index of chain elements to which the total state gets reduced
        )
        ρ = ptrace(ket, deleteat!(collect(1:N),subsystem))
        return ρ
    end    

    function reduce_to(
        ket::Ket,  # ket state of whole chain
        N::Int64,  # number of chain elements (or subsystems)
        subsystems::Vector{Int64}  # indices of chain elements to which the total state gets reduced
        )

        reducedstates = Vector{Operator}()
        for i in subsystems
            ρ = reduce_to(ket, N, i)
            push!(reducedstates, ρ)
        end
        return reducedstates
    end
 
    function reduce_to(
        ket_trajectory::Vector{<:Ket},
        N::Int64,
        subsystem::Int64
        )
        dms = Vector{Operator}(undef, length(ket_trajectory))
        for t in 1:length(ket_trajectory)
            dms[t] = reduce_to(ket_trajectory[t],N,subsystem)
        end

        return dms
    end

    function reduce_to(
        ket_trajectory::Vector{<:Ket},
        N::Int64,
        subsystems::Vector{Int64}
        )
        n_times = length(ket_trajectory)
        n_subsys = length(subsystems)
        array_dms = Array{Operator, 2}(undef, n_times,n_subsys)

        for t in 1:n_times, i in 1:n_subsys
            array_dms[t,i] = reduce_to(ket_trajectory[t],N,subsystems[i])
        end

        return array_dms
    end

    function reduce_to(
        ensemble_trajectories::Array{Ket, 2},
        N::Int64,
        subsystem::Int64)

        n_mcwf = size(ensemble_trajectories)[2]
        array_dms = 1/n_mcwf * reduce_to(ensemble_trajectories[:,1], N, subsystem)

        lock = ReentrantLock()
        @threads for i in 2:n_mcwf
            dms = 1/n_mcwf * reduce_to(ensemble_trajectories[:,i], N, subsystem) 
            @lock lock array_dms += dms
        end
        
        return array_dms
    end
   
    function reduce_to(
        ensemble_trajectories::Array{Ket, 2},
        N::Int64,
        subsystems::Vector{Int64}
        )

        n_mcwf = size(ensemble_trajectories)[2]
        array_dms = 1/n_mcwf * reduce_to(ensemble_trajectories[:,1], N, subsystems)

        lock = ReentrantLock()
        @threads for i in 2:n_mcwf
            dms = 1/n_mcwf * reduce_to(ensemble_trajectories[:,i], N, subsystems) 
            @lock lock array_dms += dms 
        end
        
        return array_dms
    end


    """
    This function reduces a sequence of kets - i.e. through time - to the density matrices for the combined subsystem given by 'subsystems'.
    Note that unlike kets_to_partialdms() the function reduce_kets_partially returns for each ket a single dm representing the state of the
    COMBINED (thus only partially reduced) subsystem and not several dms for each atomic subsystem.
    """
    function reduce_partially(
        ket_trajectory::Vector{<:Ket},
        subsystems::Vector{Int64}
        )

        array_dms = Vector{Operator}(undef, length(ket_trajectory))
        for t in 1:length(ket_trajectory)
            array_dms[t] = reduced(ket_trajectory[t],subsystems)
        end
        return array_dms
    end

    function reduce_partially(
        ensemble_trajectories::Array{Ket, 2},
        subsystems::Vector{Int64}
        )
        n_mcwf = size(ensemble_trajectories)[2]
        array_dms = 1/n_mcwf * reduce_partially(ensemble_trajectories[:,1], subsystems)

        lock = ReentrantLock()
        @threads for i in 2:n_mcwf
            dms = 1/n_mcwf * reduce_partially(ensemble_trajectories[:,i], subsystems) 
            @lock lock array_dms += dms
        end
        return array_dms
    end


    """
    Correctring supposed density matrices to be hermitian and of trace one.
    """    
    function correct_dm(dm::Operator)
        ρ = 0.5*(dm + dagger(dm))
        return ρ/tr(ρ)
    end

end # module ReduceSystems
