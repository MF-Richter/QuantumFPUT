
"""
This module contains functions to compute, for a given state of a multi-partite quantum system,
the correlations between two of its subsystems. The correlation is measured by the trace
distance between the fully correlated state and the product of its reduced states (decorrelated
state).

It also contains a function to compute the quantum and classical correlations. To do so, the
correlations of the given state are reduced to its classical ones only using a classification
procedure. The quantum correlations are then measured by the distance between the fully cor-
related state and the classified state.

The classical correlations are accordingly measured by the distance between the classified state and
the decorrelated state.
"""
module Correlations

    import QuantumOptics: Operator, dense, Ket, reduced, dm, ⊗, tracedistance, identityoperator, dagger, eigenstates
    using  Base.Threads
    export totalcorrelation, classicalcorrelation, quantumcorrelation, correlations


    """
    Correlations within a bipartite state measured by the trace distance between the fully correlated
    state and its decorrelation
    """    
    function totalcorrelation(ρ::Operator)
        return tracedistance(ρ, reduced(ρ,1)⊗reduced(ρ,2))
    end

    """
    The classification procedure: A function to reduce the correlations of a bipartite state 'ρ' to its
    classical correlations only.
    """    
    function classify_correlations(ρ::Operator)
        # separate state into its partial states
        ρA = reduced(ρ, 1)
        ρB = reduced(ρ, 2)
        dimA = length(ρA.basis_l)
        idB = identityoperator(ρB.basis_l)  

        ρA_hermitian = (ρA + dagger(ρA))/2
        eigenstatesA = eigenstates(ρA_hermitian)[2]

        Π = dm(eigenstatesA[1])                     # projector on 1st. eigenstate of ρA ⊗ the identity on B
        ρ_ClassCorr = Π⊗idB * ρ * dagger(Π⊗idB)
        
        # loop over the remaining eigenstates of ρA
        for i in 2:dimA
            Π = dm(eigenstatesA[i])
            ρ_ClassCorr += Π⊗idB * ρ * dagger(Π⊗idB)
        end

        return ρ_ClassCorr
    end


    """
    Quantum Correlation measured by the distance between the classified version of the state 'ρ' and its
    fully correlated version
    """    
    function quantumcorrelation(ρ::Operator)
        return tracedistance(ρ, classify_correlations(ρ))
    end


    """
    Classical correlations measured by the distance between the classified version of the boipartite state
    'ρ' and its fully decorrelated version
    """    
    function classicalcorrelation(ρ::Operator)
        return tracedistance( classify_correlations(ρ), reduced(ρ, 1)⊗reduced(ρ, 2) )
    end


    """
    A function that takes a trajectory of ket states through time (or a Monte-Carlo ensemble of them) and a tuple
    containing the indices of the two subsystems between which the correlations should be measured.
    It returns three vectors containing the total correlation, the classical correlation, and the quantum correlation
    through time between these two subsystems.

    If just a ket trajectory is passed, the additional kwarg 'onlytotalcorrelation' (default is 'true') can be used
    to only compute the total correlation. Since ket states cannot possess correlations different from entanglement, the
    distinction between classical and quantum correlation becomes irrelevant - every correlation is of genuine quantum
    origin (however, entanglement can trigger classically looking correlations too, although their origin must be
    quantum)
    """
    function correlations(
        trajectory::Vector{<:Ket},
        subsystems::Tuple{Int64,Int64};
        only_total_correlation::Bool = false
        )

        n_times = length(trajectory)  # number of time steps

        totalcorrelations = Vector{Float64}(undef, n_times)

        if only_total_correlation == false
            classcorrelations = Vector{Float64}(undef, n_times)
            quantcorrelations = Vector{Float64}(undef, n_times)
        end

        # ReentrantLock for multi-threading to prevent data race
        lock1 = ReentrantLock()
        lock2 = ReentrantLock()
        lock3 = ReentrantLock()

        @threads for τ in 1:n_times  # mutli-threaded loop over time steps

            # reducing the ket vector to a bipartite desnity matrix of the two subsystems
            dm = reduced(trajectory[τ], collect(subsystems))

            totalcorr_dm = totalcorrelation(dm)
            @lock lock1 totalcorrelations[τ] = totalcorr_dm

            if only_total_correlation == false
                classcorr_dm = classicalcorrelation(dm)
                quantcorr_dm = quantumcorrelation(dm)

                @lock lock2 classcorrelations[τ] = classcorr_dm
                @lock lock3 quantcorrelations[τ] = quantcorr_dm
            end
        end
        
        if only_total_correlation == true
            return totalcorrelations
        else
            return totalcorrelations, classcorrelations, quantcorrelations
        end
    end

    function correlations(
        ensemble_trajectories::Array{Ket, 2},
        subsystems::Tuple{Int64,Int64}
        )

        n_times = size(ensemble_trajectories)[1]  # number of time steps
        n_mcwf  = size(ensemble_trajectories)[2]  # number of Monte-Calro wave function trajectories

        totalcorrelations = Vector{Float64}(undef, n_times)
        classcorrelations = Vector{Float64}(undef, n_times)
        quantcorrelations = Vector{Float64}(undef, n_times)

        # ReentrantLock for multi-threading to prevent data race
        lock1 = ReentrantLock()
        lock2 = ReentrantLock()
        lock3 = ReentrantLock()

        @threads for t in 1:n_times  # mutli-threaded loop over time steps

            # computing bi-partite density matrix at current time step by averaging over the Monte-Carlo
            # wave functions
            dm = 1/n_mcwf * reduced(ensemble_trajectories[t,1], collect(subsystems))
            for i in 2:n_mcwf
                dm += 1/n_mcwf * reduced(ensemble_trajectories[t,i], collect(subsystems))
            end

            totalcorr_dm = totalcorrelation(dm)
            classcorr_dm = classicalcorrelation(dm)
            quantcorr_dm = quantumcorrelation(dm)
            @lock lock1 totalcorrelations[t] = totalcorr_dm
            @lock lock2 classcorrelations[t] = classcorr_dm
            @lock lock3 quantcorrelations[t] = quantcorr_dm
        end
        
        return totalcorrelations, classcorrelations, quantcorrelations
    end

end
