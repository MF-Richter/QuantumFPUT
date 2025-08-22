module Correlations

    import QuantumOptics: Operator, dense, Ket, reduced, dm, ⊗, tracedistance, identityoperator, dagger, eigenstates
    using  Base.Threads
    export totalcorrelation, classicalcorrelation, quantumcorrelation, correlations


    """
    Correlations measured by the distance between full state and its uncorrelated counterpart
    """    
    function totalcorrelation(ρ::Operator)
        return tracedistance(ρ, reduced(ρ,1)⊗reduced(ρ,2))
    end

    """
    This function takes a state of the full FPUT chain and computes its corresponding classically
    correlated state.
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
    Classical Correlation measured by the distance between the classified version of a state and its
    fully decorrelated version
    """    
    function classicalcorrelation(ρ::Operator)
        return tracedistance( classify_correlations(ρ), reduced(ρ, 1)⊗reduced(ρ, 2) )
    end

    """
    Quantum Correlation measured by the distance between the classified version of a state and its
    fully correlated version
    """    
    function quantumcorrelation(ρ::Operator)
        return tracedistance(ρ, classify_correlations(ρ))
    end


    """
    A function that takes a trajectory of ket states through time or an ensemble of such as well as a tupel containg
    the indices of the two subsystems between which one wants to measure correlations.
    It returns three vectors containg the total correlation, the classical correlation and the quantum correlation
    through time between these two subsystems.

    If just a ket trajectory is passed the additional kwarg 'only_total_correlation' (default is 'true') can be used
    to only compute the total correlation. Since ket states cannot posses correlation different than entanglement the
    distingution between classical and quantum correlation becomes irrelevant - every correlation is of genuine quan-
    tum origin (however, entanglement can trigger classically looking correlations too, although their origin must be
    quantum)
    """
    function correlations(
        trajectory::Vector{<:Ket},
        subsystems::Tuple{Int64,Int64};
        only_total_correlation::Bool = false
        )

        n_times = length(trajectory)

        totalcorrelations = Vector{Float64}(undef, n_times)

        if only_total_correlation == false
            classcorrelations = Vector{Float64}(undef, n_times)
            quantcorrelations = Vector{Float64}(undef, n_times)
        end

        lock1 = ReentrantLock()
        lock2 = ReentrantLock()
        lock3 = ReentrantLock()
        @threads for τ in 1:n_times

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

        n_times = size(ensemble_trajectories)[1]
        n_mcwf  = size(ensemble_trajectories)[2]

        totalcorrelations = Vector{Float64}(undef, n_times)
        classcorrelations = Vector{Float64}(undef, n_times)
        quantcorrelations = Vector{Float64}(undef, n_times)

        lock1 = ReentrantLock()
        lock2 = ReentrantLock()
        lock3 = ReentrantLock()
        @threads for t in 1:n_times
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
