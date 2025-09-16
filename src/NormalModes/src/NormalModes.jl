
"""
The module 'NormalModes' contains several functions concerning the normal modes of a harmonic chain, such as
normal mode quadrature operators, frequencies, and Hamiltonian operators. Additionally, it also provides a
function to initialize the chain in a normal mode coherent state, i.e., a collective displacement in phase
space according to the given mode.
"""
module NormalModes

    import QuantumOptics: Ket, Operator, FockBasis, fockstate, create, destroy, identityoperator, ⊗, dagger, exp, expect, eigenstates, dense, norm
    export QuadraturesNormalMode, FrequenceNormalMode, HamiltonianNormalMode, DisplacementoperatorNormalMode, coherentstateNM, means

    """
    Position operator of single oscillator
    """
    function position(basis::FockBasis)
        a = destroy(basis)
        at = create(basis)
        return 0.5( 1/sqrt(2)*(at+a) + dagger(1/sqrt(2)*(at+a)) )
    end

    """
    Momentum operator of single oscillator
    """    
    function momentum(basis::FockBasis)
        a = destroy(basis)
        at = create(basis)
        return 0.5( -1*im/sqrt(2)*(a-at) + dagger(-1*im/sqrt(2)*(a-at)) )
    end

    """
    Building a operator at one 'site' of a FPUT chain with 'N' elements. 'singlebasis' is
    the Fock basis of each single site of the chain.
    """   
    function operator_at_site(
        O::Operator,   # single oscillator operator
        N::Int64,      # number of sites in chain
        site::Int64,   # site of operator
        )

        basis = O.basis_l
        id = identityoperator(basis)
        operatorFPUT = O
        if site>1
            for i in 2:site
                operatorFPUT = id⊗operatorFPUT
            end
        end

        if site<N
            for i in (site+1):N
                operatorFPUT = operatorFPUT⊗id
            end
        end

        return operatorFPUT
    end


    """
    Position and momentum operator of the normal mode 'mode' of a harmonic chain with 'N' sites. 
    'singlebasis' is the Fock basis of each single site of the chain.
    """
    function QuadraturesNormalMode(
        N::Int64,  # number of sites in chain
        mode::Int64,  # level of normal mode
        singlebasis::FockBasis
        )

        positionNM = sqrt(2/(N+1)) * operator_at_site(position(singlebasis),N,1) * sin(pi*mode*1/(N+1))
        momentumNM = sqrt(2/(N+1)) * operator_at_site(momentum(singlebasis),N,1) * sin(pi*mode*1/(N+1))

        for i in 2:N
            positionNM += sqrt(2/(N+1)) * operator_at_site(position(singlebasis),N,i) * sin(pi*mode*i/(N+1))
            momentumNM += sqrt(2/(N+1)) * operator_at_site(momentum(singlebasis),N,i) * sin(pi*mode*i/(N+1))
        end

        return positionNM, momentumNM
    end

    """
    Frequency of the normal mode 'mode' of a harmonic chain with 'N' sites.
    """    
    function FrequenceNormalMode(
        N::Int64,         # number of sites in chain
        mode::Int64;      # index of normal mode
        κ::Float64 = 1.0  # harmonic coupling strength
        )
        return 2*sqrt(κ)*sin(0.5*pi*mode/(N+1))
    end

    """
    Hamiltonian of the normal mode 'mode' of a harmonic chain with 'N' sites. If a vector 'modes'
    is passed, a vector of normal mode Hamiltonians is returned. 'singlebasis' is the Fock basis
    of each single site of the chain.
    """
    function HamiltonianNormalMode(
        N::Int64,     # number of sites in chain
        mode::Int64,  # index of normal mode
        singlebasis::FockBasis;  
        κ::Float64 = 1.0 # harmonic coupling strength
        )

        Q, P = QuadraturesNormalMode(N,mode,singlebasis)
        ω = FrequenceNormalMode(N,mode; κ=κ)
        HamNM = 0.5 * (P*dagger(P) + ω^2*Q*dagger(Q))

        return (HamNM + dagger(HamNM))/2
    end

    function HamiltonianNormalMode(
        N::Int64,              # number of sites in chain
        modes::Vector{Int64},  # indeces of normal modes
        singlebasis::FockBasis;
        κ::Float64=1.0  # harmonic coupling strength
        )

        HamOps = Vector{Operator}()
        for k in modes
            Ham = HamiltonianNormalMode(N, k, singlebasis; κ=κ)
            push!(HamOps, Ham)
        end

        return HamOps
    end

    """
    Normal mode displacement operator, i.e. the usual displacement operator yet with normal mode
    quaratures of nomral mode 'mode' of a chain of 'N' sites. 'singlebasis' is the Fock basis of
    each single site of the chain. 'x' is the displacement vector in the phase-space of the mode,
    while 'ϕ' is the corresponding complex number for the complex plain representation of the
    phase-space
    """
    function DisplacementoperatorNormalMode(
        x::Vector{Float64},  # normal mode phase-space displacement
        N::Int64,            # number of sites in chain
        mode::Int64,         # index of normal mode
        singlebasis::FockBasis
        )
        Q,P = QuadraturesNormalMode(N,mode,singlebasis)
        Ω = [0.0 1.0; -1.0 0.0]
        return QuantumOptics.exp(-1im*transpose(Ω*x)*[Q,P])
    end

    function DisplacementoperatorNormalMode(
        ϕ::ComplexF64,  # normal mode phase-space displacement
        N::Int64,       # number of sites in chain
        mode::Int64,    # index of normal mode
        singlebasis::FockBasis
        )

        x = sqrt(2)*[real(ϕ), imag(ϕ)]
        Q,P = QuadraturesNormalMode(N,mode,singlebasis)
        Ω = [0.0 1.0; -1.0 0.0]
        return exp(-1im*transpose(Ω*x)*[Q,P])
    end

    """
    Excitation of a harmonic chain of 'N' sites in a normal mode coherent state |ϕ>.
    """
    function coherentstateNM(
        singlebasis::FockBasis,
        ϕ,         # normal mode phase-space displacement
        N::Int64;  # number of sites in chain
        mode::Int64 = 1 # index of normal mode
        )

        vacuumstateNM = fockstate(singlebasis,0)
        for i in 2:N
            vacuumstateNM = vacuumstateNM⊗fockstate(singlebasis,0)
        end
        Dop = DisplacementoperatorNormalMode(ϕ,N,mode,singlebasis)

        return Dop*vacuumstateNM
    end


    """
    computing mean values for an observable (or a collection of such) over a trajectory of
    ket states or an ensemble of Monte-Carlo wave function trajectories
    """
    function means(
        operator::Operator,
        trajectory::Vector{<:Ket}
        )
        mean(ket::Ket) = real(expect(operator,ket))
        return mean.(trajectory)
    end

    function means(
        operators::Vector{<:Operator},
        trajectory::Vector{<:Ket}
        )
        meanvalues = Array{Float64, 2}(undef, length(trajectory), length(operators))
        for i in eachindex(operators)
            op = operators[i]
            meanvalues[:,i] = means(op, trajectory)
        end
        return meanvalues
    end

    function means(
        operator::Operator,
        ensemble_trajectories::Array{Ket, 2}
        )
        mean(ket::Ket) = real(expect(operator,ket))
        n_mcwf = size(ensemble_trajectories)[2]

        meanvalues = 1/n_mcwf * mean.(ensemble_trajectories[:,1])
        for i in 2:n_mcwf
            meanvalues += 1/n_mcwf * mean.(ensemble_trajectories[:,i])
        end

        return meanvalues
    end

    function means(
        operators::Vector{<:Operator},
        ensemble_trajectories::Array{Ket, 2}
        )
        meanvalues = Array{Float64, 2}(undef, size(ensemble_trajectories)[1], length(operators))
        for i in eachindex(operators)
            op = operators[i]
            meanvalues[:,i] = means(op, ensemble_trajectories)
        end
        return meanvalues
    end    

end
