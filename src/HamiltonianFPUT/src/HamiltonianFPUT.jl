
"""
The module HamiltonianFPUT contains several functions to build the Hamiltonian operator for a FPUT chain of
arbitrary length and coupling. It also provides functions to create jump operators and jumping rates for the
Monte-Carlo wave function technique to attach a thermal bath to single sites of the chain and evolve it
by the quantum optical master equation (Lindblad equation). Finally, some functions are also given to build
a locally excited initial state of the chain.
"""
module HamiltonianFPUT

    import QuantumOptics: Ket, Operator, FockBasis, fockstate, create, destroy, identityoperator, ⊗, dagger
    export HamOpFPUT, JumpOperators, JumpRates, localKets, operator_at_site

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
    Creating an multipartite operator on a chain with 'N' sites, which acts on one 'site' as the monopartite
    'operator' and as the identity operator on the remaining sites. 
    """    
    function operator_at_site(
        O::Operator,  # monopartite operator
        N::Int64,     # number sites in the chain
        site::Int64,  # site on which the operator acts
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
    Kinetic term of the Hamiltonian of a FPUT chain with 'N' elements
    """
    function KineticOpFPUT(
        N::Int64,               # number of chain elements
        singlebasis::FockBasis  # basis of single chain element
        )

        P = momentum(singlebasis)  # momentum oprator of single oscillator

        T = 0.5*operator_at_site(P, N, 1)^2
        for i in 2:N
            T += 0.5*operator_at_site(P, N, i)^2
        end

        return T
    end

    """
    Potential term of the Hamiltonian of a FPUT chain with 'N' elements and fixed ends.
    """    
    function PotentialOpFPUT_fixed(
        N::Int64,    # number of chain elements
        κ::Float64,  # strength of harmonic interaction
        α::Float64,  # strengt of cubic interaction
        β::Float64,  # strength of tetric interaction
        singlebasis::FockBasis  # basis of single chain element
        )

        Q = position(singlebasis)  # position operator of single oscillator

        # FPUT potential
        potential(O::Operator) = κ/2*O^2 + α/6*O^3 + β/24*O^4

        V = potential(-operator_at_site(Q,N,1))
        for i in 1:(N-1)
            V += potential(operator_at_site(Q,N,i) - operator_at_site(Q,N,i+1))
        end
        V += potential(operator_at_site(Q,N,N))
        return V
    end


    """
    Returns the Hamiltonian of a FPUT chain with 'N' elements and fixed ends for the coupling coefficients
    'κ', 'α' and 'β'.
    """
    function HamOpFPUT(
        N::Int64,    # number of chain elements
        κ::Float64,  # strength of harmonic interaction
        α::Float64,  # strengt of cubic interaction
        β::Float64,  # strength of tetric interaction
        singlebasis::FockBasis  # basis of single chain element
        )

        HamT = KineticOpFPUT(N, singlebasis)
        HamV = PotentialOpFPUT_fixed(N, κ,α,β, singlebasis)
        return HamT + HamV
    end


    """
    Planck distribution of a bosonic field at temperature 'kbT'
    """
    function PlanckDistribution(kbT::Float64)
        if kbT == 0.0
            return  0.0
        else
            return 1/(exp(1/kbT)-1)
        end
    end


    """
    Returns a vector containing the jump rates for a bosonic bath at temperature 'kbT' and coupling strength
    'γ' attached to the subsystems - indexed by 'BathSites' - of a multipartite system. If no vector 'BathSites'
    is passed only a single bath is assumed and thus, the information is not needed.
    """
    function JumpRates(kbT::Float64, γ::Float64)
        η = PlanckDistribution(kbT)
        return [γ*(η+1),γ*η]
    end

    function JumpRates(kbT::Float64, γ::Float64, BathSites::Vector{Int64})
        η = PlanckDistribution(kbT)
        rates = Vector{Float64}(undef, 2*length(BathSites))
        for i in eachindex(BathSites)
            rates[i] = γ*(η+1)
            rates[length(BathSites)+i] = γ*η
        end
        return rates
    end


    """
    Returns a vector containing the jump operators to the subsystems - indexed by 'BathSites' - of a multipartite
    system. If only a single integer 'BathSite' is passed only a single bath at this side is assumed.
    """
    function JumpOperators(J::Operator, N::Int64, BathSite::Int64)
        return [operator_at_site(J,N,BathSite), dagger(operator_at_site(J,N,BathSite))]
    end

    function JumpOperators(J::Operator, N::Int64, BathSites::Vector{Int64})
        operators = Vector{Operator}(undef, 2*length(BathSites))
        for i in 1:length(BathSites)
            operators[i] = operator_at_site(J,N,BathSites[i])
            operators[length(BathSites)+i] = dagger(operator_at_site(J,N,BathSites[i]))
        end
        return operators
    end



    """
    This function builds the ket state of a FPUT chain with 'N' elements where the element at 'site' is excited
    to the state 'ket' while all other chain elements remain in their ground state (i.e. Fock state |0>).  
    """
    function localKets(
        ketS::Ket,  # state vector of initial exitation at site k
        N::Int64;  # number of chain elements
        siteS::Int64=1  # site of system oscillator
        )

        singlebasis = ketS.basis
        ketCHAIN = ketS

        # iteratively attaching all necessary ground states to the left
        if siteS>1
            for i in 2:site
                ketCHAIN = fockstate(singlebasis,0)⊗ketCHAIN
            end
        end

        # iteratively attaching all remaining ground states to the right
        if siteS<N
            for i in (siteS+1):N
                ketCHAIN = ketCHAIN⊗fockstate(singlebasis,0)
            end
        end

        return ketCHAIN
    end

end