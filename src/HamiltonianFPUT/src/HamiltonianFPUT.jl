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
    Building a operator at one 'site' of a FPUT chain with 'N' elements. 'singlebasis' is the Fock basis of each single
    chain element.
    """    
    function operator_at_site(
        O::Operator,   # single oscillator operator
        N::Int64,      # number of chain elements
        site::Int64, # site of operator
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
    Builds the kinetic term of the Hamiltonian of a FPUT chain with 'N' elements
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
    Builds the potential term of the Hamiltonian of a FPUT chain with 'N' elements and fixed ends.
    """    
    function PotentialOpFPUT_fixed(
        N::Int64,    # number of chain elements
        κ::Float64,  # strength of harmonic interaction
        α::Float64,  # strengt of cubic interaction
        β::Float64,  # strength of tetric interaction
        singlebasis::FockBasis;  # basis of single chain element

        potential_factor::Float64 = 1.0
        )

        Q = position(singlebasis)  # position operator of single oscillator

        # potentail from ONE side to each chain element
        potential(O::Operator) = potential_factor * (κ/2*O^2 + α/6*O^3 + β/24*O^4)

        V = potential_factor * potential(-operator_at_site(Q,N,1))
        for i in 1:(N-1)
            V += potential_factor * potential(operator_at_site(Q,N,i) - operator_at_site(Q,N,i+1))
        end
        V += potential_factor * potential(operator_at_site(Q,N,N))
        return V
    end


    """
    Builds the potential term of the Hamiltonian of a FPUT chain with 'N' elements and closed ends, i.e. periodic boundry conditions.
    Only for developement usage so far
    """
    function PotentialOpFPUT_periodic(
        N::Int64,    # number of chain elements
        κ::Float64,  # strength of harmonic interaction
        α::Float64,  # strengt of cubic interaction
        β::Float64,  # strength of tetric interaction
        singlebasis::FockBasis;  # basis of single chain element

        potential_factor::Float64 = 1.0
        )

        Q = position(singlebasis)  # position operator of single oscillator
        id = identityoperator(singlebasis)

        # potentail from ONE side to each chain element
        potential(O::Operator) = potential_factor * (κ/2*O^2 + α/6*O^3 + β/24*O^4)

        # build potential/interaction Hamiltonian; best try out by hand to see that it does ;-)
        Q_left  = Q⊗id  # potition operator of left most chain element
        Q_right = id⊗Q  # potition operator of right most chain element

        I_left = potential(Q⊗id-id⊗Q) # interaction to the left
        I_right = I_left                # interaction to the right
        I = I_left                      # interaction for 2 chain elements

        if N<=2
            return I
        else
            ## bulidng iteratively the interactions within the chain
            for i in 3:N
                Q_left = Q_left⊗id      # trivial extension of coupling to left fixed end
                Q_right = id⊗Q_right    # trivial extension of coupling to right fixed end
                I_left = I⊗id         # new left part of interaction is old interaction with trivial extension to the right
                I_right = id⊗I_right  # new right part of interaction is old right part with trivial extension to the left 
                I = I_left + I_right   # new interaction is new left part plus new right part
            end
            I_closing = potential(Q_right - Q_left)  # interaction between left most and right most ends to close the loop
            return  I + I_closing
        end
    end

    """
    Returns the Hamiltonian of a FPUT chain with 'N' elements and fixed ends.
    """
    function HamOpFPUT(
        N::Int64,    # number of chain elements
        κ::Float64,  # strength of harmonic interaction
        α::Float64,  # strengt of cubic interaction
        β::Float64,  # strength of tetric interaction
        singlebasis::FockBasis;  # basis of single chain element
        potential_factor::Float64 = 1.0,
        )

        HamT = KineticOpFPUT(N, singlebasis)
        HamV = PotentialOpFPUT_fixed(N, κ,α,β, singlebasis; potential_factor=potential_factor)
        return HamT + HamV
    end


    """
    Planck distribution of a bosonic field at temperature kbT
    """
    function PlanckDistribution(kbT::Float64)
        if kbT == 0.0
            return  0.0
        else
            return 1/(exp(1/kbT)-1)
        end
    end


    """
    Returns a vector containing the jump rates for a bosonic bath at temperature 'kbT' and coupling strength 'γ' attached
    to the subsystems - indexed by 'BathSites' - of a multipartite system. If no vector 'BathSites' is passed only a single
    bath is assumed and thus, the information is not needed.
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
    Returns a vector containing the jump operators to the subsystems - indexed by 'BathSites' - of a multipartite system. If
    only a single integer 'BathSite' is passed only a single bath at this side is assumed.
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
    This function builds the ket state of a FPUT chain with 'N' elements where the the element at 'site' is excited to the state 'ket'
    while all other chain elements remain in their ground state (i.e. Fock state |o>).  
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


    # """
    # This function builds the ket state of a FPUT chain with 'N' elements where the the element at 'site' is excited to the state 'ket'
    # while all other chain elements remain in their ground state (i.e. Fock state |o>).  
    # """
    # function EntangledKets(
    #     ket1::Ket,  # state vector of initial exitation at site k
    #     ket2::Ket,  # state vector of initial exitation at site k
    #     N::Int64,   # number of chain elements
    #     plus::Bool=true; # signe of superposition
    #     siteA::Int64=1,  # site of 1st oscillator in the entanglement
    #     siteB::Int64=N   # site of 2nd oscillator in the entanglement
    #     )

    #     if ket1.basis != ket2.basis
    #         error("ket1 and ket2 do not have the same basis")
    #     else
    #         singlebasis = ket1.basis
    #     end

    #     ketCHAIN1 = ket1
    #     ketCHAIN2 = ket2

    #     # iteratively attaching all necessary ground states to the left
    #     if siteA>1
    #         for i in 2:site
    #             ketCHAIN1 = fockstate(singlebasis,0)⊗ketCHAIN1
    #             ketCHAIN2 = fockstate(singlebasis,0)⊗ketCHAIN2
    #         end
    #     end

    #     # iteratively attaching all remaining ground states to the right
    #     if siteS<N
    #         for i in (siteS+1):N
    #             ketCHAIN1 = ketCHAIN1⊗fockstate(singlebasis,0)
    #         end
    #     end

    #     return ketCHAIN
    # end

end