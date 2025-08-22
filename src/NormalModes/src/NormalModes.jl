module NormalModes

    import QuantumOptics: Ket, Operator, FockBasis, fockstate, create, destroy, identityoperator, ⊗, dagger, exp, expect, eigenstates, dense, norm
    export QuadraturesNormalMode, FrequenceNormalMode, NumberoperatorNormalMode, HamiltonianNormalMode, DisplacementoperatorNormalMode, coherentstateNM, fockstateNM, means

    function position(basis::FockBasis)
        """
        Position operator of single oscillator
        """
        a = destroy(basis)
        at = create(basis)
        return 0.5( 1/sqrt(2)*(at+a) + dagger(1/sqrt(2)*(at+a)) )
    end

    function momentum(basis::FockBasis)
        """
        Momentum operator of single oscillator
        """
        a = destroy(basis)
        at = create(basis)
        return 0.5( -1*im/sqrt(2)*(a-at) + dagger(-1*im/sqrt(2)*(a-at)) )
    end

    function operator_at_site(
        O::Operator,   # single oscillator operator
        N::Int64,      # number of chain elements
        site::Int64,   # site of operator
        )
        """
        Building a operator at one 'site' of a FPUT chain with 'N' elements. 'singlebasis' is the Fock basis of each single
        chain element.
        """
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


    function QuadraturesNormalMode(
        N::Int64,  # number of chain elements
        mode::Int64,  # level of normal mode
        singlebasis::FockBasis
        )
        """
        Position and momentum operator of the 'k'-th normal mode of a harmonic FPUT chain with 'N' elements
        """
        positionNM = sqrt(2/(N+1)) * operator_at_site(position(singlebasis),N,1) * sin(pi*mode*1/(N+1))
        momentumNM = sqrt(2/(N+1)) * operator_at_site(momentum(singlebasis),N,1) * sin(pi*mode*1/(N+1))

        for i in 2:N
            positionNM += sqrt(2/(N+1)) * operator_at_site(position(singlebasis),N,i) * sin(pi*mode*i/(N+1))
            momentumNM += sqrt(2/(N+1)) * operator_at_site(momentum(singlebasis),N,i) * sin(pi*mode*i/(N+1))
        end

        return positionNM, momentumNM
    end

    function FrequenceNormalMode(
        N::Int64,
        mode::Int64;
        potential_factor::Float64 = 1.0
        )
        """
        Frequency of the 'k'-th normal mode of a harmonic FPUT chain with 'N' elements
        """
        return 2*sqrt(potential_factor)*sin(0.5*pi*mode/(N+1))
    end

    function LaddersNormalMode(
        N::Int64,  # number of chain elements
        mode::Int64,  # level of normal mode
        singlebasis::FockBasis; 
        potential_factor::Float64 = 1.0
        ) 
        """
        Creation and destruction operator of the 'k'-th normal mode of a harmonic FPUT chain with 'N' elements
        """
        Q, P = QuadraturesNormalMode(N, mode, singlebasis)
        Ω = FrequenceNormalMode(N, mode; potential_factor=potential_factor)
        destroyNM = sqrt(Ω/2)*Q + 1im/sqrt(2*Ω)*P
        return destroyNM, dagger(destroyNM)
    end

    function HamiltonianNormalMode(
        N::Int64,     # number of chain elements
        mode::Int64,  # level of normal mode
        singlebasis::FockBasis;
        potential_factor::Float64 = 1.0,
        build_by_quadratures::Bool = true
        )

        if build_by_quadratures == true
            Q, P = QuadraturesNormalMode(N,mode,singlebasis)
            ω = FrequenceNormalMode(N,mode; potential_factor=potential_factor)
            HamNM = 0.5 * (P*dagger(P) + ω^2*Q*dagger(Q))
        else
            A, At = LaddersNormalMode(N, mode, singlebasis; potential_factor=potential_factor)
            idFPUT = operator_at_site(identityoperator(singlebasis), N, 1)
            ω = FrequenceNormalMode(N,mode; potential_factor=potential_factor)
            HamNM = ω*(At*A + 0.5*idFPUT)
        end

        return (HamNM + dagger(HamNM))/2
    end

    function HamiltonianNormalMode(
        N::Int64,  # number of chain elements
        modes::Vector{Int64},  # level of normal mode
        singlebasis::FockBasis;
        potential_factor::Float64=1.0
        )

        HamOps = Vector{Operator}()
        for k in modes
            Ham = HamiltonianNormalMode(N, k, singlebasis; potential_factor=potential_factor)
            push!(HamOps, Ham)
        end

        return HamOps
    end



    function NumberoperatorNormalMode(
        N::Int64,  # ...number of chain elements
        mode::Int64,
        singlebasis::FockBasis;
        potential_factor::Float64 = 1.0
        )

        # idCHAIN = operator_at_site(identityoperator(singlebasis), N, 1)
        # Q, P = QuadraturesNormalMode(N,mode,singlebasis)
        # ω = FrequenceNormalMode(N,mode; potential_factor=potential_factor)
        # return 0.5 * ( P^2/ω + ω*Q^2 - idCHAIN )

        A, At = LaddersNormalMode(N, mode, singlebasis; potential_factor=potential_factor)
        return At*A
    end


    function NumberoperatorNormalMode(
        N::Int64,  # number of chain elements
        modes::Vector{Int64},  # level of normal mode
        singlebasis::FockBasis;
        potential_factor::Float64=1.0
        )

        Noperators = Vector{Operator}()
        for k in modes
            Nop = NumberoperatorNormalMode(N, k, singlebasis; potential_factor=potential_factor)
            push!(Noperators, Nop)
        end

        return Noperators
    end


    function DisplacementoperatorNormalMode(
        x::Vector{Float64},
        N::Int64,  # number of chain elements
        mode::Int64,  # level of normal mode
        singlebasis::FockBasis
        )
        Q,P = QuadraturesNormalMode(N,mode,singlebasis)
        Ω = [0.0 1.0; -1.0 0.0]
        return QuantumOptics.exp(-1im*transpose(Ω*x)*[Q,P])
    end

    function DisplacementoperatorNormalMode(
        ϕ::ComplexF64,
        N::Int64,  # number of chain elements
        mode::Int64,  # level of normal mode
        singlebasis::FockBasis
        )

        x = sqrt(2)*[real(ϕ), imag(ϕ)]
        Q,P = QuadraturesNormalMode(N,mode,singlebasis)
        Ω = [0.0 1.0; -1.0 0.0]
        return exp(-1im*transpose(Ω*x)*[Q,P])
    end

    function coherentstateNM(
        singlebasis::FockBasis,
        ϕ,
        N::Int64;
        mode::Int64 = 1
        )

        vacuumstateNM = fockstate(singlebasis,0)
        for i in 2:N
            vacuumstateNM = vacuumstateNM⊗fockstate(singlebasis,0)
        end
        Dop = DisplacementoperatorNormalMode(ϕ,N,mode,singlebasis)

        return Dop*vacuumstateNM
    end

    function fockstateNM(
        singlebasis::FockBasis,
        n,
        N::Int64;
        mode::Int64 = 1
        )

        CreateNM = LaddersNormalMode(N,mode,singlebasis)[1]

        vacuumstateNM = fockstate(singlebasis,0)
        for i in 2:N
            vacuumstateNM = vacuumstateNM⊗fockstate(singlebasis,0)
        end
        return CreateNM^n * vacuumstateNM/norm(CreateNM^n * vacuumstateNM)

    
        # HamNM = HamiltonianNormalMode(N,mode,singlebasis)
        # eigenvalsNM, eigenstatesNM = eigenstates(dense(HamNM))
        # groundstateNM = eigenstatesNM[1]
        # return CreateNM^n * groundstateNM/norm(CreateNM^n * groundstateNM)
    end



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
