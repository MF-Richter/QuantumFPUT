"""
The module 'StateCharacteristics' contains several quantities, such as mean values and (co)variances of
number operators or phase-space quadratures, as well as benchmarking quantities for the dimension cutoff.
"""
module StateCharacteristics

    import LinearAlgebra: det, eigvals
    import QuantumOptics: Operator, FockBasis, destroy, create, number, dagger, tr, tracenorm
    export meannumber,variancenumber,
    meanposition,meanmomentum,meandisplacement,phasespaceangle,meanPhSp,
    covariance_qq,covariance_pp,covariance_qp,covariancematrix,covariancedeterminant,covariancesqueezing,
    relative_dimcutoffvalue,tracenorm_diviation

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
    Mean value of number operator for given states
    """
    function meannumber(state::Operator)
        basis = state.basis_l
        return tr( number(basis)*state )
    end

    """
    Variance of number operator for given state
    """    
    function variancenumber(state::Operator)
        basis = state.basis_l
        return tr( number(basis)^2 * state ) - meannumber(state)^2
    end

    """
    Mean position value in optical phase space of state in its baiss
    """
    function meanposition(state::Operator)
        basis = state.basis_l
        Q = position(basis)
        return real(tr(Q*state))
    end

    """
    Mean position value in optical phase space of state in its baiss
    """    
    function meanmomentum(state::Operator)
        basis = state.basis_l
        P = momentum(basis)
        return real(tr(P*state))
    end

    """
    Vector of phase space means
    """
    function meanPhSp(state::Operator)
        return [meanposition(state), meanmomentum(state)]
    end

    """
    Mean displacement in optical phase space for given states
    """    
    function meandisplacement(state::Operator)
        q = meanposition(state)
        p = meanmomentum(state)
        return sqrt(q^2+p^2)
    end

    """
    Angle of mean displacement vector in optical phase space
    """
    function phasespaceangle(state::Operator)
        q = meanposition(state)
        p = meanmomentum(state)
        return atan(p,q)
    end

    """
    position-position entry of the covariance matrix of a 'state' given by its density matrix
    """
    function covariance_qq(state::Operator)
        basis = state.basis_l
        Q = position(basis)
        q = meanposition(state)
        return 2*real(tr(Q^2*state)) - 2*q^2
    end

    """
    momentum-momentum entry of the covariance matrix of a 'state' given by its density matrix
    """
    function covariance_pp(state::Operator)
        basis = state.basis_l
        P = momentum(basis)
        p = meanmomentum(state)
        return 2*real(tr(P^2*state)) - 2*p^2
    end

    """
    position-momentum entry of the covariance matrix of a 'state' given by its density matrix
    """
    function covariance_qp(state::Operator)
        basis = state.basis_l
        Q, P = position(basis), momentum(basis)
        q, p = meanposition(state), meanmomentum(state)
        return real(tr(Q*P*state + P*Q*state)) - 2*q*p
    end

    """
    Covariance matrix of a 'state' given by its density matrix
    """
    function covariancematrix(state::Operator)
        Σqq = covariance_qq(state)
        Σpp = covariance_pp(state)
        Σqp = covariance_qp(state)
        return [Σqq Σqp; Σqp Σpp]
    end

    """
    Determinante of the covariance matrix of a 'state' given by its density matrix
    """
    function covariancedeterminant(state::Operator)
        return det(covariancematrix(state))  
    end

    """
    Squeezing parameter of the covariance matrix of a 'state' given by its density matrix
    """
    function covariancesqueezing(state::Operator)
        σ1,σ2 = abs.(eigvals(covariancematrix(state)))
        return 1-σ1/σ2
    end


    """
    Ratio between absolute value of the entry of ρ at dim. cutoff and the maximal value of ρ.
    Can be used as measure of error done by dim. cutoff.
    """
    function relative_dimcutoffvalue(ρ::Operator)
        dm = ρ.data
        dim = size(dm)[1]
        dm_max = abs(maximum(abs.(dm)))
        return abs(dm[dim,dim])/dm_max
    end

    "This function computes the absolute deviation of the trace norm of a density operator 'ρ' from one"
    function tracenorm_diviation(ρ::Operator)
        return abs(1-tracenorm(ρ))
    end

end 
