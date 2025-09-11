# QuantumFPUT

[QuantumFPUT](https://github.com/MF-Richter/QuantumFPUT) is a julia project that together with its twin project [ClassicalFPUT](https://github.com/MF-Richter/ClassicalFPUT) is written to investigate quantum and classical information flow through FPUT chains.

The FPUT chain was firstly simulated by Enrico Fermi, John Pasta, Stanislaw Ulam and Mary Tsingou and consists of oscillating particles coupled to their nearest neighbors by the potential

$V(\Delta q) = \frac{\kappa}{2}\Delta q^2 + \frac{\alpha}{3!}\Delta q^3 + \frac{\beta}{4!}\Delta q^4$

i.e. a harmonic chain of coupling strength $\kappa$ with some anharmonic perturbation of strength $\alpha$ and $\beta$. FPUT studied this chain originally in the context of energy sharing between its normal modes to test the ergodicity hypothesis (which in fact they found not to be fulfilled for this model[1]).

[QuantumFPUT](https://github.com/MF-Richter/QuantumFPUT) allows to construct such a FPUT chain of quantum oscillators of arbitrary length and coupling strengths. The whole project is build using the package [QuantumOptics](https://github.com/qojulia/QuantumOptics.jl), the differential equations are solved under the hood using [DifferentialEquations](https://github.com/SciML/DifferentialEquations.jl). The Hilbert space effectively is the tensor product of the Hilbert spaces for the single oscillators (called sites). [QuantumFPUT](https://github.com/MF-Richter/QuantumFPUT) provides tools to construct their mutual Hamiltonian and initial states by exciting either a single site by a coherent state or a Fock state, or excite the whole chain in one of its normal modes, i.e. each site starts in a coherent state of the corresponding displacement. The chain is then evolved by its Hamiltonian according to the Schroedinger equation. Additionally, one can attach a thermal bath at single sites and evolve the chain by a Lindblad master equation using the quantum jump method with Monte-Carlo wave function technique [2,3,4]. [QuantumFPUT](https://github.com/MF-Richter/QuantumFPUT) further provides modules to reduce the chains ket states or ensembles of such to local density operators at single sites and compute further from them data like information flow, correlations, phase-space means etc.

To study the information flow through a FPUT chain and the effects of anharmonic coupling to it [QuantumFPUT](https://github.com/MF-Richter/QuantumFPUT) uses the Breuer-Laine-Piilo (BLP) measure based on the trace distance between quantum density operators [5,6,7], while [ClassicalFPUT](https://github.com/MF-Richter/ClassicalFPUT) uses the classical limit of the BLP measure based on the Kolmogorov distance between phase-space distributions [8,9]. The phase-space distributions can be approximated by initializing the chain in an ensemble as described above.



## Author
- Moritz F. Richter


## Installation & Setup

This code base is using the [Julia Language](https://julialang.org/) and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/) to make a reproducible scientific project named
To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "QuantumFPUT"
```
which auto-activate the project and enable local path handling from DrWatson.


## Acknowledgments
This work has been supported by the German Research Foundation (DFG) through FOR 5099.


## Data Availability
Raw data not included. Please request from author.


## References
[1] S. Ulam (and M. Tsingou) E. fermi, J. Pasta. Studies of nonlinear problems. Los Alamos Report, (LA-1940), 1955.

[2] Jean Dalibard, Yvan Castin, and Klaus Mølmer. Wave-function approach to dissipative processes in quantum optics. Phys. Rev. Lett., 68:580–583, Feb 1992.

[3] Klaus Mølmer, Yvan Castin, and Jean Dalibard. Monte carlo wave-function method in quantum optics. J. Opt. Soc. Am. B, 10(3):524–538, Mar 1993.

[4] Howard Carmichael. An open systems approach to quantum opticslectures presented at the Université Libre de Bruxelles, October 28 to November 4, 1991. Springer, Berlin, 1993.

[5] Heinz-Peter Breuer, Elsi-Mari Laine, and Jyrki Piilo. Measure for the degree of non-Markovian behavior of quantum processes in open systems. Phys. Rev. Lett., 103:210401,2009.

[6] Elsi-Mari Laine, Jyrki Piilo, and Heinz-Peter Breuer. Measure for the non-markovianity of quantum processes. Phys. Rev. A, 81:062115, Jun 2010.

[7] E.-M. Laine, J. Piilo, and H.-P. Breuer. Witness for initial system-environment correlations in open-system dynamics. Europhysics Letters, 92(6):60010, jan 2011.

[8] Moritz Ferdinand Richter, Raphael Wiedenmann, and Heinz-Peter Breuer. Witnessing non-markovianity by quantum quasi-probability distributions. New Journal of Physics, 2022.

[9] Moritz F. Richter and Heinz-Peter Breuer. Phase-space measures of information flow in open systems: A quantum and classical perspective of non-markovianity. Phys. Rev. A, 110:062401, Dec 2024.


## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Release Note
This is still a beta version. The Documentation and single commend in the code are still work in progress