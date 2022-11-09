#### Graph-Theoretical Representation of Multiphoton Resonance Processes

[Hossein Jooya, *et. al.*, Scientific Reports 6, 37544](https://www.nature.com/articles/srep37544)

A graph-theoretical formalism is proposed to study generic circuit quantum electrodynamics systems consisting of a two level qubit coupled with a single-mode resonator in arbitrary coupling strength regimes beyond rotating-wave approximation. The intuitive and predictive picture provided by this method, and the simplicity of the mathematical construction, are demonstrated with some numerical studies of the multiphoton resonance processes and quantum interference phenomena for the superconducting qubit systems driven by intense ac fields.

### Usage

Simply compile and run the `Transverse_Coupling.f90` to obtain the transition probabilities and corresponding quasi-energies, as reported in Fig.(1) of the above referenced publication. This algorithm is suitable to model natural atoms which couple with electromagnetic fields at transverse mode due to the well-defined inversion symmetry of the potential energy. Within the Bloch 
representation, the time-dependent Hamiltonian of such two-level atom with transverse coupling is given by `H=-(1/2)[e0 sigma_z + e(t) sigma_x]` where `e0` is the separation between the energy levels, and `e(t)` is the microwave drive. `sigma_x` and `sigma_z` represent the Pauli matrices.

![alt text](https://github.com/hjooya/Chemical-Theory-and-Computation/blob/main/Graph%20Theoretical%20Approach%20to%20Multiphoton%20Absorption%20Spectra/Transverse_Coupling.jpg)

![alt text](https://github.com/hjooya/Chemical-Theory-and-Computation/blob/main/Graph%20Theoretical%20Approach%20to%20Multiphoton%20Absorption%20Spectra/Longitudinal_Coupling.jpg)

![alt text](https://github.com/hjooya/Chemical-Theory-and-Computation/blob/main/Graph%20Theoretical%20Approach%20to%20Multiphoton%20Absorption%20Spectra/Bidirectional_Coupling.jpg)




