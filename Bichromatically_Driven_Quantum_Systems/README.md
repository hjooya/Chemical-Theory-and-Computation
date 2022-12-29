## Absorption Spectroscopy of Bichromatically Driven Quantum Systems 

[Hossein Jooya, *et. al.*, Phys. Rev. B 96, 174518](https://www.nature.com/articles/srep37544)

![](https://github.com/hjooya/Chemical-Theory-and-Computation/blob/main/Bichromatically_Driven_Quantum_Systems/Bichromatic_Driven_TLS.gif)

<img src="https://github.com/hjooya/Chemical-Theory-and-Computation/blob/main/Bichromatically_Driven_Quantum_Systems/Bichromatic_Driven_TLS.gif" width="250" height="250"/>



A graph-theoretical formalism is proposed to study generic circuit quantum electrodynamics systems consisting of a two level qubit coupled with a single-mode resonator in arbitrary coupling strength regimes beyond rotating-wave approximation. The intuitive and predictive picture provided by this method, and the simplicity of the mathematical construction, are demonstrated with some numerical studies of the multiphoton resonance processes and quantum interference phenomena for the superconducting qubit systems driven by intense ac fields.

### Usage

Simply compile the Makefile, 'Makefile_Transverse_Coupling', to run the `Code_Transverse_Coupling.f90` and obtain the transition probabilities and corresponding quasi-energies, as reported in Fig.(1) of the above referenced publication. Notice how only `odd` transitions are aloowe. This algorithm is suitable to model natural atoms which couple with electromagnetic fields at transverse mode due to the well-defined inversion symmetry of the potential energy. Within the Bloch 
representation, the time-dependent Hamiltonian of such two-level atom with transverse coupling is given by `H=-(1/2)[e0 sigma_z + e(t) sigma_x]` where `e0` is the separation between the energy levels, and `e(t)` is the microwave drive. `sigma_x` and `sigma_z` represent the Pauli matrices.


Since the potential energy for superconducting qubits can be tuned, the inversion symmetry for these artificial atoms can be broken and `even` multiphoton 
processes can be observed. The existence of the longitudinal coupling between superconducting qubits and applied magnetic fields have been observed, when the inversion symmetry of the potential energy of the superconducting qubit is broken. The time dependent Hamiltonian with the longitudinal coupling is given by `H=-(1/2)[e(t) sigma_z + D sigma_x]` where `D (Delta)` is the tunnel splitting. Simply compile and run `Code_Longitudinal_Coupling.f90` and see the results are shown below and presented in Fig.(2) of the cited publication. Notice how both `odd` and `even` transitions are observed in this model.

![alt text](https://github.com/hjooya/Chemical-Theory-and-Computation/blob/main/Graph%20Theoretical%20Approach%20to%20Multiphoton%20Absorption%20Spectra/Image_Longitudinal_Coupling.jpg)

It has been shown that superconducting qubits and external fields can have both transverse and longitudinal coupling. The `Code_Bidirectional_Coupling.f90` algorithm models such light-matter interactions. The result should look like below and as discussed in Fig.(5) of the referenced article.

![alt text](https://github.com/hjooya/Chemical-Theory-and-Computation/blob/main/Graph%20Theoretical%20Approach%20to%20Multiphoton%20Absorption%20Spectra/Image_Bidirectional_Coupling.jpg)


### Conclusion
In summary, a generalized graph theoretical method is introduced to investigate some of the characteristic multiphoton resonance processes and quantum interference phenomena for the superconducting qubit systems driven by intense ac fields. The various interacting designs at arbitrary coupling strengths are modeled by different graph products on colored-weighted graphs. The intuitive picture provided by this beyond rotating-wave approximation approach helps us to demonstrate some characteristic features of the superconducting qubit systems. 
