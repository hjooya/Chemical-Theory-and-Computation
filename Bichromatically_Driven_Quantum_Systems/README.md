## Absorption Spectroscopy of Bichromatically Driven Quantum Systems 

[Hossein Jooya, *et. al.*, Phys. Rev. B 96, 174518](https://kuscholarworks.ku.edu/bitstream/handle/1808/27228/Pan_2017.pdf?sequence=1)

<p align="center">
<img src="https://github.com/hjooya/Chemical-Theory-and-Computation/blob/main/Bichromatically_Driven_Quantum_Systems/Bichromatic_Driven_TLS.gif" width="500" height="250"/>
</p>

This example shows how to simulate the observation of two distinct quantum interference patterns in the absorption spectra
when a transmon superconducor qubit is subjected to a bichromatic microwave field with the same Rabi frequencies. The calculations are done within the two-mode Floquet formalism with no dissipation processes. Please refer to the above paper for details.

Such multiphoton interference may be used in manipulating the quantum states of various quantum systems. For example see the excitation spectra of a quantum dot in Fig.(2) of [PHYSICAL REVIEW B 89, 155305 (2014)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.89.155305).


### Usage

Simply compile the Makefile, to run the `Bichromatically_Driven_TLS.f90` and obtain the transition probabilities and corresponding quasi-energies. The animated image presented is the two-mode Floquet calculations for the evolution of a two-level-system (TLS) subject to a bichromatic field, when the intensity of the second field is gradually increased until it becomes equal to the intensity of the first one.  



