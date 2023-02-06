## Psi4-MATLAB Molecular Dynamics Simulation Workflow  


<img src="https://github.com/hjooya/Chemical-Theory-and-Computation/blob/main/Psi4-MATLAB-MD-Simulations/Penicillamine_MD_Simulation.gif" width="400" height="200"/> <img src="https://github.com/hjooya/Chemical-Theory-and-Computation/blob/main/Psi4-MATLAB-MD-Simulations/HF_Energy_Plot.jpg" width="350" height="200" />


### Introduction
In this example, we show how [Psi4](https://psicode.org/), an open-source suite of ab initio quantum chemistry programs, is used for molecular dynamics (MD) simulations. This Psi-MATLAB workflow starts with a single molecular structure input, rotates it around a C-C bond, and calculates the molecular energy at the desired theory level at each time step.  The output of the Psi4 computations is then processed in MATLAB to extract data and to build a single .mat file for further analysis.  

### Usage 
Simply run the MATLAB live script 
> MATLAB_Psi4_MD_Simulation.mlx

### Setting up the conda environment
Psi4 can be installed from this [source](https://psicode.org/installs/v17/). This example is using conda environment setup for this software.

### Input molecule in sdf format and convert to xyz
We use [penicillamine](https://pubchem.ncbi.nlm.nih.gov/compound/Penicillamine) for this example. Penicillamine is a FDA approved chelating agent used to decrease copper stores in Wilson disease. The initial configuration is the distorted molecular geometry due to 60 degrees rotation around the C5-C6 single bond. The final optimized geometry is also visualized here using [Molviewer](https://www.mathworks.com/help/bioinfo/ref/molviewer.html) for comparison. 
> 'Convert_sdf_to_xyz.m' 

is a MATLAB function that converts the Structure-Data File (sdf) format to XYZ format for further processing.

## Set up the MD simulations
Refer to Psi4 documentation for [Single-Point Energy](https://psicode.org/psi4manual/master/energy.html) setup. In this example Hartree–Fock (HF) self consistent field (SCF) theory level is used with cc-pvtz [basis set](https://psicode.org/psi4manual/master/basissets.html). The MD simulations are done when the molecule is rotated around C5-C6 single bond (z-axis). Molviewer is used to identify the indices of the rotating atoms. The generated xyz and Python input files can be stored in dedicated folders, for which the paths should be given.

## Generate the MD configurations and running the Psi4 simulations
Since the C5-C6 bond is initialy set as the intermolecular z-axis, we can use MATLAB's built-in rotation matrix function, [rotz](https://www.mathworks.com/help/phased/ref/rotz.html), to generate the new configurations for this structural transformation. 
> 'ZRotate_Molecule.m' 

is the function that takes in the initial configuration saves the rotated geometry as xyz file in the xyzpath folder. 
> 'Psi4_py_input_builder.m' 

builds the python script to run the Psi4 calculations. The output is overwritten in the 'Psi4_Output.dat' datafile for further processes.
> 'Get_Atomic_Numbers.m' 

extracts the atomic numbers and atomic symbols from the Psi4 output. "Get_Atomic_Geometries" and 'Get_Total_Energy' extracts the atomic coordinates and the computed energy at each time step, respectively.

All extracted variables can then be saved in one '.mat' data file for future work. 

## Results
The above figure (top right panel) shows the total energy vs. the rotation angle. The [steric effect](https://www.sciencedirect.com/topics/chemistry/steric-effect) in the initial configuration explains the observation of the high starting energy during the MD simulation. This energy monotonically decreases when the conformational restriction is removed during the simulation, until the structure reaches its minimum energy.   

We can store the conformational structures during the simulation by calling the 'Trajcetory' function. The generated XYZ file can be then viewed with a molecular visualization program like [VMD](https://www.ks.uiuc.edu/Research/vmd/). The above movie (top left panel) shows this result.


