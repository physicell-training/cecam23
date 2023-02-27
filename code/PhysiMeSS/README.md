# PhysiMESS Directory
PhysiMESS (PhysiCell MicroEnvironment Structure Simulation) is a PhysiCell add-on which allows users to simulate ECM components as agents. 

This source code is for the CECAM-Lorentz workshop "The Extracellular Matrix: How to model structure complexity"

## Directory files
* Source code based on PhysiCell v1.10.1 with PhysiMESS modifications (see DOCUMENTATION.txt)
* Makefile
* config directory with pre-loaded examples
* [setup guide](https://github.com/physicell-training/cecam23/tree/main/code/PhysiMeSS/setup/PhysiMESS - Model Builder.pdf)


## Pre-loaded examples

### Fibre_Initialisation
* mymodel_initialisatin.xml and initialfibres.csv files for initialising fibres in the domain

### Fibre_Degradation 
* mymodel_fibre_degradation.xml and cells_and_fibres_attractant.csv to model one cell degrading fibres to reach attractant
* mymodel_matrix_degradation.xml and cells_and_fibres.csv to model growth of cell mass degrading matrix

### Cell_Fibre_Mechanics
* mymodel_fibremaze.xml and fibre_maze.csv to model cell moving around a maze made of fibres
* mymodel_potentials.xml and snowplough.csv to model both fibre pushing and rotation by cells
