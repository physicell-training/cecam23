# PhysiMESS Modifications
The following files have been modified for PhysiMESS:

Code read between lines:

     // !!! PHYSIMESS CODE BLOCK START !!! //
     // !!! PHYSIMESS CODE BLOCK END !!! //


## PhysiMESS/modules/PhysiCell_pathology.cpp
* Code to visualise fibres as thin rectangles rather than spheres added to ```SVG_plot```
* Similarly fibres visualised as rectangles in the plot legend via ```create_plot_legend```

## PhysiMESS/custom-modules/custom.cpp
* Changes made in ```set_up_tissue``` to allow us to initialise fibres either manually or from a csv file 

## PhysiMESS/core/PhysiCell_cell.h

    #include <list> required 

### The following additional functions:
* ```std::vector<Cell*> crosslinkers;```
* ```std::vector<double> crosslink_point;```
* ```void force_update_motility_vector(double dt_);```
* ```void check_fibre_crosslinks(Cell*);```
* ```void degrade_fibre(Cell*);```
* ```std::vector<double> nearest_point_on_fibre(std::vector<double> point, Cell* , std::vector<double>& displacement);```
* ```std::vector<double> CrossProduct(std::vector<double> vector_A, std::vector<double> vector_B, std::vector<double>& C_P);```
* ```double DotProduct(std::vector<double> vector_A, std::vector<double> vector_B);```
* ```std::list<int> register_fibre_voxels( Cell* pCell );```
* ```void deregister_fibre_voxels( Cell* pCell );```
* ```std::list<int> find_agent_voxels(Cell * pCell );```
* ```void find_agent_neighbors( Cell* pCell );```
* ```void add_crosslinks( Cell* pCell );```

### The following additional parameters (XML):
* ```double mLength = 0;```
* ```double mRadius = 0;```
* ```int fail_count = 0;```
* ```double stuck_counter = 0;```
* ```double unstuck_counter = 0;```
* ```bool degradation_flag = false;```
* ```bool fibre_degradation = false;```
* ```double fibreDegradationRate = 0.0;```
* ```double stuck_threshold = 0.0;```
* ```bool fibre_rotation = false;```
* ```double mFibreStickiness = 1.0;```
* ```bool fibre_pushing = false;```
* ```int X_crosslink_count;```
* ```double mVelocityAdhesion = 0;```
* ```double mVelocityContact = 0;```
* ```double mCellVelocityMaximum= 0;```

## PhysiMESS/core/PhysiCell_cell.cpp
* In ```Cell_State: crosslinkers.resize(0);```
* Significant changes to the ```add_potentials``` function 
* Addition to ```update_position``` function to flag if a cell is stuck 
    
### New functions 
* ```check_fibre_crosslinks``` and ```add_crosslinks``` - determine fibre crosslinks
* ```CrossProduct``` and ```DotProduct``` - required for fibre crosslinks
* ```degrade_fibre``` - remove a fibre if it’s determined that it will be degraded 
* ```force_update_motility_vector``` - alters a cells non-force dictated movement if deemed to be stuck 
* ```nearest_point_on_fibre``` - calculate the displacement vector from a point to the nearest point on a fibre
* ```find_agent_neighbors``` and ```find_agent_voxels``` - determine neighbours of fibres and cells 
* ```register_fibre_voxels``` and ```deregister_fibre_voxels``` - add/remove fibres to/from the list of agents in each of its voxels


## PhysiMESS/core/PhysiCell_standard_models.cpp
* Substantial changes to ```standard_update_cell_velocity```

## PhysiMESS/core/PhysiCell_cell_container.cpp
* Substantially changes to ```update_all_cells``` the following functions are called in order: 
    
      register_fibre_voxels; find_agent_neighbors; deregister_fibre_voxels; add_crosslinks, degrade_fibre
    
* In function remove_agent_from_voxel change to the following lines:       
      
      agent_grid[voxel_index][delete_index] = agent_grid[voxel_index][agent_grid[voxel_index].size()-1 ];
      agent_grid[voxel_index].pop_back();


