/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON1 ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"


void stiffness_substrate_calculation(double dt)
{
	//wheres waldo
	static int collagen_index = microenvironment.find_density_index("collagen");
	static int crosslinks_index = microenvironment.find_density_index("crosslinks");
	static int stiffness_index = microenvironment.find_density_index("stiffness");

		//iterate thorugh voxels
			//read values 
			//calculate
			//write
		for (int n=0; n<microenvironment.mesh.voxels.size(); n++)
		{
			std::vector<double> rho = microenvironment(n);
			double c = rho[ collagen_index ];
			double cl = rho[ crosslinks_index ];
			double s = rho[ stiffness_index ];

			double dc = c * cl;

			microenvironment(n)[stiffness_index] = dc;	
		}
	return;
}

void degrade_matrix_soluble_mmp()
{
	// Get the index of the MMP and Collagen substrate
	int collagen_index = microenvironment.find_density_index("collagen");
	int mmp_index = microenvironment.find_density_index("mmp");

	//Finding the crosslinks substrate index
	int crosslinks_index = microenvironment.find_density_index("crosslinks");

	// Loop over all the voxel in the domain
	#pragma omp parallel for
	for (int voxel_index = 0; voxel_index<microenvironment.number_of_voxels(); voxel_index++)
	{	
		// Get the density of Collagen and MMP in the current voxel
		double collagen_density = microenvironment.density_vector(voxel_index)[collagen_index];
		double mmp_density = microenvironment.density_vector(voxel_index)[mmp_index];
		double crosslinking_quantity = microenvironment.density_vector(voxel_index)[crosslinks_index];

		// Compute the degradation rate depending on the user parameter and
		// the density of both Collagen and MMP
		double degradation_rate_soluble_mmp = parameters.doubles("degradation_rate_soluble_mmp") * mmp_density * (1.1 - crosslinking_quantity) * (1.1 - collagen_density);
		
		// Compute the new collagen density to simulate degradation
		collagen_density *= (1 - degradation_rate_soluble_mmp); 
		
		// Update the current voxel with the new collagen density
		microenvironment.density_vector(voxel_index)[collagen_index] = collagen_density;
	}
}

void custom_cell(Cell* pCell, Phenotype& phenotype, double dt)
{
	if (pCell->phenotype.motility.is_motile)
	{
		// Scaling our orientation to the radius of the cell
		std::vector<double> scaled_orientation = phenotype.geometry.radius * pCell->state.orientation;
		
		// Calculating the position of the membrane in the direction of the cell
		std::vector<double> position_membrane = pCell->position + scaled_orientation;
		
		// Computing the index of the voxel at that position
		int voxel_membrane = microenvironment.nearest_voxel_index( position_membrane );
		
		// Finding the collagen substrate index
		int collagen_index = microenvironment.find_density_index("collagen");
		
		// Finding out the quantity of collagen in that voxel
		double collagen_quantity = microenvironment.density_vector(voxel_membrane)[collagen_index];

		//Finding the crosslinks substrate index
		int crosslinks_index = microenvironment.find_density_index("crosslinks");

		// Finding out the crosslinking average in that voxel
		double crosslinking_quantity = microenvironment.density_vector(voxel_membrane)[crosslinks_index];
		
		// Reducing the quantity of collagen by 2% everytime we call this function (every mechanics_dt)
		microenvironment.density_vector(voxel_membrane)[collagen_index] *=  (1.0 - parameters.doubles("degradation_rate_membrane_mmp")*(1.1 - crosslinking_quantity)); 
	}
}

void custom_update_cell_velocity( Cell* pCell, Phenotype& phenotype, double dt)
{
	///////////////// Update migration speed //////////////////////////
	// Scaling our orientation to the radius of the cell
	std::vector<double> scaled_orientation = phenotype.geometry.radius * pCell->state.orientation;
	
	// Calculating the position of the membrane in the direction of the cell
	std::vector<double> position_membrane = pCell->position + scaled_orientation;
	
	// Computing the index of the voxel at that position
	int voxel_membrane = microenvironment.nearest_voxel_index( position_membrane );

	// Finding the collagen substrate index
	int collagen_index = microenvironment.find_density_index("collagen");
	
	// Finding out the quantity of collagen in that voxel
	double collagen_quantity = microenvironment.density_vector(voxel_membrane)[collagen_index];

	// Base migration speed
	double base_migration_speed = get_single_base_behavior(pCell, "migration speed") ;
	// std::cout << "Base migration speed: " << base_migration_speed << std::endl;

	// Calculate collagen dependent speed
	double collagen_dependent_speed = 4 * collagen_quantity * (1-collagen_quantity);

	pCell->phenotype.motility.migration_speed = collagen_dependent_speed * base_migration_speed;
	// std::cout << "New speed: " << pCell->phenotype.motility.migration_speed << std::endl;

	///////////////////////// Update standard velocity /////////////////////////////
	
	// Calling the standard update velocity of PhysiCell
	standard_update_cell_velocity(pCell, phenotype, dt);
	
	// Setting our orientation vector from the velocity vector
	pCell->state.orientation = normalize(pCell->velocity);
	
	return; 
}

void create_cell_types( void )
{
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
	
	cell_defaults.functions.custom_cell_rule = custom_cell; 
	cell_defaults.functions.update_velocity = custom_update_cell_velocity;
		
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 
	
	// create some of each type of cell 
	
	Cell* pC;
	
	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k]; 
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
		for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
		{
			std::vector<double> position = {0,0,0}; 
			// position[0] = Xmin + UniformRandom()*Xrange; 
			// position[1] = Ymin + UniformRandom()*Yrange; 
			// position[2] = Zmin + UniformRandom()*Zrange; 
			
			pC = create_cell( *pCD ); 
			pC->assign_position( position );
		}
	}
	std::cout << std::endl; 
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml(); 	
	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	// If the cell is dead, we don't need to update it anymore (we remove the update function)
	if (pCell->phenotype.death.dead){
		pCell->functions.update_phenotype = NULL;
		return;
	}

	// Update the cell division and apoptosis rate depending on the extracellular oxygen concentration
    update_cell_and_death_parameters_O2_based(pCell,phenotype,dt);

	// If we don't want to simulate the EMT transition we don't need to update the motility of the cell
	if (parameters.bools("emt_transition"))
	{
		// Get oxygen substrate index
		int oxygen_index = microenvironment.find_density_index("oxygen");

		// Get the index of the voxel the cell is in
		int voxel_index = pCell->get_current_voxel_index();

		// Get extra-cellular oxygen density (in the voxel)
		double oxygen_density = microenvironment.density_vector(voxel_index)[oxygen_index];

		// Update the motility and the degradation by MMPs of the cell depending on the
		// oxygen density
		if(oxygen_density < pCell->parameters.o2_hypoxic_response)
		{
			pCell->phenotype.motility.is_motile = true;
			pCell->phenotype.secretion.secretion_rate("mmp") = parameters.doubles("secretion_rate_soluble_mmp");
		} else
		{
			pCell->phenotype.motility.is_motile = false;
			pCell->phenotype.secretion.secretion_rate("mmp") = 0.0;
		}
	}

	return;
}

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 
