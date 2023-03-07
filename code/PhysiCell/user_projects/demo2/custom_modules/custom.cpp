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
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"

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

	Cell_Definition* pCD = find_cell_definition( "epithelial"); 
	pCD->functions.update_phenotype = epithelial_phenotype;
	pCD->functions.custom_cell_rule = epithelial_custom;

	pCD = find_cell_definition( "macrophage"); 
	pCD->functions.update_phenotype = macrophage_phenotype; 
	
	pCD = find_cell_definition( "fibroblast"); 
	pCD->functions.update_phenotype = fibroblast_phenotype; 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
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
			position[0] = Xmin + UniformRandom()*Xrange; 
			position[1] = Ymin + UniformRandom()*Yrange; 
			position[2] = Zmin + UniformRandom()*Zrange; 
			
			pC = create_cell( *pCD ); 
			pC->assign_position( position );
		}
	}
	std::cout << std::endl; 
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml(); 	

	for( int n =0; n < (*all_cells).size() ; n++ )
	{
		Cell* pC = (*all_cells)[n]; 
		if( fabs( pC->position[0] ) < 50 && fabs( pC->position[1] ) < 50 )
		{ set_single_behavior( pC , "necrosis" , 9e99 ); }
	}
	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 

void epithelial_phenotype( Cell* pCell , Phenotype& phenotype , double dt )
{
	// if apoptotic, release apoptotic debris and exit 
	if( get_single_signal( pCell, "apoptotic") > 0.5 )
	{
		set_single_behavior( pCell, "apoptotic debris secretion" , 1); 
		pCell->functions.update_phenotype = NULL; 
		return; 
	}

	// if necrotic, release necrotic debris and exit 
	if( get_single_signal( pCell, "necrotic") > 0.5 )
	{
		set_single_behavior( pCell, "necrotic debris secretion" , 1); 
		pCell->functions.update_phenotype = NULL; 
		return; 
	}

	// get signals 
	double p = get_single_signal( pCell , "pressure"); 
	double p_max = 1.0; 
	double nd = get_single_signal( pCell , "necrotic debris"); 
	double f = get_single_signal( pCell , "fibrosis"); 

	// get (nonzero) base values, and needed constants
	double b0 = get_single_base_behavior( pCell , "cycle entry"); 

	double dA0 = get_single_base_behavior( pCell , "apoptosis"); 
	double dAM = 100 * dA0; 

	bool attached = (bool) get_single_signal( pCell , "custom:attached");
	double rA = get_single_behavior( pCell ,"custom:attachment_rate"); 


	// calculate and set values 
		// birth
	p /= p_max; 
	if( p > 1 )
	{ p = 1; }
	double b = b0 * (1-p);
	set_single_behavior( pCell , "cycle entry" , b ); 

		// apoptosis 
	double dA = dA0 + (dAM-dA0) * Hill_response_function( nd , 0.1 , 4 ); 
	set_single_behavior( pCell , "apoptosis" , dA ); 

		// dynamic attachment to ECM 
	if( attached == false )
	{
		double prob_attach = rA * Hill_response_function( f , 0.5 , 4 ) * dt ; 
		if( UniformRandom() <= prob_attach )
		{
			set_single_behavior( pCell , "custom:attached_x" , pCell->position[0] ); 
			set_single_behavior( pCell , "custom:attached_y" , pCell->position[1] ); 
			set_single_behavior( pCell , "custom:attached_z" , pCell->position[2] ); 

			set_single_behavior( pCell , "custom:attached" , 1.0 ); 
		}
	}


}

void epithelial_custom( Cell* pCell, Phenotype& phenotype, double dt )
{
	// exit if not attached 
	if( get_single_signal( pCell , "custom:attached") < 0.5 )
	{ return; }

	// get constants 
	double k = get_single_behavior( pCell ,"custom:spring_k"); 
	double attached_x = get_single_signal( pCell , "custom:attached_x");
	double attached_y = get_single_signal( pCell , "custom:attached_y");
	double attached_z = get_single_signal( pCell , "custom:attached_z");

	// calculate stretch and dv 
	std::vector<double> dv = { attached_x , attached_y , attached_z }; // (xA) 
	dv -= pCell->position; // (xA-x)
	dv *= k; // k*(xA-x)

	pCell->velocity += dv; 

	return; 
}

void macrophage_phenotype( Cell* pCell , Phenotype& phenotype , double dt )
{
	// get signals 
	double nd = get_single_signal( pCell , "necrotic debris");
	double vol = get_single_signal( pCell , "volume"); 

	// get ref values 
	double s0 = get_single_base_behavior( pCell , "pro-inflammatory factor secretion"); 
	double sM = 10; //

	double ph0 = get_single_base_behavior( pCell , "phagocytose dead cell"); 
	double ph = ph0;  
	ph = ph0; 
	if( vol > 2500 )
	{ ph = 0; }
	set_single_behavior( pCell , "phagocytose dead cell" , ph ); 

	// calculate and set 
	double s = s0 + (sM-s0)*Hill_response_function( nd , 0.5 , 2 ); 
	set_single_behavior( pCell , "pro-inflammatory factor secretion" , s ); 

	return; 
}

void fibroblast_phenotype( Cell* pCell , Phenotype& phenotype , double dt )
{
	// get signals 
	double pif = get_single_signal( pCell , "pro-inflammatory factor");

	// get ref values 
	double s0 = get_single_base_behavior( pCell , "fibrosis secretion"); 
	double sM = 1; //

	// calculate and set 
	double s = s0 + (sM-s0)*Hill_response_function( pif , 0.5 , 3 ); // 4 
	set_single_behavior( pCell , "fibrosis secretion" , s ); 

	return; 
}