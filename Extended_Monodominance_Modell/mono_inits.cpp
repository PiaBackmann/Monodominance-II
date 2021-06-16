/*mono_inits.cpp     June 2018 - January 2020 by Pia Backmann   */

/** Initialisation of all structs, the grid etc. **/

#include "mainwindow.h"
#include "iostream" //if console output is needed
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>


using namespace std;


World_t *new_World(int worldsize, int num_species, double mortality, int seed_mass, int seed, char *pfad, int radius, double prop_factor, int density_radius, int init_option, int self_replacement, int geometric_replacement, double shade, int JC_radius, int JC_maxcount, char *ende)
/** **************************************** **/
/** default constructor of world-class       **/
/**  PARAMETERS: world_t*                    **/
/**  RETURNS: (nothing) changes everything   **/
/** 			via pointers         **/
/**  params:  1: number trees                **/
/**      2: world size  3: mortality-rate    **/
/** **************************************** **/
{
    World_t *welt;
    Patch_t ***world_grid;
    double ***startgrid;
    elem_t *start;
    int i, j;
    welt = (World_t *) malloc(sizeof(World_t));
    world_grid = (Patch_t ***) malloc(worldsize*worldsize*sizeof(Patch_t**));
    for (i = 0; i< worldsize; i++)
    {
        world_grid[i] = (Patch_t **) malloc(worldsize*sizeof(Patch_t*));
    }

    startgrid = (double***)malloc(worldsize*sizeof(double**));
    for(i = 0; i < worldsize; i++)
    {
          startgrid[i] = (double**)malloc(worldsize*sizeof(double*));
          for(j = 0; j < worldsize; j++)
   	  {
               startgrid[i][j] = (double*)malloc(num_species*sizeof(double));
    	  }
    }
    start = NULL;
    if(welt == NULL)
    {
        fprintf(stderr, "out of Memory in new_World");
        cout << "war wohl nix";
        exit(1);
    }
    welt->self_replacement = self_replacement;
    welt->geometric_replacement = geometric_replacement;
    welt->radius = radius;
    welt->init_option = init_option;
    welt->real_propagule = prop_factor*seed_mass;
    welt->pfad = pfad;
    welt->ende = ende;
    welt->empty_patches = start;
    welt->seed = seed;
    welt->size = worldsize;                           // default 512 x 512 grid-cells
    welt->tick = 1;
    welt->mortality_rate = mortality;                 // default 0.0150000
    welt->disturbance_chance =  0.000000;             // not used now
    welt-> disturbance_range = 5.000000;              // not used now
    welt->disturbance_severity = 60.000000;           // not used now
    welt->shadetolerance_neutral = 0.4;
    welt->shadetolerance_mono = 1.0;
    welt->seed_mass = seed_mass;                      // default: 9450
    welt->num_species = num_species;                  // default 8: 1 monodominant, 7  identical species
    welt->shade = shade;
    welt->JCradius = JC_radius;
    welt->JC_maxcount = JC_maxcount;
    welt->clusteroutput[0] = 1;
    welt->density_radius = density_radius;
    j = 1;
    for(i = 100; i <= 10000; i = i+100) // normalerweise in 500er schritten
    {
        welt->clusteroutput[j] = i;
        j++;
    }
    startgrid = initialize_seeds(welt, startgrid);
    // Fülle das grid der Welt mit patches: (erstmal alle gleich)
    for(i = 0; i < worldsize; i++)
    {
        for(j = 0; j < worldsize; j++)
        {
	    // setze alle patches der Welt erstmal auf default:
            world_grid[i][j] = new_Patch(i,j, startgrid[i][j], welt);
        }
    }
    welt->grid = world_grid;
    welt->maxdispersal = radius;
    for(i = 0; i < worldsize; i++)
    {
          for(j = 0; j < worldsize; j++)
   	  {
               free(startgrid[i][j]);
    	  }
    }
    for(i = 0; i < worldsize; i++)
    {
          free(startgrid[i]);
    }
    free(startgrid);
    return(welt);
}

double*** initialize_seeds(World_t* welt, double*** startgrid)
{
    int i, j, k, counter, clustermax, center;
	double factor_mono = 0.0;


	// first: all random init: [init_option == 0]
	for(i = 0; i < welt->size; i++)
	{
		for(j = 0; j < welt->size; j++)
		{
			for(k = 0; k < welt->num_species; k++)
			{
				startgrid[i][j][k] = 12.5;  // all species (without monodom) have
									   // the same proportion within the seedbanks
			}
		}
	}
	// other initializations:

	if(welt->init_option == -1) //streets/lines of monodominance every third row [init_opt == -1]
	{

		factor_mono = 4.2;  // 7*size(einheitszelle)/(8-size(einheitszelle))
		for(i = 0; i < welt->size; i++)
		{
			for(j = 0; j < welt->size; j++)
			{
				startgrid[i][j][0] = 0;
				if(j%3 == 0)
				{
					startgrid[i][j][0] = 12.5*factor_mono;
				}
			}
		}
	}
	if(welt->init_option == 1) // checkerboard configuration
	{

		factor_mono = 7.0;
		for(i = 0; i < welt->size; i++)
		{
			for(j = 0; j < welt->size; j++)
			{
				startgrid[i][j][0] = 0;
				if((i%(welt->init_option + 1) == 0) && (j%(welt->init_option + 1) == 0))
				{
					startgrid[i][j][0] = 12.5*factor_mono;
				}
			}
		}
	}
	if(welt->init_option == 3) // knight movement configuration
	{

		factor_mono = 11.666667;
		for(i = 0; i < welt->size; i++)
		{
			for(j = 0; j < welt->size; j++)
			{
				startgrid[i][j][0] = 0;

				if(i == 0 && j == 0)
				{
					startgrid[i][j][0] = 12.5*factor_mono;
				}
				else if(((j+5) % 5 == 0) && ((i+5)%5 == 0))
				{
					startgrid[i][j][0] = 12.5*factor_mono;
				}
				else if (i > 0 && j > 0)
				{
					if(((i+4)%5 == 0) && ((j+3)%5 == 0))
					{
						startgrid[i][j][0] = 12.5*factor_mono;
					}
					else if(((i+3)%5 == 0) && ((j+1)%5 == 0))
					{
						startgrid[i][j][0] = 12.5*factor_mono;
					}
					else if(((i+2)%5 == 0) && ((j+4)%5 == 0))
					{
						startgrid[i][j][0] = 12.5*factor_mono;
					}
					else if(((i+1)%5 == 0) && ((j+2)%5 == 0))
					{
						startgrid[i][j][0] = 12.5*factor_mono;
					}
				}
			}
		}
	}
	if(welt->init_option == 4) // one big cluster in the middle
    {
        for(i = 0; i < welt->size; i++)
		{
			for(j = 0; j < welt->size; j++)
			{
				startgrid[i][j][0] = 0;  // No Monodom species at all
			}
		}
		clustermax = 15;
		center = round(welt->size/2);
		j = center;
		i = center;
		//startgrid[i][j][0] = 1000000;
	//	for(int d = -clustermax; d <= clustermax; d++) // make one square
	    for(int d = 0; d< welt->size; d++)
		{
			//for(int f = -clustermax; f <= clustermax; f++)
			for(int f = 0; f < welt->size; f++)
			{
				if(calculate_distance(center, center, d,f) <= clustermax)
				{
				    startgrid[d][f][0] = 1000000;
				}

			}
		}
    }
	return(startgrid);
}


Tree_t *new_Tree(int species, int seed_mass, int radius)
{
    Tree_t* tree;
    tree = (Tree_t*)malloc(sizeof(Tree_t));
    tree->species = species;
    if( species == 0)                                     // monodominant
    {
        tree->dispersalrange = 1;
        tree->dispersalarea = 0;
        tree->seednumberperpatch = seed_mass;   // seed mass of monodominant species, default: 9450
        tree->shadetolerance = 0.9;
        tree->status = 1;                                // 1 = "juvenile"; 2 = "adult"
    }
    else if( species >= 1)                               // other species
    {
        tree->dispersalrange = radius;
        tree->dispersalarea = 0;
        tree->seednumberperpatch = 20;
        tree->shadetolerance = 0.2;
        tree->status = 1;                                // 1 = "juvenile"; 2 = "adult"
    }
    return(tree);
}


Patch_t* new_Patch(int i, int j, double* startseeds, World_t* welt)
/** **************************************** **/
/** Patch constructor                        **/
/**  PARAMETERS: positions i,j               **/
/**  RETURNS: nothing                        **/
/** 			                     **/
/**  params:                                 **/
/** **************************************** **/
{
    Patch_t* patch;
    patch = (Patch_t*) malloc(sizeof(Patch_t));
    patch->seed_array = (double*) malloc(welt->num_species*sizeof(double));
    patch->former_species = 1;
    patch->color = 0;                                 // in the beginning are all grid cells brown.
    for(int l = 0; l < welt->num_species; l++)
    {
        patch->seed_array[l] = startseeds[l];   //erste Verteilung alle gleich initialisieren
    }
    patch->pos_i = i;
    patch->pos_j = j;
    patch->current_tree = NULL;
    patch->shade = 0.0;

    return(patch);
}




int *make_local_grid(World_t *welt)
{
    int neighbourhood_count, i,j,o,p;
    int *local_density;
    local_density = (int*) malloc(welt->size*welt->size*sizeof(int));
    int density_radius = welt->density_radius;
    for(i = 0; i < welt->size; i++)
    {
        for(j = 0; j < welt->size; j++)
        {
            neighbourhood_count = 0;
            // iterate over all patches around focal grid cell
            for(int k = i-density_radius; k <= i+density_radius; k++)
            {
                for(int l = j-density_radius; l <= j+density_radius; l++)
                {
                    o = k % welt->size;
                    p = l % welt->size;
                    // 1. : Flip coordinates at the borders (periodic boundary conditions)

                                                                               // check periodic boundaries
                    if(o < 0){o = welt->size + o;}
                    if(o >= welt->size){o = o - welt->size;}
                    if(p < 0){p = welt->size + p;}
                    if(p >= welt->size){p = p - welt->size;}
                    if(welt->grid[o][p]->color == 1)                           // iterate over all monodominant trees within max-dispersal range
                    {

                       if (calculate_distance(i, j, k, l) <= density_radius)   // dispersal area is circular, not quadratic
                       {
                           neighbourhood_count++;
                       }
                    }
               }
            }
            local_density[j+i*welt->size] = neighbourhood_count;               // switch i und j if necessary
        }
    }
    return(local_density);
}
