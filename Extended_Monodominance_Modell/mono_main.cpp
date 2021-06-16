/*mono_main.cpp     June 2018 - January 2020 by Pia Backmann   */
/** Hier sind alle Fensterspezifischen Funktionen **/

#include "mainwindow.h"
#include "iostream"     //if console output is needed
#include<time.h>        /* time_t, struct tm, difftime, time, mktime */
#include<string.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

using namespace std;

int main(int argc, char *argv[])
{
    // setze default params:
    //ijmax = 512;  // Weltgröße
    //number_pfts = 8;  // hier später dynamisch ändern lassen...
    //mortality_rate = 0.015; // Sterblichkeit der Bäume
    //neighbourhood = 0; // 0 = "in-radius" 1 = "Manhattan"
    int num_cluster, ijmax = 512, number_pfts = 8, seed_mass = 9450, i = 1,  seed = 1, radius = 20, density_radius = 1, random_init = 0, self_replacement = 0, geometric_replacement = 0,  init_option = 0, JC_radius = 0;
    double prop_factor = 1.0, shade = 0.0;
    int* one_dim_grid;
    int* local_grid;
    int* next; int* cluster;
    int pos = 0;
    World_t *welt;
    double prop_radius[100];
    char *cluster_file = (char*)malloc(1000*sizeof(char));
    char *histo_file= (char*)malloc(1000*sizeof(char));
    FILE *out; FILE *histogramm;
    char *pfad_allgemein, *pfad_ende;
    pfad_allgemein = (char*)malloc(1000*sizeof(char));
    pfad_ende = (char*)malloc(1000*sizeof(char));
    prop_radius[1] =  5.0/1257.0;
    prop_radius[2] = (13.0/1257.0);
    prop_radius[3] = (29.0/1257.0);
    prop_radius[4] = (49.0/1257.0);
    prop_radius[5] = (81.0/1257.0);
    prop_radius[6] = (111.0/1257.0);
    prop_radius[7] = (149.0/1257.0);
    prop_radius[8] = (199.0/1257.0);
    prop_radius[9] = (253.0/1257.0);
    prop_radius[10] = (317.0/1257.0);
    prop_radius[11] = (377.0/1257.0);
    prop_radius[12] = (441.0/1257.0);
    prop_radius[13] = (529.0/1257.0);
    prop_radius[14] = (613.0/1257.0);
    prop_radius[15] = (709.0/1257.0);
    prop_radius[16] = (797.0/1257.0);
    prop_radius[17] = (901.0/1257.0);
    prop_radius[18] = (1009.0/1257.0);
    prop_radius[19] = 1129.0/1257.0;
    prop_radius[20] = 1257.0/1257.0;
    prop_radius[21] = 1373.0/1257.0;
    prop_radius[22] = 1515.0/1257.0;
    prop_radius[23] = 1653.0/1257.0;
    prop_radius[24] = 1793.0/1257.0;
    prop_radius[25] = 1961.0/1257.0;
    prop_radius[26] = 2121.0/1257.0;
    prop_radius[27] = 2289.0/1257.0;
    prop_radius[28] = 2453.0/1257.0;
    prop_radius[29] = 2629.0/1257.0;
    prop_radius[30] = 2821.0/1257.0;
    prop_radius[40] = 5025.0/1257.0;
    prop_radius[50] = 7845.0/1257.0;
    prop_radius[60] = 11289.0/1257.0;


    double mortality_rate = 0.015; // mortality of trees
     while( i<argc)
     {
              if(strcmp(argv[i],"--size") == 0)
              {
                  ijmax = atoi(argv[++i]);
              }
              else if(strcmp(argv[i],"--num_species") == 0)
              {
                  number_pfts = atoi(argv[++i]);
              }
              else if(strcmp(argv[i],"--mort_rate") == 0)
              {
                  mortality_rate = atof(argv[++i]);
              }
              else if(strcmp(argv[i],"--neighbourhood") == 0)
              {
                 // neighbourhood = atoi(argv[++i]);
              }
              else if(strcmp(argv[i],"--seed_mass") == 0)
              {
                 seed_mass = atoi(argv[++i]);
              }
              else if(strcmp(argv[i],"--seed") == 0)
              {
                 seed = atoi(argv[++i]);
              }
              else if(strcmp(argv[i],"--radius") == 0)
              {
                radius = atoi(argv[++i]);
              }
              else if(strcmp(argv[i],"--densrad") == 0)
              {
                density_radius = atoi(argv[++i]);
              }
              else if(strcmp(argv[i],"--random_init") == 0)
              {
                random_init = atoi(argv[++i]);
              }
              else if(strcmp(argv[i],"--init_option") == 0)
              {
                init_option = atoi(argv[++i]);
              }
              else if(strcmp(argv[i],"--self_replacement") == 0) // 0 = no, 1 = yes
              {
                self_replacement = atoi(argv[++i]);
              }
              else if(strcmp(argv[i],"--geometric_replacement") == 0) // 0 = no, 1 = yes
              {
                geometric_replacement = atoi(argv[++i]);  // 0 = all neighbours get seeds; 1 = mono seeds are distributed randomly
              }
              else if(strcmp(argv[i],"--shade") == 0) //[0.0 (no shade) ... 1.0 (max)]
              {
                shade = atof(argv[++i]);  // 0 = all neighbours get seeds; 1 = mono seeds are distributed randomly
              }
              else if(strcmp(argv[i],"--JC_radius") == 0) // 0 = no, 1 = yes
              {
                JC_radius = atoi(argv[++i]);  // 0 = all neighbours get seeds; 1 = mono seeds are distributed randomly
              }
              else
              {
                  i++;
              }
      }
      prop_factor = prop_radius[radius];
      if (self_replacement == 1)
      {
            prop_factor = prop_radius[radius]*4/5; //* 4/5 because of self replacement
      }
      int JC_maxcount =  prop_radius[JC_radius]*1257.0; // geht bis von 1-30 und 40 50 60
    /************************************* FILE Routines:  ***********************************************************************/

    char* myfileroot="";
    sprintf(pfad_allgemein, "%sData/" , myfileroot); // kein slash am anfang -> ab working directory
  //  sprintf(pfad_ende, "_size_%d_seedmass_%d_radius_%d_mort_%f_init_option_%d_seed_%d", ijmax, seed_mass, radius, mortality_rate, init_option, seed);
// without seed
    sprintf(pfad_ende, "_size_%d_seedmass_%d_radius_%d_mort_%f_init_option_%d_self_replacement_%d_geometric_replacement_%d_shade_%1.2f_JC_radius_%d", ijmax, seed_mass, radius, mortality_rate, init_option, self_replacement, geometric_replacement, shade, JC_radius);
    sprintf(cluster_file, "%sStats/Cluster_stats%s.out",pfad_allgemein, pfad_ende);
    out = fopen(cluster_file, "w");
    sprintf(histo_file, "%sClusterhistogramme/ClusterHistogramm%s.out", pfad_allgemein, pfad_ende);
    histogramm = fopen(histo_file, "a");

     /****************************************************************************************************************************/
     /*  Eingabeparameter  */
     next = (int*) malloc(4*ijmax*ijmax*sizeof(int));
     srand48(seed); //random seed setzen

    // double test = drand48();
     //printf( "seed = %d \t zufallszahl 1 = %f \n", seed, test);
     /***************      INTITIALISATION      *****************/
     welt = new_World(ijmax, number_pfts, mortality_rate, seed_mass, seed, pfad_allgemein, radius, prop_factor, density_radius, init_option, self_replacement, geometric_replacement, shade, JC_radius, JC_maxcount, pfad_ende);

     distribute_seeds(welt); // verteile bäume nach lotterie

     one_dim_grid = (int*) malloc(welt->size*welt->size*sizeof(int));
     local_grid = (int*) malloc(welt->size*welt->size*sizeof(int));
     cluster = (int*) malloc(welt->size*welt->size*sizeof(int));
    /******************************************************************************************************************************/
    /**********  Routines  **********/
    fill_neighbour_grid(next, welt->size);
    while(welt->tick <= 10000)
    {


        // if( welt->tick == 1 || welt->tick == 5000 || welt->tick == 10000)
        // {
	   //  two_to_one_dim(welt, one_dim_grid);
	  //   num_cluster = percol_cluster(4,welt->size*welt->size, next, one_dim_grid, cluster, welt, out, histogramm);
            //local_grid = make_local_grid(welt);
	 //    two_to_one_dim(welt, one_dim_grid);
          //   print_one_dim_grid_file(1, welt, cluster);
       //  }

        if (welt->tick == welt->clusteroutput[pos])
        {
	         // clusteralgorithmus
	         two_to_one_dim(welt, one_dim_grid);
	         num_cluster = percol_cluster(4,welt->size*welt->size, next, one_dim_grid, cluster, welt, out, histogramm);
	         print_one_dim_grid_file(1, welt, cluster);
	         pos++;
	         printf("Clustercount = %d\n", welt->tick);

        }
        if(welt->tick % 5000 == 0 || welt->tick == 1)
        {
	       // cout << "Tick = " << welt->tick << endl;
            print_world(welt);
            //print_shade(welt);
           // local_grid = make_local_grid(welt);
	  //  two_to_one_dim(welt, one_dim_grid);
      //      print_one_dim_grid_file(2, welt, cluster);
	   // print_recruitment_wahrscheinlichkeit(welt);
           // print_recruitment_wahrscheinlichkeit(welt);
	        two_to_one_dim(welt, one_dim_grid);
            print_one_dim_grid_file(1, welt, cluster);
        }
       /* if(welt->tick % 10 == 0 || welt->tick < 500)
        {
            print_world(welt);
        }*/
	 welt = go(welt);
    }
    /**************** FREE MEMORY **********************/
    fclose(out);
    fclose(histogramm);
    kill_World(welt);
    free(cluster_file);
    free(histo_file);
    free(pfad_allgemein);
    free(pfad_ende);
    free(one_dim_grid);
    free(local_grid);
    free(next);
    free(cluster);
    return(0);
}



