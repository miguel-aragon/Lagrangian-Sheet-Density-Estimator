/*

  TRIANGLE SOFTWARE RASTERIZER

  Triangle renderer shamelesly ripped from: http://devmaster.net/forums/topic/1145-advanced-rasterization/
  (at least I read the whole post...)


  A few notes:
   - Triangles must be in counter-clockwise order.
   - http://stackoverflow.com/questions/11184473/how-to-link-a-c-shared-library-with-gcc


  Compile (shared library), IDL complains if not compiled with openMP flags...

  PARALLEL (openMP)
   g++ -O3 -D_REENTRANT -fopenmp ltfe_soft_all_0.0.cpp -o ltfe_soft_all

Written by Miguel A. Aragon Calvo, 2013

*/
   
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
   
#include "ltfe_soft_all.h"

void *__gxx_personality_v0;

//===============================================================
//--- Global variables
//===============================================================
struct io_header_1 {
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  //--- Fills to 256 Bytes 
} header1;

struct particle_data {
  float  Pos[3];
} *P;     //--- Pointer to structure

int     nPart;
int    *tetras;
//int     ng,nx,ny,nz;
float **cube=NULL;


//===============================================================
//   MAIN
//===============================================================
int main(int argc, char **argv)
{

  int ng,nx,ny,nz;

  if (argc != 7)
    {
      printf("<-----------------------------------------------------------> \n");
      printf("Usage: \n");
      printf("       ./ltfe_soft gadget_file n_grid_particles nx ny nz file_out\n");
      printf("       \n");
      printf("<-----------------------------------------------------------> \n");
      exit(0);
    }
  ng = atoi(argv[2]);
  nx = atoi(argv[3]);
  ny = atoi(argv[4]);
  nz = atoi(argv[5]);
  printf(">>> ng = %d\n", ng);
  printf(">>> nx = %d\n", nx);
  printf(">>> ny = %d\n", ny);
  printf(">>> nz = %d\n", nz);


  //--- Read Gadget file
  load_snapshot(argv[1], 1);

  //--- Allocate memory for particles
  allocate_memory_cube(nx,ny,nz);

  //--- Compute densities
  printf("densities 0,0,0\n");
  get_densities(ng, nx,ny,nz, cube, 0,0,0);
  printf("densities 0,0,1\n");
  get_densities(ng, nx,ny,nz, cube, 0,0,1);
  printf("densities 0,1,0\n");
  get_densities(ng, nx,ny,nz, cube, 0,1,0);
  printf("densities 0,1,1\n");
  get_densities(ng, nx,ny,nz, cube, 0,1,1);
  printf("densities 1,0,0\n");
  get_densities(ng, nx,ny,nz, cube, 1,0,0);
  printf("densities 1,0,1\n");
  get_densities(ng, nx,ny,nz, cube, 1,0,1);
  printf("densities 1,1,0\n");
  get_densities(ng, nx,ny,nz, cube, 1,1,0);
  printf("densities 1,1,1\n");
  get_densities(ng, nx,ny,nz, cube, 1,1,1);

  //--- 
  save_cube(argv[6], cube, nx, nx,ny,nz);

  //--- Free memory
  free(P);
  for (int i=0; i<nz; ++i) free(&(cube[i][0]));
  free(cube);

  return(0);

}

//----------------------------------------------------------
// 
//----------------------------------------------------------
void get_densities(int ng, int nx, int ny, int nz, float **cube, int mirror_x, int mirror_y, int mirror_z){
  
  printf("N = %d\n", nPart);

  long ng2 = ng*ng;
  long ng3 = ng*ng*ng;

  //=========================================
  //   Parallel block
  //=========================================
  #pragma omp parallel
  {

    long i0,i1,j0,j1,k0,k1;

    //======================================
    // MAIN LOOP
    //======================================
    #pragma omp for

    for (long k=0; k<ng; k++) {
      //--- Allocate tetrahedron array
      float **tetra=NULL;
      tetra = (float **) malloc(sizeof(float **)*4);
      for (int q=0;q<4;q++) tetra[q] = (float *) malloc(sizeof(float *)*4);       
      
      long   ind[6][4];
      double xt[4],yt[4],zt[4];

      for (long j=0; j<ng; j++) {
	for (long i=0; i<ng; i++) {
	  
	  //--- Default case
	  i0 = i;
	  i1 = i+1;
	  j0 = j;
	  j1 = j+1;
	  k0 = k;
	  k1 = k+1;

	  //--- Mirror indices
	  if (mirror_x == 1) {
	    i0 = i+1;
	    i1 = i;
	  }
	  if (mirror_y == 1) {
	    j0 = j+1;
	    j1 = j;
	  }
	  if (mirror_z == 1) {
	    k0 = k+1;
	    k1 = k;
	  }

	  //---Check boundaries
	  if (i0 == -1) i0 = ng-1;
	  if (i1 == -1) i1 = ng-1;
	  if (j0 == -1) j0 = ng-1;
	  if (j1 == -1) j1 = ng-1;
	  if (k0 == -1) k0 = ng-1;
	  if (k1 == -1) k1 = ng-1;

	  ind[0][0] = i0 + ng*j0 + ng2*k0;
	  ind[0][1] = i0 + ng*j1 + ng2*k0;
	  ind[0][2] = i0 + ng*j0 + ng2*k1;
	  ind[0][3] = i1 + ng*j0 + ng2*k1;

	  ind[1][0] = i0 + ng*j0 + ng2*k0;
	  ind[1][1] = i0 + ng*j1 + ng2*k0;
	  ind[1][2] = i1 + ng*j0 + ng2*k1;
	  ind[1][3] = i1 + ng*j0 + ng2*k0;

	  ind[2][0] = i0 + ng*j0 + ng2*k1;
	  ind[2][1] = i0 + ng*j1 + ng2*k1;
	  ind[2][2] = i1 + ng*j0 + ng2*k1;
	  ind[2][3] = i0 + ng*j1 + ng2*k0;

	  ind[3][0] = i0 + ng*j1 + ng2*k0;
	  ind[3][1] = i1 + ng*j0 + ng2*k1;
	  ind[3][2] = i1 + ng*j1 + ng2*k1;
	  ind[3][3] = i0 + ng*j1 + ng2*k1;

	  ind[4][0] = i0 + ng*j1 + ng2*k0;
	  ind[4][1] = i1 + ng*j0 + ng2*k1;
	  ind[4][2] = i1 + ng*j1 + ng2*k1;
	  ind[4][3] = i1 + ng*j1 + ng2*k0;

	  ind[5][0] = i0 + ng*j1 + ng2*k0;
	  ind[5][1] = i1 + ng*j0 + ng2*k1;
	  ind[5][2] = i1 + ng*j1 + ng2*k0;
	  ind[5][3] = i1 + ng*j0 + ng2*k0;
	  
	  //--- Fix overflow
	  if (ind[0][0] >= nPart) ind[0][0] -= nPart;
	  if (ind[0][1] >= nPart) ind[0][1] -= nPart;
	  if (ind[0][2] >= nPart) ind[0][2] -= nPart;
	  if (ind[0][3] >= nPart) ind[0][3] -= nPart;
	  //
	  if (ind[1][0] >= nPart) ind[1][0] -= nPart;
	  if (ind[1][1] >= nPart) ind[1][1] -= nPart;
	  if (ind[1][2] >= nPart) ind[1][2] -= nPart;
	  if (ind[1][3] >= nPart) ind[1][3] -= nPart;
	  //
	  if (ind[2][0] >= nPart) ind[2][0] -= nPart;
	  if (ind[2][1] >= nPart) ind[2][1] -= nPart;
	  if (ind[2][2] >= nPart) ind[2][2] -= nPart;
	  if (ind[2][3] >= nPart) ind[2][3] -= nPart;
	  //
	  if (ind[3][0] >= nPart) ind[3][0] -= nPart;
	  if (ind[3][1] >= nPart) ind[3][1] -= nPart;
	  if (ind[3][2] >= nPart) ind[3][2] -= nPart;
	  if (ind[3][3] >= nPart) ind[3][3] -= nPart;
	  //
	  if (ind[4][0] >= nPart) ind[4][0] -= nPart;
	  if (ind[4][1] >= nPart) ind[4][1] -= nPart;
	  if (ind[4][2] >= nPart) ind[4][2] -= nPart;
	  if (ind[4][3] >= nPart) ind[4][3] -= nPart;
	  //
	  if (ind[5][0] >= nPart) ind[5][0] -= nPart;
	  if (ind[5][1] >= nPart) ind[5][1] -= nPart;
	  if (ind[5][2] >= nPart) ind[5][2] -= nPart;
	  if (ind[5][3] >= nPart) ind[5][3] -= nPart;

	  float cnt_x = ((float) nx) / ((float) header1.BoxSize);
	  float cnt_y = ((float) ny) / ((float) header1.BoxSize);
	  float cnt_z = ((float) nz) / ((float) header1.BoxSize);
	  //--- Make sure to unroll at compile time
	  for (int w=0; w<6; w++){
	    //--- Make sure to unroll at compile time
	    for (int ww=0; ww<4;ww++) {
	      xt[ww] = (double) (P[ind[w][ww]].Pos[0]);
	      yt[ww] = (double) (P[ind[w][ww]].Pos[1]);
	      zt[ww] = (double) (P[ind[w][ww]].Pos[2]);
	      tetra[ww][0] =     P[ind[w][ww]].Pos[0] * cnt_x;
	      tetra[ww][1] =     P[ind[w][ww]].Pos[1] * cnt_y;
	      tetra[ww][2] =     P[ind[w][ww]].Pos[2] * cnt_z;
	    } //--- end for ww

	    //--- Fix periodic tetrahedra, OJO: TODO: must replicate!
	    fix_periodic_tetra(&xt[0],&yt[0],&zt[0], header1.BoxSize);

	    //--- Compute tetrahedron volume
	    double tetra_vol = fabs(tetra_volu_ori(xt,yt,zt));
	    
	    //--- Process this tetrahedron...
	    float den_i = (float) (1.0/tetra_vol);
	    process_tetrahedron(&(tetra[0]), den_i, &(cube[0]), nx,ny,nz);
	  } //--- end for w
	  
	}//--- end for k
      }//--- end for j
    //--- Each thread frees its tetra
    for (int q=0; q<4; q++) free(&(tetra[q][0]));
    free(tetra);
    }//--- end for i, main loop
    
    
    //--- Make sure all threads are ready before continuing
    #pragma omp barrier
    
  } // #pragma omp parallel for private(i)
  
  printf(">>> ready with densities!\n");
  
} //--- end get_densities()  




//----------------------------------------------------------------
//
//----------------------------------------------------------------
void process_tetrahedron(float **tetra, float den, float **cube, int nx, int ny, int nz){


  int   ind_x[4], ind_y[4], ind_z[4];
  int   per_x=0, per_y=0, per_z=0, per_xyz;
  //--- Asssume particles cover all box and go from 0-nxyz
  float nx2 = nx/2.0;
  float ny2 = ny/2.0;
  float nz2 = nz/2.0;

  //--- Fix periodic boundaries
  /*
  if (tetra[0][0] >= nx) tetra[0][0] -= nx;
  if (tetra[1][0] >= nx) tetra[1][0] -= nx;
  if (tetra[2][0] >= nx) tetra[2][0] -= nx;
  if (tetra[3][0] >= nx) tetra[3][0] -= nx;
  if (tetra[0][1] >= nx) tetra[0][1] -= ny;
  if (tetra[1][1] >= nx) tetra[1][1] -= ny;
  if (tetra[2][1] >= nx) tetra[2][1] -= ny;
  if (tetra[3][1] >= nx) tetra[3][1] -= ny;
  if (tetra[0][2] >= nx) tetra[0][2] -= nz;
  if (tetra[1][2] >= nx) tetra[1][2] -= nz;
  if (tetra[2][2] >= nx) tetra[2][2] -= nz;
  if (tetra[3][2] >= nx) tetra[3][2] -= nz;
  */

  //--- Sort indexes...
  sortFour(tetra[0][0], tetra[1][0], tetra[2][0], tetra[3][0], ind_x);
  sortFour(tetra[0][1], tetra[1][1], tetra[2][1], tetra[3][1], ind_y);
  sortFour(tetra[0][2], tetra[1][2], tetra[2][2], tetra[3][2], ind_z);  
  //--- Find periodic tetrahedron
  if ( (tetra[ind_x[3]][0] - tetra[ind_x[0]][0]) > nx2 ) per_x=1; // 2^0
  if ( (tetra[ind_y[3]][1] - tetra[ind_y[0]][1]) > ny2 ) per_y=2; // 2^1
  if ( (tetra[ind_z[3]][2] - tetra[ind_z[0]][2]) > nz2 ) per_z=4; // 2^2
  //--- Encode periodicity
  per_xyz = per_x + per_y + per_z;
  
  //--- Choose 
  switch( per_xyz ) {
  case 0:
    render_tetrahedron_in_cube(&(tetra[0]), den, cube, nx,ny,nz);
    break;
  case 1:
    if (tetra[0][0] >= nx2 ) tetra[0][0] -= nx;
    if (tetra[1][0] >= nx2 ) tetra[1][0] -= nx;
    if (tetra[2][0] >= nx2 ) tetra[2][0] -= nx;
    if (tetra[3][0] >= nx2 ) tetra[3][0] -= nx;
    render_tetrahedron_in_cube(&(tetra[0]), den, cube, nx,ny,nz);
    tetra[0][0] += nx;
    tetra[1][0] += nx;
    tetra[2][0] += nx;
    tetra[3][0] += nx;
    render_tetrahedron_in_cube(&(tetra[0]), den, cube, nx,ny,nz);
    break;
  case 2:
    if (tetra[0][1] >= ny2 ) tetra[0][1] -= ny;
    if (tetra[1][1] >= ny2 ) tetra[1][1] -= ny;
    if (tetra[2][1] >= ny2 ) tetra[2][1] -= ny;
    if (tetra[3][1] >= ny2 ) tetra[3][1] -= ny;
    render_tetrahedron_in_cube(&(tetra[0]), den, cube, nx,ny,nz);
    tetra[0][1] += ny;
    tetra[1][1] += ny;
    tetra[2][1] += ny;
    tetra[3][1] += ny;
    render_tetrahedron_in_cube(&(tetra[0]), den, cube, nx,ny,nz);
    break;
  case 4:
    if (tetra[0][2] > nz2 ) tetra[0][2] -= nz;
    if (tetra[1][2] > nz2 ) tetra[1][2] -= nz;
    if (tetra[2][2] > nz2 ) tetra[2][2] -= nz;
    if (tetra[3][2] > nz2 ) tetra[3][2] -= nz;
    render_tetrahedron_in_cube(&(tetra[0]), den, cube, nx,ny,nz);
    tetra[0][2] += nz;
    tetra[1][2] += nz;
    tetra[2][2] += nz;
    tetra[3][2] += nz;
    render_tetrahedron_in_cube(&(tetra[0]), den, cube, nx,ny,nz);
    break;
  case 3:
    if (tetra[0][0] > nx2 ) tetra[0][0] -= nx;
    if (tetra[1][0] > nx2 ) tetra[1][0] -= nx;
    if (tetra[2][0] > nx2 ) tetra[2][0] -= nx;
    if (tetra[3][0] > nx2 ) tetra[3][0] -= nx;
    if (tetra[0][1] > ny2 ) tetra[0][1] -= ny;
    if (tetra[1][1] > ny2 ) tetra[1][1] -= ny;
    if (tetra[2][1] > ny2 ) tetra[2][1] -= ny;
    if (tetra[3][1] > ny2 ) tetra[3][1] -= ny;
    render_tetrahedron_in_cube(&(tetra[0]), den, cube, nx,ny,nz);
    tetra[0][0] += nx;
    tetra[1][0] += nx;
    tetra[2][0] += nx;
    tetra[3][0] += nx;
    render_tetrahedron_in_cube(&(tetra[0]), den, cube, nx,ny,nz);
    tetra[0][0] -= nx;
    tetra[1][0] -= nx;
    tetra[2][0] -= nx;
    tetra[3][0] -= nx;
    tetra[0][1] += ny;
    tetra[1][1] += ny;
    tetra[2][1] += ny;
    tetra[3][1] += ny;
    render_tetrahedron_in_cube(&(tetra[0]), den, cube, nx,ny,nz);
    tetra[0][0] += nx;
    tetra[1][0] += nx;
    tetra[2][0] += nx;
    tetra[3][0] += nx;
    render_tetrahedron_in_cube(&(tetra[0]), den, cube, nx,ny,nz);
    break;
  case 6:
    if (tetra[0][1] > ny2 ) tetra[0][1] -= ny;
    if (tetra[1][1] > ny2 ) tetra[1][1] -= ny;
    if (tetra[2][1] > ny2 ) tetra[2][1] -= ny;
    if (tetra[3][1] > ny2 ) tetra[3][1] -= ny;
    if (tetra[0][2] > nz2 ) tetra[0][2] -= nz;
    if (tetra[1][2] > nz2 ) tetra[1][2] -= nz;
    if (tetra[2][2] > nz2 ) tetra[2][2] -= nz;
    if (tetra[3][2] > nz2 ) tetra[3][2] -= nz;
    render_tetrahedron_in_cube(&(tetra[0]), den, cube, nx,ny,nz);
    tetra[0][1] += ny;
    tetra[1][1] += ny;
    tetra[2][1] += ny;
    tetra[3][1] += ny;
    render_tetrahedron_in_cube(&(tetra[0]), den, cube, nx,ny,nz);
    tetra[0][1] -= ny;
    tetra[1][1] -= ny;
    tetra[2][1] -= ny;
    tetra[3][1] -= ny;
    tetra[0][2] += nz;
    tetra[1][2] += nz;
    tetra[2][2] += nz;
    tetra[3][2] += nz;
    render_tetrahedron_in_cube(&(tetra[0]), den, cube, nx,ny,nz);
    tetra[0][1] += ny;
    tetra[1][1] += ny;
    tetra[2][1] += ny;
    tetra[3][1] += ny;
    render_tetrahedron_in_cube(&(tetra[0]), den, cube, nx,ny,nz);
    break;
  case 5:
    if (tetra[0][0] > nx2 ) tetra[0][0] -= nx;
    if (tetra[1][0] > nx2 ) tetra[1][0] -= nx;
    if (tetra[2][0] > nx2 ) tetra[2][0] -= nx;
    if (tetra[3][0] > nx2 ) tetra[3][0] -= nx;
    if (tetra[0][2] > nz2 ) tetra[0][2] -= nz;
    if (tetra[1][2] > nz2 ) tetra[1][2] -= nz;
    if (tetra[2][2] > nz2 ) tetra[2][2] -= nz;
    if (tetra[3][2] > nz2 ) tetra[3][2] -= nz;
    render_tetrahedron_in_cube(&(tetra[0]), den, cube, nx,ny,nz);
    tetra[0][0] += nx;
    tetra[1][0] += nx;
    tetra[2][0] += nx;
    tetra[3][0] += nx;
    render_tetrahedron_in_cube(&(tetra[0]), den, cube, nx,ny,nz);
    tetra[0][0] -= nx;
    tetra[1][0] -= nx;
    tetra[2][0] -= nx;
    tetra[3][0] -= nx;
    tetra[0][2] += nz;
    tetra[1][2] += nz;
    tetra[2][2] += nz;
    tetra[3][2] += nz;
    render_tetrahedron_in_cube(&(tetra[0]), den, cube, nx,ny,nz);
    tetra[0][0] += nx;
    tetra[1][0] += nx;
    tetra[2][0] += nx;
    tetra[3][0] += nx;
    render_tetrahedron_in_cube(&(tetra[0]), den, cube, nx,ny,nz);
    break;
  case 7:
    if (tetra[0][0] > nx2 ) tetra[0][0] -= nx;
    if (tetra[1][0] > nx2 ) tetra[1][0] -= nx;
    if (tetra[2][0] > nx2 ) tetra[2][0] -= nx;
    if (tetra[3][0] > nx2 ) tetra[3][0] -= nx;
    if (tetra[0][1] > ny2 ) tetra[0][1] -= ny;
    if (tetra[1][1] > ny2 ) tetra[1][1] -= ny;
    if (tetra[2][1] > ny2 ) tetra[2][1] -= ny;
    if (tetra[3][1] > ny2 ) tetra[3][1] -= ny;
    if (tetra[0][2] > nz2 ) tetra[0][2] -= nz;
    if (tetra[1][2] > nz2 ) tetra[1][2] -= nz;
    if (tetra[2][2] > nz2 ) tetra[2][2] -= nz;
    if (tetra[3][2] > nz2 ) tetra[3][2] -= nz;
    render_tetrahedron_in_cube(&(tetra[0]), den, cube, nx,ny,nz); // 0 0 0
    tetra[0][0] += nx;
    tetra[1][0] += nx;
    tetra[2][0] += nx;
    tetra[3][0] += nx;
    render_tetrahedron_in_cube(&(tetra[0]), den, cube, nx,ny,nz); // 0 0 1
    tetra[0][0] -= nx;
    tetra[1][0] -= nx;
    tetra[2][0] -= nx;
    tetra[3][0] -= nx;
    tetra[0][1] += ny;
    tetra[1][1] += ny;
    tetra[2][1] += ny;
    tetra[3][1] += ny;
    render_tetrahedron_in_cube(&(tetra[0]), den, cube, nx,ny,nz); // 0 1 0
    tetra[0][0] += nx;
    tetra[1][0] += nx;
    tetra[2][0] += nx;
    tetra[3][0] += nx;
    render_tetrahedron_in_cube(&(tetra[0]), den, cube, nx,ny,nz); // 0 1 1 
    tetra[0][0] -= nx;
    tetra[1][0] -= nx;
    tetra[2][0] -= nx;
    tetra[3][0] -= nx;
    tetra[0][1] -= ny;
    tetra[1][1] -= ny;
    tetra[2][1] -= ny;
    tetra[3][1] -= ny;
    tetra[0][2] += nz;
    tetra[1][2] += nz;
    tetra[2][2] += nz;
    tetra[3][2] += nz;
    render_tetrahedron_in_cube(&(tetra[0]), den, cube, nx,ny,nz); // 1 0 0 
    tetra[0][0] += nx;
    tetra[1][0] += nx;
    tetra[2][0] += nx;
    tetra[3][0] += nx;
    render_tetrahedron_in_cube(&(tetra[0]), den, cube, nx,ny,nz); // 1 0 1
    tetra[0][0] -= nx;
    tetra[1][0] -= nx;
    tetra[2][0] -= nx;
    tetra[3][0] -= nx;
    tetra[0][1] += ny;
    tetra[1][1] += ny;
    tetra[2][1] += ny;
    tetra[3][1] += ny;
    render_tetrahedron_in_cube(&(tetra[0]), den, cube, nx,ny,nz); // 1 1 0
    tetra[0][0] += nx;
    tetra[1][0] += nx;
    tetra[2][0] += nx;
    tetra[3][0] += nx;
    render_tetrahedron_in_cube(&(tetra[0]), den, cube, nx,ny,nz); // 1 1 1
   } //--- switch

}

//----------------------------------------------------------------
//
//----------------------------------------------------------------
void render_tetrahedron_in_cube(float **tetra, float den, float **cube, int nx, int ny, int nz){

  int   i,w;
     int   z0,z1;  //--- loop indexes
     float iso_ver[6][3];
     float area;
     float isoval;
     int   n_valid_triangles;

     //--- Sort tetras with their z-position
     float tet_val[4] = { tetra[0][2], tetra[1][2], tetra[2][2], tetra[3][2] };
     int   tet_ind[4];

     //--- Discard tetrahedron outside of rendering box...
     //--- X
     //sortFour(tetra[0][0], tetra[1][0], tetra[2][0], tetra[3][0], tet_ind);
     //if ( ((tetra[tet_ind[3]][0] < 0) || (tetra[tet_ind[0]][0] > nx)) ) return;
     //--- Y
     //sortFour(tetra[0][1], tetra[1][1], tetra[2][1], tetra[3][1], tet_ind);
     //if ( ((tetra[tet_ind[3]][1] < 0) || (tetra[tet_ind[0]][1] > ny)) ) return;
     //--- Z
     //sortFour(tetra[0][2], tetra[1][2], tetra[2][2], tetra[3][2], tet_ind); 
     //if ( ((tetra[tet_ind[3]][2] < 0) || (tetra[tet_ind[0]][2] > nz)) ) return;
     //--- Reject tetras outside of rendering volume
     //if ( !((tetra[tet_ind[0]][0] < 0) || (tetra[tet_ind[3]][0] > nx)) &&
	  // !((tetra[tet_ind[0]][1] < 0) || (tetra[tet_ind[3]][1] > ny)) &&
	  //  !((tetra[tet_ind[0]][2] < 0) || (tetra[tet_ind[3]][2] > nz)) ) return;

     //--- Tetrahedron inside pixel give its total mass to pixel (OJO: should weight..)
     //if ( (tetra[tet_ind[0]][0] == tetra[tet_ind[3]][0]) &&
     // 	  (tetra[tet_ind[0]][1] == tetra[tet_ind[3]][1]) &&
     //	  (tetra[tet_ind[0]][2] == tetra[tet_ind[3]][2]) ) {
     // cube[lround(tetra[tet_ind[0]][2])][lround(tetra[tet_ind[0]][1])*nx + lround(tetra[tet_ind[0]][0])] = tetra[tet_ind[0]][3];
     //  return;
     //}
     
     //--- Integer steps, this only works if positions are in range [0,nx]...
     //    Reuse z-sorted vertices
     sortFour(tetra[0][2], tetra[1][2], tetra[2][2], tetra[3][2], tet_ind);
     z0 = lround(tet_val[tet_ind[0]]);
     z1 = lround(tet_val[tet_ind[3]]);

     //--- Flat tetrahedra, try finer sampling?
     //if (z0 == z1) return;
     //--- Clip tetrahedron of out of box
     if (z0 <   0) z0 = 0;
     if (z1 >= nz) z1 = nz-1;
     // OJO: TODO: Discard tetrahedron outside the rendering volume

     //----------------------------------------------
     //  SLICE THE TETRAHEDRON
     //----------------------------------------------
     for (int i=z0; i<=z1;i++){  // = faster than <= ?

       //--- Define isoplane
       isoval = float(i) + 0.5;
          
       //--- CUT TETRA ALONG ISOVALUE
       n_valid_triangles = isocut_tetrahedron(isoval, tetra, iso_ver);

       //-- If no intersection then return
       if (n_valid_triangles == 0) return;

       //===================================================================================
       //  TRIANGLE INTERPOLATION
       //===================================================================================

       int FLAG_INTERPOL = 0;

       if (FLAG_INTERPOL != 0) {
	 //------------------
	 //  TRIANGLE
	 //------------------
	 if (n_valid_triangles == 3) {
	   
	   area = ((iso_ver[1][0]-iso_ver[0][0])*(iso_ver[2][1]-iso_ver[0][1]) - (iso_ver[2][0]-iso_ver[0][0])*(iso_ver[1][1]-iso_ver[0][1]));
	   if (area < 0)
	     raster_triangle_interpol( &(iso_ver[0][0]), &(iso_ver[1][0]), &(iso_ver[2][0]), den, &(cube[i][0]), nx,ny, 8);
	   else
	     raster_triangle_interpol( &(iso_ver[0][0]), &(iso_ver[2][0]), &(iso_ver[1][0]), den, &(cube[i][0]), nx,ny, 8);	 
	 } 
	 //------------------
	 //  QUAD
	 //------------------
	 else{ 
	   area = ((iso_ver[1][0]-iso_ver[0][0])*(iso_ver[2][1]-iso_ver[0][1]) - (iso_ver[2][0]-iso_ver[0][0])*(iso_ver[1][1]-iso_ver[0][1]));
	   if (area < 0)
	     raster_triangle_interpol( &(iso_ver[0][0]), &(iso_ver[1][0]), &(iso_ver[2][0]), den, &(cube[i][0]), nx,ny, 8);
	   else
	     raster_triangle_interpol( &(iso_ver[0][0]), &(iso_ver[2][0]), &(iso_ver[1][0]), den, &(cube[i][0]), nx,ny, 8);
	   
	   //=== ADD MISSING VERTICES FOR SECOND TRIANGLE [triangle repeated, fix !!!]
	   //--- VERTEX 2
	   //--- Add the two elements before the last
	   iso_ver[4][0] = iso_ver[2][0];
	   iso_ver[4][1] = iso_ver[2][1];
	   iso_ver[4][2] = iso_ver[2][2];
	   //--- VERTEX 2
	   iso_ver[5][0] = iso_ver[1][0];
	   iso_ver[5][1] = iso_ver[1][1];
	   iso_ver[5][2] = iso_ver[1][2];
	   
	   area = ((iso_ver[4][0]-iso_ver[3][0])*(iso_ver[5][1]-iso_ver[3][1]) - (iso_ver[5][0]-iso_ver[3][0])*(iso_ver[4][1]-iso_ver[3][1]));
	   if (area < 0)
	     raster_triangle_interpol( &(iso_ver[3][0]), &(iso_ver[4][0]), &(iso_ver[5][0]), den, &(cube[i][0]), nx,ny, 8);
	   else
	     raster_triangle_interpol( &(iso_ver[3][0]), &(iso_ver[5][0]), &(iso_ver[4][0]), den, &(cube[i][0]), nx,ny, 8);
	   
	 } // else

       } //--- FLAG_INTERPOL
       else {
	 float x_min,x_max,y_min,y_max,l_thres;
	 //------------------
	 //  TRIANGLE
	 //------------------
	 if (n_valid_triangles == 3) {
	      
	   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	   //--- OJO OJO OJO!!! Problem is here!!! finally found!
	   //    Some triangles are given very large sizes.
	   x_min = MIN(iso_ver[0][0], iso_ver[1][0], iso_ver[2][0]);
	   x_max = MAX(iso_ver[0][0], iso_ver[1][0], iso_ver[2][0]);
	   y_min = MIN(iso_ver[0][1], iso_ver[1][1], iso_ver[2][1]);
	   y_max = MAX(iso_ver[0][1], iso_ver[1][1], iso_ver[2][1]);
	   l_thres = 8.0; //--- This trick only works for large particle grids >> 8
	   if ((x_max-x_min) > nx/l_thres ) return;
	   if ((y_max-y_min) > ny/l_thres ) return;

	   area = ((iso_ver[1][0]-iso_ver[0][0])*(iso_ver[2][1]-iso_ver[0][1]) - (iso_ver[2][0]-iso_ver[0][0])*(iso_ver[1][1]-iso_ver[0][1]));
	   if (area < 0)
	     raster_triangle( &(iso_ver[0][0]), &(iso_ver[1][0]), &(iso_ver[2][0]), den, &(cube[i][0]), nx,ny, 8);
	   else
	     raster_triangle( &(iso_ver[0][0]), &(iso_ver[2][0]), &(iso_ver[1][0]), den, &(cube[i][0]), nx,ny, 8); 
	 } 
	 //------------------
	 //  QUAD
	 //------------------
	 else{    

	   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	   //--- OJO OJO OJO!!! Problem is here!!! finally found!
	   //    Some triangles are given very large sizes.
	   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	   x_min = MIN(iso_ver[0][0], iso_ver[1][0], iso_ver[2][0]);
	   x_max = MAX(iso_ver[0][0], iso_ver[1][0], iso_ver[2][0]);
	   y_min = MIN(iso_ver[0][1], iso_ver[1][1], iso_ver[2][1]);
	   y_max = MAX(iso_ver[0][1], iso_ver[1][1], iso_ver[2][1]);
	   l_thres = 8.0; //--- This trick only works for large particle grids >> 8
	   if ((x_max-x_min) > nx/l_thres ) return;
	   if ((y_max-y_min) > ny/l_thres ) return;

	   area = ((iso_ver[1][0]-iso_ver[0][0])*(iso_ver[2][1]-iso_ver[0][1]) - (iso_ver[2][0]-iso_ver[0][0])*(iso_ver[1][1]-iso_ver[0][1]));
	   if (area < 0)
	     raster_triangle( &(iso_ver[0][0]), &(iso_ver[1][0]), &(iso_ver[2][0]), den, &(cube[i][0]), nx,ny, 8);
	   else
	     raster_triangle( &(iso_ver[0][0]), &(iso_ver[2][0]), &(iso_ver[1][0]), den, &(cube[i][0]), nx,ny, 8);

	   //=== ADD MISSING VERTICES FOR SECOND TRIANGLE [triangle repeated, fix !!!]
	   //--- VERTEX 2
	   //--- Add the two elements before the last
	   iso_ver[4][0] = iso_ver[2][0];
	   iso_ver[4][1] = iso_ver[2][1];
	   iso_ver[4][2] = iso_ver[2][2];
	   //--- VERTEX 2
	   iso_ver[5][0] = iso_ver[1][0];
	   iso_ver[5][1] = iso_ver[1][1];
	   iso_ver[5][2] = iso_ver[1][2];
	      
	   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	   //--- OJO OJO OJO!!! Problem is here!!! finally found!
	   //    Some triangles are given very large sizes.
	   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	   x_min = MIN(iso_ver[3][0], iso_ver[4][0], iso_ver[5][0]);
	   x_max = MAX(iso_ver[3][0], iso_ver[4][0], iso_ver[5][0]);
	   y_min = MIN(iso_ver[3][1], iso_ver[4][1], iso_ver[5][1]);
	   y_max = MAX(iso_ver[3][1], iso_ver[4][1], iso_ver[5][1]);
	   l_thres = 8.0; //--- This trick only works for large particle grids >> 8
	   if ((x_max-x_min) > nx/l_thres ) return;
	   if ((y_max-y_min) > ny/l_thres ) return;

	   area = ((iso_ver[4][0]-iso_ver[3][0])*(iso_ver[5][1]-iso_ver[3][1]) - (iso_ver[5][0]-iso_ver[3][0])*(iso_ver[4][1]-iso_ver[3][1]));
	   if (area < 0)
	     raster_triangle( &(iso_ver[3][0]), &(iso_ver[4][0]), &(iso_ver[5][0]), den, &(cube[i][0]), nx,ny, 8);
	   else
	     raster_triangle( &(iso_ver[3][0]), &(iso_ver[5][0]), &(iso_ver[4][0]), den, &(cube[i][0]), nx,ny, 8);
	      
	 } // else

       }

     } //--- for i
     

} // render_tetrahedron_in_cube()


//----------------------------------------------------------------
//
//----------------------------------------------------------------
int isocut_tetrahedron(float isoval, float **tetra, float iso_ver[6][3]){
     int i,j;
     int z0,z1;  //--- loop indexes

     //--- Get density at the vertices tetrahedra and sort them  tetra[vert][coord]
     float tet_val[4] = { tetra[0][2], tetra[1][2], tetra[2][2], tetra[3][2] }; 
     int   tet_ind[4];
     sortFour(tetra[0][2], tetra[1][2], tetra[2][2], tetra[3][2], tet_ind); //--- Sort z direction

     //--- Get rid of tetrahedra in plane. Remember vertices are sorted
     //if(tet_val[tet_ind[0]]==tet_val[tet_ind[3]]) return 0;
     
     //--- Intersect tetrahedron with isovalue
     int   cnt = 0;
     if (isocut_line(isoval, tet_val[tet_ind[0]], tet_val[tet_ind[1]], &(tetra[tet_ind[0]][0]), &(tetra[tet_ind[1]][0]), &(iso_ver[cnt][0]))) cnt++;
     if (isocut_line(isoval, tet_val[tet_ind[0]], tet_val[tet_ind[2]], &(tetra[tet_ind[0]][0]), &(tetra[tet_ind[2]][0]), &(iso_ver[cnt][0]))) cnt++;
     if (isocut_line(isoval, tet_val[tet_ind[0]], tet_val[tet_ind[3]], &(tetra[tet_ind[0]][0]), &(tetra[tet_ind[3]][0]), &(iso_ver[cnt][0]))) cnt++;
     if (isocut_line(isoval, tet_val[tet_ind[1]], tet_val[tet_ind[2]], &(tetra[tet_ind[1]][0]), &(tetra[tet_ind[2]][0]), &(iso_ver[cnt][0]))) cnt++;
     if (isocut_line(isoval, tet_val[tet_ind[1]], tet_val[tet_ind[3]], &(tetra[tet_ind[1]][0]), &(tetra[tet_ind[3]][0]), &(iso_ver[cnt][0]))) cnt++;
     if (isocut_line(isoval, tet_val[tet_ind[2]], tet_val[tet_ind[3]], &(tetra[tet_ind[2]][0]), &(tetra[tet_ind[3]][0]), &(iso_ver[cnt][0]))) cnt++;
 
     return cnt;

} // isocut_tetrahedron()


//----------------------------------------------------------------
//
//----------------------------------------------------------------
bool isocut_line(const float isoval, const float val1, const float val2, const float *vert1, const float *vert2, float *vert_inter){

  //--- Check if isoval lies between the two values:
  if ( (val1 < isoval) && (isoval < val2) ) {
    //--- Two-point parametric curve of the line along the two vertices
    float t_par = (isoval - vert1[2]) / (vert2[2] - vert1[2]);
    //--- Interpolate values
    vert_inter[0] = (1-t_par)*vert1[0] + t_par*vert2[0]; //--- X
    vert_inter[1] = (1-t_par)*vert1[1] + t_par*vert2[1]; //--- Y
    vert_inter[2] = (1-t_par)*vert1[3] + t_par*vert2[3]; //--- Density
   return true;
  } else {
    return false;
  }

}// isocut_line()


//----------------------------------------------------------------
//
//  v1  ->  [x,y,den]
//  den ->  density of full triangle
//
//----------------------------------------------------------------
void raster_triangle(const float *v1, 
		     const float *v2, 
		     const float *v3, 
		     const float den, float *buffer, const int xg, const int yg, const int q){
  
  int ind;
  int stride = xg;
  int y,x,iy,ix;
  int tr;


  // Block size, standard 8x8 (must be power of two)
  //const int q = 16;

  float x1 = v1[0];
  float x2 = v2[0];
  float x3 = v3[0];
  float y1 = v1[1];
  float y2 = v2[1];
  float y3 = v3[1];

  // 28.4 fixed-point coordinates
  const int   X1 = (int) round(16.0f * x1);
  const int   X2 = (int) round(16.0f * x2);
  const int   X3 = (int) round(16.0f * x3);
  
  const int   Y1 = (int) round(16.0f * y1);
  const int   Y2 = (int) round(16.0f * y2);
  const int   Y3 = (int) round(16.0f * y3);

  //--- Densities at each vertex
  const float w1 = v1[2];
  const float w2 = v2[2];
  const float w3 = v3[2];

  // Deltas
  const int DX12 = X1 - X2;
  const int DX23 = X2 - X3;
  const int DX31 = X3 - X1;
  
  const int DY12 = Y1 - Y2;
  const int DY23 = Y2 - Y3;
  const int DY31 = Y3 - Y1;
  
  // Fixed-point deltas
  const int FDX12 = DX12 << 4;
  const int FDX23 = DX23 << 4;
  const int FDX31 = DX31 << 4;

  const int FDY12 = DY12 << 4;
  const int FDY23 = DY23 << 4;
  const int FDY31 = DY31 << 4;

  // Bounding rectangle
  int minx = (int) ceil((LOW(MIN(X1, X2, X3),    0))/16.0);
  int maxx = (int) ceil((TOP(MAX(X1, X2, X3),xg*16))/16.0);
  int miny = (int) ceil((LOW(MIN(Y1, Y2, Y3),    0))/16.0);
  int maxy = (int) ceil((TOP(MAX(Y1, Y2, Y3),yg*16))/16.0);
  
  //--- Original, no clipping...
  //int minx = ( MIN(X1, X2, X3) + 0xF) >> 4;
  //int maxx = ( MAX(X1, X2, X3) + 0xF) >> 4;
  //int miny = ( MIN(Y1, Y2, Y3) + 0xF) >> 4;
  //int maxy = ( MAX(Y1, Y2, Y3) + 0xF) >> 4;    

  // Start in corner of 8x8 block
  minx &= ~(q - 1);
  miny &= ~(q - 1);
  
  // Half-edge constants
  int C1 = DY12 * X1 - DX12 * Y1;
  int C2 = DY23 * X2 - DX23 * Y2;
  int C3 = DY31 * X3 - DX31 * Y3;
  
  // Correct for fill convention
  if(DY12 < 0 || (DY12 == 0 && DX12 > 0)) C1++;
  if(DY23 < 0 || (DY23 == 0 && DX23 > 0)) C2++;
  if(DY31 < 0 || (DY31 == 0 && DX31 > 0)) C3++;
  
  
  // Loop through blocks
  for(y = miny; y < maxy; y += q){
    for(x = minx; x < maxx; x += q){
      
      // Corners of block
      int x0 = x << 4;
      int x1 = (x + q - 1) << 4;
      int y0 = y << 4;
      int y1 = (y + q - 1) << 4;
      
      // Evaluate half-space functions
      bool a00 = C1 + DX12 * y0 - DY12 * x0 > 0;
      bool a10 = C1 + DX12 * y0 - DY12 * x1 > 0;
      bool a01 = C1 + DX12 * y1 - DY12 * x0 > 0;
      bool a11 = C1 + DX12 * y1 - DY12 * x1 > 0;
      int  a   = (a00 << 0) | (a10 << 1) | (a01 << 2) | (a11 << 3);
      
      bool b00 = C2 + DX23 * y0 - DY23 * x0 > 0;
      bool b10 = C2 + DX23 * y0 - DY23 * x1 > 0;
      bool b01 = C2 + DX23 * y1 - DY23 * x0 > 0;
      bool b11 = C2 + DX23 * y1 - DY23 * x1 > 0;
      int  b   = (b00 << 0) | (b10 << 1) | (b01 << 2) | (b11 << 3);
      
      bool c00 = C3 + DX31 * y0 - DY31 * x0 > 0;
      bool c10 = C3 + DX31 * y0 - DY31 * x1 > 0;
      bool c01 = C3 + DX31 * y1 - DY31 * x0 > 0;
      bool c11 = C3 + DX31 * y1 - DY31 * x1 > 0;
      int  c   = (c00 << 0) | (c10 << 1) | (c01 << 2) | (c11 << 3);
      
      // Skip block when outside an edge
      if(a == 0x0 || b == 0x0 || c == 0x0) continue;
      
      // Accept whole block when totally covered
      if (a == 0xF && b == 0xF && c == 0xF){

	for(iy = y; iy < y+q; iy++){
	  for(ix = x; ix < x + q; ix++){		      
	    ind = iy*xg +  ix;
	    buffer[ind] += den;
	  } // end for ix
	} // end for iy		
      } // end if
      
      else { // << // Partially covered block
	
	int CY1 = C1 + DX12 * y0 - DY12 * x0;
	int CY2 = C2 + DX23 * y0 - DY23 * x0;
	int CY3 = C3 + DX31 * y0 - DY31 * x0;
	for(iy = y; iy < y + q; iy++){
	  int CX1 = CY1;
	  int CX2 = CY2;
	  int CX3 = CY3;
	  
	  for(ix = x; ix < x + q; ix++){
	    if(CX1 > 0 && CX2 > 0 && CX3 > 0){
	      ind = iy*xg +  ix;
	      buffer[ind] += den;
	    }
	    CX1 -= FDY12;
	    CX2 -= FDY23;
	    CX3 -= FDY31;
	  } // end for ix
	  CY1 += FDX12;
	  CY2 += FDX23;
	  CY3 += FDX31;
	} // for ix
      } // else
      
    } // x
  } // y 
  
} // raster_triangle()


//---------------------------------------------------------------
//  Raster triangles using bilinear interpolation
//---------------------------------------------------------------
void raster_triangle_interpol(const float *v1, 
			      const float *v2, 
			      const float *v3, 
			      const float den, float *buffer, const int xg, const int yg, const int q){
  
  int ind;
  int stride = xg;
  int y,x,iy,ix;
  int tr;

  //--- Bad fix, intercept very large triangles. This is probably some weird thing from using inline functions...
  

  // Block size, standard 8x8 (must be power of two)
  //const int q = 16;

  float x1 = v1[0];
  float x2 = v2[0];
  float x3 = v3[0];
  float y1 = v1[1];
  float y2 = v2[1];
  float y3 = v3[1];

  // 28.4 fixed-point coordinates
  const int   X1 = (int) round(16.0f * x1);
  const int   X2 = (int) round(16.0f * x2);
  const int   X3 = (int) round(16.0f * x3);
  
  const int   Y1 = (int) round(16.0f * y1);
  const int   Y2 = (int) round(16.0f * y2);
  const int   Y3 = (int) round(16.0f * y3);

  //--- Densities at each vertex
  const float w1 = v1[2];
  const float w2 = v2[2];
  const float w3 = v3[2];

  //--- Compute bilinear interpolation coeficcients 
  //    http://www1.eonfusion.com/manual/index.php/Formulae_for_interpolation
  float A,B,C;
  float DET = x1*y2 - x2*y1 + x2*y3 - x3*y2 + x3*y1 - x1*y3;
  if (DET==0){  //--- If determinant zero use constan mean density inside triangle
    A = 0;
    B = 0;
    C = (w1+w2+w3)/3.0;
  } 
  else {        //--- Otherwise bilinear interpolation
    DET = 1.0/DET;
    A   = ((y2-y3)*w1 + (y3-y1)*w2 + (y1-y2)*w3) * DET;
    B   = ((x3-x2)*w1 + (x1-x3)*w2 + (x2-x1)*w3) * DET;
    C   = ((x2*y3-x3*y2)*w1 + (x3*y1-x1*y3)*w2 + (x1*y2-x2*y1)*w3) * DET;
  }

  // Deltas
  const int DX12 = X1 - X2;
  const int DX23 = X2 - X3;
  const int DX31 = X3 - X1;
  
  const int DY12 = Y1 - Y2;
  const int DY23 = Y2 - Y3;
  const int DY31 = Y3 - Y1;
  
  // Fixed-point deltas
  const int FDX12 = DX12 << 4;
  const int FDX23 = DX23 << 4;
  const int FDX31 = DX31 << 4;

  const int FDY12 = DY12 << 4;
  const int FDY23 = DY23 << 4;
  const int FDY31 = DY31 << 4;

  // Bounding rectangle
  int minx = (int) ceil((LOW(MIN(X1, X2, X3),    0))/16.0);
  int maxx = (int) ceil((TOP(MAX(X1, X2, X3),xg*16))/16.0);
  int miny = (int) ceil((LOW(MIN(Y1, Y2, Y3),    0))/16.0);
  int maxy = (int) ceil((TOP(MAX(Y1, Y2, Y3),yg*16))/16.0);
  
  //--- Original, no clipping...
  //int minx = ( MIN(X1, X2, X3) + 0xF) >> 4;
  //int maxx = ( MAX(X1, X2, X3) + 0xF) >> 4;
  //int miny = ( MIN(Y1, Y2, Y3) + 0xF) >> 4;
  //int maxy = ( MAX(Y1, Y2, Y3) + 0xF) >> 4;    

  // Start in corner of 8x8 block
  minx &= ~(q - 1);
  miny &= ~(q - 1);
  
  // Half-edge constants
  int C1 = DY12 * X1 - DX12 * Y1;
  int C2 = DY23 * X2 - DX23 * Y2;
  int C3 = DY31 * X3 - DX31 * Y3;
  
  // Correct for fill convention
  if(DY12 < 0 || (DY12 == 0 && DX12 > 0)) C1++;
  if(DY23 < 0 || (DY23 == 0 && DX23 > 0)) C2++;
  if(DY31 < 0 || (DY31 == 0 && DX31 > 0)) C3++;


  // Loop through blocks
  for(y = miny; y < maxy; y += q){
    for(x = minx; x < maxx; x += q){
      
      // Corners of block
      int x0 = x << 4;
      int x1 = (x + q - 1) << 4;
      int y0 = y << 4;
      int y1 = (y + q - 1) << 4;
      
      // Evaluate half-space functions
      bool a00 = C1 + DX12 * y0 - DY12 * x0 > 0;
      bool a10 = C1 + DX12 * y0 - DY12 * x1 > 0;
      bool a01 = C1 + DX12 * y1 - DY12 * x0 > 0;
      bool a11 = C1 + DX12 * y1 - DY12 * x1 > 0;
      int  a   = (a00 << 0) | (a10 << 1) | (a01 << 2) | (a11 << 3);
      
      bool b00 = C2 + DX23 * y0 - DY23 * x0 > 0;
      bool b10 = C2 + DX23 * y0 - DY23 * x1 > 0;
      bool b01 = C2 + DX23 * y1 - DY23 * x0 > 0;
      bool b11 = C2 + DX23 * y1 - DY23 * x1 > 0;
      int  b   = (b00 << 0) | (b10 << 1) | (b01 << 2) | (b11 << 3);
      
      bool c00 = C3 + DX31 * y0 - DY31 * x0 > 0;
      bool c10 = C3 + DX31 * y0 - DY31 * x1 > 0;
      bool c01 = C3 + DX31 * y1 - DY31 * x0 > 0;
      bool c11 = C3 + DX31 * y1 - DY31 * x1 > 0;
      int  c   = (c00 << 0) | (c10 << 1) | (c01 << 2) | (c11 << 3);
      
      // Skip block when outside an edge
      if(a == 0x0 || b == 0x0 || c == 0x0) continue;
      
      // Accept whole block when totally covered
      if (a == 0xF && b == 0xF && c == 0xF){

	for(iy = y; iy < y+q; iy++){
	  for(ix = x; ix < x + q; ix++){      
	    ind = iy*xg + ix;
	    //    #pragma omp critical
	    {
	      float tempf = A*ix + B*iy + C;
	      if (tempf < 0) buffer[ind] += (w1+w2+w3)/3.0;
	      else           buffer[ind] += tempf;
	    }
	  } // end for ix
	} // end for iy
      } // end if
      
      else { // << // Partially covered block
	
	int CY1 = C1 + DX12 * y0 - DY12 * x0;
	int CY2 = C2 + DX23 * y0 - DY23 * x0;
	int CY3 = C3 + DX31 * y0 - DY31 * x0;
	for(iy = y; iy < y + q; iy++){
	  int CX1 = CY1;
	  int CX2 = CY2;
	  int CX3 = CY3;
	    
	  for(ix = x; ix < x + q; ix++){
	    if(CX1 > 0 && CX2 > 0 && CX3 > 0){
	      ind = iy*xg + ix;
	      //#pragma omp critical
	      {
		float tempf = A*ix + B*iy + C;
		if (tempf < 0) buffer[ind] += (w1+w2+w3)/3.0;
		else           buffer[ind] += tempf;
	      }
	    }
	    CX1 -= FDY12;
	    CX2 -= FDY23;
	    CX3 -= FDY31;
	  } // end for ix
	  CY1 += FDX12;
	  CY2 += FDX23;
	  CY3 += FDX31;
	} // for ix
      } // else
      
    } // x
  } // y 
  
} // raster_triangle()


//----------------------------------------------------------------
//
//----------------------------------------------------------------
void min_max(float a, float b, float c, float d, int *ind_minmax){
  return;
}


//----------------------------------------------------------------
//
//----------------------------------------------------------------
void sortFour(float a, float b, float c, float d, int *inds){
  float l1, h1, l2, h2;
  float L,m,M,H;
  
  int   il1, ih1, il2, ih2, iL,im,iM, iH;

  //--- First round bin 1 (a,b)
  if (a < b){
    l1  = a;
    h1  = b;
    il1 = 0;
    ih1 = 1;
  } else {
    l1  = b;
    h1  = a;
    il1 = 1;
    ih1 = 0;
  }

  //--- First round bin 2 (b,c)
  if (c < d){
    l2  = c;
    h2  = d;
    il2 = 2;
    ih2 = 3;
  } else {
    l2  = d;
    h2  = c;
    il2 = 3;
    ih2 = 2;
  }
  //--- Second round, lowest
  if (l1 < l2){
    L  = l1;
    m  = l2;
    iL = il1;
    im = il2;
  } else{
    L  = l2;
    m  = l1;
    iL = il2;
    im = il1;
  }
  //--- Second round, highest
  if (h1 > h2){
    H  = h1;
    M  = h2;
    iH = ih1;
    iM = ih2;
  } else{
    H  = h2;
    M  = h1;
    iH = ih2;
    iM = ih1;
  }
  //--- Third round middle
  if (m < M){
    inds[0] = iL;
    inds[1] = im;
    inds[2] = iM;
    inds[3] = iH;
  } else{
    inds[0] = iL;
    inds[1] = iM;
    inds[2] = im;
    inds[3] = iH;
  }

  
} // sortFour()

//----------------------------------------------------------
//
//----------------------------------------------------------
void save_cube(char *fname, float **cube, int nxyz, int nx, int ny, int nz){
  FILE *fd;
  char buf[256];
  int  i;
  double tempd; 
  float  tempf=0;
  int    tempi=0;
  char  remaining[256-36];

  //--- Load filename
  sprintf(buf,"%s",fname);
  
  printf("Writting densities to file '%s' ...\n",buf); fflush(stdout);
  fd=fopen(buf,"wb");
  //--- Header of file (number of partices)
  fwrite(&tempi,sizeof(int),1,fd);  
  fwrite(&tempi,sizeof(int),1,fd);
  fwrite(&nxyz,sizeof(int),1,fd);
  fwrite(&nx,sizeof(int),1,fd);
  fwrite(&ny,sizeof(int),1,fd);
  fwrite(&nz,sizeof(int),1,fd);
  fwrite(&tempf,sizeof(int),1,fd);
  fwrite(&tempf,sizeof(int),1,fd);
  fwrite(&tempf,sizeof(int),1,fd);
  fwrite(&remaining,sizeof(char),256-36,fd);
  
  
  for (int i=0;i<nz;i++){
    for (int j=0; j<nx*ny; j++) {
      tempf = cube[i][j];
      fwrite(&tempf,sizeof(float),1,fd);
    }
  }
  fclose(fd);
  printf("     Ready writting to file...\n");
  
} //--- end save_density()


//==================================================================
//           FUNCTIONS
//==================================================================

//----------------------------------------------------------
//---
//----------------------------------------------------------
double tetra_volu_ori(double xt[4], double yt[4], double zt[4]) {

  double v1[3] = {xt[1]-xt[0], yt[1]-yt[0], zt[1]-zt[0]};
  double v2[3] = {xt[2]-xt[0], yt[2]-yt[0], zt[2]-zt[0]};
  double v3[3] = {xt[3]-xt[0], yt[3]-yt[0], zt[3]-zt[0]};

  return ( v1[0]*v2[1]*v3[2] + v1[1]*v2[2]*v3[0] + v1[2]*v2[0]*v3[1] - 
	   (v1[2]*v2[1]*v3[0] + v1[1]*v2[0]*v3[2] + v1[0]*v2[2]*v3[1]) )/6.0;
  
} //--- end tetra_volu_ori()


//----------------------------------------------------------
//---
//----------------------------------------------------------
void fix_periodic_tetra(double x[4], double y[4], double z[4], double box) {

  double box2 = box/2.0;
  double dis;
  double min_x = box;
  double min_y = box;
  double min_z = box;
  double max_x = 0.0;
  double max_y = 0.0;
  double max_z = 0.0;

  for (int i=0; i<4; i++) {
    if (x[i] < min_x) min_x = x[i];
    if (y[i] < min_y) min_y = y[i];
    if (z[i] < min_z) min_z = z[i];
    if (x[i] > max_x) max_x = x[i];
    if (y[i] > max_y) max_y = y[i];
    if (z[i] > max_z) max_z = z[i];
  }

  double z0 = z[0];
  double z1 = z[1];
  double z2 = z[2];
  double z3 = z[3];
  //--- Decomposed (XYZ) distances from center of halo
  for (int w=0;w<4;w++){
    if ( (max_x-min_x > box2) && (x[w] > box2) ) x[w] -= box;
    if ( (max_y-min_y > box2) && (y[w] > box2) ) y[w] -= box;
    if ( (max_z-min_z > box2) && (z[w] > box2) ) z[w] -= box;
  } //--- end for w

} //--- end fix_periodic_tretra()


//----------------------------------------------------------
//--- Loads particle data from Gadget's default file
//----------------------------------------------------------
void load_snapshot(char *fname, int files) {
  FILE *fd;
  char buf[200];
  int  i,k,dummy,ntot_withmasses;
  int  n,pc,pc_new;
  int  count_err=0;

  size_t result;
#define SKIP result = fread(&dummy, sizeof(dummy), 1, fd);
 
  printf("reading \n");
  
  sprintf(buf,"%s",fname);      
  if(!(fd=fopen(buf,"r"))){
    printf("can't open file '%s'\n",buf);
    exit(0);
  }
  
  //--- Reader header
  SKIP;
  result = fread(&header1, sizeof(header1), 1, fd);  //--- Read Header1 in one pass
  SKIP;
  nPart = header1.npart[1];

  printf("header %d\n", nPart);
  allocate_memory_particles();
 
  //--- Read Particle's positions
  SKIP;
  for(n=0;n<header1.npart[1];n++) {
    result = fread(&(P[n].Pos[0]), sizeof(float), 3, fd);
  }
  SKIP;

  fclose(fd);

} //--- end load_snapshot()

//----------------------------------------------------------
//---Allocates the memory for the particle data.
//----------------------------------------------------------
void allocate_memory_particles(void) {
  printf("   Allocating memory for %d particles...",header1.npart[1] );
  if(!( P = (particle_data *) malloc(header1.npart[1]*sizeof(struct  particle_data) ))){
    fprintf(stderr,"failed to allocate memory HIGH RES.\n");
    exit(0);
  }
} //--- end  allocate_memory_particles()


//----------------------------------------------------------
//---Allocates the memory for the particle data.
//----------------------------------------------------------
void allocate_memory_cube(int nx, int ny, int nz) {
  printf("   Allocating memory for cube [%d,%d,%d]...", nx,ny,nz);
  //--- Allocate cube
  cube = (float **) malloc(sizeof(float **)*nz);
  for (int i=0; i<nz; i++) cube[i] = (float *) malloc(sizeof(float *) *nx*ny);  //--- OJO: ++i ???

  //--- Initialize array
  printf("\n>>> Initializing array...\n");
  for (int i=0; i<nz; i++) 
    for (int j=0; j<nx*ny; j++) cube[i][j]=0;
} //--- end  allocate_memory_cube()
