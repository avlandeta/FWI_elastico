/*------------------------------------------------------------------------
 *   output of material parameters after each stage
 *
 *   Daniel Koehn
 *   last update 9.11.2007
 *
 *  ---------------------------------------------------------------------*/

#include "fd.h"
void model_freq_out_PSV(float **ppi, float **rho, float **pu, int iter,
                        float freq) {

  /*--------------------------------------------------------------------------*/
  FILE *FP1;
  /* extern variables */
  extern int POS[3], MYID, INVMAT1;
  extern char INV_MODELFILE[STRING_SIZE];

  /* local variables */
  char modfile[STRING_SIZE];

  /* parameter names according to INVMAT1 (same order as readmod_elastic_PSV) */
  char name1[8] = "vp", name2[8] = "vs", name3[8] = "rho";

  if (INVMAT1 == 2) {
    sprintf(name1, "ip");
    sprintf(name2, "is");
    sprintf(name3, "rho");
  }

  if (INVMAT1 == 3) {
    sprintf(name1, "lam");
    sprintf(name2, "mu");
    sprintf(name3, "rho");
  }

  if (INVMAT1 == 4) {
    sprintf(name1, "rps");
    sprintf(name2, "ip");
    sprintf(name3, "vp");
  }

  sprintf(modfile, "%s_%s_stage_%d.bin", INV_MODELFILE, name1, iter);
  writemod(modfile, ppi, 3);
  MPI_Barrier(MPI_COMM_WORLD);

  if (MYID == 0)
    mergemod(modfile, 3);
  MPI_Barrier(MPI_COMM_WORLD);
  sprintf(modfile, "%s_%s_stage_%d.bin.%i.%i", INV_MODELFILE, name1, iter, POS[1],
          POS[2]);
  remove(modfile);

  sprintf(modfile, "%s_%s_stage_%d.bin", INV_MODELFILE, name2, iter);
  writemod(modfile, pu, 3);
  MPI_Barrier(MPI_COMM_WORLD);

  if (MYID == 0)
    mergemod(modfile, 3);
  MPI_Barrier(MPI_COMM_WORLD);
  sprintf(modfile, "%s_%s_stage_%d.bin.%i.%i", INV_MODELFILE, name2, iter, POS[1],
          POS[2]);
  remove(modfile);

  sprintf(modfile, "%s_%s_stage_%d.bin", INV_MODELFILE, name3, iter);
  writemod(modfile, rho, 3);
  MPI_Barrier(MPI_COMM_WORLD);

  if (MYID == 0)
    mergemod(modfile, 3);
  MPI_Barrier(MPI_COMM_WORLD);
  sprintf(modfile, "%s_%s_stage_%d.bin.%i.%i", INV_MODELFILE, name3, iter, POS[1],
          POS[2]);
  remove(modfile);
}
