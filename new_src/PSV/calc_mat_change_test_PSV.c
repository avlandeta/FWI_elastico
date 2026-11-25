/*------------------------------------------------------------------------
 *   calculate test step length for material parameter update
 *
 *   Daniel Koehn
 *   last update 9.11.2007
 *
 *  ---------------------------------------------------------------------*/

#include "fd.h"
float calc_mat_change_test_PSV(float **waveconv, float **waveconv_rho,
                               float **waveconv_u, float **rho, float **rhonp1,
                               float **pi, float **pinp1, float **u,
                               float **unp1, int iter, int epstest,
                               float eps_scale, int itest) {

  /*--------------------------------------------------------------------------*/
  FILE *FP1;
  /* extern variables */
  extern float DH, DT;
  extern float EPSILON, EPSILON_u, EPSILON_rho, MUN, SCALERHO;
  extern int NX, NY, NXG, NYG, POS[3], MYID, INVMAT1;

  extern int INV_RHO_ITER, INV_VP_ITER, INV_VS_ITER;

  extern char INV_MODELFILE[STRING_SIZE];

  extern float VPUPPERLIM, VPLOWERLIM, VSUPPERLIM, VSLOWERLIM, RHOUPPERLIM,
      RHOLOWERLIM;

  /* local variables */

  float Rho, Vp, Vs, Vsnp1, Vpnp1, x, y, undf, r, pi0, K, mu, Zp, Zs, eps_true;
  float dpi, pimax, rhomax, umax, gradmax, gradmax_rho, gradmax_u, epsilon1,
      pimaxr, gradmaxr, gradmaxr_u, umaxr, gradmaxr_rho, rhomaxr;
  int i, j, ii, jj, testuplow;
  char modfile[STRING_SIZE];

  /*SCALERHO=0.5;*/

  /* invert for Zp and Zs ???  */
  /* ------------------------------------------------------------------------------------
   */

  /* find maximum of first parameter and gradient waveconv */

  pimax = 0.0;
  gradmax = 0.0;

  for (i = 1; i <= NX; i++) {
    for (j = 1; j <= NY; j++) {

      Zp = pi[j][i];

      if (Zp > pimax) {
        pimax = Zp;
      }

      if ((i * j == 1) || (gradmax == 0.0)) {
        gradmax = fabs(waveconv[j][i]);
      } else {
        if (fabs(waveconv[j][i]) > gradmax) {
          gradmax = fabs(waveconv[j][i]);
        }
      }
    }
  }

  /* find maximum of second parameter and gradient waveconv_u */

  umax = 0.0;
  gradmax_u = 0.0;

  for (i = 1; i <= NX; i++) {
    for (j = 1; j <= NY; j++) {

      Zs = u[j][i];

      if (Zs > umax) {
        umax = Zs;
      }

      if ((i * j == 1) || (gradmax_u == 0.0)) {
        gradmax_u = fabs(waveconv_u[j][i]);
      } else {
        if (fabs(waveconv_u[j][i]) > gradmax_u) {
          gradmax_u = fabs(waveconv_u[j][i]);
        }
      }
    }
  }

  /* find maximum of rho (thrid parameter) and gradient waveconv_rho */
  rhomax = 0.0;
  gradmax_rho = 0.0;

  for (i = 1; i <= NX; i++) {
    for (j = 1; j <= NY; j++) {

      if (rho[j][i] > rhomax) {
        rhomax = rho[j][i];
      }

      if ((i * j == 1) || (gradmax_rho == 0.0)) {
        gradmax_rho = fabs(waveconv_rho[j][i]);
      } else {
        if (fabs(waveconv_rho[j][i]) > gradmax_rho) {
          gradmax_rho = fabs(waveconv_rho[j][i]);
        }
      }
    }
  }

  /* calculate scaling factor for the gradients */
  /* --------------------------------------------- */

  /* parameter 1 */

  MPI_Allreduce(&pimax, &pimaxr, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&gradmax, &gradmaxr, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);

  /* INVMAT1 == 5: m = ln(v/v0) is dimensionless -> max |dm| = eps_scale
     (max relative change of the physical parameter ~ eps_scale) */
  if (INVMAT1 == 5) {
    EPSILON = eps_scale / gradmaxr;
  } else {
    EPSILON = eps_scale * (pimaxr / gradmaxr);
  }
  if (iter < INV_VP_ITER) {
    EPSILON = 0.0;
  }
  epsilon1 = EPSILON;
  MPI_Allreduce(&EPSILON, &epsilon1, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
  EPSILON = epsilon1;

  /* parameter 2 */
  MPI_Allreduce(&umax, &umaxr, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&gradmax_u, &gradmaxr_u, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);

  if (INVMAT1 == 5) {
    EPSILON_u = eps_scale / gradmaxr_u;
  } else {
    EPSILON_u = eps_scale * (umaxr / gradmaxr_u);
  }
  if (iter < INV_VS_ITER) {
    EPSILON_u = 0.0;
  }
  epsilon1 = EPSILON_u;
  MPI_Allreduce(&EPSILON_u, &epsilon1, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);

  EPSILON_u = epsilon1;

  /* parameter 3 */

  MPI_Allreduce(&rhomax, &rhomaxr, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&gradmax_rho, &gradmaxr_rho, 1, MPI_FLOAT, MPI_MAX,
                MPI_COMM_WORLD);

  if (INVMAT1 == 5) {
    EPSILON_rho = (eps_scale / gradmaxr_rho) * SCALERHO;
  } else {
    EPSILON_rho = eps_scale * (rhomaxr / gradmaxr_rho) * SCALERHO;
  }
  if (iter < INV_RHO_ITER) {
    EPSILON_rho = 0.0;
  }
  epsilon1 = EPSILON_rho;
  MPI_Allreduce(&EPSILON_rho, &epsilon1, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);

  EPSILON_rho = epsilon1;

  if (MYID == 0) {
    /* 1st, 2nd and 3rd model slots (pi, u, rho); their meaning depends on
       INVMAT1 (see readmod_elastic_PSV.c) */
    printf("MYID = %d \t 1st_Class_max = %e \t 1st_Class_gradmax = %e \n",
           MYID, pimaxr, gradmaxr);
    printf("MYID = %d \t 2nd_Class_max = %e \t 2nd_Class_gradmax = %e \n",
           MYID, umaxr, gradmaxr_u);
    printf("MYID = %d \t 3rd_Class_max = %e \t 3rd_Class_gradmax = %e \n",
           MYID, rhomaxr, gradmaxr_rho);
  }

  /* save true step length */
  eps_true = EPSILON_u;

  /* loop over local grid */
  for (i = 1; i <= NX; i++) {
    for (j = 1; j <= NY; j++) {

      /* update lambda, mu, rho */
      if ((INVMAT1 == 3) || (INVMAT1 == 2) || (INVMAT1 == 1) ||
          (INVMAT1 == 4) || (INVMAT1 == 5)) {

        testuplow = 0;

        if (INVMAT1 == 5) {
          /* m_new = m - eps*G_m  <=>  v_new = v * exp(-eps*G_m)
             (v0 cancels; v keeps its sign and Vs = 0 stays 0) */
          pinp1[j][i] = pi[j][i] * exp(-EPSILON * waveconv[j][i]);
          unp1[j][i] = u[j][i] * exp(-EPSILON_u * waveconv_u[j][i]);
          rhonp1[j][i] = rho[j][i] * exp(-EPSILON_rho * waveconv_rho[j][i]);
        } else {
          pinp1[j][i] = pi[j][i] - EPSILON * waveconv[j][i];
          unp1[j][i] = u[j][i] - EPSILON_u * waveconv_u[j][i];
          rhonp1[j][i] = rho[j][i] - EPSILON_rho * waveconv_rho[j][i];
        }

        if (INVMAT1 == 1 || INVMAT1 == 5) {
          if (pinp1[j][i] < VPLOWERLIM) {
            pinp1[j][i] = pi[j][i];
          }

          if (pinp1[j][i] > VPUPPERLIM) {
            pinp1[j][i] = pi[j][i];
          }

          if ((unp1[j][i] < VSLOWERLIM) && (unp1[j][i] > 1e-6)) {
            unp1[j][i] = u[j][i];
          }

          if (unp1[j][i] > VSUPPERLIM) {
            unp1[j][i] = u[j][i];
          }

          if (rhonp1[j][i] < RHOLOWERLIM) {
            rhonp1[j][i] = rho[j][i];
          }

          if (rhonp1[j][i] > RHOUPPERLIM) {
            rhonp1[j][i] = rho[j][i];
          }
        }

        if (INVMAT1 == 2) {
          /* pi = Ip, u = Is, rho = density
             Vp = Ip/rho, Vs = Is/rho */

          /* density limits -> revert rho (checked first, Vp and Vs depend on it) */
          if ((rhonp1[j][i] < RHOLOWERLIM) || (rhonp1[j][i] > RHOUPPERLIM)) {
            rhonp1[j][i] = rho[j][i];
          }

          if (rhonp1[j][i] > 0.0) {
            Vp = pinp1[j][i] / rhonp1[j][i];
            Vs = unp1[j][i] / rhonp1[j][i];

            /* Vp limits -> revert Ip */
            if ((Vp < VPLOWERLIM) || (Vp > VPUPPERLIM)) {
              pinp1[j][i] = pi[j][i];
            }

            /* Vs limits -> revert Is (Vs ~ 0 allowed for fluids) */
            if ((Vs < VSLOWERLIM) && (Vs > 1e-6)) {
              unp1[j][i] = u[j][i];
            }
            if (Vs > VSUPPERLIM) {
              unp1[j][i] = u[j][i];
            }
          }
        }

        if (INVMAT1 == 4) {
          /* pi = Rps, u = Ip, rho = Vp
             Vs = Vp/Rps, density = Ip/Vp */

          /* Vp limits -> revert Vp (checked first, Vs and density depend on it) */
          if ((rhonp1[j][i] < VPLOWERLIM) || (rhonp1[j][i] > VPUPPERLIM)) {
            rhonp1[j][i] = rho[j][i];
          }

          if (rhonp1[j][i] > 0.0) {
            Rho = unp1[j][i] / rhonp1[j][i];

            /* density limits -> revert Ip */
            if ((Rho < RHOLOWERLIM) || (Rho > RHOUPPERLIM)) {
              unp1[j][i] = u[j][i];
            }
          }

          /* Rps = 0 -> fluid (Vs = 0), skip Vs limits */
          if (pinp1[j][i] != 0.0) {
            Vs = rhonp1[j][i] / pinp1[j][i];

            /* Vs limits -> revert Rps (Vs ~ 0 allowed for fluids) */
            if ((Vs < VSLOWERLIM) && (Vs > 1e-6)) {
              pinp1[j][i] = pi[j][i];
            }
            if (Vs > VSUPPERLIM) {
              pinp1[j][i] = pi[j][i];
            }
          }
        }

        /* None of these parameters should be smaller than zero */
        if (pinp1[j][i] < 0.0) {
          pinp1[j][i] = pi[j][i];
        }
        if (unp1[j][i] < 0.0) {
          unp1[j][i] = u[j][i];
        }
        if (rhonp1[j][i] < 0.0) {
          rhonp1[j][i] = rho[j][i];
        }

        if (itest == 0) {
          rho[j][i] = rhonp1[j][i];
          pi[j][i] = pinp1[j][i];
          u[j][i] = unp1[j][i];
        }
      }
    }
  }

  if (itest == 0) {

    /* Write first model parameter */
    if (INVMAT1 == 1 || INVMAT1 == 5) {
        sprintf(modfile, "%s_vp.bin", INV_MODELFILE);
    } else if (INVMAT1 == 2) {
        sprintf(modfile, "%s_ip.bin", INV_MODELFILE);
    } else if (INVMAT1 == 3) {
        sprintf(modfile, "%s_lam.bin", INV_MODELFILE);
    } else if (INVMAT1 == 4) {
        sprintf(modfile, "%s_rps.bin", INV_MODELFILE);
    }

    writemod(modfile, pinp1, 3);

    MPI_Barrier(MPI_COMM_WORLD);

    if (MYID == 0)
        mergemod(modfile, 3);


    /* Write second model parameter */
    if (INVMAT1 == 1 || INVMAT1 == 5) {
        sprintf(modfile, "%s_vs.bin", INV_MODELFILE);
    } else if (INVMAT1 == 2) {
        sprintf(modfile, "%s_is.bin", INV_MODELFILE);
    } else if (INVMAT1 == 3) {
        sprintf(modfile, "%s_mu.bin", INV_MODELFILE);
    } else if (INVMAT1 == 4) {
        sprintf(modfile, "%s_ip.bin", INV_MODELFILE);
    }

    writemod(modfile, unp1, 3);

    MPI_Barrier(MPI_COMM_WORLD);

    if (MYID == 0)
        mergemod(modfile, 3);


    /* Write third model parameter */
    if (INVMAT1 == 4) {
        sprintf(modfile, "%s_vp.bin", INV_MODELFILE);
    } else {
        sprintf(modfile, "%s_rho.bin", INV_MODELFILE);
    }

    writemod(modfile, rho, 3);

    MPI_Barrier(MPI_COMM_WORLD);

    if (MYID == 0)
        mergemod(modfile, 3);
}

  return eps_true;
}
