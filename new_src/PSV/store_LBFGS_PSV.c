/*------------------------------------------------------------------------
 * Store gradients and model parameters for
 * Limited Memory - Broyden-Fletcher-Goldfarb-Shanno (l-BFGS)
 * optimization (Nocedal and Wright, 2006) in case of the PSV problem
 *
 * Daniel Koehn
 * Kiel, 25.07.2016
 *
 * ----------------------------------------------------------------------*/

#include "fd.h"

void store_LBFGS_PSV(float **taper_coeff, int nsrc, float **srcpos,
                     int **recpos, int ntr_glob, int iter, float **waveconv,
                     float **gradp, float **waveconv_u, float **gradp_u,
                     float **waveconv_rho, float **gradp_rho, float *y_LBFGS,
                     float *s_LBFGS, float *q_LBFGS, float **ppi, float **pu,
                     float **prho, int nxnyi, int LBFGS_pointer, int NLBFGS,
                     int NLBFGS_vec) {

  extern int NX, NY, IDX, IDY, SPATFILTER;
  extern int HESSIAN, SWS_TAPER_GRAD_VERT, SWS_TAPER_GRAD_HOR,
      SWS_TAPER_GRAD_SOURCES, SWS_TAPER_FILE;
  extern int POS[3], MYID;
  extern float C_vp, C_vs, C_rho;
  extern int INVMAT1;
  extern char JACOBIAN[STRING_SIZE];

  char jac[225], jac1[225];
  int i, j, k, h, h1, h2;
  float betaz, betan, gradplastiter, gradclastiter, betar, beta;
  float gamma_LBFGS, sum_nom, sum_denom;
  float LBFGSTMP, LBFGSTMP1, LBFGSTMP2, LBFGSTMP3, modellastiter, norm_fac,
      norm_fac_u, norm_fac_rho;
  float beta_LBFGS_1;
  float sc_vp, sc_vs, sc_rho;
  int ki, itershift, iter1;
  extern FILE *FP;
  FILE *FP3, *FP4, *FP6, *FP5, *FP7;

  itershift = 1;

  /* scaling of the l-BFGS variables
     INVMAT1 != 5: x = v / C_v       -> dJ/dx = C_v * dJ/dv,  s = dv / C_v
     INVMAT1 == 5: x = m = ln(v/v0)  -> already dimensionless, no scaling,
                   s = m - m_old = ln(v / v_old)  (v0 cancels out) */
  if (INVMAT1 == 5) {
    sc_vp = 1.0;
    sc_vs = 1.0;
    sc_rho = 1.0;
  } else {
    sc_vp = C_vp;
    sc_vs = C_vs;
    sc_rho = C_rho;
  }

  /* =====================================================================================================================================================
   */
  /* ===================================================== GRADIENT Vp/Zp/lambda
   * =========================================================================
   */
  /* =====================================================================================================================================================
   */

  /* Normalization of the gradient   */
  /* ------------------------------- */
  for (i = 1; i <= NX; i = i + IDX) {
    for (j = 1; j <= NY; j = j + IDY) {
      waveconv[j][i] = sc_vp * waveconv[j][i];
    }
  }

  for (i = 1; i <= NX; i = i + IDX) {
    for (j = 1; j <= NY; j = j + IDY) {
      gradp[j][i] = waveconv[j][i];
    }
  }

  /* =====================================================================================================================================================
   */
  /* ===================================================== GRADIENT Vs/Zs/mu
   * =============================================================================
   */
  /* =====================================================================================================================================================
   */

  /* Normalization of the gradient   */
  /* ------------------------------- */
  for (i = 1; i <= NX; i = i + IDX) {
    for (j = 1; j <= NY; j = j + IDY) {
      waveconv_u[j][i] = sc_vs * waveconv_u[j][i];
    }
  }

  for (i = 1; i <= NX; i = i + IDX) {
    for (j = 1; j <= NY; j = j + IDY) {
      gradp_u[j][i] = waveconv_u[j][i];
    }
  }

  /* =====================================================================================================================================================
   */
  /* ===================================================== GRADIENT rho
   * ==================================================================================
   */
  /* =====================================================================================================================================================
   */

  /* Normalization of the gradient   */
  /* ------------------------------- */
  for (i = 1; i <= NX; i = i + IDX) {
    for (j = 1; j <= NY; j = j + IDY) {
      waveconv_rho[j][i] = sc_rho * waveconv_rho[j][i];
    }
  }

  for (i = 1; i <= NX; i = i + IDX) {
    for (j = 1; j <= NY; j = j + IDY) {
      gradp_rho[j][i] = waveconv_rho[j][i];
    }
  }

  /* -------------------- */
  /* store l-BFGS vectors */
  /* -------------------- */

  if (iter > 1) {

    /* load old models and gradients - rho and store them in the LBFGS vectors
     */
    /* ------------------------------------------------------------------------
     */

    sprintf(jac, "%s_p_rho.old.%i.%i", JACOBIAN, POS[1], POS[2]);
    FP6 = fopen(jac, "rb");

    sprintf(jac1, "%s_p_mrho.old.%i.%i", JACOBIAN, POS[1], POS[2]);
    FP7 = fopen(jac1, "rb");

    /*iter1 = iter-itershift;*/ /* shift iter counter by 1 because L-BFGS method
                                   starts at iter > 1 */

    h = NLBFGS_vec * (LBFGS_pointer - 1) +
        1; /* locate current initial position in LBFGS-vector */

    for (i = 1; i <= NX; i = i + IDX) {
      for (j = 1; j <= NY; j = j + IDY) {

        /* calculate and save y, s at iteration step iter */
        fread(&gradplastiter, sizeof(float), 1, FP6);
        y_LBFGS[h] = waveconv_rho[j][i] - gradplastiter;

        fread(&modellastiter, sizeof(float), 1, FP7);
        if (INVMAT1 == 5) {
          /* s = ln(v/v_old); fluid (v = 0) -> no change */
          if ((prho[j][i] > 0.0) && (modellastiter > 0.0)) {
            s_LBFGS[h] = log(prho[j][i] / modellastiter);
          } else {
            s_LBFGS[h] = 0.0;
          }
        } else {
          s_LBFGS[h] = (prho[j][i] - modellastiter) / C_rho;
        }

        h++;
      }
    }

    fclose(FP6);
    fclose(FP7);

    /* load old models and gradients - Vs and store them in the LBFGS vectors */
    /* -----------------------------------------------------------------------
     */
    sprintf(jac, "%s_p_u.old.%i.%i", JACOBIAN, POS[1], POS[2]);
    FP6 = fopen(jac, "rb");

    sprintf(jac1, "%s_p_vs.old.%i.%i", JACOBIAN, POS[1], POS[2]);
    FP7 = fopen(jac1, "rb");

    for (i = 1; i <= NX; i = i + IDX) {
      for (j = 1; j <= NY; j = j + IDY) {

        /* calculate and save y, s at iteration step iter */
        fread(&gradplastiter, sizeof(float), 1, FP6);
        y_LBFGS[h] = waveconv_u[j][i] - gradplastiter;

        fread(&modellastiter, sizeof(float), 1, FP7);
        if (INVMAT1 == 5) {
          /* s = ln(v/v_old); fluid (v = 0) -> no change */
          if ((pu[j][i] > 0.0) && (modellastiter > 0.0)) {
            s_LBFGS[h] = log(pu[j][i] / modellastiter);
          } else {
            s_LBFGS[h] = 0.0;
          }
        } else {
          s_LBFGS[h] = (pu[j][i] - modellastiter) / C_vs;
        }

        h++;
      }
    }

    fclose(FP6);
    fclose(FP7);

    /* load old models and gradients - Vp and store them in the LBFGS vectors */
    /* -----------------------------------------------------------------------
     */
    sprintf(jac, "%s_p.old.%i.%i", JACOBIAN, POS[1], POS[2]);
    FP6 = fopen(jac, "rb");

    sprintf(jac1, "%s_p_vp.old.%i.%i", JACOBIAN, POS[1], POS[2]);
    FP7 = fopen(jac1, "rb");

    for (i = 1; i <= NX; i = i + IDX) {
      for (j = 1; j <= NY; j = j + IDY) {

        /* calculate and save y, s at iteration step iter */
        fread(&gradplastiter, sizeof(float), 1, FP6);
        y_LBFGS[h] = waveconv[j][i] - gradplastiter;

        fread(&modellastiter, sizeof(float), 1, FP7);
        if (INVMAT1 == 5) {
          /* s = ln(v/v_old); fluid (v = 0) -> no change */
          if ((ppi[j][i] > 0.0) && (modellastiter > 0.0)) {
            s_LBFGS[h] = log(ppi[j][i] / modellastiter);
          } else {
            s_LBFGS[h] = 0.0;
          }
        } else {
          s_LBFGS[h] = (ppi[j][i] - modellastiter) / C_vp;
        }

        h++;
      }
    }

    fclose(FP6);
    fclose(FP7);

    /* save q_LBFGS for all material parameters */
    h = 1;

    for (i = 1; i <= NX; i = i + IDX) {
      for (j = 1; j <= NY; j = j + IDY) {

        q_LBFGS[h] = waveconv_rho[j][i];
        h++;
      }
    }

    for (i = 1; i <= NX; i = i + IDX) {
      for (j = 1; j <= NY; j = j + IDY) {

        q_LBFGS[h] = waveconv_u[j][i];
        h++;
      }
    }

    for (i = 1; i <= NX; i = i + IDX) {
      for (j = 1; j <= NY; j = j + IDY) {

        q_LBFGS[h] = waveconv[j][i];
        h++;
      }
    }
  }
}
