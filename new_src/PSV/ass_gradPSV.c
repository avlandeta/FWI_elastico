/*
 * Assemble gradients for each shot (PSV problem)
 *
 * Daniel Koehn
 * Kiel, 24/04/2016
 */

#include "fd.h"
#include <math.h>

void ass_gradPSV(struct fwiPSV *fwiPSV, struct matPSV *matPSV, int iter) {

  /* global variables */
  extern int NX, NY, IDX, IDY, INVMAT1;
  extern int GRAD_FORM;
  extern int INV_VP_ITER, INV_VS_ITER, INV_RHO_ITER;
  extern float DT;

  /* local variables */
  int i, j;
  float muss, lamss;
  float Rps, Ip, Vp, Glam, Gmu, Grho;

  /* calculate gradient for first class */
  /* --------------------------------------- */

  /* norm(fwiPSV.waveconv_shot);
  norm(fwiPSV.waveconv_u_shot);
  norm(fwiPSV.waveconv_rho_shot); */

  for (i = 1; i <= NX; i = i + IDX) {
    for (j = 1; j <= NY; j = j + IDY) {

      /* calculate lambda gradient */
      (*fwiPSV).waveconv_lam[j][i] = -DT * (*fwiPSV).waveconv_shot[j][i];

      if (INVMAT1 == 3) {

        if (GRAD_FORM == 1) {

          muss = (*matPSV).pu[j][i];
          lamss = (*matPSV).ppi[j][i];

          if ((muss > 0.0) || (lamss > 0.0)) {
            (*fwiPSV).waveconv_lam[j][i] =
                (1.0 / (4.0 * (lamss + muss) * (lamss + muss))) *
                (*fwiPSV).waveconv_lam[j][i];
          }
        }

        (*fwiPSV).waveconv_shot[j][i] = (*fwiPSV).waveconv_lam[j][i];
      }

      if (INVMAT1 == 1 || INVMAT1 == 5) {

        /* calculate Vp gradient */
        if (GRAD_FORM == 1) {

          muss = (*matPSV).prho[j][i] * (*matPSV).pu[j][i] * (*matPSV).pu[j][i];
          lamss =
              (*matPSV).prho[j][i] * (*matPSV).ppi[j][i] * (*matPSV).ppi[j][i] -
              2.0 * muss;

          if ((muss > 0.0) || (lamss > 0.0)) {
            (*fwiPSV).waveconv_lam[j][i] =
                (1.0 / (4.0 * (lamss + muss) * (lamss + muss))) *
                (*fwiPSV).waveconv_lam[j][i];
          }

          (*fwiPSV).waveconv_shot[j][i] = 2.0 * (*matPSV).ppi[j][i] *
                                          (*matPSV).prho[j][i] *
                                          (*fwiPSV).waveconv_lam[j][i];
        }

        if (GRAD_FORM == 2) {

          muss = (*matPSV).prho[j][i] * (*matPSV).pu[j][i] * (*matPSV).pu[j][i];
          lamss =
              (*matPSV).prho[j][i] * (*matPSV).ppi[j][i] * (*matPSV).ppi[j][i] -
              2.0 * muss;

          if ((muss > 0.0) || (lamss > 0.0)) {
            (*fwiPSV).waveconv_lam[j][i] =
                (1.0 / (4.0 * (lamss + muss) * (lamss + muss))) *
                (*fwiPSV).waveconv_lam[j][i];
          }

          (*fwiPSV).waveconv_shot[j][i] = 2.0 * (*matPSV).ppi[j][i] *
                                          (*matPSV).prho[j][i] *
                                          (*fwiPSV).waveconv_lam[j][i];
        }
      }

      if (INVMAT1 == 2) {
        /* calculate Ip gradient */

        muss = (*matPSV).pu[j][i] * (*matPSV).pu[j][i] / (*matPSV).prho[j][i];
        lamss = ((*matPSV).ppi[j][i] * (*matPSV).ppi[j][i] -
                 2.0 * (*matPSV).pu[j][i] * (*matPSV).pu[j][i]) /
                (*matPSV).prho[j][i];

        if ((muss > 0.0) || (lamss > 0.0)) {
          (*fwiPSV).waveconv_lam[j][i] *=
              1.0 / (4.0 * (lamss + muss) * (lamss + muss));
        }

        (*fwiPSV).waveconv_shot[j][i] =
            (2.0 * (*matPSV).ppi[j][i] / (*matPSV).prho[j][i]) *
            (*fwiPSV).waveconv_lam[j][i];
      }

      if (INVMAT1 == 4) {
        /* ppi = Rps, pu = Ip, prho = Vp
           only scale lambda gradient here; Rps, Ip and Vp gradients need
           G_mu and G_rho, so they are assembled at the end */
        if ((*matPSV).ppi[j][i] != 0.0) {
          muss = (*matPSV).pu[j][i] * (*matPSV).prho[j][i] /
                 ((*matPSV).ppi[j][i] * (*matPSV).ppi[j][i]);
          lamss = (*matPSV).pu[j][i] * (*matPSV).prho[j][i] *
                  (1.0 - 2.0 / ((*matPSV).ppi[j][i] * (*matPSV).ppi[j][i]));
        } else {
          /* fluid: mu = 0, lambda = Ip*Vp */
          muss = 0.0;
          lamss = (*matPSV).pu[j][i] * (*matPSV).prho[j][i];
        }

        if ((muss > 0.0) || (lamss > 0.0)) {
          (*fwiPSV).waveconv_lam[j][i] *=
              1.0 / (4.0 * (lamss + muss) * (lamss + muss));
        }
      }

      if (iter < INV_VP_ITER) {
        (*fwiPSV).waveconv_shot[j][i] = 0.0;
      }
    }
  }

  /* calculate gradient for second class */
  /* ----------------------------------- */

  for (i = 1; i <= NX; i = i + IDX) {
    for (j = 1; j <= NY; j = j + IDY) {

      /* calculate mu gradient */
      (*fwiPSV).waveconv_mu[j][i] = -DT * (*fwiPSV).waveconv_u_shot[j][i];

      /* VERIFICAR valor INF en gradiente max */
      /*
      if (!isfinite((*fwiPSV).waveconv_mu[j][i])) {
        printf("INF/NAN MU: MYID=%d i=%d j=%d "
               "waveconv_mu=%e waveconv_u_shot=%e\n",
               i, j,
               (*fwiPSV).waveconv_mu[j][i],
               (*fwiPSV).waveconv_u_shot[j][i]);
      }
      */

      if (INVMAT1 == 1 || INVMAT1 == 5) {
        /* calculate Vs gradient */
        (*fwiPSV).waveconv_u_shot[j][i] =
            (-4.0 * (*matPSV).prho[j][i] * (*matPSV).pu[j][i] *
             (*fwiPSV).waveconv_lam[j][i]) +
            2.0 * (*matPSV).prho[j][i] * (*matPSV).pu[j][i] *
                (*fwiPSV).waveconv_mu[j][i];
      }

      if (INVMAT1 == 2) {
        /* calculate Is gradient */
        (*fwiPSV).waveconv_u_shot[j][i] =
            (-4.0 * (*matPSV).pu[j][i] / (*matPSV).prho[j][i]) *
                (*fwiPSV).waveconv_lam[j][i] +
            (2.0 * (*matPSV).pu[j][i] / (*matPSV).prho[j][i]) *
                (*fwiPSV).waveconv_mu[j][i];
      }

      /* VERIFICAR valor INF en gradiente max */
      /*
      if (!isfinite((*fwiPSV).waveconv_u_shot[j][i])) {
        printf("INF/NAN IS: MYID=%d i=%d j=%d "
               "Is=%e rho=%e grad_lam=%e grad_mu=%e grad_Is=%e\n",
               i, j,
               (*matPSV).pu[j][i],
               (*matPSV).prho[j][i],
               (*fwiPSV).waveconv_lam[j][i],
               (*fwiPSV).waveconv_mu[j][i],
               (*fwiPSV).waveconv_u_shot[j][i]);
      }
      */

      if (INVMAT1 == 3) {
        /* calculate u gradient */
        (*fwiPSV).waveconv_u_shot[j][i] = (*fwiPSV).waveconv_mu[j][i];
      }

      if (iter < INV_VS_ITER) {
        (*fwiPSV).waveconv_u_shot[j][i] = 0.0;
      }
    }
  }

  /* calculate gradient for third class */
  /* ------------------------------ */
  for (i = 1; i <= NX; i = i + IDX) {
    for (j = 1; j <= NY; j = j + IDY) {

      /* calculate density gradient rho' */
      (*fwiPSV).waveconv_rho_s[j][i] = -DT * (*fwiPSV).waveconv_rho_shot[j][i];

      if (INVMAT1 == 1 || INVMAT1 == 5) {
        /* calculate density gradient */
        (*fwiPSV).waveconv_rho_shot[j][i] =
            (((((*matPSV).ppi[j][i] * (*matPSV).ppi[j][i]) -
               (2.0 * (*matPSV).pu[j][i] * (*matPSV).pu[j][i])) *
              (*fwiPSV).waveconv_lam[j][i]) +
             ((*matPSV).pu[j][i] * (*matPSV).pu[j][i] *
              (*fwiPSV).waveconv_mu[j][i]) +
             (*fwiPSV).waveconv_rho_s[j][i]);
      }

      if (INVMAT1 == 3) {
        /* calculate density gradient */
        (*fwiPSV).waveconv_rho_shot[j][i] = (*fwiPSV).waveconv_rho_s[j][i];
      }

      if (INVMAT1 == 2) {
        /* calculate density gradient for impedance parametrization */
        (*fwiPSV).waveconv_rho_shot[j][i] =
            (*fwiPSV).waveconv_rho_s[j][i] -
            ((((*matPSV).ppi[j][i] * (*matPSV).ppi[j][i]) -
              (2.0 * (*matPSV).pu[j][i] * (*matPSV).pu[j][i])) /
             ((*matPSV).prho[j][i] * (*matPSV).prho[j][i])) *
                (*fwiPSV).waveconv_lam[j][i] -
            (((*matPSV).pu[j][i] * (*matPSV).pu[j][i]) /
             ((*matPSV).prho[j][i] * (*matPSV).prho[j][i])) *
                (*fwiPSV).waveconv_mu[j][i];
      }

      if (iter < INV_RHO_ITER) {
        (*fwiPSV).waveconv_rho_shot[j][i] = 0.0;
      }
    }
  }

  /* INVMAT1 == 4: Rps, Ip and Vp gradients (chain rule) */
  /* ------------------------------------------------- */
  /* rho = Ip/Vp, mu = Ip*Vp/Rps^2, lambda = Ip*Vp*(1 - 2/Rps^2)
     G_Rps = 4*Ip*Vp/Rps^3 * G_lam - 2*Ip*Vp/Rps^3 * G_mu
     G_Ip  = Vp*(1 - 2/Rps^2) * G_lam + Vp/Rps^2 * G_mu + 1/Vp * G_rho
     G_Vp  = Ip*(1 - 2/Rps^2) * G_lam + Ip/Rps^2 * G_mu - Ip/Vp^2 * G_rho */
  if (INVMAT1 == 4) {
    for (i = 1; i <= NX; i = i + IDX) {
      for (j = 1; j <= NY; j = j + IDY) {

        Rps = (*matPSV).ppi[j][i];
        Ip = (*matPSV).pu[j][i];
        Vp = (*matPSV).prho[j][i];

        Glam = (*fwiPSV).waveconv_lam[j][i];
        Gmu = (*fwiPSV).waveconv_mu[j][i];
        Grho = (*fwiPSV).waveconv_rho_s[j][i];

        if (Rps != 0.0) {
          /* Rps gradient */
          (*fwiPSV).waveconv_shot[j][i] =
              (4.0 * Ip * Vp / (Rps * Rps * Rps)) * Glam -
              (2.0 * Ip * Vp / (Rps * Rps * Rps)) * Gmu;

          /* Ip gradient */
          (*fwiPSV).waveconv_u_shot[j][i] =
              Vp * (1.0 - 2.0 / (Rps * Rps)) * Glam +
              (Vp / (Rps * Rps)) * Gmu;

          /* Vp gradient */
          (*fwiPSV).waveconv_rho_shot[j][i] =
              Ip * (1.0 - 2.0 / (Rps * Rps)) * Glam +
              (Ip / (Rps * Rps)) * Gmu;
        } else {
          /* fluid: mu = 0 fixed, Rps gradient undefined -> 0 */
          (*fwiPSV).waveconv_shot[j][i] = 0.0;
          (*fwiPSV).waveconv_u_shot[j][i] = Vp * Glam;
          (*fwiPSV).waveconv_rho_shot[j][i] = Ip * Glam;
        }

        /* density terms */
        if (Vp != 0.0) {
          (*fwiPSV).waveconv_u_shot[j][i] += Grho / Vp;
          (*fwiPSV).waveconv_rho_shot[j][i] -= (Ip / (Vp * Vp)) * Grho;
        }

        if (iter < INV_VP_ITER) {
          (*fwiPSV).waveconv_shot[j][i] = 0.0;
        }
        if (iter < INV_VS_ITER) {
          (*fwiPSV).waveconv_u_shot[j][i] = 0.0;
        }
        if (iter < INV_RHO_ITER) {
          (*fwiPSV).waveconv_rho_shot[j][i] = 0.0;
        }
      }
    }
  }

  /* INVMAT1 == 5: logarithmic (dimensionless) parameters */
  /* ---------------------------------------------------- */
  /* m_vp = ln(vp/vp0), m_vs = ln(vs/vs0), m_rho = ln(rho/rho0)
     dv/dm = v  ->  G_m = v * G_v   (v0 cancels out)
     The blocks above (shared with INVMAT1 == 1) already hold G_vp, G_vs
     and G_rho; here they are scaled by the current physical parameters:
       G_mvp  = 2*(lambda + 2*mu) * G_lam
       G_mvs  = 2*mu * (G_mu - 2*G_lam)
       G_mrho = lambda*G_lam + mu*G_mu + rho*G_rho'                       */
  if (INVMAT1 == 5) {
    for (i = 1; i <= NX; i = i + IDX) {
      for (j = 1; j <= NY; j = j + IDY) {

        (*fwiPSV).waveconv_shot[j][i] *= (*matPSV).ppi[j][i];
        (*fwiPSV).waveconv_u_shot[j][i] *= (*matPSV).pu[j][i];
        (*fwiPSV).waveconv_rho_shot[j][i] *= (*matPSV).prho[j][i];
      }
    }
  }
}
