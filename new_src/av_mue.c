
#include "fd.h"

void av_mue(float **u, float **uipjp, float **rho, float **pi) {

  extern int NX, NY, INVMAT1;
  int i, j;
  float u1, u2, u3, u4;

  if (INVMAT1 == 3) {

    for (j = 1; j <= NY; j++) {
      for (i = 1; i <= NX; i++) {

        uipjp[j][i] = 4.0 / ((1.0 / u[j][i]) + (1.0 / u[j][i + 1]) +
                             (1.0 / u[j + 1][i]) + (1.0 / u[j + 1][i + 1]));

        if ((u[j][i] == 0.0) || (u[j][i + 1] == 0.0) || (u[j + 1][i] == 0.0) ||
            (u[j + 1][i + 1] == 0.0)) {
          uipjp[j][i] = 0.0;
        }
      }
    }
  }

  if (INVMAT1 == 1 || INVMAT1 == 5) {

    for (j = 1; j <= NY; j++) {
      for (i = 1; i <= NX; i++) {

        u1 = rho[j][i] * u[j][i] * u[j][i];
        u2 = rho[j][i + 1] * u[j][i + 1] * u[j][i + 1];
        u3 = rho[j + 1][i] * u[j + 1][i] * u[j + 1][i];
        u4 = rho[j + 1][i + 1] * u[j + 1][i + 1] * u[j + 1][i + 1];

        uipjp[j][i] = 4.0 / ((1.0 / u1) + (1.0 / u2) + (1.0 / u3) + (1.0 / u4));

        if ((u1 == 0.0) || (u2 == 0.0) || (u3 == 0.0) || (u4 == 0.0)) {
          uipjp[j][i] = 0.0;
        }
      }
    }
  }

  if (INVMAT1 == 2) {
	for (j = 1; j <= NY; j++) {
     	 for (i = 1; i <= NX; i++) {

  	  u1 = (u[j][i] * u[j][i]) / rho[j][i];

 	   u2 = (u[j][i + 1] * u[j][i + 1]) / rho[j][i + 1];

  	  u3 = (u[j + 1][i] * u[j + 1][i]) / rho[j + 1][i];

	u4 = (u[j + 1][i + 1] * u[j + 1][i + 1]) / rho[j + 1][i + 1];

    uipjp[j][i] =   4.0 / ((1.0 / u1) + (1.0 / u2) +   (1.0 / u3) + (1.0 / u4));
    if ((u1 == 0.0) || (u2 == 0.0) || (u3 == 0.0) || (u4 == 0.0)) {
          uipjp[j][i] = 0.0;
        }
    
      }
    }
  }

  /* INVMAT1 == 4: pi = Rps, u = Ip, rho = Vp
     mu = rho_d * Vs^2 = (Ip/Vp) * (Vp/Rps)^2 = Ip * Vp / Rps^2 */
  if (INVMAT1 == 4) {

    for (j = 1; j <= NY; j++) {
      for (i = 1; i <= NX; i++) {

        /* Rps = 0 -> fluid: mu = 0 */
        if (pi[j][i] != 0.0) {
          u1 = u[j][i] * rho[j][i] / (pi[j][i] * pi[j][i]);
        } else {
          u1 = 0.0;
        }

        if (pi[j][i + 1] != 0.0) {
          u2 = u[j][i + 1] * rho[j][i + 1] / (pi[j][i + 1] * pi[j][i + 1]);
        } else {
          u2 = 0.0;
        }

        if (pi[j + 1][i] != 0.0) {
          u3 = u[j + 1][i] * rho[j + 1][i] / (pi[j + 1][i] * pi[j + 1][i]);
        } else {
          u3 = 0.0;
        }

        if (pi[j + 1][i + 1] != 0.0) {
          u4 = u[j + 1][i + 1] * rho[j + 1][i + 1] /
               (pi[j + 1][i + 1] * pi[j + 1][i + 1]);
        } else {
          u4 = 0.0;
        }

        uipjp[j][i] = 4.0 / ((1.0 / u1) + (1.0 / u2) + (1.0 / u3) + (1.0 / u4));

        if ((u1 == 0.0) || (u2 == 0.0) || (u3 == 0.0) || (u4 == 0.0)) {
          uipjp[j][i] = 0.0;
        }
      }
    }
  }
}
