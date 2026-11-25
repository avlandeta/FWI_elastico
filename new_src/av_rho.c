/* $Id: av_rho.c,v 1.1.1.1 2007/11/21 22:44:52 koehn Exp $*/

#include "fd.h"

void av_rho(float **rho, float **rip, float **rjp) {

  extern int NX, NY;
  int i, j;

  for (j = 1; j <= NY; j++) {
    for (i = 1; i <= NX; i++) {

      rjp[j][i] = 1.0 / (0.5 * (rho[j][i] + rho[j + 1][i]));
      rip[j][i] = 1.0 / (0.5 * (rho[j][i] + rho[j][i + 1]));

      if ((rho[j][i] < 1e-4) && (rho[j + 1][i] < 1e-4)) {
        rjp[j][i] = 0.0;
      }

      if ((rho[j][i] < 1e-4) && (rho[j][i + 1] < 1e-4)) {
        rip[j][i] = 0.0;
      }
    }
  }
}

/* PSV version: converts the third model slot to density according to
   INVMAT1 before averaging (see readmod_elastic_PSV.c for slot order)
   INVMAT1 = 1,2,3: rho slot already holds density
   INVMAT1 = 4    : pi = Rps, u = Ip, rho = Vp  ->  density = Ip / Vp */
void av_rho_PSV(float **pi, float **u, float **rho, float **rip, float **rjp) {

  extern int NX, NY, INVMAT1;
  int i, j;
  float r1, r2, r3;

  if (INVMAT1 != 4) {
    av_rho(rho, rip, rjp);
    return;
  }

  for (j = 1; j <= NY; j++) {
    for (i = 1; i <= NX; i++) {

      /* density = Ip / Vp (avoid division by zero) */
      if (rho[j][i] != 0.0) {
        r1 = u[j][i] / rho[j][i];
      } else {
        r1 = 0.0;
      }

      if (rho[j + 1][i] != 0.0) {
        r2 = u[j + 1][i] / rho[j + 1][i];
      } else {
        r2 = 0.0;
      }

      if (rho[j][i + 1] != 0.0) {
        r3 = u[j][i + 1] / rho[j][i + 1];
      } else {
        r3 = 0.0;
      }

      rjp[j][i] = 1.0 / (0.5 * (r1 + r2));
      rip[j][i] = 1.0 / (0.5 * (r1 + r3));

      if ((r1 < 1e-4) && (r2 < 1e-4)) {
        rjp[j][i] = 0.0;
      }

      if ((r1 < 1e-4) && (r3 < 1e-4)) {
        rip[j][i] = 0.0;
      }
    }
  }
}
