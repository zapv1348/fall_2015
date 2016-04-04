/*
 * cost.c
 *
 *    incremental cost - L_2 traj error (squared)
 *
 *      l(x,u,wlt)   [ w_l(t) = (x_des(t), u_des(t)) ]
 *            = || x - x_des(t) ||^2_Q / 2 +  | u - u_des(t) |^2_R / 2
 *
 *      diag Q and R are to be specified in L2_cost_params.h
 *
   cost(x,u,wlt,ders, lxu,   lxu_x_,lxu_u_, lxu_x_x_,lxu_x_u_,lxu_u_u_);
           ders:      lxu(1),  a(2),b(4),     Q(8),    S(16),   R(32)
 *
 * state             (x)
 *
 *   phi phid th thd
 *
 * input             (u)
 *
 *   u
 *
 * exogenous input   (wlt)
 *
 *   x_des u_des
 *
 * JH, Boulder, apr 12
 * checked oct 13
 * nov 15
 */

/* prontoTK spec:
dynamics(x,u,wt,ders, dx,y,  fxu_x_,fxu_u_, q,q_fxu_x_x_,q_fxu_x_u_,q_fxu_u_u_);
           ders:   dx(1),y(2), A(4),B(8),     Q(16),     S(32),     R(64)
   cost(x,u,wlt,ders, lxu,   lxu_x_,lxu_u_, lxu_x_x_,lxu_x_u_,lxu_u_u_);
           ders:      lxu(1),  a(2),b(4),     Q(8),    S(16),   R(32)
*/

#include <math.h>

#include "L2_cost_params.h"

/*
   cost(x,u,wlt,ders, lxu,   lxu_x_,lxu_u_, lxu_x_x_,lxu_x_u_,lxu_u_u_);
           ders:      lxu(1),  a(2),b(4),     Q(8),    S(16),   R(32)
*/
void cost(
  double *x, double *u, double *wlt,
  int ders,
  double *lxu,
  double *lxu_x_, double *lxu_u_,
  double *lxu_x_x_, double *lxu_x_u_, double *lxu_u_u_
  )
{

#define NX  (NS)

#define NU  (NI)

#define lxu_x_x(i,j)      lxu_x_x_[ ((i)-1) + ((j)-1)*(NX) ]  // columnwise
#define lxu_x_u(i,j)      lxu_x_u_[ ((i)-1) + ((j)-1)*(NX) ]  // columnwise
#define lxu_u_u(i,j)      lxu_u_u_[ ((i)-1) + ((j)-1)*(NU) ]  // columnwise
// #define lxu_x_x_safe(i,j) lxu_x_x_safe_[((i)-1) + ((j)-1)*(NX) ]

  double *x_des = wlt, *u_des=x_des+NX;
  double dx[NX], du[NU];

  int do_l, do_a, do_b, do_Q, do_S, do_R, do_Qsafe;

  int i, j, k;

  double ll;

  // determine what jobs to do
  /*
   *    ders - which ders? use binary code to specify
   *    1 - lxu
   *    2 - lxu_x = a
   *    4 - lxu_u = b
   *    8 - lxu_x_x = Q
   *   16 - lxu_x_u = S
   *   32 - lxu_u_u = R
   */
  do_l     =  ders & 1       ;
  do_a     = (ders & 2)  >> 1;
  do_b     = (ders & 4)  >> 2;
  do_Q     = (ders & 8)  >> 3;
  do_S     = (ders & 16) >> 4;
  do_R     = (ders & 32) >> 5;

  // make states and controls more accessible
  //   not needed here

  if (do_l || do_a || do_b) {
    for (i=0; i<NX; i++)
      dx[i] = x[i] - x_des[i];
    for (i=0; i<NU; i++)
      du[i] = u[i] - u_des[i];
  }

  if (do_a) {
    // zero out a
    for (i=0; i<NX; i++) {
      lxu_x_[i] = 0.0;
    }
  }

  if (do_b) {
    // zero out b
    for (i=0; i<NU; i++) {
      lxu_u_[i] = 0.0;
    }
  }

  if (do_Q) {
    // zero out Q = lxu_x_x
    for (i=0; i<NX*NX; i++) {
      lxu_x_x_[i] = 0.0;
    }
  }

  if (do_S) {
    // zero out S = lxu_x_u
    for (i=0; i<NX*NU; i++) {
      lxu_x_u_[i] = 0.0;
    }
  }

  if (do_R) {
    // zero out R = lxu_u_u
    for (i=0; i<NU*NU; i++) {
      lxu_u_u_[i] = 0.0;
    }
  }
  
  if (do_l) {
    ll = 0.0;
    for (i=0; i<NX; i++)
      ll += QD[i]*dx[i]*dx[i]/2.0;
    for (i=0; i<NU; i++)
      ll += RD[i]*du[i]*du[i]/2.0;
    lxu[0] = ll;
  }

  if (do_a) {
    for (i=0; i<NX; i++)
      lxu_x_[i] = QD[i]*dx[i];
  }

  if (do_b) {
    for (i=0; i<NU; i++)
      lxu_u_[i] = RD[i]*du[i];
  }

  if (do_Q) {
    for (i=1; i<=NX; i++)
      lxu_x_x(i,i) = QD[i-1];

#if 0
    printf("lxu_x_x:\n");
    for (i=1; i<=NX; i++) {
      for (j=1; j<=NX; j++) {
        printf("  %g", lxu_x_x(i,j));
      }
      printf("\n");
    }
#endif

  }

  if (do_S) {
    // no lxu_x_u_ terms --- YET!
  }

  if (do_R) {
    for (i=1; i<=NU; i++)
      lxu_u_u(i,i) = RD[i-1];
  }
  
#undef lxu_x_x
#undef lxu_x_u
#undef lxu_u_u
// #undef lxu_x_x_safe

#undef NX
#undef NU

}
