/*
 * dynamics.c
 *
 *   Thrust Vectored Wing:
 *      Inputs: u1 = fx, u2=fz
 *      Dynamics:
 *      d v     = -D(v,alpha)/m -g*sin(theta-alpha)+(1/m)*cos(alpha)*u1 + (1/m)*sin(alpha)*u2   
 *      d alpha  = -L(v,alpha)/(m*v) + (g/m)*cos(theta-alpha)+omega - (1/m*v)*sin(alpha)*u1 + (1/m*v)*cos(alpha)*u2 
 *      d omega = (1/J)*M(v,alpha)+ (lz/J)*u2
 *      d theta = omega
 *
 *    Drag:
 *    D(v,alpha) = (1/2)*p*v^2*S*Cd(alpha)
 *    D(v,alpha) = (1/2)*p*v^2*S*(Cd0+Cd2*alpha^2)
 *    
 *    Lift:
 *    L(v,alpha) = (1/2)*p*v^2*S*Cl(alpha)
 *    L(v,alpha) = (1/2)*p*v^2*S*(cl1*alpha)
 *
 *    Moment:
 *    M(v,alpha) = (1/2)*p*v^2*c_bar*S*Cm(alpha)
 *    M(v,alpha) = (1/2)*p*v^2*c_bar*S*(cm1*alpha)
 *
 *      
 *      
 

dynamics(x,u,wt,ders, dx,y,  fxu_x_,fxu_u_, q,q_fxu_x_x_,q_fxu_x_u_,q_fxu_u_u_);
           ders:   dx(1),y(2), A(4),B(8),     Q(16),     S(32),     R(64)
 *       q  is the 'projection operator stabilized adjoint q'
 *
 * state             (x)
 *
 *   v alpha omega theta
 *
 * input             (u)
 *
 *   u1=fx u2=fz
 *
 * exogenous input   (wt)
 *
 *   none
 *
 * output            (y)
 *
 *   none
 *
 *
 * JH
 * nov 15 boulder
 */

/* prontoTK spec:
dynamics(x,u,wt,ders, dx,y,  fxu_x_,fxu_u_, q,q_fxu_x_x_,q_fxu_x_u_,q_fxu_u_u_);
           ders:   dx(1),y(2), A(4),B(8),     Q(16),     S(32),     R(64)
   cost(x,u,wlt,ders, lxu,   lxu_x_,lxu_u_, lxu_x_x_,lxu_x_u_,lxu_u_u_);
           ders:      lxu(1),  a(2),b(4),     Q(8),    S(16),   R(32)
*/

#include <math.h>

/*
dynamics(x,u,wt,ders, dx,y,  fxu_x_,fxu_u_, q,q_fxu_x_x_,q_fxu_x_u_,q_fxu_u_u_);
           ders:   dx(1),y(2), A(4),B(8),     Q(16),     S(32),     R(64)
*/
void dynamics(
  double *x, double *u, double *wt,
  int ders,
  double *dx, double *y,
  double *fxu_x_, double *fxu_u_,
  double *q,
  double *q_fxu_x_x_, double *q_fxu_x_u_, double *q_fxu_u_u_
  )
{

#define NX  (NS)

#define NU  (NI)

#define A(i,j)   fxu_x_[ ((i)-1) + ((j)-1)*(NX) ]  // columnwise
#define B(i,j)   fxu_u_[ ((i)-1) + ((j)-1)*(NX) ]  // columnwise

#define q_fxu_x_x(i,j)  q_fxu_x_x_[ ((i)-1) + ((j)-1)*(NX) ]  // columnwise
#define q_fxu_x_u(i,j)  q_fxu_x_u_[ ((i)-1) + ((j)-1)*(NX) ]  // columnwise
#define q_fxu_u_u(i,j)  q_fxu_u_u_[ ((i)-1) + ((j)-1)*(NU) ]  // columnwise

  int do_dx, do_A, do_B, do_Q, do_S, do_R, do_y;

  int i, j, k;

  // state and control
  double v, alpha, omega, gamma;
  double u1, u2;

  // d/dt x = f(x,u)
  double dv, dalpha, domega, dgamma;
  
  // partials: d/d xi and d/d ui of f(x,u)
  double  dv_v, dv_alpha, dv_gamma, dv_u1, dv_u2;
  double  dalpha_v, dalpha_alpha, dalpha_gamma, dalpha_omega, dalpha_u1,dalpha_u2;
  double domega_v, domega_alpha,domega_u2;
  double dgamma_v, dgamma_alpha, dgamma_gamma, dgamma_u1, dgamma_u2;
  
  
  // supporting cast

  double v2, v3;
  double alpha2;

  double sin_alpha, cos_alpha;
  double cos_gamma, sin_gamma;
  
  double L_ov,D,M;
  double L_ov_v, L_ov_alpha,L_ov_v_alpha;
  double D_v, D_alpha, D_v_alpha;
  double M_v, M_alpha, M_v_alpha;
 



  // determine what jobs to do
  /*
   *    ders - which ders? use binary code to specify
   *    1 - fxu
   *    2 - y
   *    4 - fxu_x = A
   *    8 - fxu_u = B
   *   16 - q_fxu_x_x = Q
   *   32 - q_fxu_x_u = S
   *   64 - q_fxu_u_u = R
   */
  do_dx =  ders & 1       ;
  do_y  = (ders & 2)  >> 1;
  do_A  = (ders & 4)  >> 2;
  do_B  = (ders & 8)  >> 3;
  do_Q  = (ders & 16) >> 4;
  do_S  = (ders & 32) >> 5;
  do_R  = (ders & 64) >> 6;

  // make states and controls more accessible

  v     = x[0];
  gamma = x[1];
  alpha = x[2];
  omega = x[3];

  u1   = u[0];
  u2   = u[1];

  v2 = v*v;
  v3 = v2*v;
  alpha2 = alpha*alpha;

  cos_gamma = cos(gamma);
  sin_gamma = sin(gamma);
  sin_alpha = sin(alpha);
  cos_alpha = cos(alpha);
  
#define Cd0   (0.1716)
#define Cd2   (2.395)
#define Cl1   (3.256)
#define Cm1   (-0.1)
#define g     (0.6)  
#define p     (1.20)
#define SS     (0.61)
#define mass  (12.0)
#define o_m   (1.0/mass)
#define c_bar (0.5)
#define J     (0.24)
#define lz    (0.31)

  D = 0.5*p*SS*v2*(Cd0+Cd2*alpha2);
  L_ov = 0.5*p*SS*Cl1*v*alpha; // L/v
  M = 0.5*p*SS*c_bar*Cm1*v2*alpha;
  
  D_v = p*SS*v*(Cd0+Cd2*alpha2);
  D_alpha = p*SS*v2*Cd2*alpha;
  D_v_alpha = 2*p*SS*v*Cd2*alpha;
  //D_alpha_v = D_v_alpha
  
  L_ov_v = 0.5*p*SS*Cl1*alpha;
  L_ov_alpha = 0.5*p*SS*Cl1*v;
  L_ov_v_alpha = 0.5*p*SS*Cl1;
  //L_ov_alpha_v = L_ov_v_alpha
  
  M_v = p*SS*c_bar*Cm1*v*alpha;
  M_alpha = 0.5*p*SS*c_bar*Cm1*v2;
  M_v_alpha = p*SS*c_bar*Cm1*v;
  //M_alpha_v = M_v_alpha
 
  
  ///////////////////// Dynamics/////////////////////////////
  dv = -D*o_m -g*sin_gamma + o_m*(u1*cos_alpha+u2*sin_alpha);
  dalpha = -o_m*L_ov + (g/v)*cos_gamma + omega +(o_m/v)*(-u1*sin_alpha+u2*cos_alpha);
  domega = M/J + (lz/J)*u2;
  dgamma = o_m*L_ov - (g/v)*cos_gamma + (o_m/v)*(u1*sin_alpha - u2*cos_alpha);

  ///////////////// derivatives //////////////////////////////
  
  ///////////////////dv////////////////////////////////////////
  dv_v = -o_m*D_v;
  dv_alpha = -o_m*D_alpha + (o_m)*(-u1*sin_alpha+u2*cos_alpha);
  //dv_omega = 0
  dv_gamma = -g*cos_gamma;
  
  dv_u1 = (o_m)*cos_alpha;
  dv_u2 = (o_m)*sin_alpha;
  
  ///////////////////dalpha /////////////////////////////////////////       
  dalpha_v = -o_m*L_ov_v - (g/v2)*cos_gamma - (o_m/v2)*(-u1*sin_alpha + u2*cos_alpha);
  dalpha_alpha = -o_m*L_ov_alpha + (o_m/v)*(-u1*cos_alpha - u2*sin_alpha);
  dalpha_omega = 1;
  dalpha_gamma = -(g/v)*sin_gamma;
  
  dalpha_u1 = -(o_m/v)*sin_alpha;
  dalpha_u2 = (o_m/v)*cos_alpha;
  
  /////////////////////domega////////////////////////////////////////////
  domega_v = M_v/J;
  domega_alpha = M_alpha/J;
  //domega_omega = 0;
  //domega_gamma = 0;
  
  //domega_u1 = 0;
  domega_u2 = lz/J;
  
  /////////// dgamma /////////////////////////////////////////////////
  
  dgamma_v = o_m*L_ov_v + (g/v2)*cos_gamma + (o_m/v2)*(-u1*sin_alpha + u2*cos_alpha);
  dgamma_alpha = o_m*L_ov_alpha + (o_m/v)*(u1*cos_alpha + u2*sin_alpha);
  //dgamma_omega = 0;
  dgamma_gamma = (g/v)*sin_gamma;
  
  dgamma_u1 = (o_m/v)*sin_alpha;
  dgamma_u2 = -(o_m/v)*cos_alpha;
  
  ///////////////// Second Derivatives //////////////////////////////

  
          
 
  if (do_A || do_Q || do_S) {
    // nothing here
  }

  if (do_A) {
    // zero out A
    for (i=0; i<NX*NX; i++) {
      fxu_x_[i] = 0.0;
    }
  }

  if (do_B) {
    // zero out B
    for (i=0; i<NX*NU; i++) {
      fxu_u_[i] = 0.0;
    }
  }

  if (do_Q) {
    // zero out Q = q_fxu_x_x
    for (i=0; i<NX*NX; i++) {
      q_fxu_x_x_[i] = 0.0;
    }
  }

  if (do_S) {
    // zero out S = q_fxu_x_u
    for (i=0; i<NX*NU; i++) {
      q_fxu_x_u_[i] = 0.0;
    }
  }

  if (do_R) {
    // zero out R = q_fxu_u_u
    for (i=0; i<NU*NU; i++) {
      q_fxu_u_u_[i] = 0.0;
    }
  }

  if (do_dx) {    
    dx[0] = dv;
    dx[1] = dgamma;
    dx[2] = dalpha;
    dx[3] = domega;
  }

  if (do_y) {

  }

  if (do_A) {
    
    A(1,1) = dv_v;
    A(1,2) = dv_gamma;
    A(1,3) = dv_alpha;
    //A(1,4) = 0;
    
    A(2,1) = dgamma_v;
    A(2,2) = dgamma_gamma;
    A(2,3) = dgamma_alpha;
    //A(2,4) = 0;
    
    A(3,1) = dalpha_v;
    A(3,2) = dalpha_gamma;
    A(3,3) = dalpha_alpha;
    A(3,4) = dalpha_omega;
    
    A(4,1) = domega_v;
    //A(4,2) = 0;
    A(4,3) = domega_alpha;
    //A(4,4) = 0;
  }

  if (do_B) {
    B(1,1) = dv_u1;
    B(1,2) = dv_u2;
    
    B(2,1) = dgamma_u1;
    B(2,2) = dgamma_u2;
    
    B(3,1) = dalpha_u1;
    B(3,2) = dalpha_u2;
    
    //B(4,1) = 0;
    B(4,2)=domega_u2;
  }

  if (do_Q) { // f_x_x terms
#if 0
    q_fxu_x_x(2,2) += q[0]*dv_beta_beta;

    q_fxu_x_x(1,1) += q[1]*dbeta_v_v;
    q_fxu_x_x(1,2) += q[1]*dbeta_v_beta;

    q_fxu_x_x(2,2) += q[1]*dbeta_beta_beta;

    // symmetrize Q
    q_fxu_x_x(2,1) = q_fxu_x_x(1,2);
#endif
  }

  if (do_S) { // f_x_u terms
#if 0
    q_fxu_x_u(2,1) += q[0]*dv_beta_u1;

    q_fxu_x_u(1,1) += q[1]*dbeta_v_u1;

    q_fxu_x_u(2,1) += q[1]*dbeta_beta_u1;
#endif
  }

  if (do_R) { // f_u_u terms
    // none
  }

#undef A
#undef B

#undef q_fxu_x_x
#undef q_fxu_x_u
#undef q_fxu_u_u

#undef NX
#undef NU
}



/*
 *  dynamics_unbdd(x,u,wt, dx);
 *
 *  trap to freeze a system that is becoming unbounded
 */
void dynamics_unbdd(
  double *x, double *u, double *wt,
  double *dx
  )
{
#define NX (NS)
#define NU (NI)

  int i;

// #if 0
//   can be used to "comment out" a section of code
#if 0
  // if pendulum angle exceeds 30 radians, just quit! 
  if ( x[0]*x[0] > 900.0 ) {
    for (i=0; i<NX; i++) {
      dx[i] = 0.0;
    }
  }
#endif

#undef NX
#undef NU
}
