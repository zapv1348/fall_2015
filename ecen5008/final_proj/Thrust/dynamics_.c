/*
 * dynamics.c
 *
 *   Thrust Vectored Wing:
 *      Inputs: u1 = fx, u2=fz
 *      Dynamics:
 *      d v     = -D(v,alpha)/m -g*sin(theta-alpha)+(1/m)*cos(alpha)*u1 + (1/m)*sin(alpha)*u2   
 *      d alpha  = -L(v,alpha)/(m*v) + (g/m)*cos(theta-alpha) - (1/m*v)*sin(alpha)*u1 + (1/m*v)*cos(alpha)*u2 
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

// Model Parameters

#define Cd0   (0.1716)
#define Cd2   (2.395)
#define Cl1   (3.256)
#define Cm1   (-0.1)
#define g     (9.81)
#define p     (1.20)
#define S     (0.61)
#define m     (12)
#define c_bar (0.5)
#define J     (0.24)
#define lz    (0.31)
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
  double v, alpha, omega, theta;
  double u1, u2;

  // d/dt x = f(x,u)
  double dv, dalpha, domega, dtheta;
  
  // partials: d/d xi and d/d ui of f(x,u)
  double          dv_beta;
  double dbeta_v, dbeta_beta, dbeta_omega;
  
  double dv_u1;
  double dbeta_u1;
  double           domega_u2;

  double dv_beta_beta, dv_beta_u1;
  double dbeta_v_v, dbeta_v_beta,    dbeta_v_u1;
  double            dbeta_beta_beta, dbeta_beta_u1;

  // supporting cast

  double v2, v3;

  double sin_beta, cos_beta;

  double sin_2beta, cos_2beta;

  double phi, phi_beta, phi_beta_beta;
  double sin_phi, cos_phi;
  double sin_phi_beta, sin_phi_beta_beta;

  double Fy, Fy_beta, Fy_beta_beta;



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
  alpha = x[1];
  omega = x[2];
  theta = x[3];

  u1   = u[0];
  u2   = u[1];

  v2 = v*v;
  v3 = v2*v;
  alpha2 = alpha*alpha

  cos_theta_alpha = cos(theta-alpha);
  sin_theta_alpha = sin(theta-alpha);
  sin_alpha = sin(alpha);
  cos_alpha = cos(alpha);
  
  D = 0.5*p*S*v2*(Cd0+Cd2*alpha2)
  L = 0.5*p*S*Cl1*v2*alpha;
  Lv = 0.5*p*S*Cl1*v*alpha;
  M = 0.5*p*S*c_bar*Cm1*v2*alpha;
  
  // Dynamics
  dv = -D/m -g*sin_theta_alpha + (1/m)*(u1*cos_alpha+u2*sin_alpha);
  dalpha = -L/(m*v) + (g/v)*cos_theta_alpha+(1/(m*v))*(-u1*sin_alpha+u2*cos_alpha)
  domega = M/J + (lz/J)*u2;
  dtheta = omega;

  k0 = -p*(S/m);

  // derivatives
  //dv
  
  dv_v = k0*(Cd0*v+Cd2*v*alpha2);
  dv_alpha = k0*Cd2*v2*alpha + g*cos_theta_alpha + (1/m)*(-u1*sin_alpha+u2*cos_alpha)
  //dv_omega = 0
  dv_theta = -g*cos_theta_alpha
  
  dv_u1 = (1/m)*cos_alpha;
  dv_u2 = (1/m)*sin_alpha;
  //dalpha
          
  dalpha_v = 0.5*k0*Cl1*alpha - (g/v2)*cos_theta_alpha + (1/(m*v2))*(u1*sin_alpha-u2*cos_alpha);
  dalpha_alpha = 0.5*k0*Cl1*v + (g/v)*sin_theta_alpha + (1/(m*v))*(-u1*cos_alpha - u2*sin_alpha);
  //dalpha_omega = 0;
  dalpha_theta = (-g/v)*sin_theta_alpha;
  
  dalpha_u1 = -(1/(m*v))*sin(alpha)
  dalpha_u2 = (1/(m*v))*cos(alpha)
  
  //domega
  
  domega_v = (p*S*c_bar/J)*Cm1*v*alpha;
  domega_alpha = (0.5*p*S*c_bar/J)*Cm1*v2;
  //domega_omega = 0;
  //domega_theta = 0;
  
  //domega_u1 = 0;
  domega_u2 = lz/J;
  
  //dtheta
  
  //dtheta_v = 0;
  //dtheta_alpha = 0;
  //dtheta_theta = 0;
  dtheta_omega = 1;
  
  //dtheta_u1 = 0;
  //dtheta_u2 = 0;
  
  /* **************** Second Derivatives fxu_x_x **************/
  
  // ## dv ## 
  dv_v_v = k0*(Cd0+Cd2*alpha2);
  dv_v_alpha = k0*(2*Cd2*v*alpha);
  //dv_v_omega = 0;
  dv_v_theta = -g*sin_theta_alpha;
  
  dv_alpha_v = k0*(2*Cd2*v*alpha);
  dv_alpha_alpha = k0*(2*v2*Cd2)+g*sin_theta_alpha + (1/m)*(-u1*sin_alpha-u2*cos_alpha);
  //dv_alpha_omega= 0;
  dv_alpha_theta = -g*sin_theta_alpha;
  
  //dv_omega_v = 0;
  //dv_omega_alpha = 0;
  //dv_omega_omega=0;
  //dv_omega_theta =0;
  
  //dv_theta_v = 0;
  dv_theta_alpha = -g*sin_theta_alpha;
  //dv_theta_omega=0;
  dv_theta_theta = g*sin_theta_alpha;
          
  // ## dalpha ##
  dalpha_v_v = (2*g/v3)*sin_theta_alpha - (2/(m*v3))*(-u1*cos_alpha - u2*sin_alpha)
  dalpha_v_alpha = 0.5*k0*Cl1 - (g/v2)*sin_theta_alpha + (1/(m*v2))*(u1*cos_alpha+u2*sin_alpha);
  //dalpha_v_omega = 0;
  dalpha_v_theta =  (g/v2)*sin_theta_alpha;
          
  dalpha_alpha_v = 0.5*k0*Cl1 - (g/v2)*sin_theta_alpha + (1/(m*v2))*(u1*cos_alpha + u2*sin_alpha);
  dalpha_alpha_alpha = -(g/v)*cos_theta_alpha + (1/(m*v))*(u1*sin_alpha - u2*cos_alpha);
  //dalpha_alpha_omega=0;
  dalpha_alpha_theta = (g/v)*cos_theta_alpha;
  
  //dalpha_omega_v = 0;
  //dalpha_omega_alpha = 0;
  //dalpha_omega_omega=0;
  //dalpha_omega_theta =0;
  
  dalpha_theta_v = (g/v2)*sin_theta_alpha;
  dalpha_theta_alpha = (g/v)*cos_theta_alpha;
  //dalpha_theta_omega= 0;
  dalpha_theta_theta = (-g/v)*cos_theta_alpha;
 
   // ## domega ##
  domega_v_v = (p*S*c_bar/J)*Cm1*alpha;
  domega_v_alpha = (p*S*c_bar/J)*Cm1*v;
  //domega_v_omega = 0;
  //domega_v_theta = 0;
          
  domega_alpha_v = (p*S*c_bar/J)*Cm1*v;
  //domega_alpha_alpha = 0;
  //domega_alpha_omega=0;
  //domega_alpha_theta =0;
  
  //domega_omega_v = 0;
  //domega_omega_alpha = 0;
  //domega_omega_omega=0;
  //domega_omega_theta =0;
  
  //domega_theta_v = 0;
  //domega_theta_alpha = 0;
  //domega_theta_omega=0;
  //domega_theta_theta =0;
          
  // ## dtheta ##
  /* all zero */
  /*  ******** Second Derivatives fxu_x_u ******** */
  //dv_v_u1 =0
  //dv_v_u2 =0
  dv_alpha_u1 = (1/m)*(-sin_alpha)
  dv_alpha_u2 = (1/m)*(cos_alpha)
  //dv_omega_u1 = 0;
  //dv_omega_u2 = 0;
  //dv_theta_u1 = 0;
  //dv_theta_u2 = 0;
  
  dalpha_v_u1 = (1/(m*v2))*(sin_alpha);
  dalpha_v_u2 = (1/(m*v2))*(-cos_alpha);
  dalpha_alpha_u1 = (1/(m*v))*(-cos_alpha);
  dalpha_alpha_u2 = (1/(m*v))*(-sin_alpha);
  //dalpha_omega_u1 = 0;
  //dalpha_omega_u2 = 0;
  //dalpha_theta_u1 = 0;
  //dalpha_theta_u2 = 0;
  
//   domega_v_u1 =0;
//   domega_v_u2 = 0;
//   domega_alpha_u1 =0;
//   domega_alpha_u2 = 0;
//   domega_omega_u1 = 0;
//   domega_omega_u2 = 0;
  //domega_theta_u1 = 0;
  //domega_theta_u2 = 0;
  
//   dtheta_v_u1 =0;
//   dtheta_v_u2 = 0;
//   dtheta_alpha_u1 =0;
//   dtheta_alpha_u2 = 0;
//   dtheta_omega_u1 =0;
//   dtheta_omega_u2 = 0;
  //dtheta_theta_u1 =0
  //dtheta_theta_u2 = 0
          
 
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
    dx[1] = dalpha;
    dx[2] = domega;
    dx[3] = dtheta;
  }

  if (do_y) {

  }

  if (do_A) {
    
    A(1,2) = dv_beta;
    
    A(2,1) = dbeta_v;
    A(2,2) = dbeta_beta;
    A(2,3) = dbeta_omega;
    
  }

  if (do_B) {
    B(1,1) = dv_u1;

    B(2,1) = dbeta_u1;

    B(3,2) = domega_u2;
  }

  if (do_Q) { // f_x_x terms
    q_fxu_x_x(2,2) += q[0]*dv_beta_beta;

    q_fxu_x_x(1,1) += q[1]*dbeta_v_v;
    q_fxu_x_x(1,2) += q[1]*dbeta_v_beta;

    q_fxu_x_x(2,2) += q[1]*dbeta_beta_beta;

    // symmetrize Q
    q_fxu_x_x(2,1) = q_fxu_x_x(1,2);
  }

  if (do_S) { // f_x_u terms
    q_fxu_x_u(2,1) += q[0]*dv_beta_u1;

    q_fxu_x_u(1,1) += q[1]*dbeta_v_u1;

    q_fxu_x_u(2,1) += q[1]*dbeta_beta_u1;
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
