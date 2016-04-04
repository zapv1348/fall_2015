/*
 * dynamics.c
 *
 *    
 *
 *
 *
 *
 Thrust Vectored UAV
 *
 *      d_v     =   - D(v,alpha)/m - g*sin(theta-alpha) + ( cos(alpha)/m )*f_x + ( sin(alpha)/m )*f_z   
 *      d_alpha =   - L(v,alpha)/(m*v) + (g/v)*cos(theta-alpha) - ( sin(alpha)/(mv) )*f_x + ( cos(alpha)/(mv) )*f_z
 *      d_omega =     M(v,Alpha)/J + (l_tau/J)*f_z                  
 *      d_theta =     omega
 * 
 sliding car
 *
 *      d v     =    u1 cos(beta) + Fy(beta) sin(beta)
 *      d beta  = ( -u1 sin(beta) + Fy(beta) cos(beta) )/v - omega
 *      d omega = u2
 *
 *    the sideslip force (accel) model has parameters
 *
 *      Aslip - max sideslip accel
 *      Cfy   - cornering stiffness
 *    with
 *      Aslip <= Cfy <= 2 Aslip
 *
 *    model      
 *      Fy(beta) = -Aslip*sin( beta + Bslip/2*sin(2*beta) )
 *    with  Bslip = Cfy/Aslip - 1  so that  0 <= Bslip <= 1
 *
 *    note that
 *      Fy'(0) = -Cfy = -Aslip*(1 + Bslip)
 *
 *    the constraint on Cfy (hence Bslip) ensures that
 *      Fy'(beta) >= 0
 *    on (-pi/2,pi/2)
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
 *   f_x f_y
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

// Fy model parameters
//   Fy(beta) = -Aslip*sin( beta + Bslip/2*sin(2*beta) )

#define C_l         3.256
#define C_d_0       0.1716
#define C_d_alhpa   0.09999
// #define Aslip   (15.0)
// #define Cfy     (30.0)  // choose Cfy st Aslip <= Cfy <= 2 Aslip
// #define Bslip   (Cfy/Aslip - 1.0)
// #define Bslip2  (Bslip/2.0)

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
  double v, beta, omega;
  double u1, u2;

  // d/dt x = f(x,u)
  double dv, dbeta, domega;
  
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
  beta  = x[1];
  omega = x[2];

  u1   = u[0];
  u2   = u[1];

  v2 = v*v;
  v3 = v*v2;

  cos_beta = cos(beta);
  sin_beta = sin(beta);

  sin_2beta = sin(2.0*beta);
  cos_2beta = cos(2.0*beta);

  phi           = beta + Bslip2*sin_2beta;
  phi_beta      = 1.0  + Bslip *cos_2beta;
  phi_beta_beta =   -2.0*Bslip *sin_2beta;

  sin_phi = sin(phi);
  cos_phi = cos(phi);

  sin_phi_beta      = cos_phi*phi_beta;
  sin_phi_beta_beta = cos_phi*phi_beta_beta - sin_phi*phi_beta*phi_beta;

  // Fy(beta) = -Aslip*sin( beta + Bslip/2*sin(2*beta) )
  //          = -Aslip*sin( phi(beta) )
  Fy           = -Aslip*sin_phi;
  Fy_beta      = -Aslip*sin_phi_beta;
  Fy_beta_beta = -Aslip*sin_phi_beta_beta;


  // dynamics
  


  d_v     =   - D(v,alpha)/m - g*sin(theta-alpha) + ( cos(alpha)/m )*f_x + ( sin(alpha)/m )*f_z   
  d_alpha =   - L(v,alpha)/(m*v) + (g/v)*cos(theta-alpha) - ( sin(alpha)/(mv) )*f_x + ( cos(alpha)/(mv) )*f_z
  d_omega =     M(v,Alpha)/J + (l_tau/J)*f_z                  
  d_theta =     omega
  // dv     =    u1*cos_beta + Fy*sin_beta;
  // dbeta  = ( -u1*sin_beta + Fy*cos_beta )/v - omega;
  // domega = u2;

  // derivatives

  //dv         =  u1*cos_beta + Fy*sin_beta;
  dv_beta      = -u1*sin_beta + Fy*cos_beta + Fy_beta*sin_beta;
  dv_u1        =     cos_beta;

  dv_beta_beta = -u1*cos_beta - Fy*sin_beta + 2.0*Fy_beta*cos_beta
                              + Fy_beta_beta*sin_beta;
  dv_beta_u1   = -   sin_beta;

  //dbeta     =   ( -u1*sin_beta + Fy*cos_beta )/v - omega;
  dbeta_v     = - ( -u1*sin_beta + Fy*cos_beta )/v2;
  dbeta_beta  =   ( -u1*cos_beta - Fy*sin_beta + Fy_beta*cos_beta )/v;
  dbeta_omega =                                    - 1.0;
  dbeta_u1    =   ( -   sin_beta               )/v;

  //dbeta_v    = - ( -u1*sin_beta + Fy*cos_beta )/v2;
  dbeta_v_v    =   ( -u1*sin_beta + Fy*cos_beta )*2.0/v3;
  dbeta_v_beta = - ( -u1*cos_beta - Fy*sin_beta + Fy_beta*cos_beta )/v2;
  dbeta_v_u1   = - ( -   sin_beta               )/v2;

  //dbeta_beta    =   ( -u1*cos_beta - Fy*sin_beta + Fy_beta*cos_beta )/v;
  dbeta_beta_beta =   (  u1*sin_beta - Fy*cos_beta - 2.0*Fy_beta*sin_beta
			 + Fy_beta_beta*cos_beta )/v;
  dbeta_beta_u1   =   ( -   cos_beta                                  )/v;

  //domega  = u2;
  domega_u2 = 1.0;


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
    dx[1] = dbeta;
    dx[2] = domega;
  }

  if (do_y) {
    // example output: lateral acceleration
    // yy = v*(dbeta + omega);
    // y[0] = yy;

    // output
    //   side slip force  Fy(beta)
    //   lateral accel = v*(dbeta + omega)
    y[0] = Fy;
    y[1] = v * (dbeta + omega);
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
