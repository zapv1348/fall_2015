#define S_FUNCTION_NAME LQR_Kr

/*
 * LQR_Kr.c
 *
 *  S-function to compute
 *    the time-varying Linear Quadratic Regulator Kr(.)
 *    for the dynamics described in
 *
 *      dynamics.c
 *
 *  ** runs BACKWARD in time
 *  (in matlab, we reverse the time trajectories ... :(  ))
 *  (and eat the '-' sign in -Pdot = A' P ... etc. )
 *  [latest matlab can probably go backwards, but stay the same]
 *
 *  John Hauser
 *  Nov/Dec 06
 *  Dec 11
 *
 *  states
 *	P - quadratic part of cost to go	n(n+1)/2
 *          lower triangle stored P(1,1) P(2,1) P(2,2) P(3,1) ...
 *          we use an array of pointers & a macro to access
 *
 *  inputs (in time reversed order)
 *	xi(NS+NI)   :	alf(NS)   mu(NI)    - trajectory for linearization
 *
 *  outputs
 *	Kr(NI*NS) :		LQR feedback gain
 *
 *  parameters
 *	none
 */

/* prontoTK spec:
dynamics(x,u,wt,ders, dx,y,  fxu_x_,fxu_u_, q,q_fxu_x_x_,q_fxu_x_u_,q_fxu_u_u_);
           ders:   dx(1),y(2), A(4),B(8),     Q(16),     S(32),     R(64)
   cost(x,u,wlt,ders, lxu,   lxu_x_,lxu_u_, lxu_x_x_,lxu_x_u_,lxu_u_u_);
           ders:      lxu(1),  a(2),b(4),     Q(8),    S(16),   R(32)
*/

#define S_FUNCTION_LEVEL 2

#include "simstruc.h"

#include <math.h>

// system description: dynamics with derivatives
#include "sys_sizes.h"
#include "dynamics.c"

// get regulator Q and R (and params)
#include "QR_params.h"

// "sys_sizes.h" defines
// NS   - number of system states
// NI   - number of system inputs
// NW   - number of exogenous inputs for dynamics
// NWL  - number of exogenous inputs for cost and constraints

#define NS2		(NS*NS)
#define NN12		((NS*(NS+1))/2)

#define N_STATES	(NN12)
#define N_INPUTS	(NS+NI+NW)
#define N_OUTPUTS	(NS*NI)
#define N_PARAMS	(0)

static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumSFcnParams(S, N_PARAMS);  /* Number of expected parameters */
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        /* Return if number of expected != number of actual parameters */
        return;
    }

    ssSetNumContStates(S, N_STATES);
    ssSetNumDiscStates(S, 0);

    if (!ssSetNumInputPorts(S, 1)) return;
    ssSetInputPortWidth(S, 0, N_INPUTS);
    ssSetInputPortDirectFeedThrough(S, 0, 1);

    if (!ssSetNumOutputPorts(S, 1)) return;
    ssSetOutputPortWidth(S, 0, N_OUTPUTS);

    ssSetNumSampleTimes(S, 1);
    ssSetNumRWork(S, 0);
    ssSetNumIWork(S, 0);
    ssSetNumPWork(S, 0);
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);

    ssSetOptions(S, 0);
}

static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
}

/* set initial conditions with  set options */
#if defined(MDL_INITIALIZE_CONDITIONS)
#undef MDL_INITIALIZE_CONDITIONS
#endif

/*
 * most matrices are stored COLUMNwise
 *  THE exception is the feedback gain matrix K
 *  which is stored ROWwise
 *
 * in each case, we set up macros to allow matlab-like indexing
 *   [ (1..n) vs (0..n-1) ]
 */

static void optBK(double *B_, double *Ko_,
		  double *x, double *u)
{

#define B(i,j)	B_[ ((i)-1) + ((j)-1)*(NS) ]	// columnwise
#define Ko(i,j)	Ko_[ ((i)-1)*NS + ((j)-1) ]	// rowwise
#define P(i,j)	(( (j<=i) ? P_[i-1][j-1] : P_[j-1][i-1] ))	// low tri

  double *P__ = x;

  // pointers to 'input structure'
  double *alf = u, *mu = alf+NS, *wt = mu+NI;

  double *P_[NS];       // array of pointers for lower triangular P

  int i, j, k;
  double ko;

  j = 0;	// set up lower triangular structure for P
  for (i=0; i<NS; i++) {
    P_[i] = P__ + j;
    j += i+1;
  }

  // get B of linearization about (alf(t),mu(t))

/*
dynamics(x,u,wt,ders, dx,y,  fxu_x_,fxu_u_, q,q_fxu_x_x_,q_fxu_x_u_,q_fxu_u_u_);
         ders:     dx(1),y(2), A(4),B(8),     Q(16),     S(32),     R(64)
  cost(x,u,wlt, ders, lxu,   lxu_x_,lxu_u_, lxu_x_x_,lxu_x_u_,lxu_u_u_);
                ders: lxu(1),  a(2),b(4),     Q(8),    S(16),   R(32)
*/

  // get B(8)
  dynamics(alf,mu,wt, 8, NULL,NULL,  NULL,B_,  NULL, NULL,NULL,NULL);

  // Kopt = R^{-1} B^T P
  //     this calc requires R to be diagonal
  for (i=1; i<=NI; i++) {
    for (j=1; j<=NS; j++) {
      ko = 0.0;
      for (k=1; k<=NS; k++) {
	ko += B(k,i)*P(k,j);	// B^T P
      }
      Ko(i,j) = Rri[i-1]*ko;	// R^{-1} B^T P
    }
  }

#undef P
#undef Ko
#undef B

}

static void mdlOutputs(SimStruct *S, int_T tid)
{
  real_T	*y	= ssGetOutputPortRealSignal(S,0);
  real_T	*x	= ssGetContStates(S);	
  InputRealPtrsType uPtrs   = ssGetInputPortRealSignalPtrs(S,0);
  // double        *params0 = mxGetPr(ssGetSFcnParam(S,0));

  int i;
  double u[N_INPUTS];
  double B_[NS*NI], *Ko_ = y;

  /* copy inputs to u[] to ensure contiguity ... and easy access*/
  for (i=0; i<N_INPUTS; i++)
    u[i] = *uPtrs[i];

  /* output new optimal gain */
  optBK(B_, Ko_, x, u);
}

/* no model update */
#if defined(MDL_UPDATE)
#undef MDL_UPDATE
#endif

#define MDL_DERIVATIVES
static void mdlDerivatives(SimStruct *S)
{
  real_T            *dx     = ssGetdX(S);
  real_T            *x      = ssGetContStates(S);
  InputRealPtrsType uPtrs   = ssGetInputPortRealSignalPtrs(S,0);
  // double         *param0 = mxGetPr(ssGetSFcnParam(S,0));

  int i, j, k;
  double u[N_INPUTS];

  // pointers to 'input structure'
  double *alf = u, *mu = alf+NS, *wt = mu+NI;

  // double a[NS], b[NI], Ko_b[NS], K_b[NS];

  double *P__ = x;
  double *dP__ = dx;
  double dp;
  double *P_[NS], *dP_[NS];	// array of pointers for lower triangular P, dP

  double A_[NS2], A_BKo_[NS2];
  double B_[NS*NI], Ko_[NI*NS];
  double Q_[NS2];

  /* index A & P with 1:NS  (NOT 0:NS-1!) */
#define A(i,j)		A_[((i)-1) + NS*((j)-1)]
#define Q(i,j)		Q_[((i)-1) + NS*((j)-1)]
#define P(i,j)	(( (j<=i) ? P_[i-1][j-1]  : P_[j-1][i-1] ))	// low tri
#define dP(i,j)	( dP_[i-1][j-1] )	// *** i <= j *** REQUIRED!
#define B(i,j)	B_[ ((i)-1) + ((j)-1)*(NS) ]	// columnwise
#define Ko(i,j)		Ko_[ ((i)-1)*NS + ((j)-1) ]	// rowwise

  // set up lower triangular structure for P and dP
  j = 0;
  for (i=0; i<NS; i++) {
    P_[i]  = P__ + j;
    dP_[i] = dP__ + j;
    j += i+1;
  }

  /* copy inputs to u[] to ensure contiguity ... and easy access*/
  for (i=0; i<N_INPUTS; i++)
    u[i] = *uPtrs[i];

/* prontoTK spec:
dynamics(x,u,wt,ders, dx,y,  fxu_x_,fxu_u_, q,q_fxu_x_x_,q_fxu_x_u_,q_fxu_u_u_);
           ders:   dx(1),y(2), A(4),B(8),     Q(16),     S(32),     R(64)
   cost(x,u,wlt,ders, lxu,   lxu_x_,lxu_u_, lxu_x_x_,lxu_x_u_,lxu_u_u_);
           ders:      lxu(1),  a(2),b(4),     Q(8),    S(16),   R(32)
*/

  // get A of linearization about (alf(t),mu(t))

  // get A(4)
  dynamics(alf,mu,wt, 4, NULL,NULL,  A_,NULL,  NULL, NULL,NULL,NULL);

  /* set up Q matrix --- start with zero */
  for (i=0; i<NS2; i++){
    Q_[i] = 0.0;
  }
  
  /* add diagonal cost fcn terms  */
  for (i=1; i<=NS; i++){
    Q(i,i) = Qr[i-1];
  }

  optBK(B_, Ko_, x, u);

  /* only *lower triangle* of dP */
  for (i=1; i<=NS; i++) {
    for (j=1; j<=i; j++) {
      dp = Q(i,j);
      for (k=1; k<=NS; k++) {
        dp += A(k,i)*P(k,j) + P(i,k)*A(k,j);    // A^T P + P A
      }
      for (k=1; k<=NI; k++) {
        dp -= Ko(k,i)*Rr[k-1]*Ko(k,j);  // -P B R^-1 B^T P = -Ko^T R Ko
      }
      dP(i,j) /* = dP(j,i) */  = dp;
    }
  }

#undef B
#undef K
#undef Ko
#undef dP
#undef P
#undef A
}

static void mdlTerminate(SimStruct *S)
{
}

#ifdef	MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
