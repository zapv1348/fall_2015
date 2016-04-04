/*
 *
 * QR_params_m.c
 *
 *   .MEX function to obtain the
 *
 *     cost paramaters  
 *
 *   in matlab.
 *
 * The calling syntax is:
 *
 *   [Qreg, Rreg] = QR_params_m();
 *
 *
 * JH, Boulder
 * jan 15
 */

#include "mex.h"

#include "sys_sizes.h"
#include "QR_params.h"

/* "QR_params.h" defines */
/* Qr - regulator state cost */
/* Rr - regulator input cost */

/* no input arguments */

/* Output Arguments */

#define	_Qr   plhs[0]
#define _Rr   plhs[1]

void mexFunction(
                 int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]
		 )
{
  double *qr, *rr;
  int i;


  // create and assign output variables

  // the cryptic code below fills in the diagonal elements of Q and R

  _Qr = mxCreateDoubleMatrix((mwSize)NS, (mwSize)NS, mxREAL);
  qr = mxGetPr(_Qr);
  for(i=0;i<NS;i++)
  {
    *qr = Qr[i];
    qr += NS+1;
  }

  if (nlhs > 1) {
    _Rr = mxCreateDoubleMatrix((mwSize)NI, (mwSize)NI, mxREAL);
    rr = mxGetPr(_Rr);
    for(i=0;i<NI;i++)
    {
      *rr = Rr[i];
      rr += NI+1;
    }
  }

}

