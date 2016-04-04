/*
 * sys_sizes.h
 *
 * specify the system sizes (NS, NI, etc.)
 *
 * JH
 * nov 15 boulder
 */

// flying wing
//   x: v alpha omega theta
//   u: u1 u2
//   ## y: Fy(beta)

#define NS    (4)        // number of system states
#define NI    (2)        // number of system inputs
#define NO    (0)        // number of system outputs
#define NW    (0)        // number of exogenous inputs for dynamics
#define NWL   (6)        // number of exogenous inputs for cost l(x,u,w_l(t))   
#define NWC   (0)        // number of constraints (eventually)
