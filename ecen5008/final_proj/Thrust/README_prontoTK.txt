brief guide to the pronto toolkit zoo of files
  [ JH nov 15 boulder ]

QUICK START
  run

    mex_newt_files

    newt_slcar

  to do trajectory optimization for the sliding car model

  change "verbose" in newt_slcar.m  to control plotting and pausing


USER needs to supply

  vector field & derivs
    dynamics.c

  incremental cost
    cost.c
    (if L2, set Q, R in cost_params.h)

  system size description
    sys_sizes.h

  TV LQR weights for projection operator feedback
    QR_params.h


USE

  mex_newt_files

to 'make' the S-functions
  (for simulations, there is an associated mdl, e.g., nonl_KS.mdl)

  proj operator + cost
    nonl_K.c

  regulator
    LQR_Kr.c

  descent direction calcs
    Prq1_Kv.c
    Prq2_Kv.c
    linK_zv.c

and the MEX functions (for evaluation at command line or in scripts)
  [ the suffix "_m"  means accessible from M file scripts, etc. ]

  evaluation of system vector field & derivs
    dynamics_m.c

  evaluation of incremental cost & derivs
    cost_m.c

  access to system sizeinfo
    sys_sizes_m.c

  access to QR weights
    QR_params_m.c


You can obtain calling sequences for the  _m  MEX functions using

  help dynamics_m
  help cost_m
  help sys_sizes_m
  help QR_params_m

which uses the M version comments, e.g.,   dynamics_m.m




