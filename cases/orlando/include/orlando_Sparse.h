/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/* Sparse Data Header File                                          */
/*                                                                  */
/* Generated by KPP-2.2.3 symbolic chemistry Kinetics PreProcessor  */
/*       (http://www.cs.vt.edu/~asandu/Software/KPP)                */
/* KPP is distributed under GPL, the general public licence         */
/*       (http://www.gnu.org/copyleft/gpl.html)                     */
/* (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa           */
/* (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech            */
/*     With important contributions from:                           */
/*        M. Damian, Villanova University, USA                      */
/*        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany */
/*                                                                  */
/* File                 : orlando_Sparse.h                          */
/* Time                 : Wed Aug 19 11:17:37 2020                  */
/* Working directory    : /home/WUR/krol005/kpp/examples            */
/* Equation file        : orlando.kpp                               */
/* Output root filename : orlando                                   */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */




/*  ----------> Sparse Jacobian Data                                */

extern int LU_IROW[LU_NONZERO];                 /* Row indexes of the LU Jacobian of variables */
extern int LU_ICOL[LU_NONZERO];                 /* Column indexes of the LU Jacobian of variables */
extern int LU_CROW[CNVAR];                      /* Compressed row indexes of the LU Jacobian of variables */
extern int LU_DIAG[CNVAR];                      /* Diagonal indexes of the LU Jacobian of variables */


/*  ----------> Sparse Hessian Data                                 */

extern int IHESS_I[NHESS];                      /* Index i of Hessian element d^2 f_i/dv_j.dv_k */
extern int IHESS_J[NHESS];                      /* Index j of Hessian element d^2 f_i/dv_j.dv_k */
extern int IHESS_K[NHESS];                      /* Index k of Hessian element d^2 f_i/dv_j.dv_k */


/*  ----------> Sparse Stoichiometric Matrix                        */

extern double STOICM[NSTOICM];                  /* Stoichiometric Matrix in compressed column format */
extern int IROW_STOICM[NSTOICM];                /* Row indices in STOICM */
extern int CCOL_STOICM[CNEQN];                  /* Beginning of columns in STOICM */
extern int ICOL_STOICM[NSTOICM];                /* Column indices in STOICM */


/*  ----------> Sparse Data for Jacobian of Reactant Products       */

extern int ICOL_JVRP[NJVRP];                    /* Column indices in JVRP */
extern int IROW_JVRP[NJVRP];                    /* Row indices in JVRP */
extern int CROW_JVRP[CNEQN];                    /* Beginning of rows in JVRP */

