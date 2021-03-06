/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/* The Stoichiometric Chemical Model File                           */
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
/* File                 : ifs_Stoichiom.c                           */
/* Time                 : Mon Apr 12 11:55:23 2021                  */
/* Working directory    : /home/WUR/krol005/kpp/examples            */
/* Equation file        : ifs.kpp                                   */
/* Output root filename : ifs                                       */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "ifs_Parameters.h"
#include "ifs_Global.h"
#include "ifs_Sparse.h"


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/* ReactantProd - Reactant Products in each equation                */
/*   Arguments :                                                    */
/*      V         - Concentrations of variable species (local)      */
/*      F         - Concentrations of fixed species (local)         */
/*      ARP       - Reactant product in each equation               */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void ReactantProd( 
  double V[],                            /* Concentrations of variable species (local) */
  double F[],                            /* Concentrations of fixed species (local) */
  double ARP[]                           /* Reactant product in each equation */
)
{

/* Reactant Products in each equation are useful in the             */
/*     stoichiometric formulation of mass action law                */
  ARP[0] = V[11]*F[2];
  ARP[1] = V[9]*V[11];
  ARP[2] = V[9]*F[2];
  ARP[3] = V[9]*V[9];
  ARP[4] = F[2];
  ARP[5] = V[1]*F[2];
  ARP[6] = V[8]*V[11];
  ARP[7] = V[7]*V[11];
  ARP[8] = V[8]*V[9];
  ARP[9] = V[7]*F[2];
  ARP[10] = F[0]*F[2];
  ARP[11] = V[8]*V[10];
  ARP[12] = V[9]*V[10];
  ARP[13] = V[4]*F[2];
  ARP[14] = V[2]*F[2];
  ARP[15] = V[3]*F[2];
  ARP[16] = V[5]*V[11];
  ARP[17] = V[5]*F[2];
  ARP[18] = V[5]*V[6];
  ARP[19] = V[11];
  ARP[20] = V[7];
  ARP[21] = V[3];
  ARP[22] = V[3];
  ARP[23] = V[4];
  ARP[24] = V[6];
  ARP[25] = V[6];
  ARP[26] = V[1];
  ARP[27] = V[6];
  ARP[28] = F[3];
  ARP[29] = F[3];
  ARP[30] = V[11];
  ARP[31] = V[8];
  ARP[32] = V[7];
  ARP[33] = V[5];
  ARP[34] = V[0];
  ARP[35] = V[4];
  ARP[36] = V[1];
}

/* End of ReactantProd function                                     */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/* JacReactantProd - Jacobian of Reactant Products vector           */
/*   Arguments :                                                    */
/*      V         - Concentrations of variable species (local)      */
/*      F         - Concentrations of fixed species (local)         */
/*      JVRP      - d ARP(1:NREACT)/d VAR (1:NVAR)                  */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void JacReactantProd( 
  double V[],                            /* Concentrations of variable species (local) */
  double F[],                            /* Concentrations of fixed species (local) */
  double JVRP[]                          /* d ARP(1:NREACT)/d VAR (1:NVAR) */
)
{

/* Reactant Products in each equation are useful in the             */
/*    stoichiometric formulation of mass action law                 */
/* Below we compute the Jacobian of the Reactant Products vector    */
/*    w.r.t. variable species: d ARP(1:NREACT) / d Var(1:NVAR)      */

/* JVRP(0) = dARP(0)/dV(11)                                         */
  JVRP[0] = F[2];
/* JVRP(1) = dARP(1)/dV(9)                                          */
  JVRP[1] = V[11];
/* JVRP(2) = dARP(1)/dV(11)                                         */
  JVRP[2] = V[9];
/* JVRP(3) = dARP(2)/dV(9)                                          */
  JVRP[3] = F[2];
/* JVRP(4) = dARP(3)/dV(9)                                          */
  JVRP[4] = 2*V[9];
/* JVRP(5) = dARP(5)/dV(1)                                          */
  JVRP[5] = F[2];
/* JVRP(6) = dARP(6)/dV(8)                                          */
  JVRP[6] = V[11];
/* JVRP(7) = dARP(6)/dV(11)                                         */
  JVRP[7] = V[8];
/* JVRP(8) = dARP(7)/dV(7)                                          */
  JVRP[8] = V[11];
/* JVRP(9) = dARP(7)/dV(11)                                         */
  JVRP[9] = V[7];
/* JVRP(10) = dARP(8)/dV(8)                                         */
  JVRP[10] = V[9];
/* JVRP(11) = dARP(8)/dV(9)                                         */
  JVRP[11] = V[8];
/* JVRP(12) = dARP(9)/dV(7)                                         */
  JVRP[12] = F[2];
/* JVRP(13) = dARP(11)/dV(8)                                        */
  JVRP[13] = V[10];
/* JVRP(14) = dARP(11)/dV(10)                                       */
  JVRP[14] = V[8];
/* JVRP(15) = dARP(12)/dV(9)                                        */
  JVRP[15] = V[10];
/* JVRP(16) = dARP(12)/dV(10)                                       */
  JVRP[16] = V[9];
/* JVRP(17) = dARP(13)/dV(4)                                        */
  JVRP[17] = F[2];
/* JVRP(18) = dARP(14)/dV(2)                                        */
  JVRP[18] = F[2];
/* JVRP(19) = dARP(15)/dV(3)                                        */
  JVRP[19] = F[2];
/* JVRP(20) = dARP(16)/dV(5)                                        */
  JVRP[20] = V[11];
/* JVRP(21) = dARP(16)/dV(11)                                       */
  JVRP[21] = V[5];
/* JVRP(22) = dARP(17)/dV(5)                                        */
  JVRP[22] = F[2];
/* JVRP(23) = dARP(18)/dV(5)                                        */
  JVRP[23] = V[6];
/* JVRP(24) = dARP(18)/dV(6)                                        */
  JVRP[24] = V[5];
/* JVRP(25) = dARP(19)/dV(11)                                       */
  JVRP[25] = 1;
/* JVRP(26) = dARP(20)/dV(7)                                        */
  JVRP[26] = 1;
/* JVRP(27) = dARP(21)/dV(3)                                        */
  JVRP[27] = 1;
/* JVRP(28) = dARP(22)/dV(3)                                        */
  JVRP[28] = 1;
/* JVRP(29) = dARP(23)/dV(4)                                        */
  JVRP[29] = 1;
/* JVRP(30) = dARP(24)/dV(6)                                        */
  JVRP[30] = 1;
/* JVRP(31) = dARP(25)/dV(6)                                        */
  JVRP[31] = 1;
/* JVRP(32) = dARP(26)/dV(1)                                        */
  JVRP[32] = 1;
/* JVRP(33) = dARP(27)/dV(6)                                        */
  JVRP[33] = 1;
/* JVRP(34) = dARP(30)/dV(11)                                       */
  JVRP[34] = 1;
/* JVRP(35) = dARP(31)/dV(8)                                        */
  JVRP[35] = 1;
/* JVRP(36) = dARP(32)/dV(7)                                        */
  JVRP[36] = 1;
/* JVRP(37) = dARP(33)/dV(5)                                        */
  JVRP[37] = 1;
/* JVRP(38) = dARP(34)/dV(0)                                        */
  JVRP[38] = 1;
/* JVRP(39) = dARP(35)/dV(4)                                        */
  JVRP[39] = 1;
/* JVRP(40) = dARP(36)/dV(1)                                        */
  JVRP[40] = 1;
}

/* End of JacReactantProd function                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */



/* Begin Derivative w.r.t. Rate Coefficients                        */


/* End Derivative w.r.t. Rate Coefficients                          */


/* Begin Jacobian Derivative w.r.t. Rate Coefficients               */


/* End Jacobian Derivative w.r.t. Rate Coefficients                 */

