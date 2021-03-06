/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/* Hessian File                                                     */
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
/* File                 : orlando_ohss_Hessian.c                    */
/* Time                 : Sun Aug  9 15:18:54 2020                  */
/* Working directory    : /home/WUR/krol005/kpp/examples            */
/* Equation file        : orlando_ohss.kpp                          */
/* Output root filename : orlando_ohss                              */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "orlando_ohss_Parameters.h"
#include "orlando_ohss_Global.h"
#include "orlando_ohss_Sparse.h"


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/* Hessian - function for Hessian (Jac derivative w.r.t. variables) */
/*   Arguments :                                                    */
/*      V         - Concentrations of variable species (local)      */
/*      F         - Concentrations of fixed species (local)         */
/*      RCT       - Rate constants (local)                          */
/*      HESS      - Hessian of Var (i.e. the 3-tensor d Jac / d Var) */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void Hessian( 
  double V[],                            /* Concentrations of variable species (local) */
  double F[],                            /* Concentrations of fixed species (local) */
  double RCT[],                          /* Rate constants (local) */
  double HESS[]                          /* Hessian of Var (i.e. the 3-tensor d Jac / d Var) */
)
{
/* --------------------------------------------------------         */
/* Note: HESS is represented in coordinate sparse format:           */
/*       HESS(m) = d^2 f_i / dv_j dv_k = d Jac_{i,j} / dv_k         */
/*       where i = IHESS_I(m), j = IHESS_J(m), k = IHESS_K(m).      */
/* --------------------------------------------------------         */
/* Note: d^2 f_i / dv_j dv_k = d^2 f_i / dv_k dv_j,                 */
/*       therefore only the terms d^2 f_i / dv_j dv_k               */
/*       with j <= k are computed and stored in HESS.               */
/* --------------------------------------------------------         */

/* Local variables                                                  */
double D2A[20];                          /* Second derivatives of equation rates */

/* Computation of the second derivatives of equation rates          */
/* D2A(0) = d^2 A(6)/{dV(17)dV(17)}                                 */
  D2A[0] = RCT[6]*2;
/* D2A(1) = d^2 A(7) / dV(12)dV(17)                                 */
  D2A[1] = RCT[7];
/* D2A(2) = d^2 A(8) / dV(12)dV(13)                                 */
  D2A[2] = RCT[8];
/* D2A(3) = d^2 A(10) / dV(11)dV(12)                                */
  D2A[3] = RCT[10];
/* D2A(4) = d^2 A(11) / dV(13)dV(17)                                */
  D2A[4] = RCT[11];
/* D2A(5) = d^2 A(12) / dV(11)dV(16)                                */
  D2A[5] = RCT[12];
/* D2A(6) = d^2 A(13) / dV(13)dV(16)                                */
  D2A[6] = RCT[13];
/* D2A(7) = d^2 A(14) / dV(16)dV(17)                                */
  D2A[7] = RCT[14];
/* D2A(8) = d^2 A(18) / dV(13)dV(14)                                */
  D2A[8] = RCT[18];
/* D2A(9) = d^2 A(19) / dV(14)dV(17)                                */
  D2A[9] = RCT[19];
/* D2A(10) = d^2 A(24) / dV(10)dV(14)                               */
  D2A[10] = RCT[24];
/* D2A(11) = d^2 A(25) / dV(10)dV(17)                               */
  D2A[11] = RCT[25];
/* D2A(12) = d^2 A(26) / dV(10)dV(13)                               */
  D2A[12] = RCT[26];
/* D2A(13) = d^2 A(27) / dV(14)dV(15)                               */
  D2A[13] = RCT[27];
/* D2A(14) = d^2 A(28) / dV(15)dV(17)                               */
  D2A[14] = RCT[28];
/* D2A(15) = d^2 A(30) / dV(13)dV(15)                               */
  D2A[15] = RCT[30];
/* D2A(16) = d^2 A(35) / dV(9)dV(17)                                */
  D2A[16] = RCT[35];
/* D2A(17) = d^2 A(36) / dV(9)dV(13)                                */
  D2A[17] = RCT[36];
/* D2A(18) = d^2 A(38) / dV(7)dV(12)                                */
  D2A[18] = RCT[38];
/* D2A(19) = d^2 A(39) / dV(7)dV(16)                                */
  D2A[19] = RCT[39];

/* Computation of the Jacobian derivative                           */
/* HESS(0) = d^2 Vdot(0)/{dV(7)dV(12)} = d^2 Vdot(0)/{dV(12)dV(7)}  */
  HESS[0] = 0.29*D2A[18];
/* HESS(1) = d^2 Vdot(0)/{dV(17)dV(17)} = d^2 Vdot(0)/{dV(17)dV(17)} */
  HESS[1] = D2A[0];
/* HESS(2) = d^2 Vdot(2)/{dV(11)dV(16)} = d^2 Vdot(2)/{dV(16)dV(11)} */
  HESS[2] = D2A[5];
/* HESS(3) = d^2 Vdot(3)/{dV(13)dV(15)} = d^2 Vdot(3)/{dV(15)dV(13)} */
  HESS[3] = 0.87*D2A[15];
/* HESS(4) = d^2 Vdot(3)/{dV(14)dV(15)} = d^2 Vdot(3)/{dV(15)dV(14)} */
  HESS[4] = 0.75*D2A[13];
/* HESS(5) = d^2 Vdot(4)/{dV(7)dV(12)} = d^2 Vdot(4)/{dV(12)dV(7)}  */
  HESS[5] = 0.363*D2A[18];
/* HESS(6) = d^2 Vdot(4)/{dV(9)dV(13)} = d^2 Vdot(4)/{dV(13)dV(9)}  */
  HESS[6] = 0.25*D2A[17];
/* HESS(7) = d^2 Vdot(4)/{dV(10)dV(13)} = d^2 Vdot(4)/{dV(13)dV(10)} */
  HESS[7] = 0.9*D2A[12];
/* HESS(8) = d^2 Vdot(4)/{dV(13)dV(15)} = d^2 Vdot(4)/{dV(15)dV(13)} */
  HESS[8] = 0.15*D2A[15];
/* HESS(9) = d^2 Vdot(5)/{dV(7)dV(12)} = d^2 Vdot(5)/{dV(12)dV(7)}  */
  HESS[9] = 0.785*D2A[18];
/* HESS(10) = d^2 Vdot(5)/{dV(9)dV(13)} = d^2 Vdot(5)/{dV(13)dV(9)} */
  HESS[10] = 0.25*D2A[17];
/* HESS(11) = d^2 Vdot(5)/{dV(10)dV(13)} = d^2 Vdot(5)/{dV(13)dV(10)} */
  HESS[11] = 0.9*D2A[12];
/* HESS(12) = d^2 Vdot(5)/{dV(10)dV(14)} = d^2 Vdot(5)/{dV(14)dV(10)} */
  HESS[12] = 1.5*D2A[10];
/* HESS(13) = d^2 Vdot(5)/{dV(13)dV(14)} = d^2 Vdot(5)/{dV(14)dV(13)} */
  HESS[13] = D2A[8];
/* HESS(14) = d^2 Vdot(5)/{dV(13)dV(15)} = d^2 Vdot(5)/{dV(15)dV(13)} */
  HESS[14] = 0.1*D2A[15];
/* HESS(15) = d^2 Vdot(5)/{dV(14)dV(15)} = d^2 Vdot(5)/{dV(15)dV(14)} */
  HESS[15] = 0.75*D2A[13];
/* HESS(16) = d^2 Vdot(6)/{dV(10)dV(17)} = d^2 Vdot(6)/{dV(17)dV(10)} */
  HESS[16] = D2A[11];
/* HESS(17) = d^2 Vdot(6)/{dV(15)dV(17)} = d^2 Vdot(6)/{dV(17)dV(15)} */
  HESS[17] = D2A[14];
/* HESS(18) = d^2 Vdot(7)/{dV(7)dV(12)} = d^2 Vdot(7)/{dV(12)dV(7)} */
  HESS[18] = -D2A[18];
/* HESS(19) = d^2 Vdot(7)/{dV(7)dV(16)} = d^2 Vdot(7)/{dV(16)dV(7)} */
  HESS[19] = -D2A[19];
/* HESS(20) = d^2 Vdot(8)/{dV(7)dV(12)} = d^2 Vdot(8)/{dV(12)dV(7)} */
  HESS[20] = 0.675*D2A[18];
/* HESS(21) = d^2 Vdot(8)/{dV(10)dV(13)} = d^2 Vdot(8)/{dV(13)dV(10)} */
  HESS[21] = 0.9*D2A[12];
/* HESS(22) = d^2 Vdot(8)/{dV(10)dV(14)} = d^2 Vdot(8)/{dV(14)dV(10)} */
  HESS[22] = 0.75*D2A[10];
/* HESS(23) = d^2 Vdot(9)/{dV(9)dV(13)} = d^2 Vdot(9)/{dV(13)dV(9)} */
  HESS[23] = -D2A[17];
/* HESS(24) = d^2 Vdot(9)/{dV(9)dV(17)} = d^2 Vdot(9)/{dV(17)dV(9)} */
  HESS[24] = -D2A[16];
/* HESS(25) = d^2 Vdot(10)/{dV(10)dV(13)} = d^2 Vdot(10)/{dV(13)dV(10)} */
  HESS[25] = -D2A[12];
/* HESS(26) = d^2 Vdot(10)/{dV(10)dV(14)} = d^2 Vdot(10)/{dV(14)dV(10)} */
  HESS[26] = -D2A[10];
/* HESS(27) = d^2 Vdot(10)/{dV(10)dV(17)} = d^2 Vdot(10)/{dV(17)dV(10)} */
  HESS[27] = -D2A[11];
/* HESS(28) = d^2 Vdot(11)/{dV(9)dV(13)} = d^2 Vdot(11)/{dV(13)dV(9)} */
  HESS[28] = 0.9*D2A[17];
/* HESS(29) = d^2 Vdot(11)/{dV(10)dV(13)} = d^2 Vdot(11)/{dV(13)dV(10)} */
  HESS[29] = 0.9*D2A[12];
/* HESS(30) = d^2 Vdot(11)/{dV(11)dV(12)} = d^2 Vdot(11)/{dV(12)dV(11)} */
  HESS[30] = -D2A[3];
/* HESS(31) = d^2 Vdot(11)/{dV(11)dV(16)} = d^2 Vdot(11)/{dV(16)dV(11)} */
  HESS[31] = -D2A[5];
/* HESS(32) = d^2 Vdot(11)/{dV(12)dV(13)} = d^2 Vdot(11)/{dV(13)dV(12)} */
  HESS[32] = D2A[2];
/* HESS(33) = d^2 Vdot(11)/{dV(13)dV(14)} = d^2 Vdot(11)/{dV(14)dV(13)} */
  HESS[33] = 0.99*D2A[8];
/* HESS(34) = d^2 Vdot(11)/{dV(13)dV(15)} = d^2 Vdot(11)/{dV(15)dV(13)} */
  HESS[34] = 0.9*D2A[15];
/* HESS(35) = d^2 Vdot(11)/{dV(13)dV(16)} = d^2 Vdot(11)/{dV(16)dV(13)} */
  HESS[35] = 2*D2A[6];
/* HESS(36) = d^2 Vdot(11)/{dV(13)dV(17)} = d^2 Vdot(11)/{dV(17)dV(13)} */
  HESS[36] = D2A[4];
/* HESS(37) = d^2 Vdot(11)/{dV(16)dV(17)} = d^2 Vdot(11)/{dV(17)dV(16)} */
  HESS[37] = D2A[7];
/* HESS(38) = d^2 Vdot(12)/{dV(7)dV(12)} = d^2 Vdot(12)/{dV(12)dV(7)} */
  HESS[38] = -D2A[18];
/* HESS(39) = d^2 Vdot(12)/{dV(11)dV(12)} = d^2 Vdot(12)/{dV(12)dV(11)} */
  HESS[39] = -D2A[3];
/* HESS(40) = d^2 Vdot(12)/{dV(12)dV(13)} = d^2 Vdot(12)/{dV(13)dV(12)} */
  HESS[40] = -D2A[2];
/* HESS(41) = d^2 Vdot(12)/{dV(12)dV(17)} = d^2 Vdot(12)/{dV(17)dV(12)} */
  HESS[41] = -D2A[1];
/* HESS(42) = d^2 Vdot(13)/{dV(9)dV(13)} = d^2 Vdot(13)/{dV(13)dV(9)} */
  HESS[42] = -D2A[17];
/* HESS(43) = d^2 Vdot(13)/{dV(10)dV(13)} = d^2 Vdot(13)/{dV(13)dV(10)} */
  HESS[43] = -D2A[12];
/* HESS(44) = d^2 Vdot(13)/{dV(12)dV(13)} = d^2 Vdot(13)/{dV(13)dV(12)} */
  HESS[44] = -D2A[2];
/* HESS(45) = d^2 Vdot(13)/{dV(13)dV(14)} = d^2 Vdot(13)/{dV(14)dV(13)} */
  HESS[45] = -D2A[8];
/* HESS(46) = d^2 Vdot(13)/{dV(13)dV(15)} = d^2 Vdot(13)/{dV(15)dV(13)} */
  HESS[46] = -D2A[15];
/* HESS(47) = d^2 Vdot(13)/{dV(13)dV(16)} = d^2 Vdot(13)/{dV(16)dV(13)} */
  HESS[47] = -D2A[6];
/* HESS(48) = d^2 Vdot(13)/{dV(13)dV(17)} = d^2 Vdot(13)/{dV(17)dV(13)} */
  HESS[48] = -D2A[4];
/* HESS(49) = d^2 Vdot(14)/{dV(7)dV(12)} = d^2 Vdot(14)/{dV(12)dV(7)} */
  HESS[49] = 0.181*D2A[18];
/* HESS(50) = d^2 Vdot(14)/{dV(10)dV(14)} = d^2 Vdot(14)/{dV(14)dV(10)} */
  HESS[50] = -D2A[10];
/* HESS(51) = d^2 Vdot(14)/{dV(13)dV(14)} = d^2 Vdot(14)/{dV(14)dV(13)} */
  HESS[51] = -D2A[8];
/* HESS(52) = d^2 Vdot(14)/{dV(14)dV(15)} = d^2 Vdot(14)/{dV(15)dV(14)} */
  HESS[52] = -D2A[13];
/* HESS(53) = d^2 Vdot(14)/{dV(14)dV(17)} = d^2 Vdot(14)/{dV(17)dV(14)} */
  HESS[53] = -D2A[9];
/* HESS(54) = d^2 Vdot(15)/{dV(13)dV(15)} = d^2 Vdot(15)/{dV(15)dV(13)} */
  HESS[54] = -D2A[15];
/* HESS(55) = d^2 Vdot(15)/{dV(14)dV(15)} = d^2 Vdot(15)/{dV(15)dV(14)} */
  HESS[55] = -D2A[13];
/* HESS(56) = d^2 Vdot(15)/{dV(15)dV(17)} = d^2 Vdot(15)/{dV(17)dV(15)} */
  HESS[56] = -D2A[14];
/* HESS(57) = d^2 Vdot(16)/{dV(7)dV(16)} = d^2 Vdot(16)/{dV(16)dV(7)} */
  HESS[57] = -D2A[19];
/* HESS(58) = d^2 Vdot(16)/{dV(11)dV(12)} = d^2 Vdot(16)/{dV(12)dV(11)} */
  HESS[58] = D2A[3];
/* HESS(59) = d^2 Vdot(16)/{dV(11)dV(16)} = d^2 Vdot(16)/{dV(16)dV(11)} */
  HESS[59] = -D2A[5];
/* HESS(60) = d^2 Vdot(16)/{dV(13)dV(16)} = d^2 Vdot(16)/{dV(16)dV(13)} */
  HESS[60] = -D2A[6];
/* HESS(61) = d^2 Vdot(16)/{dV(16)dV(17)} = d^2 Vdot(16)/{dV(17)dV(16)} */
  HESS[61] = -D2A[7];
/* HESS(62) = d^2 Vdot(17)/{dV(7)dV(12)} = d^2 Vdot(17)/{dV(12)dV(7)} */
  HESS[62] = 0.125*D2A[18];
/* HESS(63) = d^2 Vdot(17)/{dV(9)dV(13)} = d^2 Vdot(17)/{dV(13)dV(9)} */
  HESS[63] = D2A[17];
/* HESS(64) = d^2 Vdot(17)/{dV(9)dV(17)} = d^2 Vdot(17)/{dV(17)dV(9)} */
  HESS[64] = -D2A[16];
/* HESS(65) = d^2 Vdot(17)/{dV(10)dV(13)} = d^2 Vdot(17)/{dV(13)dV(10)} */
  HESS[65] = 0.9*D2A[12];
/* HESS(66) = d^2 Vdot(17)/{dV(10)dV(14)} = d^2 Vdot(17)/{dV(14)dV(10)} */
  HESS[66] = D2A[10];
/* HESS(67) = d^2 Vdot(17)/{dV(10)dV(17)} = d^2 Vdot(17)/{dV(17)dV(10)} */
  HESS[67] = -D2A[11];
/* HESS(68) = d^2 Vdot(17)/{dV(12)dV(17)} = d^2 Vdot(17)/{dV(17)dV(12)} */
  HESS[68] = -D2A[1];
/* HESS(69) = d^2 Vdot(17)/{dV(13)dV(14)} = d^2 Vdot(17)/{dV(14)dV(13)} */
  HESS[69] = D2A[8];
/* HESS(70) = d^2 Vdot(17)/{dV(13)dV(15)} = d^2 Vdot(17)/{dV(15)dV(13)} */
  HESS[70] = 0.97*D2A[15];
/* HESS(71) = d^2 Vdot(17)/{dV(13)dV(17)} = d^2 Vdot(17)/{dV(17)dV(13)} */
  HESS[71] = -D2A[4];
/* HESS(72) = d^2 Vdot(17)/{dV(14)dV(15)} = d^2 Vdot(17)/{dV(15)dV(14)} */
  HESS[72] = D2A[13];
/* HESS(73) = d^2 Vdot(17)/{dV(14)dV(17)} = d^2 Vdot(17)/{dV(17)dV(14)} */
  HESS[73] = -D2A[9];
/* HESS(74) = d^2 Vdot(17)/{dV(15)dV(17)} = d^2 Vdot(17)/{dV(17)dV(15)} */
  HESS[74] = -D2A[14];
/* HESS(75) = d^2 Vdot(17)/{dV(16)dV(17)} = d^2 Vdot(17)/{dV(17)dV(16)} */
  HESS[75] = -D2A[7];
/* HESS(76) = d^2 Vdot(17)/{dV(17)dV(17)} = d^2 Vdot(17)/{dV(17)dV(17)} */
  HESS[76] = -2*D2A[0];
}

/* End of Hessian function                                          */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/* HessTR_Vec - Hessian transposed times user vectors               */
/*   Arguments :                                                    */
/*      HESS      - Hessian of Var (i.e. the 3-tensor d Jac / d Var) */
/*      U1        - User vector                                     */
/*      U2        - User vector                                     */
/*      HTU       - Transposed Hessian times user vectors: (Hess x U2)^T * U1 = [d (Jac^T*U1)/d Var] * U2  */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void HessTR_Vec( 
  double HESS[],                         /* Hessian of Var (i.e. the 3-tensor d Jac / d Var) */
  double U1[],                           /* User vector */
  double U2[],                           /* User vector */
  double HTU[]                           /* Transposed Hessian times user vectors: (Hess x U2)^T * U1 = [d (Jac^T*U1)/d Var] * U2  */
)
{
/* Compute the vector HTU =(Hess x U2)^T * U1 = d (Jac^T*U1)/d Var * U2  */
  HTU[0] = 0;
  HTU[1] = 0;
  HTU[2] = 0;
  HTU[3] = 0;
  HTU[4] = 0;
  HTU[5] = 0;
  HTU[6] = 0;
  HTU[7] = HESS[0]*(U1[0]*U2[12])+HESS[5]*(U1[4]*U2[12])+HESS[9]*(U1[5]
          *U2[12])+HESS[18]*(U1[7]*U2[12])+HESS[19]*(U1[7]*U2[16])
          +HESS[20]*(U1[8]*U2[12])+HESS[38]*(U1[12]*U2[12])+HESS[49]
          *(U1[14]*U2[12])+HESS[57]*(U1[16]*U2[16])+HESS[62]*(U1[17]
          *U2[12]);
  HTU[8] = 0;
  HTU[9] = HESS[6]*(U1[4]*U2[13])+HESS[10]*(U1[5]*U2[13])+HESS[23]
          *(U1[9]*U2[13])+HESS[24]*(U1[9]*U2[17])+HESS[28]*(U1[11]
          *U2[13])+HESS[42]*(U1[13]*U2[13])+HESS[63]*(U1[17]*U2[13])
          +HESS[64]*(U1[17]*U2[17]);
  HTU[10] = HESS[7]*(U1[4]*U2[13])+HESS[11]*(U1[5]*U2[13])+HESS[12]
           *(U1[5]*U2[14])+HESS[16]*(U1[6]*U2[17])+HESS[21]*(U1[8]
           *U2[13])+HESS[22]*(U1[8]*U2[14])+HESS[25]*(U1[10]*U2[13])
           +HESS[26]*(U1[10]*U2[14])+HESS[27]*(U1[10]*U2[17])+HESS[29]
           *(U1[11]*U2[13])+HESS[43]*(U1[13]*U2[13])+HESS[50]*(U1[14]
           *U2[14])+HESS[65]*(U1[17]*U2[13])+HESS[66]*(U1[17]*U2[14])
           +HESS[67]*(U1[17]*U2[17]);
  HTU[11] = HESS[2]*(U1[2]*U2[16])+HESS[30]*(U1[11]*U2[12])+HESS[31]
           *(U1[11]*U2[16])+HESS[39]*(U1[12]*U2[12])+HESS[58]*(U1[16]
           *U2[12])+HESS[59]*(U1[16]*U2[16]);
  HTU[12] = HESS[0]*(U1[0]*U2[7])+HESS[5]*(U1[4]*U2[7])+HESS[9]*(U1[5]
           *U2[7])+HESS[18]*(U1[7]*U2[7])+HESS[20]*(U1[8]*U2[7])
           +HESS[30]*(U1[11]*U2[11])+HESS[32]*(U1[11]*U2[13])+HESS[38]
           *(U1[12]*U2[7])+HESS[39]*(U1[12]*U2[11])+HESS[40]*(U1[12]
           *U2[13])+HESS[41]*(U1[12]*U2[17])+HESS[44]*(U1[13]*U2[13])
           +HESS[49]*(U1[14]*U2[7])+HESS[58]*(U1[16]*U2[11])+HESS[62]
           *(U1[17]*U2[7])+HESS[68]*(U1[17]*U2[17]);
  HTU[13] = HESS[3]*(U1[3]*U2[15])+HESS[6]*(U1[4]*U2[9])+HESS[7]*(U1[4]
           *U2[10])+HESS[8]*(U1[4]*U2[15])+HESS[10]*(U1[5]*U2[9])
           +HESS[11]*(U1[5]*U2[10])+HESS[13]*(U1[5]*U2[14])+HESS[14]
           *(U1[5]*U2[15])+HESS[21]*(U1[8]*U2[10])+HESS[23]*(U1[9]
           *U2[9])+HESS[25]*(U1[10]*U2[10])+HESS[28]*(U1[11]*U2[9])
           +HESS[29]*(U1[11]*U2[10])+HESS[32]*(U1[11]*U2[12])+HESS[33]
           *(U1[11]*U2[14])+HESS[34]*(U1[11]*U2[15])+HESS[35]*(U1[11]
           *U2[16])+HESS[36]*(U1[11]*U2[17])+HESS[40]*(U1[12]*U2[12])
           +HESS[42]*(U1[13]*U2[9])+HESS[43]*(U1[13]*U2[10])+HESS[44]
           *(U1[13]*U2[12])+HESS[45]*(U1[13]*U2[14])+HESS[46]*(U1[13]
           *U2[15])+HESS[47]*(U1[13]*U2[16])+HESS[48]*(U1[13]*U2[17])
           +HESS[51]*(U1[14]*U2[14])+HESS[54]*(U1[15]*U2[15])+HESS[60]
           *(U1[16]*U2[16])+HESS[63]*(U1[17]*U2[9])+HESS[65]*(U1[17]
           *U2[10])+HESS[69]*(U1[17]*U2[14])+HESS[70]*(U1[17]*U2[15])
           +HESS[71]*(U1[17]*U2[17]);
  HTU[14] = HESS[4]*(U1[3]*U2[15])+HESS[12]*(U1[5]*U2[10])+HESS[13]
           *(U1[5]*U2[13])+HESS[15]*(U1[5]*U2[15])+HESS[22]*(U1[8]
           *U2[10])+HESS[26]*(U1[10]*U2[10])+HESS[33]*(U1[11]*U2[13])
           +HESS[45]*(U1[13]*U2[13])+HESS[50]*(U1[14]*U2[10])+HESS[51]
           *(U1[14]*U2[13])+HESS[52]*(U1[14]*U2[15])+HESS[53]*(U1[14]
           *U2[17])+HESS[55]*(U1[15]*U2[15])+HESS[66]*(U1[17]*U2[10])
           +HESS[69]*(U1[17]*U2[13])+HESS[72]*(U1[17]*U2[15])+HESS[73]
           *(U1[17]*U2[17]);
  HTU[15] = HESS[3]*(U1[3]*U2[13])+HESS[4]*(U1[3]*U2[14])+HESS[8]
           *(U1[4]*U2[13])+HESS[14]*(U1[5]*U2[13])+HESS[15]*(U1[5]
           *U2[14])+HESS[17]*(U1[6]*U2[17])+HESS[34]*(U1[11]*U2[13])
           +HESS[46]*(U1[13]*U2[13])+HESS[52]*(U1[14]*U2[14])+HESS[54]
           *(U1[15]*U2[13])+HESS[55]*(U1[15]*U2[14])+HESS[56]*(U1[15]
           *U2[17])+HESS[70]*(U1[17]*U2[13])+HESS[72]*(U1[17]*U2[14])
           +HESS[74]*(U1[17]*U2[17]);
  HTU[16] = HESS[2]*(U1[2]*U2[11])+HESS[19]*(U1[7]*U2[7])+HESS[31]
           *(U1[11]*U2[11])+HESS[35]*(U1[11]*U2[13])+HESS[37]*(U1[11]
           *U2[17])+HESS[47]*(U1[13]*U2[13])+HESS[57]*(U1[16]*U2[7])
           +HESS[59]*(U1[16]*U2[11])+HESS[60]*(U1[16]*U2[13])+HESS[61]
           *(U1[16]*U2[17])+HESS[75]*(U1[17]*U2[17]);
  HTU[17] = HESS[1]*(U1[0]*U2[17])+HESS[16]*(U1[6]*U2[10])+HESS[17]
           *(U1[6]*U2[15])+HESS[24]*(U1[9]*U2[9])+HESS[27]*(U1[10]
           *U2[10])+HESS[36]*(U1[11]*U2[13])+HESS[37]*(U1[11]*U2[16])
           +HESS[41]*(U1[12]*U2[12])+HESS[48]*(U1[13]*U2[13])+HESS[53]
           *(U1[14]*U2[14])+HESS[56]*(U1[15]*U2[15])+HESS[61]*(U1[16]
           *U2[16])+HESS[64]*(U1[17]*U2[9])+HESS[67]*(U1[17]*U2[10])
           +HESS[68]*(U1[17]*U2[12])+HESS[71]*(U1[17]*U2[13])+HESS[73]
           *(U1[17]*U2[14])+HESS[74]*(U1[17]*U2[15])+HESS[75]*(U1[17]
           *U2[16])+HESS[76]*(U1[17]*U2[17]);
}

/* End of HessTR_Vec function                                       */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/* Hess_Vec - Hessian times user vectors                            */
/*   Arguments :                                                    */
/*      HESS      - Hessian of Var (i.e. the 3-tensor d Jac / d Var) */
/*      U1        - User vector                                     */
/*      U2        - User vector                                     */
/*      HU        - Hessian times user vectors: (Hess x U2) * U1 = [d (Jac*U1)/d Var] * U2 */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void Hess_Vec( 
  double HESS[],                         /* Hessian of Var (i.e. the 3-tensor d Jac / d Var) */
  double U1[],                           /* User vector */
  double U2[],                           /* User vector */
  double HU[]                            /* Hessian times user vectors: (Hess x U2) * U1 = [d (Jac*U1)/d Var] * U2 */
)
{
/* Compute the vector HU =(Hess x U2) * U1 = d (Jac*U1)/d Var * U2  */
  HU[0] = HESS[0]*(U1[7]*U2[12])+HESS[0]*(U1[12]*U2[7])+HESS[1]*(U1[17]
         *U2[17]);
  HU[1] = 0;
  HU[2] = HESS[2]*(U1[11]*U2[16])+HESS[2]*(U1[16]*U2[11]);
  HU[3] = HESS[3]*(U1[13]*U2[15])+HESS[3]*(U1[15]*U2[13])+HESS[4]
         *(U1[14]*U2[15])+HESS[4]*(U1[15]*U2[14]);
  HU[4] = HESS[5]*(U1[7]*U2[12])+HESS[5]*(U1[12]*U2[7])+HESS[6]*(U1[9]
         *U2[13])+HESS[6]*(U1[13]*U2[9])+HESS[7]*(U1[10]*U2[13])
         +HESS[7]*(U1[13]*U2[10])+HESS[8]*(U1[13]*U2[15])+HESS[8]
         *(U1[15]*U2[13]);
  HU[5] = HESS[9]*(U1[7]*U2[12])+HESS[9]*(U1[12]*U2[7])+HESS[10]*(U1[9]
         *U2[13])+HESS[10]*(U1[13]*U2[9])+HESS[11]*(U1[10]*U2[13])
         +HESS[11]*(U1[13]*U2[10])+HESS[12]*(U1[10]*U2[14])+HESS[12]
         *(U1[14]*U2[10])+HESS[13]*(U1[13]*U2[14])+HESS[13]*(U1[14]
         *U2[13])+HESS[14]*(U1[13]*U2[15])+HESS[14]*(U1[15]*U2[13])
         +HESS[15]*(U1[14]*U2[15])+HESS[15]*(U1[15]*U2[14]);
  HU[6] = HESS[16]*(U1[10]*U2[17])+HESS[16]*(U1[17]*U2[10])+HESS[17]
         *(U1[15]*U2[17])+HESS[17]*(U1[17]*U2[15]);
  HU[7] = HESS[18]*(U1[7]*U2[12])+HESS[18]*(U1[12]*U2[7])+HESS[19]
         *(U1[7]*U2[16])+HESS[19]*(U1[16]*U2[7]);
  HU[8] = HESS[20]*(U1[7]*U2[12])+HESS[20]*(U1[12]*U2[7])+HESS[21]
         *(U1[10]*U2[13])+HESS[21]*(U1[13]*U2[10])+HESS[22]*(U1[10]
         *U2[14])+HESS[22]*(U1[14]*U2[10]);
  HU[9] = HESS[23]*(U1[9]*U2[13])+HESS[23]*(U1[13]*U2[9])+HESS[24]
         *(U1[9]*U2[17])+HESS[24]*(U1[17]*U2[9]);
  HU[10] = HESS[25]*(U1[10]*U2[13])+HESS[25]*(U1[13]*U2[10])+HESS[26]
          *(U1[10]*U2[14])+HESS[26]*(U1[14]*U2[10])+HESS[27]*(U1[10]
          *U2[17])+HESS[27]*(U1[17]*U2[10]);
  HU[11] = HESS[28]*(U1[9]*U2[13])+HESS[28]*(U1[13]*U2[9])+HESS[29]
          *(U1[10]*U2[13])+HESS[29]*(U1[13]*U2[10])+HESS[30]*(U1[11]
          *U2[12])+HESS[30]*(U1[12]*U2[11])+HESS[31]*(U1[11]*U2[16])
          +HESS[31]*(U1[16]*U2[11])+HESS[32]*(U1[12]*U2[13])+HESS[32]
          *(U1[13]*U2[12])+HESS[33]*(U1[13]*U2[14])+HESS[33]*(U1[14]
          *U2[13])+HESS[34]*(U1[13]*U2[15])+HESS[34]*(U1[15]*U2[13])
          +HESS[35]*(U1[13]*U2[16])+HESS[35]*(U1[16]*U2[13])+HESS[36]
          *(U1[13]*U2[17])+HESS[36]*(U1[17]*U2[13])+HESS[37]*(U1[16]
          *U2[17])+HESS[37]*(U1[17]*U2[16]);
  HU[12] = HESS[38]*(U1[7]*U2[12])+HESS[38]*(U1[12]*U2[7])+HESS[39]
          *(U1[11]*U2[12])+HESS[39]*(U1[12]*U2[11])+HESS[40]*(U1[12]
          *U2[13])+HESS[40]*(U1[13]*U2[12])+HESS[41]*(U1[12]*U2[17])
          +HESS[41]*(U1[17]*U2[12]);
  HU[13] = HESS[42]*(U1[9]*U2[13])+HESS[42]*(U1[13]*U2[9])+HESS[43]
          *(U1[10]*U2[13])+HESS[43]*(U1[13]*U2[10])+HESS[44]*(U1[12]
          *U2[13])+HESS[44]*(U1[13]*U2[12])+HESS[45]*(U1[13]*U2[14])
          +HESS[45]*(U1[14]*U2[13])+HESS[46]*(U1[13]*U2[15])+HESS[46]
          *(U1[15]*U2[13])+HESS[47]*(U1[13]*U2[16])+HESS[47]*(U1[16]
          *U2[13])+HESS[48]*(U1[13]*U2[17])+HESS[48]*(U1[17]*U2[13]);
  HU[14] = HESS[49]*(U1[7]*U2[12])+HESS[49]*(U1[12]*U2[7])+HESS[50]
          *(U1[10]*U2[14])+HESS[50]*(U1[14]*U2[10])+HESS[51]*(U1[13]
          *U2[14])+HESS[51]*(U1[14]*U2[13])+HESS[52]*(U1[14]*U2[15])
          +HESS[52]*(U1[15]*U2[14])+HESS[53]*(U1[14]*U2[17])+HESS[53]
          *(U1[17]*U2[14]);
  HU[15] = HESS[54]*(U1[13]*U2[15])+HESS[54]*(U1[15]*U2[13])+HESS[55]
          *(U1[14]*U2[15])+HESS[55]*(U1[15]*U2[14])+HESS[56]*(U1[15]
          *U2[17])+HESS[56]*(U1[17]*U2[15]);
  HU[16] = HESS[57]*(U1[7]*U2[16])+HESS[57]*(U1[16]*U2[7])+HESS[58]
          *(U1[11]*U2[12])+HESS[58]*(U1[12]*U2[11])+HESS[59]*(U1[11]
          *U2[16])+HESS[59]*(U1[16]*U2[11])+HESS[60]*(U1[13]*U2[16])
          +HESS[60]*(U1[16]*U2[13])+HESS[61]*(U1[16]*U2[17])+HESS[61]
          *(U1[17]*U2[16]);
  HU[17] = HESS[62]*(U1[7]*U2[12])+HESS[62]*(U1[12]*U2[7])+HESS[63]
          *(U1[9]*U2[13])+HESS[63]*(U1[13]*U2[9])+HESS[64]*(U1[9]
          *U2[17])+HESS[64]*(U1[17]*U2[9])+HESS[65]*(U1[10]*U2[13])
          +HESS[65]*(U1[13]*U2[10])+HESS[66]*(U1[10]*U2[14])+HESS[66]
          *(U1[14]*U2[10])+HESS[67]*(U1[10]*U2[17])+HESS[67]*(U1[17]
          *U2[10])+HESS[68]*(U1[12]*U2[17])+HESS[68]*(U1[17]*U2[12])
          +HESS[69]*(U1[13]*U2[14])+HESS[69]*(U1[14]*U2[13])+HESS[70]
          *(U1[13]*U2[15])+HESS[70]*(U1[15]*U2[13])+HESS[71]*(U1[13]
          *U2[17])+HESS[71]*(U1[17]*U2[13])+HESS[72]*(U1[14]*U2[15])
          +HESS[72]*(U1[15]*U2[14])+HESS[73]*(U1[14]*U2[17])+HESS[73]
          *(U1[17]*U2[14])+HESS[74]*(U1[15]*U2[17])+HESS[74]*(U1[17]
          *U2[15])+HESS[75]*(U1[16]*U2[17])+HESS[75]*(U1[17]*U2[16])
          +HESS[76]*(U1[17]*U2[17]);
}

/* End of Hess_Vec function                                         */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


