/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/* The ODE Jacobian of Chemical Model File                          */
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
/* File                 : orlando_Jacobian.c                        */
/* Time                 : Wed Aug 19 11:17:37 2020                  */
/* Working directory    : /home/WUR/krol005/kpp/examples            */
/* Equation file        : orlando.kpp                               */
/* Output root filename : orlando                                   */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "orlando_Parameters.h"
#include "orlando_Global.h"
#include "orlando_Sparse.h"


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/* Jac_SP - the Jacobian of Variables in sparse matrix representation */
/*   Arguments :                                                    */
/*      V         - Concentrations of variable species (local)      */
/*      F         - Concentrations of fixed species (local)         */
/*      RCT       - Rate constants (local)                          */
/*      JVS       - sparse Jacobian of variables                    */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void Jac_SP( 
  double V[],                            /* Concentrations of variable species (local) */
  double F[],                            /* Concentrations of fixed species (local) */
  double RCT[],                          /* Rate constants (local) */
  double JVS[]                           /* sparse Jacobian of variables */
)
{

/* Local variables                                                  */
double B[80];                            /* Temporary array */

/* B(0) = dA(0)/dV(17)                                              */
  B[0] = RCT[0];
/* B(1) = dA(1)/dV(11)                                              */
  B[1] = RCT[1];
/* B(2) = dA(2)/dV(5)                                               */
  B[2] = RCT[2];
/* B(3) = dA(3)/dV(5)                                               */
  B[3] = RCT[3];
/* B(4) = dA(4)/dV(12)                                              */
  B[4] = RCT[4]*V[18];
/* B(5) = dA(4)/dV(18)                                              */
  B[5] = RCT[4]*V[12];
/* B(6) = dA(5)/dV(17)                                              */
  B[6] = RCT[5]*V[18];
/* B(7) = dA(5)/dV(18)                                              */
  B[7] = RCT[5]*V[17];
/* B(8) = dA(6)/dV(12)                                              */
  B[8] = RCT[6]*2*V[12];
/* B(9) = dA(7)/dV(12)                                              */
  B[9] = RCT[7]*V[17];
/* B(10) = dA(7)/dV(17)                                             */
  B[10] = RCT[7]*V[12];
/* B(11) = dA(8)/dV(13)                                             */
  B[11] = RCT[8]*V[17];
/* B(12) = dA(8)/dV(17)                                             */
  B[12] = RCT[8]*V[13];
/* B(13) = dA(9)/dV(11)                                             */
  B[13] = RCT[9]*V[18];
/* B(14) = dA(9)/dV(18)                                             */
  B[14] = RCT[9]*V[11];
/* B(15) = dA(10)/dV(11)                                            */
  B[15] = RCT[10]*V[17];
/* B(16) = dA(10)/dV(17)                                            */
  B[16] = RCT[10]*V[11];
/* B(17) = dA(11)/dV(12)                                            */
  B[17] = RCT[11]*V[13];
/* B(18) = dA(11)/dV(13)                                            */
  B[18] = RCT[11]*V[12];
/* B(19) = dA(12)/dV(11)                                            */
  B[19] = RCT[12]*V[16];
/* B(20) = dA(12)/dV(16)                                            */
  B[20] = RCT[12]*V[11];
/* B(21) = dA(13)/dV(13)                                            */
  B[21] = RCT[13]*V[16];
/* B(22) = dA(13)/dV(16)                                            */
  B[22] = RCT[13]*V[13];
/* B(23) = dA(14)/dV(12)                                            */
  B[23] = RCT[14]*V[16];
/* B(24) = dA(14)/dV(16)                                            */
  B[24] = RCT[14]*V[12];
/* B(25) = dA(15)/dV(2)                                             */
  B[25] = RCT[15];
/* B(26) = dA(16)/dV(2)                                             */
  B[26] = RCT[16];
/* B(27) = dA(17)/dV(16)                                            */
  B[27] = RCT[17];
/* B(28) = dA(18)/dV(13)                                            */
  B[28] = RCT[18]*V[15];
/* B(29) = dA(18)/dV(15)                                            */
  B[29] = RCT[18]*V[13];
/* B(30) = dA(19)/dV(12)                                            */
  B[30] = RCT[19]*V[15];
/* B(31) = dA(19)/dV(15)                                            */
  B[31] = RCT[19]*V[12];
/* B(32) = dA(20)/dV(4)                                             */
  B[32] = RCT[20]*V[18];
/* B(33) = dA(20)/dV(18)                                            */
  B[33] = RCT[20]*V[4];
/* B(34) = dA(21)/dV(0)                                             */
  B[34] = RCT[21]*V[18];
/* B(35) = dA(21)/dV(18)                                            */
  B[35] = RCT[21]*V[0];
/* B(36) = dA(22)/dV(5)                                             */
  B[36] = RCT[22]*V[18];
/* B(37) = dA(22)/dV(18)                                            */
  B[37] = RCT[22]*V[5];
/* B(38) = dA(23)/dV(7)                                             */
  B[38] = RCT[23]*V[18];
/* B(39) = dA(23)/dV(18)                                            */
  B[39] = RCT[23]*V[7];
/* B(40) = dA(24)/dV(10)                                            */
  B[40] = RCT[24]*V[15];
/* B(41) = dA(24)/dV(15)                                            */
  B[41] = RCT[24]*V[10];
/* B(42) = dA(25)/dV(10)                                            */
  B[42] = RCT[25]*V[12];
/* B(43) = dA(25)/dV(12)                                            */
  B[43] = RCT[25]*V[10];
/* B(44) = dA(26)/dV(10)                                            */
  B[44] = RCT[26]*V[13];
/* B(45) = dA(26)/dV(13)                                            */
  B[45] = RCT[26]*V[10];
/* B(46) = dA(27)/dV(14)                                            */
  B[46] = RCT[27]*V[15];
/* B(47) = dA(27)/dV(15)                                            */
  B[47] = RCT[27]*V[14];
/* B(48) = dA(28)/dV(12)                                            */
  B[48] = RCT[28]*V[14];
/* B(49) = dA(28)/dV(14)                                            */
  B[49] = RCT[28]*V[12];
/* B(50) = dA(29)/dV(14)                                            */
  B[50] = RCT[29];
/* B(51) = dA(30)/dV(13)                                            */
  B[51] = RCT[30]*V[14];
/* B(52) = dA(30)/dV(14)                                            */
  B[52] = RCT[30]*V[13];
/* B(53) = dA(31)/dV(6)                                             */
  B[53] = RCT[31]*V[18];
/* B(54) = dA(31)/dV(18)                                            */
  B[54] = RCT[31]*V[6];
/* B(55) = dA(32)/dV(6)                                             */
  B[55] = RCT[32];
/* B(56) = dA(33)/dV(3)                                             */
  B[56] = RCT[33]*V[18];
/* B(57) = dA(33)/dV(18)                                            */
  B[57] = RCT[33]*V[3];
/* B(58) = dA(34)/dV(3)                                             */
  B[58] = RCT[34];
/* B(59) = dA(35)/dV(8)                                             */
  B[59] = RCT[35]*V[12];
/* B(60) = dA(35)/dV(12)                                            */
  B[60] = RCT[35]*V[8];
/* B(61) = dA(36)/dV(8)                                             */
  B[61] = RCT[36]*V[13];
/* B(62) = dA(36)/dV(13)                                            */
  B[62] = RCT[36]*V[8];
/* B(63) = dA(37)/dV(9)                                             */
  B[63] = RCT[37]*V[18];
/* B(64) = dA(37)/dV(18)                                            */
  B[64] = RCT[37]*V[9];
/* B(65) = dA(38)/dV(7)                                             */
  B[65] = RCT[38]*V[17];
/* B(66) = dA(38)/dV(17)                                            */
  B[66] = RCT[38]*V[7];
/* B(67) = dA(39)/dV(7)                                             */
  B[67] = RCT[39]*V[16];
/* B(68) = dA(39)/dV(16)                                            */
  B[68] = RCT[39]*V[7];
/* B(69) = dA(40)/dV(1)                                             */
  B[69] = RCT[40];
/* B(70) = dA(41)/dV(1)                                             */
  B[70] = RCT[41];
/* B(71) = dA(42)/dV(17)                                            */
  B[71] = RCT[42];
/* B(72) = dA(43)/dV(13)                                            */
  B[72] = RCT[43];
/* B(73) = dA(44)/dV(11)                                            */
  B[73] = RCT[44];
/* B(74) = dA(45)/dV(5)                                             */
  B[74] = RCT[45];
/* B(75) = dA(46)/dV(6)                                             */
  B[75] = RCT[46];
/* B(76) = dA(47)/dV(3)                                             */
  B[76] = RCT[47];

/* Construct the Jacobian terms from B's                            */
/* JVS(0) = Jac_FULL(0,0)                                           */
  JVS[0] = -B[34];
/* JVS(1) = Jac_FULL(0,18)                                          */
  JVS[1] = -B[35];
/* JVS(2) = Jac_FULL(1,1)                                           */
  JVS[2] = -B[69]-B[70];
/* JVS(3) = Jac_FULL(1,7)                                           */
  JVS[3] = 0.29*B[65];
/* JVS(4) = Jac_FULL(1,12)                                          */
  JVS[4] = B[8];
/* JVS(5) = Jac_FULL(1,17)                                          */
  JVS[5] = 0.29*B[66];
/* JVS(6) = Jac_FULL(2,2)                                           */
  JVS[6] = -B[25]-B[26];
/* JVS(7) = Jac_FULL(2,11)                                          */
  JVS[7] = B[19];
/* JVS(8) = Jac_FULL(2,16)                                          */
  JVS[8] = B[20];
/* JVS(9) = Jac_FULL(3,3)                                           */
  JVS[9] = -B[56]-B[58]-B[76];
/* JVS(10) = Jac_FULL(3,13)                                         */
  JVS[10] = 0.87*B[51];
/* JVS(11) = Jac_FULL(3,14)                                         */
  JVS[11] = 0.75*B[46]+B[50]+0.87*B[52];
/* JVS(12) = Jac_FULL(3,15)                                         */
  JVS[12] = 0.75*B[47];
/* JVS(13) = Jac_FULL(3,18)                                         */
  JVS[13] = -B[57];
/* JVS(14) = Jac_FULL(4,4)                                          */
  JVS[14] = -B[32];
/* JVS(15) = Jac_FULL(4,5)                                          */
  JVS[15] = B[2]+B[3]+B[36];
/* JVS(16) = Jac_FULL(4,7)                                          */
  JVS[16] = 0.363*B[65];
/* JVS(17) = Jac_FULL(4,8)                                          */
  JVS[17] = 0.25*B[61];
/* JVS(18) = Jac_FULL(4,10)                                         */
  JVS[18] = 0.9*B[44];
/* JVS(19) = Jac_FULL(4,13)                                         */
  JVS[19] = 0.9*B[45]+0.15*B[51]+0.25*B[62];
/* JVS(20) = Jac_FULL(4,14)                                         */
  JVS[20] = 0.15*B[52];
/* JVS(21) = Jac_FULL(4,17)                                         */
  JVS[21] = 0.363*B[66];
/* JVS(22) = Jac_FULL(4,18)                                         */
  JVS[22] = -B[33]+B[37];
/* JVS(23) = Jac_FULL(5,5)                                          */
  JVS[23] = -B[2]-B[3]-B[36]-B[74];
/* JVS(24) = Jac_FULL(5,6)                                          */
  JVS[24] = 0.69*B[55];
/* JVS(25) = Jac_FULL(5,7)                                          */
  JVS[25] = 0.785*B[65];
/* JVS(26) = Jac_FULL(5,8)                                          */
  JVS[26] = 0.25*B[61];
/* JVS(27) = Jac_FULL(5,9)                                          */
  JVS[27] = B[63];
/* JVS(28) = Jac_FULL(5,10)                                         */
  JVS[28] = 1.5*B[40]+0.9*B[44];
/* JVS(29) = Jac_FULL(5,13)                                         */
  JVS[29] = B[28]+0.9*B[45]+0.1*B[51]+0.25*B[62];
/* JVS(30) = Jac_FULL(5,14)                                         */
  JVS[30] = 0.75*B[46]+0.1*B[52];
/* JVS(31) = Jac_FULL(5,15)                                         */
  JVS[31] = B[29]+1.5*B[41]+0.75*B[47];
/* JVS(32) = Jac_FULL(5,17)                                         */
  JVS[32] = 0.785*B[66];
/* JVS(33) = Jac_FULL(5,18)                                         */
  JVS[33] = -B[37]+B[64];
/* JVS(34) = Jac_FULL(6,6)                                          */
  JVS[34] = -B[53]-B[55]-B[75];
/* JVS(35) = Jac_FULL(6,10)                                         */
  JVS[35] = B[42];
/* JVS(36) = Jac_FULL(6,12)                                         */
  JVS[36] = B[43]+B[48];
/* JVS(37) = Jac_FULL(6,14)                                         */
  JVS[37] = B[49];
/* JVS(38) = Jac_FULL(6,18)                                         */
  JVS[38] = -B[54];
/* JVS(39) = Jac_FULL(7,7)                                          */
  JVS[39] = -B[38]-B[65]-B[67];
/* JVS(40) = Jac_FULL(7,16)                                         */
  JVS[40] = -B[68];
/* JVS(41) = Jac_FULL(7,17)                                         */
  JVS[41] = -B[66];
/* JVS(42) = Jac_FULL(7,18)                                         */
  JVS[42] = -B[39];
/* JVS(43) = Jac_FULL(8,3)                                          */
  JVS[43] = B[56];
/* JVS(44) = Jac_FULL(8,6)                                          */
  JVS[44] = 0.4*B[53];
/* JVS(45) = Jac_FULL(8,8)                                          */
  JVS[45] = -B[59]-B[61];
/* JVS(46) = Jac_FULL(8,10)                                         */
  JVS[46] = 0;
/* JVS(47) = Jac_FULL(8,12)                                         */
  JVS[47] = -B[60];
/* JVS(48) = Jac_FULL(8,13)                                         */
  JVS[48] = -B[62];
/* JVS(49) = Jac_FULL(8,14)                                         */
  JVS[49] = 0;
/* JVS(50) = Jac_FULL(8,15)                                         */
  JVS[50] = 0;
/* JVS(51) = Jac_FULL(8,18)                                         */
  JVS[51] = 0.4*B[54]+B[57];
/* JVS(52) = Jac_FULL(9,6)                                          */
  JVS[52] = 0.69*B[55];
/* JVS(53) = Jac_FULL(9,7)                                          */
  JVS[53] = 0.675*B[65];
/* JVS(54) = Jac_FULL(9,9)                                          */
  JVS[54] = -B[63];
/* JVS(55) = Jac_FULL(9,10)                                         */
  JVS[55] = 0.75*B[40]+0.9*B[44];
/* JVS(56) = Jac_FULL(9,12)                                         */
  JVS[56] = 0;
/* JVS(57) = Jac_FULL(9,13)                                         */
  JVS[57] = 0.9*B[45];
/* JVS(58) = Jac_FULL(9,14)                                         */
  JVS[58] = 0;
/* JVS(59) = Jac_FULL(9,15)                                         */
  JVS[59] = 0.75*B[41];
/* JVS(60) = Jac_FULL(9,16)                                         */
  JVS[60] = 0;
/* JVS(61) = Jac_FULL(9,17)                                         */
  JVS[61] = 0.675*B[66];
/* JVS(62) = Jac_FULL(9,18)                                         */
  JVS[62] = -B[64];
/* JVS(63) = Jac_FULL(10,7)                                         */
  JVS[63] = 0.6*B[38];
/* JVS(64) = Jac_FULL(10,10)                                        */
  JVS[64] = -B[40]-B[42]-B[44];
/* JVS(65) = Jac_FULL(10,12)                                        */
  JVS[65] = -B[43];
/* JVS(66) = Jac_FULL(10,13)                                        */
  JVS[66] = -B[45];
/* JVS(67) = Jac_FULL(10,15)                                        */
  JVS[67] = -B[41];
/* JVS(68) = Jac_FULL(10,16)                                        */
  JVS[68] = 0;
/* JVS(69) = Jac_FULL(10,17)                                        */
  JVS[69] = 0;
/* JVS(70) = Jac_FULL(10,18)                                        */
  JVS[70] = 0.6*B[39];
/* JVS(71) = Jac_FULL(11,2)                                         */
  JVS[71] = B[25];
/* JVS(72) = Jac_FULL(11,8)                                         */
  JVS[72] = 0.9*B[61];
/* JVS(73) = Jac_FULL(11,10)                                        */
  JVS[73] = 0.9*B[44];
/* JVS(74) = Jac_FULL(11,11)                                        */
  JVS[74] = -B[1]-B[13]-B[15]-B[19]-B[73];
/* JVS(75) = Jac_FULL(11,12)                                        */
  JVS[75] = B[17]+B[23];
/* JVS(76) = Jac_FULL(11,13)                                        */
  JVS[76] = B[11]+B[18]+2*B[21]+0.99*B[28]+0.9*B[45]+0.9*B[51]+0.9
           *B[62];
/* JVS(77) = Jac_FULL(11,14)                                        */
  JVS[77] = 0.9*B[52];
/* JVS(78) = Jac_FULL(11,15)                                        */
  JVS[78] = 0.99*B[29];
/* JVS(79) = Jac_FULL(11,16)                                        */
  JVS[79] = -B[20]+2*B[22]+B[24]+0.89*B[27];
/* JVS(80) = Jac_FULL(11,17)                                        */
  JVS[80] = B[12]-B[16];
/* JVS(81) = Jac_FULL(11,18)                                        */
  JVS[81] = -B[14];
/* JVS(82) = Jac_FULL(12,3)                                         */
  JVS[82] = B[58];
/* JVS(83) = Jac_FULL(12,4)                                         */
  JVS[83] = B[32];
/* JVS(84) = Jac_FULL(12,5)                                         */
  JVS[84] = 2*B[2]+B[36];
/* JVS(85) = Jac_FULL(12,6)                                         */
  JVS[85] = B[55];
/* JVS(86) = Jac_FULL(12,7)                                         */
  JVS[86] = 0.125*B[65];
/* JVS(87) = Jac_FULL(12,8)                                         */
  JVS[87] = -B[59]+B[61];
/* JVS(88) = Jac_FULL(12,9)                                         */
  JVS[88] = B[63];
/* JVS(89) = Jac_FULL(12,10)                                        */
  JVS[89] = B[40]-B[42]+0.9*B[44];
/* JVS(90) = Jac_FULL(12,12)                                        */
  JVS[90] = -B[4]-2*B[8]-B[9]-B[17]-B[23]-B[30]-B[43]-B[48]-B[60];
/* JVS(91) = Jac_FULL(12,13)                                        */
  JVS[91] = -B[18]+B[28]+0.9*B[45]+0.97*B[51]+B[62];
/* JVS(92) = Jac_FULL(12,14)                                        */
  JVS[92] = B[46]-B[49]+B[50]+0.97*B[52];
/* JVS(93) = Jac_FULL(12,15)                                        */
  JVS[93] = B[29]-B[31]+B[41]+B[47];
/* JVS(94) = Jac_FULL(12,16)                                        */
  JVS[94] = -B[24];
/* JVS(95) = Jac_FULL(12,17)                                        */
  JVS[95] = B[6]-B[10]+0.125*B[66];
/* JVS(96) = Jac_FULL(12,18)                                        */
  JVS[96] = -B[5]+B[7]+B[33]+B[37]+B[64];
/* JVS(97) = Jac_FULL(13,8)                                         */
  JVS[97] = -B[61];
/* JVS(98) = Jac_FULL(13,10)                                        */
  JVS[98] = -B[44];
/* JVS(99) = Jac_FULL(13,11)                                        */
  JVS[99] = B[1];
/* JVS(100) = Jac_FULL(13,12)                                       */
  JVS[100] = -B[17];
/* JVS(101) = Jac_FULL(13,13)                                       */
  JVS[101] = -B[11]-B[18]-B[21]-B[28]-B[45]-B[51]-B[62]-B[72];
/* JVS(102) = Jac_FULL(13,14)                                       */
  JVS[102] = -B[52];
/* JVS(103) = Jac_FULL(13,15)                                       */
  JVS[103] = -B[29];
/* JVS(104) = Jac_FULL(13,16)                                       */
  JVS[104] = -B[22]+0.11*B[27];
/* JVS(105) = Jac_FULL(13,17)                                       */
  JVS[105] = -B[12];
/* JVS(106) = Jac_FULL(13,18)                                       */
  JVS[106] = 0;
/* JVS(107) = Jac_FULL(14,7)                                        */
  JVS[107] = 0.4*B[38];
/* JVS(108) = Jac_FULL(14,12)                                       */
  JVS[108] = -B[48];
/* JVS(109) = Jac_FULL(14,13)                                       */
  JVS[109] = -B[51];
/* JVS(110) = Jac_FULL(14,14)                                       */
  JVS[110] = -B[46]-B[49]-B[50]-B[52];
/* JVS(111) = Jac_FULL(14,15)                                       */
  JVS[111] = -B[47];
/* JVS(112) = Jac_FULL(14,16)                                       */
  JVS[112] = 0;
/* JVS(113) = Jac_FULL(14,17)                                       */
  JVS[113] = 0;
/* JVS(114) = Jac_FULL(14,18)                                       */
  JVS[114] = 0.4*B[39];
/* JVS(115) = Jac_FULL(15,0)                                        */
  JVS[115] = B[34];
/* JVS(116) = Jac_FULL(15,7)                                        */
  JVS[116] = 0.181*B[65];
/* JVS(117) = Jac_FULL(15,9)                                        */
  JVS[117] = 1.5*B[63];
/* JVS(118) = Jac_FULL(15,10)                                       */
  JVS[118] = -B[40];
/* JVS(119) = Jac_FULL(15,12)                                       */
  JVS[119] = -B[30];
/* JVS(120) = Jac_FULL(15,13)                                       */
  JVS[120] = -B[28];
/* JVS(121) = Jac_FULL(15,14)                                       */
  JVS[121] = -B[46];
/* JVS(122) = Jac_FULL(15,15)                                       */
  JVS[122] = -B[29]-B[31]-B[41]-B[47];
/* JVS(123) = Jac_FULL(15,16)                                       */
  JVS[123] = 0;
/* JVS(124) = Jac_FULL(15,17)                                       */
  JVS[124] = 0.181*B[66];
/* JVS(125) = Jac_FULL(15,18)                                       */
  JVS[125] = B[35]+1.5*B[64];
/* JVS(126) = Jac_FULL(16,2)                                        */
  JVS[126] = B[25];
/* JVS(127) = Jac_FULL(16,7)                                        */
  JVS[127] = -B[67];
/* JVS(128) = Jac_FULL(16,11)                                       */
  JVS[128] = B[15]-B[19];
/* JVS(129) = Jac_FULL(16,12)                                       */
  JVS[129] = -B[23];
/* JVS(130) = Jac_FULL(16,13)                                       */
  JVS[130] = -B[21];
/* JVS(131) = Jac_FULL(16,14)                                       */
  JVS[131] = 0;
/* JVS(132) = Jac_FULL(16,15)                                       */
  JVS[132] = 0;
/* JVS(133) = Jac_FULL(16,16)                                       */
  JVS[133] = -B[20]-B[22]-B[24]-B[27]-B[68];
/* JVS(134) = Jac_FULL(16,17)                                       */
  JVS[134] = B[16];
/* JVS(135) = Jac_FULL(16,18)                                       */
  JVS[135] = 0;
/* JVS(136) = Jac_FULL(17,7)                                        */
  JVS[136] = -B[65];
/* JVS(137) = Jac_FULL(17,11)                                       */
  JVS[137] = B[1]-B[15];
/* JVS(138) = Jac_FULL(17,12)                                       */
  JVS[138] = -B[9];
/* JVS(139) = Jac_FULL(17,13)                                       */
  JVS[139] = -B[11];
/* JVS(140) = Jac_FULL(17,14)                                       */
  JVS[140] = 0;
/* JVS(141) = Jac_FULL(17,15)                                       */
  JVS[141] = 0;
/* JVS(142) = Jac_FULL(17,16)                                       */
  JVS[142] = 0.89*B[27];
/* JVS(143) = Jac_FULL(17,17)                                       */
  JVS[143] = -B[0]-B[6]-B[10]-B[12]-B[16]-B[66]-B[71];
/* JVS(144) = Jac_FULL(17,18)                                       */
  JVS[144] = -B[7];
/* JVS(145) = Jac_FULL(18,0)                                        */
  JVS[145] = -B[34];
/* JVS(146) = Jac_FULL(18,1)                                        */
  JVS[146] = 2*B[69];
/* JVS(147) = Jac_FULL(18,3)                                        */
  JVS[147] = -B[56]+B[58];
/* JVS(148) = Jac_FULL(18,4)                                        */
  JVS[148] = -B[32];
/* JVS(149) = Jac_FULL(18,5)                                        */
  JVS[149] = -B[36];
/* JVS(150) = Jac_FULL(18,6)                                        */
  JVS[150] = -0.4*B[53];
/* JVS(151) = Jac_FULL(18,7)                                        */
  JVS[151] = -B[38]+0.205*B[65];
/* JVS(152) = Jac_FULL(18,8)                                        */
  JVS[152] = 0;
/* JVS(153) = Jac_FULL(18,9)                                        */
  JVS[153] = -B[63];
/* JVS(154) = Jac_FULL(18,10)                                       */
  JVS[154] = 0;
/* JVS(155) = Jac_FULL(18,11)                                       */
  JVS[155] = -B[13];
/* JVS(156) = Jac_FULL(18,12)                                       */
  JVS[156] = -B[4]+B[9]+B[17]+B[23];
/* JVS(157) = Jac_FULL(18,13)                                       */
  JVS[157] = B[18];
/* JVS(158) = Jac_FULL(18,14)                                       */
  JVS[158] = 0;
/* JVS(159) = Jac_FULL(18,15)                                       */
  JVS[159] = 0;
/* JVS(160) = Jac_FULL(18,16)                                       */
  JVS[160] = B[24];
/* JVS(161) = Jac_FULL(18,17)                                       */
  JVS[161] = 2*B[0]-B[6]+B[10]+0.205*B[66];
/* JVS(162) = Jac_FULL(18,18)                                       */
  JVS[162] = -B[5]-B[7]-B[14]-B[33]-B[35]-B[37]-B[39]-0.4*B[54]-B[57]
            -B[64];
}

/* End of Jac_SP function                                           */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/* Jac_SP_Vec - function for sparse multiplication: sparse Jacobian times vector */
/*   Arguments :                                                    */
/*      JVS       - sparse Jacobian of variables                    */
/*      UV        - User vector for variables                       */
/*      JUV       - Jacobian times user vector                      */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void Jac_SP_Vec( 
  double JVS[],                          /* sparse Jacobian of variables */
  double UV[],                           /* User vector for variables */
  double JUV[]                           /* Jacobian times user vector */
)
{
  JUV[0] = JVS[0]*UV[0]+JVS[1]*UV[18];
  JUV[1] = JVS[2]*UV[1]+JVS[3]*UV[7]+JVS[4]*UV[12]+JVS[5]*UV[17];
  JUV[2] = JVS[6]*UV[2]+JVS[7]*UV[11]+JVS[8]*UV[16];
  JUV[3] = JVS[9]*UV[3]+JVS[10]*UV[13]+JVS[11]*UV[14]+JVS[12]*UV[15]
          +JVS[13]*UV[18];
  JUV[4] = JVS[14]*UV[4]+JVS[15]*UV[5]+JVS[16]*UV[7]+JVS[17]*UV[8]
          +JVS[18]*UV[10]+JVS[19]*UV[13]+JVS[20]*UV[14]+JVS[21]*UV[17]
          +JVS[22]*UV[18];
  JUV[5] = JVS[23]*UV[5]+JVS[24]*UV[6]+JVS[25]*UV[7]+JVS[26]*UV[8]
          +JVS[27]*UV[9]+JVS[28]*UV[10]+JVS[29]*UV[13]+JVS[30]*UV[14]
          +JVS[31]*UV[15]+JVS[32]*UV[17]+JVS[33]*UV[18];
  JUV[6] = JVS[34]*UV[6]+JVS[35]*UV[10]+JVS[36]*UV[12]+JVS[37]*UV[14]
          +JVS[38]*UV[18];
  JUV[7] = JVS[39]*UV[7]+JVS[40]*UV[16]+JVS[41]*UV[17]+JVS[42]*UV[18];
  JUV[8] = JVS[43]*UV[3]+JVS[44]*UV[6]+JVS[45]*UV[8]+JVS[47]*UV[12]
          +JVS[48]*UV[13]+JVS[51]*UV[18];
  JUV[9] = JVS[52]*UV[6]+JVS[53]*UV[7]+JVS[54]*UV[9]+JVS[55]*UV[10]
          +JVS[57]*UV[13]+JVS[59]*UV[15]+JVS[61]*UV[17]+JVS[62]*UV[18];
  JUV[10] = JVS[63]*UV[7]+JVS[64]*UV[10]+JVS[65]*UV[12]+JVS[66]*UV[13]
           +JVS[67]*UV[15]+JVS[70]*UV[18];
  JUV[11] = JVS[71]*UV[2]+JVS[72]*UV[8]+JVS[73]*UV[10]+JVS[74]*UV[11]
           +JVS[75]*UV[12]+JVS[76]*UV[13]+JVS[77]*UV[14]+JVS[78]*UV[15]
           +JVS[79]*UV[16]+JVS[80]*UV[17]+JVS[81]*UV[18];
  JUV[12] = JVS[82]*UV[3]+JVS[83]*UV[4]+JVS[84]*UV[5]+JVS[85]*UV[6]
           +JVS[86]*UV[7]+JVS[87]*UV[8]+JVS[88]*UV[9]+JVS[89]*UV[10]
           +JVS[90]*UV[12]+JVS[91]*UV[13]+JVS[92]*UV[14]+JVS[93]*UV[15]
           +JVS[94]*UV[16]+JVS[95]*UV[17]+JVS[96]*UV[18];
  JUV[13] = JVS[97]*UV[8]+JVS[98]*UV[10]+JVS[99]*UV[11]+JVS[100]*UV[12]
           +JVS[101]*UV[13]+JVS[102]*UV[14]+JVS[103]*UV[15]+JVS[104]
           *UV[16]+JVS[105]*UV[17];
  JUV[14] = JVS[107]*UV[7]+JVS[108]*UV[12]+JVS[109]*UV[13]+JVS[110]
           *UV[14]+JVS[111]*UV[15]+JVS[114]*UV[18];
  JUV[15] = JVS[115]*UV[0]+JVS[116]*UV[7]+JVS[117]*UV[9]+JVS[118]
           *UV[10]+JVS[119]*UV[12]+JVS[120]*UV[13]+JVS[121]*UV[14]
           +JVS[122]*UV[15]+JVS[124]*UV[17]+JVS[125]*UV[18];
  JUV[16] = JVS[126]*UV[2]+JVS[127]*UV[7]+JVS[128]*UV[11]+JVS[129]
           *UV[12]+JVS[130]*UV[13]+JVS[133]*UV[16]+JVS[134]*UV[17];
  JUV[17] = JVS[136]*UV[7]+JVS[137]*UV[11]+JVS[138]*UV[12]+JVS[139]
           *UV[13]+JVS[142]*UV[16]+JVS[143]*UV[17]+JVS[144]*UV[18];
  JUV[18] = JVS[145]*UV[0]+JVS[146]*UV[1]+JVS[147]*UV[3]+JVS[148]*UV[4]
           +JVS[149]*UV[5]+JVS[150]*UV[6]+JVS[151]*UV[7]+JVS[153]*UV[9]
           +JVS[155]*UV[11]+JVS[156]*UV[12]+JVS[157]*UV[13]+JVS[160]
           *UV[16]+JVS[161]*UV[17]+JVS[162]*UV[18];
}

/* End of Jac_SP_Vec function                                       */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/* JacTR_SP_Vec - sparse multiplication: sparse Jacobian transposed times vector */
/*   Arguments :                                                    */
/*      JVS       - sparse Jacobian of variables                    */
/*      UV        - User vector for variables                       */
/*      JTUV      - Jacobian transposed times user vector           */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void JacTR_SP_Vec( 
  double JVS[],                          /* sparse Jacobian of variables */
  double UV[],                           /* User vector for variables */
  double JTUV[]                          /* Jacobian transposed times user vector */
)
{
  JTUV[0] = JVS[0]*UV[0]+JVS[115]*UV[15]+JVS[145]*UV[18];
  JTUV[1] = JVS[2]*UV[1]+JVS[146]*UV[18];
  JTUV[2] = JVS[6]*UV[2]+JVS[71]*UV[11]+JVS[126]*UV[16];
  JTUV[3] = JVS[9]*UV[3]+JVS[43]*UV[8]+JVS[82]*UV[12]+JVS[147]*UV[18];
  JTUV[4] = JVS[14]*UV[4]+JVS[83]*UV[12]+JVS[148]*UV[18];
  JTUV[5] = JVS[15]*UV[4]+JVS[23]*UV[5]+JVS[84]*UV[12]+JVS[149]*UV[18];
  JTUV[6] = JVS[24]*UV[5]+JVS[34]*UV[6]+JVS[44]*UV[8]+JVS[52]*UV[9]
           +JVS[85]*UV[12]+JVS[150]*UV[18];
  JTUV[7] = JVS[3]*UV[1]+JVS[16]*UV[4]+JVS[25]*UV[5]+JVS[39]*UV[7]
           +JVS[53]*UV[9]+JVS[63]*UV[10]+JVS[86]*UV[12]+JVS[107]*UV[14]
           +JVS[116]*UV[15]+JVS[127]*UV[16]+JVS[136]*UV[17]+JVS[151]
           *UV[18];
  JTUV[8] = JVS[17]*UV[4]+JVS[26]*UV[5]+JVS[45]*UV[8]+JVS[72]*UV[11]
           +JVS[87]*UV[12]+JVS[97]*UV[13];
  JTUV[9] = JVS[27]*UV[5]+JVS[54]*UV[9]+JVS[88]*UV[12]+JVS[117]*UV[15]
           +JVS[153]*UV[18];
  JTUV[10] = JVS[18]*UV[4]+JVS[28]*UV[5]+JVS[35]*UV[6]+JVS[55]*UV[9]
            +JVS[64]*UV[10]+JVS[73]*UV[11]+JVS[89]*UV[12]+JVS[98]
            *UV[13]+JVS[118]*UV[15];
  JTUV[11] = JVS[7]*UV[2]+JVS[74]*UV[11]+JVS[99]*UV[13]+JVS[128]*UV[16]
            +JVS[137]*UV[17]+JVS[155]*UV[18];
  JTUV[12] = JVS[4]*UV[1]+JVS[36]*UV[6]+JVS[47]*UV[8]+JVS[65]*UV[10]
            +JVS[75]*UV[11]+JVS[90]*UV[12]+JVS[100]*UV[13]+JVS[108]
            *UV[14]+JVS[119]*UV[15]+JVS[129]*UV[16]+JVS[138]*UV[17]
            +JVS[156]*UV[18];
  JTUV[13] = JVS[10]*UV[3]+JVS[19]*UV[4]+JVS[29]*UV[5]+JVS[48]*UV[8]
            +JVS[57]*UV[9]+JVS[66]*UV[10]+JVS[76]*UV[11]+JVS[91]*UV[12]
            +JVS[101]*UV[13]+JVS[109]*UV[14]+JVS[120]*UV[15]+JVS[130]
            *UV[16]+JVS[139]*UV[17]+JVS[157]*UV[18];
  JTUV[14] = JVS[11]*UV[3]+JVS[20]*UV[4]+JVS[30]*UV[5]+JVS[37]*UV[6]
            +JVS[77]*UV[11]+JVS[92]*UV[12]+JVS[102]*UV[13]+JVS[110]
            *UV[14]+JVS[121]*UV[15];
  JTUV[15] = JVS[12]*UV[3]+JVS[31]*UV[5]+JVS[59]*UV[9]+JVS[67]*UV[10]
            +JVS[78]*UV[11]+JVS[93]*UV[12]+JVS[103]*UV[13]+JVS[111]
            *UV[14]+JVS[122]*UV[15];
  JTUV[16] = JVS[8]*UV[2]+JVS[40]*UV[7]+JVS[79]*UV[11]+JVS[94]*UV[12]
            +JVS[104]*UV[13]+JVS[133]*UV[16]+JVS[142]*UV[17]+JVS[160]
            *UV[18];
  JTUV[17] = JVS[5]*UV[1]+JVS[21]*UV[4]+JVS[32]*UV[5]+JVS[41]*UV[7]
            +JVS[61]*UV[9]+JVS[80]*UV[11]+JVS[95]*UV[12]+JVS[105]
            *UV[13]+JVS[124]*UV[15]+JVS[134]*UV[16]+JVS[143]*UV[17]
            +JVS[161]*UV[18];
  JTUV[18] = JVS[1]*UV[0]+JVS[13]*UV[3]+JVS[22]*UV[4]+JVS[33]*UV[5]
            +JVS[38]*UV[6]+JVS[42]*UV[7]+JVS[51]*UV[8]+JVS[62]*UV[9]
            +JVS[70]*UV[10]+JVS[81]*UV[11]+JVS[96]*UV[12]+JVS[114]
            *UV[14]+JVS[125]*UV[15]+JVS[144]*UV[17]+JVS[162]*UV[18];
}

/* End of JacTR_SP_Vec function                                     */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


