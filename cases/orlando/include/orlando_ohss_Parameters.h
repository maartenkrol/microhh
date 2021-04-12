/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/* Parameter Header File                                            */
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
/* File                 : orlando_ohss_Parameters.h                 */
/* Time                 : Sun Aug  9 15:18:54 2020                  */
/* Working directory    : /home/WUR/krol005/kpp/examples            */
/* Equation file        : orlando_ohss.kpp                          */
/* Output root filename : orlando_ohss                              */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */




#define NSPEC                22          /* Number of chemical species */
#define NVAR                 18          /* Number of Variable species */
#define NVARACT              18          /* Number of Active species */
#define NFIX                 4           /* Number of Fixed species */
#define NREACT               51          /* Number of reactions */
#define NVARST               0           /* Starting of variables in conc. vect. */
#define NFIXST               18          /* Starting of fixed in conc. vect. */
#define NONZERO              114         /* Number of nonzero entries in Jacobian */
#define LU_NONZERO           129         /* Number of nonzero entries in LU factoriz. of Jacobian */
#define CNVAR                19          /* (NVAR+1) Number of elements in compressed row format */
#define CNEQN                52          /* (NREACT+1) Number stoicm elements in compressed col format */
#define NHESS                77          /* Length of Sparse Hessian */
#define NLOOKAT              22          /* Number of species to look at */
#define NMONITOR             1           /* Number of species to monitor */
#define NMASS                1           /* Number of atoms to check mass balance */

/* Index declaration for variable species in C and VAR              */
/*   VAR(ind_spc) = C(ind_spc)                                      */

#define ind_H2O2             0          
#define ind_CH4              1          
#define ind_N2O5             2          
#define ind_HALD             3          
#define ind_CO               4          
#define ind_HCHO             5          
#define ind_ISOPOOH          6          
#define ind_ISOP             7          
#define ind_MVKMACR          8          
#define ind_XO2              9          
#define ind_ISOPAO2          10         
#define ind_NO2              11         
#define ind_O3               12         
#define ind_NO               13         
#define ind_CH3O2            14         
#define ind_ISOPBO2          15         
#define ind_NO3              16         
#define ind_HO2              17         

/* Index declaration for fixed species in C                         */
/*   C(ind_spc)                                                     */

#define ind_OH               18         
#define ind_H2O              19         
#define ind_M                20         
#define ind_DUMMY            21         

/* Index declaration for fixed species in FIX                       */
/*    FIX(indf_spc) = C(ind_spc) = C(NVAR+indf_spc)                 */

#define indf_OH              0          
#define indf_H2O             1          
#define indf_M               2          
#define indf_DUMMY           3          

#define NJVRP                67          /* Length of sparse Jacobian JVRP */

#define NSTOICM              136         /* Length of Sparse Stoichiometric Matrix */