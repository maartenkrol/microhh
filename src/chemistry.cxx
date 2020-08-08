
/*
 * MicroHH
 * Copyright (c) 2011-2020 Chiel van Heerwaarden
 * Copyright (c) 2011-2020 Thijs Heus
 * Copyright (c) 2014-2020 Bart van Stratum
 *
 * This file is part of MicroHH
 *
 * MicroHH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * MicroHH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
 */

//#include <cstdio>
#include <algorithm>
#include <iostream>
#include <math.h>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "thermo.h"
#include "stats.h"
#include "chemistry.h"
#include "constants.h" 


namespace 
{
#include "/Users/krol/kpp/examples/orlando_ohss_Parameters.h"
#include "/Users/krol/kpp/examples/orlando_ohss_Global.h"
#include "/Users/krol/kpp/examples/orlando_ohss_Sparse.h"
#include "/Users/krol/kpp/examples/orlando_ohss_Integrator.c"
#include "/Users/krol/kpp/examples/orlando_ohss_Function.c"
#include "/Users/krol/kpp/examples/orlando_ohss_LinearAlgebra.c"
#include "/Users/krol/kpp/examples/orlando_ohss_JacobianSP.c"
#include "/Users/krol/kpp/examples/orlando_ohss_Jacobian.c"
	

double C[NSPEC];                         /* Concentration of all species */
double * VAR = & C[0];
double * FIX = & C[19];
double RCONST[NREACT];                   /* Rate constants (global) */
double TIME;                             /* Current integration time */
double SUN;                              /* Sunlight intensity between [0,1] */
double TEMP;                             /* Temperature */
double RTOLS;                            /* (scalar) Relative tolerance */
double TSTART;                           /* Integration start time */
double TEND;                             /* Integration end time */
double DT;                               /* Integration step */
double ATOL[NVAR];                       /* Absolute tolerance */
double RTOL[NVAR];                       /* Relative tolerance */
double STEPMIN;                          /* Lower bound for integration step */
double STEPMAX;                          /* Upper bound for integration step */
double CFACTOR;                          /* Conversion factor for concentration units */

	template<typename TF>
	TF  ARRM( TF A0, TF B0, TF C0, TF TEMP )
	{
	return A0 * exp( -B0/TEMP ) 
		* pow( (TEMP/300.0), C0 );   
	}           

	template<typename TF>
	TF usr_O3_hv_H2O( TF TEMP, TF C_M, TF C_H2O, TF J_O1D)
	{
	TF KH2O;
	TF KN2;
	TF KO2;
	KH2O = ((TF)1.63e-10 *C_H2O * exp((TF)60.0/TEMP)  )  ;
	KN2  = ((TF)2.15e-11 * exp((TF)110.0/TEMP) *(TF)0.79*C_M) ;
	KO2  = ((TF)3.30e-11 * exp((TF)55.0 /TEMP) *(TF)0.21*C_M) ;
	return (KH2O *J_O1D) / (KH2O + KN2 + KO2);
	}           

	template<typename TF>
	TF usr_HO2_HO2( TF TEMP, TF C_M, TF C_H2O )
	/* for cesm-consistent reaction labels, equivalent to usr9 */
	/* HO2+HO2 -> H2O2+O2 */
	/* H2O included in fc calculation, not as reactant */
	{
	TF ko;
	TF kinf;
	TF fc;
	TF kr;

	if( C_H2O > (TF)0.0 ) 
	  {
	  ko   = (TF)2.3e-13 * exp( (TF)600./TEMP );
	  kinf = (TF)1.7e-33 * C_M * exp( (TF)1000./TEMP );
	  fc = (TF)1.0 + (TF)1.4e-21 *C_H2O* exp( (TF)2200./TEMP );
	  kr = (ko + kinf) * fc; 
	  }
	else
	  {
	  kr = (TF)0.0 ;
	  }
	return kr;
	}           

	template<typename TF>
	TF JPL_TROE( TF k0_300K, TF n, TF kinf_300K, TF m, 
		 TF base, TF temp, TF cair )

	{
	/* !------------------------------------------------------------ */
	/* ! ... local variables */
	/* !------------------------------------------------------------ */
	TF zt_help;
	TF k0_T;
	TF kinf_T;
	TF k_ratio;

	zt_help = (TF)300./temp;
	k0_T    = k0_300K   * pow(zt_help,n) * cair ; /* ! k_0   at current T */
	kinf_T  = kinf_300K * pow(zt_help,m)        ; /* ! k_inf at current T */
	k_ratio = k0_T/kinf_T;

	return  k0_T/((TF)1.+k_ratio) * 
	      pow(base,(TF)1.0/
	      ((TF)1.0+pow(log10(k_ratio),(TF)2.0)));
	}           

	template<typename TF>
	TF usr_CO_OH_a( TF temp, TF c_m )
	/* ! for cesm-consistent reaction labels, equivalent to usr8 */
	/* ! CO+OH -> CO2+HO2 */
	{
	TF boltz = 1.38044e-16;
	return (TF)1.5e-13 * ((TF)1.+ (TF)6.e-7*boltz*c_m*temp);
	}

	template<typename TF>
	TF usr_N2O5_H2O( TF k, TF c_h2o )
	{
	return k * c_h2o;
	}

	template<typename TF>
	TF ARR2M( TF A0, TF B0, TF TEMP)
	{
	return A0 * exp( -B0/TEMP );
	}

    template<typename TF>
    void pss(
            TF* restrict tch4, const TF* ch4, 
            TF* restrict th2o2, const TF* h2o2, 
            TF* restrict tn2o5, const TF* n2o5, 
            TF* restrict thald, const TF* hald, 
            TF* restrict tco, const TF* co, 
            TF* restrict thcho, const TF* hcho, 
            TF* restrict tisopooh, const TF* isopooh, 
            TF* restrict tisop, const TF* isop, 
            TF* restrict txo2, const TF* xo2, 
            TF* restrict tmvkmacr, const TF* mvkmacr, 
            TF* restrict tisopao2, const TF* isopao2, 
            TF* restrict tno2, const TF* no2, 
            TF* restrict tho2, const TF* ho2, 
            TF* restrict tno, const TF* no, 
            TF* restrict tisopbo2, const TF* isopbo2, 
            TF* restrict tch3o2, const TF* ch3o2, 
            TF* restrict tno3, const TF* no3, 
            TF* restrict to3, const TF* o3, 
            TF*oh,
	    const TF* restrict qt,
	    const TF* restrict Temp, const TF rkdt, const TF switch_dt,
	    const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk, const TF* const restrict dz, TF* const restrict rhoref)
    {
        const TF jval[] = {4.964E-05, 8.881E-06, 1.172E-02, 2.620E-02, 2.075E-01,
      	                 3.848E-05, 5.544E-05, 6.732E-06}; 
	const int Pj_o31d = 0;
	const int Pj_h2o2 = 1;
	const int Pj_no2  = 2;
	const int Pj_no3a = 3;
	const int Pj_no3b = 4;
	const int Pj_ch2or = 5;
	const int Pj_ch2om = 6;
	const int Pj_ch3o2h = 7;
	const TF RTOLS = 1e-3;
	const TF xmh2o = 18.015265;
        const TF xmair = 28.9647;       // Molar mass of dry air  [kg kmol-1]
        const TF Na    = 6.02214086e23; // Avogadros number [molecules mol-1]
	TF C_M = 2.55e19;
	TF C_H2O = 0.01*C_M;
	TF tscale[NVAR] ;
	TF deriv[NVAR] ;
	TF VAR0[NVAR] ;

      
	for( int i = 0; i < NVAR; i++ ) {
	  RTOL[i] = RTOLS;
	  ATOL[i] = 1.0;
	  tscale[i] = 1e20;
	}
	TF STEPMIN = 0.01;
	TF STEPMAX = 90;
	
	TF A[41];
	TF V[19];
	TF TEMP = 0.0;

	TF vdo3 = 0.0;
	TF vdno = 0.0;
	TF vdno2 = 0.0;
	TF vdho2 = 0.0;
	TF vdoh = 0.0;
	VAR = &C[0];
	FIX = &C[18];
        int nkpp = 0;
	int nderiv = 0;
        for (int k=kstart; k<kend; ++k)
            {	
	    C_M = (TF)1e-3*rhoref[k]*Na/xmair;   // molecules/cm3 for chemistry!
	    const TF CFACTOR = C_M*(TF)1e-9 ;               // from ppb (units mixing ratio) to molecules/cm3
            if (k==kstart) {
	        vdo3  = TF(0.005)/dz[k];   // 1/s
	        vdno  = TF(0.002)/dz[k];   // 1/s
	        vdno2 = TF(0.005)/dz[k];   // 1/s
	        vdho2 = TF(0.010)/dz[k];   // 1/s
	        vdoh  = TF(0.010)/dz[k];   // 1/s
	    }
            else {
	        vdo3  = TF(0.0);
	        vdno  = TF(0.0);
	        vdno2 = TF(0.0);
	        vdho2 = TF(0.0);
	        vdoh  = TF(0.0);
	    }
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
		    C_H2O = qt[ijk]*xmair*C_M/xmh2o;                   // kg/kg --> molH2O/molAir --*C_M--> molecules/cm3
		    const TF SUN = 1.0; 
		    TEMP = Temp[ijk];
		    // convert to molecules per cm3:
                    VAR[0] = std::max(h2o2[ijk]*CFACTOR,(TF)0.0);
		    VAR[1] = std::max(ch4[ijk]*CFACTOR,(TF)0.0);
                    VAR[2] = std::max(n2o5[ijk]*CFACTOR,(TF)0.0);
                    VAR[3] = std::max(hald[ijk]*CFACTOR,(TF)0.0);
                    VAR[4] = std::max(co[ijk]*CFACTOR,(TF)0.0);
                    VAR[5] = std::max(hcho[ijk]*CFACTOR,(TF)0.0);
                    VAR[6] = std::max(isopooh[ijk]*CFACTOR,(TF)0.0);
                    VAR[7] = std::max(isop[ijk]*CFACTOR,(TF)0.0);
                    VAR[8] = std::max(mvkmacr[ijk]*CFACTOR,(TF)0.0);
                    VAR[9] = std::max(xo2[ijk]*CFACTOR,(TF)0.0);
                    VAR[10] = std::max(isopao2[ijk]*CFACTOR,(TF)0.0);
                    VAR[11] = std::max(no2[ijk]*CFACTOR,(TF)0.0);
                    VAR[12] = std::max(o3[ijk]*CFACTOR,(TF)0.0);
                    VAR[13] = std::max(no[ijk]*CFACTOR,(TF)0.0);
                    VAR[14] = std::max(ch3o2[ijk]*CFACTOR,(TF)0.0);
                    VAR[15] = std::max(isopbo2[ijk]*CFACTOR,(TF)0.0);
                    VAR[16] = std::max(no3[ijk]*CFACTOR,(TF)0.0);
                    VAR[17] = std::max(ho2[ijk]*CFACTOR,(TF)0.0);
		    
                    RCONST[0] = (usr_O3_hv_H2O(TEMP,C_M,C_H2O,SUN*jval[Pj_o31d]));
                    RCONST[1] = (SUN*jval[Pj_no2]);
                    RCONST[2] = (SUN*jval[Pj_ch2or]);
                    RCONST[3] = (SUN*jval[Pj_ch2om]);
                    RCONST[4] = (ARR2M((TF)4.8e-11,-(TF)250.0,TEMP));
                    RCONST[5] = (ARR2M((TF)1.70e-12,(TF)940.0,TEMP));
                    RCONST[6] = (usr_HO2_HO2(TEMP,C_M,C_H2O));
                    RCONST[7] = (ARRM((TF)2.03e-16,-(TF)693.0,(TF)4.57,TEMP));
                    RCONST[8] = (ARR2M((TF)1.4e-12,(TF)1310.0,TEMP));
                    RCONST[9] = (JPL_TROE((TF)1.8e-30,(TF)3.,(TF)2.8e-11,(TF)0.0,(TF)0.6,TEMP,C_M));
                    RCONST[10] = (ARR2M((TF)1.4e-13,(TF)2470.0,TEMP));
                    RCONST[11] = (ARR2M((TF)3.50e-12,-(TF)270.0,TEMP));
                    RCONST[12] = (JPL_TROE((TF)2.4e-30,(TF)3.0,(TF)1.6e-12,-(TF)0.10,
                                (TF)0.60,TEMP,C_M));
                    RCONST[13] = (ARR2M((TF)1.8E-11,-(TF)110.0,TEMP));
                    RCONST[14] = ((TF)4.0e-12);
		    RCONST[15] = RCONST[12]*(TF)1.724e+26 * exp(-(TF)10840.0 / TEMP);
                    RCONST[16] = (usr_N2O5_H2O((TF)2.5e-22,C_H2O));
                    RCONST[17] = (SUN*jval[Pj_no3a]+SUN*jval[Pj_no3b]);
                    RCONST[18] = (ARR2M((TF)2.80e-12,-(TF)300.0,TEMP));
                    RCONST[19] = (ARR2M((TF)4.10e-13,-(TF)750.0,TEMP));
                    RCONST[20] = (usr_CO_OH_a(TEMP,C_M));
                    RCONST[21] = (ARR2M((TF)2.45e-12,(TF)1775.0,TEMP));
                    RCONST[22] = (ARR2M((TF)5.50e-12,-(TF)125.0,TEMP));
                    RCONST[23] = (ARR2M((TF)2.54e-11,-(TF)410.0,TEMP));
                    RCONST[24] = (ARR2M((TF)5.0e-13,-(TF)400.0,TEMP));
                    RCONST[25] = (ARR2M((TF)8.0e-13,-(TF)700.0,TEMP));
                    RCONST[26] = (ARR2M((TF)4.4e-12,-(TF)180.0,TEMP));
                    RCONST[27] = (ARR2M((TF)5.0e-13,-(TF)400.0,TEMP));
                    RCONST[28] = (ARR2M((TF)8.0e-13,-(TF)700.0,TEMP));
                    RCONST[29] = (ARR2M((TF)1.6e+09,(TF)8300.0,TEMP));
                    RCONST[30] = (ARR2M((TF)4.4e-12,-(TF)180.0,TEMP));
                    RCONST[31] = (ARR2M((TF)1.52e-11,-(TF)200.0,TEMP));
                    RCONST[32] = (SUN*jval[Pj_ch3o2h]);
                    RCONST[33] = (ARR2M((TF)1.86e-11,-(TF)175.0,TEMP));
                    RCONST[34] = ((TF)0.004*SUN*jval[Pj_no2]);
                    RCONST[35] = (ARR2M((TF)8.0e-13,-(TF)700.0,TEMP));
                    RCONST[36] = (ARR2M((TF)2.7e-12,-(TF)360.0,TEMP));
                    RCONST[37] = ((TF)2.4e-11);
                    RCONST[38] = (ARR2M((TF)1.03E-14,(TF)1995.0,TEMP));
                    RCONST[39] = (ARR2M((TF)3.15e-12,(TF)450.0,TEMP));
                    RCONST[40] = (SUN*jval[Pj_h2o2]);

		    FIX[0] = ((TF)2.0*RCONST[0]*VAR[12] + RCONST[7]*VAR[12]*VAR[17] + RCONST[11]*VAR[13]*VAR[17] + RCONST[14]*VAR[16]*VAR[17] + (TF)2.0*RCONST[40]*VAR[0])/
			        ( RCONST[4]*VAR[17] +  RCONST[5]*VAR[12] +     RCONST[9]*VAR[11] +  RCONST[20]*VAR[4] + RCONST[21]*VAR[1] + 
		       	          RCONST[22]*VAR[5] +  RCONST[23]*VAR[7] + 0.4*RCONST[31]*VAR[6] +  RCONST[33]*VAR[3] + RCONST[37]*VAR[8]);
		    FIX[1] = C_H2O;
                    FIX[2] = C_M;
		    oh[ijk] = FIX[0]/CFACTOR;

		    Fun( VAR, FIX, RCONST, deriv);
		    // for (int l=0; l<NVAR; ++l) printf (" %i %13.3e %13.3e %13.3e \n", l,VAR[l],deriv[l],VAR[l]/ABS(deriv[l]));
		    TF mint = (TF)1e20;
		    for (int l=0; l<NVAR; ++l)
			    if (ABS(deriv[l]) > (TF)1e-5 && VAR[l]> (TF)1e-5) mint = std::min(mint,VAR[l]/ABS(deriv[l])); 
			//printf (" %i %13.3e \n ", ijk, mint);

		    if (mint < switch_dt)
		    {
			    nkpp += 1;

			    WCOPY(NVAR,VAR,1,VAR0,1);
			    INTEGRATE(  (TF)0.0 , rkdt );
			    
			    th2o2[ijk] =    (VAR[0]-VAR0[0])/(rkdt*CFACTOR);
			    tch4[ijk] =     (VAR[1]-VAR0[1])/(rkdt*CFACTOR);
			    tn2o5[ijk] =    (VAR[2]-VAR0[2])/(rkdt*CFACTOR);
			    thald[ijk] =    (VAR[3]-VAR0[3])/(rkdt*CFACTOR);
			    tco[ijk] =      (VAR[4]-VAR0[4])/(rkdt*CFACTOR);
			    thcho[ijk] =    (VAR[5]-VAR0[5])/(rkdt*CFACTOR);
			    tisopooh[ijk] = (VAR[6]-VAR0[6])/(rkdt*CFACTOR);
			    tisop[ijk] =    (VAR[7]-VAR0[7])/(rkdt*CFACTOR);
			    tmvkmacr[ijk] = (VAR[8]-VAR0[8])/(rkdt*CFACTOR);
			    txo2[ijk] =     (VAR[9]-VAR0[9])/(rkdt*CFACTOR);
			    tisopao2[ijk] = (VAR[10]-VAR0[10])/(rkdt*CFACTOR);
			    tno2[ijk] =     (VAR[11]-VAR0[11])/(rkdt*CFACTOR);
			    to3[ijk] =      (VAR[12]-VAR0[12])/(rkdt*CFACTOR);
			    tno[ijk] =      (VAR[13]-VAR0[13])/(rkdt*CFACTOR);
			    tch3o2[ijk] =   (VAR[14]-VAR0[14])/(rkdt*CFACTOR);
			    tisopbo2[ijk] = (VAR[15]-VAR0[15])/(rkdt*CFACTOR);
			    tno3[ijk] =     (VAR[16]-VAR0[16])/(rkdt*CFACTOR);
			    tho2[ijk] =     (VAR[17]-VAR0[17])/(rkdt*CFACTOR);
		    }
		    else
		    {
			    nderiv += 1;
			    th2o2[ijk] =    deriv[0]/CFACTOR;
   			    tch4[ijk] =     deriv[1]/CFACTOR;
			    tn2o5[ijk] =    deriv[2]/CFACTOR;
			    thald[ijk] =    deriv[3]/CFACTOR;
			    tco[ijk] =      deriv[4]/CFACTOR;
			    thcho[ijk] =    deriv[5]/CFACTOR;
			    tisopooh[ijk] = deriv[6]/CFACTOR;
			    tisop[ijk] =    deriv[7]/CFACTOR;
			    tmvkmacr[ijk] = deriv[8]/CFACTOR;
			    txo2[ijk] =     deriv[9]/CFACTOR;
			    tisopao2[ijk] = deriv[10]/CFACTOR;
			    tno2[ijk] =     deriv[11]/CFACTOR;
			    to3[ijk] =      deriv[12]/CFACTOR;
			    tno[ijk] =      deriv[13]/CFACTOR;
			    tch3o2[ijk] =   deriv[14]/CFACTOR;
			    tisopbo2[ijk] = deriv[15]/CFACTOR;
			    tno3[ijk] =     deriv[16]/CFACTOR;
			    tho2[ijk] =     deriv[17]/CFACTOR;
		    }
		    tscale[0] = std::min(tscale[0],h2o2[ijk]/ABS(th2o2[ijk]));
		    tscale[1] = std::min(tscale[1],ch4[ijk]/ABS(tch4[ijk]));
		    tscale[2] = std::min(tscale[2],n2o5[ijk]/ABS(tn2o5[ijk]));
		    tscale[3] = std::min(tscale[3],hald[ijk]/ABS(thald[ijk]));
		    tscale[4] = std::min(tscale[4],co[ijk]/ABS(tco[ijk]));
		    tscale[5] = std::min(tscale[5],hcho[ijk]/ABS(thcho[ijk]));
		    tscale[6] = std::min(tscale[6],isopooh[ijk]/ABS(tisopooh[ijk]));
		    tscale[7] = std::min(tscale[7],isop[ijk]/ABS(tisop[ijk]));
		    tscale[8] = std::min(tscale[8],mvkmacr[ijk]/ABS(tmvkmacr[ijk]));
		    tscale[9] = std::min(tscale[9],xo2[ijk]/ABS(txo2[ijk]));
		    tscale[10] = std::min(tscale[10],isopao2[ijk]/ABS(tisopao2[ijk]));
		    tscale[11] = std::min(tscale[11],no2[ijk]/ABS(tno2[ijk]));
		    tscale[12] = std::min(tscale[12],o3[ijk]/ABS(to3[ijk]));
		    tscale[13] = std::min(tscale[13],no[ijk]/ABS(tno[ijk]));
		    tscale[14] = std::min(tscale[14],ch3o2[ijk]/ABS(tch3o2[ijk]));
		    tscale[15] = std::min(tscale[15],isopbo2[ijk]/ABS(tisopbo2[ijk]));
		    tscale[16] = std::min(tscale[16],no3[ijk]/ABS(tno3[ijk]));
		    tscale[17] = std::min(tscale[17],ho2[ijk]/ABS(tho2[ijk]));
                }
	}
        printf ("N2O5:  %12.3e ", tscale[2] );
        printf ("NO3:  %12.3e ", tscale[16] );
        printf (" %i  %i  ", nkpp, nderiv);
        printf ("\n");
    }

}

template<typename TF>
Chemistry<TF>::Chemistry(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin)
{
    const std::string group_name = "default";
    auto& gd = grid.get_grid_data();
    fields.init_diagnostic_field("oh","oh","ppb", group_name, gd.sloc);
}

template <typename TF>
Chemistry<TF>::~Chemistry()
{
}

template <typename TF>
void Chemistry<TF>::create_stats(Stats<TF>& stats)
{
   const std::string group_name = "default";

   const std::vector<std::string> stat_op_def = {"mean", "2", "3", "4", "w", "grad", "diff", "flux", "path"};
   const std::vector<std::string> stat_op_w = {"mean", "2", "3", "4"};
   const std::vector<std::string> stat_op_p = {"mean", "2", "w", "grad"};

    // Add the profiles to te statistics
   if (stats.get_switch())
   {
       stats.add_profs(*fields.sd.at("oh"), "z", stat_op_w, group_name);

   }
}
template<typename TF>
void Chemistry<TF>::exec_stats(Stats<TF>& stats)
{
    const TF no_offset = 0.;
    const TF no_threshold = 0.;

    stats.calc_stats("oh", *fields.sd.at("oh"), no_offset, no_threshold);
}

template <typename TF>
void Chemistry<TF>::init(Input& inputin)
{
    for (auto& it : fields.st)
    {
        const std::string type = inputin.get_item<std::string>("chemistry", "swchemistry", it.first, "0");
        if (type == "0")
        {
            // Cycle to avoid reading unneeded namelist options.
            continue;
        }
        else if (type == "enabled")
        {
            cmap[it.first].type = Chemistry_type::enabled;
        }
        else if (type == "disabled")
        {
            cmap[it.first].type = Chemistry_type::disabled;
        }
        else
            throw std::runtime_error("Invalid option for \"Chemistry type\"");
    }
    switch_dt = inputin.get_item<TF>("chemistry", "switch_dt","", (TF)1e5);

}

template <typename TF>
void Chemistry<TF>::create(Input& inputin, Stats<TF>& stats, Thermo<TF>&)
{
}

#ifndef USECUDA
template <typename TF>
void Chemistry<TF>::exec(Thermo<TF>& thermo,double sdt,double dt)
{
    auto& gd = grid.get_grid_data();
    
    auto Temp = fields.get_tmp();
    thermo.get_thermo_field(*Temp, "T", false, false);
    
    // determine sub time step:
    TF rkdt = 0.0;
    if (abs(sdt/dt - 1./3.) < 1e-5)
	    rkdt = dt*(TF)1./(TF)3.;
    else if (abs(sdt/dt - 15./16.) < 1e-5)
	    rkdt = dt*(TF)5./(TF)12.;
    else if (abs(sdt/dt - 8./15.) < 1e-5)
	    rkdt = dt*(TF)1./(TF)4.;
    else
	    throw std::runtime_error("Invalid time step in RK3");
 
    pss<TF>(
	    fields.st.at("ch4")    ->fld.data(), fields.sp.at("ch4")->fld.data(), 
	    fields.st.at("h2o2")   ->fld.data(), fields.sp.at("h2o2")->fld.data(), 
	    fields.st.at("n2o5")   ->fld.data(), fields.sp.at("n2o5")->fld.data(), 
	    fields.st.at("hald")   ->fld.data(), fields.sp.at("hald")->fld.data(), 
	    fields.st.at("co")     ->fld.data(), fields.sp.at("co")->fld.data(), 
	    fields.st.at("hcho")   ->fld.data(), fields.sp.at("hcho")->fld.data(), 
	    fields.st.at("isopooh")->fld.data(), fields.sp.at("isopooh")->fld.data(), 
	    fields.st.at("isop")   ->fld.data(), fields.sp.at("isop")->fld.data(), 
	    fields.st.at("xo2")    ->fld.data(), fields.sp.at("xo2")->fld.data(), 
	    fields.st.at("mvkmacr")->fld.data(), fields.sp.at("mvkmacr")->fld.data(), 
	    fields.st.at("isopao2")->fld.data(), fields.sp.at("isopao2")->fld.data(), 
	    fields.st.at("no2")    ->fld.data(), fields.sp.at("no2")->fld.data(), 
	    fields.st.at("ho2")    ->fld.data(), fields.sp.at("ho2")->fld.data(), 
	    fields.st.at("no")     ->fld.data(), fields.sp.at("no")->fld.data(), 
	    fields.st.at("isopbo2")->fld.data(), fields.sp.at("isopbo2")->fld.data(), 
	    fields.st.at("ch3o2")  ->fld.data(), fields.sp.at("ch3o2")->fld.data(), 
	    fields.st.at("no3")    ->fld.data(), fields.sp.at("no3")->fld.data(), 
	    fields.st.at("o3")     ->fld.data(), fields.sp.at("o3")->fld.data(), 
	    fields.sd.at("oh")     ->fld.data(),
	    fields.sp.at("qt")     ->fld.data(),
	    Temp ->fld.data(), rkdt, switch_dt,
	    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
	    gd.icells, gd.ijcells, gd.dz.data(), fields.rhoref.data());
    fields.release_tmp(Temp);
}
#endif

template class Chemistry<double>;
//:template class Chemistry<float>;
