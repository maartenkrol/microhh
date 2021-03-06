
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
#include "stats.h"
#include "chemistry.h"

namespace
{
    template<typename TF>
    void pss(
            TF* restrict to3, const TF* o3, 
            TF* restrict tno, const TF* no, 
            TF* restrict tno2, const TF* no2, 
            TF* restrict trh, const TF* rh, 
            TF* restrict tho2, const TF* ho2, 
            TF* restrict oh, 
	    const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk, const TF* const restrict dz)
    {
	const TF ko3no = 4.75E-4;
	const TF jno2   = 8.9E-3;
	const TF jo3    = 2.7E-6;
	const TF kohco  = 6E-3;
	const TF fco    = 100;   // scaled reaction rate
	const TF kho2no = 0.21; 
	const TF kho2o3  = 5E-5;
	const TF kho2ho2 = 7.25E-2;
	const TF kohno2  = 0.275;
	const TF koho3   = 1.75E-3;
	const TF kohho2  = 2.75;
	const TF co   = 100.0    ;  // Krol, 2000
	TF po3 = 0.0;
	TF pno = 0.0;
	TF pno2 = 0.0;
	TF pho2 = 0.0;
	TF prh = 0.0;
	TF poh = 0.0;
	TF vdo3 = 0.0;
	TF vdno = 0.0;
	TF vdno2 = 0.0;
	TF vdho2 = 0.0;
	TF vdrh = 0.0;
	TF vdoh = 0.0;


        for (int k=kstart; k<kend; ++k)
	{
            if (k==kstart) {
	        vdo3  = TF(0.005)/dz[k];   // 1/s
	        vdno  = TF(0.002)/dz[k];   // 1/s
	        vdno2 = TF(0.005)/dz[k];   // 1/s
	        vdrh  = TF(0.001)/dz[k];   // 1/s
	        vdho2 = TF(0.010)/dz[k];   // 1/s
	        vdoh  = TF(0.010)/dz[k];   // 1/s
	    }
            else {
	        vdo3  = TF(0.0);
	        vdno  = TF(0.0);
	        vdno2 = TF(0.0);
	        vdrh  = TF(0.0);
	        vdho2 = TF(0.0);
	        vdoh  = TF(0.0);
	    }
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
		            prh = std::max(rh[ijk],TF(0.0));
		            pno = std::max(no[ijk],TF(0.0));
		            pno2 = std::max(no2[ijk],TF(0.0));
		            pho2 = std::max(ho2[ijk],TF(0.0));
		            po3  = std::max(o3[ijk],TF(0.0));
                    //prh = rh[ijk];
                    //pno = no[ijk];
                    //pno2 = no2[ijk];
                    //pho2 = ho2[ijk];
                    //po3  = o3[ijk];
                    //poh  = oh[ijk];
                    oh[ijk] = ( (TF(2.0)*jo3 + kho2o3*pho2)*po3 + kho2no*pho2*pno)/
                         (kohco*(co+fco*prh) + kohno2*pno2 + koho3*po3 + kohho2*pho2 + vdoh);
                    //poh = ( (TF(2.0)*jo3 + kho2o3*pho2)*po3 + kho2no*pho2*pno)/
                    //    (kohco*(co+fco*prh) + kohno2*pno2 + koho3*po3 + kohho2*pho2 + vdoh);

                    poh  = std::max(oh[ijk],TF(0.0));

                    const TF fo3no = ko3no*po3*pno;
                    const TF fjno2 = jno2*pno2;
                    const TF fjo3  = jo3*po3;
                    const TF fohco = kohco*co*poh;
                    const TF fohrh = fco*kohco*poh*prh;
                    const TF fho2no = kho2no*pho2*pno;
                    const TF fho2o3 = kho2o3*pho2*po3;
                    const TF fho2ho2 = kho2ho2*pho2*pho2;
                    const TF fohno2 = kohno2*poh*pno2;
                    const TF foho3  = koho3*poh*po3;
                    const TF fohho2 = kohho2*poh*pho2;


                    to3[ijk] += fjno2 - fo3no - fjo3 - fho2o3 - foho3 - vdo3*po3;
                    tno[ijk] += fjno2 - fo3no - fho2no - vdno*pno;
                    tno2[ijk] += fo3no - fjno2 + fho2no - fohno2 - vdno2*pno2;
                    tho2[ijk] += fohco + fohrh - fho2no - fho2o3 -2*fho2ho2 - fohho2 + foho3 - vdho2*pho2;
                    trh[ijk] += -fohrh - vdrh*prh;
                }
	}
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


}

template <typename TF>
void Chemistry<TF>::create(Input& inputin, Stats<TF>& stats)
{
}

#ifndef USECUDA
template <typename TF>
void Chemistry<TF>::exec(double sdt,double dt)
{
    auto& gd = grid.get_grid_data();
    
    
    pss<TF>(
	    fields.st.at("o3")->fld.data(), fields.sp.at("o3")->fld.data(), 
	    fields.st.at("no")->fld.data(), fields.sp.at("no")->fld.data(), 
	    fields.st.at("no2")->fld.data(), fields.sp.at("no2")->fld.data(), 
	    fields.st.at("rh")->fld.data(), fields.sp.at("rh")->fld.data(), 
	    fields.st.at("ho2")->fld.data(), fields.sp.at("ho2")->fld.data(), 
	    fields.sd.at("oh")->fld.data(), 
	    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
	    gd.icells, gd.ijcells, gd.dz.data());
}
#endif

template class Chemistry<double>;
template class Chemistry<float>;
