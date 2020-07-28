
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
	    const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
	const TF rko3no = 4.75E-4;
	const TF jno2   = 8.89E-3;
	TF rate = 0.0;

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
		    rate = (rko3no*o3[ijk]*no[ijk] - jno2*no2[ijk]);
                    to3[ijk] -= rate;
                    tno[ijk] -= rate;
                    tno2[ijk] += rate;
                }
    }

}
template<typename TF>
Chemistry<TF>::Chemistry(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin)
{
}

template <typename TF>
Chemistry<TF>::~Chemistry()
{
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
	    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
	    gd.icells, gd.ijcells);
}
#endif

template class Chemistry<double>;
template class Chemistry<float>;
