
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
#include <cstdio>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <math.h>
#include <iomanip>
#include <utility>

#include "master.h"
#include "grid.h"
#include "fields.h"
#include "thermo.h"
#include "stats.h"
#include "netcdf_interface.h"
#include "chemistry.h"
#include "constants.h"
#include "timeloop.h"
#include "deposition.h"
#include "boundary.h"
#include "cross.h"

namespace
{
    std::pair<std::string, int> check_for_unique_time_dim(const std::map<std::string, int>& dims)
    {
        // Check for the existence of a unique time dimension.
        bool only_one_time_dim = false;
        std::string time_dim;
        int time_dim_length = 0;

        for (auto i : dims)
        {
            if (i.first.substr(0, 9) == "time_chem")
            {
                if (only_one_time_dim)
                    throw std::runtime_error("More than one time dimensions in input");
                else
                {
                    only_one_time_dim = true;
                    time_dim = i.first;
                    time_dim_length = i.second;
                }
            }
        }

        return std::make_pair(time_dim, time_dim_length);
    }


    template<typename TF>
    void pss(
            TF* restrict tnh3, const TF* const restrict nh3,
            const TF* const restrict vdnh3,
            const TF* const restrict tprof,
            const TF* const restrict qprof,
            const TF* const restrict dzi,
            const TF* const restrict rhoref,
            const TF dt,
            const TF sdt,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride)
    {


        const TF xmh2o = 18.015265;
        const TF xmh2o_i = TF(1) / xmh2o;
        const TF xmair = 28.9647;       // Molar mass of dry air  [kg kmol-1]
        const TF xmair_i = TF(1) / xmair;
        const TF Na = 6.02214086e23; // Avogadros number [molecules mol-1]
        TF C_M;
        // Update the time integration of the reaction fluxes with the full timestep on first RK3 step

        for (int k=kstart; k<kend; ++k)
        {
            C_M = TF(1e-3) * rhoref[k] * Na * xmair_i;   // molecules/cm3 for chmistry!

            // From mol/mol (units mixing ratio) to molecules/cm3
            const TF CFACTOR = C_M;

            // rate constants on horizontal average quantities, 
            const TF TEMP = tprof[k];
            const TF C_H2O = std::max(qprof[k] * xmair * C_M * xmh2o_i, TF(1));
            //printf("Temp and H2o on layer = %13.5e %13.5e %4i \n", TEMP, C_H2O, k);

          
            TF nh3_l;
            TF dep_nh3;
	
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jstride + k*kstride;
                    const int ij = i + j*jstride;

                    // kg/kg --> molH2O/molAir --*C_M--> molecules/cm3 limit to 1 molecule/cm3 to avoid error usr_HO2_HO2
                    // const TF C_H2O = std::max(qt[ijk] * xmair * C_M * xmh2o_i, TF(1));
                    // const TF TEMP = temp[ijk];

                    // Convert to molecules per cm3 and add tendenccies of other processes:
                    nh3_l  = std::max((nh3[ijk] + tnh3[ijk] * sdt) * CFACTOR, TF(0));

                    const TF sdt_cfac_i = TF(1) / (sdt * CFACTOR); 
                    if (k==kstart)
                        {
                            dep_nh3 = vdnh3[ij]   * dzi[k] * nh3_l * sdt; // [molecules cm-3]
                        }
                    else
                        {
                            dep_nh3 = 0.0;
                        }

                    //  Calculate tendency and add to the tendency of the transported tracers:
                    tnh3[ijk] -=  dep_nh3 * sdt_cfac_i; // [mol/mol s-1]

 
                } // i,j
        } // k
    }
}

template<typename TF>
Chemistry<TF>::Chemistry(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin), field3d_operators(master, grid, fields)
{
    const std::string group_name = "default";
    auto& gd = grid.get_grid_data();

    sw_chemistry = inputin.get_item<bool>("chemistry", "swchemistry", "", false);

    if (!sw_chemistry)
        return;

    deposition = std::make_shared<Deposition <TF>>(masterin, gridin, fieldsin, inputin);
}

template <typename TF>
Chemistry<TF>::~Chemistry()
{
}

template<typename TF>
void Chemistry<TF>::exec_stats(const int iteration, const double time, Stats<TF>& stats)
{
    if (!sw_chemistry or stats.get_switch())
        return;

    const TF no_offset = 0.;
    const TF no_threshold = 0.;
    auto& gd = grid.get_grid_data();

    if (iteration != 0)   // this does not make sense for first step = t=0.
    {
        // add deposition velocities to statistics:
        stats.calc_stats_2d("vdnh3"   , vdnh3,   no_offset);

        // Increment the statistics index.
        ++statistics_counter;

    }

}

template <typename TF>
void Chemistry<TF>::init(Input& inputin)
{
    if (!sw_chemistry)
        return;

    auto& gd = grid.get_grid_data();

    statistics_counter = 0;

    // initialize 2D deposition arrays:
    vdnh3.resize(gd.ijcells);

    // initialize deposition routine:
    deposition-> init(inputin);

    // fill deposition with standard values:
    std::fill(vdnh3.begin(), vdnh3.end(), deposition-> get_vd("nh3"));

    master.print_message("Deposition arrays initialized, e.g. with vdnh3 = %13.5e m/s \n", deposition-> get_vd("nh3"));
    
}

template <typename TF>
void Chemistry<TF>::create(
        const Timeloop<TF>& timeloop, std::string sim_name, Netcdf_handle& input_nc,
        Stats<TF>& stats, Cross<TF>& cross)
{
    if (!sw_chemistry)
        return;

    auto& gd = grid.get_grid_data();
    int iotime = timeloop.get_iotime();

    qprof.resize(gd.kcells);
    tprof.resize(gd.kcells);

    if (stats.get_switch())
    {
        // Stats:
        const std::string group_name = "default";
        const std::vector<std::string> stat_op_def = {"mean", "2", "3", "4", "w", "grad", "diff", "flux", "path"};
        const std::vector<std::string> stat_op_w = {"mean", "2", "3", "4"};
        const std::vector<std::string> stat_op_p = {"mean", "2", "w", "grad"};



        // add the deposition-velocity timeseries in deposition group statistics
        const std::string group_named = "deposition";

        // used in chemistry:
        stats.add_time_series("vdnh3", "NH3 deposition velocity", "m s-1", group_named);
    }

    // add cross-sections
    if (cross.get_switch())
    {
        std::vector<std::string> allowed_crossvars = {"vdnh3"};
        cross_list = cross.get_enabled_variables(allowed_crossvars);

        // `deposition->create()` only creates cross-sections.
        deposition->create(stats, cross);
    }
}

template<typename TF>
void Chemistry<TF>::exec_cross(Cross<TF>& cross, unsigned long iotime)
{
    if (!sw_chemistry)
        return;

    auto& gd = grid.get_grid_data();

    for (auto& it : cross_list)
    {
        if (it == "vdnh3")
            cross.cross_plane(vdnh3.data(), "vdnh3", iotime);
    }

    // see if to write per tile:
    deposition->exec_cross(cross, iotime);
}

template <typename TF>
void Chemistry<TF>::update_time_dependent(Timeloop<TF>& timeloop, Boundary<TF>& boundary)
{
    if (!sw_chemistry)
        return;

    deposition->update_time_dependent(
            timeloop,
            boundary,
            vdnh3.data());
}


#ifndef USECUDA
template <typename TF>
void Chemistry<TF>::exec(Thermo<TF>& thermo,double sdt,double dt)
{
    if (!sw_chemistry)
        return;

    auto& gd = grid.get_grid_data();

    auto tmp = fields.get_tmp();
    thermo.get_thermo_field(*tmp, "T", true, false);

    // Calculate the mean temperature and water vapor mixing ratio.
    field3d_operators.calc_mean_profile(tprof.data(), tmp->fld.data());
    qprof = fields.sp.at("qt")->fld_mean;
    //field3d_operators.calc_mean_profile(qprof.data(), fields.sp.at("qt")->fld.data());
 

    pss<TF>(
        fields.st.at("nh3")->fld.data(), fields.sp.at("nh3")->fld.data(),
        vdnh3.data(),
        tprof.data(),
        qprof.data(),
        gd.dzi.data(),
        fields.rhoref.data(),
        dt, sdt, 
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kstart, gd.kend,
        gd.icells, gd.ijcells);

    fields.release_tmp(tmp);

    //isop_stat<TF>(
    //      fields.st.at("isop")->fld.data(), fields.sp.at("isop")->fld.data(),
    //      gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
    //      gd.icells, gd.ijcells);
}
#endif

template class Chemistry<double>;
//:template class Chemistry<float>;
