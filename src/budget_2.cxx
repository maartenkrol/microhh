/*
 * MicroHH
 * Copyright (c) 2011-2019 Chiel van Heerwaarden
 * Copyright (c) 2011-2019 Thijs Heus
 * Copyright (c) 2014-2019 Bart van Stratum
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

#include <cstdio>
#include <cmath>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "finite_difference.h"
#include "fast_math.h"
#include "model.h"
#include "thermo.h"
#include "diff.h"
#include "advec.h"
#include "force.h"
#include "stats.h"

#include "budget.h"
#include "budget_2.h"

namespace
{
    using namespace Finite_difference::O2;

    /**
     * Calculate the kinetic and turbulence kinetic energy
     */
    template<typename TF>
    void calc_kinetic_energy(
            TF* const restrict ke, TF* const restrict tke,
            const TF* const restrict u, const TF* const restrict v, const TF* const restrict w,
            const TF* const restrict umodel, const TF* const restrict vmodel,
            const TF utrans, const TF vtrans,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        using Fast_math::pow2;

        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;
    
        for (int k=kstart; k<kend; ++k)
        {
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
    
                    const TF u2 = pow2(interp2(u[ijk]+utrans, u[ijk+ii]+utrans));
                    const TF v2 = pow2(interp2(v[ijk]+vtrans, v[ijk+jj]+vtrans));
                    const TF w2 = pow2(interp2(w[ijk]       , w[ijk+kk]       ));
    
                    ke[ijk] = 0.5 * (u2 + v2 + w2);
                }
    
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
    
                    const TF u2 = pow2(interp2(u[ijk]-umodel[k], u[ijk+ii]-umodel[k]));
                    const TF v2 = pow2(interp2(v[ijk]-vmodel[k], v[ijk+jj]-vmodel[k]));
                    const TF w2 = pow2(interp2(w[ijk]          , w[ijk+kk]          ));
    
                    tke[ijk] = TF(0.5) * (u2 + v2 + w2);
                }
        }
    }
}

template<typename TF>
Budget_2<TF>::Budget_2(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin,
        Thermo<TF>& thermoin, Diff<TF>& diffin, Advec<TF>& advecin, Force<TF>& forcein, Input& inputin) :
    Budget<TF>(masterin, gridin, fieldsin, thermoin, diffin, advecin, forcein, inputin),
    field3d_operators(masterin, gridin, fieldsin)
{
    // The LES flux budget requires one additional ghost cell in the horizontal
    if (diff.get_switch() == Diffusion_type::Diff_smag2)
    {
        const int igc = 2;
        const int jgc = 2;
        const int kgc = 1;

        grid.set_minimum_ghost_cells(igc, jgc, kgc);
    }
}

template<typename TF>
Budget_2<TF>::~Budget_2()
{
}

template<typename TF>
void Budget_2<TF>::init()
{
    auto& gd = grid.get_grid_data();

    umodel.resize(gd.kcells);
    vmodel.resize(gd.kcells);
}

template<typename TF>
void Budget_2<TF>::create(Stats<TF>& stats)
{
    const std::string group_name = "budget";

    // Add the profiles for the kinetic energy to the statistics.
    stats.add_prof("ke" , "Kinetic energy" , "m2 s-2", "z", group_name);
    stats.add_prof("tke", "Turbulent kinetic energy" , "m2 s-2", "z", group_name);

    /*
    // Add the profiles for the kinetic energy budget to the statistics.
    stats.add_prof("u2_shear" , "Shear production term in U2 budget" , "m2 s-3", "z" , group_name);
    stats.add_prof("v2_shear" , "Shear production term in V2 budget" , "m2 s-3", "z" , group_name);
    stats.add_prof("tke_shear", "Shear production term in TKE budget", "m2 s-3", "z" , group_name);
    stats.add_prof("uw_shear" , "Shear production term in UW budget" , "m2 s-3", "zh", group_name);
    stats.add_prof("vw_shear" , "Shear production term in VW budget" , "m2 s-3", "zh", group_name);

    stats.add_prof("u2_turb" , "Turbulent transport term in U2 budget" , "m2 s-3", "z" , group_name);
    stats.add_prof("v2_turb" , "Turbulent transport term in V2 budget" , "m2 s-3", "z" , group_name);
    stats.add_prof("w2_turb" , "Turbulent transport term in W2 budget" , "m2 s-3", "zh", group_name);
    stats.add_prof("tke_turb", "Turbulent transport term in TKE budget", "m2 s-3", "z" , group_name);
    stats.add_prof("uw_turb" , "Turbulent transport term in UW budget" , "m2 s-3", "zh", group_name);
    stats.add_prof("vw_turb" , "Turbulent transport term in VW budget" , "m2 s-3", "zh", group_name);

    if (diff.get_switch() != Diffusion_type::Disabled)
    {
        stats.add_prof("u2_diss" , "Dissipation term in U2 budget" , "m2 s-3", "z" , group_name);
        stats.add_prof("v2_diss" , "Dissipation term in V2 budget" , "m2 s-3", "z" , group_name);
        stats.add_prof("w2_diss" , "Dissipation term in W2 budget" , "m2 s-3", "zh", group_name);
        stats.add_prof("tke_diss", "Dissipation term in TKE budget", "m2 s-3", "z" , group_name);
        stats.add_prof("uw_diss" , "Dissipation term in UW budget" , "m2 s-3", "zh", group_name);
        stats.add_prof("vw_diss" , "Dissipation term in VW budget" , "m2 s-3", "zh", group_name);

        stats.add_prof("u2_visc" , "Viscous transport term in U2 budget" , "m2 s-3", "z" , group_name);
        stats.add_prof("v2_visc" , "Viscous transport term in V2 budget" , "m2 s-3", "z" , group_name);
        stats.add_prof("w2_visc" , "Viscous transport term in W2 budget" , "m2 s-3", "zh", group_name);
        stats.add_prof("tke_visc", "Viscous transport term in TKE budget", "m2 s-3", "z" , group_name);
        stats.add_prof("uw_visc" , "Viscous transport term in UW budget" , "m2 s-3", "zh", group_name);
        stats.add_prof("vw_visc" , "Viscous transport term in VW budget" , "m2 s-3", "zh", group_name);

        // For LES, add the total diffusive budget terms, which (unlike diss + visc) close
        if (diff.get_switch() == Diffusion_type::Diff_smag2)
        {
            stats.add_prof("u2_diff" , "Total diffusive term in U2 budget" , "m2 s-3", "z" , group_name);
            stats.add_prof("v2_diff" , "Total diffusive term in V2 budget" , "m2 s-3", "z" , group_name);
            stats.add_prof("w2_diff" , "Total diffusive term in W2 budget" , "m2 s-3", "zh", group_name);
            stats.add_prof("tke_diff", "Total diffusive term in TKE budget", "m2 s-3", "z" , group_name);
            stats.add_prof("uw_diff" , "Total diffusive term in UW budget" , "m2 s-3", "zh", group_name);
            stats.add_prof("vw_diff" , "Total diffusive term in VW budget" , "m2 s-3", "zh", group_name);
        }
    }

    if (force.get_switch_lspres() == Large_scale_pressure_type::Geo_wind)
    {
        stats.add_prof("u2_cor", "Coriolis term in U2 budget", "m2 s-3", "z" , group_name);
        stats.add_prof("v2_cor", "Coriolis term in V2 budget", "m2 s-3", "z" , group_name);
        stats.add_prof("uw_cor", "Coriolis term in UW budget", "m2 s-3", "zh", group_name);
        stats.add_prof("vw_cor", "Coriolis term in VW budget", "m2 s-3", "zh", group_name);
    }

    if (thermo.get_switch() != "0")
    {
        stats.add_prof("w2_buoy" , "Buoyancy production/destruction term in W2 budget" , "m2 s-3", "zh", group_name);
        stats.add_prof("tke_buoy", "Buoyancy production/destruction term in TKE budget", "m2 s-3", "z" , group_name);
        stats.add_prof("uw_buoy" , "Buoyancy production/destruction term in UW budget" , "m2 s-3", "zh", group_name);
        stats.add_prof("vw_buoy" , "Buoyancy production/destruction term in VW budget" , "m2 s-3", "zh", group_name);

        stats.add_prof("b2_shear", "Shear production term in B2 budget"   , "m2 s-5", "z", group_name);
        stats.add_prof("b2_turb" , "Turbulent transport term in B2 budget", "m2 s-5", "z", group_name);

        stats.add_prof("bw_shear", "Shear production term in B2 budget"   , "m2 s-4", "zh", group_name);
        stats.add_prof("bw_turb" , "Turbulent transport term in B2 budget", "m2 s-4", "zh", group_name);

        if (diff.get_switch() != Diffusion_type::Disabled)
        {
            stats.add_prof("b2_visc" , "Viscous transport term in B2 budget", "m2 s-5", "z" , group_name);
            stats.add_prof("b2_diss" , "Dissipation term in B2 budget"      , "m2 s-5", "z" , group_name);
            stats.add_prof("bw_visc" , "Viscous transport term in BW budget", "m2 s-4", "zh", group_name);
            stats.add_prof("bw_diss" , "Dissipation term in BW budget"      , "m2 s-4", "zh", group_name);
        }

        stats.add_prof("bw_rdstr", "Redistribution term in BW budget"     , "m2 s-4", "zh", group_name);
        stats.add_prof("bw_buoy" , "Buoyancy term in BW budget"           , "m2 s-4", "zh", group_name);
        stats.add_prof("bw_pres" , "Pressure transport term in BW budget" , "m2 s-4", "zh", group_name);
    }

    stats.add_prof("w2_pres" , "Pressure transport term in W2 budget" , "m2 s-3", "zh", group_name);
    stats.add_prof("tke_pres", "Pressure transport term in TKE budget", "m2 s-3", "z" , group_name);
    stats.add_prof("uw_pres" , "Pressure transport term in UW budget" , "m2 s-3", "zh", group_name);
    stats.add_prof("vw_pres" , "Pressure transport term in VW budget" , "m2 s-3", "zh", group_name);

    stats.add_prof("u2_rdstr", "Pressure redistribution term in U2 budget", "m2 s-3", "z" , group_name);
    stats.add_prof("v2_rdstr", "Pressure redistribution term in V2 budget", "m2 s-3", "z" , group_name);
    stats.add_prof("w2_rdstr", "Pressure redistribution term in W2 budget", "m2 s-3", "zh", group_name);
    stats.add_prof("uw_rdstr", "Pressure redistribution term in UW budget", "m2 s-3", "zh", group_name);
    stats.add_prof("vw_rdstr", "Pressure redistribution term in VW budget", "m2 s-3", "zh", group_name);
    */
}

template<typename TF>
void Budget_2<TF>::exec_stats(Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    // Calculate the mean of the fields.
    field3d_operators.calc_mean_profile(umodel.data(), fields.mp.at("u")->fld.data());
    field3d_operators.calc_mean_profile(vmodel.data(), fields.mp.at("v")->fld.data());

    // Calculate kinetic and turbulent kinetic energy
    auto ke  = fields.get_tmp();
    auto tke = fields.get_tmp();

    const TF no_offset = 0.;
    const TF no_threshold = 0.;

    calc_kinetic_energy(
            ke->fld.data(), tke->fld.data(),
            fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
            umodel.data(), vmodel.data(),
            grid.utrans, grid.vtrans,
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    stats.calc_stats("ke" , *ke , no_offset, no_threshold);
    stats.calc_stats("tke", *tke, no_offset, no_threshold);

    fields.release_tmp(ke );
    fields.release_tmp(tke);

    /*
    // Calculate the shear production and turbulent transport terms
    calc_advection_terms(m->profs["u2_shear"].data, m->profs["v2_shear"].data, m->profs["tke_shear"].data,
                         m->profs["uw_shear"].data, m->profs["vw_shear"].data,
                         m->profs["u2_turb"].data,  m->profs["v2_turb"].data,  m->profs["w2_turb"].data, m->profs["tke_turb"].data,
                         m->profs["uw_turb"].data, m->profs["vw_turb"].data,
                         fields.u->data, fields.v->data, fields.w->data, umodel, vmodel,
                         fields.atmp["tmp1"]->data, fields.atmp["tmp2"]->data, grid.dzi, grid.dzhi);

    if(diff.get_switch() != "0")
    {
        // Calculate the diffusive transport and dissipation terms
        if(diff.get_switch() == "2" || diff.get_switch() == "4")
            calc_diffusion_terms_DNS(m->profs["u2_visc"].data, m->profs["v2_visc"].data, m->profs["w2_visc"].data, m->profs["tke_visc"].data, m->profs["uw_visc"].data,
                                     m->profs["u2_diss"].data, m->profs["v2_diss"].data, m->profs["w2_diss"].data, m->profs["tke_diss"].data, m->profs["uw_diss"].data,
                                     fields.atmp["tmp1"]->data, fields.atmp["tmp2"]->data, fields.atmp["tmp3"]->data, fields.u->data, fields.v->data, fields.w->data, umodel, vmodel,
                                     grid.dzi, grid.dzhi, grid.dxi, grid.dyi, fields.visc);
        else if(diff.get_switch() == "smag2")
            calc_diffusion_terms_LES(m->profs["u2_diss"].data,  m->profs["v2_diss"].data, m->profs["w2_diss"].data,
                                     m->profs["tke_diss"].data, m->profs["uw_diss"].data, m->profs["vw_diss"].data,
                                     m->profs["u2_visc"].data,  m->profs["v2_visc"].data, m->profs["w2_visc"].data,
                                     m->profs["tke_visc"].data, m->profs["uw_visc"].data, m->profs["vw_visc"].data,
                                     m->profs["u2_diff"].data,  m->profs["v2_diff"].data, m->profs["w2_diff"].data,
                                     m->profs["tke_diff"].data, m->profs["uw_diff"].data, m->profs["vw_diff"].data,
                                     fields.atmp["tmp1"]->data, fields.atmp["tmp2"]->data, fields.atmp["tmp3"]->data,
                                     fields.u->data, fields.v->data, fields.w->data,
                                     fields.u->datafluxbot, fields.v->datafluxbot,
                                     fields.sd.at("evisc")->data, umodel, vmodel,
                                     grid.dzi, grid.dzhi, grid.dxi, grid.dyi);
    }

    if(thermo.get_switch() != "0")
    {
        // Get the buoyancy diffusivity from the thermo class
        const double diff_b = thermo.get_buoyancy_diffusivity();

        // Store the buoyancy in the tmp1 field
        thermo.get_thermo_field(fields.atmp["tmp1"], fields.atmp["tmp2"], "b", true);

        // Calculate mean fields
        grid.calc_mean(fields.atmp["tmp1"]->datamean, fields.atmp["tmp1"]->data, grid.kcells);
        grid.calc_mean(fields.sd.at("p")->datamean, fields.sd.at("p")->data, grid.kcells);

        // Calculate buoyancy terms
        calc_buoyancy_terms(m->profs["w2_buoy"].data, m->profs["tke_buoy"].data,
                            m->profs["uw_buoy"].data, m->profs["vw_buoy"].data,
                            fields.u->data, fields.v->data, fields.w->data, fields.atmp["tmp1"]->data,
                            umodel, vmodel, fields.atmp["tmp1"]->datamean);

        // Buoyancy variance and flux budgets
        calc_buoyancy_terms_scalar(m->profs["bw_buoy"].data,
                                   fields.atmp["tmp1"]->data, fields.atmp["tmp1"]->data,
                                   fields.atmp["tmp1"]->datamean, fields.atmp["tmp1"]->datamean);

        if (advec.get_switch() != "0")
            calc_advection_terms_scalar(m->profs["b2_shear"].data, m->profs["b2_turb"].data,
                                        m->profs["bw_shear"].data, m->profs["bw_turb"].data,
                                        fields.atmp["tmp1"]->data, fields.w->data, fields.atmp["tmp1"]->datamean,
                                        grid.dzi, grid.dzhi);

        if (diff.get_switch() == "2" || diff.get_switch() == "4")
            calc_diffusion_terms_scalar_DNS(m->profs["b2_visc"].data, m->profs["b2_diss"].data,
                                            m->profs["bw_visc"].data, m->profs["bw_diss"].data,
                                            fields.atmp["tmp1"]->data, fields.w->data,
                                            fields.atmp["tmp1"]->datamean,
                                            grid.dzi, grid.dzhi, grid.dxi, grid.dyi, fields.visc, diff_b);

        calc_pressure_terms_scalar(m->profs["bw_pres"].data,  m->profs["bw_rdstr"].data,
                                   fields.atmp["tmp1"]->data, fields.sd.at("p")->data,
                                   fields.atmp["tmp1"]->datamean, fields.sd.at("p")->datamean,
                                   grid.dzi, grid.dzhi);
    }

    if(force.get_switch_lspres() == "geo")
    {
        const double fc = force.get_coriolis_parameter();
        calc_coriolis_terms(m->profs["u2_cor"].data, m->profs["v2_cor"].data,
                            m->profs["uw_cor"].data, m->profs["vw_cor"].data,
                            fields.u->data, fields.v->data, fields.w->data,
                            umodel, vmodel, fc);
    }

    // Calculate the pressure transport and redistribution terms
    calc_pressure_terms(m->profs["w2_pres"].data,     m->profs["tke_pres"].data,
                        m->profs["uw_pres"].data,     m->profs["vw_pres"].data,
                        m->profs["u2_rdstr"].data,    m->profs["v2_rdstr"].data, m->profs["w2_rdstr"].data,
                        m->profs["uw_rdstr"].data,    m->profs["vw_rdstr"].data,
                        fields.u->data, fields.v->data, fields.w->data, fields.sd.at("p")->data, umodel, vmodel,
                        grid.dzi, grid.dzhi, grid.dxi, grid.dyi);
    */
}

namespace
{
    // Double linear interpolation
    inline double interp2_4(const double a, const double b, const double c, const double d)
    {
        return 0.25 * (a + b + c + d);
    }
}

template class Budget_2<double>;
template class Budget_2<float>;
