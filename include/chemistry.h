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

#ifndef CHEMISTRY_H
#define CHEMISTRY_H

#include <vector>
#include <string>
#include <map>

class Master;
class Input;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Stats;

/**
 * Class that creates a chemistry term for scalars.
 */

enum class Chemistry_type {disabled, enabled, simple};

template<typename TF>
class Chemistry
{
    public:
        Chemistry(Master&, Grid<TF>&, Fields<TF>&, Input&); ///< Constructor of the chemistry class.
        ~Chemistry();                                       ///< Destructor  of the chemistry class.

        void init(Input&);                 ///< Initialize the arrays that contain the profiles.
        void create(Input&, Stats<TF>&);   ///< Read the profiles of the forces from the input.
        void exec(double,double);     ///< Add the tendencies belonging to the chemistry processes.

    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;

	// internal variable
	struct Chemistry_var
	{
	     Chemistry_type type;
	};

        typedef std::map<std::string, Chemistry_var> Chemistry_map;
        Chemistry_map cmap;

        const std::string tend_name = "chemistry";
        const std::string tend_longname = "Chemistry";

};
#endif
