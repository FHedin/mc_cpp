/*
 *  mc_cpp : A basic Monte Carlo simulations software.
 *  Copyright (C) 2013  Florent Hedin
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef IO_MDBAS_H
#define	IO_MDBAS_H

#include "IO.h"

class IO_MDBAS : public IO
{
public:
    IO_MDBAS(std::string configf_name, std::string forfieldf_name,
            std::vector<Atom>& _at_List, PerConditions& _pbc, Ensemble& _ens);
    ~IO_MDBAS();

private:
    void read_coord();
    void read_ff();

    FILE* conff;
    FILE* forff;
};

#endif	/* IO_MDBAS_H */
