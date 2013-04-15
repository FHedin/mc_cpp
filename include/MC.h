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

#ifndef MC_H
#define MC_H

#include <vector>

#include "Atom.h"
#include "Ensemble.h"
#include "PerConditions.h"
#include "FField.h"

class MC
{
    public:
        MC(std::vector<Atom>& _at_List, PerConditions& _pbc, Ensemble& _ens, FField& _ff, int _steps, double _dmax);
        ~MC();

        void run();
        void write_traj() const;

    protected:

    private:
        std::vector<Atom>& at_List;
        PerConditions& pbc;
        Ensemble& ens;
        FField& ff;

        int nsteps;
        double dmax;
        bool isAccepted;

        FILE *xyz;
        FILE *efile;
        FILE *pfile;

        //assign initial random coordinates for atoms
        void randInit();
        void move(Atom& newAt) const;
        void move();
        void apply_criterion(Atom const& oldAt, Atom const& newAt, int candidate);
        void adj_dmax(double acc, double each);
};

#endif // MC_H
