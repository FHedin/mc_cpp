/*
 *  mc_cpp : A Molecular Monte Carlo simulations software.
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

#ifndef MC_SPAV_H
#define	MC_SPAV_H

#ifdef SPAV_EXPERIMENTAL

#include "Global_include.hpp"

#include "MC.hpp"

class MC_spav : public MC
{
public:
    MC_spav(std::vector<Atom>& _at_List, PerConditions& _pbc,
            Ensemble& _ens, FField& _ff, int _steps = 1000000,
            double _dmax = 0.25, int _update_frequency = 1000,
            double _we = 0.5, int _me = 5, int _ne = 5);
    virtual ~MC_spav();

    //overriding class MC
    virtual void run();

private:
    //normal distribution
    std::normal_distribution<double> distributionNormal;

    // spatial averaging parameters
    double we;
    int me;
    int ne;
    std::vector<Atom> oldAtList, newAtList;
    std::vector<double> oldEList, newEList;
    //            std::vector<double> Sold , Snew;
    std::vector<double> deltaM;

    // spatial averaging methods
    void buildSpavConfigs(Atom const& oldAt, Atom const& newAt);

    //overriding class MC
    virtual void apply_criterion(Atom const& oldAt, Atom const& newAt, int candidate); //one atom
    virtual void apply_criterion(std::vector<Atom>& candidateVector); //all atoms

    double rndNorm();
    void rndNorm(double _crd[3]);
};

#endif //SPAV_EXPERIMENTAL

#endif	// MC_SPAV_H

