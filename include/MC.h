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

#ifndef MC_H
#define MC_H

#include <vector>
#include <random>

#include <cstdio>
#include <cstdint>

#include "Global_include.h"

#include "Atom.h"
#include "Ensemble.h"
#include "PerConditions.h"
#include "FField.h"
#include "List_Moves.h"

class MC
{
public:
    MC(std::vector<Atom>& _at_List, PerConditions& _pbc, Ensemble& _ens, FField& _ff, List_Moves& _mvlist);
    virtual ~MC();

    virtual void run() = 0;

protected:
    // the random number generator stuff    
    std::random_device seed;
    std::mt19937 generator;
    std::uniform_real_distribution<double> distributionAlpha;
    std::uniform_real_distribution<double> distributionMove;

    void rndInit();
    void rndInit(uint64_t _seed);
    double rndUnifMove();
    double rndUnifAlpha();
    int rndIntCandidate(int _n);
    void rndSphere(double rnd[3]);

    void scaleVec(double r[3], double dmax);

    std::vector<Atom>& at_List;
    PerConditions& pbc;
    Ensemble& ens;
    FField& ff;
    List_Moves& mvlist;

    int nsteps;
    double dmax;
    bool isAccepted;
    int svFreq;

    FILE *xyz;
    FILE *efile;

    virtual void apply_criterion(double de) = 0;
    
    double E_moving_set(int natom, int nmvtyp, int imvtyp, int imvatm);
    void write_traj(int st) const;

};

#endif // MC_H
