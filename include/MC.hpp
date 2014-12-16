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

#include "Global_include.hpp"

#include "Atom.hpp"
#include "Ensemble.hpp"
#include "PerConditions.hpp"
#include "FField.hpp"
#include "List_Moves.hpp"

class MC
{
public:
    MC(std::vector<Atom>& _at_List, PerConditions& _pbc, Ensemble& _ens, FField& _ff, List_Moves& _mvlist,
       int _steps, int _save_freq, double _dmax_value, double _dmax_target, int _dmax_each);
    virtual ~MC();

    // has to be overriden
    virtual void run() = 0;

protected:
    // the random number generator stuff    
    std::random_device seed;
    std::mt19937 generator;
    std::uniform_real_distribution<double> distributionAlpha;
    std::uniform_real_distribution<double> distributionMove;

    void write_traj(int st) const;
    
    //random generator definition
    void rndInit();
    void rndInit(uint64_t _seed);
    double rndUnifMove();
    double rndUnifAlpha();
    int rndIntCandidate(int _n);
    void rndSphere(double rnd[3]);
    
    //rescales random displacement vector
    void scaleVec(double r[3], double dmax);
    
    //can be overriden
    virtual void adjust_dmax(int acc, int currentStep);
    virtual bool initial_checks_before_running();
    
    // has to be overriden
    virtual void apply_criterion(double de) = 0;

    std::vector<Atom>& at_List;
    PerConditions& pbc;
    Ensemble& ens;
    FField& ff;
    List_Moves& mvlist;

    int nsteps;
    
    double dmax;
    double target;
    int each;
    
    bool isAccepted;
    
    int svFreq;

    FILE *xyz;
    FILE *efile;

};

#endif // MC_H
