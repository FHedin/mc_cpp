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
#include <random>
#include <cstdint>

#include "Atom.h"
#include "Ensemble.h"
#include "PerConditions.h"
#include "FField.h"

class MC
{
    public:
        MC(std::vector<Atom>& _at_List, PerConditions& _pbc, Ensemble& _ens, FField& _ff);
        virtual ~MC();

        virtual void run()=0;

    protected:
        // the random number generator stuff    
        std::random_device seed;
        std::mt19937 generator;
        std::uniform_real_distribution<double> distributionAlpha;
        std::uniform_real_distribution<double> distributionMove;
        void rndInit();
        void rndInit(uint64_t _seed);
        double rndUnifMove(double scale=1.0);
        double rndUnifAlpha();
        int    rndCandidate(int _nat);
//        int rndCandidate(int candidate);
//        double rndNorm();
        
        std::vector<Atom>& at_List;
        PerConditions& pbc;
        Ensemble& ens;
        FField& ff;

        int nsteps;
        double dmax;
        bool isAccepted;
        int upFreq;

        FILE *xyz;
        FILE *efile;
        FILE *pfile; 
        
        //assign initial random coordinates for atoms
        void Init();
        
        void move(Atom& newAt); //one atom
        void move(std::vector<Atom>& candidateVector); //all atoms
        
        virtual void apply_criterion(Atom const& oldAt, Atom const& newAt, int candidate)=0; //one atom
        virtual void apply_criterion(std::vector<Atom>& candidateVector)=0; //all atoms
        
        void adj_dmax(double acc, double each);
        
        void write_traj() const;
        
        void recentre();
};

#endif // MC_H
