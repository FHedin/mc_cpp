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

#ifndef FFIELD_H
#define FFIELD_H

#include "Atom.h"
#include "Atom.h"
#include "Ensemble.h"
#include "PerConditions.h"

class FField
{
    public:
        FField(std::vector<Atom>& _at_List, PerConditions& _pbc, Ensemble& _ens);
        ~FField();

        void resetE();
        void resetW();

        void getLJV(bool dV);
        void getLJV(Atom const& newAt, int candidate, bool dV);
        double PressFromVirial(int each);

        double get_E();
        double get_W();
        
        double tail_energy();
        double tail_pressure();
        
        static const double epsilon;
        static const double sigma;
        static const double sigma3;
        static const double sigma6;
        static const double sigma12;
        static const double kb_ch;
        static const double kb_si;
        static const double NA;
        static const double rc;
        
    protected:

    private:
        std::vector<Atom>& at_List;
        PerConditions& pbc;
        Ensemble& ens;

        double e,w;   //energy and virial
};

#endif // FFIELD_H
