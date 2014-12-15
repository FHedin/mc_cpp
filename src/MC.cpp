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

#include <cstdlib>
#include <iostream>
#include <functional>

#include "MC.hpp"
#include "Tools.hpp"
#include "FField.hpp"

using namespace std;

MC::MC(vector<Atom>& _at_List, PerConditions& _pbc, Ensemble& _ens, FField& _ff, List_Moves& _mvlist)
: at_List(_at_List), pbc(_pbc), ens(_ens), ff(_ff), mvlist(_mvlist)
{
    //    rndInit(1566636691);
    rndInit();
    
    xyz = nullptr;
    xyz = fopen("tr.xyz", "w");

    efile = nullptr;
    efile = fopen("ener.dat", "w");
    
}

MC::~MC()
{
    fclose(xyz);
    fclose(efile);
}

void MC::write_traj(int st) const
{
    double crd[3] = {0.};
    int n = ens.getN();
    const char *symb = nullptr;

    fprintf(xyz, "%d\n", n);
    fprintf(xyz, "#At step %d\n", st);

    for ( auto& it : at_List )
    {
        it.getCoords(crd);
        symb = it.getSymbol().c_str();
        fprintf(xyz, "%s\t%10.5lf\t%10.5lf\t%10.5lf\n", symb, crd[0], crd[1], crd[2]);
    }
}

void MC::rndInit()
{
    unsigned int lseed = seed();
    cout << "SEED : " << lseed << endl;
    generator.seed(lseed);
    distributionAlpha = uniform_real_distribution<double>(0.0, 1.0);
    distributionMove = uniform_real_distribution<double>(-1.0, 1.0);
}

void MC::rndInit(uint64_t _seed)
{
    cout << "SEED : " << _seed << endl;
    generator.seed(_seed);
    distributionAlpha = uniform_real_distribution<double>(0.0, 1.0);
    distributionMove = uniform_real_distribution<double>(-1.0, 1.0);
}

double MC::rndUnifMove()
{
    return distributionMove(generator);
}

double MC::rndUnifAlpha()
{
    return distributionAlpha(generator);
}

int MC::rndIntCandidate(int _n)
{
    return (int) _n * (rndUnifAlpha());
}

void MC::rndSphere(double rnd[3])
{
    double rx, ry, rz;

    do
    {
        rx = rndUnifMove();
        ry = rndUnifMove();
        rz = rndUnifMove();
    }
    while ( (rx * rx + ry * ry + rz * rz) > 1.0 );

    rnd[0] = rx;
    rnd[1] = ry;
    rnd[2] = rz;
}

void MC::scaleVec(double r[3], double dmax)
{
    r[0] *= dmax;
    r[1] *= dmax;
    r[2] *= dmax;
}

void MC::adjust_dmax(int acc, int currentStep)
{   
    if(currentStep!=0 && currentStep%each==0)
    {
        cout << "dmax adjusted at step " << currentStep  ;
        double ratio = (double)acc/(double)each;
        cout << " ratio is " << ratio << " target is "<< target  << " \% ; dmax " << dmax ;
        
        if(ratio > target/100)
            dmax *= 1.10;
        else
            dmax *= 0.9;
        
        cout << " --> " << dmax << endl;
    }
}


