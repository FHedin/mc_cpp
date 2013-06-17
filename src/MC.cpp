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

#include <cstdio>
#include <iostream>
#include <functional>

#include "MC.h"

MC::MC(std::vector<Atom>& _at_List, PerConditions& _pbc, Ensemble& _ens, FField& _ff) : at_List(_at_List), pbc(_pbc), ens(_ens), ff(_ff)
{
    //    rndInit(1566636691);
    rndInit();
}

MC::~MC()
{
}

void MC::Init()
{
    double crd[3];
    double pbv[3];

    switch (pbc.getType())
    {
        case NONE:
            pbv[0] = ens.getN() / 8.0;
            pbv[1] = ens.getN() / 8.0;
            pbv[2] = ens.getN() / 8.0;
            break;
        default:
            pbc.get_pbc_vectors(pbv);
            break;
    }

    for (auto& it : at_List)
    {
        crd[0] = pbv[0] * rndUnifMove();
        crd[1] = pbv[1] * rndUnifMove();
        crd[2] = pbv[2] * rndUnifMove();

        it.setCoords(crd);
        //        pbc.applyPBC(*it); 
    }

    recentre();

    for (auto& it : at_List)
        pbc.applyPBC(it);

}

void MC::move(Atom& newAt)
{

    double initial[3] = {0.}, trial[3] = {0.};

    newAt.getCoords(initial);

    trial[0] = initial[0] + rndUnifMove(dmax);
    trial[1] = initial[1] + rndUnifMove(dmax);
    trial[2] = initial[2] + rndUnifMove(dmax);

    newAt.setCoords(trial);

    pbc.applyPBC(newAt);
}

//random move for all atoms or a list of atoms

void MC::move(std::vector<Atom>& candidateVector)
{

    double initial[3] = {0.}, trial[3] = {0.};

    for (auto& it : candidateVector)
    {
        it.getCoords(initial);

        trial[0] = initial[0] + rndUnifMove(dmax);
        trial[1] = initial[1] + rndUnifMove(dmax);
        trial[2] = initial[2] + rndUnifMove(dmax);

        it.setCoords(trial);

        pbc.applyPBC(it);
    }
}

void MC::write_traj() const
{
    double crd[3] = {0.};
    int n = ens.getN();
    const char *symb = nullptr;

    fprintf(xyz, "%d\n", n);
    fprintf(xyz, "\n");

    for (auto& it : at_List)
    {
        it.getCoords(crd);
        symb = it.getSymbol().c_str();
        fprintf(xyz, "%s\t%10.5lf\t%10.5lf\t%10.5lf\n", symb, crd[0], crd[1], crd[2]);
    }
}

void MC::adj_dmax(double acc, double each)
{
    //    std::cout << "dmax update : " << dmax << " --> "; 
    (acc / each) <= 0.5 ? dmax *= 0.95 : dmax *= 1.05;
    //    std::cout << dmax << " : targeting acceptance of 50 % " << std::endl;

    //    double pbv[3];
    //    pbc.get_pbc_vectors(pbv);

    //    dmax > (pbv[0]/2.0 - 1.0) ? dmax = pbv[0]/2.0 - 1.0  : dmax;
    //    dmax < 0.05 ? dmax = 0.05 : dmax;
}

void MC::recentre()
{
    double crd[3];
    double cmass[3] = {0., 0., 0.};

    Atom::getCentreOfMass(at_List, cmass, ens.getN());

    for (auto& it : at_List)
    {
        it.getCoords(crd);
        crd[0] -= cmass[0];
        crd[1] -= cmass[1];
        crd[2] -= cmass[2];

        it.setCoords(crd);
    }
}

void MC::rndInit()
{
    unsigned int lseed = seed();
    std::cout << "SEED : " << lseed << std::endl;
    generator.seed(lseed);
    distributionAlpha = std::uniform_real_distribution<double>(0.0, 1.0);
    distributionMove = std::uniform_real_distribution<double>(-1.0, 1.0);
}

void MC::rndInit(uint64_t _seed)
{
    std::cout << "SEED : " << _seed << std::endl;
    generator.seed(_seed);
    distributionAlpha = std::uniform_real_distribution<double>(0.0, 1.0);
    distributionMove = std::uniform_real_distribution<double>(-1.0, 1.0);
}

double MC::rndUnifMove(double scale)
{
    return scale * distributionMove(generator);
}

double MC::rndUnifAlpha()
{
    return distributionAlpha(generator);
}

int MC::rndCandidate(int _nat)
{
    return  _nat * (int) (rndUnifAlpha());
}


