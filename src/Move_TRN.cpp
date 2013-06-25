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

#include "Move_TRN.h"

using namespace std;

Move_TRN::Move_TRN()
{
}

Move_TRN::~Move_TRN()
{
}

void Move_TRN::translate_set(vector<Atom>& at_List, PerConditions& pbc, int moveAtomList[],
                                  double rndx, double rndy, double rndz)
{
    int ng = moveAtomList[0];
    int endng = ng + 2;
    int nn;
    int iaf, ial;

    for ( int it1 = 1; it1 <= ng; it1++ )
    {
        nn = moveAtomList[it1];
        for ( int it2 = endng; it2 <= nn; it2 += 2 )
        {
            iaf = moveAtomList[it2 - 1];
            ial = moveAtomList[it2];
            tr_a_b(at_List, pbc, iaf, ial, rndx, rndy, rndz);
        }
        endng = nn + 2;
    }
}

void Move_TRN::tr_a_b(vector<Atom>& at_List, PerConditions& pbc, int atFirst, int atLast,
                      double rndx, double rndy, double rndz)
{
    for ( int i = atFirst; i <= atLast; i++ )
    {
        at_List[i].addCoords(rndx, rndy, rndz);
        pbc.applyPBC(at_List[i]);
    }
}
