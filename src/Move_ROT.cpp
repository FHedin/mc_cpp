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
#include <cmath>

#include "Move_ROT.h"
#include "Tools.h"
#include "FField.h"

using namespace std;

Move_ROT::Move_ROT()
{
}

Move_ROT::~Move_ROT()
{
}

void Move_ROT::rotate_set(vector<Atom>& at_List, PerConditions& pbc, int moveAtomList[],
                          int Pivot, double rndAngle, double rndvec[3])
{
    double pivcrd[3];

    double rotmat[3][3] = {
        {1., 0., 0.},
        {0., 1., 0.},
        {0., 0., 1.}
    };

    find_rot_mat(rotmat, rndvec, rndAngle);

    if ( Pivot >= 0 )
    {
        at_List[Pivot].getCoords(pivcrd);
    }
    else
    {
        Tools::getCentreOfMass(at_List, moveAtomList, pivcrd);
    }

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
            apply_rotation(at_List, pbc, pivcrd, rotmat, iaf, ial);
        }
        endng = nn + 2;
    }

}

void Move_ROT::find_rot_mat(double mat[3][3], double vec[3], double angle)
{
    double angleRad = angle * FField::PI / 180.0;

    double rnormfact = 1.0 / sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    vec[0] *= rnormfact;
    vec[1] *= rnormfact;
    vec[2] *= rnormfact;

    double T[3][3];
    double T2[3][3] = {
        {0., 0., 0.},
        {0., 0., 0.},
        {0., 0., 0.}
    };

    T[0][0] = T[1][1] = T[2][2] = 0.0;
    T[2][1] = vec[0];
    T[1][2] = -vec[0];
    T[0][2] = vec[1];
    T[2][0] = -vec[1];
    T[1][0] = vec[3];
    T[0][1] = -vec[3];

    double v1 = sin(angleRad);
    double v2 = 1.0 - cos(angleRad);

    for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
            mat[i][j] += T[i][j] * v1;

    for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
            for ( int k = 0; k < 3; k++ )
                T2[i][j] += T[i][k] * T[k][j];

    for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
            mat[i][j] += T2[i][j] * v2;

}

void Move_ROT::apply_rotation(vector<Atom>& at_List, PerConditions& pbc, double pivot[3], double rotmat[3][3],
                              int first, int last)
{
    double x, y, z;
    double atcrd[3];

    for ( int i = first; i <= last; i++ )
    {
        at_List[i].getCoords(atcrd);

        x = atcrd[0] - pivot[0];
        y = atcrd[1] - pivot[1];
        z = atcrd[2] - pivot[2];

        atcrd[0] = pivot[0] + rotmat[0][0] * x + rotmat[1][0] * y + rotmat[2][0] * z;
        atcrd[1] = pivot[1] + rotmat[0][1] * x + rotmat[1][1] * y + rotmat[2][1] * z;
        atcrd[2] = pivot[2] + rotmat[0][2] * x + rotmat[1][2] * y + rotmat[2][2] * z;

        at_List[i].setCoords(atcrd);
        pbc.applyPBC(at_List[i]);
    }
}
