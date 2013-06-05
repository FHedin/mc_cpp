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

#include <iostream>
#include <cmath>
#include <limits>

#include "PerConditions.h"

PerConditions::PerConditions(pbcond _pbtype)
{
    switch (_pbtype)
    {
        //no periodic boundaries conditions
    case NONE:
        std::cout << "No PBC." << std::endl;
        pbtype = NONE;
        break;
    case CUBIC:
        std::cout << "Cubic PBC." << std::endl;
        pbtype = CUBIC;
        break;
    default:
        std::cout << "Problem : " << _pbtype << " is not a valid PBC type. Set by default to 0 (no PBC)." << std::endl;
        pbtype = NONE;
        break;
    }
}

PerConditions::~PerConditions()
{
}

pbcond PerConditions::getType() const
{
    return pbtype;
}

void PerConditions::set_pbc_vectors(double _pbx, double _pby, double _pbz)
{
    switch (pbtype)
    {
        case CUBIC:
            pbx = _pbx;
            pby = _pbx;
            pbz = _pbx;
            rpbx=1.0/pbx;
            rpby=1.0/pby;
            rpbz=1.0/pbz;
            break;
        default:
            break;
    }
}

void PerConditions::set_pbc_angles(double _alpha, double _beta, double _gamma)
{
    switch (pbtype)
    {
        case CUBIC:
            alpha = _alpha;
            beta  = alpha;
            gamma = alpha;
            break;
        default:
            break;
    }
}

void PerConditions::get_pbc_vectors(double _pbv[3]) const
{
    _pbv[0] = pbx;
    _pbv[1] = pby;
    _pbv[2] = pbz;
}

void PerConditions::get_pbc_angles(double _pba[3]) const
{
    _pba[0] = alpha;
    _pba[1] = beta;
    _pba[2] = gamma;
}

double PerConditions::computeVol() const
{
    switch(pbtype)
    {
        case NONE:
            return std::numeric_limits<float>::infinity();
            break;
        case CUBIC:
            return (pbx*pby*pbz);
            break;
        default:
            return 0;
            break;
    }
}

void PerConditions::applyPBC(Atom& _at) const
{
    switch(pbtype)
    {
        case CUBIC:
        {
            double tmp[3];
            _at.getCoords(tmp);
//            tmp[0] -= pbx*PerConditions::rint(rpbx*tmp[0]) ;
//            tmp[1] -= pby*PerConditions::rint(rpby*tmp[1]) ;
//            tmp[2] -= pbz*PerConditions::rint(rpbz*tmp[2]) ;
            tmp[0] -= pbx*rint(rpbx*tmp[0]) ;
            tmp[1] -= pby*rint(rpby*tmp[1]) ;
            tmp[2] -= pbz*rint(rpbz*tmp[2]) ;
            _at.setCoords(tmp);
            break;
        }
        default:
            break;
    }
}

void PerConditions::applyPBC(double& dx, double& dy, double& dz) const
{
    switch(pbtype)
    {
        case CUBIC:
//            dx -= pbx*PerConditions::rint(rpbx*dx) ;
//            dy -= pby*PerConditions::rint(rpby*dy) ;
//            dz -= pbz*PerConditions::rint(rpbz*dz) ;
            dx -= pbx*rint(rpbx*dx) ;
            dy -= pby*rint(rpby*dy) ;
            dz -= pbz*rint(rpbz*dz) ;
            break;
        default:
            break;
    }
}

//double PerConditions::rint(double x)
//{
//    int temp = (x >= 0. ? (int)(x + 0.5) : (int)(x - 0.5));
//    return (double)temp;
//}
