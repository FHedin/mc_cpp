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
#include <limits>
#include <string>
#include <algorithm>

//#include <cstdlib>
#include <cctype>
#include <cmath>

#include "PerConditions.h"
#include "Tools.h"

PerConditions::PerConditions(pbcond _pbtype, double _pbx, double _pby, double _pbz,
        double _alpha, double _beta, double _gamma) {
    switch (_pbtype) {
            //no periodic boundaries conditions
        case NONE:
            std::cout << "No PBC." << std::endl;
            pbtype = NONE;
            break;
        case CUBIC:
            std::cout << "Cubic PBC." << std::endl;
            pbtype = CUBIC;
            break;
        case ORBIC:
            std::cout << "Orthorhombic PBC." << std::endl;
            pbtype = ORBIC;
            break;
//        case TCLIN:
//            std::cout << "Triclinic PBC." << std::endl;
//            pbtype = TCLIN;
//            break;
        default:
            std::cout << "Warning : " << _pbtype << " is not a valid PBC type. Set by default to 0 (no PBC)." << std::endl;
            pbtype = NONE;
            break;
    }

    set_pbc_vectors(_pbx, _pby, _pbz);
    set_pbc_angles(_alpha, _beta, _gamma);
}

PerConditions::PerConditions(std::string str, double _pbx, double _pby, double _pbz,
        double _alpha, double _beta, double _gamma) {
    Tools::str_rm_blank_spaces(str);
    Tools::str_to_lower_case(str);

    if (!str.compare("none")) {
        std::cout << "No PBC." << std::endl;
        pbtype = NONE;
    } else if (!str.compare("cubic")) {
        std::cout << "Cubic PBC." << std::endl;
        pbtype = CUBIC;
    } else if (!str.compare("orbic")) {
        std::cout << "Orthorhombic PBC." << std::endl;
        pbtype = ORBIC;
    } /*else if (!str.compare("tclin")) {
        std::cout << "Triclinic PBC." << std::endl;
        pbtype = TCLIN;
    }*/ else {
        std::cout << "Warning : " << str << " is not a valid PBC type. Set by default to NONE (no PBC)." << std::endl;
        pbtype = NONE;
    }

    set_pbc_vectors(_pbx, _pby, _pbz);
    set_pbc_angles(_alpha, _beta, _gamma);
}

PerConditions::~PerConditions() {
}

pbcond PerConditions::getType() const {
    return pbtype;
}

void PerConditions::set_pbc_vectors(double _pbx, double _pby, double _pbz) {
    switch (pbtype) {
        case CUBIC:
            pbx = _pbx;
            pby = _pbx;
            pbz = _pbx;
            rpbx = 1.0 / pbx;
            rpby = 1.0 / pby;
            rpbz = 1.0 / pbz;
            break;
        case ORBIC:
            pbx = _pbx;
            pby = _pby;
            pbz = _pbz;
            rpbx = 1.0 / pbx;
            rpby = 1.0 / pby;
            rpbz = 1.0 / pbz;
            break;
        case TCLIN:
            pbx = _pbx;
            pby = _pby;
            pbz = _pbz;
            rpbx = 1.0 / pbx;
            rpby = 1.0 / pby;
            rpbz = 1.0 / pbz;
            break;
        default:
            break;
    }
}

void PerConditions::set_pbc_angles(double _alpha, double _beta, double _gamma) {
    switch (pbtype) {
        case CUBIC:
            alpha = _alpha;
            beta = alpha;
            gamma = alpha;
            break;
        case ORBIC:
            alpha = _alpha;
            beta = alpha;
            gamma = alpha;
            break;
        case TCLIN:
            alpha = _alpha;
            beta = _beta;
            gamma = _gamma;
            break;
        default:
            break;
    }
}

void PerConditions::get_pbc_vectors(double _pbv[3]) const {
    _pbv[0] = pbx;
    _pbv[1] = pby;
    _pbv[2] = pbz;
}

void PerConditions::get_pbc_angles(double _pba[3]) const {
    _pba[0] = alpha;
    _pba[1] = beta;
    _pba[2] = gamma;
}

/*
 Volume of triclinic system is : 
 * V = abc (1- cos2 α - cos2 β - cos2 γ) + 2(cos(α) cos(β) cos(γ))½
 */
double PerConditions::computeVol() const {
    switch (pbtype) {
        case NONE:
            return std::numeric_limits<float>::infinity();
            break;
        case CUBIC:
            return (pbx * pby * pbz);
            break;
        case ORBIC:
            return (pbx * pby * pbz);
            break;
        case TCLIN:
        {
            double cos2, cossq;
            cos2 = 1.0 - cos(alpha)*cos(alpha) - cos(beta)*cos(beta) - cos(gamma)*cos(gamma);
            cossq= sqrt( cos(alpha)*cos(beta)*cos(gamma) );
            return( (pbx * pby * pbz) * cos2 + 2.0 * cossq );
        }
            break;
        default:
            return std::numeric_limits<float>::infinity();
            break;
    }
}

void PerConditions::applyPBC(Atom& _at) const {
    switch (pbtype) {
        case CUBIC:
        {
            double tmp[3];
            _at.getCoords(tmp);
            tmp[0] -= pbx * rint(rpbx * tmp[0]);
            tmp[1] -= pby * rint(rpby * tmp[1]);
            tmp[2] -= pbz * rint(rpbz * tmp[2]);
            _at.setCoords(tmp);
            break;
        }
        case ORBIC:
        {
            double tmp[3];
            _at.getCoords(tmp);
            tmp[0] -= pbx * rint(rpbx * tmp[0]);
            tmp[1] -= pby * rint(rpby * tmp[1]);
            tmp[2] -= pbz * rint(rpbz * tmp[2]);
            _at.setCoords(tmp);
            break;
        }
        case TCLIN:
        {
            double tmp[3];
            _at.getCoords(tmp);
            tmp[0] -= pbx * rint(rpbx * tmp[0]);
            tmp[1] -= pby * rint(rpby * tmp[1]);
            tmp[2] -= pbz * rint(rpbz * tmp[2]);
            _at.setCoords(tmp);
            break;
        }
        default:
            break;
    }
}

void PerConditions::applyPBC(double tmp[3]) const {
    switch (pbtype) {
        case CUBIC:
        {
            tmp[0] -= pbx * rint(rpbx * tmp[0]);
            tmp[1] -= pby * rint(rpby * tmp[1]);
            tmp[2] -= pbz * rint(rpbz * tmp[2]);
            break;
        }
        case ORBIC:
        {
            tmp[0] -= pbx * rint(rpbx * tmp[0]);
            tmp[1] -= pby * rint(rpby * tmp[1]);
            tmp[2] -= pbz * rint(rpbz * tmp[2]);
            break;
        }
        default:
            break;
    }
}

void PerConditions::applyPBC(double& dx, double& dy, double& dz) const {
    switch (pbtype) {
        case CUBIC:
            dx -= pbx * rint(rpbx * dx);
            dy -= pby * rint(rpby * dy);
            dz -= pbz * rint(rpbz * dz);
            break;
        case ORBIC:
            dx -= pbx * rint(rpbx * dx);
            dy -= pby * rint(rpby * dy);
            dz -= pbz * rint(rpbz * dz);
            break;
        default:
            break;
    }
}

double PerConditions::rint(double x)
{
    int temp = (x >= 0. ? (int)(x + 0.5) : (int)(x - 0.5));
    return (double)temp;
}
