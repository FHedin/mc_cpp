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

#ifndef PERCONDITIONS_H
#define PERCONDITIONS_H

#include "Atom.h"

enum pbcond {NONE=0,CUBIC=1};

class PerConditions
{
public:
    PerConditions(pbcond _pbtype);
    ~PerConditions();

    pbcond getType() const;

    void set_pbc_vectors(double _pbx, double _pby=0.0, double _pbz=0.0);
    void set_pbc_angles(double _alpha, double _beta=0.0, double _gamma=0.0);
    void get_pbc_vectors(double pbv[3]) const;
    void get_pbc_angles(double pba[3]) const;

    double computeVol() const;
    void applyPBC(Atom& _at) const;
    void applyPBC(double& dx, double& dy, double& dz) const;

	static double rint(double x);

private:
    pbcond pbtype;
    double alpha, beta, gamma;
    double pbx, pby, pbz;
    double rpbx, rpby, rpbz; // it is 1/pbx , 1/pby ...
    
};

#endif // PERCONDITIONS_H

