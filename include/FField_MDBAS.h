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

#include "FField.h"

#ifndef FFIELD_MDBAS_H
#define	FFIELD_MDBAS_H

class FField_MDBAS : public FField
{
public:
    FField_MDBAS(std::vector<Atom>& _at_List, PerConditions& _pbc, Ensemble& _ens);
    virtual ~FField_MDBAS();

    virtual double getEtot();
    
    virtual double computeNonBonded_full_range(int first, int last);
    virtual double computeNonBonded14_full_range(int first, int last);


private:

    virtual void computeNonBonded_full();
    virtual void computeNonBonded14_full();
    virtual double computeEelec(const double qi, const double qj, const double rt);
    virtual double computeEvdw(const double epsi, const double epsj, const double sigi,
                               const double sigj, const double r);

    virtual void computeEbond();
    virtual void computeEang();
    virtual void computeEub();
    virtual void computeEdihe();
    virtual void computeEimpr();

};

#endif	/* FFIELD_MDBAS_H */

