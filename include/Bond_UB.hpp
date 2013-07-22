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

#ifndef BOND_UB_H
#define	BOND_UB_H

#include "Global_include.hpp"

#include "Bond.hpp"

class Bond_UB : public Bond
{
    friend std::ostream& operator<<( std::ostream& overloadStream, const Bond_UB& bnd_ub );

public:
    Bond_UB();
    Bond_UB(int _a1, int _a2, int _typ, double _k, double _r);

    virtual ~Bond_UB();

private:

    virtual void toString(std::ostream& stream) const;

};

#endif	/* BOND_UB_H */

