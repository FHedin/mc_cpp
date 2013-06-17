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

#ifndef BOND_H
#define	BOND_H

#include <iostream>

enum BOND_TYPE
{
    BHARM = 0,
    BMORSE = 1,
};

class Bond
{
    friend std::ostream& operator<<( std::ostream& overloadStream, const Bond& bnd );

public:
    Bond();
    Bond(int _a1, int _a2, int _typ, double _k, double _r, double _beta);

    virtual ~Bond();
    double getBeta() const;
    double getR0() const;
    double getK() const;
    int getType() const;
    int getAt2() const;
    int getAt1() const;

protected:
    int at1, at2; // Id of atoms 1 and 2 of the bond

    int type; // type of the bond : harmonic, Morse ... forcefield dependent

    // for harmonic bond
    double k; //force constant, unit is forcefield dependent (usually kcal/mol)
    double r0; //equil. distance, unit is ff dependent (usually Angstroem)

    // for other type of bonds . i.e. Morse
    double beta;

    virtual void toString(std::ostream& stream) const;
};

#endif	/* BOND_H */

