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

#ifndef BOND_H
#define	BOND_H

class Bond
{
public:
    Bond(int _a1, int _a2, int _typ, double _k, double _r, double _beta);
    ~Bond();

private:
    int at1, at2; // Id of atoms 1 and 2 of the bond

    int type; // type of the bond : harmonic, Morse ... forcefield dependent

    // for harmonic bond
    int k; //force constant, unit is forcefield dependent (usually kcal/mol)
    int r0; //equil. distance, unit is ff dependent (usually Angstroem)

    // for other type of bonds . i.e. Morse
    int beta;
};

#endif	/* BOND_H */

