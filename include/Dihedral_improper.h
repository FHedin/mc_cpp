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

#ifndef DIHEDRAL_IMPROPER_H
#define	DIHEDRAL_IMPROPER_H

#include "Dihedral.h"

class Dihedral_improper : public Dihedral
{
public:
    Dihedral_improper();
    Dihedral_improper(int _a1, int _a2, int _a3, int _a4,
             int _typ, int _ord, double _k, double _phi0, double _mult);
    
    ~Dihedral_improper();
    
private:

};

#endif	/* DIHEDRAL_IMPROPER_H */

