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

#ifndef DIHEDRAL_H
#define	DIHEDRAL_H

#include <iostream>

enum DIHE_TYPE{
  DCOS  = 1,
  DHARM = 2,
};

class Dihedral
{
    friend std::ostream& operator<<(std::ostream& overloadStream, const Dihedral& dihe);
    
public:
    Dihedral();
    Dihedral(int _a1, int _a2, int _a3, int _a4,
             int _typ, int _ord, double _k, double _phi0, double _mult);
    
    virtual ~Dihedral();
    double getMult() const;
    double getPhi0() const;
    double getK() const;
    int getOrder() const;
    int getType() const;
    int getAt4() const;
    int getAt3() const;
    int getAt2() const;
    int getAt1() const;
    
protected:
    int at1, at2, at3, at4;
    int type;
    int order;
    
    double k;
    double phi0;
    double mult;
    
    virtual void toString(std::ostream& stream) const;
};

#endif	/* DIHEDRAL_H */

