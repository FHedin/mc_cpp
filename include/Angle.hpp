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

#ifndef ANGLE_H
#define	ANGLE_H

#include <iostream>

#include "Global_include.hpp"

class Angle
{
    friend std::ostream& operator<<( std::ostream& out, const Angle& angle );

public:
    Angle();
    Angle(int _a1, int _a2, int _a3, int _typ, double _k, double _theta0);

    virtual ~Angle();
    double getTheta0() const;
    double getK() const;
    int getType() const;
    int getAt3() const;
    int getAt2() const;
    int getAt1() const;

private:
    int at1, at2, at3;
    int type;

    double k;
    double theta0;

    virtual void toString(std::ostream& stream) const;

};

#endif	/* ANGLE_H */

