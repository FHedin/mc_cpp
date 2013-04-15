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

#ifndef ATOM_H
#define ATOM_H

#include <vector>

class Atom
{
public:
    /** Public Methods**/
    Atom();
    Atom(int _id, int _type, double _mass=0.0, double _charge=0.0);
    ~Atom();

    // return the unique atom ID and the type
    int getID() const;
    int getType() const;

    // different ways of setting/getting coordinates
    void setCoords(double _x, double _y, double _z);
    void setCoords(double _crd[3]);
    void setX(double _x);
    void setY(double _y);
    void setZ(double _z);
    void getCoords(double _crd[3]) const;
    double getX() const;
    double getY() const;
    double getZ() const;

    // setting/getting charge and mass
    void setMass(double _mass);
    void setCharge(double _charge);
    double getMass() const;
    double getCharge() const;

    /** Public Attributes **/
    // nothing for the moment

protected:

private:
    /** Private Attributes **/
    int id;  // unique identifier for this atom
    int type;   // atom type
    double mass;
    double charge;
    double x,y,z;

//    std::vector<Atom> Bonded_List;     // List of other atoms linked to this one by a bond
//    std::vector<Atom> NBonded_List;    // List of other atoms interacting with this one via non-bonded interaction
};

#endif // ATOM_H
