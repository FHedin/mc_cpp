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

#include "Atom.h"

Atom::Atom()
{
}

Atom::Atom(int _id, int _type, double _mass, double _charge)
{
    id = _id;
    type = _type;
    mass = _mass;
    charge = _charge;
}

Atom::~Atom()
{
}

/** Get the unique id and type of this atom **/
int Atom::getID() const
{
    return id;
}

int  Atom::getType() const
{
    return type;
}

/** Coordinates manipulation methods **/
void Atom::setCoords(double _x, double _y, double _z)
{
    x = _x;
    y = _y;
    z = _z;
}

void Atom::setCoords(double _crd[3])
{
    x = _crd[0];
    y = _crd[1];
    z = _crd[2];
}

void Atom::setX(double _x)
{
    x = _x;
}

void Atom::setY(double _y)
{
    y = _y;
}

void Atom::setZ(double _z)
{
    z = _z;
}

void Atom::getCoords(double _crd[3]) const
{
    _crd[0] = x;
    _crd[1] = y;
    _crd[2] = z;
}

double Atom::getX() const
{
    return x;
}

double Atom::getY() const
{
    return y;
}

double Atom::getZ() const
{
    return z;
}

/** Mass and charge manips **/
void Atom::setCharge(double _charge)
{
    charge = _charge;
}

void Atom::setMass(double _mass)
{
    mass = _mass;
}

double Atom::getMass() const
{
    return mass;
}

double Atom::getCharge() const
{
    return charge;
}


