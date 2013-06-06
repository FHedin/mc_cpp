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


Atom::Atom(int _id, std::string _symbol)
{
    id = _id;
    symbol = _symbol;
    x = y = z = 0.0;
    charge = 0.0;
    epsilon = sigma = 1.0;
}

Atom::Atom()
{
}

Atom::~Atom()
{
}

/** Get the unique id and type of this atom **/
int Atom::getID() const
{
    return id;
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

void Atom::addX(double _x)
{
    x += _x;
}

void Atom::addY(double _y)
{
    y += _y; 
}

void Atom::addZ(double _z)
{
    z += _z;
}

void Atom::addCoords(double _x, double _y, double _z)
{
    x += _x;
    y += _y;
    z += _z;
}

void Atom::addCoords(double _crd[3])
{
    x += _crd[0];
    y += _crd[1];
    z += _crd[2];
}
    

/** Mass and charge manips **/
void Atom::setCharge(double _charge)
{
    charge = _charge;
}

//void Atom::setMass(double _mass)
//{
//    mass = _mass;
//}
//
//double Atom::getMass() const
//{
//    return mass;
//}

double Atom::getCharge() const
{
    return charge;
}

void Atom::setSigma(double _sigma) {
    this->sigma = _sigma;
}

double Atom::getSigma() const {
    return sigma;
}

void Atom::setEpsilon(double _epsilon) {
    this->epsilon = _epsilon;
}

double Atom::getEpsilon() const {
    return epsilon;
}

std::string Atom::getSymbol() const {
    return symbol;
}

int Atom::getId() const {
    return id;
}

void Atom::getCentreOfMass(std::vector<Atom>& at_List, double cmass[3], int n)
{
    double crd[3];
    cmass[0]=cmass[1]=cmass[2]=0.0;
    
    for(std::vector<Atom>::iterator it = at_List.begin() ; it != at_List.end() ; ++it)
    {
        it->getCoords(crd);
        cmass[0]+=crd[0];
        cmass[1]+=crd[1];
        cmass[2]+=crd[2];
    }
    cmass[0]/=n;
    cmass[1]/=n;
    cmass[2]/=n;
}

void Atom::toString()
{
    std::cout << "Symbol " << symbol << " | id " << id << " | charge " << charge ;
    std::cout << " | epsilon " << epsilon << " | sigma " << sigma << std::endl;
}