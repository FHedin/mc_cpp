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
    id = 0;
    symbol = res_label = seg_label = "";
    x = y = z = 0.0;
    charge = 0.0;
    epsilon = sigma = 0.0;
    residue_id_global = residue_id_seg = 0;
}

Atom::~Atom()
{
}

/** Get the unique id and type of this atom **/
int Atom::getID() const
{
    return id;
}

void Atom::setId(int id)
{
    this->id = id;
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

void Atom::setSigma(double _sigma)
{
    this->sigma = _sigma;
}

double Atom::getSigma() const
{
    return sigma;
}

void Atom::setEpsilon(double _epsilon)
{
    this->epsilon = _epsilon;
}

double Atom::getEpsilon() const
{
    return epsilon;
}

void Atom::setSymbol(const char symbol[])
{
    this->symbol = symbol;
}

std::string Atom::getSymbol() const
{
    return symbol;
}

void Atom::setMass(double mass)
{
    this->mass = mass;
}

double Atom::getMass() const
{
    return mass;
}

//int Atom::getId() const
//{
//    return id;
//}

void Atom::setIs_frozen(bool is_frozen)
{
    this->is_frozen = is_frozen;
}

bool Atom::isIs_frozen() const
{
    return is_frozen;
}

void Atom::setType(int type)
{
    this->type = type;
}

int Atom::getType() const
{
    return type;
}

void Atom::getCentreOfMass(std::vector<Atom>& at_List, double cmass[3], int n)
{
    double crd[3];
    cmass[0] = cmass[1] = cmass[2] = 0.0;

    for (std::vector<Atom>::iterator it = at_List.begin(); it != at_List.end(); ++it)
    {
        it->getCoords(crd);
        cmass[0] += crd[0];
        cmass[1] += crd[1];
        cmass[2] += crd[2];
    }
    cmass[0] /= n;
    cmass[1] /= n;
    cmass[2] /= n;
}

void Atom::setResidue_id_seg(int residue_id_seg)
{
    this->residue_id_seg = residue_id_seg;
}

int Atom::getResidue_id_seg() const
{
    return residue_id_seg;
}

void Atom::setResidue_id_global(int residue_id_global)
{
    this->residue_id_global = residue_id_global;
}

int Atom::getResidue_id_global() const
{
    return residue_id_global;
}

void Atom::setSeg_label(const char seg_label[])
{
    this->seg_label = seg_label;
}

std::string Atom::getSeg_label() const
{
    return seg_label;
}

void Atom::setRes_label(const char res_label[])
{
    this->res_label = res_label;
}

std::string Atom::getRes_label() const
{
    return res_label;
}

void Atom::toString()
{
    std::cout << id << '\t' << residue_id_global << '\t' << res_label << '\t';
    std::cout << symbol << '\t' << x << '\t' << y << '\t' << z << '\t' << seg_label << '\t' << residue_id_seg << '\t';
    std::cout << std::endl;
}