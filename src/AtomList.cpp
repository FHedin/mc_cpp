/*
 * mc_cpp : A Molecular Monte Carlo simulations software.
 * Copyright (C) 2015  Florent Hédin <hedin@fhedin.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "AtomList.hpp"

AtomList::AtomList()
{

}

AtomList::~AtomList()
{

}

AtomList& AtomList::operator[](size_t index)
{
    this->arrayIndex = index;
    return *this;
}

// const AtomList& AtomList::operator[](size_t index) const
// {
//     this->arrayIndex = index;
//     return *this;
// }

AtomList& AtomList::at(size_t index)
{
    this->arrayIndex = index;
    return *this;
}

// const AtomList& AtomList::at(size_t index) const
// {
//     this->arrayIndex = index;
//     return *this;
// }

/** Get the unique id and type of this atom **/
int AtomList::getID() const
{
    return this->id[arrayIndex];
}

void AtomList::setId(int _id)
{
    this->id[arrayIndex] = _id;
}

/** Coordinates manipulation methods **/
void AtomList::setCoords(double _x, double _y, double _z)
{
    this->x[arrayIndex] = _x;
    this->y[arrayIndex] = _y;
    this->z[arrayIndex] = _z;
}

void AtomList::setCoords(double _crd[3])
{
    this->x[arrayIndex] = _crd[0];
    this->y[arrayIndex] = _crd[1];
    this->z[arrayIndex] = _crd[2];
}

void AtomList::setX(double _x)
{
    this->x[arrayIndex] = _x;
}

void AtomList::setY(double _y)
{
    this->y[arrayIndex] = _y;
}

void AtomList::setZ(double _z)
{
    this->z[arrayIndex] = _z;
}

void AtomList::getCoords(double _crd[3]) const
{
    _crd[0] = this->x[arrayIndex];
    _crd[1] = this->y[arrayIndex];
    _crd[2] = this->z[arrayIndex];
}

double AtomList::getX() const
{
    return this->x[arrayIndex];
}

double AtomList::getY() const
{
    return this->y[arrayIndex];
}

double AtomList::getZ() const
{
    return this->z[arrayIndex];
}

void AtomList::addX(double _x)
{
    this->x[arrayIndex] += _x;
}

void AtomList::addY(double _y)
{
    this->y[arrayIndex] += _y;
}

void AtomList::addZ(double _z)
{
    this->z[arrayIndex] += _z;
}

void AtomList::addCoords(double _x, double _y, double _z)
{
    this->x[arrayIndex] += _x;
    this->y[arrayIndex] += _y;
    this->z[arrayIndex] += _z;
}

void AtomList::addCoords(double _crd[3])
{
    this->x[arrayIndex] += _crd[0];
    this->y[arrayIndex] += _crd[1];
    this->z[arrayIndex] += _crd[2];
}

/** Mass and charge manips **/
void AtomList::setCharge(double _charge)
{
    this->charge[arrayIndex] = _charge;
}

void AtomList::setMass(double _mass)
{
   this->mass[arrayIndex] = _mass;
}

double AtomList::getMass() const
{
   return this->mass[arrayIndex];
}

double AtomList::getCharge() const
{
    return this->charge[arrayIndex];
}

void AtomList::setSigma(double _sigma)
{
    this->sigma[arrayIndex] = _sigma;
}

double AtomList::getSigma() const
{
    return this->sigma[arrayIndex];
}

void AtomList::setEpsilon(double _epsilon)
{
    this->epsilon[arrayIndex] = _epsilon;
}

double AtomList::getEpsilon() const
{
    return this->epsilon[arrayIndex];
}

void AtomList::setSymbol(const char symbol[])
{
    this->symbol[arrayIndex] = string(symbol);
}

string AtomList::getSymbol() const
{
    return this->symbol[arrayIndex];
}

void AtomList::setIs_frozen(bool is_frozen)
{
    this->is_frozen = is_frozen;
}

bool AtomList::Is_frozen() const
{
    return is_frozen;
}

void AtomList::setType(int type)
{
    this->type = type;
}

int AtomList::getType() const
{
    return type;
}

void AtomList::setSeg_label(string seg_label)
{
    this->seg_label = seg_label;
}

void AtomList::setRes_label(string res_label)
{
    this->res_label = res_label;
}

void AtomList::setSymbol(string symbol)
{
    this->symbol = symbol;
}

void AtomList::setBeta(double beta)
{
    this->beta = beta;
}

double AtomList::getBeta() const
{
    return beta;
}

void AtomList::setSigma14(double sigma14)
{
    this->sigma14 = sigma14;
}

double AtomList::getSigma14() const
{
    return sigma14;
}

void AtomList::setEpsilon14(double epsilon14)
{
    this->epsilon14 = epsilon14;
}

double AtomList::getEpsilon14() const
{
    return epsilon14;
}

void AtomList::setResidue_id_seg(int residue_id_seg)
{
    this->residue_id_seg = residue_id_seg;
}

int AtomList::getResidue_id_seg() const
{
    return residue_id_seg;
}

void AtomList::setResidue_id_global(int residue_id_global)
{
    this->residue_id_global = residue_id_global;
}

int AtomList::getResidue_id_global() const
{
    return residue_id_global;
}

void AtomList::setSeg_label(const char seg_label[])
{
    this->seg_label = string(seg_label);
}

string AtomList::getSeg_label() const
{
    return seg_label;
}

void AtomList::setRes_label(const char res_label[])
{
    this->res_label = string(res_label);
}

string AtomList::getRes_label() const
{
    return res_label;
}
