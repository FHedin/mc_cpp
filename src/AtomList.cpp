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

using namespace std;

AtomList::AtomList()
{

}

AtomList::~AtomList()
{

}

void AtomList::resize(size_t siz)
{
    id.resize(siz);
    type.resize(siz);
    charge.resize(siz);
    mass.resize(siz);
    epsilon.resize(siz);
    sigma.resize(siz);
    epsilon14.resize(siz); 
    sigma14.resize(siz);
    beta.resize(siz);
    x.resize(siz);
    y.resize(siz);
    z.resize(siz);
    symbol.resize(siz);
    residue_id_global.resize(siz);
    residue_id_seg.resize(siz);
    res_label.resize(siz);
    seg_label.resize(siz);
    is_frozen.resize(siz);
}

/** Get the unique id and type of this atom **/
int AtomList::getID(size_t which) const
{
    return this->id[which];
}

void AtomList::setId(size_t which, int _id)
{
    this->id[which] = _id;
}

/** Coordinates manipulation methods **/
void AtomList::setCoords(size_t which, double _x, double _y, double _z)
{
    this->x[which] = _x;
    this->y[which] = _y;
    this->z[which] = _z;
}

void AtomList::setCoords(size_t which, double _crd[3])
{
    this->x[which] = _crd[0];
    this->y[which] = _crd[1];
    this->z[which] = _crd[2];
}

void AtomList::setX(size_t which, double _x)
{
    this->x[which] = _x;
}

void AtomList::setY(size_t which, double _y)
{
    this->y[which] = _y;
}

void AtomList::setZ(size_t which, double _z)
{
    this->z[which] = _z;
}

void AtomList::getCoords(size_t which, double _crd[3]) const
{
    _crd[0] = this->x[which];
    _crd[1] = this->y[which];
    _crd[2] = this->z[which];
}

double AtomList::getX(size_t which) const
{
    return this->x[which];
}

double AtomList::getY(size_t which) const
{
    return this->y[which];
}

double AtomList::getZ(size_t which) const
{
    return this->z[which];
}

double& AtomList::getX(size_t which)
{
    return this->x[which];
}

double& AtomList::getY(size_t which)
{
    return this->y[which];
}

double& AtomList::getZ(size_t which)
{
    return this->z[which];
}

void AtomList::addX(size_t which, double _x)
{
    this->x[which] += _x;
}

void AtomList::addY(size_t which, double _y)
{
    this->y[which] += _y;
}

void AtomList::addZ(size_t which, double _z)
{
    this->z[which] += _z;
}

void AtomList::addCoords(size_t which, double _x, double _y, double _z)
{
    this->x[which] += _x;
    this->y[which] += _y;
    this->z[which] += _z;
}

void AtomList::addCoords(size_t which, double _crd[3])
{
    this->x[which] += _crd[0];
    this->y[which] += _crd[1];
    this->z[which] += _crd[2];
}

/** Mass and charge manips **/
void AtomList::setCharge(size_t which, double _charge)
{
    this->charge[which] = _charge;
}

void AtomList::setMass(size_t which, double _mass)
{
   this->mass[which] = _mass;
}

double AtomList::getMass(size_t which) const
{
   return this->mass[which];
}

double AtomList::getCharge(size_t which) const
{
    return this->charge[which];
}

void AtomList::setSigma(size_t which, double _sigma)
{
    this->sigma[which] = _sigma;
}

double AtomList::getSigma(size_t which) const
{
    return this->sigma[which];
}

void AtomList::setEpsilon(size_t which, double _epsilon)
{
    this->epsilon[which] = _epsilon;
}

double AtomList::getEpsilon(size_t which) const
{
    return this->epsilon[which];
}

void AtomList::setSymbol(size_t which, const char _symbol[])
{
    this->symbol[which] = string(_symbol);
}

string AtomList::getSymbol(size_t which) const
{
    return this->symbol[which];
}

void AtomList::setIs_frozen(size_t which, bool _is_frozen)
{
    this->is_frozen[which] = _is_frozen;
}

bool AtomList::Is_frozen(size_t which) const
{
    return this->is_frozen[which];
}

void AtomList::setType(size_t which, int _type)
{
    this->type[which]= _type;
}

int AtomList::getType(size_t which) const
{
    return this->type[which];
}

void AtomList::setSeg_label(size_t which, string _seg_label)
{
    this->seg_label[which] = _seg_label;
}

void AtomList::setRes_label(size_t which, string _res_label)
{
    this->res_label[which] = _res_label;
}

void AtomList::setSymbol(size_t which, string _symbol)
{
    this->symbol[which] = _symbol;
}

void AtomList::setBeta(size_t which, double _beta)
{
    this->beta[which] = _beta;
}

double AtomList::getBeta(size_t which) const
{
    return this->beta[which];
}

void AtomList::setSigma14(size_t which, double _sigma14)
{
    this->sigma14[which] = _sigma14;
}

double AtomList::getSigma14(size_t which) const
{
    return this->sigma14[which];
}

void AtomList::setEpsilon14(size_t which, double _epsilon14)
{
    this->epsilon14[which] = _epsilon14;
}

double AtomList::getEpsilon14(size_t which) const
{
    return this->epsilon14[which];
}

void AtomList::setResidue_id_seg(size_t which, int _residue_id_seg)
{
    this->residue_id_seg[which] = _residue_id_seg;
}

int AtomList::getResidue_id_seg(size_t which) const
{
    return this->residue_id_seg[which];
}

void AtomList::setResidue_id_global(size_t which, int _residue_id_global)
{
    this->residue_id_global[which] = _residue_id_global;
}

int AtomList::getResidue_id_global(size_t which) const
{
    return this->residue_id_global[which];
}

void AtomList::setSeg_label(size_t which, const char _seg_label[])
{
    this->seg_label[which] = string(_seg_label);
}

string AtomList::getSeg_label(size_t which) const
{
    return this->seg_label[which];
}

void AtomList::setRes_label(size_t which, const char _res_label[])
{
    this->res_label[which] = string(_res_label);
}

string AtomList::getRes_label(size_t which) const
{
    return this->res_label[which];
}

const vector<double>& AtomList::getXvect() const
{
    return this->x;
}

const vector<double>& AtomList::getYvect() const
{
    return this->y;
}

const vector<double>& AtomList::getZvect() const
{
    return this->z;
}

const vector<double>& AtomList::getChargevect() const
{
    return this->charge;
}

const vector<double>& AtomList::getEpsilonvect() const
{
    return this->epsilon;
}

const vector<double>& AtomList::getSigmavect() const
{
    return this->sigma;
}

const vector<bool>& AtomList::getFrozenList() const{
  return this->is_frozen;
}

ostream& operator<<(ostream& overloadStream, const AtomList& atomlist)
{
    atomlist.toString(overloadStream);
    
    return overloadStream;
}

void AtomList::toString(ostream& stream) const
{
//     stream << fixed << setprecision(6);
    for(size_t which=0; which < id.size(); which++)
    {
    stream << "Atom" << '\t';
    stream << id[which] << '\t' << residue_id_global[which] << '\t' << res_label[which] << '\t';
    stream << symbol[which] << '\t' << x[which] << '\t' << y[which] << '\t' << z[which];
    stream << '\t' << seg_label[which] << '\t' << residue_id_seg[which] << '\t';
    //    stream << epsilon << '\t' << sigma << '\t' << epsilon14 << '\t' << sigma14 ;
    stream << charge[which] << '\t' << mass[which];
    }
}

void AtomList::crd_backup_save(vector<tuple<double, double, double >> &crdbackup, AtomList& at_List, int moveAtomList[])
{
    int ng = moveAtomList[0];
    int endng = ng + 2;
    int nn;
    int iaf, ial;
    
    double crd[3];
    
    for ( int it1 = 1; it1 <= ng; it1++ )
    {
        nn = moveAtomList[it1];
        for ( int it2 = endng; it2 <= nn; it2 += 2 )
        {
            iaf = moveAtomList[it2 - 1];
            ial = moveAtomList[it2];
            
            for ( int it3 = iaf; it3 <= ial; it3++ )
            {
                //                 cout << "Backup of crd for atom : " << it3 << endl;
                at_List.getCoords(it3,crd);
                crdbackup[it3] = tuple<double, double, double >(crd[0], crd[1], crd[2]);
            }
        }
        endng = nn + 2;
    }
}

void AtomList::crd_backup_load(vector<tuple<double, double, double >> &crdbackup, AtomList& at_List, int moveAtomList[])
{
    int ng = moveAtomList[0];
    int endng = ng + 2;
    int nn;
    int iaf, ial;
    
    double crd[3];
    
    for ( int it1 = 1; it1 <= ng; it1++ )
    {
        nn = moveAtomList[it1];
        for ( int it2 = endng; it2 <= nn; it2 += 2 )
        {
            iaf = moveAtomList[it2 - 1];
            ial = moveAtomList[it2];
            
            for ( int it3 = iaf; it3 <= ial; it3++ )
            {
                //                 cout << "Restore of crd for atom : " << it3 << endl;
                
                crd[0] = get<0>(crdbackup[it3]);
                crd[1] = get<1>(crdbackup[it3]);
                crd[2] = get<2>(crdbackup[it3]);
                
                at_List.setCoords(it3,crd);
            }
        }
        endng = nn + 2;
    }
}