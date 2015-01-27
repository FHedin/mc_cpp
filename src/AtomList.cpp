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

#include <cstdlib>

#include "AtomList.hpp"

using namespace std;

AtomList::AtomList()
{

}

AtomList::~AtomList()
{
    free(charge);
    free(epsilon);
    free(sigma);
    free(x);
    free(y);
    free(z);
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

void AtomList::resize(size_t siz)
{
    id.resize(siz);
    type.resize(siz);
//     charge.resize(siz);
    mass.resize(siz);
//     epsilon.resize(siz);
//     sigma.resize(siz);
    epsilon14.resize(siz);
    sigma14.resize(siz);
    beta.resize(siz);
//     x.resize(siz);
//     y.resize(siz);
//     z.resize(siz);
    symbol.resize(siz);
    residue_id_global.resize(siz);
    residue_id_seg.resize(siz);
    res_label.resize(siz);
    seg_label.resize(siz);
    is_frozen.resize(siz);
}

void AtomList::alloc(size_t siz)
{
    posix_memalign((void**)&charge,32,siz*sizeof(double));
    posix_memalign((void**)&epsilon,32,siz*sizeof(double));
    posix_memalign((void**)&sigma,32,siz*sizeof(double));
    posix_memalign((void**)&x,32,siz*sizeof(double));
    posix_memalign((void**)&y,32,siz*sizeof(double));
    posix_memalign((void**)&z,32,siz*sizeof(double));
}

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

double& AtomList::getX()
{
    return this->x[arrayIndex];
}

double& AtomList::getY()
{
    return this->y[arrayIndex];
}

double& AtomList::getZ()
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

void AtomList::setSymbol(const char _symbol[])
{
    this->symbol[arrayIndex] = string(_symbol);
}

string AtomList::getSymbol() const
{
    return this->symbol[arrayIndex];
}

void AtomList::setIs_frozen(bool _is_frozen)
{
    this->is_frozen[arrayIndex] = _is_frozen;
}

bool AtomList::Is_frozen() const
{
    return this->is_frozen[arrayIndex];
}

void AtomList::setType(int _type)
{
    this->type[arrayIndex]= _type;
}

int AtomList::getType() const
{
    return this->type[arrayIndex];
}

void AtomList::setSeg_label(string _seg_label)
{
    this->seg_label[arrayIndex] = _seg_label;
}

void AtomList::setRes_label(string _res_label)
{
    this->res_label[arrayIndex] = _res_label;
}

void AtomList::setSymbol(string _symbol)
{
    this->symbol[arrayIndex] = _symbol;
}

void AtomList::setBeta(double _beta)
{
    this->beta[arrayIndex] = _beta;
}

double AtomList::getBeta() const
{
    return this->beta[arrayIndex];
}

void AtomList::setSigma14(double _sigma14)
{
    this->sigma14[arrayIndex] = _sigma14;
}

double AtomList::getSigma14() const
{
    return this->sigma14[arrayIndex];
}

void AtomList::setEpsilon14(double _epsilon14)
{
    this->epsilon14[arrayIndex] = _epsilon14;
}

double AtomList::getEpsilon14() const
{
    return this->epsilon14[arrayIndex];
}

void AtomList::setResidue_id_seg(int _residue_id_seg)
{
    this->residue_id_seg[arrayIndex] = _residue_id_seg;
}

int AtomList::getResidue_id_seg() const
{
    return this->residue_id_seg[arrayIndex];
}

void AtomList::setResidue_id_global(int _residue_id_global)
{
    this->residue_id_global[arrayIndex] = _residue_id_global;
}

int AtomList::getResidue_id_global() const
{
    return this->residue_id_global[arrayIndex];
}

void AtomList::setSeg_label(const char _seg_label[])
{
    this->seg_label[arrayIndex] = string(_seg_label);
}

string AtomList::getSeg_label() const
{
    return this->seg_label[arrayIndex];
}

void AtomList::setRes_label(const char _res_label[])
{
    this->res_label[arrayIndex] = string(_res_label);
}

string AtomList::getRes_label() const
{
    return this->res_label[arrayIndex];
}

// const vector<double>& AtomList::getXvect() const
// {
//     return this->x;
// }
// 
// const vector<double>& AtomList::getYvect() const
// {
//     return this->y;
// }
// 
// const vector<double>& AtomList::getZvect() const
// {
//     return this->z;
// }
// 
// const vector<double>& AtomList::getChargevect() const
// {
//     return this->charge;
// }
// 
// const vector<double>& AtomList::getEpsilonvect() const
// {
//     return this->epsilon;
// }
// 
// const vector<double>& AtomList::getSigmavect() const
// {
//     return this->sigma;
// }

const double* AtomList::getXvect() const
{
    return this->x;
}

const double* AtomList::getYvect() const
{
    return this->y;
}

const double* AtomList::getZvect() const
{
    return this->z;
}

const double* AtomList::getChargevect() const
{
    return this->charge;
}

const double* AtomList::getEpsilonvect() const
{
    return this->epsilon;
}

const double* AtomList::getSigmavect() const
{
    return this->sigma;
}

ostream& operator<<(ostream& overloadStream, const AtomList& atomlist)
{
    atomlist.toString(overloadStream);

    return overloadStream;
}

void AtomList::toString(ostream& stream) const
{
//     stream << fixed << setprecision(6);
    stream << "Atom" << '\t';
    stream << id[arrayIndex] << '\t' << residue_id_global[arrayIndex] << '\t' << res_label[arrayIndex] << '\t';
    stream << symbol[arrayIndex] << '\t' << x[arrayIndex] << '\t' << y[arrayIndex] << '\t' << z[arrayIndex];
    stream << '\t' << seg_label[arrayIndex] << '\t' << residue_id_seg[arrayIndex] << '\t';
    //    stream << epsilon << '\t' << sigma << '\t' << epsilon14 << '\t' << sigma14 ;
    stream << charge[arrayIndex] << '\t' << mass[arrayIndex];
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
                at_List[it3].getCoords(crd);
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

                at_List[it3].setCoords(crd);
            }
        }
        endng = nn + 2;
    }
}
