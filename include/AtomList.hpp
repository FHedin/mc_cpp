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

#ifndef ATOMLIST_H
#define ATOMLIST_H

class PerConditions;

#include <vector>
#include <string>
#include <iostream>
#include <tuple>

#include "PerConditions.hpp"

class AtomList
{
    friend std::ostream& operator<<( std::ostream& out, const AtomList& atomlist );

public:
    AtomList();
    ~AtomList();

    //for array like access
    AtomList& operator[](size_t index);
//     const AtomList& operator[](size_t index) const;
    AtomList& at(size_t index);
//     const AtomList& at(size_t index) const;

    void resize(size_t siz);
    void alloc(size_t siz);

    // return the unique atom ID and the type
    int getID() const;
    void setId(int _id);

    // different ways of setting/getting/adding
    void setCoords(double _x, double _y, double _z);
    void setCoords(double _crd[3]);
    void getCoords(double _crd[3]) const;

    void setX(double _x);
    void setY(double _y);
    void setZ(double _z);
    double getX() const;
    double getY() const;
    double getZ() const;
    double& getX();
    double& getY();
    double& getZ();

    void setCharge(double _charge);
    void setSigma(double _sigma);
    void setEpsilon(double _epsilon);
    double getCharge() const;
    double getSigma() const;
    double getEpsilon() const;

    void setResidue_id_seg(int residue_id_seg);
    int getResidue_id_seg() const;

    void setResidue_id_global(int residue_id_global);
    int getResidue_id_global() const;

    void setSeg_label(const char[]);
    std::string getSeg_label() const;

    void setRes_label(const char[]);
    std::string getRes_label() const;

    void setSymbol(const char[]);
    std::string getSymbol() const;

    void setMass(double _mass);
    double getMass() const;

    void setIs_frozen(bool _is_frozen);
    bool Is_frozen() const;

    void setType(int _type);
    int getType() const;

    void setSeg_label(std::string _seg_label);
    void setRes_label(std::string _res_label);
    void setSymbol(std::string _symbol);
    void setBeta(double _beta);
    double getBeta() const;
    void setSigma14(double _sigma14);
    double getSigma14() const;
    void setEpsilon14(double _epsilon14);
    double getEpsilon14() const;

    void addX(double _x);
    void addY(double _y);
    void addZ(double _z);
    void addCoords(double _x, double _y, double _z);
    void addCoords(double _crd[3]);

    //vector accessors
//     const std::vector<double>& getXvect() const;
//     const std::vector<double>& getYvect() const;
//     const std::vector<double>& getZvect() const;
//     const std::vector<double>& getChargevect() const;
//     const std::vector<double>& getSigmavect() const;
//     const std::vector<double>& getEpsilonvect() const;
    const double* getXvect() const;
    const double* getYvect() const;
    const double* getZvect() const;
    const double* getChargevect() const;
    const double* getSigmavect() const;
    const double* getEpsilonvect() const;

    static void crd_backup_save(std::vector<std::tuple<double, double, double >> &crdbackup, AtomList& at_List, int moveAtomList[]);
    static void crd_backup_load(std::vector<std::tuple<double, double, double >> &crdbackup, AtomList& at_List, int moveAtomList[]);

private:
    /** For accessing this object as an array **/
    size_t arrayIndex;

    // unique identifier for this atom
    // internal type, forcefield dependent
    std::vector<int> id,type;

    // coordinates
//     std::vector<double> x, y, z;
    // electrostatic and atomic mass
//     std::vector<double> charge,mass;
    // Lennard-Jones parameters
//     std::vector<double> epsilon, sigma;
    
    // coordinates
    double *x, *y, *z;
    // electrostatic
    double *charge;
    // Lennard-Jones parameters
    double *epsilon, *sigma;
    
    //the mass
    std::vector<double> mass;
    
    // Lennard-Jones 1-4 parameters
    std::vector<double> epsilon14, sigma14;
    std::vector<double> beta; //Buckingham potential parameter
    
    // strings and ids for cor and psf
    std::vector<std::string> symbol;
    std::vector<int> residue_id_global;
    std::vector<int> residue_id_seg;
    std::vector<std::string> res_label;
    std::vector<std::string> seg_label;

    // is the atom frozen ?
    std::vector<bool> is_frozen;

    virtual void toString(std::ostream& stream) const;
};

#endif // ATOMLIST_H
