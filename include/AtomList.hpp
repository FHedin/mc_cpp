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
    
    void resize(size_t siz);
    
    // return the unique atom ID and the type
    int getID(size_t which) const;
    void setId(size_t which, int _id);
    
    // different ways of setting/getting/adding
    void setCoords(size_t which, double _x, double _y, double _z);
    void setCoords(size_t which, double _crd[3]);
    void getCoords(size_t which, double _crd[3]) const;
    
    void setX(size_t which, double _x);
    void setY(size_t which, double _y);
    void setZ(size_t which, double _z);
    double getX(size_t which) const;
    double getY(size_t which) const;
    double getZ(size_t which) const;
    double& getX(size_t which);
    double& getY(size_t which);
    double& getZ(size_t which);
    
    void setCharge(size_t which, double _charge);
    void setSigma(size_t which, double _sigma);
    void setEpsilon(size_t which, double _epsilon);
    double getCharge(size_t which) const;
    double getSigma(size_t which) const;
    double getEpsilon(size_t which) const;
    
    void setResidue_id_seg(size_t which, int residue_id_seg);
    int getResidue_id_seg(size_t which) const;
    
    void setResidue_id_global(size_t which, int residue_id_global);
    int getResidue_id_global(size_t which) const;
    
    void setSeg_label(size_t which, const char[]);
    std::string getSeg_label(size_t which) const;
    
    void setRes_label(size_t which, const char[]);
    std::string getRes_label(size_t which) const;
    
    void setSymbol(size_t which, const char[]);
    std::string getSymbol(size_t which) const;
    
    void setMass(size_t which, double _mass);
    double getMass(size_t which) const;
    
    void setIs_frozen(size_t which, bool _is_frozen);
    bool Is_frozen(size_t which) const;
    
    void setType(size_t which, int _type);
    int getType(size_t which) const;
    
    void setSeg_label(size_t which, std::string _seg_label);
    void setRes_label(size_t which, std::string _res_label);
    void setSymbol(size_t which, std::string _symbol);
    void setBeta(size_t which, double _beta);
    double getBeta(size_t which) const;
    void setSigma14(size_t which, double _sigma14);
    double getSigma14(size_t which) const;
    void setEpsilon14(size_t which, double _epsilon14);
    double getEpsilon14(size_t which) const;
    
    void addX(size_t which, double _x);
    void addY(size_t which, double _y);
    void addZ(size_t which, double _z);
    void addCoords(size_t which, double _x, double _y, double _z);
    void addCoords(size_t which, double _crd[3]);
    
    //vector accessors
    const std::vector<double>& getXvect() const;
    const std::vector<double>& getYvect() const;
    const std::vector<double>& getZvect() const;
    const std::vector<double>& getChargevect() const;
    const std::vector<double>& getSigmavect() const;
    const std::vector<double>& getEpsilonvect() const;
    
    static void crd_backup_save(std::vector<std::tuple<double, double, double >> &crdbackup, AtomList& at_List, int moveAtomList[]);
    static void crd_backup_load(std::vector<std::tuple<double, double, double >> &crdbackup, AtomList& at_List, int moveAtomList[]);
    
private:
    /** For accessing this object as an array **/
//     size_t arrayIndex;
    
    // unique identifier for this atom
    // internal type, forcefield dependent
    std::vector<int> id,type;
    
    // electrostatic and atomic mass
    std::vector<double> charge,mass;
    
    // Lennard-Jones parameters
    std::vector<double> epsilon, sigma, epsilon14, sigma14;
    std::vector<double> beta; //Buckingham potential parameter
    std::vector<double> x, y, z; // coordinates
    
    std::vector<std::string> symbol;
    std::vector<int> residue_id_global;
    std::vector<int> residue_id_seg;
    std::vector<std::string> res_label;
    std::vector<std::string> seg_label;
    
    std::vector<bool> is_frozen;
    
    virtual void toString(std::ostream& stream) const;
};

#endif // ATOMLIST_H
