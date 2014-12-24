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

#ifndef ATOM_H
#define ATOM_H

class PerConditions;

#include <iostream>
#include <vector>
#include <string>
#include <tuple>

//#include "Global_include.hpp"

#include "PerConditions.hpp"

class Atom
{
    friend std::ostream& operator<<( std::ostream& out, const Atom& atom );

public:
    /** Public Methods**/
    Atom(int _id, std::string symbol);
    Atom();
    virtual ~Atom();

    // return the unique atom ID and the type
    int getID() const;
    void setId(int id);

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

    void setMass(double mass);
    double getMass() const;

    //    int getId() const;
    void setIs_frozen(bool is_frozen);
    bool Is_frozen() const;

    void setType(int type);
    int getType() const;

    void setSeg_label(std::string seg_label);
    void setRes_label(std::string res_label);
    void setSymbol(std::string symbol);
    void setBeta(double beta);
    double getBeta() const;
    void setSigma14(double sigma14);
    double getSigma14() const;
    void setEpsilon14(double epsilon14);
    double getEpsilon14() const;

    void addX(double _x);
    void addY(double _y);
    void addZ(double _z);
    void addCoords(double _x, double _y, double _z);
    void addCoords(double _crd[3]);

    static void crd_backup_save(std::vector<std::tuple<double, double, double >> &crdbackup, std::vector<Atom>& at_List, int moveAtomList[]);
    static void crd_backup_load(std::vector<std::tuple<double, double, double >> &crdbackup, std::vector<Atom>& at_List, int moveAtomList[]);

    /** Public Attributes **/
    // nothing for the moment

private:
    /** Private Attributes **/
    int id; // unique identifier for this atom
    int type; // internal type, forcefield dependent

    double charge; // electrostatic
    double mass; // useless for mc but needed for forcefields
    double epsilon, sigma, epsilon14, sigma14; // Lennard-Jones parameters
    double beta; //Buckingham potential parameter
    double x, y, z; // coordinates

    std::string symbol;
    int residue_id_global;
    int residue_id_seg;
    std::string res_label;
    std::string seg_label;

    bool is_frozen;

    virtual void toString(std::ostream& stream) const;

};

#endif // ATOM_H