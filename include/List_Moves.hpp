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

#ifndef LIST_MOVES_H
#define	LIST_MOVES_H

#include <string>
#include <iostream>
#include <vector>

//#include "Global_include.hpp"

#include "AtomList.hpp"
#include "FField.hpp"
#include "Move.hpp"

struct BOND_UPDATE
{
    bool bonds;
    bool angles;
    bool dihe;
    bool impr;
};

/*CHARMM MODE in movead.src*/
enum MOVEMODE
{
    RESIDUE = 1,
    ALL = 2,
    //    HEAVY = 3,
    ATOM = 4
};

class List_Moves
{
    friend std::ostream& operator<<( std::ostream& overloadStream, const List_Moves& lst );
public:
    List_Moves(AtomList& _at_List, FField& _ff, int _natom);
    virtual ~List_Moves();

    const std::vector<int**>& getMovePivotList() const;
    std::vector<double>& getMoveLimitsList();
    const std::vector<double>& getTargetAcceptanceList() const;
    const std::vector<int>& getMoveUpdateFreqList() const;
    const std::vector<int**>& getMoveAtomList() const;
    const std::vector<BOND_UPDATE>& getMoveBondUpdate() const;
    const std::vector<int**>& getMoveBondList() const;
    const std::vector<MOVEMODE>& getMoveModeList() const;
    const std::vector<MOVETYPE>& getMoveTypeList() const;
    const std::vector<int>& getNMoveAtm() const;
    int getNMoveTypes() const;

    void addNewMoveType(std::string mvtypName, std::string modeName,
                        std::string selMode, std::string selName,
                        double dmax_value, double dmax_target, int dmax_each);
    
    void addNewMoveType(std::string mvtypName, std::string modeName,
                        std::string selMode, std::string selName);
    
    //void addNewMoveType(std::string mvtypName, std::string modeName,
    //                    std::vector<std::tuple<std::string,std::string>> seleList);

private:

    AtomList& at_List;
    FField& ff;
    int natom;

    //    MCMBND       Maximum number of bonds     for a single atom.
    //    MCMTHT       Maximum number of angles    for a single atom.
    //    MCMPHI       Maximum number of dihedrals for a single atom.
    //    MCMIMP       Maximum number of impropers for a single atom.
    static const int MCMBND, MCMTHT, MCMPHI, MCMIMP, MMVTYP;
    // number of move types for this given simulation ; NMVTYP in CHARMM's MC.
    int nMoveTypes; // max is 50

    /* number of moving instances for a given nMoveTypes type.
     * One instance can contain several or all atoms. Or just one.
     * in CHARMM it is NMVATM
     */
    std::vector<int> nMoveAtm;

    // vector of size nMoveTypes : stores the type of the move. MVTYPE in CHARMM.
    std::vector<MOVETYPE> moveTypeList;
    // selection mode for a given move instance. MODE in CHARMM.
    std::vector<MOVEMODE> moveModeList;

    /* List of pointers to another list of bonded terms altered by the current move.
     * IBLSTP in CHARMM MC
     * Size is nMoveTypes
     */
    std::vector<int**> moveBondList;

    /* This is for storing a boolean which indicates if bonded terms are altered by a given movetype.
     * Vector is of size nMoveTypes, and for each it stores a structure of 4 booleans
     * in CHARMM it is a 2d array called QBND
     */
    std::vector<BOND_UPDATE> moveBondUpdate;

    /* natom length vectors with pointers to lists of bond,
     * angles, dihedrals and impropers affected by that atom
     */
    std::vector<int*> IABNDP; // bonds
    std::vector<int*> IATHTP; // angles
    std::vector<int*> IAPHIP; // dihedrals
    std::vector<int*> IAIMPP; // impropers

    /*
     * Contains a compacted list of moving atoms.
     * See IMVNGP in CHARMM.
     */
    std::vector<int**> moveAtomList;

    /*
     * Vector of maximal random moves (in Angstroems or degrees),
     * and vector of desired acceptance target, and vector of update frequency of the first vector for reaching target value
     */
    std::vector<double> moveLimitsList;
    std::vector<double> targetAcceptanceList;
    std::vector<int> moveUpdateFreqList;

    // list of pivots, i.e. for rotations
    std::vector<int**> movePivotList;

    void addDmaxValues(double dmax_value, double dmax_target, int dmax_each);

    // move methods all private by default
    bool NewMove_TRN_ROT(std::string modeName, std::string selMode, std::string selName);
    //    bool NewMove_TORS(std::string modeName, std::string selMode, std::string selName);

    void makeBondList();
    void fillLists(int iic, int* ilist, int size) const;
    void freeBondList();
    void makeMoveList(int*& list, int natom, const std::vector<int>& sele) const;
    void gtrsfl(int*& listp, int atomidx, int& a1, int& a2, int natom) const;
    void gnbndl(int atomidx);
    int nbtf(int* list, bool isActive) const;
    void assibl(int& ne, int n, int *orig, int* dest) const;

    void toString(std::ostream& stream) const;

    void fillNull(int** array, int size) const;
};

#endif	/* LIST_MOVES_H */

