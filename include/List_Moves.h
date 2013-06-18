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
#include <vector>

#include "Atom.h"
#include "Ensemble.h"
#include "PerConditions.h"
#include "FField.h"
#include "List_Exclude.h"

enum MVTYP
{
    TRN = 1,
    ROT = 2,
    TORS = 4
};

enum SELEMODE
{
    RESIDUE = 1,
    ALL = 2,
    //    HEAVY = 3,
    ATOM = 4
};

struct BOND_UPDATE
{
    bool bonds;
    bool angles;
    bool dihe;
    bool impr;
};

class List_Moves
{
public:
    List_Moves(std::string mvtypName, std::string modeName,
               std::vector<Atom>& _at_List, PerConditions& _pbc,
               Ensemble& _ens, FField& _ff, List_Exclude& _excl);

    virtual ~List_Moves();

    void addNewMoveType(const std::string mvtypName, const std::string modeName);

private:

    std::vector<Atom>& at_List;
    PerConditions& pbc;
    Ensemble& ens;
    FField& ff;
    List_Exclude& excl;

    //    MCMBND       Maximum number of bonds     for a single atom.
    //    MCMTHT       Maximum number of angles    for a single atom.
    //    MCMPHI       Maximum number of dihedrals for a single atom.
    //    MCMIMP       Maximum number of impropers for a single atom.
    const int MCMBND = 10, MCMTHT = 20, MCMPHI = 35, MCMIMP = 20;
    const int MMVTYP = 50;
    // number of move types for this given simulation ; NMVTYP in CHARMM's MC.
    int nMoveTypes; // max is 50

    /* number of moving instances for a given nMoveTypes type.
     * One instance can contain several or all atoms. Or just one.
     * in CHARMM it is NMVATM
     */
    std::vector<int> nMoveAtm;

    // vector of size nMoveTypes : stores the type of the move. MVTYPE in CHARMM.
    std::vector<MVTYP> moveTypeList;
    // selection mode for a given move instance. MODE in CHARMM.
    std::vector<SELEMODE> moveSeleList;

    /* List of pointers to another list of bonded terms altered by the current move.
     * IBLSTP in CHARMM MC
     * Size is nMoveTypes
     */
    std::vector<int*> moveBondList;

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
    std::vector<int*> moveAtomList;

    /*
     * A vector of dmax values
     */
    std::vector<double> moveLimitsList;

    // list of pivots, i.e. for rotations
    std::vector<int*> movePivotList;

    // move methods all private by default
    bool NewMove_TRN_ROT(std::string modeName);
    void makeBondList();
    void fillLists(int iic, int* ilist, int size);
    void freeBondList();
    void makeMoveList(int* list, int natom, std::vector<int>& sele);
};

#endif	/* LIST_MOVES_H */

