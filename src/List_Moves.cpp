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
#include <string>

#include "List_Moves.h"
#include "Atom.h"
#include "Tools.h"

using namespace std;

List_Moves::List_Moves(std::vector<Atom>& _at_List, PerConditions& _pbc,
        Ensemble& _ens, FField& _ff, List_Exclude& _excl)
: at_List(_at_List), pbc(_pbc), ens(_ens), ff(_ff), excl(_excl)
{
    this->nMoveTypes = 0;
    this->nMoveAtm.resize(MMVTYP, 0);
    this->moveTypeList.resize(MMVTYP);
    this->moveSeleList.resize(MMVTYP);
    this->moveBondList.resize(MMVTYP, nullptr);
    this->moveBondUpdate.resize(MMVTYP);
}

List_Moves::~List_Moves()
{
    int natom = ens.getN();

    for (int it = 0; it < natom; it++)
    {
        delete[] IABNDP[it];
        delete[] IATHTP[it];
        delete[] IAPHIP[it];
        delete[] IAIMPP[it];
    }
}

void List_Moves::addNewMoveType(string mvtypName, string modeName)
{
    bool success = false;
    Tools::str_rm_blank_spaces(mvtypName);
    Tools::str_to_lower_case(mvtypName);

    if (!mvtypName.compare("trn")) // charmm MVTYPE 1
    {
        moveTypeList.at(nMoveTypes) = TRN;
        success = NewMove_TRN_ROT(modeName);
        /* ... */
        nMoveTypes++;
    }
    else if (!mvtypName.compare("rot")) // MVTYPE 2
    {
        moveTypeList.at(nMoveTypes) = ROT;
        success = NewMove_TRN_ROT(modeName);
        /* ... */
        nMoveTypes++;
    }
    else if (!mvtypName.compare("tors")) // MVTYPE ?
    {
        moveTypeList.at(nMoveTypes) = TORS;
        /* ... */
        nMoveTypes++;
    }
    else
    {
        cout << "Warning : " << mvtypName << " is not a valid type of move ; skipping ...";
    }

    if (!success)
    {
        nMoveTypes--;
    }
}

bool List_Moves::NewMove_TRN_ROT(string modeName)
{
    Tools::str_rm_blank_spaces(modeName);
    Tools::str_to_lower_case(modeName);

    if (!modeName.compare("residue")) // charmm MODE 1
    {
        moveSeleList.at(nMoveTypes) = RESIDUE;
        BOND_UPDATE bndulTmp = {false, false, false, false};
        moveBondUpdate.at(nMoveTypes) = bndulTmp;
    }
    else if (!modeName.compare("all"))
    {
        moveSeleList.at(nMoveTypes) = ALL; // charmm MODE 2
        BOND_UPDATE bndulTmp = {false, false, false, false};
        moveBondUpdate.at(nMoveTypes) = bndulTmp;
    }
    else if (!modeName.compare("atom")) // charmm MODE 4
    {
        moveSeleList.at(nMoveTypes) = ATOM;
        BOND_UPDATE bndulTmp = {true, true, true, true};
        moveBondUpdate.at(nMoveTypes) = bndulTmp;
        makeBondList();
    }
    else
    {
        cout << "Warning : the following selection mode is not available : "
                << modeName << " skipping ..." << endl;

        return false;
    }

    const int natom = ens.getN();
    vector<bool> selection(natom, false);

    /* TODO : fill selection vector */

    bool fewer = true;
    /* TODO : get fewer from XML parsing */

    nMoveAtm.at(nMoveTypes) = 0;
    for (int i = 0; i < natom; i++)
    {
        if (selection.at(i))
        {
            nMoveAtm.at(nMoveTypes)++;
        }
    }

    /* ... */
}

void List_Moves::makeBondList()
{
    int* tmpPtr = nullptr;

    const int natom = ens.getN();
    const int nbond = ff.getNBond();
    const int nangl = ff.getNAngle();
    const int ndihe = ff.getNDihedral();
    const int nimpr = ff.getNImproper();

    IABNDP.resize(natom, nullptr);
    IATHTP.resize(natom, nullptr);
    IAPHIP.resize(natom, nullptr);
    IAIMPP.resize(natom, nullptr);

    for (int it = 0; it < natom; it++)
    {
        IABNDP[it] = new int(MCMBND);
        IATHTP[it] = new int(MCMTHT);
        IAPHIP[it] = new int(MCMPHI);
        IAIMPP[it] = new int(MCMIMP);
    }

    // work on bonds
    const vector<Bond>& bndl = ff.getBndList();
    for (int it = 0; it < nbond; it++)
    {
        tmpPtr = IABNDP[ bndl[it].getAt1() ];
        fillLists(it, tmpPtr, MCMBND);

        tmpPtr = IABNDP[ bndl[it].getAt2() ];
        fillLists(it, tmpPtr, MCMBND);
    }
    tmpPtr = nullptr;

    // work on angles
    const vector<Angle>& angl = ff.getAngList();
    for (int it = 0; it < nangl; it++)
    {
        tmpPtr = IATHTP[ angl[it].getAt1() ];
        fillLists(it, tmpPtr, MCMTHT);

        tmpPtr = IATHTP[ angl[it].getAt2() ];
        fillLists(it, tmpPtr, MCMTHT);

        tmpPtr = IATHTP[ angl[it].getAt3() ];
        fillLists(it, tmpPtr, MCMTHT);
    }
    tmpPtr = nullptr;

    // work on dihedrals
    const vector<Dihedral>& dihel = ff.getDiheList();
    for (int it = 0; it < ndihe; it++)
    {
        tmpPtr = IAPHIP[ dihel[it].getAt1() ];
        fillLists(it, tmpPtr, MCMPHI);

        tmpPtr = IAPHIP[ dihel[it].getAt2() ];
        fillLists(it, tmpPtr, MCMPHI);

        tmpPtr = IAPHIP[ dihel[it].getAt3() ];
        fillLists(it, tmpPtr, MCMPHI);

        tmpPtr = IAPHIP[ dihel[it].getAt4() ];
        fillLists(it, tmpPtr, MCMPHI);
    }
    tmpPtr = nullptr;


    // work on impropers
    const vector<Dihedral_improper>& imprl = ff.getImprList();
    for (int it = 0; it < nimpr; it++)
    {
        tmpPtr = IAIMPP[ imprl[it].getAt1() ];
        fillLists(it, tmpPtr, MCMIMP);

        tmpPtr = IAIMPP[ imprl[it].getAt2() ];
        fillLists(it, tmpPtr, MCMIMP);

        tmpPtr = IAIMPP[ imprl[it].getAt3() ];
        fillLists(it, tmpPtr, MCMIMP);

        tmpPtr = IAIMPP[ imprl[it].getAt4() ];
        fillLists(it, tmpPtr, MCMIMP);
    }
    tmpPtr = nullptr;
}

void List_Moves::fillLists(int iic, int* ilist, int size)
{
    ilist[0] += 1;
    int n = ilist[0] + 1;
    if (n > size)
    {
        cerr << "Error : exceeded maximum number of bonded terms." << endl;
        cerr << "From file '" << __FILE__ << "'." << endl;
    }
    ilist[n] = iic;
}

