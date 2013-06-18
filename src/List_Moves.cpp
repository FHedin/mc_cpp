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
#include <cstdlib>
#include <string>

#include "List_Moves.h"
#include "Atom.h"
#include "Tools.h"
#include "Selection.h"

using namespace std;

List_Moves::List_Moves(string mvtypName, string modeName,
                       std::vector<Atom>& _at_List, PerConditions& _pbc,
                       Ensemble& _ens, FField& _ff, List_Exclude& _excl)
: at_List(_at_List), pbc(_pbc), ens(_ens), ff(_ff), excl(_excl)
{
    this->nMoveTypes = 0;

    this->nMoveAtm.resize(MMVTYP, 0);
    this->moveTypeList.resize(MMVTYP);
    this->moveSeleList.resize(MMVTYP);

    this->moveBondList.resize(MMVTYP, nullptr);
    this->moveBondUpdate.resize(MMVTYP);

    this->moveAtomList.resize(MMVTYP, nullptr);

    this->moveLimitsList.resize(MMVTYP, 0.0);

    this->movePivotList.resize(MMVTYP, nullptr);

    addNewMoveType(mvtypName, modeName);

}

List_Moves::~List_Moves() { }

void List_Moves::addNewMoveType(string mvtypName, string modeName)
{
    bool success = false;
    Tools::str_rm_blank_spaces(mvtypName);
    Tools::str_to_lower_case(mvtypName);

    if ( !mvtypName.compare("trn") ) // charmm MVTYPE 1
    {
        moveTypeList.at(nMoveTypes) = TRN;
        success = NewMove_TRN_ROT(modeName);
        /* ... */
        nMoveTypes++;
    }
    else if ( !mvtypName.compare("rot") ) // MVTYPE 2
    {
        moveTypeList.at(nMoveTypes) = ROT;
        success = NewMove_TRN_ROT(modeName);
        /* ... */
        nMoveTypes++;
    }
    else if ( !mvtypName.compare("tors") ) // MVTYPE ?
    {
        moveTypeList.at(nMoveTypes) = TORS;
        /* ... */
        nMoveTypes++;
    }
    else
    {
        cout << "Warning : " << mvtypName << " is not a valid type of move ; skipping ...";
    }

    if ( !success )
    {
        nMoveTypes--;
    }
}

bool List_Moves::NewMove_TRN_ROT(string modeName)
{
    Tools::str_rm_blank_spaces(modeName);
    Tools::str_to_lower_case(modeName);

    if ( !modeName.compare("residue") ) // charmm MODE 1
    {
        moveSeleList.at(nMoveTypes) = RESIDUE;
        BOND_UPDATE bndulTmp = {false, false, false, false};
        moveBondUpdate.at(nMoveTypes) = bndulTmp;
    }
    else if ( !modeName.compare("all") )
    {
        moveSeleList.at(nMoveTypes) = ALL; // charmm MODE 2
        BOND_UPDATE bndulTmp = {false, false, false, false};
        moveBondUpdate.at(nMoveTypes) = bndulTmp;
    }
    else if ( !modeName.compare("atom") ) // charmm MODE 4
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

    vector<int> selection(natom, 0);

    /* TODO : fill selection vector */

    bool fewer = true;
    /* TODO : get fewer from XML parsing */

    /* Here we count the number of moving atoms for a given move instance
     *
     */
    nMoveAtm.at(nMoveTypes) = 0;
    for ( int i = 0; i < natom; i++ )
    {
        if ( selection.at(i) == 1 )
        {
            nMoveAtm.at(nMoveTypes)++;
        }
    }

    if ( nMoveAtm.at(nMoveTypes) == 0 )
    {
        if ( !modeName.compare("atom") )
            freeBondList();

        return false;
    }
    /* ... */

    // cofmas : for rotation, should rotation occur around a center of mass (true)
    // or a selection (false)
    // only centre of mass available now
    bool cofmas = false;

    // pointer to list of moving atoms
    int* listp = nullptr;
    int* tempp = nullptr;

    if ( moveSeleList.at(nMoveTypes) == ALL ) //mode 2
    {
        if ( moveTypeList.at(nMoveTypes) == TRN ) // mvtyp 1
        {
            makeMoveList(listp, natom, selection);
            nMoveAtm.at(nMoveTypes) = 1;
        }
        else if ( moveTypeList.at(nMoveTypes) == ROT ) // mvtyp 2
        {
            makeMoveList(listp, natom, selection);
            nMoveAtm.at(nMoveTypes) = 0;
            /* future second selection work here */
            cofmas = (nMoveAtm.at(nMoveTypes) == 0);
            if ( cofmas )
                nMoveAtm.at(nMoveTypes) = 1;
        }
    }

    // Now allocate space for the atom list according nMoveAtm size
    moveAtomList.at(nMoveTypes) = new int[ nMoveAtm.at(nMoveTypes) ];
    if ( moveSeleList.at(nMoveTypes) == ATOM ) //mode 4
    {
        moveBondList.at(nMoveTypes) = new int[ nMoveAtm.at(nMoveTypes) ];
    }

    // anisotropic move, not available now
    bool aniso = false;

    if ( aniso )
    {
        /* TODO : anisotropic moves*/
        moveLimitsList.resize(9 * natom, 0.0);
    }
    else
    {
        moveLimitsList.resize(natom, 0.0);
    }

    nMoveAtm.at(nMoveTypes) = 0;
    for ( int i = 0; i < natom; i++ )
    {
        if ( selection.at(i) == 1 && cofmas )
        {
            nMoveAtm.at(nMoveTypes)++;

            if ( moveSeleList.at(nMoveTypes) == RESIDUE ) //mode 1 by residue
            {
                gtrslf();
                moveAtomList.at(nMoveTypes) = tempp; // IMVNGP%A(NMVATM)%A => TEMPP%A
            }
            else if ( moveSeleList.at(nMoveTypes) == ALL ) // mode 2 by all
            {
                moveAtomList.at(nMoveTypes) = listp;
            }
            else if ( moveSeleList.at(nMoveTypes) == ATOM ) // mode 4 by atom
            {
                tempp = new int[4];
                tempp[0] = 2;
                tempp[1] = 4;
                tempp[2] = i;
                tempp[3] = i;
                moveAtomList.at(nMoveTypes) = tempp;
                gnbndl();
            }

            if ( aniso )
            {
                // flanis(...);
            }
            else
            {
                // MDXP%A(NMVATM) = RMDX
            }

            if ( moveTypeList.at(nMoveTypes) == TORS ) //mvtyp 2
            {
                tempp = new int[4];

                if ( cofmas )
                    tempp[0] = -1;
                else
                    tempp[0] = i;

                movePivotList.at(nMoveTypes) = tempp;
            }

            if ( (moveSeleList.at(nMoveTypes) == ALL && moveTypeList.at(nMoveTypes) == TRN)
                 || cofmas ) // mode_2_all && mvtype_1_trn
            {
                break;
            }

        }// long if ...
    }// long for ...

    /* ... */

    if ( !modeName.compare("atom") )
        freeBondList();

    return true;
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

    for ( int it = 0; it < natom; it++ )
    {
        IABNDP[it] = new int[MCMBND];
        IATHTP[it] = new int[MCMTHT];
        IAPHIP[it] = new int[MCMPHI];
        IAIMPP[it] = new int[MCMIMP];

        *(IABNDP[it]) = 0;
        *(IATHTP[it]) = 0;
        *(IAPHIP[it]) = 0;
        *(IAIMPP[it]) = 0;
    }

    // work on bonds
    const vector<Bond>& bndl = ff.getBndList();
    for ( int it = 0; it < nbond; it++ )
    {
        tmpPtr = IABNDP[ bndl[it].getAt1() ];
        fillLists(it, tmpPtr, MCMBND);

        tmpPtr = IABNDP[ bndl[it].getAt2() ];
        fillLists(it, tmpPtr, MCMBND);
    }
    tmpPtr = nullptr;

    // work on angles
    const vector<Angle>& angl = ff.getAngList();
    for ( int it = 0; it < nangl; it++ )
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
    for ( int it = 0; it < ndihe; it++ )
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
    for ( int it = 0; it < nimpr; it++ )
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
    int n = ilist[0];
    if ( n >= size )
    {
        cerr << "Error : exceeded maximum number of bonded terms." << endl;
        cerr << "From file '" << __FILE__ << "'." << endl;
        exit(-14);
    }
    ilist[n] = iic;
}

void List_Moves::freeBondList()
{
    int natom = ens.getN();

    for ( int it = 0; it < natom; it++ )
    {
        delete[] IABNDP[it];
        delete[] IATHTP[it];
        delete[] IAPHIP[it];
        delete[] IAIMPP[it];
    }
}

void List_Moves::makeMoveList(int* list, int natom, vector<int>& sele)
{
    int i, j, n, iprev = 0, npair = 0;

    // count how many pairs
    for ( i = 0; i < natom; i++ )
    {
        if ( iprev == 0 && sele[i] >= 1 )
        {
            npair++;
        }
        iprev = sele[i];
    }

    // fill the list
    n = 2 * npair + 2;
    list = new int[n];
    list[0] = 2;
    list[1] = n;

    j = 2;
    iprev = 0;
    for ( i = 0; i < natom; i++ )
    {
        if ( iprev == 0 && sele[i] >= 1 )
        {
            list[j] = i;
            j++;
        }
        else if ( iprev >= 1 && sele[i] == 0 )
        {
            list[j] = i - 1;
            j++;
        }
        iprev = sele[i];
    }

    if ( sele[natom - 1] >= 1 )
    {
        list[j] = natom;
    }
}

