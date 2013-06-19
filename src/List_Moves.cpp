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
#include <algorithm>

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
    this->moveModeList.resize(MMVTYP);

    this->moveBondList.resize(MMVTYP, nullptr);
    this->moveBondUpdate.resize(MMVTYP);

    this->moveAtomList.resize(MMVTYP, nullptr);

    this->moveLimitsList.resize(MMVTYP, 0.0);

    this->movePivotList.resize(MMVTYP, nullptr);

    addNewMoveType(mvtypName, modeName);

}

List_Moves::~List_Moves()
{
    for ( int i = 0; i < MMVTYP; i++ )
    {
        delete[] moveBondList[i];
        delete[] moveAtomList[i];
        delete[] movePivotList[i];
    }
}

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
        moveModeList.at(nMoveTypes) = RESIDUE;
        BOND_UPDATE bndulTmp = {false, false, false, false};
        moveBondUpdate.at(nMoveTypes) = bndulTmp;
        moveBondList.at(nMoveTypes) = nullptr;
    }
    else if ( !modeName.compare("all") )
    {
        moveModeList.at(nMoveTypes) = ALL; // charmm MODE 2
        BOND_UPDATE bndulTmp = {false, false, false, false};
        moveBondUpdate.at(nMoveTypes) = bndulTmp;
        moveBondList.at(nMoveTypes) = nullptr;
    }
    else if ( !modeName.compare("atom") ) // charmm MODE 4
    {
        moveModeList.at(nMoveTypes) = ATOM;
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

    Selection selec("RESIDUE_NAME", "NMA", at_List, natom);
    const vector<int>& seleList = selec.getSelection();
    
    //    for ( auto iter : seleList )
    //    {
    //        cout << iter << '\t';
    //    }
    //    cout << endl;
    //    exit(0);


    //    bool fewer = true;
    /* TODO : get fewer from XML parsing */

    /* Here we count the number of moving atoms for a given move instance
     *
     */
    nMoveAtm.at(nMoveTypes) = 0;
    for ( int i = 0; i < natom; i++ )
    {
        if ( seleList.at(i) == 1 )
        {
            nMoveAtm.at(nMoveTypes)++;
        }
    }

    //    cout << "On line " << __LINE__ << '\t' << *this;

    if ( nMoveAtm.at(nMoveTypes) == 0 )
    {
        if ( moveModeList.at(nMoveTypes) == ATOM )
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
    int ia1, ia2;

    if ( moveModeList.at(nMoveTypes) == ALL ) //mode 2
    {
        if ( moveTypeList.at(nMoveTypes) == TRN ) // mvtyp 1
        {
            makeMoveList(listp, natom, seleList);
            nMoveAtm.at(nMoveTypes) = 1;
        }
        else if ( moveTypeList.at(nMoveTypes) == ROT ) // mvtyp 2
        {
            makeMoveList(listp, natom, seleList);
            nMoveAtm.at(nMoveTypes) = 0;
            /* future second selection work here */
            cofmas = (nMoveAtm.at(nMoveTypes) == 0);
            if ( cofmas )
                nMoveAtm.at(nMoveTypes) = 1;
        }
    }

    // Now allocate space for the atom list according nMoveAtm size
    moveAtomList.at(nMoveTypes) = new int *[ nMoveAtm.at(nMoveTypes) ];

    if ( moveTypeList.at(nMoveTypes) == ROT ) //mvtyp 2
    {
        movePivotList.at(nMoveTypes) = new int *[nMoveAtm.at(nMoveTypes)];
    }

    if ( moveModeList.at(nMoveTypes) == ATOM ) //mode 4
    {
        //        cout << "On line " << __LINE__ << '\t' << *this;
        moveBondList.at(nMoveTypes) = new int *[ nMoveAtm.at(nMoveTypes) ];
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
        if ( seleList.at(i) == 1 || cofmas )
        {
            nMoveAtm.at(nMoveTypes)++;

            if ( moveModeList.at(nMoveTypes) == RESIDUE ) //mode 1 by residue
            {
                gtrsfl(tempp, i, ia1, ia2, natom);
                moveAtomList.at(nMoveTypes)[nMoveAtm.at(nMoveTypes)] = tempp; // IMVNGP%A(NMVATM)%A => TEMPP%A
            }
            else if ( moveModeList.at(nMoveTypes) == ALL ) // mode 2 by all
            {
                moveAtomList.at(nMoveTypes)[nMoveAtm.at(nMoveTypes)] = listp;
            }
            else if ( moveModeList.at(nMoveTypes) == ATOM ) // mode 4 by atom
            {
                //                cout << "On line " << __LINE__ << '\t' << *this;
                tempp = new int[4];
                tempp[0] = 2 - 1;
                tempp[1] = 4 - 1;
                tempp[2] = i;
                tempp[3] = i;
                moveAtomList.at(nMoveTypes)[nMoveAtm.at(nMoveTypes)] = tempp;
                gnbndl(i);
            }

            if ( aniso )
            {
                // flanis(...);
            }
            else
            {
                // MDXP%A(NMVATM) = RMDX
            }

            if ( moveTypeList.at(nMoveTypes) == ROT ) //mvtyp 2
            {
                tempp = new int[4];

                if ( cofmas )
                    tempp[0] = -1;
                else
                    tempp[0] = i;

                movePivotList.at(nMoveTypes)[nMoveAtm.at(nMoveTypes)] = tempp;
            }

            if ( (moveModeList.at(nMoveTypes) == ALL && moveTypeList.at(nMoveTypes) == TRN)
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

void List_Moves::fillLists(int iic, int* ilist, int size) const
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

void List_Moves::makeMoveList(int* list, int natom, const vector<int>& sele) const
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

// finds the first and last atom of residue, for RESIDUE selection mode

void List_Moves::gtrsfl(int *listp, int atomidx, int& a1, int& a2, int natom) const
{
    int iresf, iresl, mid;

    iresf = 1;

    // get the total number of residues
    int* resList = new int[natom];
    for ( int it = 0; it < natom; it++ )
        resList[it] = at_List.at(it).getResidue_id_global();
    iresl = *max_element(resList, resList + natom);


    // now build a list of size iresl+1 where iresl[0]=-1
    // and iresl[i] returns the last atom of residue i.
    int* resLast = new int[iresl + 1];
    resLast[0] = -1;
    int last = 1, idx = 0;
    for ( int it = 0; it < natom; it++ )
    {
        if ( resList[it] > last )
        {
            last++;
            resLast[++idx] = it - 1;
        }
    }
    resLast[iresl] = natom - 1;

    // iterate over residues
    while ( true )
    {
        if ( iresf <= iresl )
        {
            mid = (iresl + iresf) / 2;
            if ( resLast[mid - 1] <= atomidx && resLast[mid] >= atomidx )
            {
                a1 = resLast[mid - 1] + 1;
                a2 = resLast[mid];

                listp = new int[4];
                listp[0] = 2 - 1;
                listp[1] = 4 - 1;
                listp[2] = a1;
                listp[3] = a2;

                break;
            }
            else if ( resLast[mid] <= atomidx )
            {
                iresf = mid + 1;
            }
            else
            {
                iresl = mid - 1;
            }
        }
        else
        {
            cerr << "Error with residue bounds, check line " << __LINE__
                    << "of file " << __FILE__ << endl;
            exit(-16);
        }
    }

    /* */
    delete[] resList;
    delete[] resLast;
}

void List_Moves::gnbndl(int atomidx)
{
    int nb, nt, ni, np, ne;

    nb = nbtf(IABNDP[atomidx], moveBondUpdate.at(nMoveTypes).bonds);
    nt = nbtf(IATHTP[atomidx], moveBondUpdate.at(nMoveTypes).angles);
    np = nbtf(IATHTP[atomidx], moveBondUpdate.at(nMoveTypes).dihe);
    ni = nbtf(IAIMPP[atomidx], moveBondUpdate.at(nMoveTypes).impr);

    ne = nb + nt + np + ni;

    int nmvat = nMoveAtm.at(nMoveTypes);
    if ( ne == 0 )
        moveBondList.at(nMoveTypes)[nmvat] = nullptr;
    else
        moveBondList.at(nMoveTypes)[nmvat] = new int[ne];

    ne = 0;
    if ( moveBondUpdate.at(nMoveTypes).bonds )
    {
        ne++;
        assibl(ne, nb, IABNDP[atomidx], moveBondList[nMoveTypes][nmvat]);
    }
    if ( moveBondUpdate.at(nMoveTypes).angles )
    {
        ne++;
        assibl(ne, nt, IATHTP[atomidx], moveBondList[nMoveTypes][nmvat]);
    }
    if ( moveBondUpdate.at(nMoveTypes).dihe )
    {
        ne++;
        assibl(ne, np, IATHTP[atomidx], moveBondList[nMoveTypes][nmvat]);
    }
    if ( moveBondUpdate.at(nMoveTypes).impr )
    {
        ne++;
        assibl(ne, ni, IAIMPP[atomidx], moveBondList[nMoveTypes][nmvat]);
    }
}

int List_Moves::nbtf(int* list, bool isActive) const
{
    return ((isActive) ? list[0] : 0);
}

void List_Moves::assibl(int& ne, int n, int *orig, int* dest) const
{
    int i, is;

    is = ne;
    for ( i = 1; i < n; i++ )
    {
        ne++;
        dest[ne - 1] = orig[i];
    }
    dest[is - 1] = ne;
}
std::ostream& operator<<(std::ostream& overloadStream, const List_Moves& lst)
{
    lst.toString(overloadStream);

    return overloadStream;
}

void List_Moves::toString(std::ostream& stream) const
{
    int nmvat = nMoveAtm.at(nMoveTypes);
    stream << "NMVAT is : \t" << nmvat << endl;
    //    stream << "nMoveTypes : \t" << nMoveTypes << endl << endl;
    //    for ( int i = 0; i < nMoveTypes; i++ )
    //    {
    //        int nmvat = nMoveAtm.at(i);
    //        stream << "Dump of moveAtomList[" << i << "] \t";
    //        stream << "NMVAT is : \t" << nmvat << endl;
    //
    //        //        int** ptr = moveAtomList.at(i);
    //        //        int* idx = ptr[nMoveAtm.at(nMoveTypes)];
    //
    //    }
}
