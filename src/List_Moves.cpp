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
#include "Tools.h"
#include "Selection.h"

const int List_Moves::MCMBND = 10;
const int List_Moves::MCMTHT = 100;
const int List_Moves::MCMPHI = 200;
const int List_Moves::MCMIMP = 100;

const int List_Moves::MMVTYP = 50;

using namespace std;

List_Moves::List_Moves(std::vector<Atom>& _at_List, FField& _ff, int _natom)
: at_List(_at_List), ff(_ff)
{
    natom = _natom;

    this->nMoveTypes = 0;

    this->nMoveAtm.resize(MMVTYP, 0);
    this->moveTypeList.resize(MMVTYP);
    this->moveModeList.resize(MMVTYP);

    this->moveBondList.resize(MMVTYP, nullptr);
    this->moveBondUpdate.resize(MMVTYP);

    this->moveAtomList.resize(MMVTYP, nullptr);

    this->moveLimitsList.resize(MMVTYP, 0.0);

    this->movePivotList.resize(MMVTYP, nullptr);
}

List_Moves::~List_Moves()
{
    int* tmp = nullptr;
    int i = 0, j = 0;

    //    cout << "Hello from Dtor List_Moves::~List_Moves()" << endl;

    for ( i = 0; i < MMVTYP; i++ )
    {
        if ( moveBondList[i] != nullptr )
        {
            //            cout << "Dtor is cleaning  moveBondList[" << i << "]" << endl;
            for ( j = 0; j < nMoveAtm[i]; j++ )
            {
                tmp = moveBondList[i][j];
                if ( tmp != nullptr )
                {
                    //                    cout << "Dtor is cleaning  moveBondList[" << i << "][" << j << "]" << endl;
                    delete[] moveBondList[i][j];
                }
            }
            delete[] moveBondList[i];
        }

        if ( moveAtomList[i] != nullptr )
        {
            //            cout << "Dtor is cleaning  moveAtomList[" << i << "]" << endl;
            for ( j = 0; j < nMoveAtm[i]; j++ )
            {
                tmp = moveAtomList[i][j];
                if ( tmp != nullptr )
                {
                    //                    cout << "Dtor is cleaning  moveAtomList[" << i << "][" << j << "]" << endl;
                    delete[] moveAtomList[i][j];
                }
            }
            delete[] moveAtomList[i];
        }

        if ( movePivotList[i] != nullptr )
        {
            //            cout << "Dtor is cleaning  movePivotList[" << i << "]" << endl;
            for ( j = 0; j < nMoveAtm[i]; j++ )
            {
                tmp = movePivotList[i][j];
                if ( tmp != nullptr )
                {
                    //                    cout << "Dtor is cleaning  movePivotList[" << i << "][" << j << "]" << endl;
                    delete[] movePivotList[i][j];
                }
            }
            delete[] movePivotList[i];
        }
    }
}

void List_Moves::addNewMoveType(string mvtypName, string modeName, string selMode, string selName)
{
    bool success = true;

    Tools::str_rm_blank_spaces(mvtypName);
    Tools::str_to_lower_case(mvtypName);

    if ( !mvtypName.compare("trn") ) // charmm MVTYPE 1
    {
        moveTypeList.at(nMoveTypes) = TRN;
        success = NewMove_TRN_ROT(modeName, selMode, selName);
        /* ... */
        nMoveTypes++;
    }
    else if ( !mvtypName.compare("rot") ) // MVTYPE 2
    {
        moveTypeList.at(nMoveTypes) = ROT;
        success = NewMove_TRN_ROT(modeName, selMode, selName);
        /* ... */
        nMoveTypes++;
    }
        //    else if ( !mvtypName.compare("tors") ) // MVTYPE 4
        //    {
        //        moveTypeList.at(nMoveTypes) = TORS;
        //        success = NewMove_TORS(modeName, selMode, selName);
        //        /* ... */
        //        nMoveTypes++;
        //    }
    else
    {
        cout << "Warning : " << mvtypName << " is not a valid type of move ; skipping ...";
    }

    if ( !success )
    {
        nMoveTypes--;
    }
}

bool List_Moves::NewMove_TRN_ROT(string modeName, string selMode, string selName)
{
    Selection selec(selMode, selName, at_List, natom);
    const vector<int>& seleList = selec.getSelection();

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
        if ( moveTypeList.at(nMoveTypes) == ROT )
        {
            cerr << "Error : rotation with single atom mode is useless ; please fix the input file" << endl;
            exit(-20);
        }
        moveModeList.at(nMoveTypes) = ATOM;
        BOND_UPDATE bndulTmp = {true, true, true, true};
        moveBondUpdate.at(nMoveTypes) = bndulTmp;
        makeBondList();
    }
    else
    {
        cout << "Warning : the following move mode is not available : "
                << modeName << " skipping ..." << endl;

        return false;
    }

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
    // only centre of mass available now, or one atom selection (automatic if ROT RESIDUE)
    bool cofmas = false;

    // pointer to list of moving atoms
    int* listp = nullptr;
    int* tempp = nullptr;
    int ia1, ia2;

    if ( moveModeList.at(nMoveTypes) == ALL ) //mode 2
    {
        if ( moveTypeList.at(nMoveTypes) == TRN ) // mvtyp 1
        {
            //            cout << listp << endl;
            makeMoveList(listp, natom, seleList);
            //            cout << listp << endl;
            nMoveAtm.at(nMoveTypes) = 1;
        }
        else if ( moveTypeList.at(nMoveTypes) == ROT ) // mvtyp 2
        {
            //            cout << listp << endl;
            makeMoveList(listp, natom, seleList);
            //            cout << listp << endl;
            nMoveAtm.at(nMoveTypes) = 0;
            /* future second selection work here */
            cofmas = (nMoveAtm.at(nMoveTypes) == 0);
            if ( cofmas )
                nMoveAtm.at(nMoveTypes) = 1;
        }
    }

    // Now allocate space for the atom list according nMoveAtm size
    //    cout << "On line " << __LINE__ << '\t' << *this;
    moveAtomList.at(nMoveTypes) = new int *[ nMoveAtm.at(nMoveTypes) ];
    fillNull(moveAtomList.at(nMoveTypes), nMoveAtm.at(nMoveTypes));

    if ( moveTypeList.at(nMoveTypes) == ROT ) //mvtyp 2
    {
        movePivotList.at(nMoveTypes) = new int *[nMoveAtm.at(nMoveTypes)];
        fillNull(movePivotList.at(nMoveTypes), nMoveAtm.at(nMoveTypes));
    }

    if ( moveModeList.at(nMoveTypes) == ATOM ) //mode 4
    {
        //        cout << "On line " << __LINE__ << '\t' << *this;
        moveBondList.at(nMoveTypes) = new int *[ nMoveAtm.at(nMoveTypes) ];
        fillNull(moveBondList.at(nMoveTypes), nMoveAtm.at(nMoveTypes));
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
            if ( moveModeList.at(nMoveTypes) == RESIDUE ) //mode 1 by residue
            {
                gtrsfl(tempp, i, ia1, ia2, natom);
                moveAtomList.at(nMoveTypes)[nMoveAtm.at(nMoveTypes)] = tempp; // IMVNGP%A(NMVATM)%A => TEMPP%A
            }
            else if ( moveModeList.at(nMoveTypes) == ALL ) // mode 2 by all
            {
                //                cout << nMoveTypes << '\t' << nMoveAtm.at(nMoveTypes) << '\t' << moveAtomList.at(nMoveTypes)[nMoveAtm.at(nMoveTypes)] << '\t' << listp << endl;
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
                tempp = new int[1];

                if ( cofmas )
                    tempp[0] = -1;
                else
                    tempp[0] = i;

                movePivotList.at(nMoveTypes)[nMoveAtm.at(nMoveTypes)] = tempp;
            }

            nMoveAtm.at(nMoveTypes)++;

            if ( (moveModeList.at(nMoveTypes) == ALL && moveTypeList.at(nMoveTypes) == TRN)
                 || cofmas ) // mode_2_all && mvtype_1_trn
            {
                break;
            }

        }// long if ...
    }// long for ...

    /* ... */

    if ( moveModeList.at(nMoveTypes) == ATOM )
        freeBondList();

    return true;
}

bool List_Moves::NewMove_TORS(string modeName, string selMode, string selName)
{
    return false;
}

void List_Moves::makeBondList()
{
    int* tmpPtr = nullptr;

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
        //        cout << "makeBondList() works on bond " << it << endl;
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
        //        cout << "makeBondList() works on angle " << it << endl;
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
        //        cout << "makeBondList() works on dihedral " << it << endl;
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
        //        cout << "makeBondList() works on improper " << it << endl;
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
        cerr << "iic = " << iic << "\t n = " << n << "\t size = " << size << endl;
        cerr << "From file '" << __FILE__ << "'." << endl;
        exit(-14);
    }
    ilist[n] = iic;
}

void List_Moves::freeBondList()
{
    for ( int it = 0; it < natom; it++ )
    {
        delete[] IABNDP[it];
        delete[] IATHTP[it];
        delete[] IAPHIP[it];
        delete[] IAIMPP[it];
    }
}

void List_Moves::makeMoveList(int*& list, int natom, const vector<int>& sele) const
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
    list[0] = 2 - 1;
    list[1] = n - 1;

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

void List_Moves::gtrsfl(int*& listp, int atomidx, int& a1, int& a2, int natom) const
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

    //    cout << "From gnbndl : atomid = " << atomidx << " \t ne = " << ne << endl;
    //    cout << "Booleans : " << moveBondUpdate.at(nMoveTypes).bonds << moveBondUpdate.at(nMoveTypes).angles;
    //    cout << moveBondUpdate.at(nMoveTypes).dihe << moveBondUpdate.at(nMoveTypes).impr << endl;

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
    return ((isActive) ? list[0] + 1 : 0);
}

void List_Moves::assibl(int& ne, int n, int *orig, int* dest) const
{
    int i, is;
    //        cout << "From assibl : ne = " << ne << "\t n = " << n << endl;
    is = ne;
    for ( i = 1; i < n; i++ )
    {
        ne++;
        dest[ne - 1] = orig[i];
    }
    dest[is - 1] = ne - 1;
}

std::ostream& operator<<(std::ostream& overloadStream, const List_Moves& lst)
{
    lst.toString(overloadStream);

    return overloadStream;
}

void List_Moves::toString(std::ostream& stream) const
{
    int** ptr;
    stream << "nMoveTypes : \t" << nMoveTypes << endl << endl;
    for ( int i = 0; i < nMoveTypes; i++ )
    {
        int nmvat = nMoveAtm.at(i);
        stream << "NMVAT is : \t" << nmvat << endl << endl;

        ptr = moveAtomList.at(i);
        if ( ptr != nullptr )
        {
            /* Dump of moveAtomList, a vector of pointers to pointers of int containing
             * a list of moving atoms */
            stream << "Dump of moveAtomList[" << i << "]" << endl;

            for ( int j = 0; j < nmvat; j++ )
            {
                int* idx = ptr[j];
                if ( idx != nullptr )
                {
                    int ng = idx[0];
                    int endng = ng + 2;

                    stream << "For nmvat " << j << " : \t NG = " << ng << endl;
                    for ( int it1 = 1; it1 <= ng; it1++ )
                    {
                        int nn = idx[it1];
                        //                stream << "For mvgroup " << nn << endl;
                        for ( int it2 = endng; it2 <= nn; it2 += 2 )
                        {
                            int iaf = idx[it2 - 1];
                            int ial = idx[it2];
                            stream << "First and Last atom to translate/rotate are : " << iaf << '\t' << ial << endl;
                        }
                        endng = nn + 2;
                    }
                } // idx nullptr test
                stream << endl;
            } // loop nmvat
        } // end of moveAtomList dump

        ptr = moveBondList.at(i);
        if ( ptr != nullptr )
        {
            /* Dump of moveBondList a vector of pointers to pointers of int containing
             * a list of bonds/angles/dihedrals/impropers affected by a move
             */
            bool bbond = moveBondUpdate.at(i).bonds;
            bool bangl = moveBondUpdate.at(i).angles;
            bool bdihe = moveBondUpdate.at(i).dihe;
            bool bimpr = moveBondUpdate.at(i).impr;

            stream << "Dump of moveBondList[" << i << "]" << endl;
            stream << "Booleans (bonds,angles,dihe,impr) are : " << bbond << '\t' << bangl << '\t' << bdihe << '\t' << bimpr << endl;
            for ( int j = 0; j < nmvat; j++ )
            {
                int* idx = ptr[j];
                if ( idx != nullptr )
                {
                    int current_idx = 0, end_list;

                    stream << "For nmvat " << j << " : " << endl;

                    if ( bbond )
                    {
                        stream << "Bonds List : " << '\t';
                        end_list = idx[current_idx++];
                        for ( int it = current_idx; it <= end_list; it++ )
                        {
                            stream << idx[it] << '\t';
                        }
                        stream << endl;
                        current_idx = end_list + 1;
                    } // bbond

                    if ( bangl )
                    {
                        stream << "Angles List : " << '\t';
                        end_list = idx[current_idx++];
                        for ( int it = current_idx; it <= end_list; it++ )
                        {
                            stream << idx[it] << '\t';
                        }
                        stream << endl;
                        current_idx = end_list + 1;
                    } // bangl

                    if ( bdihe )
                    {
                        stream << "Dihedrals List : " << '\t';
                        end_list = idx[current_idx++];
                        for ( int it = current_idx; it <= end_list; it++ )
                        {
                            stream << idx[it] << '\t';
                        }
                        stream << endl;
                        current_idx = end_list + 1;
                    } // bdihe

                    if ( bimpr )
                    {
                        stream << "Impropers List : " << '\t';
                        end_list = idx[current_idx++];
                        for ( int it = current_idx; it <= end_list; it++ )
                        {
                            stream << idx[it] << '\t';
                        }
                        stream << endl;
                        current_idx = end_list + 1;
                    } // bimpr

                    stream << endl;
                } // idx nullptr test
            }// loop nmvat

        } // end of moveBondList dump

        ptr = movePivotList.at(i);
        if ( ptr != nullptr )
        {
            stream << "Dump of movePivotList[" << i << "]" << endl;
            for ( int j = 0; j < nmvat; j++ )
            {
                int* idx = ptr[j];
                if ( idx != nullptr )
                {
                    stream << "For nmvat " << j << " : " << endl;
                    stream << "Pivot is : " << idx[0] << endl;
                }

            } // loop nmvat

        } // end of movePivotList dump


    } // end of loop on nmvtyp



}// end of to string function

void List_Moves::fillNull(int** array, int size) const
{
    for ( int i = 0; i < size; i++ )
        array[i] = nullptr;
}

