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

#include <iostream>
#include <sstream>
#include <algorithm> // for std::fill

#include "Selection.h"
#include "Tools.h"

using namespace std;

Selection::Selection(std::string _selectMode, std::string _selectionString,
                     std::vector<Atom>& _at_List, const int& _natom)
: at_List(_at_List), natom(_natom)
{
    selectionString = _selectionString;
    selection.resize(natom, 0);

    getMode(_selectMode);

    select_main();
}

Selection::~Selection()
{
}

void Selection::getMode(std::string _selectMode)
{
    Tools::str_rm_blank_spaces(_selectMode);
    Tools::str_to_lower_case(_selectMode);

    if ( !_selectMode.compare("all") )
        selectMode = SELEMODE::ALL;
    else if ( !_selectMode.compare("none") )
        selectMode = SELEMODE::NONE;
    else if ( !_selectMode.compare("residue_id") )
        selectMode = SELEMODE::RESIDUE_ID;
    else if ( !_selectMode.compare("residue_name") )
        selectMode = SELEMODE::RESIDUE_NAME;
        //    else if ( !_selectMode.compare("segment_id") )
        //        selectMode = SEGMENT_ID;
        //    else if ( !_selectMode.compare("segment_name") )
        //        selectMode = SEGMENT_NAME;
    else if ( !_selectMode.compare("atom_idx") )
        selectMode = SELEMODE::ATOM_IDX;
    else
    {
        cerr << "Error : the following 'sel_mode' is unknown : " << _selectMode << endl;
        exit(-17);
    }
}

void Selection::select_main()
{
    switch ( this->selectMode )
    {
        case SELEMODE::ALL:
            fill(selection.begin(), selection.begin() + natom, 1);
            break;
        case SELEMODE::NONE:
            break;
        case SELEMODE::RESIDUE_ID:
            select_resId();
            break;
        case SELEMODE::RESIDUE_NAME:
            select_resName();
            break;
            //        case SEGMENT_ID:
            //            select_segId();
            //            break;
            //        case SEGMENT_NAME:
            //            select_segName();
            //            break;
        case SELEMODE::ATOM_IDX:
            select_atomIdx();
            break;
        default:
            break;
    }
}

void Selection::select_resId()
{
    int num;
    istringstream ss(selectionString);
    ss >> num;
    for ( int it = 0; it < natom; it++ )
    {
        const int& residx = at_List.at(it).getResidue_id_global();
        if ( residx == num )
        {
            selection.at(it) = 1;
        }
    }
}

void Selection::select_resName()
{
    for ( int it = 0; it < natom; it++ )
    {
        const string& str = at_List.at(it).getRes_label();
        if ( !selectionString.compare(str) )
        {
            selection.at(it) = 1;
        }
    }
}

//void Selection::select_segId()
//{
//
//}
//
//void Selection::select_segName()
//{
//
//}

void Selection::select_atomIdx()
{
    int num;
    istringstream ss(selectionString);
    ss >> num;
    for ( int it = 0; it < natom; it++ )
    {
        const int& idx = at_List.at(it).getID();
        if ( idx == num )
        {
            selection.at(it) = 1;
        }
    }
}

const std::vector<int>& Selection::getSelection() const
{
    return selection;
}

