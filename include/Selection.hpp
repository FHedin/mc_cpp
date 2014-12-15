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

#ifndef SELECTION_H
#define	SELECTION_H

#include <string>
#include <vector>
#include <iostream>

#include "Global_include.hpp"

#include "Atom.hpp"

class Selection
{
    friend std::ostream& operator<<( std::ostream& overloadStream, const Selection& sele );

public:
    Selection(std::string _selectMode, std::string _selectionString, std::vector<Atom>& _at_List, const int& _natom);
    virtual ~Selection();

    const std::vector<int>& getSelection() const;

private:

    enum SELEMODE
    {
        ALL,
        NONE,
        RESIDUE_ID,
        RESIDUE_NAME,
        //    SEGMENT_ID,
        //    SEGMENT_NAME,
        ATOM_IDX
    };

    std::string selectionString;
    std::vector<Atom>& at_List;
    const int& natom;

    SELEMODE selectMode;
    std::vector<int> selection;

    void getMode(std::string _selectMode);

    void select_main();
    void select_resId();
    void select_resName();
    //    void select_segId();
    //    void select_segName();
    void select_atomIdx();

    void toString(std::ostream& stream) const;

};

#endif	/* SELECTION_H */

