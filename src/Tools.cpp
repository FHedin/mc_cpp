/*
 *  mc_cpp : A basic Monte Carlo simulations software.
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
#include <algorithm>

#include "Tools.h"

//Tools::Tools() {
//}

//Tools::~Tools() {
//}

//removes white spaces from input string
void Tools::str_rm_blank_spaces(std::string& str)
{
    
    str.erase(std::remove_if(str.begin(), str.end(),std::ptr_fun(isspace)), str.end());
}

// for comparisons, be sure that the string str is in low caps
void Tools::str_to_lower_case(std::string& str)
{
    std::transform(str.begin(), str.end(), str.begin(), std::ptr_fun(tolower));
}