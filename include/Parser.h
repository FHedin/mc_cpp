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

#ifndef PARSER_H
#define	PARSER_H

#include <map>
#include <iostream>       
#include <stdexcept>      
#include <string>
#include <sstream>

#include "rapidxml-1.13/rapidxml.hpp"
#include "rapidxml-1.13/rapidxml_utils.hpp"

class Parser_XML
{
    
public:
    Parser_XML(const char inpfileName[], bool verbose=false);
    ~Parser_XML();
    
    template <typename T>
    T val_from_attr(const std::string& str);
    
private:
    void Dump();
    void node_processing(rapidxml::xml_node<> *src);
    int attribute_processing(rapidxml::xml_node<> *src);
    void check_has_son(rapidxml::xml_node<> *src);
    
    std::map<std::string,int>           nodes_list;   // the int is the number of attributes for a node
    std::map<std::string,std::string>   attrs_list;   // attribute name and value stored has string and processed later
    
    template <typename T>
    T string_to_T(const std::string& str);

};

template <typename T>
T Parser_XML::string_to_T(const std::string& str)
{
    std::stringstream st(str);
    T type;
    st >> type;
    return type;
}

template <typename T>
T Parser_XML::val_from_attr(const std::string& str)
{
    std::string value;
    
    try
    {
        value = attrs_list.at(str);
    }
    catch (const std::out_of_range& oor)
    {
        std::cerr << "Out of Range error when accessing to std::map<std::string,std::string> attrs_list at rank ['" << str << "'] : \n" << oor.what() << std::endl;
    }
    
    return string_to_T<T>(value);
}

#endif	/* PARSER_H */

