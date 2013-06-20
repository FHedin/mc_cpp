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
    Parser_XML(const char inpfileName[], bool verbose = false);
    ~Parser_XML();

    template <typename T>
    T val_from_attr(const std::string& str, bool verbose = false);

    template <typename T>
    T val_from_node(const std::string& str, bool verbose = false);

private:
    bool verbose;
    rapidxml::file<>* xmlFile;
    rapidxml::xml_document<>* doc;

    std::map<std::string, int> nodes_list; // the int is the number of attributes for a node
    std::map<std::string, std::string> vals_list;
    std::map<std::string, std::string> attrs_list; // attribute name and value stored has string and processed later

    void Dump();
    void node_processing(rapidxml::xml_node<> *src);
    int attribute_processing(rapidxml::xml_node<> *src);
    rapidxml::xml_node<>* check_has_son(rapidxml::xml_node<> *src);

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
T Parser_XML::val_from_attr(const std::string& str, bool verbose)
{
    std::string value;

    try
    {
        value = attrs_list.at(str);
    }
    catch ( const std::out_of_range& oor )
    {
        std::cerr << "Warning : keyword " << str << " not found in the input XML file." << std::endl;
        if ( verbose )
        {
            std::cerr << "Out of Range error when accessing to "
                    "std::map<std::string,std::string> attrs_list at rank ['" << str << "'] : \n" << oor.what() << std::endl;
        }
        //        exit(-21);
    }

    T toreturn = string_to_T<T>( value );

    if ( verbose )
    {
        std::cout << "Attribute name : " << str << "\t Value : " << toreturn << std::endl;
    }

    return toreturn;
}

template <typename T>
T Parser_XML::val_from_node(const std::string& str, bool verbose)
{
    std::string value;

    try
    {
        value = vals_list.at(str);
    }
    catch ( const std::out_of_range& oor )
    {
        std::cerr << "Warning : node " << str << " not found in the input XML file." << std::endl;
        if ( verbose )
        {
            std::cerr << "Out of Range error when accessing to "
                    "std::map<std::string,std::string> vals_list at rank ['" << str << "'] : \n" << oor.what() << std::endl;
        }
        //        exit(-22);
    }

    T toreturn = string_to_T<T>( value );

    if ( verbose )
    {
        std::cout << "Value : " << toreturn << std::endl;
    }

    return toreturn;
}

#endif	/* PARSER_H */

