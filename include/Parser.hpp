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

#include <iostream>       
#include <stdexcept>      
#include <string>
#include <sstream>
#include <map> 
#include <exception> 

//#include "Global_include.hpp"

#include "rapidxml.hpp"
#include "rapidxml_utils.hpp"

#include "PerConditions.hpp"

#include "Ensemble.hpp"
#include "Ens_NVT.hpp"

#ifdef NPT_EXPERIMENTAL
#include "Ens_NPT.hpp"
#endif

#include "Atom.hpp"

#include "FField.hpp"

#include "List_nonBonded.hpp"
#include "List_Moves.hpp"

#include "MC.hpp"
#include "MC_metropolis.hpp"

#ifdef SPAV_EXPERIMENTAL
#include "MC_spav.hpp"
#endif 

#ifdef CHARMM_EXPERIMENTAL
#include "IO_CHARMM.hpp"
#endif

#include "IO_MDBAS.hpp"

#include "FField_MDBAS.hpp"

// template <typename T>
// class useDefaultValue: public std::exception
// {
// private:
//     std::string paramname;
// public:
//     //ctor
// //     useDefaultValue(std::string name, T _value);
//     useDefaultValue(std::string _name):paramname(_name){};
//     
//     
//     virtual const char* what() const throw()
//     {
//         std::string message = "Note : no value in xml file for the required parameter '" + paramname + "' ; using default value ..." ;
//         return message.c_str();
//     }
// };

class Parser_XML
{
public:
    Parser_XML(const char inpfileName[], PerConditions** pbc, Ensemble** ens,
                                std::vector<Atom>& atomList, FField** ff, List_nonBonded** exlst,
                                List_Moves** mvlist, MC** simulation, bool _verbose = true);
    
    virtual ~Parser_XML();

	bool requiredBenchmark() const;

private:
    //parsing related variables
    bool verbose;
    rapidxml::file<>* xmlFile;
    rapidxml::xml_document<>* doc;
    
    std::map<std::string, std::string> attrs_list; // attribute name and value stored has string and processed later

    int attribute_processing(rapidxml::xml_node<> *src);
    
    template <typename T>
    T string_to_T(const std::string& str);
    
//     template <typename T>
//     T val_from_attr(const std::string& str);
    
    template <typename T>
    void val_from_attr(const std::string& str, T& var);
    
    //here all the parameters parsed from xml file, with default values
    std::string input_title="";
    
    std::string pbtype="none";
    double a=0.,b=0.,c=0.,alpha=0.,beta=0.,gamma=0.;
    
    std::string enstype="nvt";
    int natom=0;
    double Temp=0.;
    
    std::string ffmode="", fffile="";
    bool is_mdbas=false;
    #ifdef CHARMM_EXPERIMENTAL
    bool is_charmm=false;
    #endif
    std::string cutMode="full";
    double ctoff=12., cton=10., dcut=2.;
    
    IO* io = nullptr;
    std::string atMode="file";
    std::string symb="";
    double q=0., lj_eps=1., lj_sig=1.;
    std::string corname="";
    
    std::string mvtyp="trn", mvmode="all";
    std::string smode="all", selec="";
//     int nsele=1;
//     std::vector<std::tuple<std::string,std::string>> seleList;

	bool benchmark=false;
    
    int nsteps=10000;
    double dmax_value=0.5, dmax_target=50.0;
    int dmax_each=100;
    int save_frequency=100;
	uint64_t seed = 0;

};

template <typename T>
T Parser_XML::string_to_T(const std::string& str)
{
    std::stringstream st(str);
    T type;
    st >> type;
    return type;
}

// template <typename T>
// T Parser_XML::val_from_attr(const std::string& str)
// {
//     std::string value;
// 
//     try
//     {
//         value = attrs_list.at(str);
//     }
//     catch ( const std::out_of_range& oor )
//     {
// //         std::cerr << std::endl << "Warning : keyword " << str << " not found in the input XML file." << std::endl;
// //         if ( verbose )
// //         {
// //             std::cerr << "Out of Range error when accessing to "
// //                     "std::map<std::string,std::string> attrs_list at rank ['" << str << "'] : \n" << oor.what() << std::endl << std::endl;
// //         }
//         value = "NA"; // NA for Not Available
//     }
// 
//     T toreturn;
//     
//     // no value specified so keep default one
//     if(!value.compare("NA"))
//     {
//         throw new useDefaultValue(str);
//     }
//     else
//     {
//         toreturn = string_to_T<T>( value );
//     }
// 
//     if ( verbose )
//     {
//         std::cout << "Attribute name : " << str << "\t Value : " << toreturn << std::endl;
//     }
// 
//     return toreturn;
// }

template <typename T>
void Parser_XML::val_from_attr(const std::string& str, T& var)
{
    std::string value;
    
    try
    {
        value = attrs_list.at(str);
    }
    catch ( const std::out_of_range& oor )
    {
        //         std::cerr << std::endl << "Warning : keyword " << str << " not found in the input XML file." << std::endl;
        //         if ( verbose )
        //         {
        //             std::cerr << "Out of Range error when accessing to "
        //                     "std::map<std::string,std::string> attrs_list at rank ['" << str << "'] : \n" << oor.what() << std::endl << std::endl;
        //         }
        value = "NA"; // NA for Not Available
    }
    
//     T toreturn;
    
    // no value specified so keep default one from var
    if(!value.compare("NA"))
    {
        std::cerr << "Note : no value in xml file for the required parameter '" << str << "' ; using default value of " << var << std::endl;
        return;
    }
    else
    {
        var = string_to_T<T>( value );
    }
    
    if ( verbose )
    {
        std::cout << "Attribute name : " << str << "\t Value : " << var << std::endl;
    }
    
}

#endif	/* PARSER_H */

