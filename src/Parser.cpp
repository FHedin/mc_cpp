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
#include <map>
#include <string>

#include "Parser.hpp"

using namespace std;
using namespace rapidxml;

Parser_XML::Parser_XML(const char inpfileName[], PerConditions** pbc, Ensemble** ens,
                       std::vector<Atom>& atomList, FField** ff, List_nonBonded** exlst,
                       List_Moves** mvlist, MC** simulation, bool _verbose) : verbose(_verbose)
{
    string input_title;

    string pbtype;
    double a,b,c,alpha,beta,gamma;

    string enstype;
    int natom;
    double T;

    string ffmode, fffile;
    bool is_mdbas=false;
#ifdef CHARMM_EXPERIMENTAL
    bool is_charmm=false;
#endif
    string cutMode;
    double ctoff, cton, dcut;

    IO* io = nullptr;
    string atMode;
    string symb;
    double q, lj_eps, lj_sig;
    string corname;

    string mvtyp, mvmode;
    int nsele=1;
    string smode, selec;
    vector<tuple<string,string>> seleList;
    
    int nsteps;
    double dmax_value, dmax_target;
    int dmax_each;
    int save_frequency;
    
    //--------------------------------------------------

    xmlFile = new file<>(inpfileName);
    doc = new xml_document<>();
    doc->parse<0>(xmlFile->data());

    //case insensitive strings used for node names and attribute names
    xml_node<> *inpf_root = doc->first_node("inputFile");
    //we expect only one title attribute
    input_title = inpf_root->first_attribute("title")->value();
    cout << "Simulation title is : " << input_title << endl ;

    //get nodes sons of inpf_root
    xml_node<> *son_of_root = nullptr;
    string son_name = "";
    for(son_of_root = inpf_root->first_node(); son_of_root; son_of_root = son_of_root->next_sibling())
    {
        son_name = son_of_root->name();
        cout << "RootNode " << inpf_root->name() <<  " has a son named " << son_name << endl;

        // Periodic boundaries conditions parsing goes here
        if(!son_name.compare("pbc"))
        {
            //read all attributes of node pbc from xml
            attribute_processing(son_of_root);

            //convert each of the attribute to the required value
            pbtype = val_from_attr<string>("pbc_type");
            a = val_from_attr<double>("a");
            b = val_from_attr<double>("b");
            c = val_from_attr<double>("c");
            alpha = val_from_attr<double>("alpha");
            beta = val_from_attr<double>("beta");
            gamma = val_from_attr<double>("gamma");
            // and instantiate a pbc object
            *pbc = new PerConditions(pbtype, a, b, c, alpha, beta, gamma);

            // now we clear the attributes list
            attrs_list.clear();

        }
        // Ensemble parsing goes here
        else if(!son_name.compare("ensemble"))
        {
            //read all attributes of node ensemble from xml
            attribute_processing(son_of_root);

            enstype = val_from_attr<string>("ens_type");
            natom = val_from_attr<int>("N");
            T = val_from_attr<double>("T");
            Tools::str_rm_blank_spaces(enstype);
            Tools::str_to_lower_case(enstype);
            if ( !enstype.compare("nvt") )
            {
                *ens = new Ens_NVT(natom, (*pbc)->computeVol(), T);
            }
            else
            {
                cerr << "Error : the current version only supports the 'NVT' ensemble, "
                     "with a possibly infinite volume (i.e. no pbc)." << std::endl;
                exit(-3);
            }

            // now we clear the attributes list
            attrs_list.clear();

//             xml_node<> *son_of_ens = son_of_root->first_node();
//             if (son_of_ens != nullptr)
//                 cout << "Node " << son_name << " has a son called : " << son_of_ens->name() << endl;
        }
        // forcefield parsing including type, param file, and cutoff parameters goes here
        else if(!son_name.compare("forcefield"))
        {
            //read all attributes of node forcefield from xml
            attribute_processing(son_of_root);

            ffmode = val_from_attr<string>("ff_type");
            fffile = val_from_attr<string>("ff_file");
            Tools::str_rm_blank_spaces(ffmode);
            Tools::str_to_lower_case(ffmode);
            Tools::str_rm_blank_spaces(fffile);

            is_mdbas = !ffmode.compare("mdbas");
#ifdef CHARMM_EXPERIMENTAL
            is_charmm = !ffmode.compare("charmm");
#endif

            cutMode = val_from_attr<string>("cut_type");

            Tools::str_rm_blank_spaces(cutMode);
            Tools::str_to_lower_case(cutMode);
            if ( !cutMode.compare("switch"))
            {
                ctoff = val_from_attr<double>("cutoff");
                cton = val_from_attr<double>("cuton");
                dcut = val_from_attr<double>("delta_cut");
                *ff = new FField_MDBAS(atomList, **pbc, **ens, cutMode, ctoff, cton, dcut);
            }
            else if (!cutMode.compare("full"))
            {
                *ff = new FField_MDBAS(atomList, **pbc, **ens, cutMode);
            }
            else
            {
                cerr << "Error : the current version only supports the 'switch' cut function and 'full'"
                     " for non bonded interactions : please fix your input file." << std::endl;
                exit(-50);
            }

            attrs_list.clear();

//             xml_node<> *son_of_ff = son_of_root->first_node();
//             if (son_of_ff != nullptr)
//                 cout << "Node " << son_name << " has a son called : " << son_of_ff->name() << endl;
        }
        // definition of atoms list (from file or automatically generated) goes here
        else if(!son_name.compare("atomlist"))
        {
            attribute_processing(son_of_root);

            // Atom list + coordinates reading
            atMode = val_from_attr<string>("list_mode");
            Tools::str_rm_blank_spaces(atMode);
            Tools::str_to_lower_case(atMode);
            if ( !atMode.compare("manual") )
            {
                symb = val_from_attr<string>("symbol");
                q = val_from_attr<double>("charge");
                lj_eps = val_from_attr<double>("lj_epsilon");
                lj_sig = val_from_attr<double>("lj_sigma");
                for ( int i = 0; i < natom; i++ )
                {
                    atomList.push_back(Atom(i, symb));
                    atomList.at(i).setCharge(q);
                    atomList.at(i).setEpsilon(lj_eps);
                    atomList.at(i).setSigma(lj_sig);
                }
            }
            else if ( !atMode.compare("file") )
            {
                corname = val_from_attr<string>("cor_path");
                if ( is_mdbas )
                {
                    io = new IO_MDBAS(corname, fffile, **ff, atomList, **pbc, **ens);
                }
                else
                {
                    cerr << "Error : only MDBAS forcefield and coordinates (COR) supported ." << std::endl;
                    exit(-7);
                }
            }
            else
            {
                cerr << "Error : the current version only supports "
                     "the 'manual' mode, or 'file' mode for the atomlist" << std::endl;
                exit(-4);
            }

            //build exclude list and link it to ff
            *exlst = new List_nonBonded(atomList, **ff, **pbc, **ens);
            (*ff)->setExcl(**exlst);

            attrs_list.clear();
        }
        // definition of MC move types and the associated selection(s) goes here
        else if(!son_name.compare("movelist"))
        {
            attribute_processing(son_of_root);

            //first create an empty list of mc moves
            *mvlist = new List_Moves(atomList, **ff, natom);

//             int nmoves=1;
//             for(int itmv=0; itmv<nmoves; itmv++)
//             {

            string mvtyp = val_from_attr<string>("move_type");
            string mvmode = val_from_attr<string>("move_mode");

            xml_node<> *son_of_move = son_of_root->first_node("selection");
            if (son_of_move != nullptr)
                cout << "\tNode " << son_name << " has a son called : " << son_of_move->name() << endl;

// 	    int nsele=1;
// 	    string smode, selec;
// 	    vector<tuple<string,string>> seleList;
// 	    for(int itsel=0; itsel<nsele; itsel++)
// 	    {
            attribute_processing(son_of_move);
            smode = val_from_attr<string>("sel_mode");
            selec = val_from_attr<string>("sele");
            seleList.push_back( tuple<string,string>(smode,selec) );
            attrs_list.clear();
// 	    }
            (*mvlist)->addNewMoveType(mvtyp, mvmode, smode, selec);
            //         (*mvlist)->addNewMoveType(mvtyp, mvmode, seleList);

//             }

            (*ff)->setMcmvlst(**mvlist);

	    attrs_list.clear();
	    
//             xml_node<> *son_of_move = son_of_root->first_node();
//             if (son_of_move != nullptr)
//                 cout << "\tNode " << son_name << " has a son called : " << son_of_move->name() << endl;
        }
        // type of MC simulation, and the associated parameters go here
        else if(!son_name.compare("mc"))
        {
	    attribute_processing(son_of_root);
	    
            nsteps = val_from_attr<int>("nsteps");

            dmax_value = val_from_attr<double>("dmax_value");
            dmax_target = val_from_attr<double>("dmax_target");
            dmax_each = val_from_attr<int>("dmax_each");

            save_frequency = val_from_attr<int>("save_each");

            *simulation = new MC_metropolis(atomList, **pbc, **ens, **ff, **mvlist, nsteps, save_frequency, dmax_value, dmax_target, dmax_each);

	    attrs_list.clear();
	    
//             xml_node<> *son_of_mc = son_of_root->first_node();
//             if (son_of_mc != nullptr)
//                 cout << "Node " << son_name << " has a son called : " << son_of_mc->name() << endl;
        }
        else
        {
            cerr << "Warning : unknown section \"<" << son_name << ">\" will be ignored." << endl;
        }

    }//end iteration on all sons of root node inputFile

    if (is_mdbas)
        delete io;

    delete xmlFile;
    delete doc;

}

// Parser_XML::Parser_XML(const char inpfileName[], bool _verbose) : verbose(_verbose)
// {
//     xmlFile = new file<>(inpfileName);
//     doc = new xml_document<>();
//     doc->parse<0>(xmlFile->data());
//
//     //case insensitive strings used for node names and attribute names
//     xml_node<> *root = doc->first_node("inputFile");
//
//     node_processing(root);
//
//     if ( verbose )
//         Dump();
//
//     delete xmlFile;
//     delete doc;
// }

Parser_XML::~Parser_XML()
{
}

// void Parser_XML::Dump()
// {
//     cerr << "Dump of nodes_list [node_name => number_of_attributes] : " << endl;
//     for ( auto it_nodes = nodes_list.begin(); it_nodes != nodes_list.end(); ++it_nodes )
//     {
//         cerr << it_nodes->first << " => " << it_nodes->second << endl;
//     }
//     cerr << endl;
//
//     cerr << "Dump of vals_list [node_name => value] : " << endl;
//     for ( auto it_attrs = vals_list.begin(); it_attrs != vals_list.end(); ++it_attrs )
//     {
//         cerr << it_attrs->first << " => " << it_attrs->second << endl;
//     }
//     cerr << endl;
//
//     cerr << "Dump of attrs_list [attr_name => attr_value] : " << endl;
//     for ( auto it_attrs = attrs_list.begin(); it_attrs != attrs_list.end(); ++it_attrs )
//     {
//         cerr << it_attrs->first << " => " << it_attrs->second << endl;
//     }
//     cerr << endl;
// }

// void Parser_XML::node_processing(xml_node<> *src)
// {
//     int n_attr;
//     for ( xml_node<> *node = src->first_node(); node; node = node->next_sibling() )
//     {
//         if ( node->name_size() > 0 && node->value_size() > 0 )
//             vals_list.insert(pair<string, string>(node->name(), node->value()));
//
//         n_attr = attribute_processing(node);
//
//         if ( node->name_size() > 0 && n_attr > 0 )
//             nodes_list.insert(pair<string, int>(node->name(), n_attr));
//
//         check_has_son(node);
//     }
// }

int Parser_XML::attribute_processing(xml_node<> *src)
{
    int n_attr = 0;
    for ( xml_attribute<> *attr = src->first_attribute(); attr; attr = attr->next_attribute() )
    {
        n_attr++;
        attrs_list.insert(pair<string, string>(attr->name(), attr->value()));
    }
    return n_attr;
}

// xml_node<>* Parser_XML::check_has_son(xml_node<> *src)
// {
//     xml_node<> *son = src->first_node();
//     int n_attr;
//     if ( son != 0 )
//     {
//         if ( son->name_size() > 0 && son->value_size() > 0 )
//             vals_list.insert(pair<string, string>(son->name(), son->value()));
//
//         n_attr = attribute_processing(son);
//         nodes_list.insert(pair<string, int>(son->name(), n_attr));
//         node_processing(son);
//         if ( verbose )
//             cerr << "Node " << src->name() << " has a son called : " << son->name() << endl;
//     }
//     return son;
//
// 	    xml_node<> *son_of_pbc = son_of_root->first_node();
//             if (son_of_pbc != nullptr)
//                 cout << "Node " << son_name << " has a son called : " << son_of_pbc->name() << endl;
// }




