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
    xmlFile = new file<>(inpfileName);
    doc = new xml_document<>();
    doc->parse<0>(xmlFile->data());
    
    //case insensitive strings used for node names and attribute names
    xml_node<> *inpf_root = doc->first_node("inputFile");
    //we expect only one title attribute
    string input_title = inpf_root->first_attribute("title")->value();
    cout << "title is : " << input_title << endl ;

    /*
     *    // box and periodic boundary conditions
    string pbtype = xmlfp->val_from_attr<string>("pbctype");
    double a = xmlfp->val_from_attr<double>("a");
    double b = xmlfp->val_from_attr<double>("b");
    double c = xmlfp->val_from_attr<double>("c");
    double alpha = xmlfp->val_from_attr<double>("alpha");
    double beta = xmlfp->val_from_attr<double>("beta");
    double gamma = xmlfp->val_from_attr<double>("gamma");
    *pbc = new PerConditions(pbtype, a, b, c, alpha, beta, gamma);

    // informations about the Ensemble
    string enstype = xmlfp->val_from_attr<string>("enstype");
    int natom = xmlfp->val_from_attr<int>("N");
    double T = xmlfp->val_from_attr<double>("T");
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

    // Forcefield parameters  + opening of coordinates files
    IO* io = nullptr;
    string ffmode = xmlfp->val_from_attr<string>("ff_mode");
    string fffile = xmlfp->val_from_attr<string>("ff_file");
    Tools::str_rm_blank_spaces(ffmode);
    Tools::str_to_lower_case(ffmode);
    Tools::str_rm_blank_spaces(fffile);

    bool is_mdbas = !ffmode.compare("mdbas");
#ifdef CHARMM_EXPERIMENTAL
    bool is_charmm = !ffmode.compare("charmm");
#endif 

    string cutMode = xmlfp->val_from_attr<string>("cut");
    
    Tools::str_rm_blank_spaces(cutMode);
    Tools::str_to_lower_case(cutMode);
    if ( !cutMode.compare("switch"))
    {
        double ctoff = xmlfp->val_from_attr<double>("cutoff");
        double cton = xmlfp->val_from_attr<double>("cuton");
        double dcut = xmlfp->val_from_attr<double>("delta_cut");
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

    // Atom list + coordinates reading
    string atMode = xmlfp->val_from_attr<string>("at_list");
    Tools::str_rm_blank_spaces(atMode);
    Tools::str_to_lower_case(atMode);
    if ( !atMode.compare("repeat") )
    {
        string symb = xmlfp->val_from_attr<string>("symbol");
        double q = xmlfp->val_from_attr<double>("charge");
        double lj_eps = xmlfp->val_from_attr<double>("lj_epsilon");
        double lj_sig = xmlfp->val_from_attr<double>("lj_sigma");
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
        string corname = xmlfp->val_from_attr<string>("at_file");
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
                "the 'repeat' mode, or 'file' mode for the atomlist" << std::endl;
        exit(-4);
    }

    //build exclude list and link it to ff
    *exlst = new List_nonBonded(atomList, **ff, **pbc, **ens);
    (*ff)->setExcl(**exlst);

    // selection list and moves list
    // <move move_type="TRN" move_mode="ATOM" >
    *mvlist = new List_Moves(atomList, **ff, (*ens)->getN());

    int nmoves=1;
    for(int itmv=0;itmv<nmoves;itmv++)
    {
        string mvtyp = xmlfp->val_from_attr<string>("move_type");
        string mvmode = xmlfp->val_from_attr<string>("move_mode");
        
        int nsele=1;
        string smode, selec;
        vector<tuple<string,string>> seleList;
        for(int itsel=0;itsel<nsele;itsel++)
        {
            smode = xmlfp->val_from_attr<string>("sel_mode");
            selec = xmlfp->val_from_node<string>("selection");
            seleList.push_back( tuple<string,string>(smode,selec) );
        }
        (*mvlist)->addNewMoveType(mvtyp, mvmode, smode, selec);
//         (*mvlist)->addNewMoveType(mvtyp, mvmode, seleList);
    }
    
    (*ff)->setMcmvlst(**mvlist);

    int nsteps = xmlfp->val_from_attr<int>("nsteps");
    
    double dmax_value = xmlfp->val_from_attr<double>("dmax_value");
    double dmax_target = xmlfp->val_from_attr<double>("dmax_target");
    int dmax_each = xmlfp->val_from_attr<int>("dmax_each");
    
    int save_frequency = xmlfp->val_from_attr<int>("save_each");

    *simulation = new MC_metropolis(atomList, **pbc, **ens, **ff, **mvlist, nsteps, save_frequency, dmax_value, dmax_target, dmax_each);

    delete io;
     */
    
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
            xml_node<> *son_of_pbc = son_of_root->first_node();
            if (son_of_pbc != nullptr)
                cout << "Node " << son_name << " has a son called : " << son_of_pbc->name() << endl;
        }
        // Ensemble parsing goes here
        else if(!son_name.compare("ensemble"))
        {
            xml_node<> *son_of_ens = son_of_root->first_node();
            if (son_of_ens != nullptr)
                cout << "Node " << son_name << " has a son called : " << son_of_ens->name() << endl;
        }
        // forcefield parsing including type, param file, and cutoff parameters goes here
        else if(!son_name.compare("ff"))
        {
            xml_node<> *son_of_ff = son_of_root->first_node();
            if (son_of_ff != nullptr)
                cout << "Node " << son_name << " has a son called : " << son_of_ff->name() << endl;
        }
        // definition of atoms list (from file or automatically generated) goes here
        else if(!son_name.compare("atomlist"))
        {
            xml_node<> *son_of_atlist = son_of_root->first_node();
            if (son_of_atlist != nullptr)
                cout << "Node " << son_name << " has a son called : " << son_of_atlist->name() << endl;
        }
        // definition of move types and the associated selection(s) goes here
        else if(!son_name.compare("move"))
        {
            xml_node<> *son_of_move = son_of_root->first_node();
            if (son_of_move != nullptr)
                cout << "\tNode " << son_name << " has a son called : " << son_of_move->name() << endl;
        }
        // type of MC simulation, and the associated parameters go here 
        else if(!son_name.compare("mc"))
        {
            xml_node<> *son_of_mc = son_of_root->first_node();
            if (son_of_mc != nullptr)
            cout << "Node " << son_name << " has a son called : " << son_of_mc->name() << endl;
        }
        else
        {
            cerr << "Warning : unknown section \"<" << son_name << ">\" will be ignored." << endl;
        }
        
    }//end iteration on all sons of root node inputFile
    
    delete xmlFile;
    delete doc;
    
}

Parser_XML::Parser_XML(const char inpfileName[], bool _verbose) : verbose(_verbose)
{
    xmlFile = new file<>(inpfileName);
    doc = new xml_document<>();
    doc->parse<0>(xmlFile->data());
    
    //case insensitive strings used for node names and attribute names
    xml_node<> *root = doc->first_node("inputFile");

    node_processing(root);

    if ( verbose )
        Dump();

    delete xmlFile;
    delete doc;
}

Parser_XML::~Parser_XML()
{
}

void Parser_XML::Dump()
{
    cerr << "Dump of nodes_list [node_name => number_of_attributes] : " << endl;
    for ( auto it_nodes = nodes_list.begin(); it_nodes != nodes_list.end(); ++it_nodes )
    {
        cerr << it_nodes->first << " => " << it_nodes->second << endl;
    }
    cerr << endl;

    cerr << "Dump of vals_list [node_name => value] : " << endl;
    for ( auto it_attrs = vals_list.begin(); it_attrs != vals_list.end(); ++it_attrs )
    {
        cerr << it_attrs->first << " => " << it_attrs->second << endl;
    }
    cerr << endl;

    cerr << "Dump of attrs_list [attr_name => attr_value] : " << endl;
    for ( auto it_attrs = attrs_list.begin(); it_attrs != attrs_list.end(); ++it_attrs )
    {
        cerr << it_attrs->first << " => " << it_attrs->second << endl;
    }
    cerr << endl;
}

void Parser_XML::node_processing(xml_node<> *src)
{
    int n_attr;
    for ( xml_node<> *node = src->first_node(); node; node = node->next_sibling() )
    {
        if ( node->name_size() > 0 && node->value_size() > 0 )
            vals_list.insert(pair<string, string>(node->name(), node->value()));

        n_attr = attribute_processing(node);

        if ( node->name_size() > 0 && n_attr > 0 )
            nodes_list.insert(pair<string, int>(node->name(), n_attr));

        check_has_son(node);
    }
}

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

xml_node<>* Parser_XML::check_has_son(xml_node<> *src)
{
    xml_node<> *son = src->first_node();
    int n_attr;
    if ( son != 0 )
    {
        if ( son->name_size() > 0 && son->value_size() > 0 )
            vals_list.insert(pair<string, string>(son->name(), son->value()));

        n_attr = attribute_processing(son);
        nodes_list.insert(pair<string, int>(son->name(), n_attr));
        node_processing(son);
        if ( verbose )
            cerr << "Node " << src->name() << " has a son called : " << son->name() << endl;
    }
    return son;
}



