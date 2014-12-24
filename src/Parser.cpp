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
#include <exception>

#include "Parser.hpp"

using namespace std;
using namespace rapidxml;

Parser_XML::Parser_XML(const char inpfileName[], PerConditions** pbc, Ensemble** ens,
                       std::vector<Atom>& atomList, FField** ff, List_nonBonded** exlst,
                       List_Moves** mvlist, MC** simulation, bool _verbose) : verbose(_verbose)
{
    //--------------------------------------------------

    xmlFile = new file<char>(inpfileName);
    doc = new xml_document<char>();
    
    try{
        doc->parse<0>(xmlFile->data());
    }catch(parse_error e)
    {
        cerr << "Error detected when parsing file "  << inpfileName << endl;
        cerr << "What : " << e.what() << endl;
        cerr << "Where : " << e.where<char*>() << endl;
    }

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
        //cout << "RootNode " << inpf_root->name() <<  " has a son named " << son_name << endl;

        // Periodic boundaries conditions parsing goes here
        if(!son_name.compare("pbc"))
        {
            //read all attributes of node pbc from xml
            attribute_processing(son_of_root);

            //convert each of the attribute to the required value
            val_from_attr<string>("pbc_type",pbtype);
            // if type none not used we need lattice parameters
            if(pbtype.compare("none")!=0)
            {
                val_from_attr<double>("a",a);
                val_from_attr<double>("b",b);
                val_from_attr<double>("c",c);
                val_from_attr<double>("alpha",alpha);
                val_from_attr<double>("beta",beta);
                val_from_attr<double>("gamma",gamma);
            }
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

            val_from_attr<string>("ens_type",enstype);
            val_from_attr<int>("N",natom);
            val_from_attr<double>("T",Temp);
            Tools::str_rm_blank_spaces(enstype);
            Tools::str_to_lower_case(enstype);
            if ( !enstype.compare("nvt") )
            {
                *ens = new Ens_NVT(natom, (*pbc)->computeVol(), Temp);
            }
            else
            {
                cerr << "Error : the current version only supports the 'NVT' ensemble, "
                     "with a possibly infinite volume (i.e. no pbc)." << std::endl;
                exit(-3);
            }

            // now we clear the attributes list
            attrs_list.clear();

        }
        // forcefield parsing including type, param file, and cutoff parameters goes here
        else if(!son_name.compare("forcefield"))
        {
            //read all attributes of node forcefield from xml
            attribute_processing(son_of_root);

            val_from_attr<string>("ff_type",ffmode);
            val_from_attr<string>("ff_file",fffile);
            Tools::str_rm_blank_spaces(ffmode);
            Tools::str_to_lower_case(ffmode);
            Tools::str_rm_blank_spaces(fffile);

            is_mdbas = !ffmode.compare("mdbas");
#ifdef CHARMM_EXPERIMENTAL
            is_charmm = !ffmode.compare("charmm");
#endif

            val_from_attr<string>("cut_type",cutMode);

            Tools::str_rm_blank_spaces(cutMode);
            Tools::str_to_lower_case(cutMode);
            if ( !cutMode.compare("switch"))
            {
                val_from_attr<double>("cutoff",ctoff);
                val_from_attr<double>("cuton",cton);
                val_from_attr<double>("delta_cut",dcut);
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

        }
        // definition of atoms list (from file or automatically generated) goes here
        else if(!son_name.compare("atomlist"))
        {
            attribute_processing(son_of_root);

            // Atom list + coordinates reading
            val_from_attr<string>("list_mode",atMode);
            Tools::str_rm_blank_spaces(atMode);
            Tools::str_to_lower_case(atMode);
            if ( !atMode.compare("manual") )
            {
                val_from_attr<string>("symbol",symb);
                val_from_attr<double>("charge",q);
                val_from_attr<double>("lj_epsilon",lj_eps);
                val_from_attr<double>("lj_sigma",lj_sig);
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
                val_from_attr<string>("cor_path",corname);
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
            //attribute_processing(son_of_root);

            //first create an empty list of mc moves
            *mvlist = new List_Moves(atomList, **ff, natom);
            
            xml_node<> *moves = nullptr;
            for(moves=son_of_root->first_node("move");moves;moves=moves->next_sibling("move"))
            {
                cout << "\tNode " << son_name << " has a son called : " << moves->name() << endl;
                attribute_processing(moves);
                val_from_attr<string>("move_type",mvtyp);
                val_from_attr<string>("move_mode",mvmode);
                val_from_attr<string>("sel_mode",smode);
                // if sele mode all or none we don't need the xtra selection attribute
                if(smode.compare("none")!=0 && smode.compare("all")!=0)
                    val_from_attr<string>("sele",selec);
                (*mvlist)->addNewMoveType(mvtyp, mvmode, smode, selec);
                attrs_list.clear();
            }
            
            (*ff)->setMcmvlst(**mvlist);
            
//             int nmoves=1;
//             for(int itmv=0; itmv<nmoves; itmv++)
//             {

//             val_from_attr<string>("move_type",mvtyp);
//             val_from_attr<string>("move_mode",mvmode);

//             xml_node<> *son_of_move = son_of_root->first_node("selection");
//             if (son_of_move != nullptr)
//                 cout << "\tNode " << son_name << " has a son called : " << son_of_move->name() << endl;

// 	    int nsele=1;
// 	    string smode, selec;
// 	    vector<tuple<string,string>> seleList;
// 	    for(int itsel=0; itsel<nsele; itsel++)
// 	    {
//             attribute_processing(son_of_move);
//             val_from_attr<string>("sel_mode",smode);
//             cerr << smode << '\t' << smode.compare("none") << '\t' << smode.compare("all") << endl;
            // if sele mode all or none we don't need the xtra selection attribute
//             if(smode.compare("none")!=0 && smode.compare("all")!=0)
//             {
//                 val_from_attr<string>("sele",selec);
//             }
//             seleList.push_back( tuple<string,string>(smode,selec) );
//             attrs_list.clear();
// 	    }
//             (*mvlist)->addNewMoveType(mvtyp, mvmode, smode, selec);
            //         (*mvlist)->addNewMoveType(mvtyp, mvmode, seleList);

//             }

//             (*ff)->setMcmvlst(**mvlist);

//             attrs_list.clear();
	    
        }
        // type of MC simulation, and the associated parameters go here
        else if(!son_name.compare("mc"))
        {
            attribute_processing(son_of_root);
	    
            val_from_attr<int>("nsteps",nsteps);
            if(nsteps!=0)
            {
                val_from_attr<double>("dmax_value",dmax_value);
                val_from_attr<double>("dmax_target",dmax_target);
                val_from_attr<int>("dmax_each",dmax_each);

                val_from_attr<int>("save_each",save_frequency);
            }

			val_from_attr<uint64_t>("seed", seed);

            *simulation = new MC_metropolis(atomList, **pbc, **ens, **ff, **mvlist, nsteps, save_frequency, dmax_value, dmax_target, dmax_each, seed);

            attrs_list.clear();
	    

        }
		else if (!son_name.compare("benchmark"))
		{
			this->benchmark = true;
		}
        else
        {
            cerr << "Warning : unknown section \"<" << son_name << ">\" will be ignored." << endl;
        }

    }//end iteration on all sons of root node inputFile

}

Parser_XML::~Parser_XML()
{
    if (is_mdbas)
        delete io;
    
    delete xmlFile;
    delete doc; 
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

bool Parser_XML::requiredBenchmark() const
{ 
	return benchmark; 
}

