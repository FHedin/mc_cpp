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

#include "Parser.h"

using namespace std;
using namespace rapidxml;

Parser_XML::Parser_XML(const char inpfileName[], bool _verbose) : verbose(_verbose)
{
    xmlFile = new file<>(inpfileName);
    doc = new xml_document<>();

    doc->parse<0>(xmlFile->data());
    xml_node<> *root = doc->first_node();

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



