#include <iostream>
#include <map>
#include <string>

#include "Parser.h"

using namespace rapidxml;

Parser_XML::Parser_XML(const char inpfileName[], bool verbose)
{
    file<> xmlFile(inpfileName);
    xml_document<> doc;
    doc.parse<0>(xmlFile.data());
    
    node_processing(&doc);
    
    if(verbose)
        Dump();
}

Parser_XML::~Parser_XML()
{
}

void Parser_XML::Dump()
{
    std::cout << "Dump of nodes_list [node_name => number_of_attributes] : " << std::endl;
    std::map<std::string,int>::iterator it_nodes;
    for( it_nodes=nodes_list.begin(); it_nodes!=nodes_list.end(); ++it_nodes )
    {
        std::cout << it_nodes->first << " => " << it_nodes->second << std::endl;
    }
    
    std::cout << std::endl << "Dump of attrs_list [attr_name => attr_value] : " << std::endl;
    std::map<std::string,std::string>::iterator it_attrs;
    for( it_attrs=attrs_list.begin(); it_attrs!=attrs_list.end(); ++it_attrs )
    {
        std::cout << it_attrs->first << " => " << it_attrs->second << std::endl;
    }
}


void Parser_XML::node_processing(xml_node<> *src)
{
    int n_attr;
    for(xml_node<> *node = src->first_node();node; node = node->next_sibling())
    {
        n_attr = attribute_processing(node);
        nodes_list.insert( std::pair<std::string,int>(node->name(),n_attr) );
        
        check_has_son(node);
    }
}

int Parser_XML::attribute_processing(xml_node<> *src)
{
    int n_attr=0;
    for(xml_attribute<> *attr = src->first_attribute();attr; attr = attr->next_attribute())
    {
        n_attr++;
        attrs_list.insert( std::pair<std::string,std::string>(attr->name(),attr->value()) );
    }
    return n_attr;
}

void Parser_XML::check_has_son(xml_node<> *src)
{
    xml_node<> *son = NULL;
    son = src->first_node();
    int n_attr;
    if(son != NULL)
    {
        n_attr = attribute_processing(son);
        nodes_list.insert( std::pair<std::string,int>(son->name(),n_attr) );
    }
}



