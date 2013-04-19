#include <iostream>

#include "rapidxml.hpp"
#include "rapidxml_utils.hpp"

using namespace std;
using namespace rapidxml;

void node_processing(xml_node<> *src);
void attribute_processing(xml_node<> *src);
void check_has_son(xml_node<> *src);

int main()
{
  file<> xmlFile("test.xml");
  xml_document<> doc;
  doc.parse<0>(xmlFile.data());
  
  node_processing(&doc);

  return 0;

}

void node_processing(xml_node<> *src)
{
  for(xml_node<> *node = src->first_node();node; node = node->next_sibling())
  {
    cout << "Name of node is: '" << node->name() << "'" << endl;
    attribute_processing(node);
    
    check_has_son(node);
  }
  
}

void attribute_processing(xml_node<> *src)
{
  for(xml_attribute<> *attr = src->first_attribute();attr; attr = attr->next_attribute())
  {
    cout << "Node has attribute '" << attr->name() << "' ";
    cout << "with value '" << attr->value() << "'" << endl;
  }
  
  cout << endl;
  
}

void check_has_son(xml_node<> *src)
{
  xml_node<> *son = NULL;
  son = src->first_node();
  if(son != NULL)
  {
    cout << "Node has a son with name '" << son->name() << "'" << endl;
    attribute_processing(son);
  }
}

