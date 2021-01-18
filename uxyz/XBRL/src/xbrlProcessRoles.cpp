// Copyright (C) 2014-2016 Roberto Bertolusso
//
// This file is part of XBRL.
//
// XBRL is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// XBRL is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with XBRL. If not, see <http://www.gnu.org/licenses/>.

#include "XBRL.h"


RcppExport SEXP xbrlProcessRoles(SEXP epaDoc) {
  xmlDocPtr doc = (xmlDocPtr) R_ExternalPtrAddr(epaDoc);

  xmlXPathContextPtr context = xmlXPathNewContext(doc);
  xmlXPathObjectPtr roleType_res = xmlXPathEvalExpression((xmlChar*) "//*[local-name()='roleType']", context);
  xmlNodeSetPtr roleType_nodeset = roleType_res->nodesetval;
  xmlXPathFreeContext(context);

  int roleType_nodeset_ln = roleType_nodeset->nodeNr;

  CharacterVector roleId(roleType_nodeset_ln);
  CharacterVector order(roleType_nodeset_ln);
  CharacterVector type(roleType_nodeset_ln);
  CharacterVector description(roleType_nodeset_ln);
  CharacterVector definition(roleType_nodeset_ln);

  for (int i=0; i < roleType_nodeset_ln; i++) {
    xmlNodePtr roleType_node = roleType_nodeset->nodeTab[i];

    xmlChar *tmp_str;
    if ((tmp_str = xmlGetProp(roleType_node, (xmlChar*) "roleURI"))) { 
      roleId[i] = (char *) tmp_str;
      xmlFree(tmp_str);
    } else {
      roleId[i] = NA_STRING;
    }
    xmlNodePtr child_node = roleType_node->xmlChildrenNode;
    while (child_node) {
      if (xmlStrcmp(child_node->name, (xmlChar*) "definition")) {
	child_node = child_node->next;
	continue;
      }
      if ((tmp_str = xmlNodeListGetString(doc, child_node->xmlChildrenNode, 1))) {
	definition[i] = (char *) tmp_str;
	string str = (char *) tmp_str;
	xmlFree(tmp_str);
	size_t found1 = str.find(" - ");
	if (found1 != string::npos) {
	  order[i] = str.substr(0, found1);
	  size_t found2 = str.find(" - ", found1+3);
	  if (found2 != string::npos) {
	    type[i] = str.substr(found1+3, found2-found1-3);
	    description[i] = str.substr(found2+3);
	  }
	}
      }

      break;
    }
  }
  xmlXPathFreeObject(roleType_res);

  return DataFrame::create(Named("roleId")=roleId,
			   Named("order")=order,
			   Named("type")=type,
			   Named("description")=description,
			   Named("definition")=definition);
}
