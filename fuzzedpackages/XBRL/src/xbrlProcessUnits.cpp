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


RcppExport SEXP xbrlProcessUnits(SEXP epaDoc) {
  xmlDocPtr doc = (xmlDocPtr) R_ExternalPtrAddr(epaDoc);

  xmlXPathContextPtr context = xmlXPathNewContext(doc);
  xmlXPathObjectPtr unit_res = xmlXPathEvalExpression((xmlChar*) "//*[local-name()='unit']", context);
  xmlNodeSetPtr unit_nodeset = unit_res->nodesetval;
  int unit_nodeset_ln = unit_nodeset->nodeNr;
  xmlXPathFreeContext(context);

  CharacterVector unitId(unit_nodeset_ln);
  CharacterVector measure(unit_nodeset_ln);
  CharacterVector unitNumerator(unit_nodeset_ln);
  CharacterVector unitDenominator(unit_nodeset_ln);

  for (int i=0; i < unit_nodeset_ln; i++) {
    xmlNodePtr unit_node = unit_nodeset->nodeTab[i];
    xmlChar *tmp_str;
    if ((tmp_str = xmlGetProp(unit_node, (xmlChar*) "id"))) { 
      unitId[i] = (char *) tmp_str;
      xmlFree(tmp_str);
    } else {
      unitId[i] = NA_STRING;
    }
    measure[i] = unitNumerator[i] = unitDenominator[i] = NA_STRING;
    xmlNodePtr child_node = unit_node->xmlChildrenNode;
    while (child_node) {
      if (!xmlStrcmp(child_node->name, (xmlChar*) "measure")) {
	if ((tmp_str = xmlNodeListGetString(doc, child_node->xmlChildrenNode, 1))) {
	  measure[i] = (char *) tmp_str;
	  xmlFree(tmp_str);
	}
      } else if (!xmlStrcmp(child_node->name, (xmlChar*) "divide")) {
	xmlNodePtr gchild_node = child_node->xmlChildrenNode;
	while (gchild_node) {
	  if (!xmlStrcmp(gchild_node->name, (xmlChar*) "unitNumerator")) {
	    xmlNodePtr ggchild_node = gchild_node->xmlChildrenNode;
	    while (ggchild_node) {
	      if (!xmlStrcmp(ggchild_node->name, (xmlChar*) "measure")) {
		if ((tmp_str = xmlNodeListGetString(doc, ggchild_node->xmlChildrenNode, 1))) {
		  unitNumerator[i] = (char *) tmp_str;
		  xmlFree(tmp_str);
		}
	      }
	      ggchild_node = ggchild_node->next;
	    }
	  } else if (!xmlStrcmp(gchild_node->name, (xmlChar*) "unitDenominator")) {
	    xmlNodePtr ggchild_node = gchild_node->xmlChildrenNode;
	    while (ggchild_node) {
	      if (!xmlStrcmp(ggchild_node->name, (xmlChar*) "measure")) {
		if ((tmp_str = xmlNodeListGetString(doc, ggchild_node->xmlChildrenNode, 1))) {
		  unitDenominator[i] = (char *) tmp_str;
		  xmlFree(tmp_str);
		}
	      }
	      ggchild_node = ggchild_node->next;
	    }
	  }
	  gchild_node = gchild_node->next;
	}
      }
      child_node = child_node->next;
    }
  }
  xmlXPathFreeObject(unit_res);

  return DataFrame::create(Named("unitId")=unitId,
			   Named("measure")=measure,
			   Named("unitNumerator")=unitNumerator,
			   Named("unitDenominator")=unitDenominator);
}
