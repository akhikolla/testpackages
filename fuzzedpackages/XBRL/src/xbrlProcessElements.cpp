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


RcppExport SEXP xbrlProcessElements(SEXP epaDoc) {
  xmlDocPtr doc = (xmlDocPtr) R_ExternalPtrAddr(epaDoc);

  xmlXPathContextPtr context = xmlXPathNewContext(doc);
  xmlXPathObjectPtr schema_res = xmlXPathEvalExpression((xmlChar*) "//*[local-name()='schema']", context);
  xmlNodeSetPtr schema_nodeset = schema_res->nodesetval;
  xmlChar *ns_txt;
  ns_txt = xmlGetProp(schema_nodeset->nodeTab[0], (xmlChar*) "targetNamespace");

  xmlXPathObjectPtr element_res = xmlXPathEvalExpression((xmlChar*) "//*[local-name()='element'][@*[local-name()='periodType']]", context);
  xmlNodeSetPtr element_nodeset = element_res->nodesetval;
  xmlXPathFreeContext(context);

  int element_nodeset_ln = element_nodeset->nodeNr;

  CharacterVector elementId(element_nodeset_ln);
  CharacterVector type(element_nodeset_ln);
  CharacterVector substitutionGroup(element_nodeset_ln);
  CharacterVector periodType(element_nodeset_ln);
  CharacterVector abstract(element_nodeset_ln);
  CharacterVector nillable(element_nodeset_ln);
  CharacterVector balance(element_nodeset_ln);
  CharacterVector ns(element_nodeset_ln);

  for (int i=0; i < element_nodeset_ln; i++) {
    xmlNodePtr element_node = element_nodeset->nodeTab[i];

    xmlChar *tmp_str;
    if ((tmp_str = xmlGetProp(element_node, (xmlChar*) "id"))) { 
      elementId[i] = (char *) tmp_str;
      xmlFree(tmp_str);
    } else {
      elementId[i] = NA_STRING;
    }
    if ((tmp_str = xmlGetProp(element_node, (xmlChar*) "type"))) { 
      type[i] = (char *) tmp_str;
      xmlFree(tmp_str);
    } else {
      type[i] = NA_STRING;
    }
    if ((tmp_str = xmlGetProp(element_node, (xmlChar*) "substitutionGroup"))) { 
      substitutionGroup[i] = (char *) tmp_str;
      xmlFree(tmp_str);
    } else {
      substitutionGroup[i] = NA_STRING;
    }
    if ((tmp_str = xmlGetProp(element_node, (xmlChar*) "periodType"))) { 
      periodType[i] = (char *) tmp_str;
      xmlFree(tmp_str);
    } else {
      periodType[i] = NA_STRING;
    }
    if ((tmp_str = xmlGetProp(element_node, (xmlChar*) "abstract"))) { 
      abstract[i] = (char *) tmp_str;
      xmlFree(tmp_str);
    } else {
      abstract[i] = NA_STRING;
    }
    if ((tmp_str = xmlGetProp(element_node, (xmlChar*) "nillable"))) { 
      nillable[i] = (char *) tmp_str;
      xmlFree(tmp_str);
    } else {
      nillable[i] = NA_STRING;
    }
    if ((tmp_str = xmlGetProp(element_node, (xmlChar*) "balance"))) { 
      balance[i] = (char *) tmp_str;
      xmlFree(tmp_str);
    } else {
      balance[i] = NA_STRING;
    }
    ns[i] = (char *) ns_txt;
  }
  xmlFree(ns_txt);
  xmlXPathFreeObject(element_res);
  xmlXPathFreeObject(schema_res);

  return DataFrame::create(Named("elementId")=elementId,
			   Named("type")=type,
			   Named("substitutionGroup")=substitutionGroup,
			   Named("periodType")=periodType,
			   Named("abstract")=abstract,
			   Named("nillable")=nillable,
			   Named("balance")=balance,
			   Named("ns")=ns);
}
