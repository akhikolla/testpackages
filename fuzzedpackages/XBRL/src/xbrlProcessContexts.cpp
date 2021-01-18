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


RcppExport SEXP xbrlProcessContexts(SEXP epaDoc) {
  xmlDocPtr doc = (xmlDocPtr) R_ExternalPtrAddr(epaDoc);

  xmlXPathContextPtr context = xmlXPathNewContext(doc);
  xmlXPathObjectPtr context_res = xmlXPathEvalExpression((xmlChar*) "//*[local-name()='context']", context);
  xmlNodeSetPtr context_nodeset = context_res->nodesetval;
  int context_nodeset_ln = context_nodeset->nodeNr;
  xmlXPathFreeContext(context);

  CharacterVector contextId(context_nodeset_ln);
  CharacterVector scheme(context_nodeset_ln);
  CharacterVector identifier(context_nodeset_ln);
  CharacterVector startDate(context_nodeset_ln);
  CharacterVector endDate(context_nodeset_ln);
  CharacterVector dimension1(context_nodeset_ln);
  CharacterVector value1(context_nodeset_ln);
  CharacterVector dimension2(context_nodeset_ln);
  CharacterVector value2(context_nodeset_ln);
  CharacterVector dimension3(context_nodeset_ln);
  CharacterVector value3(context_nodeset_ln);
  CharacterVector dimension4(context_nodeset_ln);
  CharacterVector value4(context_nodeset_ln);

  for (int i=0; i < context_nodeset_ln; i++) {
    xmlNodePtr context_node = context_nodeset->nodeTab[i];
    xmlChar *tmp_str;
    if ((tmp_str = xmlGetProp(context_node, (xmlChar*) "id"))) { 
      contextId[i] = (char *) tmp_str;
      xmlFree(tmp_str);
    } else {
      contextId[i] = NA_STRING;
    }
    scheme[i] = identifier[i] = startDate[i] = endDate[i] =
      dimension1[i] = value1[i] = dimension2[i] = value2[i] =
      dimension3[i] = value3[i] = dimension4[i] = value4[i] = NA_STRING;
    xmlNodePtr child_node = context_node->xmlChildrenNode;
    while (child_node) {
      if (!xmlStrcmp(child_node->name, (xmlChar*) "entity")) {
	xmlNodePtr gchild_node = child_node->xmlChildrenNode;
	while (gchild_node) {
	  if (!xmlStrcmp(gchild_node->name, (xmlChar*) "identifier")) {
	    if ((tmp_str = xmlGetProp(gchild_node, (xmlChar*) "scheme"))) { 
	      scheme[i] = (char *) tmp_str;
	      xmlFree(tmp_str);
	    }
	    if ((tmp_str = xmlNodeListGetString(doc, gchild_node->xmlChildrenNode, 1))) {
	      identifier[i] = (char *) tmp_str;
	      xmlFree(tmp_str);
	    }
	  } else if (!xmlStrcmp(gchild_node->name, (xmlChar*) "segment")) {
	    xmlNodePtr ggchild_node = gchild_node->xmlChildrenNode;
	    int dimn = 1;
	    while (ggchild_node) {
	      if (!xmlStrcmp(ggchild_node->name, (xmlChar*) "explicitMember")) {
		if ((tmp_str = xmlGetProp(ggchild_node, (xmlChar*) "dimension"))) {
		  if (dimn == 1)
		    dimension1[i] = (char *) tmp_str;
		  else if (dimn == 2)
		    dimension2[i] = (char *) tmp_str;
		  else if (dimn == 3)
		    dimension3[i] = (char *) tmp_str;
		  else if (dimn == 4)
		    dimension4[i] = (char *) tmp_str;
		  xmlFree(tmp_str);
		}
		if ((tmp_str = xmlNodeListGetString(doc, ggchild_node->xmlChildrenNode, 1))) {
		  if (dimn == 1)
		    value1[i] = (char *) tmp_str;
		  else if (dimn == 2)
		    value2[i] = (char *) tmp_str;
		  else if (dimn == 3)
		    value3[i] = (char *) tmp_str;
		  else if (dimn == 4)
		    value4[i] = (char *) tmp_str;
		  xmlFree(tmp_str);
		}
		dimn++;
	      }
	      ggchild_node = ggchild_node->next;
	    }
	  }
	  gchild_node = gchild_node->next;
	}
      } else if (!xmlStrcmp(child_node->name, (xmlChar*) "period")) {
	xmlNodePtr gchild_node = child_node->xmlChildrenNode;
	while (gchild_node) {
	  if (!xmlStrcmp(gchild_node->name, (xmlChar*) "startDate")) {
	    if ((tmp_str = xmlNodeListGetString(doc, gchild_node->xmlChildrenNode, 1))) {
	      startDate[i] = (char *) tmp_str;
	      xmlFree(tmp_str);
	    }
	  } else if (!xmlStrcmp(gchild_node->name, (xmlChar*) "endDate")) {
	    if ((tmp_str = xmlNodeListGetString(doc, gchild_node->xmlChildrenNode, 1))) {
	      endDate[i] = (char *) tmp_str;
	      xmlFree(tmp_str);
	    }
	  } else if (!xmlStrcmp(gchild_node->name, (xmlChar*) "instant")) {
	    if ((tmp_str = xmlNodeListGetString(doc, gchild_node->xmlChildrenNode, 1))) {
	      endDate[i] = (char *) tmp_str;
	      xmlFree(tmp_str);
	    }
	  }
	  gchild_node = gchild_node->next;
	}
      }
      child_node = child_node->next;
    }
  }
  xmlXPathFreeObject(context_res);

  return DataFrame::create(Named("contextId")=contextId,
			   Named("scheme")=scheme,
			   Named("identifier")=identifier,
			   Named("startDate")=startDate,
			   Named("endDate")=endDate,
			   Named("dimension1")=dimension1,
			   Named("value1")=value1,
			   Named("dimension2")=dimension2,
			   Named("value2")=value2,
			   Named("dimension3")=dimension3,
			   Named("value3")=value3,
			   Named("dimension4")=dimension4,
			   Named("value4")=value4);
}
