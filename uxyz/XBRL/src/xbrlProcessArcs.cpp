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


RcppExport SEXP xbrlProcessArcs(SEXP epaDoc, SEXP arcTypeS) {
  xmlDocPtr doc = (xmlDocPtr) R_ExternalPtrAddr(epaDoc);
  CharacterVector arcType(arcTypeS);

  string arcType_str = (string) arcType[0];
  xmlXPathContextPtr context = xmlXPathNewContext(doc);
  string typeArc_str = "//*[local-name()='" + arcType_str + "Arc']";
  xmlXPathObjectPtr typeArc_res = xmlXPathEvalExpression((xmlChar*) typeArc_str.data(), context);
  int typeArc_nodeset_ln = typeArc_res->nodesetval->nodeNr;
  xmlXPathFreeObject(typeArc_res);

  if (typeArc_nodeset_ln == 0) {
    xmlXPathFreeContext(context);
    return R_NilValue;
  }

  CharacterVector roleId(typeArc_nodeset_ln);
  CharacterVector fromElementId(typeArc_nodeset_ln);
  CharacterVector toElementId(typeArc_nodeset_ln);
  CharacterVector arcrole(typeArc_nodeset_ln);
  CharacterVector order(typeArc_nodeset_ln);
  CharacterVector closed(typeArc_nodeset_ln);
  CharacterVector usable(typeArc_nodeset_ln);
  CharacterVector contextElement(typeArc_nodeset_ln);
  CharacterVector preferredLabel(typeArc_nodeset_ln);
  CharacterVector fromHref(typeArc_nodeset_ln);
  CharacterVector toHref(typeArc_nodeset_ln);

  xmlXPathObjectPtr typeLink_res = xmlXPathEvalExpression((xmlChar*) ("//*[local-name()='" + arcType_str + "Link']").data(), context);
  xmlNodeSetPtr typeLink_nodeset = typeLink_res->nodesetval;

  xmlChar *tmp_str;
  int r=0;
  for (int i=0; i < typeLink_nodeset->nodeNr; i++) {
    xmlNodePtr typeLink_node = typeLink_nodeset->nodeTab[i];
    context->node = typeLink_node;  // Search inside this typeLink_node.
    xmlXPathObjectPtr typeArc_res = xmlXPathEvalExpression((xmlChar*) ("*[local-name()='" + arcType_str + "Arc']").data(), context);
    xmlNodeSetPtr typeArc_nodeset = typeArc_res->nodesetval;
    xmlXPathObjectPtr loc_res = xmlXPathEvalExpression((xmlChar*) "*[local-name()='loc']", context);
    xmlNodeSetPtr loc_nodeset = loc_res->nodesetval;

    for (int j=0; j < typeArc_nodeset->nodeNr; j++) {
      xmlNodePtr typeArc_node = typeArc_nodeset->nodeTab[j];
      xmlChar *typeArc_from = xmlGetProp(typeArc_node, (xmlChar*) "from");
      xmlChar *typeArc_to = xmlGetProp(typeArc_node, (xmlChar*) "to");

      int matches = 0;
      for (int k=0; k < loc_nodeset->nodeNr; k++) {
	xmlNodePtr loc_node = loc_nodeset->nodeTab[k];
	xmlChar *loc_label = xmlGetProp(loc_node, (xmlChar*) "label");

	if (!xmlStrcmp(loc_label, typeArc_from)) {
	  if ((tmp_str = xmlGetProp(loc_node, (xmlChar*) "href"))) { 
	    fromHref[r] = (char *) tmp_str;
	    string str = (char *) tmp_str;
	    xmlFree(tmp_str);
	    size_t found = str.find("#");
	    if (found != string::npos) {
	      str.replace(0, found+1, "");
	      fromElementId[r] = str;
	    }
	    matches++;
	  }
	} else if (!xmlStrcmp(loc_label, typeArc_to)) {
	  if ((tmp_str = xmlGetProp(loc_node, (xmlChar*) "href"))) { 
	    toHref[r] = (char *) tmp_str;
	    string str = (char *) tmp_str;
	    xmlFree(tmp_str);
	    size_t found = str.find("#");
	    if (found != string::npos) {
	      str.replace(0, found+1, "");
	      toElementId[r] = str;
	    }
	    matches++;
	  }
	}
	xmlFree(loc_label);
	if (matches == 2)
	  break;
      }
      xmlFree(typeArc_from);
      xmlFree(typeArc_to);

      if ((tmp_str = xmlGetProp(typeLink_node, (xmlChar*) "role"))) { 
	roleId[r] = (char *) tmp_str;
	xmlFree(tmp_str);
      } else {
	roleId[r] = NA_STRING;
      }
      if ((tmp_str = xmlGetProp(typeArc_node, (xmlChar*) "arcrole"))) { 
	arcrole[r] = (char *) tmp_str;
	xmlFree(tmp_str);
      } else {
	arcrole[r] = NA_STRING;
      }
      if ((tmp_str = xmlGetProp(typeArc_node, (xmlChar*) "order"))) { 
	order[r] = (char *) tmp_str;
	xmlFree(tmp_str);
      } else {
	order[r] = NA_STRING;
      }
      if ((tmp_str = xmlGetProp(typeArc_node, (xmlChar*) "closed"))) { 
	closed[r] = (char *) tmp_str;
	xmlFree(tmp_str);
      } else {
	closed[r] = NA_STRING;
      }
      if ((tmp_str = xmlGetProp(typeArc_node, (xmlChar*) "usable"))) { 
	usable[r] = (char *) tmp_str;
	xmlFree(tmp_str);
      } else {
	usable[r] = NA_STRING;
      }
      if ((tmp_str = xmlGetProp(typeArc_node, (xmlChar*) "contextElement"))) { 
	contextElement[r] = (char *) tmp_str;
	xmlFree(tmp_str);
      } else {
	contextElement[r] = NA_STRING;
      }
      if ((tmp_str = xmlGetProp(typeArc_node, (xmlChar*) "preferredLabel"))) { 
	preferredLabel[r] = (char *) tmp_str;
	xmlFree(tmp_str);
      } else {
	preferredLabel[r] = NA_STRING;
      }

      r++;
    }
    xmlXPathFreeObject(loc_res);
    xmlXPathFreeObject(typeArc_res);
  }
  xmlXPathFreeObject(typeLink_res);
  xmlXPathFreeContext(context);

  return DataFrame::create(Named("roleId")=roleId,
			   Named("fromElementId")=fromElementId,
			   Named("toElementId")=toElementId,
			   Named("arcrole")=arcrole,
			   Named("order")=order,
			   Named("closed")=closed,
			   Named("usable")=usable,
			   Named("contextElement")=contextElement,
			   Named("preferredLabel")=preferredLabel,
			   Named("fromHref")=fromHref,
			   Named("toHref")=toHref);
}
