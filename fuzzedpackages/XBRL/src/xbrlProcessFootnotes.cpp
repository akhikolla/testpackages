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


RcppExport SEXP xbrlProcessFootnotes(SEXP epaDoc) {
  xmlDocPtr doc = (xmlDocPtr) R_ExternalPtrAddr(epaDoc);

  xmlXPathContextPtr context = xmlXPathNewContext(doc);
  xmlXPathObjectPtr footnoteArc_res = xmlXPathEvalExpression((xmlChar*) "//*[local-name()='footnoteArc']", context);
  xmlNodeSetPtr footnoteArc_nodeset = footnoteArc_res->nodesetval;
  xmlXPathObjectPtr loc_res = xmlXPathEvalExpression((xmlChar*) "//*[local-name()='loc']", context);
  xmlNodeSetPtr loc_nodeset = loc_res->nodesetval;
  xmlXPathObjectPtr footnote_res = xmlXPathEvalExpression((xmlChar*) "//*[local-name()='footnote']", context);
  xmlNodeSetPtr footnote_nodeset = footnote_res->nodesetval;
  xmlXPathFreeContext(context);

  int footnoteArc_nodeset_ln = footnoteArc_nodeset->nodeNr;

  CharacterVector factId(footnoteArc_nodeset_ln);
  CharacterVector lang(footnoteArc_nodeset_ln);
  CharacterVector footnoteRole(footnoteArc_nodeset_ln);
  CharacterVector footnoteString(footnoteArc_nodeset_ln);
  CharacterVector href(footnoteArc_nodeset_ln);

  for (int i=0; i < footnoteArc_nodeset_ln; i++) {
    xmlNodePtr footnoteArc_node = footnoteArc_nodeset->nodeTab[i];
    xmlChar *footnoteArc_to = xmlGetProp(footnoteArc_node, (xmlChar*) "to");
    xmlChar *footnoteArc_from = xmlGetProp(footnoteArc_node, (xmlChar*) "from");
    for (int j=0; j < loc_nodeset->nodeNr; j++) {
      xmlNodePtr loc_node = loc_nodeset->nodeTab[j];
      xmlChar *loc_label = xmlGetProp(loc_node, (xmlChar*) "label");

      int nomatch = xmlStrcmp(loc_label, footnoteArc_from);
      xmlFree(loc_label);
      if (nomatch)
	continue;
      xmlFree(footnoteArc_from);  // There is a match. Not needed anymore.
      // Footnotes can be reused, and they are not necessarily ordered.
      // So we need to check them all.
      for (int k=0; k < footnote_nodeset->nodeNr; k++) {
	xmlNodePtr footnote_node = footnote_nodeset->nodeTab[k];
	xmlChar *footnote_label = xmlGetProp(footnote_node, (xmlChar*) "label");

	int nomatch = xmlStrcmp(footnoteArc_to, footnote_label);
	xmlFree(footnote_label);
	if (nomatch)
	  continue;
	xmlFree(footnoteArc_to);  // There is a match. Not needed anymore.

	xmlChar *tmp_str;
	if ((tmp_str = xmlGetProp(loc_node, (xmlChar*) "href"))) {
	  href[i] = (char *) tmp_str;
	  string str = (char *) tmp_str;
	  xmlFree(tmp_str);
	  size_t found = str.find("#");
	  if (found != string::npos) {
	    str.replace(0, found+1, "");
	    factId[i] = str;
	  }
	} else {
	  href[i] = NA_STRING;
	  factId[i] = NA_STRING;
	}
	if ((tmp_str = xmlGetProp(footnote_node, (xmlChar*) "lang"))) {
	  lang[i] = (char *) tmp_str;
	  xmlFree(tmp_str);
	} else {
	  lang[i] = NA_STRING;
	}
	if ((tmp_str = xmlGetProp(footnote_node, (xmlChar*) "role"))) { 
	  footnoteRole[i] = (char *) tmp_str;
	  xmlFree(tmp_str);
	} else {
	  footnoteRole[i] = NA_STRING;
	}
	if ((tmp_str = xmlNodeListGetString(doc, footnote_node->xmlChildrenNode, 1))) {
	  footnoteString[i] = (char *) tmp_str;
	  xmlFree(tmp_str);
	} else {
	  footnoteString[i] = NA_STRING;
	}
	break;
      }
      break;
    }
  }
  xmlXPathFreeObject(footnote_res);
  xmlXPathFreeObject(loc_res);
  xmlXPathFreeObject(footnoteArc_res);

  if (footnoteArc_nodeset_ln == 0)
    return R_NilValue;

  return DataFrame::create(Named("factId")=factId,
			   Named("lang")=lang,
			   Named("footnoteRole")=footnoteRole,
			   Named("footnoteString")=footnoteString,
			   Named("href")=href);
}
