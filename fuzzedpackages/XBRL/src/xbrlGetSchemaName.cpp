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


RcppExport SEXP xbrlGetSchemaName(SEXP epaDoc) {
  xmlDocPtr doc = (xmlDocPtr) R_ExternalPtrAddr(epaDoc);

  xmlXPathContextPtr context = xmlXPathNewContext(doc);
  xmlXPathObjectPtr schemaRef_res = xmlXPathEvalExpression((xmlChar*) "//*[local-name()='schemaRef']", context);
  xmlNodeSetPtr schemaRef_nodeset = schemaRef_res->nodesetval;
  xmlXPathFreeContext(context);

  int schemaRef_nodeset_ln = schemaRef_nodeset->nodeNr;

  CharacterVector href(schemaRef_nodeset_ln);

  for (int i=0; i < schemaRef_nodeset_ln; i++) {
    xmlNodePtr schemaRef_node = schemaRef_nodeset->nodeTab[i];

    xmlChar *tmp_str;
    if ((tmp_str = xmlGetProp(schemaRef_node, (xmlChar*) "href"))) { 
      href[i] = (char *) tmp_str;
      xmlFree(tmp_str);
    } else {
      href[i] = NA_STRING;
    }
  }
  xmlXPathFreeObject(schemaRef_res);

  return href;
}
