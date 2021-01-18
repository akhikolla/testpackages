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


RcppExport SEXP xbrlGetImportNames(SEXP epaDoc) {
  xmlDocPtr doc = (xmlDocPtr) R_ExternalPtrAddr(epaDoc);

  xmlXPathContextPtr context = xmlXPathNewContext(doc);
  xmlXPathObjectPtr import_res = xmlXPathEvalExpression((xmlChar*) "//*[local-name()='import']", context);
  xmlNodeSetPtr import_nodeset = import_res->nodesetval;
  xmlXPathFreeContext(context);

  int import_nodeset_ln = import_nodeset->nodeNr;

  CharacterVector schemaLocation(import_nodeset_ln);

  for (int i=0; i < import_nodeset_ln; i++) {
    xmlNodePtr import_node = import_nodeset->nodeTab[i];

    xmlChar *tmp_str;
    if ((tmp_str = xmlGetProp(import_node, (xmlChar*) "schemaLocation"))) { 
      schemaLocation[i] = (char *) tmp_str;
      xmlFree(tmp_str);
    } else {
      schemaLocation[i] = NA_STRING;
    }
  }
  xmlXPathFreeObject(import_res);

  return schemaLocation;
}
