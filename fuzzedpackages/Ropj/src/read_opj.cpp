/*
Copyright (C) 2019 Ivan Krylov

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#include <algorithm> // max, all_of
#include <iterator> // distance
#include <string> // string

#include <Rcpp.h>
using namespace Rcpp;

#include "OriginFile.h"

#include "R_ext/Riconv.h"
class decoder {
	void * cd;
public:
	decoder(const char * from) {
		cd = Riconv_open("", from);
		if (cd == (void*)(-1))
			throw std::invalid_argument(std::string("Cannot decode from ") + from);
	}
	~decoder() {
		Riconv_close(cd);
	}
	String operator()(const std::string & s) {
		std::string out(s.size(), 0);

		const char * inbuf = s.c_str();
		char * outbuf = &out[0]; // ick
		size_t inbytesleft = s.size(), outbytesleft = out.size();

		// this is what happens when you bring C to an STL fight
		while (Riconv(cd, &inbuf, &inbytesleft, &outbuf, &outbytesleft) == (size_t)-1) {
			if (errno != E2BIG) throw std::runtime_error("Cannot decode string");
			ptrdiff_t pos = outbuf - &out[0];
			outbytesleft += out.size();
			out.resize(out.size() * 2);
			outbuf = &out[pos];
		}
		// get rid of trailing \0, if any resulted from over-allocation
		out.resize(out.size() - outbytesleft);
		// there may have been NULs in the input string, too
		auto nulpos = out.find('\0');
		if (nulpos != std::string::npos) out.erase(nulpos);

		return String(out, CE_NATIVE);
	}
};

static List import_spreadsheet(const Origin::SpreadSheet & osp, decoder & dec) {
	List rsp(osp.columns.size());
	StringVector names(rsp.size()), comments(rsp.size()), commands(rsp.size());

	size_t maxRows = osp.maxRows;
	for (const Origin::SpreadColumn & osc : osp.columns)
		maxRows = std::max(osc.data.size(), maxRows);

	for (unsigned int c = 0; c < osp.columns.size(); c++) {
		const Origin::SpreadColumn & ocol = osp.columns[c];
		names[c] = dec(ocol.name);
		comments[c] = dec(ocol.comment);
		commands[c] = dec(ocol.command);
		if (
			std::all_of(
				ocol.data.begin(), ocol.data.end(),
				[](const Origin::variant & v){
					return v.type() == Origin::variant::V_DOUBLE;
				}
			)
		){
			NumericVector ncol(maxRows, NA_REAL);
			for (size_t row = 0; row < ocol.data.size(); row++) {
				ncol[row] = ocol.data[row].as_double();
				if (ncol[row] == _ONAN) ncol[row] = R_NaN;
			}
			rsp[c] = ncol;
		} else {
			StringVector ccol(maxRows, NA_STRING);
			for (size_t row = 0; row < ocol.data.size(); row++) {
				const Origin::variant & v = ocol.data[row];
				if (v.type() == Origin::variant::V_DOUBLE) {
					if (v.as_double() != _ONAN) ccol[row] = std::to_string(v.as_double()); // yuck
				} else {
					ccol[row] = dec(v.as_string());
				}
			}
			rsp[c] = ccol;
		}
	}

	rsp.attr("names") = std::move(names);
	rsp.attr("comments") = comments; // XXX: Ropj <= 0.2-2
	rsp.attr("comment") = std::move(comments);
	rsp.attr("commands") = std::move(commands);
	rsp.attr("type") = "spreadsheet";
	return rsp;
}

static List import_matrix(const Origin::Matrix & omt, decoder & dec) {
	List ret(omt.sheets.size());
	StringVector names(ret.size()), commands(ret.size());
	for (unsigned int i = 0; i < omt.sheets.size(); i++) {
		NumericMatrix rms(Dimension(omt.sheets[i].columnCount, omt.sheets[i].rowCount));
		std::copy_n(omt.sheets[i].data.begin(), rms.size(), rms.begin());
		std::replace(rms.begin(), rms.end(), _ONAN, R_NaN);

		rms.attr("dimensions") = NumericVector(omt.sheets[i].coordinates.begin(), omt.sheets[i].coordinates.end());
		ret[i] = std::move(rms);
		names[i] = dec(omt.sheets[i].name);
		commands[i] = dec(omt.sheets[i].command);
	}
	ret.attr("names") = names;
	ret.attr("commands") = commands;
	ret.attr("type") = "matrix";
	return ret;
}

static List import_tree(
	tree<Origin::ProjectNode>::sibling_iterator cur,
	tree<Origin::ProjectNode>::sibling_iterator end,
	decoder & dec
) {
	unsigned int i = 0;
	List ret(std::distance(cur, end));
	StringVector names(ret.size());

	for (; cur != end; ++cur, ++i) {
		String nm = dec(cur->name);
		names[i] = nm;
		switch(cur->type) {
		case Origin::ProjectNode::Folder:
			ret[i] = import_tree(cur.begin(), cur.end(), dec);
			break;
		case Origin::ProjectNode::SpreadSheet:
		case Origin::ProjectNode::Matrix:
		case Origin::ProjectNode::Excel:
		case Origin::ProjectNode::Note:
			ret[i] = nm;
			break;
		default: // ignore types we don't understand
			break;
		}
	}

	ret.attr("names") = names;
	return ret;
}

// [[Rcpp::export(name="read_opj")]]
List read_opj(const std::string & file, const char * encoding, bool tree) {
	decoder dec(encoding);
	OriginFile opj(file);

	if (!opj.parse()) stop("Failed to open and/or parse " + file); // throws

	unsigned int j = 0,
		items = opj.spreadCount() + opj.excelCount() + opj.matrixCount() + opj.noteCount();
	List ret(items);
	StringVector retn(items), retl(items);


	for (unsigned int i = 0; i < opj.spreadCount(); i++, j++) {
		const Origin::SpreadSheet & osp = opj.spread(i);
		retn[j] = dec(osp.name);
		retl[j] = dec(osp.label);
		ret[j] = import_spreadsheet(osp, dec);
	}

	for (unsigned int i = 0; i < opj.excelCount(); i++, j++) {
		const Origin::Excel & oex = opj.excel(i);
		retn[j] = dec(oex.name);
		retl[j] = dec(oex.label);

		List exl(oex.sheets.size());
		StringVector exln(oex.sheets.size());

		for (size_t sp = 0; sp < oex.sheets.size(); sp++) {
			exl[sp] = import_spreadsheet(oex.sheets[sp], dec);
			exln[sp] = dec(oex.sheets[sp].name);
		}

		exl.attr("names") = exln;
		exl.attr("type") = "excel";
		ret[j] = exl;
	}

	for (unsigned int i = 0; i < opj.matrixCount(); i++, j++) {
		const Origin::Matrix & omt = opj.matrix(i);
		retn[j] = dec(omt.name);
		retl[j] = dec(omt.label);
		ret[j] = import_matrix(omt, dec);
	}

	for (unsigned int i = 0; i < opj.noteCount(); i++, j++) {
		const Origin::Note & ont = opj.note(i);
		retn[j] = dec(ont.name);
		retl[j] = dec(ont.label);
		ret[j] = dec(ont.text);
	}

	ret.attr("names") = retn;
	ret.attr("comment") = retl;
	if (tree)
		ret.attr("tree") = import_tree(
			// must skip the root of the tree when iterating over it
			opj.project()->begin().begin(), opj.project()->begin().end(),
			dec
		);
	return ret;
}
