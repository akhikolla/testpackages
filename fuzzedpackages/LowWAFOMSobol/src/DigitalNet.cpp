/**
 * @file DigitalNet.cpp
 *
 * @brief DigitalNet class for Quasi Monte-Carlo Method.
 *
 * @note Currently only 64-bit DigitalNet is tested.
 *
 * @author Shinsuke Mori (Hiroshima University)
 * @author Makoto Matsumoto (Hiroshima University)
 * @author Mutsuo Saito
 *
 * Copyright (C) 2017 Shinsuke Mori, Makoto Matsumoto, Mutsuo Saito
 * and Hiroshima University.
 * All rights reserved.
 *
 * The GPL ver.3 is applied to this software, see
 * COPYING
 */
#include "config.h"
#include "bit_operator.h"
#include "DigitalNet.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <stdlib.h>
#include <cstring>
#include <cstdio>
#include <cerrno>

#if !defined(IN_RCPP)
#include "digital.h"
#include <sqlite3.h>
#else
using namespace Rcpp;
#endif

#if defined(BINARY_DF)
#include "nb64.h"
#endif
// [[Rcpp::plugins(cpp11)]]

using namespace std;

/*
 * Unnamed NameSpace for file scope things.
 */
namespace {
    using namespace DigitalNetNS;
#if !defined(IN_RCPP)
    const string digital_net_path = "DIGITAL_NET_PATH";
    struct digital_net_name {
        std::string name;
        std::string abb;
        std::string construction;
    };

    const digital_net_name digital_net_name_data[] = {
        {"NX", "nx", "Niederreiter-Xing"},
        //{"Sobol", "so", "Sobol"},
        {"Sobol", "sobolbase", "Sobol"},
        {"Old_Sobol", "oldso", "Old Sobol"},
        {"NX_LowWAFOM", "nxlw", "NX+LowWAFOM, CV = (max(CV) + min(CV))/2"},
        {"Sobol_LowWAFOM", "solw", "Sobol+LowWAFOM, CV = (max(CV) + min(CV))/2"}
    };

    const uint32_t digital_net_name_data_size = 5;
#endif

#if !defined(IN_RCPP)
    ostream& errs = cerr;
    void msgout(ostream&)
    {
    }
#endif

#if !defined(IN_RCPP)
    const char * getDataPath()
    {
        return getenv(digital_net_path.c_str());
    }
#endif
#if !defined(IN_RCPP)
    const string makePath(const string& name, const string& ext)
    {
        string path;
        const char * cpath = getDataPath();
        if (cpath == NULL) {
            path = "../data";
        } else {
            path = cpath;
        }
#if __cplusplus >= 201103L
        if (path.back() != '/') {
            path += "/";
        }
#else
        path += "/";
#endif
        path += name;
        path += ext;
        return path;
    }
#endif // !IN_RCPP

#if defined(USE_FILE)
    template<typename U>
    int read_sobol_base(const string& path, uint32_t s, uint32_t m, U base[])
    {
        //ifstream ifs(path, ios::in | ios::binary);
        ifstream ifs;
        ifs.open(path.c_str(), ios::in);
        if (!ifs) {
            errs << "can't open:" << path << endl;
            msgout(errs);
            return -1;
        }
#if defined(IN_CRAN)
        uint64_t *data = new uint64_t[s * m];
#else
        uint64_t data[s * m];
#endif
        bool r = get_sobol_base(ifs, s, m, data);
        if (!r) {
            return -1;
        }
        if (sizeof(U) * 8  == 32) {
            for (size_t i = 0; i < s * m; i++) {
                base[i] = static_cast<U>((data[i] >> 32)
                                         & UINT32_C(0xffffffff));
            }
        } else {
            for (size_t i = 0; i < s * m; i++) {
                base[i] = data[i];
            }
        }
#if defined(IN_CRAN)
        delete[] data;
#endif
        return 0;
    }
#endif // USE_FILE

#if defined(IN_RCPP) && defined(USE_DF) && defined(USE_SOBOL)
    template<typename U>
    int read_sobol_base(DataFrame df, uint32_t s, uint32_t m, U base[])
    {
        uint64_t data[s * m];
        bool r = read_sobol_base(df, s, m, data);
        if (!r) {
            return -1;
        }
        if (sizeof(U) * 8  == 32) {
            for (size_t i = 0; i < s * m; i++) {
                base[i] = static_cast<U>((data[i] >> 32)
                                         & UINT32_C(0xffffffff));
            }
        } else {
            for (size_t i = 0; i < s * m; i++) {
                base[i] = data[i];
            }
        }
        return 0;
    }
#endif

#if !defined(IN_RCPP)
    template<typename U>
    int selectSobolBase(const string& path, uint32_t s, uint32_t m, U base[])
    {
        uint64_t data[s * m];
        //int bitsize = sizeof(U) * 8;
        //bool r = select_sobol_base(path, bitsize, s, m, data);
        bool r = select_sobol_base(path, s, m, data);
        if (!r) {
            return -1;
        }
        if (sizeof(U) * 8  == 32) {
            for (size_t i = 0; i < s * m; i++) {
                base[i] = static_cast<U>((data[i] >> 32)
                                         & UINT32_C(0xffffffff));
            }
        } else {
            for (size_t i = 0; i < s * m; i++) {
                base[i] = data[i];
            }
        }
        return 0;
    }
#endif // IN_RCPP

#if !defined(IN_RCPP)
    template<typename U>
    int read_digital_net_data(std::istream& is, int n,
                              uint32_t s, uint32_t m,
                              U base[],
                              int * tvalue, double * wafom)
    {
        uint64_t data[s * m];
        uint64_t tmp;
        uint32_t i = 0;
        uint32_t j = 0;
        for (i = 0; i < m; i++) {
            for (j = 0; j < s; j++) {
                if (!is) {
                    errs << "too less data i = " << dec << i
                         << " j = " << j
                         << " s = " << s << " m = " << m
                         << endl;
                    msgout(errs);
                    return -1;
                }
                is >> tmp;
                data[i * s + j] = tmp;
            }
        }
        if (i * s + j < s * m) {
            errs << "too less data i = " << dec << i
                 << " s = " << s << " m = " << m
                 << endl;
            msgout(errs);
            return -1;
        }
        if (is) {
            is >> *wafom;
        }
        if (is) {
            is >> *tvalue;
        }
        if (n == sizeof(U) * 8) {
            for (size_t i = 0; i < s * m; i++) {
                base[i] = static_cast<U>(data[i]);
            }
        } else if (n == 32) { // U is uint64_t
            for (size_t i = 0; i < s * m; i++) {
                base[i] = data[i] << 32;
            }
        } else { // n == 64 U is uint32_t
            for (size_t i = 0; i < s * m; i++) {
                base[i] = static_cast<U>((data[i] >> 32)
                                        & UINT32_C(0xffffffff));
            }
        }
        return 0;
    }
#endif

#if defined(IN_RCPP)
#if defined(USE_DF)
    template<typename U>
    int read_digital_net_data(DataFrame df, digital_net_id id,
                              uint32_t s, uint32_t m,
                              U base[],
                              int * tvalue, double * wafom) {
#if defined(USE_SOBOL)
        if (id == SOBOL) {
            return read_sobol_base(df, s, m, base);
        }
#endif
#if defined(DEBUG)
        Rcout << "in read_digital_net_data" << endl;
#endif
        //stop("step 1");
        //int bit = sizeof(U) * 8;
        //int r = 0;
        NumericVector wafom_v = df["wafom"];
        if (wafom_v.length() < 1) {
            *wafom = NAN;
        } else {
            *wafom = wafom_v[0];
        }
        NumericVector tvalue_v = df["tvalue"];
        if (tvalue_v.length() < 1) {
            *tvalue = -1;
        } else {
            *tvalue = tvalue_v[0];
        }
#if defined(DEBUG)
        Rcout << "wafom = " << *wafom << endl;
        Rcout << "tvalue = " << *tvalue << endl;
#endif
        StringVector data = df["data"];
        if (data.length() < 1) {
            stop("not found");
            return -1;
        }
        //stop("step 2");
#if defined(BINARY_DF)
        stringstream ssbase;
        ssbase << data[0];
        string str;
        ssbase >> str;
#if defined(DEBUG)
        Rcout << "data = " << str << endl;
        Rcout << "str[0] = "<< str[0] << endl;
#endif
        //stop("step 3");
        uint64_t * b = new uint64_t[s * m];
        hexchar64tohost(reinterpret_cast<const uint8_t *>(str.c_str()+2),
                        str.size()-3, b);
        //nb64tohost(reinterpret_cast<const uint8_t *>(str.c_str()),
        //           str.size(), b);
        //nb64tohost(str.c_str(), str.size(), base);
        for (size_t i = 0; i < s * m; i++) {
            base[i] = b[i];
        }
        delete[] b;
#else
        stringstream ssbase;
        ssbase << data[0];
        for (size_t i = 0; i < s * m; i++) {
            ssbase >> base[i];
        }
#endif
#if defined(DEBUG)
        Rcout << "out select_digital_net_data" << endl;
        Rcout << "base:" << endl;
        for (size_t i = 0; i < s; i++) {
            Rcout << base[i] << " ";
        }
        Rcout << endl;
#endif
        return 0;
    }
#endif // USE_DF
#endif

#if !defined(IN_RCPP)
    template<typename U>
    int read_digital_net_data(digital_net_id id, uint32_t s, uint32_t m,
                              U base[],
                              int * tvalue, double * wafom) {
#if defined(DEBUG)
        cout << "in read_digital_net_data" << endl;
#endif
        string name = digital_net_name_data[id].abb;
        string path = makePath(name, ".dat");
        if (id == SOBOL) {
            return read_sobol_base(path, s, m, base);
        }
#if defined(DEBUG)
        cout << "fname = " << path << endl;
#endif
        const char * mode = "rb";
        size_t count;
        errno = 0;
        FILE *fp = fopen(path.c_str(), mode);
        if (fp == NULL) {
            cerr << path.c_str() << ":" << strerror(errno) << endl;
            //throw runtime_error("can't open");
            throw "can't open";
        }
        uint64_t dmy;
        count = fread(&dmy, sizeof(uint64_t), 1, fp);
        if (count != 1) {
            cerr << "fail to read magic number" << endl;
            //throw runtime_error("fail to read magic number");
            throw "fail to read magic number";
        }
        if (dmy != DIGITAL_MAGIC) {
            cerr << "magic number mismatch" << endl;
            //throw runtime_error("magic number mismatch");
            throw "magic number mismatch";
        }
        digital_net_header_t header;
        for (;;) {
            count = fread(&header, sizeof(digital_net_header_t), 1, fp);
            if (count != 1) {
                cerr << "fail to read header s = "
                     << dec << s << " m = " << m << endl;
                //throw runtime_error("can't read header");
                throw "can't read header";
            }
            if (header.s > s && header.m > m) {
                cerr << "header s = " << dec << header.s << " m = " << header.m
                     << endl;
                //throw runtime_error("can't find s and m in header");
                throw "can't find s and m in header";
            }
            if (header.s == s && header.m == m) {
                break;
            }
        }
        fseek(fp, header.pos, SEEK_SET);
        size_t size = s * m * sizeof(uint64_t);
        digital_net_data_t dn;
        count = fread(&dn, sizeof(digital_net_data_t), 1, fp);
        if (count != 1) {
            cerr << "fail to read digital net data" << endl;
            //throw runtime_error("fail to read digital net data");
            throw "fail to read digital net data";
        }
        uint64_t data[s * m];
        count = fread(data, size, 1, fp);
        if (count != 1) {
            cout << "fail to read digtal net array" << endl;
            //throw runtime_error("fail to read digtal net array");
            throw "fail to read digtal net array";
        }
#if defined(DEBUG)
        cout << "read_digital_net_data after read" << endl;
        cout << "data:" << endl;
        for (size_t i = 0; i < s; i++) {
            cout << data[i] << " ";
        }
        cout << endl;
#endif
        *wafom = dn.wafom;
        *tvalue = dn.tvalue;
        //uint64_t * data = dn.data;
        if (sizeof(U) * 8  == 32) {
            for (size_t i = 0; i < s * m; i++) {
                base[i] = static_cast<U>((data[i] >> 32)
                                         & UINT32_C(0xffffffff));
            }
        } else {
            for (size_t i = 0; i < s * m; i++) {
                base[i] = data[i];
            }
        }
#if defined(DEBUG)
        cout << "out read_digital_net_data" << endl;
        cout << "base:" << endl;
        for (size_t i = 0; i < s; i++) {
            cout << base[i] << " ";
        }
        cout << endl;
#endif
        return 0;
    }
#endif

#if !defined(IN_RCPP)
    template<typename U>
    int select_digital_net_data(digital_net_id id, uint32_t s, uint32_t m,
                                U base[],
                                int * tvalue, double * wafom) {
#if defined(DEBUG)
        cout << "in select_digital_net_data" << endl;
#endif
        string name = digital_net_name_data[id].abb;
        string path = makePath("digitalnet", ".sqlite3");
        if (id == SOBOL) {
            return selectSobolBase(path, s, m, base);
        }
#if defined(DEBUG)
        cout << "dbname = " << path << endl;
#endif
        int bit = sizeof(U) * 8;
        // prepare sql
        int r = 0;
        sqlite3 *db;
        r = sqlite3_open_v2(path.c_str(), &db, SQLITE_OPEN_READONLY, NULL);
        if (r != SQLITE_OK) {
            cout << "sqlite3_open error code = " << dec << r << endl;
            cout << sqlite3_errmsg(db) << endl;
            return -1;
        }
        string strsql = "select wafom, tvalue, data";
        strsql += " from digitalnet ";
        strsql += " where netname = ? "; // 1
        strsql += "and bitsize = ? ";    // 2
        strsql += "and dimr = ? ";       // 3
        strsql += "and dimf2 = (select min(dimf2) from digitalnet ";
        strsql += "where netname = ? "; // 4
        strsql += "and bitsize = ? ";   // 5
        strsql += "and dimr = ? ";      // 6
        strsql += "and dimf2 >= ?);";   // 7
        sqlite3_stmt* select_sql = NULL;
        stringstream ssbase;
        r = sqlite3_prepare_v2(db, strsql.c_str(), -1, &select_sql, NULL);
        if (r != SQLITE_OK) {
            cout << "sqlite3_prepare error code = " << dec << r << endl;
            cout << sqlite3_errmsg(db) << endl;
            r = sqlite3_close_v2(db);
            return -2;
        }
        if (select_sql == NULL) {
            cout << "sqlite3_prepare null statement" << endl;
            r = sqlite3_close_v2(db);
            return -3;
        }
        do {
            r = sqlite3_bind_text(select_sql, 1, name.c_str(),
                                  -1, SQLITE_STATIC);
            if (r != SQLITE_OK) {
                cout << "error bind netname r = " << dec << r << endl;
                cout << sqlite3_errmsg(db) << endl;
                return r;
            }
            r = sqlite3_bind_int(select_sql, 2, bit);
            if (r != SQLITE_OK) {
                cout << "error bind bitsize r = " << dec << r << endl;
                cout << sqlite3_errmsg(db) << endl;
                return r;
            }
            r = sqlite3_bind_int(select_sql, 3, s);
            if (r != SQLITE_OK) {
                cout << "error bind dimr r = " << dec << r << endl;
                cout << sqlite3_errmsg(db) << endl;
                return r;
            }
            r = sqlite3_bind_text(select_sql, 4, name.c_str(),
                                  -1, SQLITE_STATIC);
            if (r != SQLITE_OK) {
                cout << "error bind netname r = " << dec << r << endl;
                cout << sqlite3_errmsg(db) << endl;
                return r;
            }
            r = sqlite3_bind_int(select_sql, 5, bit);
            if (r != SQLITE_OK) {
                cout << "error bind bitsize r = " << dec << r << endl;
                cout << sqlite3_errmsg(db) << endl;
                return r;
            }
            r = sqlite3_bind_int(select_sql, 6, s);
            if (r != SQLITE_OK) {
                cout << "error bind dimr r = " << dec << r << endl;
                cout << sqlite3_errmsg(db) << endl;
                return r;
            }
            r = sqlite3_bind_int(select_sql, 7, m);
            if (r != SQLITE_OK) {
                cout << "error bind dimf2 r = " << dec << r << endl;
                cout << sqlite3_errmsg(db) << endl;
                return r;
            }
            r = sqlite3_step(select_sql);
            if (r != SQLITE_ROW) {
                cout << "not found" << endl;
                cout << "netname = " << name << endl;
                return r;
            }
            if (sqlite3_column_type(select_sql, 0) == SQLITE_NULL) {
                *wafom = NAN;
            } else {
                *wafom = sqlite3_column_double(select_sql, 0);
            }
            if (sqlite3_column_type(select_sql, 1) == SQLITE_NULL) {
                *tvalue = -1;
            } else {
                *tvalue = sqlite3_column_int(select_sql, 1);
            }
            char * tmp = (char *)sqlite3_column_text(select_sql, 2);
            ssbase << tmp;
        } while (false);
        // release sql
        r = sqlite3_finalize(select_sql);
        if (r != SQLITE_OK) {
            cout << "error finalize r = " << dec << r << endl;
            cout << sqlite3_errmsg(db) << endl;
        }
        sqlite3_close_v2(db);
        if (r != SQLITE_OK) {
            return r;
        }
        for (size_t i = 0; i < s * m; i++) {
            ssbase >> base[i];
        }
#if defined(DEBUG)
        cout << "out select_digital_net_data" << endl;
        cout << "base:" << endl;
        for (size_t i = 0; i < s; i++) {
            cout << base[i] << " ";
        }
        cout << endl;
#endif
        return 0;
    }
#endif

#if !defined(IN_RCPP)
    int get_s_max(const string& path, digital_net_id id)
    {
        // db open
        string name = digital_net_name_data[id].abb;
        sqlite3 *db;
        int r = 0;
        r = sqlite3_open_v2(path.c_str(), &db, SQLITE_OPEN_READONLY, NULL);
        if (r != SQLITE_OK) {
            cout << "sqlite3_open error code = " << dec << r << endl;
            cout << sqlite3_errmsg(db) << endl;
            return -1;
        }
        sqlite3_stmt *select_sql = NULL;
        string strsql = "select max(dimr) from digitalnet ";
        strsql += "where netname = ?;";
        stringstream ssbase;
        r = sqlite3_prepare_v2(db, strsql.c_str(), -1, &select_sql, NULL);
        if (r != SQLITE_OK) {
            cout << "sqlite3_prepare error code = " << dec << r << endl;
            cout << sqlite3_errmsg(db) << endl;
            r = sqlite3_close_v2(db);
            return -2;
        }
        if (select_sql == NULL) {
            cout << "sqlite3_prepare null statement" << endl;
            r = sqlite3_close_v2(db);
            return -3;
        }
        r = sqlite3_bind_text(select_sql, 1, name.c_str(),
                              -1, SQLITE_STATIC);
        if (r != SQLITE_OK) {
            cout << "error bind netname r = " << dec << r << endl;
            cout << sqlite3_errmsg(db) << endl;
            r = sqlite3_close_v2(db);
            return -4;
        }
        r = sqlite3_step(select_sql);
        if (r != SQLITE_ROW) {
            cout << "not found" << endl;
            cout << "netname = " << name << endl;
            r = sqlite3_close_v2(db);
            return -5;
        }
        int s_max = sqlite3_column_int(select_sql, 0);
        r = sqlite3_finalize(select_sql);
        if (r != SQLITE_OK) {
            cout << "error finalize r = " << dec << r << endl;
            cout << sqlite3_errmsg(db) << endl;
        }
        sqlite3_close_v2(db);
        return s_max;
    }

    int get_s_min(const string& path, digital_net_id id)
    {
        // db open
        string name = digital_net_name_data[id].abb;
        sqlite3 *db;
        int r = 0;
        r = sqlite3_open_v2(path.c_str(), &db, SQLITE_OPEN_READONLY, NULL);
        if (r != SQLITE_OK) {
            cout << "sqlite3_open error code = " << dec << r << endl;
            cout << sqlite3_errmsg(db) << endl;
            return -1;
        }
        sqlite3_stmt *select_sql = NULL;
        //r = select_bind(db, &select_sql);
        string strsql = "select min(dimr) from digitalnet ";
        strsql += "where netname = ?;";
        stringstream ssbase;
        r = sqlite3_prepare_v2(db, strsql.c_str(), -1, &select_sql, NULL);
        if (r != SQLITE_OK) {
            cout << "sqlite3_prepare error code = " << dec << r << endl;
            cout << sqlite3_errmsg(db) << endl;
            r = sqlite3_close_v2(db);
            return -2;
        }
        if (select_sql == NULL) {
            cout << "sqlite3_prepare null statement" << endl;
            r = sqlite3_close_v2(db);
            return -3;
        }
        r = sqlite3_bind_text(select_sql, 1, name.c_str(),
                              -1, SQLITE_STATIC);
        if (r != SQLITE_OK) {
            cout << "error bind netname r = " << dec << r << endl;
            cout << sqlite3_errmsg(db) << endl;
            r = sqlite3_close_v2(db);
            return -4;
        }
        r = sqlite3_step(select_sql);
        if (r != SQLITE_ROW) {
            cout << "not found" << endl;
            cout << "netname = " << name << endl;
            r = sqlite3_close_v2(db);
            return -5;
        }
        int s_min = sqlite3_column_int(select_sql, 0);
        r = sqlite3_finalize(select_sql);
        if (r != SQLITE_OK) {
            cout << "error finalize r = " << dec << r << endl;
            cout << sqlite3_errmsg(db) << endl;
        }
        sqlite3_close_v2(db);
        return s_min;
    }

    int get_m_max(const string& path, digital_net_id id, int s)
    {
        // db open
        string name = digital_net_name_data[id].abb;
        sqlite3 *db;
        int r = 0;
        r = sqlite3_open_v2(path.c_str(), &db, SQLITE_OPEN_READONLY, NULL);
        if (r != SQLITE_OK) {
            cout << "sqlite3_open error code = " << dec << r << endl;
            cout << sqlite3_errmsg(db) << endl;
            return -1;
        }
        sqlite3_stmt *select_sql = NULL;
        //r = select_bind(db, &select_sql);
        string strsql = "select max(dimf2) from digitalnet ";
        strsql += "where netname = ? and dimr = ?;";
        stringstream ssbase;
        r = sqlite3_prepare_v2(db, strsql.c_str(), -1, &select_sql, NULL);
        if (r != SQLITE_OK) {
            cout << "sqlite3_prepare error code = " << dec << r << endl;
            cout << sqlite3_errmsg(db) << endl;
            r = sqlite3_close_v2(db);
            return -2;
        }
        if (select_sql == NULL) {
            cout << "sqlite3_prepare null statement" << endl;
            r = sqlite3_close_v2(db);
            return -3;
        }
        r = sqlite3_bind_text(select_sql, 1, name.c_str(),
                              -1, SQLITE_STATIC);
        if (r != SQLITE_OK) {
            cout << "error bind netname r = " << dec << r << endl;
            cout << sqlite3_errmsg(db) << endl;
            r = sqlite3_close_v2(db);
            return -4;
        }
        r = sqlite3_bind_int(select_sql, 2, s);
        if (r != SQLITE_OK) {
            cout << "error bind netname r = " << dec << r << endl;
            cout << sqlite3_errmsg(db) << endl;
            r = sqlite3_close_v2(db);
            return -5;
        }
        r = sqlite3_step(select_sql);
        if (r != SQLITE_ROW) {
            cout << "not found" << endl;
            cout << "netname = " << name << endl;
            r = sqlite3_close_v2(db);
            return -6;
        }
        int m_max = sqlite3_column_int(select_sql, 0);
        r = sqlite3_finalize(select_sql);
        if (r != SQLITE_OK) {
            cout << "error finalize r = " << dec << r << endl;
            cout << sqlite3_errmsg(db) << endl;
        }
        sqlite3_close_v2(db);
        return m_max;
    }

    int get_m_min(const string& path, digital_net_id id, int s)
    {
        // db open
        string name = digital_net_name_data[id].abb;
        sqlite3 *db;
        int r = 0;
        r = sqlite3_open_v2(path.c_str(), &db, SQLITE_OPEN_READONLY, NULL);
        if (r != SQLITE_OK) {
            cout << "sqlite3_open error code = " << dec << r << endl;
            cout << sqlite3_errmsg(db) << endl;
            return -1;
        }
        sqlite3_stmt *select_sql = NULL;
        //r = select_bind(db, &select_sql);
        string strsql = "select min(dimf2) from digitalnet ";
        strsql += "where netname = ? and dimr = ?;";
        stringstream ssbase;
        r = sqlite3_prepare_v2(db, strsql.c_str(), -1, &select_sql, NULL);
        if (r != SQLITE_OK) {
            cout << "sqlite3_prepare error code = " << dec << r << endl;
            cout << sqlite3_errmsg(db) << endl;
            r = sqlite3_close_v2(db);
            return -2;
        }
        if (select_sql == NULL) {
            cout << "sqlite3_prepare null statement" << endl;
            r = sqlite3_close_v2(db);
            return -3;
        }
        r = sqlite3_bind_text(select_sql, 1, name.c_str(),
                              -1, SQLITE_STATIC);
        if (r != SQLITE_OK) {
            cout << "error bind netname r = " << dec << r << endl;
            cout << sqlite3_errmsg(db) << endl;
            r = sqlite3_close_v2(db);
            return -4;
        }
        r = sqlite3_bind_int(select_sql, 2, s);
        if (r != SQLITE_OK) {
            cout << "error bind netname r = " << dec << r << endl;
            cout << sqlite3_errmsg(db) << endl;
            r = sqlite3_close_v2(db);
            return -5;
        }
        r = sqlite3_step(select_sql);
        if (r != SQLITE_ROW) {
            cout << "not found" << endl;
            cout << "netname = " << name << endl;
            r = sqlite3_close_v2(db);
            return -6;
        }
        int m_min = sqlite3_column_int(select_sql, 0);
        r = sqlite3_finalize(select_sql);
        if (r != SQLITE_OK) {
            cout << "error finalize r = " << dec << r << endl;
            cout << sqlite3_errmsg(db) << endl;
        }
        sqlite3_close_v2(db);
        return m_min;
    }
#endif // IN_RCPP
}


namespace DigitalNetNS {
#if !defined(IN_RCPP)
    int getSMax(digital_net_id id)
    {
        string path = makePath("digitalnet", ".sqlite3");
        if (id == SOBOL) {
            return get_sobol_s_max(path);
        } else {
            return get_s_max(path, id);
        }
    }

    int getSMin(digital_net_id id)
    {
        string path = makePath("digitalnet", ".sqlite3");
        if (id == SOBOL) {
            return get_sobol_s_min(path);
        } else {
            return get_s_min(path, id);
        }
    }

    int getMMax(digital_net_id id, int s)
    {
        string path = makePath("digitalnet", ".sqlite3");
        if (id == SOBOL) {
            return get_sobol_m_max(path, s);
        } else {
            return get_m_max(path, id, s);
        }
    }

    int getMMin(digital_net_id id, int s)
    {
        string path = makePath("digitalnet", ".sqlite3");
        if (id == SOBOL) {
            return get_sobol_m_min(path, s);
        } else {
            return get_m_min(path, id, s);
        }
    }

    const string getDigitalNetName(uint32_t index)
    {
        if (index < digital_net_name_data_size) {
            return digital_net_name_data[index].name;
        } else {
            return "";
        }
    }

    int readDigitalNetHeader(std::istream& is, int * n,
                             uint32_t * s, uint32_t * m)
    {
        is >> *n;
        is >> *s;
        is >> *m;
        if (is.good()) {
            return 0;
        } else {
            return -1;
        }
    }
#endif
/**
 * Constructor from input stream
 *
 * File Format:
 * separator: white space, blank char, tab char, cr, lf, etc.
 * the first item: bit size, integer fixed to 64, currently.
 * the second item: s, unsigned integer.
 * the third item: m, unsigned integer.
 * from fourth: s * m number of 64-bit unsigned integers.
 * the last but one: wafom double precision number, optional.
 * the last: t-value, integer, optional.
 * @param is input stream, from where digital net data are read.
 * @exception runtime_error, when can't read data from is.
 */
#if !defined(IN_RCPP)
    int readDigitalNetData(std::istream& is, int n,
                           uint32_t s, uint32_t m,
                           uint64_t base[],
                           int * tvalue, double * wafom)
    {
        return read_digital_net_data(is, n, s, m, base, tvalue, wafom);
    }
    int readDigitalNetData(std::istream& is, int n,
                           uint32_t s, uint32_t m,
                           uint32_t base[],
                           int * tvalue, double * wafom)
    {
        return read_digital_net_data(is, n, s, m, base, tvalue, wafom);
    }
#endif
/**
 * Constructor from reserved data
 *
 * file are searched from environment variable DIGITAL_NET_PATH
 *
 * @param name name of digital net
 * @param s s value
 * @param m m value
 * @exception runtime_error, when can't read data from is.
 */
#if defined(IN_RCPP)
#if defined(USE_FILE)
    int readSobolBase(string& filename, uint32_t s, uint32_t m,
                      uint64_t base[])
    {
        return read_sobol_base(filename, s, m, base);
    }

    int readSobolBase(string& filename, uint32_t s, uint32_t m,
                      uint32_t base[])
    {
        return read_sobol_base(filename, s, m, base);
    }
#endif
#if defined(USE_DF)
    int readDigitalNetData(DataFrame df, digital_net_id id,
                           uint32_t s, uint32_t m,
                           uint64_t base[],
                           int * tvalue, double * wafom)
    {
        //stop("before read_digital_net_data");
        return read_digital_net_data(df, id, s, m, base, tvalue, wafom);
    }
#endif
#else // not IN_RCPP
    int readDigitalNetData(digital_net_id id, uint32_t s, uint32_t m,
                           uint64_t base[],
                           int * tvalue, double * wafom)
    {
        //return read_digital_net_data(id, s, m, base, tvalue, wafom);
        return select_digital_net_data(id, s, m, base, tvalue, wafom);
    }

    int readDigitalNetData(digital_net_id id, uint32_t s, uint32_t m,
                           uint32_t base[],
                           int * tvalue, double * wafom)
    {
        //return read_digital_net_data(id, s, m, base, tvalue, wafom);
        return select_digital_net_data(id, s, m, base, tvalue, wafom);
    }
#endif // IN_RCPP
}
