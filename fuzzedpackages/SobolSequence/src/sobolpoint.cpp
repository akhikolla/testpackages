// Frances Y. Kuo
//
// Email: <f.kuo@unsw.edu.au>
// School of Mathematics and Statistics
// University of New South Wales
// Sydney NSW 2052, Australia
//
// Last updated: 21 October 2008
//
//   You may incorporate this source code into your own program
//   provided that you
//   1) acknowledge the copyright owner in your program and publication
//   2) notify the copyright owner by email
//   3) offer feedback regarding your experience with different direction numbers
//
//
// -----------------------------------------------------------------------------
// Licence pertaining to sobol.cc and the accompanying sets of direction numbers
// -----------------------------------------------------------------------------
// Copyright (c) 2008, Frances Y. Kuo and Stephen Joe
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//
//     * Neither the names of the copyright holders nor the names of the
//       University of New South Wales and the University of Waikato
//       and its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// -----------------------------------------------------------------------------

/*
 * This file is based on sobol.cc written by Frances Y. Kuo.
 * http://web.maths.unsw.edu.au/~fkuo/sobol/
 * http://web.maths.unsw.edu.au/~fkuo/sobol/sobol.cc
 *
 * And modified much by Mutsuo Saito <saito@manieth.com> for R package.
 */

#include <inttypes.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdio>
#include <cerrno>
#include <cstring>
#include <stdexcept>
#if defined(IN_RCPP)
#include <Rcpp.h>
#endif
#if defined(USE_SQL)
#include <sqlite3.h>
#endif
#include "sobolpoint.h"

// [[Rcpp::plugins(cpp11)]]

//#define DEBUG 1
using namespace std;
//#define DEBUG_STEP(x) do { cerr << "debug step " << x << endl;} while(0)
#define DEBUG_STEP(x)

//static const int s_min = 2;
static const int max_data = 50;

#if defined(IN_RCPP)
using namespace Rcpp;
#endif

namespace {
    //int select_bind(sqlite3 *db, sqlite3_stmt** select_sql);
#if defined(USE_FILE)
    bool read_data(istream& is, uint32_t data[]);
#endif
#if defined(USE_DF)
    bool read_data(Rcpp::DataFrame df, int count, uint32_t data[]);
#endif
#if defined(USE_SQL)
    bool read_data(sqlite3_stmt* select_sql, uint32_t data[]);
#endif
}

namespace DigitalNetNS {
#if defined(USE_DF)
    bool read_sobol_base(DataFrame df,
                         uint32_t s, uint32_t m,  uint64_t base[])
    {
        uint32_t D = s + 1;
        uint32_t L = m;
        //uint32_t N = UINT32_C(1) << (m - 1);
        uint32_t col = 0;
        uint64_t V[L + 1];
        for (unsigned i=1;i<=L;i++) {
            V[i] = UINT64_C(1) << (64 - i); // all m's = 1
        }
#if defined(DEBUG)
        cout << "col = " << dec << col << endl;
        for (uint32_t i = 1; i <= L; i++) {
            cout << "V[" << dec << i << "] = " << hex << V[i] << endl;
        }
#endif
        for (uint32_t i = 1; i <= L; i++) {
            base[(i - 1) * s + col] = V[i];
        }
        uint32_t data[max_data];
        for (uint32_t c = 1; c < D - 1; c++) {
            bool success = read_data(df, c - 1, data);
            if (! success) {
                warning("data format error");
                //throw runtime_error("data format error");
                return false;
            }
            col++;
            //uint32_t d_sobol = data[0];
            uint32_t s_sobol = data[0];
            uint32_t a_sobol = data[1];
            uint32_t *m_sobol = &data[1]; // index from 1
#if defined(DEBUG)
            //cout << "d = " << dec << d_sobol << endl;
            cout << "s = " << dec << s_sobol << endl;
            cout << "a = " << dec << a_sobol << endl;
            cout << "L = " << dec << L << endl;
#endif
            if (L <= s_sobol) {
                for (unsigned i=1;i<=L;i++) {
                    V[i] = static_cast<uint64_t>(m_sobol[i]) << (64 - i);
                }
            } else {
                for (unsigned i = 1; i <= s_sobol; i++) {
                    V[i] = static_cast<uint64_t>(m_sobol[i]) << (64 - i);
                }
                for (unsigned i = s_sobol + 1; i <= L; i++) {
                    V[i] = V[i - s_sobol] ^ (V[i - s_sobol] >> s_sobol);
                    // s
                    for (unsigned k=1; k <= s_sobol-1; k++) {
                        V[i] ^= (((a_sobol >> (s_sobol-1-k)) & 1) * V[i-k]);
                    }
                }
            }
#if defined(DEBUG)
            cout << "col = " << dec << col << endl;
            for (uint32_t i = 1; i <= L; i++) {
                cout << "V[" << dec << i << "] = " << hex << V[i] << endl;
            }
#endif
            for (uint32_t i = 1; i <= L; i++) {
                base[(i - 1) * s + col] = V[i];
            }
        }
        for (uint32_t i = L - 1; i >= 1; i--) {
            for (uint32_t j = 0; j < s; j++) {
                base[i * s + j] ^= base[(i - 1) * s + j];
            }
        }
        return true;
    }
#endif

#if defined(USE_SQL)
    bool select_sobol_base(const std::string& path,
                           uint32_t s, uint32_t m,  uint64_t base[])
    {
        // db open
        sqlite3 *db;
        int r = 0;
        DEBUG_STEP(1);
        r = sqlite3_open_v2(path.c_str(), &db, SQLITE_OPEN_READONLY, NULL);
        if (r != SQLITE_OK) {
            cout << "sqlite3_open error code = " << dec << r << endl;
            cout << sqlite3_errmsg(db) << endl;
            return false;
        }
        DEBUG_STEP(2);
        sqlite3_stmt *select_sql = NULL;
        //r = select_bind(db, &select_sql);
        string strsql = "select d, s, a, mi from sobolbase ";
        strsql += "order by d;";
        r = sqlite3_prepare_v2(db, strsql.c_str(), -1, &select_sql, NULL);
        if (r != SQLITE_OK) {
            cout << "sqlite3_prepare error code = " << dec << r << endl;
            cout << sqlite3_errmsg(db) << endl;
            r = sqlite3_close_v2(db);
            return false;
        }
        if (select_sql == NULL) {
            cout << "sqlite3_prepare null statement" << endl;
            r = sqlite3_close_v2(db);
            return false;
        }

        uint32_t D = s + 1;
        uint32_t L = m;
        //uint32_t N = UINT32_C(1) << (m - 1);
        uint32_t col = 0;
        uint64_t V[L + 1];
        for (unsigned i=1;i<=L;i++) {
            V[i] = UINT64_C(1) << (64 - i); // all m's = 1
        }
#if defined(DEBUG)
        cout << "col = " << dec << col << endl;
        for (uint32_t i = 1; i <= L; i++) {
            cout << "V[" << dec << i << "] = " << hex << V[i] << endl;
        }
#endif
        for (uint32_t i = 1; i <= L; i++) {
            base[(i - 1) * s + col] = V[i];
        }
        uint32_t data[max_data];
        for (uint32_t c = 1; c < D - 1; c++) {
            bool success = read_data(select_sql, data);
            if (! success) {
                cerr << "data format error" << endl;
                //throw runtime_error("data format error");
                return false;
            }
            col++;
            //uint32_t d_sobol = data[0];
            uint32_t s_sobol = data[0];
            uint32_t a_sobol = data[1];
            uint32_t *m_sobol = &data[1]; // index from 1
#if defined(DEBUG)
            //cout << "d = " << dec << d_sobol << endl;
            cout << "s = " << dec << s_sobol << endl;
            cout << "a = " << dec << a_sobol << endl;
            cout << "L = " << dec << L << endl;
#endif
            if (L <= s_sobol) {
                for (unsigned i=1;i<=L;i++) {
                    V[i] = static_cast<uint64_t>(m_sobol[i]) << (64 - i);
                }
            } else {
                for (unsigned i = 1; i <= s_sobol; i++) {
                    V[i] = static_cast<uint64_t>(m_sobol[i]) << (64 - i);
                }
                for (unsigned i = s_sobol + 1; i <= L; i++) {
                    V[i] = V[i - s_sobol] ^ (V[i - s_sobol] >> s_sobol);
                    // s
                    for (unsigned k=1; k <= s_sobol-1; k++) {
                        V[i] ^= (((a_sobol >> (s_sobol-1-k)) & 1) * V[i-k]);
                    }
                }
            }
#if defined(DEBUG)
            cout << "col = " << dec << col << endl;
            for (uint32_t i = 1; i <= L; i++) {
                cout << "V[" << dec << i << "] = " << hex << V[i] << endl;
            }
#endif
            for (uint32_t i = 1; i <= L; i++) {
                base[(i - 1) * s + col] = V[i];
            }
        }
        r = sqlite3_finalize(select_sql);
        if (r != SQLITE_OK) {
            cout << "error finalize r = " << dec << r << endl;
            cout << sqlite3_errmsg(db) << endl;
        }
        sqlite3_close_v2(db);
        if (r != SQLITE_OK) {
            return false;
        }
        for (uint32_t i = L - 1; i >= 1; i--) {
            for (uint32_t j = 0; j < s; j++) {
                base[i * s + j] ^= base[(i - 1) * s + j];
            }
        }
        return true;
    }
#endif

#if defined(USE_FILE)
    bool get_sobol_base(std::istream& is,
                        uint32_t s, uint32_t m,  uint64_t base[])
    {
        uint32_t D = s + 1;
        uint32_t L = m;
        //uint32_t N = UINT32_C(1) << (m - 1);
        uint32_t col = 0;
#if defined(IN_CRAN)
        uint64_t * V = new uint64_t[L + 1];
#else
        uint64_t V[L + 1];
#endif
        uint64_t one = 1;
        for (unsigned i=1;i<=L;i++) {
            V[i] = one << (64 - i); // all m's = 1
        }
#if defined(DEBUG)
        cout << "col = " << dec << col << endl;
        for (uint32_t i = 1; i <= L; i++) {
            cout << "V[" << dec << i << "] = " << hex << V[i] << endl;
        }
#endif
        for (uint32_t i = 1; i <= L; i++) {
            base[(i - 1) * s + col] = V[i];
        }
        uint32_t data[max_data];
        for (uint32_t c = 1; c < D - 1; c++) {
            bool success = read_data(is, data);
            if (! success) {
#if defined(IN_RCPP)
                Rcout << "data format error" << endl;
#else
                cerr << "data format error" << endl;
#endif
                //throw runtime_error("data format error");
#if defined(IN_CRAN)
                delete[] V;
#endif
                return false;
            }
            col++;
            //uint32_t d_sobol = data[0];
            uint32_t s_sobol = data[1];
            uint32_t a_sobol = data[2];
            uint32_t *m_sobol = &data[2]; // index from 1
#if defined(DEBUG)
            //cout << "d = " << dec << d_sobol << endl;
            cout << "s = " << dec << s_sobol << endl;
            cout << "a = " << dec << a_sobol << endl;
            cout << "L = " << dec << L << endl;
#endif
            if (L <= s_sobol) {
                for (unsigned i=1;i<=L;i++) {
                    V[i] = static_cast<uint64_t>(m_sobol[i]) << (64 - i);
                }
            } else {
                for (unsigned i = 1; i <= s_sobol; i++) {
                    V[i] = static_cast<uint64_t>(m_sobol[i]) << (64 - i);
                }
                for (unsigned i = s_sobol + 1; i <= L; i++) {
                    V[i] = V[i - s_sobol] ^ (V[i - s_sobol] >> s_sobol);
                    for (unsigned k=1; k <= s_sobol-1; k++) {
                        V[i] ^= (((a_sobol >> (s_sobol-1-k)) & 1) * V[i-k]);
                    }
                }
            }
#if defined(DEBUG)
            cout << "col = " << dec << col << endl;
            for (uint32_t i = 1; i <= L; i++) {
                cout << "V[" << dec << i << "] = " << hex << V[i] << endl;
            }
#endif
            for (uint32_t i = 1; i <= L; i++) {
                base[(i - 1) * s + col] = V[i];
            }
        }
#if defined(IN_CRAN)
        delete[] V;
#endif
        for (uint32_t i = L - 1; i >= 1; i--) {
            for (uint32_t j = 0; j < s; j++) {
                base[i * s + j] ^= base[(i - 1) * s + j];
            }
        }
        return true;
    }

#endif

}

namespace {
#if defined(USE_DF)
    bool read_data(DataFrame df, int count, uint32_t data[])
    {
        NumericVector vs = df["s"];
        NumericVector va = df["a"];
        StringVector vmi = df["mi"];
        if (count >= vs.length()) {
            return false;
        }
        data[0] = vs[count];
        data[1] = va[count];
        //string tmp = vmi[count];
        stringstream ss;
        ss << vmi[count];
        for (int i = 2; ss.good(); i++) {
            ss >> data[i];
        }
        return true;
    }
#endif

#if defined(USE_SQL)
    bool read_data(sqlite3_stmt* select_sql, uint32_t data[])
    {
        DEBUG_STEP(4.37);
        int r = sqlite3_step(select_sql);
        if (r != SQLITE_ROW) {
            return false;
        }
        //sobolbase.d = sqlite3_column_int(select_sql, 0); // d
        data[0] = sqlite3_column_int(select_sql, 1); // s
        data[1] = sqlite3_column_int(select_sql, 2); // a
        char * tmp = (char *)sqlite3_column_text(select_sql, 3); // mi
        stringstream ss;
        ss << tmp;
        for (int i = 2; ss.good(); i++) {
            ss >> data[i];
        }
        return true;
    }
#endif

#if defined(USE_FILE)
#if defined(USE_BINFILE)
    /* read from binary file */
    bool read_data(istream& is, uint32_t data[])
    {
        is.read(reinterpret_cast<char *>(data), sizeof(uint32_t) * 3);
        if (!is) {
            return false;
        }
        is.read(reinterpret_cast<char *>(&data[3]), sizeof(uint32_t) * data[1]);
        if (!is) {
            return false;
        } else {
            return true;
        }
    }
#else // ! USE_BINFILE
    /* read from text file */
    bool read_data(istream& is, uint32_t data[])
    {
        is >> data[0];
        is >> data[1];
        is >> data[2];
        if (!is) {
            return false;
        }
        for (size_t i = 0; i < data[1]; i++) {
            is >> data[i + 3];
        }
        if (!is) {
            return false;
        } else {
            return true;
        }
    }
#endif // USE_BINFILE
#endif
}
