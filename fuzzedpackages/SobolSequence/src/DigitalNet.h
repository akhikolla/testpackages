#pragma once
#ifndef DIGITAL_NET_H
#define DIGITAL_NET_H
/**
 * @file DigitalNet.h
 *
 * @brief DigitalNet class for Quasi Monte-Carlo Method.
 *
 * @note Currently only 64-bit DigitalNet is implemented.
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
#include "grayindex.h"
#include <stdint.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <cerrno>
#include <cmath>
#if defined(IN_RCPP)
#include <Rcpp.h>
#else
// can't use 64bit in CRAN
#include "MersenneTwister64.h"
#endif
// [[Rcpp::plugins(cpp11)]]

namespace DigitalNetNS {

    enum digital_net_id {
        NX = 0,
        SOBOL = 1,
        OLDSO = 2,
        NXLW = 3,
        SOLW = 4,
        RANDOM = -1
    };

    uint32_t getParameterSize();
    const std::string getDigitalNetName(uint32_t index);
    const std::string getDigitalNetConstruction(uint32_t index);
#if defined(IN_RCPP)
#if defined(USE_FILE)
    int readSobolBase(std::string& fname,
                      uint32_t s, uint32_t m,
                      uint64_t base[]);
    int readSobolBase(std::string& fname,
                      uint32_t s, uint32_t m,
                      uint32_t base[]);
#endif
#if defined(USE_DF)
    int readDigitalNetData(Rcpp::DataFrame df, digital_net_id id,
                           uint32_t s, uint32_t m,
                           uint64_t base[],
                           int * tvalue, double * wafom);
    int readDigitalNetData(Rcpp::DataFrame df, digital_net_id id,
                           uint32_t s, uint32_t m,
                           uint32_t base[],
                           int * tvalue, double * wafom);
#endif
#else // not IN_RCPP
    int readDigitalNetData(digital_net_id id, uint32_t s, uint32_t m,
                           uint64_t base[],
                           int * tvalue, double * wafom);
    int readDigitalNetData(digital_net_id id, uint32_t s, uint32_t m,
                           uint32_t base[],
                           int * tvalue, double * wafom);
    int readDigitalNetHeader(std::istream& is, int * n,
                             uint32_t * s, uint32_t * m);
    int readDigitalNetData(std::istream& is, int n,
                           uint32_t s, uint32_t m,
                           uint64_t base[],
                           int * tvalue, double * wafom);
    int readDigitalNetData(std::istream& is, int n,
                           uint32_t s, uint32_t m,
                           uint32_t base[],
                           int * tvalue, double * wafom);
    int getSMax(digital_net_id id);
    int getSMin(digital_net_id id);
    int getMMax(digital_net_id id, int s);
    int getMMin(digital_net_id id, int s);
#endif // IN_RCPP
    template<typename U>
    class DigitalNet {
    private:
        // First of all, forbid copy and assign.
        DigitalNet(const DigitalNet<U>& that);
        DigitalNet<U>& operator=(const DigitalNet<U>&);

    public:
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
        DigitalNet(std::istream& is) {
            //using namespace std;
            int n;
            id = -100;
            int r = readDigitalNetHeader(is, &n, &s, &m);
            if (r != 0) {
                //throw std::runtime_error("data type mismatch!");
                throw "data type mismatch!";
            }
            base = new U[s * m]();
            r = readDigitalNetData(is, n, s, m, base,
                                   &tvalue, &wafom);
            if (r != 0) {
                //throw std::runtime_error("data type mismatch!");
                throw "data type mismatch!";
            }
            shift = NULL;
            point_base = NULL;
            point = NULL;
            count = 0;
            digitalShift = false;
            shiftVector = false;
        }
#endif
#if defined(IN_RCPP)
#if defined(USE_FILE)
        DigitalNet(std::string& filename, digital_net_id& id,
                   uint32_t s, uint32_t m) {
            using namespace std;
            this->s = s;
            this->m = m;
            this->id = static_cast<int>(id);
            base = new U[s * m]();
            int r;
            if (id == SOBOL) {
                r = readSobolBase(filename, s, m, base);
            } else {
                delete[] base;
                base = NULL;
                Rcpp::stop("id mismatch!");
            }
            if (r != 0) {
                //throw runtime_error("data type mismatch!");
                delete[] base;
                base = NULL;
                Rcpp::stop("data type mismatch!");
            }
            shift = NULL;
            point_base = NULL;
            point = NULL;
            count = 0;
            digitalShift = false;
            shiftVector = false;
        }
#endif
#if defined(USE_DF)
        DigitalNet(Rcpp::DataFrame df, digital_net_id& id,
                   uint32_t s, uint32_t m) {
            using namespace std;
            this->s = s;
            this->m = m;
            this->id = static_cast<int>(id);
            base = new U[s * m]();
            int r = readDigitalNetData(df, id, s, m, base,
                                       &tvalue, &wafom);
            if (r != 0) {
                //throw runtime_error("data type mismatch!");
                throw "data type mismatch!";
            }
            shift = NULL;
            point_base = NULL;
            point = NULL;
            count = 0;
            digitalShift = false;
            shiftVector = false;
        }
#endif //
#else // not IN_RCPP
        DigitalNet(const digital_net_id& id, uint32_t s, uint32_t m) {
            using namespace std;
            this->s = s;
            this->m = m;
            this->id = static_cast<int>(id);
            base = new U[s * m]();
            int r = readDigitalNetData(id, s, m, base,
                                       &tvalue, &wafom);
            if (r != 0) {
                //throw runtime_error("data type mismatch!");
                throw "data type mismatch!";
            }
            shift = NULL;
            point_base = NULL;
            point = NULL;
            count = 0;
            digitalShift = false;
            shiftVector = false;
        }
#endif // IN_RCPP

        ~DigitalNet() {
            if (base != NULL) {
                delete[] base;
            }
            if (shift != NULL) {
                delete[] shift;
            }
            if (point_base != NULL) {
                delete[] point_base;
            }
            if (point != NULL) {
                delete[] point;
            }
        }

        U getBase(int i, int j) const {
            return base[i * s + j];
        }

        // このあとpoint initialize せよ
        void saveBase(U save[], size_t size) const {
            for (size_t i = 0; (i < s * m) && (i < size); i++) {
                save[i] = base[i];
            }
        }

        void restoreBase(U save[], size_t size) const {
            for (size_t i = 0; (i < s * m) && (i < size); i++) {
                base[i] = save[i];
            }
        }

        double getPoint(int i) const {
            return point[i];
        }

        void setDigitalShift(bool value) {
            digitalShift = value;
        }

        void setDigitalShift(U value[]) {
            shiftVector = true;
            if (shift == NULL) {
                shift = new U[s];
            }
            for (size_t i = 0; i < s; i++) {
                shift[i] = value[i];
            }
        }

        const double * getPoint() const {
            return point;
        }

        const U * getPointBase() const {
            return point_base;
        }

        uint32_t getS() const {
            return s;
        }
        uint32_t getM() const {
            return m;
        }
        const std::string getName() {
            if (id >= 0) {
                return getDigitalNetName(id);
            } else {
                return "no name";
            }
        }

#if !defined(IN_RCPP)
        // Random Linear Scramble
        // Base を変えてしまう => いいのかも。
        void scramble() {
            const size_t N = sizeof(U) * 8;
            U LowTriMat[N];
            U tmp;
            const U one = 1;
            for (size_t i = 0; i < s; i++) {
                // 正則な下三角行列を作る
                for (size_t j = 0; j < N; j++) {
                    U p2 = one << (N - j - 1);
                    LowTriMat[j] = (mt() << (N - j - 1)) | p2;
                }
                for (size_t k = 0; k < m; k++) {
                    tmp = 0;
                    for (size_t j = 0; j < N; j++) {
                        U bit = innerProduct(LowTriMat[j], getBase(k, i));
                        tmp ^= bit << (N - j - 1);
                    }
                    setBase(k, i, tmp);
                }
            }
        }

        /** Hill Climb Linear Scramble.
         *
         *
         指定されたi(0<=i<s)に対し, C_iのみにL_iをかける操作をする
         ただし, L_iは下三角行列かつ正則で, 指定されたj, l(n>j>l>=0)に対し
         第(j, l)成分(と対角成分)が1, その他が0の行列である.
         もう一度同じL_iをC_iにかけることで元のC_iに戻る.
         L_i(j, l)をC_iにかけると, C_iのj行にl行を足した(XOR)ものとなる.

         * upos1 > upos2
         * @param idx C_i
         * @param upos1 j count from MSB 0 is MSB
         * @param upos2 l count from MSB 0 is MSB
         */
        void hc_scramble(int idx, int upos1, int upos2) {
            // getBit(a, b) : a の下から b 番目のbit
            const U one = 1;
            const int N = sizeof(U) * 8;
            //int diff = upos2 - upos1;
            int bpos1 = N - 1 - upos1;
            U umask1 = one << bpos1;
            int bpos2 = N - 1 - upos2;
            U umask2 = one << bpos2;
            for (size_t i = 0; i < m; i++) {
                int index = getIndex(i, idx);
                if (base[index] & umask2) {
                    base[index] ^= umask1;
                }
            }
        }
#endif
        void pointInitialize() {
#if defined(DEBUG)
            using namespace std;
            cout << "in pointInitialize" << endl;
#endif
            if (sizeof(U) * 8 == 64) {
                get_max = 64 - 53;
                factor = exp2(-53);
                eps = exp2(-64);
            } else {
                get_max = 0;
                factor = exp2(-32);
                eps = exp2(-33);
            }
            if (shift == NULL) {
                shift = new U[s]();
            }
            if (point_base == NULL) {
                point_base = new U[s]();
            }
            if (point == NULL) {
                point = new double[s]();
            }
            for (uint32_t i = 0; i < s; ++i) {
                point_base[i] = 0;
            }
#if !defined(IN_RCPP)
            if (digitalShift) {
                for (uint32_t i = 0; i < s; ++i) {
                    shift[i] = mt();
                }
            } else {
                for (uint32_t i = 0; i < s; ++i) {
                    shift[i] = 0;
                }
            }
#endif
            gray.clear();
            count = 0;
            count++;
            convertPoint();
#if defined(DEBUG)
            cout << "out pointInitialize" << endl;
#endif
        }

        void nextPoint() {
#if defined(DEBUG)
            using namespace std;
            cout << "in nextPoint" << endl;
#endif
            uint64_t count_max = 1;
            count_max = count_max << m;
            if (count == count_max) {
                pointInitialize();
                //return;
            }
            int bit = gray.index();
#if defined(DEBUG)
            cout << "bit = " << bit << endl;
            cout << "before boint_base:" << endl;
            for (size_t i = 0; i < s; i++) {
                cout << point_base[i] << " ";
            }
            cout << endl;
#endif
            for (uint32_t i = 0; i < s; ++i) {
                point_base[i] ^= getBase(bit, i);
            }
            convertPoint();
            if (count == count_max) {
                count = 0;
                gray.clear();
            } else {
                gray.next();
                count++;
            }
#if defined(DEBUG)
            cout << "after boint_base:" << endl;
            for (size_t i = 0; i < s; i++) {
                cout << point_base[i] << " ";
            }
            cout << endl;
            cout << "out nextPoint" << endl;
#endif
        }
#if !defined(IN_RCPP)
        //void showStatus(std::ostream& os);
        void setSeed(U seed) {
            mt.seed(seed);
        }
#endif
        double getWAFOM() {
            return wafom;
        }

        int64_t getTvalue() {
            return tvalue;
        }
    private:
        void setBase(int i, int j, U value) {
            base[i * s + j] = value;
        }
        int getIndex(int i, int j) const {
            return i * s + j;
        }
        void convertPoint() {
            for (uint32_t i = 0; i < s; i++) {
                // shift して1を立てている
                uint64_t tmp = (point_base[i] ^ shift[i]) >> get_max;
                point[i] = static_cast<double>(tmp) * factor + eps;
            }
        }
        int id;
        uint32_t s;
        uint32_t m;
        uint64_t count;
        double wafom;
        int tvalue;
        bool digitalShift;
        bool shiftVector;
        int get_max;
        double factor;
        double eps;
        GrayIndex gray;
#if !defined(IN_RCPP)
        MersenneTwister64 mt;
#endif
        U * base;
        U * point_base;
        U * shift;
        double * point;
    };


    template<typename T>
    void print(std::ostream& os, const DigitalNet<T>& dn,
               bool verbose = true, char delim = '\n')
    {
        using namespace std;
        int s = dn.getS();
        int m = dn.getM();
        os << dec;
        if (verbose) {
            os << "n = " << (sizeof(T) * 8) << endl;
            os << "s = " << s << endl;
            os << "m = " << m << endl;
            for (int k = 0; k < m; k++) {
                for (int i = 0; i < s; i++) {
                    os << "base[" << dec << k << "][" << i << "]"
                       << hex << dn.getBase(k, i) << endl;
                }
            }
        } else {
            os << (sizeof(T) * 8) << delim;
            os << s << delim;
            os << m << delim;
            for (int k = 0; k < m; k++) {
                for (int i = 0; i < s; i++) {
                    os << dec << dn.getBase(k, i) << delim;
                }
                if (delim != ' ') {
                    os << delim;
                }
            }
        }
    }
}

#endif // DIGITAL_NET_HPP
