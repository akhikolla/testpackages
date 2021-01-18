// yakmo -- yet another k-means via orthogonalization
//  $Id: yakmo.h 1866 2015-01-21 10:25:43Z ynaga $
// Copyright (c) 2012-2015 Naoki Yoshinaga <ynaga@tkl.iis.u-tokyo.ac.jp>



#include <getopt.h>
#include <stdint.h>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <limits>
#include <vector>
#include <algorithm>
#include <iomanip>


#include <sys/types.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>


#include "errx_dummy.h"
#include <Rcpp.h>
using namespace Rcpp;


// fix missing posix getline :/
#if defined(WIN32) or defined(__sun)

/* getline.c -- Replacement for GNU C library function getline

Copyright (C) 1993 Free Software Foundation, Inc.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as
published by the Free Software Foundation; either version 2 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. */

/* Written by Jan Brittenson, bson@gnu.ai.mit.edu.  */

/* Read up to (and including) a TERMINATOR from STREAM into *LINEPTR
   + OFFSET (and null-terminate it). *LINEPTR is a pointer returned from
   malloc (or NULL), pointing to *N characters of space.  It is realloc'd
   as necessary.  Return the number of characters read (not including the
   null terminator), or -1 on error or EOF.  */

int getstr (char ** lineptr, unsigned long int *n, FILE * stream, char terminator, int offset) {
    int nchars_avail;       /* Allocated but unused chars in *LINEPTR.  */
    char *read_pos;     /* Where we're reading into *LINEPTR. */
    int ret;

    if (!lineptr || !n || !stream)
        return -1;

    if (!*lineptr) {
        *n = 64;
        *lineptr = (char *) malloc (*n);

        if (!*lineptr)
            return -1;
    }

    nchars_avail = *n - offset;
    read_pos = *lineptr + offset;

    for (;;) {
        register int c = getc (stream);

        /* We always want at least one char left in the buffer, since we
        always (unless we get an error while reading the first char)
         NUL-terminate the line buffer.  */

        assert (*n - nchars_avail == read_pos - *lineptr);

        if (nchars_avail < 1) {
            if (*n > 64)
                *n *= 2;
            else
                *n += 64;

            nchars_avail = *n + *lineptr - read_pos;
            *lineptr = (char *) realloc (*lineptr, *n);

            if (!*lineptr)
                return -1;

            read_pos = *n - nchars_avail + *lineptr;
            assert (*n - nchars_avail == read_pos - *lineptr);
        }

        if (c == EOF || ferror (stream)) {
            /* Return partial line, if any.  */
            if (read_pos == *lineptr)
                return -1;
            else
                break;
        }

        *read_pos++ = c;
        nchars_avail--;

        if (c == terminator)
            /* Return the line.  */
            break;
    }

    /* Done - NUL terminate and return the number of chars read.  */
    *read_pos = '\0';

    ret = read_pos - (*lineptr + offset);
    return ret;
}

signed long int getline (char **lineptr, unsigned long int *n, FILE *stream) {
    return getstr (lineptr, n, stream, '\n', 0);
}
#endif


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#ifndef USE_HASH
#include <set>
#endif

#define YAKMO_COPYRIGHT  "yakmo - yet another k-means via orthogonalization\n\
Copyright (c) 2012-2013 Naoki Yoshinaga, All rights reserved.\n\
\n\
Usage: %s [options] train model test\n\
\n\
train     train file        set '-' to skip clustering\n\
model     model file        set '-' to train/test w/o saving centroids\n\
test      test  file        set '-' to output clustering results of train\n\
\n"

#define YAKMO_OPT  "Optional parameters in training and testing:\n\
  -t, --dist-type=TYPE      select type of distance function\n\
                            * 0 - Euclidean\n\
  -c, --init-centroid=TYPE  select method of choosing initial centroids\n\
                              0 - random\n\
                            * 1 - k-means++\n\
  -k, --num-cluster=NUM     k in k-means (3)\n\
  -m, --num-result=NUM      number of alternative results (1)\n\
  -i, --iteration=NUM       maximum number of iterations per clustering (100)\n\
  -r, --init-random-seed    initialize random seed in initializing centroids\n\
  -n, --normalize           normalize L2-norm of data points\n\
  -O, --output=TYPE         select output type of testing\n\
                            * 0 - no output\n\
                              1 - report assignment per cluster\n\
                              2 - report assignment per item\n\
  -v, --verbose             verbosity level (1)\n\
  -h, --help                show this help and exit\n"

static const  char*  yakmo_short_options = "t:c:k:m:i:rnO:v:h";

static struct option yakmo_long_options[] = {
    {"dist-type",        required_argument, NULL, 't'},
    {"kmeans++",         required_argument, NULL, 'c'},
    {"num-cluster",      required_argument, NULL, 'k'},
    {"num-result",       required_argument, NULL, 'm'},
    {"iteration",        required_argument, NULL, 'i'},
    {"normalize",        no_argument,       NULL, 'n'},
    {"init-random-seed", no_argument,       NULL, 'r'},
    {"output",           required_argument, NULL, 'O'},
    {"verbose",          no_argument,       NULL, 'v'},
    {"help",             no_argument,       NULL, 'h'},
    {NULL, 0, NULL, 0}
};

extern char* optarg;
extern int   optind;

namespace yakmo {
    typedef long unsigned int uint;

    typedef double fl_t;

    static inline bool getLine (FILE*& fp, char*& line, unsigned long int& read) {
#ifdef __APPLE__

        if ( (line = fgetln (fp, &read)) == NULL) return false;
#elif __sun
        // problems with 64bit and size_t :/
        static signed long int read_ = 0;
        static unsigned long int size = 0; // static helps inlining

        if ( (read_ = getline (&line, &size, fp)) == -1) return false;

        read = read_;

#elif WIN32
        // problems with 64bit and size_t :/
        static signed long int read_ = 0;
        static unsigned long int size = 0; // static helps inlining

        if ( (read_ = getline (&line, &size, fp)) == -1) return false;

        read = read_;
#else
        static signed long int read_ = 0;
        static size_t size = 0; // static helps inlining

        if ( (read_ = getline (&line, &size, fp)) == -1) return false;

        read = read_;
#endif
        * (line + read - 1) = '\0';
        return true;
    }
    static bool isspace (const char p) {
        return p == ' ' || p == '\t';
    }
    template <typename T> T strton (const char* s, char** error) {
        const signed long int  ret  = (signed long int) (std::strtoll (s, error, 10));
        const unsigned long int retu = (unsigned long int) (std::strtoull (s, error, 10));

        if (std::numeric_limits <T>::is_specialized &&
                (ret  < static_cast <signed long int> (std::numeric_limits <T>::min ()) ||
                 retu > static_cast <unsigned long int> (std::numeric_limits <T>::max ())))
            errx (1, "overflow: %s", s);

        return static_cast <T> (ret);
    }
    enum dist_t { EUCLIDEAN };
    enum init_t { RANDOM, KMEANSPP };
    struct option { // option handler
        enum mode_t { BOTH, TRAIN, TEST };
        const char* com, *train, *model, *test;
        //
        dist_t   dist;  // dist-type
        init_t   init;  //
        mutable uint  k;
        mutable uint  m;
        uint     iter;
        bool     random;
        bool     normalize;
        uint         output;
        uint     verbosity;
        mode_t   mode;
        option (int argc, char** argv) : com (argc ? argv[0] : "--"), train ("-"), model ("-"), test ("-"), dist (EUCLIDEAN), init (KMEANSPP), k (3), m (1), iter (100), random (false), normalize (false), output (0), verbosity (0), mode (BOTH) {
            set (argc, argv);
        }
        void set (int argc, char** argv) { // getOpt
            if (argc == 0) return;

            optind = 1;

            while (1) {
                int opt = getopt_long (argc, argv,
                                       yakmo_short_options, yakmo_long_options, NULL);

                if (opt == -1) break;

                char* err = NULL;

                switch (opt) {
                    case 't':
                        dist      = strton <dist_t> (optarg, &err);
                        break;

                    case 'c':
                        init      = strton <init_t> (optarg, &err);
                        break;

                    case 'k':
                        k         = strton <uint> (optarg, &err);
                        break;

                    case 'm':
                        m         = strton <uint> (optarg, &err);
                        break;

                    case 'i':
                        iter      = strton <uint> (optarg, &err);
                        break;

                    case 'r':
                        random    = true;
                        break;

                    case 'n':
                        normalize = true;
                        break;

                    case 'O':
                        output    = strton <uint16_t> (optarg, &err);
                        break;

                    case 'v':
                        verbosity = strton <uint> (optarg, &err);
                        break;

                    // misc
                    case 'h':
                        printCredit ();
                        printHelp ();
                        Rcpp::stop ("Printed help, exiting.");

                    default:
                        printCredit ();
                        Rcpp::stop ("Printed credits, exiting.");
                }

                if (err && *err)
                    errx (1, "unrecognized option value: %s", optarg);
            }

            if (dist != EUCLIDEAN)
                errx (1, "only euclidean distance is supported.");

            if (init != RANDOM && init != KMEANSPP)
                errx (1, "unsupported centroid initialization.");

            if (argc < optind + 3) {
                printCredit ();
                errx (1, "Type `%s --help' for option details.", com);
            }

            train = argv[optind];
            model = argv[++optind];
            test  = argv[++optind];
            setMode (); // induce appropriate mode
        }
        void setMode () {
            if (std::strcmp (train, "-") == 0 && std::strcmp (test, "-") == 0)
                errx (1, "specify at least training or test file.");
            else if (std::strcmp (test,  "-") == 0) mode = TRAIN;
            else if (std::strcmp (train, "-") == 0) mode = TEST;
            else                                    mode = BOTH;

            if (std::strcmp (model, "-") == 0 && mode == TEST)
                errx (1, "instant mode needs training files.");

            const char* mode0 [] = {"BOTH", "TRAIN", "TEST"};
            Rcout << "mode:" << mode0[mode] << "\n";
        }
        void printCredit () {
            Rcout << YAKMO_COPYRIGHT << com;
        }
        void printHelp () {
            Rcout << YAKMO_OPT;
        }
    };
// implementation of space-efficient k-means using triangle inequality:
//   G. Hamerly. Making k-means even faster (SDM 2010)
    class kmeans {
        public:
#pragma pack(1)
            struct node_t {
                uint idx;
                fl_t val;
                node_t () : idx (0), val (0) {}
                node_t (uint idx_, fl_t val_) : idx (idx_), val (val_) {}
                bool operator< (const node_t &n) const {
                    return idx < n.idx;
                }
            };
#pragma pack()
            class centroid_t;
            class point_t {
                public:
                    point_t (const node_t* n, const uint size, const fl_t norm)
                        : up_d (0), lo_d (0), id (), _size (size), _body (new node_t[_size]), _norm (norm) {
                        std::copy (n, n + size, body ());
                    }
                    point_t& operator= (const point_t& p) {
                        up_d = p.up_d;
                        lo_d = p.lo_d;
                        id = p.id;
                        _size = p._size;
                        _body = p._body;
                        _norm = p._norm;
                        return *this;
                    }
                    fl_t calc_ip (const centroid_t& c) const {
                        // return inner product between this point and the given centroid
                        fl_t ret = 0;

                        for (const node_t* n = begin (); n != end (); ++n)
                            ret += n->val * c[n->idx];

                        return ret;
                    }
                    fl_t calc_dist (const centroid_t& c, const dist_t dist) const {
                        // return distance from this point to the given centroid
                        fl_t ret = 0;

                        switch (dist) {
                            case EUCLIDEAN:
                                ret += _norm + c.norm ();

                                for (const node_t* n = begin (); n != end (); ++n)
                                    ret -= 2 * n->val * c[n->idx];
                        }

                        return  ret;
                    }
                    void set_closest (const std::vector <centroid_t> &cs, const dist_t dist) {
                        uint i   = id == 0 ? 1 : 0; // second closest (cand)
                        uint id0 = id;
                        fl_t d0 (calc_dist (cs[id0], dist)), d1 (calc_dist (cs[i], dist));

                        if (d1 < d0) {
                            id = i;
                            std::swap (d0, d1);
                        }

                        for (++i; i < cs.size (); ++i) { // for all other centers
                            if (i == id0) continue;

                            const fl_t di = calc_dist (cs[i], dist);

                            if (di < d0) {
                                d1 = d0;
                                d0 = di;
                                id = i;
                            } else if (di < d1) {
                                d1 = di;
                            }
                        }

                        up_d = std::sqrt (d0);
                        lo_d = std::sqrt (d1);
                    }
                    void shrink (const uint nf) {
                        while (! empty () && back ().idx > nf) --_size;
                    }
                    void project (const centroid_t& c) {
                        const double norm_ip = calc_ip (c) / c.norm ();
                        up_d = lo_d = id = 0;
                        _norm = 0; // reset

                        for (uint i = 0; i < _size; ++i) {
                            _body[i].val -= c[_body[i].idx] * norm_ip;
                            _norm += _body[i].val * _body[i].val;
                        }

                        _norm = floor (_norm * 1000000000000) / 1000000000000;
                    }
                    const node_t* begin () const {
                        return _body;
                    }
                    const node_t* end () const {
                        return _body + _size;
                    }
                    fl_t    norm () const {
                        return _norm;
                    }
                    uint    size () const {
                        return _size;
                    }
                    bool    empty () const {
                        return _size == 0;
                    }
                    node_t& back () const {
                        return _body[_size - 1];
                    }
                    node_t* body () const {
                        return _body;
                    }
                    void    clear () const {
                        if (_body) delete [] _body;
                    }
                    fl_t    up_d;  // distance to the closest centroid
                    fl_t    lo_d;  // distance to the second closest centroid
                    uint    id;    // cluster id
                private:
                    uint    _size;
                    node_t* _body;
                    fl_t    _norm;
            };
            class centroid_t {
                public:
                    centroid_t (point_t& p, const uint nf, const bool delegate = false) :
                        delta (0), next_d (0), _norm (p.norm ()), _dv (0), _sum (0), _body (0), _nelm (0), _nf (nf), _size (0) {
                        if (delegate) {
                            _size = p.size ();
                            _body = p.body (); // delegate
                        } else { // workaround for a bug in value initialization in gcc 4.0
                            _dv  = new fl_t[_nf + 1];
                            std::fill_n (_dv,  _nf + 1, 0);
                            _sum = new fl_t[_nf + 1];
                            std::fill_n (_sum, _nf + 1, 0);

                            for (const node_t* n = p.begin (); n != p.end (); ++n)
                                _dv[n->idx] = n->val;
                        }
                    }
                    fl_t operator[] (const uint i) const {
                        return _dv[i];
                    }
                    void pop (const point_t& p) {
                        for (const node_t* n = p.begin (); n != p.end (); ++n)
                            _sum[n->idx] -= n->val;

                        --_nelm;
                    }
                    void push (const point_t& p) {
                        for (const node_t* n = p.begin (); n != p.end (); ++n)
                            _sum[n->idx] += n->val;

                        ++_nelm;
                    }
                    fl_t calc_dist (const centroid_t& c, const dist_t dist, const bool skip = true) const {
                        // return distance from this centroid to the given centroid
                        fl_t ret = 0;

                        switch (dist) {
                            case EUCLIDEAN:
                                if (skip) {
                                    const fl_t cand = next_d * next_d;

                                    for (uint d = 0; d <= _nf; ++d)
                                        if ( (ret += (_dv[d] - c[d]) * (_dv[d] - c[d])) > cand) break;
                                } else
                                    for (uint d = 0; d <= _nf; ++d)
                                        ret += (_dv[d] - c[d]) * (_dv[d] - c[d]);
                        }

                        return ret;
                    }
                    void set_closest (const std::vector <centroid_t>& centroid, const dist_t dist) {
                        uint i = (this == &centroid[0]) ? 1 : 0;
                        next_d = calc_dist (centroid[i], dist, false);

                        for (++i; i < centroid.size (); ++i) {
                            if (this == &centroid[i]) continue;

                            const fl_t di = calc_dist (centroid[i], dist);

                            if (di < next_d) next_d = di;
                        }

                        next_d = std::sqrt (next_d);
                    }
                    void reset (const dist_t dist) { // move center
                        delta = _norm = 0;

                        switch (dist) {
                            case EUCLIDEAN:
                                for (uint i = 0; i <= _nf; ++i) {
                                    const fl_t v = _sum[i] / (fl_t) (_nelm);
                                    delta += (v - _dv[i]) * (v - _dv[i]);
                                    _norm += v * v;
                                    _dv[i] = v;
                                }
                        }

                        delta = std::sqrt (delta);
                    }
                    void compress () {
                        _size = 0;

                        for (uint i = 0; i <= _nf; ++i)
                            if (std::fpclassify (_dv[i]) != FP_ZERO)
                                ++_size;

                        _body = new node_t[_size];

                        for (uint i (0), j (0); i <= _nf; ++i)
                            if (std::fpclassify (_dv[i]) != FP_ZERO)
                                _body[j].idx = i, _body[j].val = _dv[i], ++j;

                        delete [] _dv;
                        _dv  = 0;
                        delete [] _sum;
                        _sum = 0;
                    }
                    void decompress () {
                        _dv = new fl_t[_nf + 1];
                        std::fill_n (_dv, _nf + 1, 0);

                        for (uint i = 0; i < _size; ++i)
                            _dv[_body[i].idx] = _body[i].val;

                        delete [] _body;
                        _body = 0;
                    }
                    void print (FILE* fp, const uint j) const {
                        std::fprintf (fp, "%d", j);

                        for (uint i = 0; i < _size; ++i)
                            std::fprintf (fp, " %d:%.16g", _body[i].idx, _body[i].val);

                        std::fprintf (fp, "\n");
                    }
                    /// ---
                    std::vector<fl_t> print () const {
                        std::vector<fl_t> r (_size);

                        for (uint i = 0; i < _size; ++i) {
                            r[_body[i].idx] = _body[i].val;
                        }

                        return (r);
                    }
                    /// ---
                    fl_t norm () const {
                        return _norm;
                    }
                    void clear () {
                        if (_dv)   delete [] _dv;

                        if (_sum)  delete [] _sum;

                        if (_body) delete [] _body;
                    }
                    //
                    fl_t     delta;  // moved distance
                    fl_t     next_d; // distance to neighbouring centroind
                private:
                    fl_t     _norm;  // norm
                    fl_t*    _dv;
                    fl_t*    _sum;
                    node_t*  _body;
                    uint     _nelm;  // # elements belonging to the cluster
                    uint     _nf;    // # features
                    uint     _size;  // # nozero features
            };
            kmeans (const option &opt) : _opt (opt), _point (), _centroid (), _body (), _nf (0) {
                _centroid.reserve (_opt.k);
            }
            ~kmeans () {
                clear_point ();
                clear_centroid ();
            }
            void clear_point () {
                for (std::vector <point_t>::iterator it = _point.begin ();
                        it != _point.end (); ++it)
                    it->clear ();

                std::vector <point_t> ().swap (_point);
            }
            void clear_centroid () {
                for (std::vector <centroid_t>::iterator it = _centroid.begin ();
                        it != _centroid.end (); ++it)
                    it->clear ();

                std::vector <centroid_t> ().swap (_centroid);
            }
            std::vector <point_t>&    point () {
                return _point;
            }
            std::vector <centroid_t>& centroid () {
                return _centroid;
            }
            static point_t read_point (char* const ex, const char* const ex_end, std::vector <node_t>& tmp, const bool normalize = false) {
                tmp.clear ();
                fl_t norm = 0;
                char* p = ex;

                while (p != ex_end) {
                    signed long int fi = 0;

                    for (; *p >= '0' && *p <= '9'; ++p) {
                        fi *= 10, fi += *p, fi -= '0';

                        if (fi > std::numeric_limits <signed long int>::max ())
                            errx (1, "overflow: %s", ex);
                    }

                    if (*p != ':') errx (1, "illegal feature index: %s", ex);

                    ++p;
                    const fl_t v = (fl_t) (std::strtod (p, &p));
                    tmp.push_back (node_t ( (uint) (fi), v));
                    norm += v * v;

                    while (isspace (*p)) ++p;
                }

                std::sort (tmp.begin (), tmp.end ());

                // this seems completly ok  on win vs linux.
                // Rcout << std::setprecision(16) << norm << ",";
                if (normalize) { // normalize
                    norm = std::sqrt (norm);

                    for (std::vector <node_t>::iterator it = tmp.begin ();
                            it != tmp.end (); ++it)
                        it->val /= norm;

                    norm = 1.0;
                }

                return point_t (&tmp[0], static_cast <uint> (tmp.size ()), norm); // expect RVO
            }
            void set_point (char* ex, char* ex_end, const bool normalize) {
                _point.push_back (read_point (ex, ex_end, _body, normalize));

                if (! _point.back ().empty ())
                    _nf = std::max (_point.back ().back(). idx, _nf);
            }
            void delegate (kmeans* km) {
                std::swap (_point, km->point ());
                km->nf () = _nf;
            }
            void compress () {
                for (std::vector <centroid_t>::iterator it = _centroid.begin ();
                        it != _centroid.end (); ++it)
                    it->compress ();
            }
            void decompress () {
                for (std::vector <centroid_t>::iterator it = _centroid.begin ();
                        it != _centroid.end (); ++it)
                    it->decompress ();
            }
            void push_centroid (point_t& p, const bool delegate = false) {
                _centroid.push_back (centroid_t (p, _nf, delegate));
            }
            // implementation of fast k-means:
            //   D. Arthur and S. Vassilvitskii. k-means++: the advantages of careful seeding. SODA (2007)
            void init () {
                struct rng_t {
                    // for compatibility we keep that random thing
                    rng_t (const bool r) {}

                    fl_t operator () () {
                        fl_t p = R::runif (0, 1);
                        // random numbers: 100% ok between 32/64bit and linux/windows
                        // Rcout << std::setprecision(17) << p << ", ";
                        return (p);
                    }
                } rng (_opt.random);
                //
#ifdef USE_HASH
                std::unordered_set <uint> chosen;
#else
                std::set <uint> chosen;
#endif
                std::vector <fl_t> r;
                fl_t obj = 0;

                if (_opt.init == KMEANSPP) r.resize (_point.size (), 0);

                for (uint i = 0; i < _opt.k; ++i) {
                    uint c = 0;

                    do {
                        switch (_opt.init) {
                            case RANDOM:
                                c = (std::floor (rng () * _point.size ()));
                                break;

                            case KMEANSPP:
                                c = (i == 0 ?
                                     std::floor (rng () * _point.size ()) :
                                     std::distance (r.begin (),
                                                    std::lower_bound (r.begin (), r.end (), obj * rng ())));
                                break;
                        }

                        // skip chosen centroids; fix a bug reported by Gleb
                        while (chosen.find (c) != chosen.end ()) {
                            c = c < _point.size () - 1 ? c + 1 : 0;
						}
                    } while (chosen.find (c) != chosen.end ());

                    if (c >= _point.size()) {
                        // something did not work out
                        c = 0;
                    }

                    push_centroid (_point[c]);
                    obj = 0;
                    chosen.insert (c);

                    for (uint j = 0; j < _point.size (); ++j) {
                        point_t& p = _point[j];
                        const fl_t di = p.calc_dist (_centroid[i], _opt.dist);

                        if (i == 0 || di < p.up_d) {    // closest
                            p.lo_d = p.up_d;
                            p.up_d = di;
                            p.id = i;
                        } else if (i == 1 || di < p.lo_d) { // second closest
                            p.lo_d = di;
                        }

                        if (i < _opt.k - 1) {
                            if (_opt.init == KMEANSPP) {
                                obj += p.up_d;
                                r[j] = obj;
                            }
                        } else { // i == _k - 1
                            p.up_d = std::sqrt (p.up_d);
                            p.lo_d = std::sqrt (p.lo_d);
                            _centroid[p.id].push (p);
                        }
                    }
                }

                if (_opt.verbosity > 1) Rcout << "\n";
            }
            void update_bounds () {
                uint id0 (0), id1 (1);

                if (_centroid[id1].delta > _centroid[id0].delta) {
                    std::swap (id0, id1);
                }

                for (uint j = 2; j < _opt.k; ++j)
                    if (_centroid[j].delta > _centroid[id1].delta) {
                        id1 = j;

                        if (_centroid[j].delta > _centroid[id0].delta) std::swap (id0, id1);
                    }

                for (std::vector <point_t>::iterator it = _point.begin ();
                        it != _point.end (); ++it) {
                    it->up_d += _centroid[it->id].delta;
                    it->lo_d -= _centroid[it->id == id0 ? id1 : id0].delta;
                }
            }
            uint& nf () {
                return _nf;
            }
            fl_t getObj () const { // const
                fl_t obj = 0;

                for (std::vector <point_t>::const_iterator it = _point.begin ();
                        it != _point.end (); ++it) {
                    obj += it->calc_dist (_centroid[it->id], _opt.dist);
                }

                return obj;
            }
            void run () {
                init ();
                uint moved = (_point.size ());

                for (uint i = 0; i <= _opt.iter; ++i) { // find neighbour center
                    if (moved) {
                        for (uint j = 0; j < _opt.k; ++j) // move center
                            _centroid[j].reset (_opt.dist);

                        update_bounds ();
                    }

                    if (i > 0) {
                        if (_opt.verbosity == 1) {
                            Rcout << i << std::setprecision (16) << "  obj = " << getObj () << " #moved = " << moved << "\n";
                        }

//          else
//          Rcout << ".";
                    }

                    if (! moved) break;

                    for (uint j = 0; j < _opt.k; ++j)
                        _centroid[j].set_closest (_centroid, _opt.dist);

                    moved = 0;

                    for (uint j = 0; j < _point.size (); ++j) { // for all points
                        point_t&   p  = _point[j];
                        const uint id0 = p.id;
                        const fl_t m   = std::max (_centroid[id0].next_d / 2, p.lo_d);

                        if (p.up_d > m) {
                            p.up_d = std::sqrt (p.calc_dist (_centroid[id0], _opt.dist));

                            if (p.up_d > m) {
                                p.set_closest (_centroid, _opt.dist);

                                if (p.id != id0) {
                                    ++moved;
                                    _centroid[id0].pop (p);
                                    _centroid[p.id].push (p);
                                }
                            }
                        }
                    }
                }

                if (_opt.verbosity == 1) {
                    if (moved > 0)
                        Rcout <<  "break";
                    else
                        Rcout <<  "done";

                    if (_opt.verbosity == 1)
                        Rcout << "; obj = " << getObj () << "\n";
                    else
                        Rcout << ".\n";
                }
            }
        private:
            const option _opt;
            std::vector <point_t>     _point;
            std::vector <centroid_t>  _centroid;
            std::vector <node_t>      _body;
            uint  _nf;
    };
// implementation of orthogonal k-means:
//   Y. Cui et al. Non-redundant multi-view clustering via orthogonalization (ICDM 2007)
    class orthogonal_kmeans {
        public:
            orthogonal_kmeans (const option &opt) : _opt (opt), _kms () {}
            ~orthogonal_kmeans () {
                for (std::vector <kmeans*>::iterator it = _kms.begin ();
                        it != _kms.end (); ++it)
                    delete *it;
            }
            void print (const uint i,
                        const std::vector <const char*>& label,
                        const std::vector <std::vector <uint> >& c2p) {
                for (uint j = 0; j < c2p.size (); ++j) {
                    if (_opt.m == 1)
                        Rcout << "c" << j;
                    else
                        Rcout << "c" << i << "_" << j;

                    for (std::vector <uint>::const_iterator it = c2p[j].begin ();
                            it != c2p[j].end (); ++it)
                        Rcout << label[*it];

                    Rcout << "\n";
                }
            }
            void print (const std::vector <const char*>& label,
                        const std::vector <std::vector <uint> >& p2c) {
                for (uint j = 0; j < label.size (); ++j) {
                    Rcout << label[j];

                    for (std::vector <uint>::const_iterator it = p2c[j].begin ();
                            it != p2c[j].end (); ++it)
                        Rcout << *it;

                    Rcout << "\n";
                }
            }
            void train_from_file (const char* train, const uint iter, const uint output = 0, const bool test_on_other_data = false, const bool instant = false) {
                std::vector <const char*> label;
                std::vector <std::vector <uint> > p2c; // point id to cluster id
                std::vector <std::vector <uint> > c2p (_opt.k); // cluster id to point id
                kmeans* km = new kmeans (_opt);
                FILE* fp = std::fopen (train, "r");

                if (! fp)
                    errx (1, "no such file: %s", train);

                char*  line = 0;
                unsigned long int read = 0;

                while (getLine (fp, line, read)) {
                    char* ex (line), *ex_end (line + read - 1);

                    while (ex != ex_end && ! isspace (*ex)) ++ex;

                    if (! test_on_other_data) {
                        char* copy = new char[ex - line + 1];
                        std::memcpy (copy, line, (ex - line));
                        copy[ex - line] = '\0';
                        label.push_back (copy);

                        if (output == 2) p2c.push_back (std::vector <uint> ());
                    }

                    while (isspace (*ex)) ++ex;

                    km->set_point (ex, ex_end, _opt.normalize);
                }

                if (km->point ().size () <= _opt.k)
                    errx (1, "# points (=%ld) <= k (=%d); done.",
                          km->point ().size (), _opt.k);

                _kms.push_back (km);
                std::fclose (fp);

                for (uint i = 1; i <= iter; ++i) {
                    Rcout << "iter=" << i << " k-means (k=" << _opt.k << "): ";

                    if (i >= 2) {
                        kmeans* km_ = _kms.back (); // last of mohikans
                        // project
                        std::vector <kmeans::point_t>& point_ = km_->point ();

                        for (std::vector <kmeans::point_t>::iterator it = point_.begin ();
                                it != point_.end (); ++it)
                            it->project (km_->centroid () [it->id]);

                        if (test_on_other_data || ! instant) {
                            km = new kmeans (_opt);
                            km_->delegate (km);
                            km_->compress ();
                            _kms.push_back (km);
                        } else {
                            km_->clear_centroid ();
                        }
                    }

                    km->run ();

                    if (! test_on_other_data) {
                        if (output == 1)
                            for (uint j = 0; j < _opt.k; ++j) c2p[j].clear ();

                        std::vector <kmeans::point_t> &point = km->point ();

                        for (uint j = 0; j < point.size (); ++j) {
                            if (output == 1) c2p[point[j].id].push_back (j);

                            if (output == 2) p2c[j].push_back (point[j].id);
                        }

                        if (output == 1) print (i, label, c2p);
                    }
                }

                if (! test_on_other_data)
                    if (output == 2) print (label, p2c);

                if (test_on_other_data || ! instant) {
                    _kms.back ()->compress ();
                    _kms.back ()->clear_point ();
                }

                for (std::vector <const char*>::iterator it = label.begin (); it != label.end (); ++it)
                    delete [] *it;
            }
            void save (const char* model) {
                FILE* fp = std::fopen (model, "w");
                std::fprintf (fp, "%d # m\n", _opt.m);
                std::fprintf (fp, "%d # k\n", _opt.k);
                std::fprintf (fp, "%d # number of features\n", _kms.back ()->nf ());

                for (uint i = 0; i < _kms.size (); ++i) {
                    const std::vector <kmeans::centroid_t>& centroid = _kms[i]->centroid ();

                    for (uint j = 0; j < centroid.size (); ++j) {
                        if (_opt.m == 1)
                            std::fprintf (fp, "c");
                        else
                            std::fprintf (fp, "c%d_", i);

                        centroid[j].print (fp, j);
                    }
                }

                std::fclose (fp);
            }
            void load (const char* model) {
                FILE* fp = std::fopen (model, "r");

                if (! fp)
                    errx (1, "no such file: %s", model);

                char*  line = 0;
                unsigned long int read = 0;

                if (! getLine (fp, line, read)) errx (1, "premature model: %s", model);

                _opt.m  = static_cast <uint> (std::strtol (line, NULL, 10));

                if (! getLine (fp, line, read)) errx (1, "premature model: %s", model);

                _opt.k  = static_cast <uint> (std::strtol (line, NULL, 10));

                if (! getLine (fp, line, read)) errx (1, "premature model: %s", model);

                const uint nf = static_cast <uint> (std::strtol (line, NULL, 10));
                std::vector <kmeans::node_t> body;

                for (uint i = 0; i < _opt.m; ++i) {
                    kmeans* km = new kmeans (_opt);
                    km->nf () = nf;

                    for (uint j = 0; j < _opt.k; ++j) {
                        if (! getLine (fp, line, read))
                            errx (1, "premature model: %s", model);

                        char* ex (line), *ex_end (line + read - 1);

                        while (ex != ex_end && ! isspace (*ex)) ++ex;

                        while (isspace (*ex)) ++ex;

                        kmeans::point_t p = kmeans::read_point (ex, ex_end, body, _opt.normalize);
                        km->push_centroid (p, true); // delegated
                    }

                    _kms.push_back (km);
                }

                std::fclose (fp);
            }
            void test_on_file (const char* test, const uint output = 0) {
                std::vector <kmeans::point_t>     point;
                std::vector <const char*>         label;
                std::vector <std::vector <uint> > p2c; // point id to cluster id
                std::vector <std::vector <uint> > c2p (_opt.k); // cluster id to point id
                std::vector <kmeans::node_t> body;
                FILE* fp = std::fopen (test, "r");

                if (! fp)
                    errx (1, "no such file: %s", test);

                char*  line = 0;
                unsigned long int read = 0;

                while (getLine (fp, line, read)) {
                    char* ex (line), *ex_end (line + read - 1);

                    while (ex != ex_end && ! isspace (*ex)) ++ex;

                    char* copy = new char[ex - line + 1];
                    std::memcpy (copy, line, (ex - line));
                    copy[ex - line] = '\0';
                    label.push_back (copy);

                    while (isspace (*ex)) ++ex;

                    point.push_back (kmeans::read_point (ex, ex_end, body, _opt.normalize));

                    if (output == 2) p2c.push_back (std::vector <uint> ());
                }

                std::fclose (fp);

                for (uint i = 0; i < _kms.size (); ++i) {
                    if (output == 1)
                        for (uint j = 0; j < _opt.k; ++j) c2p[j].clear ();

                    kmeans* km = _kms[i];
                    km->decompress ();

                    for (uint j = 0; j < point.size (); ++j) {
                        kmeans::point_t& p = point[j];
                        p.shrink (km->nf ());
                        p.set_closest (km->centroid (), _opt.dist);

                        if (output == 1) c2p[p.id].push_back (j);

                        if (output == 2) p2c[j].push_back (p.id);

                        p.project (km->centroid () [p.id]);
                    }

                    if (output == 1) print (i + 1, label, c2p);
                }

                if (output == 2) print (label, p2c);

                for (std::vector <kmeans::point_t>::iterator it = point.begin (); it != point.end (); ++it)
                    it->clear ();

                for (std::vector <const char*>::iterator it = label.begin (); it != label.end (); ++it)
                    delete [] *it;
            }
        private:
            const option          _opt;
            std::vector <kmeans*> _kms;
    };
}
