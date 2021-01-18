// grid.h
#ifndef GRID_H_
#define GRID_H_

#include <cstddef>
#include <vector>
#include <tuple>
#include <numeric>    // accumulate
#include <functional> // multiplies
#include <iterator>
#include <type_traits>

namespace grid {
////////////////////////////////////////////////////////////////
// interface
////////////////////////////////////////////////////////////////

////////////////////////////////
// forward declarations
//
// GridIter == DataIter + CoordIter in lockstep
// - DataIter knows nothing about the coord
// - CoordIter knows noting about the data
// - GridIter syncs the two
template<typename T, typename...CoordSpec> class GridIter;
template<typename T, typename...CoordSpec> class GridConstIter;
template<typename T> using DataIter = typename std::vector<T>::iterator;
template<typename T> using DataConstIter = typename std::vector<T>::const_iterator;
template<typename...CoordSpec> class CoordIter;

////////////////////////////////
// Main data structure
// - store datum (an array of T * blocksize) at each point on a grid
// - storage is a simple vector of T (size = blocksize * ngrid_points)
// - coordinate axes (vectors) are stored in a tuple (coords)
template<typename T, typename...CoordSpec>
class Grid {
public:
  // ctor from blocksize and a list of coordinate vectors
  Grid(std::size_t blocksize, std::vector<CoordSpec> const&...);
  Grid(std::size_t blocksize, std::tuple<std::vector<CoordSpec> const&...>&&);
  // ctor with pre-allocated storage
  Grid(std::size_t blocksize, std::vector<T>&&, std::vector<CoordSpec> const&...);

  Grid()                                     = delete;
  Grid(Grid const&)                          = delete;
  Grid& operator=(Grid const&)               = delete;
  Grid(Grid&&)                               = delete;
  Grid& operator=(Grid&&)                    = delete;
  ~Grid()                                    = default;

  inline std::size_t                           size()    const; // in the unit of T
  inline std::vector<std::size_t>       const& dimspec() const; // npoints on each coord

  // Iterator interface
  using  const_iterator                      = GridConstIter<T, CoordSpec...>;
  inline const_iterator                        cbegin() const;
  inline const_iterator                        begin()  const;
  inline const_iterator                        cend()   const;
  inline const_iterator                        end()    const;
  using  iterator                            = GridIter<T, CoordSpec...>;
  inline iterator                              begin();
  inline iterator                              end();
private:
  std::vector<std::size_t>              const  dimspec_;
  std::size_t                           const  size_;
public:
  // we copy coords as the lifetime of a grid can be longer than the
  // coords from which it is constructed
  std::tuple<std::vector<CoordSpec> const...>
                                        const  coords;
  std::size_t                           const  blocksize;
  std::vector<T>                               data;
};

////////////////////////////////
// Const iterator to traverse over the grids
template<typename T, typename...CoordSpec>
class GridConstIter {
public:
  // ctor from a grid instance
  explicit GridConstIter(Grid<T, CoordSpec...> const&);
  // cast from non-const iterator
  GridConstIter(GridIter<T, CoordSpec...> const&);

  GridConstIter()                                = delete;
  GridConstIter(GridConstIter const&)            = delete;
  GridConstIter& operator=(GridConstIter const&) = delete;
  GridConstIter(GridConstIter &&)                = default;
  GridConstIter& operator=(GridConstIter &&)     = default;
  ~GridConstIter()                               = default;

  friend class Grid<T, CoordSpec...>;

  // data on the current gridpoint
  inline DataConstIter<T>                          cbegin() const;
  inline DataConstIter<T>                          cend()   const;
  inline DataConstIter<T>                          begin()  const;
  inline DataConstIter<T>                          end()    const;

  // coord label of the current gridpoint (valid until next inc)
  // note: next() instead of ++? or operator* for coord()?
  inline std::tuple<CoordSpec...>           const& coord();

  inline GridConstIter&                            operator++(); // inc to next gridpoint
  // GridConstIter&                                   operator++(int) { throw std::logic_error("GridConstIter: use prefix ++."); }
  inline void                                      reset();

  // friend hack for template function instanciations..
  friend bool operator==(GridConstIter<T, CoordSpec...> const& i,
                         GridConstIter<T, CoordSpec...> const& j) {
    return i.equal(j);
  }
  friend bool operator!=(GridConstIter<T, CoordSpec...> const& i,
                         GridConstIter<T, CoordSpec...> const& j) {
    return !i.equal(j);
  }
  // user is responsible for ensuring that i and j point to a same grid
  friend bool operator< (GridConstIter<T, CoordSpec...> const& i,
                         GridConstIter<T, CoordSpec...> const& j) {
    return i.less(j);
  }
private:
  inline bool                                      equal(GridConstIter const&) const;
  inline bool                                      less (GridConstIter const&) const;

  std::size_t                               const  blocksize_;
  DataConstIter<T>                          const  dbeg_;
  DataConstIter<T>                          const  dend_;
  DataConstIter<T>                                 diter_;
  CoordIter<CoordSpec...>                          citer_;
};

////////////////////////////////
// Read-write iterator over grid
// (Could have merged GridConstIter and GridIter, but that leads to very ugly codes..)
template<typename T, typename...CoordSpec>
class GridIter {
public:
  // ctor from a grid instance
  explicit GridIter(Grid<T, CoordSpec...>&);
  GridIter()                                 = delete;
  GridIter(GridIter const&)                  = delete;
  GridIter& operator=(GridIter const&)       = delete;
  GridIter(GridIter &&)                      = default;
  GridIter& operator=(GridIter &&)           = default;
  ~GridIter()                                = default;

  friend class Grid<T, CoordSpec...>;
  friend class GridConstIter<T, CoordSpec...>;
  // cast to const_iterator
  // operator GridConstIter<T, CoordSpec...>();

  // data on the current gridpoint
  inline DataConstIter<T>                      cbegin() const;
  inline DataConstIter<T>                      cend()   const;
  inline DataConstIter<T>                      begin()  const;
  inline DataConstIter<T>                      end()    const;
  inline DataIter<T>                           begin();
  inline DataIter<T>                           end();

  // coord label of the current gridpoint (valid until next inc)
  inline std::tuple<CoordSpec...>       const& coord();

  inline GridIter&                             operator++(); // inc to next gridpoint
  // GridIter&                                    operator++(int) { throw std::logic_error("GridIter: use prefix ++."); }
  inline void                                  reset();

  // friend hack for template function instanciations
  friend bool operator==(GridIter<T, CoordSpec...> const& i,
                         GridIter<T, CoordSpec...> const& j) {
    return i.equal(j);
  }
  friend bool operator==(GridIter<T, CoordSpec...>      const& i,
                         GridConstIter<T, CoordSpec...> const& j) {
    return j == i; // note the order
  }
  friend bool operator!=(GridIter<T, CoordSpec...> const& i,
                         GridIter<T, CoordSpec...> const& j) {
    return !j.equal(i);
  }
  friend bool operator!=(GridIter<T, CoordSpec...>      const& i,
                         GridConstIter<T, CoordSpec...> const& j) {
    return !(j == i); // note the order
  }
  // user is responsible for ensuring that i and j point to a same grid
  friend bool operator< (GridIter<T, CoordSpec...> const& i,
                         GridIter<T, CoordSpec...> const& j) {
    return i.less(j);
  }
  friend bool operator< (GridIter<T, CoordSpec...>      const& i,
                         GridConstIter<T, CoordSpec...> const& j) {
    // i < j <=> !(j <= i) <=> (j != i) && !(j < i)
    return (j != i) && !(j < i);
  }
private:
  inline bool                                  equal(GridIter const&) const;
  inline bool                                  less (GridIter const&) const;

  std::size_t                            const blocksize_;
  DataIter<T>                            const dbeg_;
  DataIter<T>                            const dend_;
  DataIter<T>                                  diter_;
  CoordIter<CoordSpec...>                      citer_;
};

////////////////////////////////////////////////////////////////
// implementation
////////////////////////////////////////////////////////////////

////////////////////////////////
// support routines for the variadic ctor
template<std::size_t...>
struct Index {};

template<std::size_t End, std::size_t...i>
struct MakeIndex : public MakeIndex<End-1, End-1, i...> {};

template<std::size_t...i>
struct MakeIndex<0, i...> {
  using type = Index<i...>;
};

template<typename...Ts>
using make_index = typename MakeIndex<sizeof...(Ts)>::type;

////////////////////////////////
template<typename...CoordSpec, std::size_t...i>
std::vector<std::size_t>
make_dimspec(Index<i...>, std::tuple<std::vector<CoordSpec> const&...>&& coords) {
  return { std::get<i>(coords).size()... };
}

////////////////////////////////
// makes no effort to check if there are no duplicate values in each coord vector
template<typename T, typename...CoordSpec>
Grid<T, CoordSpec...>::Grid(std::size_t blocksize,
                            std::vector<CoordSpec> const&...coords)
  : dimspec_(make_dimspec(make_index<CoordSpec...>(), std::tie(coords...)))
  , size_(std::accumulate(dimspec_.begin(), dimspec_.end(),
                          blocksize, std::multiplies<std::size_t>()))
  , coords(coords...)
  , blocksize(blocksize)
  , data(std::vector<T>(size_))
{}

template<typename T, typename...CoordSpec>
Grid<T, CoordSpec...>::Grid(std::size_t blocksize,
                            std::tuple<std::vector<CoordSpec> const&...>&& coords)
  : dimspec_(make_dimspec(make_index<CoordSpec...>(), std::move(coords)))
  , size_(std::accumulate(dimspec_.begin(), dimspec_.end(),
                          blocksize, std::multiplies<std::size_t>()))
  , coords(coords)
  , blocksize(blocksize)
  , data(std::vector<T>(size_))
{}

////////////////////////////////
template<typename T, typename...CoordSpec>
inline std::size_t Grid<T, CoordSpec...>::size() const {
  return size_;
}

template<typename T, typename...CoordSpec>
inline std::vector<std::size_t> const& Grid<T, CoordSpec...>::dimspec() const {
  return dimspec_;
}

template<typename T, typename...CoordSpec>
inline typename Grid<T, CoordSpec...>::iterator
Grid<T, CoordSpec...>::begin() {
  return GridIter<T, CoordSpec...> { *this };
}

template<typename T, typename...CoordSpec>
inline typename Grid<T, CoordSpec...>::const_iterator
Grid<T, CoordSpec...>::begin() const {
  return GridConstIter<T, CoordSpec...> { *this };
}

template<typename T, typename...CoordSpec>
inline typename Grid<T, CoordSpec...>::const_iterator
Grid<T, CoordSpec...>::cbegin() const {
  return GridConstIter<T, CoordSpec...> { *this };
}

template<typename T, typename...CoordSpec>
inline typename Grid<T, CoordSpec...>::iterator
Grid<T, CoordSpec...>::end() {
  GridIter<T, CoordSpec...> it { *this };
  it.diter_ = it.dend_; // only DataIter is used for equality check
  return it;
}

template<typename T, typename...CoordSpec>
inline typename Grid<T, CoordSpec...>::const_iterator
Grid<T, CoordSpec...>::end() const {
  GridConstIter<T, CoordSpec...> it { *this };
  it.diter_ = it.dend_; // only DataIter is used for equality check
  return it;
}

template<typename T, typename...CoordSpec>
inline typename Grid<T, CoordSpec...>::const_iterator
Grid<T, CoordSpec...>::cend() const {
  GridConstIter<T, CoordSpec...> it { *this };
  it.diter_ = it.dend_; // only DataIter is used for equality check
  return it;
}

////////////////////////////////////////////////////////////////

template<typename...CoordSpec, std::size_t...i>
std::tuple<typename std::vector<CoordSpec>::const_iterator...>
make_begs(Index<i...>, std::tuple<std::vector<CoordSpec> const&...> const& cs) {
  return std::make_tuple(std::get<i>(cs).begin()...);
  // an llvm extension allows the following literal notation
  // return { std::get<i>(cs).begin()... };
}

template<typename...CoordSpec, std::size_t...i>
std::tuple<typename std::vector<CoordSpec>::const_iterator...>
make_ends(Index<i...>, std::tuple<std::vector<CoordSpec> const&...> const& cs) {
  return std::make_tuple(std::get<i>(cs).end()...);
}

////////////////////////////////
// Read-only forward iterator over coordinate labels
// - intended to be used in sync with DataIter
template<typename...CoordSpec>
class CoordIter {
public:
  using iter_tuple =
    typename std::tuple<typename std::vector<CoordSpec>::const_iterator...>;
  using val_tuple = typename std::tuple<CoordSpec...>;

  explicit CoordIter(std::tuple<std::vector<CoordSpec> const&...> const& cs)
    : begs_{make_begs(make_index<CoordSpec...>(), cs)}
    , ends_{make_ends(make_index<CoordSpec...>(), cs)}
    , iters_{begs_}
    {}
  CoordIter()                            = delete;
  // constructing GridConstIter from grid.begin() needs GridIter to GridConstIter
  // conversion, which in turn uses copy construction of CoordIter.
  CoordIter(CoordIter const&)            = default;
  CoordIter& operator=(CoordIter const&) = delete;
  CoordIter(CoordIter&&)                 = default;
  CoordIter& operator=(CoordIter&&)      = default;
  ~CoordIter()                           = default;

  CoordIter& operator++() {
    inc<sizeof...(CoordSpec)-1>();
    return *this;
  }
  // CoordIter& operator++(int) { throw std::logic_error("CoordIter: use prefix ++."); }

  val_tuple const& operator*() const {
    setvals_<sizeof...(CoordSpec)-1>();
    return vals_;
  }

  inline bool isend() { return std::get<0>(iters_) == std::get<0>(ends_); }
  inline void reset() { setiters_<sizeof...(CoordSpec)-1>(begs_); }

private:
  template <std::size_t i>
  inline typename std::enable_if<i!=0, void>::type inc() {
    // inv: iters_<0> points at the end
    //      iff the iterator has gone through all the points.
    // We make no promises on the states of iters_<1>, iters_<2>, ..
    // when iters_<0> is pointing at the end.
    ++std::get<i>(iters_); // try inc the ith slot
    // if iters_<i> is now pointing at the end, reset and carry
    if (std::get<i>(iters_) == std::get<i>(ends_)) {
      std::get<i>(iters_) = std::get<i>(begs_);
      inc<i-1>();
    }
  }
  template <std::size_t i>
  inline typename std::enable_if<i==0, void>::type inc() {
    // ++ on end (one past the last) is undefined
    // we choose iters_<0> to stay at the end (++end == end)
    if (std::get<i>(iters_) != std::get<i>(ends_))
      ++std::get<i>(iters_);
  }

  template <std::size_t i>
  inline typename std::enable_if<i!=0, void>::type setiters_(iter_tuple const& to) {
    std::get<i>(iters_) = std::get<i>(to);
    setiters_<i-1>(to);
  }
  template <std::size_t i>
  inline typename std::enable_if<i==0, void>::type setiters_(iter_tuple const& to) {
    std::get<i>(iters_) = std::get<i>(to);
  }

  // vals_ is a mutable cache, hence we define setvals_() as a const function
  template <std::size_t i>
  inline typename std::enable_if<i!=0, void>::type setvals_() const {
    std::get<i>(vals_) = *std::get<i>(iters_);
    setvals_<i-1>();
  }
  template <std::size_t i>
  inline typename std::enable_if<i==0, void>::type setvals_() const {
    std::get<i>(vals_) = *std::get<i>(iters_);
  }

  iter_tuple const begs_; // cache begin() and end() as Grid is unmovable
  iter_tuple const ends_;
  iter_tuple iters_;
  val_tuple mutable vals_; // cache *iters_ (updated when derefed)
};

////////////////////////////////////////////////////////////////
// T can come with const qualifier
template<typename T, typename...CoordSpec>
GridConstIter<T, CoordSpec...>::GridConstIter(Grid<T, CoordSpec...> const& grid)
  : blocksize_{grid.blocksize}
  , dbeg_{grid.data.begin()}
  , dend_{grid.data.end()}
  , diter_{grid.data.begin()}
  , citer_{grid.coords}
 {}

template<typename T, typename...CoordSpec>
GridConstIter<T, CoordSpec...>::GridConstIter(GridIter<T, CoordSpec...> const& i)
  : blocksize_{i.blocksize_}
  , dbeg_{i.dbeg_}
  , dend_{i.dend_}
  , diter_{i.diter_}
  , citer_{i.citer_}
 {}

template<typename T, typename...CoordSpec>
inline GridConstIter<T, CoordSpec...>&
GridConstIter<T, CoordSpec...>::operator++() {
  diter_ += blocksize_;
  ++citer_;
  return *this;
}

template<typename T, typename...CoordSpec>
inline void GridConstIter<T, CoordSpec...>::reset() {
  diter_ = dbeg_;
  citer_.reset();
}

template<typename T, typename...CoordSpec>
inline DataConstIter<T>
GridConstIter<T, CoordSpec...>::cbegin() const { return diter_; }

template<typename T, typename...CoordSpec>
inline DataConstIter<T>
GridConstIter<T, CoordSpec...>::begin() const { return diter_; }

template<typename T, typename...CoordSpec>
inline DataConstIter<T>
GridConstIter<T, CoordSpec...>::cend() const { return diter_ + blocksize_; }

template<typename T, typename...CoordSpec>
inline DataConstIter<T>
GridConstIter<T, CoordSpec...>::end() const { return diter_ + blocksize_; }

template<typename T, typename...CoordSpec>
inline std::tuple<CoordSpec...> const& GridConstIter<T, CoordSpec...>::coord() {
  return *citer_;
}

template<typename T, typename...CoordSpec>
inline bool GridConstIter<T, CoordSpec...>::equal(GridConstIter const& that) const {
  // return cur() == that.cur();
  return diter_ == that.diter_;
}

template<typename T, typename...CoordSpec>
inline bool GridConstIter<T, CoordSpec...>::less(GridConstIter const& that) const {
  // return cur() < that.cur();
  return diter_ < that.diter_;
}

////////////////////////////////
template<typename T, typename...CoordSpec>
GridIter<T, CoordSpec...>::GridIter(Grid<T, CoordSpec...>& grid)
  : blocksize_(grid.blocksize)
  , dbeg_(grid.data.begin())
  , dend_(grid.data.end())
  , diter_(grid.data.begin())
  , citer_(grid.coords)
{}

template<typename T, typename...CoordSpec>
inline GridIter<T, CoordSpec...>&
GridIter<T, CoordSpec...>::operator++() {
  diter_ += blocksize_;
  ++citer_;
  return *this;
}

template<typename T, typename...CoordSpec>
inline void GridIter<T, CoordSpec...>::reset() {
  diter_ = dbeg_;
  citer_.reset();
}

template<typename T, typename...CoordSpec>
inline DataConstIter<T> GridIter<T, CoordSpec...>::cbegin() const { return diter_; }

template<typename T, typename...CoordSpec>
inline DataConstIter<T> GridIter<T, CoordSpec...>::begin() const { return diter_; }

template<typename T, typename...CoordSpec>
inline DataIter<T> GridIter<T, CoordSpec...>::begin() { return diter_; }

template<typename T, typename...CoordSpec>
inline DataConstIter<T> GridIter<T, CoordSpec...>::cend() const { return diter_ + blocksize_; }

template<typename T, typename...CoordSpec>
inline DataConstIter<T> GridIter<T, CoordSpec...>::end() const { return diter_ + blocksize_; }

template<typename T, typename...CoordSpec>
inline DataIter<T> GridIter<T, CoordSpec...>::end() { return diter_ + blocksize_; }

template<typename T, typename...CoordSpec>
inline std::tuple<CoordSpec...> const & GridIter<T, CoordSpec...>::coord() {
  return *citer_;
}

template<typename T, typename...CoordSpec>
inline bool GridIter<T, CoordSpec...>::equal(GridIter const& that) const {
  // return cur() == that.cur();
  return diter_ == that.diter_;
}

template<typename T, typename...CoordSpec>
inline bool GridIter<T, CoordSpec...>::less(GridIter const& that) const {
  // return cur() < that.cur();
  return diter_ < that.diter_;
}

////////////////////////////////////////////////////////////////
} // namespace grid

#endif // GRID_H_
