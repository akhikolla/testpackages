#v 0.2.6.00x

#v 0.2.6

Major changes:
- Added new `dodgr_centrality` function, and associated helper functions.
- Added new `dodgr_flows_si` function
- Added new parameter `norm_sums` to `dodgr_flows_aggregate`
- `merge_directed_flows` renamed to `merge_directed_graph`, with added option
  of specifying columns to merge.
- Added new `pairwise` parameter to `dodgr_distances`; see issue #127
- Added new function `dodgr_insert_vertex` to add new vertices to graph; see #40
- Removed "radix" heap option

Minor changes:
- switch off examples that caused previous CRAN failures
- fix bug in `dodgr_dists` when number of from points >> number of to points
- fix bug in `weight_streetnet.sc` that prevented `keep_cols` from working

# v0.2.5

- bug fixes from previous versions

# v0.2.4

Major changes:
- Remove benchmark vignette

Minor changes:
- bug fixes in `dodgr_paths`, thanks to @agila5


# v0.2.1

Major changes:
- Add `dodgr_isochrones`, `dodgr_isodistances`, and `dodgr_isoverts` functions
- Considerable speed-ups for `dodgr_flows_aggregate` and `dodgr_flows_disperse`

Minor changes:
- `dodgr_flows_disperse` allows `k` to be a vector, with different coefficients
  for each `from` point.
- Add "highway:pedestrian" to weighting profiles
- `weight_streetnet` for `sc` objects automatically adds component column
- bug fix in `weight_streetnet.sc(..., wt_profile = 1)`
- bug fix in `dodgr_full_cycles` for `SC` class objects

# v0.2.0

Major changes:
- Lots of intermediate processes now executed and cached as background
  processes (via `callr` package).
- new `dodgr_cache_off` function added to suppress primary caching, for cases
  where immediate usage is critical.
- `dodgr_contract_graph` returns the contracted graph only, instead of former
  version which return list of `graph` and `edge_map` (the `edge_map` is cached
  and re-loaded when needed for graph uncontraction).

Minor changes:
- 'turn_angle' parameter of `weight_streetnet` renamed to `turn_penalty`
- Test coverage now complete (100%)
* Update internal `hampi` data to remove factor columns
- Fix some bugs in max speed calculations for weight_streetnet
* Fix bug with polygonal bbox in dodgr_streetnet()

# v0.1.4
Major changes:
* New vignette on street networks and time-based routing
- `weight_streetnet` function now returns edge times for all Open Street Map
  networks extracted with the `osmdata` package.
- weight_streetnet now accepts `SC` format data from `osmdata::osmdata_sc()`
- New `dodgr_times` function to calculate journey times, including differential
  speeds and penalties for intersections and traffic lights.
- `dodgr::weighting_profiles` data changed from single `data.frame` to list with
  additional parameters determining speeds and time penalties for `dodgr_times`
  function; former `data.frame` is now
  `dodgr::weighting_profiles$weighting_profiles`.
- New function `write_dodgr_wt_profile` writes full profile to local `.json`
  file for editing and subsequent use via
  `weight_streetnet(wt_profile_file=<local_file_name.json>)`.
- `dodgr_dists()`, `dodgr_paths()`, and `dodgr_flows()` can no longer be used
  to automatically download street networks, thus former parameters
  `wt_profile` and `expand` have been removed; networks must be explicitly
  downloaded with `dodgr_streetnet()`.

Minor changes:
- Bug fix with dodgr_to_igraph to create proper *weighted* igraph
- Add "footway" to weighting_profiles

# v0.1.3
Major changes:
- New functions `dodgr_fundamental_cycles` and `dodgr_full_cycles`
- New function `dodgr_sflines_to_poly` to convert `sf` collections of
  `LINESTRING` object into corresponding enclosed `POLYGON` objects.
- New function `dodgr_to_sf` creates full `sf` objects, extending `dodgr_to_sfc`
- New function `igraph_to_dodgr` converts `igraph` objects into `dodgr` format
- New function `dodgr_uncontract_graph` to convert from contracted back into
  original, uncontracted from, including any additional data appended on to
  contracted graph.

Minor changes:
- Bug fix with vignette caused by updates to `tinytex` rendering of `svg`
- Bug fix for `dodgr_dists (heap = "set")` with integer distances


# v0.1.2
Major changes:
- New function `dodgr_to_igraph`
- `weight_streetnet` is now a method, with implementations for objects of
  classes `.sf` and `.sc`.
- New function `weight_railway` to weight a network for railway routing.
- `dodgr_dists` implements Dijkstra paths with std::set sorting through new
  option `dodgr_dists(..., heap = "set")` (It's slower than others, but good for
  sake of completeness).

Minor changes:
- Various modifications that should result in notable speed gains
- `dodgr_streetnet` now accepts polygonal `bbox` argument, and uses
  `osmdata::trim_osmdata` to trim resultant network to within that polygon
  (issue #50).
- Extended examples for `weight_streetnet` and dodgr_flows_aggregate` to include
  a non-OSM example from `stplanr::routes_fast` (issue #45).


# v0.1.1

Major changes:
- Crucial fix of previous typo that made all `dodgr_dist` calculations wrong
  (Earth's radius is 6371, not 3671!) - thanks to @chrijo
- `weight_streetnet` function now accepts `data.frame` objects defining
  `wt_profile`, enabling modification and direct re-submission of
  `dodgr::weighting_profiles`
- `weighting_profiles$value` modified to 0-1 scores rather than previous
  percentage values.
- `weight_streetnet` flags any highway types not present in nominated or
  submitted weighting profile.
- `dodgr_paths` now has additional `pairwise` parameter to enable paths only
  between matched pairs of `from` and `to` points (so returning `n` paths rather
  than `n^2`), thanks to @mem48.
- `dodgr_to_sf` deprecated to `dodgr_to_sfc` (#43)

Minor changes:
- Added Malcolm Morgan (@mem48; bug-finder extraordinare) as contributor 
- Bug fix with `dodgr_paths` and simple `data.frame`s, thanks to James Smith.
- Bug fix of former improper handling of one-way edges, thanks to @chrijo.
- `match_pts_to_graph` has additional `connected` parameter to allow points to
  be matched only to largest connected component.

# v0.1.0

Major changes:
- New function `dodgr_flowmap` plots maps of flows. Currently only writes .png
  files, because large networks can not be effectively plotted on graphic
  devices.
- `dodgr_flows` has option to routes flows from a set of source origins to all
  points in a network, attenuated by distance from those origins.
- `dodgr_to_sf` converts a spatially-explicit `dodgr` graph into Simple Features
  (`sf`) format.

Minor changes:
- `match_pts_to_graph` now accepts Simple Features (sf) collections of
  `sfc_POINT` objects to be matched.

# v0.0.3

Tidy C++ code that flagged errors on CRAN solaris machine. Nothing else.

# v0.0.2

Major changes:
- New function, `dodgr_paths`, for returning explicit shortest path routes.
- New function, `dodgr_flows`, for aggregting flows across a network from
  multiple origin and destination points.
- New function, `merge_directed_flows`, to reduce aggregated directional flows
  to non-directional equivalent values useful for visualisation.

Minor changes:
- `weight_streetnet` now accepts arbitrary `sf`-formatted networks via
  specification of custom weighting profiles, along with highway type and ID
  columns in data.frame.
- Distance matrices from `dodgr_dists` inherit the names of routing points
  (`from` and `to` parameters).

# v0.0.1

Initial CRAN release.
