#include "csa.h"

//CSA_Multi_Return csa_matrix::csa_matrix_scan (const CSA_Parameters &csa_pars,
CSA_Multi_Return csa_matrix_scan (const CSA_Parameters &csa_pars,
        const std::unordered_set <size_t> &start_stations_set,
        std::unordered_set <size_t> &end_stations_set,
        const CSA_Inputs &csa_in,
        CSA_Outputs &csa_out)
{
    CSA_Return csa_ret;
    csa_ret.earliest_time = INFINITE_INT;
    csa_ret.end_station = INFINITE_INT;

    std::vector <bool> is_connected (csa_pars.ntrips, false);

    // trip connections:
    for (size_t i = 0; i < csa_pars.timetable_size; i++)
    {
        if (csa_in.departure_time [i] < csa_pars.start_time)
            continue; // # nocov - these lines already removed in R fn.

        // add all departures from start_stations_set:
        if (start_stations_set.find (csa_in.departure_station [i]) !=
                start_stations_set.end () &&
                csa_in.arrival_time [i] < csa_out.earliest_connection [csa_in.arrival_station [i] ])
        {
            is_connected [csa_in.trip_id [i] ] = true;
            csa::fill_one_csa_out (csa_out, csa_in, csa_in.arrival_station [i], i);
        }

        // main connection scan:
        if (((csa_out.earliest_connection [csa_in.departure_station [i] ] <= csa_in.departure_time [i]) &&
                    csa_out.n_transfers [csa_in.departure_station [i] ] < csa_pars.max_transfers) ||
                is_connected [csa_in.trip_id [i]])
        {
            if (csa_in.arrival_time [i] < csa_out.earliest_connection [csa_in.arrival_station [i] ])
            {
                csa::fill_one_csa_out (csa_out, csa_in, csa_in.arrival_station [i], i);

                csa_out.n_transfers [csa_in.arrival_station [i] ] =
                    csa_out.n_transfers [csa_in.departure_station [i] ];
            }
            csa::check_end_stations (end_stations_set, csa_in.arrival_station [i],
                    csa_in.arrival_time [i], csa_ret);

            if (csa_in.transfer_map.find (csa_in.arrival_station [i]) != csa_in.transfer_map.end ())
            {
                for (auto t: csa_in.transfer_map.at (csa_in.arrival_station [i]))
                {
                    size_t trans_dest = t.first;
                    int ttime = csa_in.arrival_time [i] + t.second;
                    if (ttime < csa_out.earliest_connection [trans_dest] &&
                            csa_out.n_transfers [trans_dest] <= csa_pars.max_transfers)
                    {
                        // modified version of fill_one_csa_out:
                        csa_out.earliest_connection [trans_dest] = ttime;
                        csa_out.prev_stn [trans_dest] = csa_in.arrival_station [i];
                        csa_out.prev_time [trans_dest] = csa_in.arrival_time [i];
                        csa_out.n_transfers [trans_dest]++;

                        csa::check_end_stations (end_stations_set,
                                trans_dest, ttime, csa_ret);

                    }
                }
            }
            is_connected [csa_in.trip_id [i]] = true;
        }
        if (end_stations_set.size () == 0)
            break;
    }
    //return csa_ret;
    CSA_Multi_Return csa_multi_ret;

    return csa_multi_ret;
}

/* *****  ALGORITHM  *****
 *
 * Scan entire timetable in one go, keeping paired list of all start_stns and
 * end_stns, each having `n` entries. The main inputs are:
 * 1. The timetable of (departure_stn, arrival_stn, departure_time,
 * arrival_time, trip_id)
 * 2. The transfers table of (from_stn, to_stn, min_transfer_time).
 *
 * Main data are:
 * 1. A square n-by-n matrix, `scan_mat`, of (from, to) times to INFINITE_INT
 * 2. An unordered_map from `trip_id` values to all start_stns connecting with
 * those, `trip_id_map`.
 *
 * The main `scan_mat` records the earliest arrival time between each station
 * pair. This is successively filled by keeping track of all previous stations
 * connecting with current one. At each step, the station list of `trip_id_map`
 * is extracted, and for each of those, `s`, the value of `scan_mat (s,
 * arrival_stn)` is updated if listed arrival time is less than current value.
 *
 * In pseudo-code format:
 *
 * 1. trip_stns = trip_id_map.at (trip_id [i])
 * 2. for s in trip_stns:
 *      if (scan_mat (s, arrival_stn [i]) > arrival_time [i]):
 *          scan_mat (s, arrival_stn [i]) = arrival_time [i]
 * 3. transfers = transfer_map.at (arrival_stn [i])
 * 4. trip_id_map.at (trip_id [i]) += transfers
 *
 * The whole thing scans the entire timetable from the specified `start_time` to
 * the very end.  Then just return `scan_mat` and job done.
 */

//' rcpp_csa_matrix
//'
//' Not really CSA, rather a full-timetable scan for all connections from
//' everywhere to everywhere which returns only start and end times for each
//' pair of (start, end) points.
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::DataFrame rcpp_csa_matrix (Rcpp::DataFrame timetable,
        Rcpp::DataFrame transfers,
        const size_t nstations,
        const size_t ntrips,
        const std::vector <size_t> start_stations,
        const std::vector <size_t> end_stations,
        const int start_time,
        const int max_transfers)
{
    CSA_Parameters csa_pars;
    csa::fill_csa_pars (csa_pars, max_transfers, start_time,
            static_cast <size_t> (timetable.nrow ()), ntrips, nstations);

    std::unordered_set <size_t> start_stations_set, end_stations_set;
    csa::make_station_sets (start_stations, end_stations,
            start_stations_set, end_stations_set);

    CSA_Inputs csa_in;
    csa::make_transfer_map (csa_in.transfer_map, transfers);

    // The csa_out vectors use nstations + 1 because it's 1-indexed throughout,
    // and the first element is ignored.
    CSA_Outputs csa_out;
    csa_out.earliest_connection.resize (csa_pars.nstations + 1, INFINITE_INT);
    csa_out.n_transfers.resize (csa_pars.nstations + 1, 0);
    csa_out.prev_time.resize (csa_pars.nstations + 1, INFINITE_INT);
    csa_out.prev_stn.resize (csa_pars.nstations + 1, INFINITE_INT);
    csa_out.current_trip.resize (csa_pars.nstations + 1, INFINITE_INT);

    csa::get_earliest_connection (start_stations, csa_pars.start_time,
            csa_in.transfer_map, csa_out.earliest_connection);

    csa::csa_in_from_df (timetable, csa_in);

    //CSA_Multi_Return csa_ret = csa_matrix::csa_matrix_scan (csa_pars, start_stations_set,
    //        end_stations_set, csa_in, csa_out);
    CSA_Multi_Return csa_ret = csa_matrix_scan (csa_pars, start_stations_set,
            end_stations_set, csa_in, csa_out);

    Rcpp::DataFrame res;

    return res;
}

