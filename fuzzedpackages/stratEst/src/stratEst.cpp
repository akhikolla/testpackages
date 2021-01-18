//#define ARMA_NO_DEBUG
#define ARMA_USE_CXX14
#include <RcppArmadillo.h>
using namespace Rcpp;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// stratEst_EM
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

arma::field<arma::mat> stratEst_EM(arma::cube& output_cube, arma::cube& sum_outputs_cube, arma::vec& strat_id, arma::vec shares, arma::vec responses, arma::vec trembles, arma::mat share_mat, arma::mat response_mat, arma::mat tremble_mat, arma::mat indices_shares, arma::mat indices_responses, arma::mat indices_trembles, arma::mat& responses_to_sum, std::string& response, int eval_pre , double tol_eval, int max_eval, arma::vec& sample_of_ids_shares, arma::vec& num_ids_sample_shares, arma::vec& sample_of_ids_responses, arma::vec& num_ids_sample_responses, arma::vec& sample_of_ids_trembles, arma::vec& num_ids_sample_trembles ) {

  arma::field<arma::mat> F(26,1);
  int num_ids = output_cube.n_slices;
  int num_rows_response_mat = response_mat.n_rows;
  int num_cols_response_mat = response_mat.n_cols;
  int num_rows_tremble_mat = tremble_mat.n_rows;
  int num_cols_tremble_mat = tremble_mat.n_cols;
  arma::uvec shares_to_est = find( indices_shares != 0 );
  arma::uvec responses_to_est = find( indices_responses != 0 );
  int num_shares_to_est = indices_shares.max();
  int num_samples_shares = sample_of_ids_shares.max();
  int num_responses_to_est = responses.n_elem;
  int num_trembles_to_est = trembles.n_elem;
  int k = share_mat.n_rows;
  int num_samples = share_mat.n_cols;
  int num_outputs = num_cols_response_mat;
  int num_samples_responses = sample_of_ids_responses.max();
  int num_samples_trembles = sample_of_ids_trembles.max();
  arma::mat remaining_shares_mat( k , num_samples , arma::fill::zeros );
  arma::vec num_fixed_shares_col( num_samples , arma::fill::zeros );
  int free_shares = 0;
  for (int i = 0; i < num_samples; i++) {
    arma::vec shares_col = share_mat( arma::span::all , i );
    int num_shares_col = shares_col.n_elem;
    arma::vec indices_shares_col = indices_shares( arma::span::all , i );
    arma::uvec fixed_shares_col = find( indices_shares_col == 0 );
    double remaining_shares = 1 - accu( shares_col( fixed_shares_col ) );
    remaining_shares_mat( arma::span::all , i ).fill( remaining_shares );
    int num_fixed_shares = fixed_shares_col.n_elem;
    num_fixed_shares_col(i) = num_fixed_shares;
    if( i == 0 ){
      if( num_fixed_shares != num_shares_col ){
        free_shares = free_shares + indices_shares_col.n_elem - fixed_shares_col.n_elem - 1;
      }
    }
    else{
      if( num_samples_shares > 1 ){
        free_shares = free_shares + indices_shares_col.n_elem - fixed_shares_col.n_elem - 1;
      }
    }
  }
  arma::umat indices_responses_to_sum = find( responses_to_sum == 1 );
  arma::umat indices_non_fixed_responses = find( indices_responses == 0 && responses_to_sum == 0 );
  int free_responses = 0;
  if( response == "pure" ){
    free_responses = num_responses_to_est/num_cols_response_mat;
  }else{
    free_responses = num_responses_to_est - num_responses_to_est/num_cols_response_mat;
  }
  int free_params = free_shares + free_responses + num_trembles_to_est;
  double eps = 1;
  double eps_now = 1;
  arma::vec ll_val( 1 , 1 , arma::fill::zeros );
  arma::vec entropy_k( 1 , arma::fill::ones);
  arma::vec crit_vals( 3 , arma::fill::ones);
  arma::vec new_ll_val( 1 , 1 , arma::fill::zeros );
  arma::vec new_entropy_k( 1 , arma::fill::ones);
  arma::mat i_shares_mat( num_ids , k , arma::fill::zeros );
  arma::cube i_shares_cube( num_rows_response_mat , num_cols_response_mat , num_ids , arma::fill::zeros );
  arma::cube weigthed_sums_cube( num_rows_response_mat , num_cols_response_mat , num_ids , arma::fill::zeros );
  //
  arma::cube weigthed_output_cube = output_cube;
  //
  arma::vec sample_of_ids_max = sample_of_ids_responses;
  if( max(sample_of_ids_max) < max(sample_of_ids_trembles ) ){
    sample_of_ids_max = sample_of_ids_trembles;
  }
  if( ( num_samples_responses == 1 || num_responses_to_est == 0 ) & ( num_samples_trembles == 1 || num_trembles_to_est == 0 ) ){
    sample_of_ids_max.fill(1);
  }
  int num_sample_of_ids_max = max( sample_of_ids_max );
  arma::mat state_obs( num_rows_response_mat*num_sample_of_ids_max , num_cols_response_mat , arma::fill::zeros );
  arma::vec new_shares = shares;
  arma::mat new_share_mat = share_mat;
  arma::mat new_responses = responses;
  arma::vec new_trembles = trembles;
  int eval = eval_pre;

  // create indices responses cube from mat
  int indices_responses_max = indices_responses.max();
  arma::cube indices_responses_cube( num_rows_response_mat , num_cols_response_mat , num_samples_responses , arma::fill::zeros );
  for (int i = 0; i < num_samples_responses; i++) {
    arma::mat indices_slice = indices_responses;
    indices_slice( find( indices_responses > 0 ) ) += i*(indices_responses_max);
    indices_responses_cube.slice(i) = indices_slice;
  }

  // create and fill response cube
  arma::cube response_cube( num_rows_response_mat , num_cols_response_mat , num_samples_responses , arma::fill::zeros );
  for (int i = 0; i < num_samples_responses; i++) {
    response_cube.slice(i) = response_mat;
  }
  for (int i = 0; i < num_responses_to_est ; i++) {
    response_cube( find( indices_responses_cube == i+1 ) ).fill( new_responses(i) );
  }

  // normalize response rows
  arma::mat normalized_response_slice = response_mat;
  for (int j = 0; j < num_samples_responses; j++) {
    arma::mat response_slice = response_cube.slice(j);
    normalized_response_slice = response_slice;
    for (int i = 0; i < num_rows_response_mat; i++) { //
      arma::rowvec row_to_normalize = normalized_response_slice.row(i);
      arma::rowvec row_indices_responses = indices_responses.row(i);
      arma::uvec fixed_indices_row = find( row_indices_responses == 0 );
      arma::uvec target_indices_row = find( row_indices_responses > 0 );
      if ( target_indices_row.n_elem > 0 ){
        double sum_density = sum( row_to_normalize( target_indices_row ) );
        arma::rowvec sum_density_vec = row_to_normalize;
        sum_density_vec.fill( sum_density );
        row_to_normalize( target_indices_row ) /= sum_density_vec( target_indices_row );
        if( fixed_indices_row.n_elem > 0 ){
          double fixed_density = sum( row_to_normalize( fixed_indices_row ) );
          arma::rowvec scale_row_vec = row_to_normalize;
          scale_row_vec.fill( 1 - fixed_density );
          row_to_normalize( target_indices_row ) %= scale_row_vec( target_indices_row );
        }
      }
      normalized_response_slice.row(i) = row_to_normalize;
    }
    response_cube.slice(j) = normalized_response_slice;
  }

  for (int i = 0; i < num_responses_to_est; i++) {
    arma::mat new_discrete_response_value = unique( response_cube( find( indices_responses_cube == i+1 ) ) );
    new_responses(i) = new_discrete_response_value(0,0);
  }

  // create indices trembles cube from mat
  int indices_trembles_max = indices_trembles.max();
  arma::cube indices_trembles_cube( num_rows_tremble_mat , num_cols_tremble_mat , num_samples_trembles , arma::fill::zeros );
  for (int i = 0; i < num_samples_trembles; i++) {
    arma::mat indices_slice = indices_trembles;
    indices_slice( find( indices_trembles > 0 ) ) += i*(indices_trembles_max);
    indices_trembles_cube.slice(i) = indices_slice;
  }

  // create and fill tremble cube
  arma::cube tremble_cube( tremble_mat.n_rows , tremble_mat.n_cols , num_samples_trembles , arma::fill::zeros );
  for (int i = 0; i < num_samples_trembles; i++) {
    tremble_cube.slice(i) = tremble_mat;
  }
  for (int i = 0; i < num_trembles_to_est ; i++) {
    tremble_cube( find( indices_trembles_cube == i+1 ) ).fill( new_trembles(i) );
  }

  // create indices response id cube & indices trembles tube
  arma::cube indices_responses_id_cube( num_rows_response_mat , num_cols_response_mat , num_ids , arma::fill::zeros );
  if( num_samples_responses > 1 ){
    for (int i = 0; i < num_ids; i++) {
      arma::mat indices_responses_id = indices_responses;
      indices_responses_id( find( indices_responses > 0 ) ) += indices_responses_max*(sample_of_ids_responses(i)-1);
      indices_responses_id_cube.slice(i) = indices_responses_id;
    }
  }else{
    for (int i = 0; i < num_ids; i++) {
      indices_responses_id_cube.slice(i) = indices_responses;
    }
  }
  arma::cube indices_trembles_tube( num_rows_response_mat , 1 , num_ids , arma::fill::zeros );
  if( num_samples_trembles > 1 ){
    for (int i = 0; i < num_ids; i++) {
      arma::vec indices_trembles_id = indices_trembles.col(0);
      indices_trembles_id( find( indices_trembles_id > 0 ) ) += indices_trembles_max*(sample_of_ids_trembles(i)-1);
      indices_trembles_tube.slice(i) = indices_trembles_id;
    }
  }else{
    for (int i = 0; i < num_ids; i++) {
      indices_trembles_tube.slice(i) = indices_trembles.col(0);
    }
  }

  // calculate remaining non-fixed response to distribute
  arma::mat incomplete_response_mat = response_mat;
  incomplete_response_mat( indices_non_fixed_responses ).fill(0);
  arma::vec remaining_response_vec = 1 - sum( incomplete_response_mat , 1 );
  arma::mat remaining_response_mat = repmat( remaining_response_vec , 1 , num_cols_response_mat );

  while (  eval < max_eval+eval_pre && ( eps < 0 || eps >= tol_eval ) && eps != arma::datum::nan ) {
    eval++;
    Rcpp::checkUserInterrupt();

    // parameters are assigned to updated parameter values
    share_mat = new_share_mat;
    trembles = new_trembles;
    responses = new_responses;
    ll_val = new_ll_val;
    entropy_k = new_entropy_k;

    // new parameters are zero
    new_share_mat.elem( shares_to_est ).fill(0);
    new_responses.fill(0);
    new_trembles.fill(0);
    new_ll_val(0) = 0;
    new_entropy_k(0) = 0;

    // loop through ids
    for (int i = 0; i < num_ids; i++) {

      arma::mat response_mat_id = response_cube.slice( sample_of_ids_responses(i) - 1 );
      arma::mat tremble_mat_id = tremble_cube.slice( sample_of_ids_trembles(i) - 1 );

      // calculate emission probabilities with trembles
      arma::mat pr_mat = response_mat_id % (1 - tremble_mat_id) + ( 1 - response_mat_id ) % ( tremble_mat_id / (num_outputs - 1) );   // probability for emission with trembles

      // calculate the probability for each outcome in each state of each strategy
      arma::mat pr_outcomes_states_strategies_entities_mat( num_rows_response_mat , num_cols_response_mat , arma::fill::zeros );
      for ( int l = 0; l < num_cols_response_mat; l++){
        for (int j = 0; j < num_rows_response_mat; j++){
          pr_outcomes_states_strategies_entities_mat(j,l) = pow( pr_mat(j,l) , output_cube(j,l,i) );
        }
      }
      arma::vec pr_states_strategies_entities_mat = prod( pr_outcomes_states_strategies_entities_mat , 1 );
      arma::vec pr_entity_k( k , arma::fill::zeros );
      for ( int j = 0; j < k; j++){
        pr_entity_k(j) = prod( pr_states_strategies_entities_mat( find( strat_id == j+1 ) ) );
      }

      // log likelihood contribution of subject
      pr_entity_k %= share_mat( arma::span::all , sample_of_ids_shares(i) - 1 );
      new_ll_val += -log( sum( pr_entity_k , 0 ) );

      // share contribution of subject as posterior probability of i to use k / N
      arma::vec i_shares = pr_entity_k.each_row() / sum( pr_entity_k , 0 );
      arma::uvec shares_to_est_col = find( indices_shares( arma::span::all , sample_of_ids_shares(i) - 1 ) != 0 );
      arma::vec new_shares_col = new_share_mat( arma::span::all , sample_of_ids_shares(i) - 1 );
      new_shares_col( shares_to_est_col ) += ( i_shares( shares_to_est_col )  / num_ids_sample_shares( sample_of_ids_shares(i) - 1 ) );
      new_share_mat( arma::span::all , sample_of_ids_shares(i) - 1 ) = new_shares_col;

      i_shares_mat.row(i) = i_shares.t();
      arma::mat lines_entity_k = repmat( i_shares , 1 , num_cols_response_mat );
      arma::mat entity_slice( num_rows_response_mat , num_cols_response_mat );
      for ( int l = 0; l < num_rows_response_mat; l++){
        entity_slice.row(l) = lines_entity_k.row( strat_id( l ) - 1 );
      }
      i_shares_cube.slice(i) = entity_slice;

      // entropy k contribution
      arma::rowvec i_entropy_k = i_shares.t() % log( i_shares.t() );
      i_entropy_k.replace(arma::datum::nan, 0);
      new_entropy_k -= sum( i_entropy_k , 1 );

    }

    // correct shares for remainder
    for (int j = 0; j < num_samples; j++) {
      if( num_fixed_shares_col(j) > 0 && num_fixed_shares_col(j) != k ){
        arma::uvec shares_to_est_col = find( indices_shares( arma::span::all , j ) != 0 );
        arma::vec new_shares_col = new_share_mat( arma::span::all , j );
        arma::vec remaining_shares_col = remaining_shares_mat( arma::span::all , j );
        arma::vec shares_of_remaining_shares = new_shares_col( shares_to_est_col );
        new_shares_col( shares_to_est_col ) = remaining_shares_col( shares_to_est_col ) %  ( shares_of_remaining_shares / accu( shares_of_remaining_shares ) );
        new_share_mat( arma::span::all , j ) = new_shares_col;
      }
    }

    // update responses
    weigthed_output_cube = i_shares_cube % output_cube;
    weigthed_sums_cube = i_shares_cube % sum_outputs_cube;
    for ( int j = 0; j < num_responses_to_est; j++){
      new_responses(j) = accu( weigthed_output_cube( find( indices_responses_id_cube == j+1 ) ) ) / accu( weigthed_sums_cube( find( indices_responses_id_cube == j+1 ) ) );
    }
    new_responses.replace(arma::datum::nan, -1 );                             // clean responses (-1 indicates no obs)

    // fill responses cube with new values
    for (int i = 0; i < num_responses_to_est ; i++) {
      response_cube( find( indices_responses_cube == i+1 ) ).fill( new_responses(i) );
    }

    // fill responses cube with responses to sum
    for (int i = 0; i < num_samples_responses ; i++) {
      arma::mat response_cube_slice = response_cube.slice(i);
      response_cube_slice( indices_responses_to_sum ).fill(0);
      arma::vec value_to_sum_vec = 1 - sum( response_cube_slice , 1 );
      arma::mat value_to_sum_mat = repmat( value_to_sum_vec , 1 , num_cols_response_mat );
      response_cube_slice( indices_responses_to_sum ) = value_to_sum_mat( indices_responses_to_sum );
      response_cube.slice(i) = response_cube_slice;
    }

    // transform into discrete responses
    if ( response == "pure" && num_responses_to_est > 0 ){
      for (int j = 0; j < num_samples_responses; j++) {
        arma::umat indices_rows_with_estimated = find( sum( indices_responses_cube.slice(j) , 1 )  > 0 );
        arma::mat response_cube_slice = response_cube.slice(j);
        arma::mat target_mat = response_cube_slice.rows( indices_rows_with_estimated );
        arma::mat response_mat_with_estimated = response_cube_slice.rows( indices_rows_with_estimated );
        int num_response_mat_with_estimated = response_mat_with_estimated.n_rows;
        for (int i = 0; i < num_response_mat_with_estimated; i++) {
          arma::rowvec target_row = response_mat_with_estimated.row( i );
          arma::uword index_of_max = index_max( target_row );
          target_row.fill(0);
          target_row( index_of_max ) = 1;
          target_mat.row(i) = target_row;
        }
        response_cube_slice.rows( indices_rows_with_estimated ) = target_mat;
        response_cube.slice(j) = response_cube_slice;
        for (int i = 0; i < num_responses_to_est; i++) {
          arma::mat new_discrete_response_value = unique( response_cube( find( indices_responses_cube == i+1 ) ) );
          new_responses(i) = new_discrete_response_value(0,0);
        }
      }
    }
    response_mat = response_cube.slice(0);

    //update trembles
    arma::cube response_cube_id( num_rows_response_mat , num_cols_response_mat , num_ids , arma::fill::zeros );
    arma::cube tremble_correction_factor_cube_id( num_rows_response_mat , num_cols_response_mat , num_ids, arma::fill::ones );
    if( num_samples_responses > 1 ){
      for (int i = 0; i < num_ids; i++) {
        arma::mat response_cube_id_slice = response_cube.slice( sample_of_ids_responses(i) - 1 );
        arma::mat tremble_correction_factor_mat( num_rows_response_mat , num_cols_response_mat , arma::fill::ones );
        tremble_correction_factor_mat( find( response_cube_id_slice == 0 ) ).fill( num_cols_response_mat-1 );
        tremble_correction_factor_mat( find( response_cube_id_slice == 1 ) ).fill( -1 );
        tremble_correction_factor_cube_id.slice(i) = tremble_correction_factor_mat;
        response_cube_id.slice(i) = response_cube_id_slice;
      }
    }else{
      response_cube_id.each_slice() = response_mat;
      arma::mat tremble_correction_factor_mat( num_rows_response_mat , num_cols_response_mat , arma::fill::ones );
      tremble_correction_factor_mat( find( response_mat == 0 ) ).fill( num_cols_response_mat-1 );
      tremble_correction_factor_mat( find( response_mat == 1 ) ).fill( -1 );
      tremble_correction_factor_cube_id.each_slice() = tremble_correction_factor_mat;
    }
    arma::cube weigthed_output_sum_diffs_cube = tremble_correction_factor_cube_id % i_shares_cube % ( output_cube - ( sum_outputs_cube % response_cube_id )  );
    arma::cube weigthed_output_sum_diffs_tube = sum( weigthed_output_sum_diffs_cube , 1 );
    arma::cube weigthed_sums_tube = sum( weigthed_sums_cube , 1 );
    for ( int j = 0; j < num_trembles_to_est; j++){
      new_trembles(j) = accu( weigthed_output_sum_diffs_tube( find( indices_trembles_tube == j+1 ) ) ) / accu( weigthed_sums_tube( find( indices_trembles_tube == j+1 ) ) );
    }
    new_trembles.replace(arma::datum::nan, -1 );

    // fill tremble cube with new values
    for (int i = 0; i < num_trembles_to_est ; i++) {
      tremble_cube( find( indices_trembles_cube == i+1 ) ).fill( new_trembles(i) );
    }
    tremble_mat = tremble_cube.slice(0);

    // check overshooting and calculate eps for tolerance
    if (eval > eval_pre+1 ) { eps_now = (1 - (new_ll_val(0) / ll_val(0))); }          // current epsilon
    if ( new_ll_val(0) == 0 ){ eps_now = 0; }
    if ( new_ll_val.is_finite() ){  eps = eps_now; }                                  // only continue if no overshoot
    else { eps = arma::datum::nan; }                                                  // if overshooting occured report results from last eval

  } // end while

  // update new_shares
  for (int i = 0; i < num_shares_to_est; i++) {
    arma::vec updated_share = new_share_mat( find( indices_shares.t() == i+1 ) );
    if( updated_share.is_finite() ){
      arma::vec unique_updates_share = unique( updated_share );
      new_shares(i) = updated_share(0);
    }
  }

  // calculate state observations
  for (int i = 0; i < num_ids; i++) {
    arma::mat weigthed_output_cube_slice = weigthed_output_cube.slice(i);
    int sample_of_id = sample_of_ids_max(i);
    state_obs( arma::span( (sample_of_id-1)*num_rows_response_mat , sample_of_id*num_rows_response_mat - 1 ) , arma::span::all ) += weigthed_output_cube_slice( arma::span( 0 , num_rows_response_mat - 1 ) , arma::span::all );
  }

  // calculate selection criteria
  crit_vals(0) = 2*new_ll_val(0) + 2*free_params;                                                 // update AIC value
  crit_vals(1) = 2*new_ll_val(0) + log( (double) num_ids )*free_params;                           // update BIC value
  crit_vals(2) = crit_vals(1) + 2*new_entropy_k(0);                                                // update ICL value

  // prepare output
  double LL = new_ll_val(0);
  if( LL == arma::datum::nan ){
    LL = ll_val(0);
  }
  double E = entropy_k(0);

  F(0,0) = new_shares;
  F(1,0) = new_responses;
  F(2,0) = new_trembles;
  F(3,0) = new_share_mat;
  F(4,0) = response_mat;
  F(5,0) = tremble_mat;
  F(6,0) = LL;
  F(7,0) = crit_vals;
  F(10,0) = eval;
  F(11,0) = eps;
  F(12,0) = E;
  F(13,0) = i_shares_mat;
  F(21,0) = state_obs;
  F(22,0) = indices_responses;
  F(23,0) = indices_trembles;
  F(24,0) = indices_shares;
  F(25,0) = free_params;
  return(F);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// stratEst_LCR_EM
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

arma::field<arma::mat> stratEst_LCR_EM(arma::cube& output_cube, arma::cube& sum_outputs_cube, arma::vec& strat_id, arma::mat& covariate_mat, arma::vec shares, arma::vec responses, arma::vec trembles, arma::vec coefficients, arma::mat share_mat, arma::mat response_mat, arma::mat tremble_mat, arma::mat coefficients_mat , arma::uvec shares_to_est, arma::mat indices_responses, arma::mat indices_trembles, arma::mat indices_coefficients, bool estimate_coefficients, arma::uvec coefficients_to_est, arma::mat& responses_to_sum, std::string& response, int eval_pre , double tol_eval, int max_eval, double newton_stepsize, bool penalty, arma::vec& sample_of_ids_shares, arma::vec& num_ids_sample_shares, arma::vec& sample_of_ids_responses, arma::vec& num_ids_sample_responses, arma::vec& sample_of_ids_trembles, arma::vec& num_ids_sample_trembles ) {

  arma::field<arma::mat> F(26,1);
  int num_ids = output_cube.n_slices;
  int num_rows_response_mat = response_mat.n_rows;
  int num_cols_response_mat = response_mat.n_cols;
  int num_rows_tremble_mat = tremble_mat.n_rows;
  int num_cols_tremble_mat = tremble_mat.n_cols;
  int num_rows_coefficients_mat = coefficients_mat.n_rows;
  int num_cols_coefficients_mat = coefficients_mat.n_cols;
  // arma::vec unique_indices_responses = unique(indices_responses( find( indices_responses != 0 ) ));
  // int num_responses_to_est = unique_indices_responses.n_elem;
  // arma::vec unique_indices_trembles = unique(indices_trembles( find( indices_trembles != 0 ) ) );
  // int num_trembles_to_est = unique_indices_trembles.n_elem;
  int num_responses_to_est = responses.n_elem;
  int num_trembles_to_est = trembles.n_elem;
  int k = share_mat.n_rows;
  int num_outputs = num_cols_response_mat;
  int num_samples_responses = sample_of_ids_responses.max();
  int num_samples_trembles = sample_of_ids_trembles.max();
  int num_coefficients_to_est = coefficients.n_elem;
  int num_coefficients = coefficients_mat.n_rows*num_cols_coefficients_mat;
  int num_coefficients_short_vec = num_coefficients - coefficients_mat.n_rows;
  arma::umat indices_responses_to_sum = find( responses_to_sum == 1 );
  arma::umat indices_non_fixed_responses = find( indices_responses == 0 && responses_to_sum == 0 );
  int free_responses = 0;
  if( response == "pure" ){
    free_responses = ( num_responses_to_est/num_cols_response_mat )*num_samples_responses;
  }else{
    free_responses = ( num_responses_to_est - num_responses_to_est/num_cols_response_mat )*num_samples_responses;
  }
  int free_params = free_responses + num_trembles_to_est + num_coefficients_to_est;
  double eps = arma::datum::inf;
  double eps_now = arma::datum::inf;
  arma::vec ll_val( 1 , 1 , arma::fill::zeros );
  arma::vec entropy_k( 1 , arma::fill::ones);
  arma::vec crit_vals( 3 , arma::fill::ones);
  arma::vec new_ll_val( 1 , 1 , arma::fill::zeros );
  arma::vec next_ll_val( 1 , 1 , arma::fill::zeros );
  arma::vec new_entropy_k( 1 , arma::fill::ones);
  arma::mat i_shares_mat( num_ids , k , arma::fill::zeros );
  arma::cube i_shares_cube( num_rows_response_mat , num_cols_response_mat , num_ids , arma::fill::zeros );
  arma::vec sample_of_ids_max = sample_of_ids_responses;
  if( max(sample_of_ids_max) < max(sample_of_ids_trembles ) ){
    sample_of_ids_max = sample_of_ids_trembles;
  }
  if( ( num_samples_responses == 1 || num_responses_to_est == 0 ) & ( num_samples_trembles == 1 || num_trembles_to_est == 0 ) ){
    sample_of_ids_max.fill(1);
  }
  int num_sample_of_ids_max = max( sample_of_ids_max );
  arma::mat state_obs( num_rows_response_mat*num_sample_of_ids_max , num_cols_response_mat , arma::fill::zeros );
  //arma::mat state_obs( num_rows_response_mat*num_sample_of_ids_max , 1 , arma::fill::zeros );
  arma::vec new_shares = shares;
  arma::mat new_share_mat = share_mat;
  arma::mat new_responses = responses;
  arma::vec new_trembles = trembles;
  arma::vec new_coefficients = coefficients;
  arma::vec changes_coefficients( num_coefficients_to_est , arma::fill::zeros );
  arma::vec score_vec( num_coefficients , arma::fill::zeros );
  arma::vec short_score_vec( num_coefficients_to_est , arma::fill::zeros );
  arma::mat hessian_mat( num_coefficients , num_coefficients , arma::fill::zeros );
  arma::mat fisher_info( num_coefficients , num_coefficients , arma::fill::zeros );
  arma::mat penalized_fisher_info( num_coefficients_to_est , num_coefficients_to_est , arma::fill::zeros );
  arma::cube d_hessian_cube( num_coefficients , num_coefficients , num_coefficients , arma::fill::zeros );
  arma::cube delta_lk( num_coefficients , num_coefficients , num_coefficients , arma::fill::zeros );
  arma::cube delta_lh( num_coefficients , num_coefficients , num_coefficients , arma::fill::zeros );
  arma::cube delta_kh( num_coefficients , num_coefficients , num_coefficients , arma::fill::zeros );
  arma::mat score_contribution_mat( num_ids , num_coefficients , arma::fill::zeros );
  arma::mat penalized_score_contribution_mat( num_ids , num_coefficients_to_est , arma::fill::zeros );
  arma::mat priors_entities_mat( num_ids , k , arma::fill::zeros );
  arma::mat pr_entities_k_mat( num_ids , k , arma::fill::zeros );
  arma::vec stepsize_vec( num_coefficients_short_vec , arma::fill::ones );
  arma::vec penalty_vec( num_coefficients_short_vec , arma::fill::zeros );
  arma::cube weigthed_sums_cube( num_rows_response_mat , num_cols_response_mat , num_ids , arma::fill::zeros );
  //
  arma::cube weigthed_output_cube = output_cube;
  //
  stepsize_vec.fill(newton_stepsize);
  int eval = eval_pre;
  bool coefficients_changed = true;

  // create indices responses cube from mat
  int indices_responses_max = indices_responses.max();
  arma::cube indices_responses_cube( num_rows_response_mat , num_cols_response_mat , num_samples_responses , arma::fill::zeros );
  for (int i = 0; i < num_samples_responses; i++) {
    arma::mat indices_slice = indices_responses;
    indices_slice( find( indices_responses > 0 ) ) += i*(indices_responses_max);
    indices_responses_cube.slice(i) = indices_slice;
  }

  // create and fill response cube
  arma::cube response_cube( num_rows_response_mat , num_cols_response_mat , num_samples_responses , arma::fill::zeros );
  for (int i = 0; i < num_samples_responses; i++) {
    response_cube.slice(i) = response_mat;
  }
  for (int i = 0; i < num_responses_to_est ; i++) {
    response_cube( find( indices_responses_cube == i+1 ) ).fill( new_responses(i) );
  }

  // create indices trembles cube from mat
  int indices_trembles_max = indices_trembles.max();
  arma::cube indices_trembles_cube( num_rows_tremble_mat , num_cols_tremble_mat , num_samples_trembles , arma::fill::zeros );
  for (int i = 0; i < num_samples_trembles; i++) {
    indices_trembles_cube.slice(i) = indices_trembles + i*(indices_trembles_max);
  }

  // create and fill tremble cube
  arma::cube tremble_cube( tremble_mat.n_rows , tremble_mat.n_cols , num_samples_trembles , arma::fill::zeros );
  for (int i = 0; i < num_samples_trembles; i++) {
    tremble_cube.slice(i) = tremble_mat;
  }
  for (int i = 0; i < num_trembles_to_est ; i++) {
    tremble_cube( find( indices_trembles_cube == i+1 ) ).fill( new_trembles(i) );
  }

  // create indices response id cube & indidices trembles tube
  arma::cube indices_responses_id_cube( num_rows_response_mat , num_cols_response_mat , num_ids , arma::fill::zeros );
  if( num_samples_responses > 1 ){
    for (int i = 0; i < num_ids; i++) {
      arma::mat indices_responses_id = indices_responses;
      indices_responses_id += indices_responses_max*(sample_of_ids_responses(i)-1);
      indices_responses_id_cube.slice(i) = indices_responses_id;
    }
  }else{
    for (int i = 0; i < num_ids; i++) {
      indices_responses_id_cube.slice(i) = indices_responses;
    }
  }
  arma::cube indices_trembles_tube( num_rows_response_mat , 1 , num_ids , arma::fill::zeros );
  if( num_samples_trembles > 1 ){
    for (int i = 0; i < num_ids; i++) {
      arma::vec indices_trembles_id = indices_trembles.col(0);
      indices_trembles_id += indices_trembles_max*(sample_of_ids_trembles(i)-1);
      indices_trembles_tube.slice(i) = indices_trembles_id;
    }
  }else{
    for (int i = 0; i < num_ids; i++) {
      indices_trembles_tube.slice(i) = indices_trembles.col(0);
    }
  }

  for (int i = 0; i < num_coefficients_to_est ; i++) {
    coefficients_mat( find( indices_coefficients == i+1 ) ).fill( new_coefficients(i) );
  }

  arma::vec updated_coefficients = vectorise( coefficients_mat( arma::span::all , arma::span( 1 , indices_coefficients.n_cols - 1 ) ) );

  //fill prior entities mat
  priors_entities_mat = exp( covariate_mat * coefficients_mat );
  priors_entities_mat = priors_entities_mat.each_col() / sum( priors_entities_mat , 1 );


  while (  eval < max_eval+eval_pre && (eps < 0 || eps >= tol_eval ) && eps != arma::datum::nan && coefficients_changed ) {
    eval++;
    Rcpp::checkUserInterrupt();

    // parameters are assigned to updated parameter values
    share_mat = new_share_mat;
    trembles = new_trembles;
    responses = new_responses;
    coefficients = new_coefficients;
    ll_val = new_ll_val;
    entropy_k = new_entropy_k;

    // new parameters are zero
    new_share_mat.fill(0);
    new_responses.fill(0);
    new_trembles.fill(0);
    new_ll_val(0) = 0;
    new_entropy_k(0) = 0;
    score_vec.fill(0);
    hessian_mat.fill(0);
    fisher_info.fill(0);
    d_hessian_cube.fill(0);

      // loop through ids
      for (int i = 0; i < num_ids; i++) {

        arma::mat response_mat_id = response_cube.slice( sample_of_ids_responses(i) - 1 );
        arma::mat tremble_mat_id = tremble_cube.slice( sample_of_ids_trembles(i) - 1 );

         // calculate emission probabilities with trembles
        arma::mat pr_mat = response_mat_id % (1 - tremble_mat_id) + ( 1 - response_mat_id ) % ( tremble_mat_id / (num_outputs - 1) );   // probability for emission with trembles

        // retrieve prior of indvidual
        arma::rowvec priors_entities_row = priors_entities_mat.row(i);

        // calculate the probability for each outcome in each state of each strategy
        arma::mat pr_outcomes_states_strategies_entities_mat( num_rows_response_mat , num_cols_response_mat , arma::fill::zeros );
        for ( int l = 0; l < num_cols_response_mat; l++){
          for (int j = 0; j < num_rows_response_mat; j++){
            pr_outcomes_states_strategies_entities_mat(j,l) = pow( pr_mat(j,l) , output_cube(j,l,i) );
          }
        }
        arma::vec pr_states_strategies_entities_mat = prod( pr_outcomes_states_strategies_entities_mat , 1 );
        arma::vec pr_entity_k( k , arma::fill::zeros );
        for ( int j = 0; j < k; j++){
          pr_entity_k(j) = prod( pr_states_strategies_entities_mat( find( strat_id == j+1 ) ) );
        }

        // log likelihood contribution of subject
        pr_entities_k_mat.row(i) = pr_entity_k.t();
        pr_entity_k %= priors_entities_row.t();
        new_ll_val += -log( accu( pr_entity_k ) );

        // share contribution of subject as posterior probability of i to use k / N
        arma::vec i_shares = pr_entity_k.each_row() / sum( pr_entity_k , 0 );
        i_shares_mat.row(i) = i_shares.t();
        arma::mat lines_entity_k = repmat( i_shares , 1 , num_cols_response_mat );
        arma::mat entity_slice( num_rows_response_mat , num_cols_response_mat );
        for ( int l = 0; l < num_rows_response_mat; l++){
          entity_slice.row(l) = lines_entity_k.row( strat_id( l ) - 1 );
        }
        i_shares_cube.slice(i) = entity_slice;

        // entropy k contribution
        arma::rowvec i_entropy_k = i_shares.t() % log( i_shares.t() );
        i_entropy_k.replace(arma::datum::nan, 0);
        new_entropy_k -= sum( i_entropy_k , 1 );

        // score and hessian contributions
        arma::cube h_covariate_cube( num_coefficients , num_coefficients , num_coefficients , arma::fill::zeros  );
        arma::vec score_contribution( num_coefficients , arma::fill::zeros );
        arma::vec individual_prior = priors_entities_row.t();
        arma::vec s_covariate_vec = repmat( covariate_mat.row(i).t() , k , 1 );
        arma::mat h_covariate_mat = repmat( covariate_mat.row(i).t() * covariate_mat.row(i) , k , k );
        arma::mat h_i_shares( num_coefficients , num_coefficients , arma::fill::zeros );
        arma::vec h_i_shares_vec( num_coefficients , arma::fill::zeros );
        arma::mat h_i_prior( num_coefficients , num_coefficients , arma::fill::zeros );
        arma::vec h_i_prior_vec( num_coefficients , arma::fill::zeros );
        arma::vec h_ones_zeros_vec( num_coefficients , arma::fill::zeros );
        arma::mat h_ones_zeros( num_coefficients , num_coefficients , arma::fill::zeros );
        h_covariate_cube.each_slice() = h_covariate_mat;
        arma::cube third_dim_cube( num_coefficients , num_coefficients , num_coefficients , arma::fill::zeros  );
        arma::mat third_dim_mat( num_coefficients , num_coefficients , arma::fill::zeros  );
        for (int s = 0; s < num_coefficients; s++){
          third_dim_mat.fill(s_covariate_vec(s));
          third_dim_cube.slice(s) = third_dim_mat;
        }
        h_covariate_cube %= third_dim_cube;
        for (int m = 0; m < k; m++){
          for (int j = 0; j < k; j++){
            for ( int l = 0; l < num_rows_coefficients_mat ; l++){
              if( m == 0 ){
                h_i_shares_vec( l + j*num_rows_coefficients_mat ) = i_shares(j);
                h_i_prior_vec( l + j*num_rows_coefficients_mat ) = individual_prior(j);
              }
              if( j == m ){
                h_ones_zeros_vec( l + j*num_rows_coefficients_mat ) = 1;
              }
              else{
                h_ones_zeros_vec( l + j*num_rows_coefficients_mat ) = 0;
              }
            }
            if( m == 0 ){
              score_contribution = s_covariate_vec % ( h_i_shares_vec - h_i_prior_vec );
            }
          }
          h_i_shares( arma::span::all , arma::span(m*num_rows_coefficients_mat , (num_rows_coefficients_mat-1) + m*num_rows_coefficients_mat ) ) = repmat( h_i_shares_vec , 1 , num_rows_coefficients_mat );
          h_i_prior( arma::span::all , arma::span(m*num_rows_coefficients_mat , (num_rows_coefficients_mat-1) + m*num_rows_coefficients_mat ) ) = repmat( h_i_prior_vec , 1 , num_rows_coefficients_mat );
          h_ones_zeros( arma::span::all , arma::span(m*num_rows_coefficients_mat , (num_rows_coefficients_mat-1) + m*num_rows_coefficients_mat ) ) = repmat( h_ones_zeros_vec , 1 , num_rows_coefficients_mat );
        }

        // objects for derivative of hessian
        arma::cube l_i_shares_cube( num_coefficients , num_coefficients , num_coefficients , arma::fill::zeros  );
        arma::cube k_i_shares_cube( num_coefficients , num_coefficients , num_coefficients , arma::fill::zeros  );
        arma::cube h_i_shares_cube( num_coefficients , num_coefficients , num_coefficients , arma::fill::zeros  );
        arma::cube l_i_prior_cube( num_coefficients , num_coefficients , num_coefficients , arma::fill::zeros  );
        arma::cube k_i_prior_cube( num_coefficients , num_coefficients , num_coefficients , arma::fill::zeros  );
        arma::cube h_i_prior_cube( num_coefficients , num_coefficients , num_coefficients , arma::fill::zeros  );

        arma::mat l_i_shares_mat( num_coefficients , num_coefficients , arma::fill::zeros  );
        arma::mat k_i_shares_mat( num_coefficients , num_coefficients , arma::fill::zeros  );
        arma::mat h_i_shares_mat( num_coefficients , num_coefficients , arma::fill::zeros  );
        arma::mat l_i_prior_mat( num_coefficients , num_coefficients , arma::fill::zeros  );
        arma::mat k_i_prior_mat( num_coefficients , num_coefficients , arma::fill::zeros  );
        arma::mat h_i_prior_mat( num_coefficients , num_coefficients , arma::fill::zeros  );

        arma::rowvec h_i_shares_row = h_i_shares.col(0).t();
        arma::rowvec h_i_prior_row = h_i_prior.col(0).t();

        l_i_shares_mat.each_col() = h_i_shares_row.t();
        k_i_shares_mat.each_row() = h_i_shares_row;
        l_i_prior_mat.each_col() = h_i_prior_row.t();
        k_i_prior_mat.each_row() = h_i_prior_row;
        l_i_shares_cube.each_slice() = l_i_shares_mat;
        k_i_shares_cube.each_slice() = k_i_shares_mat;
        l_i_prior_cube.each_slice() = l_i_prior_mat;
        k_i_prior_cube.each_slice() = k_i_prior_mat;

        if( i == 0 ){
          delta_lk.each_slice() = h_ones_zeros;
          for (int l = 0; l < num_coefficients; l++){
            for (int k = 0; k < num_coefficients; k++){
              for (int h = 0; h < num_coefficients; h++){
                double elem = delta_lk( l , k , h );
                delta_lh( l , h , k ) = elem;
                delta_kh( h , l ,  k ) = elem;
              }
            }
          }
        }

        for (int s = 0; s < num_coefficients; s++){
          h_i_shares_mat.fill(h_i_shares_row(s));
          h_i_shares_cube.slice(s) = h_i_shares_mat;
          h_i_prior_mat.fill(h_i_prior_row(s));
          h_i_prior_cube.slice(s) = h_i_prior_mat;
        }

        // individual's hessian, score and fisher cotributions
        arma::mat hessian_contribution = h_covariate_mat % ( h_i_shares % ( h_ones_zeros - h_i_shares.t() )  - h_i_prior % ( h_ones_zeros - h_i_prior.t() ) );
        hessian_mat += hessian_contribution;
        score_vec += score_contribution;
        score_contribution_mat.row(i) = score_contribution.t();
        fisher_info += score_contribution * score_contribution.t();

        // individual's derivative of hessian contribution
        if( penalty ){
          arma::cube d_hessian_i_shares = l_i_shares_cube % ( delta_lk - k_i_shares_cube ) % ( delta_lh - h_i_shares_cube) - l_i_shares_cube % k_i_shares_cube % ( delta_kh - h_i_shares_cube );
          arma::cube d_hessian_i_prior = l_i_prior_cube % ( delta_lk - k_i_prior_cube ) % ( delta_lh - h_i_prior_cube) - l_i_prior_cube % k_i_prior_cube % ( delta_kh - h_i_prior_cube );
          arma::cube d_hessian_contribution = h_covariate_cube % ( d_hessian_i_shares + d_hessian_i_prior );
          d_hessian_cube += d_hessian_contribution;

        }

      }

      // update responses
      weigthed_output_cube = i_shares_cube % output_cube;
      weigthed_sums_cube = i_shares_cube % sum_outputs_cube;
      for ( int j = 0; j < num_responses_to_est; j++){
        new_responses(j) = accu( weigthed_output_cube( find( indices_responses_id_cube == j+1 ) ) ) / accu( weigthed_sums_cube( find( indices_responses_id_cube == j+1 ) ) );
      }
      new_responses.replace(arma::datum::nan, -1 );                             // clean responses (-1 indicates no obs)

      // fill responses cube with new values
      for (int i = 0; i < num_responses_to_est ; i++) {
        response_cube( find( indices_responses_cube == i+1 ) ).fill( new_responses(i) );
      }

      // fill responses cube with responses to sum
      for (int i = 0; i < num_samples_responses ; i++) {
        arma::mat response_cube_slice = response_cube.slice(i);
        response_cube_slice( indices_responses_to_sum ).fill(0);
        arma::vec value_to_sum_vec = 1 - sum( response_cube_slice , 1 );
        arma::mat value_to_sum_mat = repmat( value_to_sum_vec , 1 , num_cols_response_mat );
        response_cube_slice( indices_responses_to_sum ) = value_to_sum_mat( indices_responses_to_sum );
        response_cube.slice(i) = response_cube_slice;
      }

      // transform into discrete responses
      if ( response == "pure" && num_responses_to_est > 0 ){
        for (int j = 0; j < num_samples_responses; j++) {
          arma::umat indices_rows_with_estimated = find( sum( indices_responses_cube.slice(j) , 1 )  > 0 );
          arma::mat response_cube_slice = response_cube.slice(j);
          arma::mat target_mat = response_cube_slice.rows( indices_rows_with_estimated );
          arma::mat response_mat_with_estimated = response_cube_slice.rows( indices_rows_with_estimated );
          int num_response_mat_with_estimated = response_mat_with_estimated.n_rows;
          for (int i = 0; i < num_response_mat_with_estimated; i++) {
            arma::rowvec target_row = response_mat_with_estimated.row( i );
            arma::uword index_of_max = index_max( target_row );
            target_row.fill(0);
            target_row( index_of_max ) = 1;
            target_mat.row(i) = target_row;
          }
          response_cube_slice.rows( indices_rows_with_estimated ) = target_mat;
          response_cube.slice(j) = response_cube_slice;
          for (int i = 0; i < num_responses_to_est; i++) {
            arma::mat new_discrete_response_value = unique( response_cube( find( indices_responses_cube == i+1 ) ) );
            new_responses(i) = new_discrete_response_value(0,0);
          }
        }
      }
      response_mat = response_cube.slice(0);

      //update trembles
      arma::cube response_cube_id( num_rows_response_mat , num_cols_response_mat , num_ids , arma::fill::zeros );
      arma::cube tremble_correction_factor_cube_id( num_rows_response_mat , num_cols_response_mat , num_ids, arma::fill::ones );
      if( num_samples_responses > 1 ){
        for (int i = 0; i < num_ids; i++) {
          arma::mat response_cube_id_slice = response_cube.slice( sample_of_ids_responses(i) - 1 );
          arma::mat tremble_correction_factor_mat( num_rows_response_mat , num_cols_response_mat , arma::fill::ones );
          tremble_correction_factor_mat( find( response_cube_id_slice == 0 ) ).fill( num_cols_response_mat-1 );
          tremble_correction_factor_mat( find( response_cube_id_slice == 1 ) ).fill( -1 );
          tremble_correction_factor_cube_id.slice(i) = tremble_correction_factor_mat;
          response_cube_id.slice(i) = response_cube_id_slice;
        }
      }else{
        response_cube_id.each_slice() = response_mat;
        arma::mat tremble_correction_factor_mat( num_rows_response_mat , num_cols_response_mat , arma::fill::ones );
        tremble_correction_factor_mat( find( response_mat == 0 ) ).fill( num_cols_response_mat-1 );
        tremble_correction_factor_mat( find( response_mat == 1 ) ).fill( -1 );
        tremble_correction_factor_cube_id.each_slice() = tremble_correction_factor_mat;
      }
      arma::cube weigthed_output_sum_diffs_cube = tremble_correction_factor_cube_id % i_shares_cube % ( output_cube - ( sum_outputs_cube % response_cube_id )  );
      arma::cube weigthed_output_sum_diffs_tube = sum( weigthed_output_sum_diffs_cube , 1 );
      arma::cube weigthed_sums_tube = sum( weigthed_sums_cube , 1 );
      for ( int j = 0; j < num_trembles_to_est; j++){
        new_trembles(j) = accu( weigthed_output_sum_diffs_tube( find( indices_trembles_tube == j+1 ) ) ) / accu( weigthed_sums_tube( find( indices_trembles_tube == j+1 ) ) );
      }
      new_trembles.replace(arma::datum::nan, -1 );

      // fill tremble cube with new values
      for (int i = 0; i < num_trembles_to_est ; i++) {
        tremble_cube( find( indices_trembles_cube == i+1 ) ).fill( new_trembles(i) );
      }
      tremble_mat = tremble_cube.slice(0);


    // update coefficients
    if( estimate_coefficients && num_coefficients_to_est > 0 ){
          short_score_vec = score_vec( arma::span( num_rows_coefficients_mat , num_coefficients-1 ) );
          arma::mat lower_hessian_mat = hessian_mat( arma::span( num_rows_coefficients_mat , num_coefficients-1 ) , arma::span( num_rows_coefficients_mat , num_coefficients-1 ) );
          arma::mat lower_fisher = fisher_info( arma::span( num_rows_coefficients_mat , num_coefficients-1 ) , arma::span( num_rows_coefficients_mat , num_coefficients-1 ) );
          arma::mat inverted_mat( num_coefficients_to_est , num_coefficients_to_est , arma::fill::none );
          arma::mat to_invert_mat = -lower_hessian_mat;                       // observed fisher is negative hesiian since we are minimizing the log likelihood
      if( pinv( inverted_mat , to_invert_mat ) ){
              if( penalty == true ){
                for(int c = 0; c < num_coefficients_to_est; c++){
                  arma::mat d_hessian_slice = d_hessian_cube.slice(num_rows_coefficients_mat+c);
                  arma::mat lower_d_hessian_mat = d_hessian_slice( arma::span( num_rows_coefficients_mat , num_coefficients-1 ) , arma::span( num_rows_coefficients_mat , num_coefficients-1 ) );
                  penalty_vec(c) = trace(inverted_mat*lower_d_hessian_mat)/2;
                }
                short_score_vec += penalty_vec;      // minus because the negative log likelihood is minimized
              }
        int correction_step = 0;
        do{
          correction_step = correction_step + 1;
          changes_coefficients = stepsize_vec % ( inverted_mat*short_score_vec );
          updated_coefficients = updated_coefficients + changes_coefficients;
          arma::mat updated_coefficients_mat = coefficients_mat;
          updated_coefficients_mat( arma::span::all , arma::span( 1 , num_cols_coefficients_mat - 1 ) ) = reshape( updated_coefficients , num_rows_coefficients_mat , num_cols_coefficients_mat - 1  );
          for (int i = 0; i < num_coefficients_to_est ; i++) {
            arma::vec updated_coefficient_value = unique( updated_coefficients_mat( find( indices_coefficients == i+1 ) ) );
            coefficients_mat( find( indices_coefficients == i+1 ) ).fill( updated_coefficient_value(0) );
            new_coefficients(i) =  updated_coefficient_value(0);
          }

          // update prior entities mat
          priors_entities_mat = exp( covariate_mat * coefficients_mat );
          priors_entities_mat = priors_entities_mat.each_col() / sum( priors_entities_mat , 1 );

          //update prior ll value
          next_ll_val(0) = 0;
          for (int i = 0; i < num_ids; i++) {
            next_ll_val += -log( sum( pr_entities_k_mat.row(i) % priors_entities_mat.row(i) ) );
          }
          if( next_ll_val(0) > new_ll_val(0) ){
            stepsize_vec.fill(0.5*stepsize_vec(0));
          }
          else{
            stepsize_vec.fill(pow(stepsize_vec(0),0.5));
          }
        }
        while ( next_ll_val(0) > new_ll_val(0) && correction_step < 2 );
        if( penalty ){
          new_ll_val = new_ll_val - log( det(-lower_hessian_mat) )/2;     // minus because ll_val is negative log likelihood
        }
        // Rcout<< "ll val: \n"  << next_ll_val << "\n";
        // Rcout<< "score vec: \n"  << short_score_vec << "\n";
        // Rcout<< "stepsize: \n"  << stepsize_vec << "\n";
      }
      else{
        eval = max_eval+eval_pre;
        new_ll_val.fill(arma::datum::nan);
      }
    }
    else{
      new_coefficients = coefficients;
    }

    // update share mat
    arma::mat new_i_shares_mat = new_share_mat;
    for (int i = 0; i < num_ids; i++) {
      arma::rowvec i_shares_vec = i_shares_mat.row(i);
      arma::rowvec pr_entities_vec = priors_entities_mat.row(i);
      new_share_mat( arma::span::all , sample_of_ids_shares(i) - 1 ) += ( pr_entities_vec.t()  / num_ids_sample_shares( sample_of_ids_shares(i) - 1 ) );
      new_i_shares_mat( arma::span::all , sample_of_ids_shares(i) - 1 ) += ( i_shares_vec.t()  / num_ids_sample_shares( sample_of_ids_shares(i) - 1 ) );
    }

    // check overshooting and calculate eps for tolerance
    if( estimate_coefficients && num_coefficients_to_est > 0 ){
      if( eval > eval_pre+1 ){
        if( max(abs(new_coefficients-coefficients)/abs(coefficients)) <= tol_eval ){
          coefficients_changed = false;
        }
      }
      if (eval > eval_pre+1 ) { eps_now = (1 - (new_ll_val(0) / ll_val(0))); }          // current epsilon
      //if ( eps_now < tol_eval ){ eps_now = tol_eval; }
      if ( new_ll_val(0) == 0 ){ eps_now = 0; }
      if ( new_ll_val.is_finite() ){  eps = eps_now; }                                  // only continue if no overshoot
      else { eps = arma::datum::nan; }
    }
    else{
      if (eval > eval_pre+1 ) { eps_now = (1 - (new_ll_val(0) / ll_val(0))); }          // current epsilon
      // if ( eps_now < tol_eval ){ eps_now = tol_eval; }
      if ( new_ll_val(0) == 0 ){ eps_now = 0; }
      if ( new_ll_val.is_finite() ){  eps = eps_now; }                                  // only continue if no overshoot
      else { eps = arma::datum::nan; }
    }                                                                                   // if overshooting occured report results from last eval

  } // end while

  // update new_shares
  new_shares = vectorise( new_share_mat );

  // calculate state observations
  for (int i = 0; i < num_ids; i++) {
    arma::mat weigthed_output_cube_slice = weigthed_output_cube.slice(i);
    int sample_of_id = sample_of_ids_max(i);
    state_obs( arma::span( (sample_of_id-1)*num_rows_response_mat , sample_of_id*num_rows_response_mat - 1 ) , arma::span::all ) += weigthed_output_cube_slice( arma::span( 0 , num_rows_response_mat - 1 ) , arma::span::all );

    // arma::mat weigthed_sums_cube_slice = weigthed_sums_cube.slice(i);
    // int sample_of_id = sample_of_ids_max(i);
    // state_obs( arma::span( (sample_of_id-1)*num_rows_response_mat , sample_of_id*num_rows_response_mat - 1 ) , arma::span::all ) += weigthed_sums_cube_slice( arma::span( 0 , num_rows_response_mat - 1 ) , 0 );
  }

  // calculate selection criteria
  crit_vals(0) = 2*new_ll_val(0) + 2*free_params;                                                 // update AIC value
  crit_vals(1) = 2*new_ll_val(0) + log( (double) num_ids )*free_params;                           // update BIC value
  crit_vals(2) = crit_vals(1) + 2*new_entropy_k(0);                                                // update ICL value

  // prepare output
  double LL = new_ll_val(0);
  if( LL == arma::datum::nan ){
    LL = ll_val(0);
  }
  double E = entropy_k(0);

  F(0,0) = new_shares;
  F(1,0) = new_responses;
  F(2,0) = new_trembles;
  F(3,0) = new_share_mat;
  F(4,0) = response_mat;
  F(5,0) = tremble_mat;
  F(6,0) = LL;
  F(7,0) = crit_vals;
  F(10,0) = eval;
  F(11,0) = eps;
  F(12,0) = E;
  F(13,0) = i_shares_mat;
  F(14,0) = new_coefficients;
  F(15,0) = coefficients_mat;
  F(16,0) = priors_entities_mat;
  F(17,0) = short_score_vec;
  F(18,0) = hessian_mat;
  F(19,0) = fisher_info;
  F(20,0) = score_contribution_mat;
  F(21,0) = state_obs;
  F(22,0) = indices_responses;
  F(23,0) = indices_trembles;
  F(24,0) = indices_responses;
  F(25,0) = free_params;
  return(F);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// stratEst_SE
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

arma::field<arma::mat> stratEst_SE(arma::cube& output_cube, arma::cube& sum_outputs_cube, arma::vec& strat_id, arma::mat& covariate_mat, arma::mat share_mat, arma::vec responses, arma::vec trembles, arma::vec coefficients, arma::mat response_mat, arma::mat tremble_mat , arma::mat coefficient_mat , arma::mat i_shares_mat, arma::mat fisher_info_coefficients, arma::vec score_coefficients, arma::mat hessian_mat, arma::mat individual_priors_mat, arma::mat indices_shares, arma::mat indices_responses, arma::mat indices_trembles, std::string response, arma::mat score_contribution_mat, bool CL, arma::vec cluster_id_vec, bool LCR, arma::vec& sample_of_ids_shares, arma::vec& num_ids_sample_shares, arma::vec& sample_of_ids_responses, arma::vec& num_ids_sample_responses, arma::vec& sample_of_ids_trembles, arma::vec& num_ids_sample_trembles ) {

  arma::field<arma::mat> F(17,1);
  int num_ids = output_cube.n_slices;
  int num_responses_to_est = responses.n_elem;
  int num_trembles_to_est = trembles.n_elem;
  int num_cols_response_mat = response_mat.n_cols;
  int num_rows_response_mat = response_mat.n_rows;
  int num_rows_tremble_mat = tremble_mat.n_rows;
  int num_cols_tremble_mat = tremble_mat.n_cols;
  int num_rows_coefficient_mat = coefficient_mat.n_rows;
  int num_cols_coefficient_mat = coefficient_mat.n_cols;
  int num_coefficients_to_est = coefficients.n_elem;
  int num_coefficients = num_rows_coefficient_mat*num_cols_coefficient_mat;
  int k = share_mat.n_rows;
  int num_samples = share_mat.n_cols;
  int num_samples_responses = sample_of_ids_responses.max();
  int num_samples_trembles = sample_of_ids_trembles.max();
  int num_outputs = num_cols_response_mat;
  arma::vec shares_vec = vectorise( share_mat );
  arma::vec indices_shares_vec = vectorise( indices_shares );
  arma::uvec shares_to_est = find( indices_shares_vec > 0 );
  int num_shares_to_est = shares_to_est.n_elem;
  arma::mat convergence( 1 , 4 , arma::fill::ones );
  convergence.fill(-1);
  arma::mat remaining_shares_mat( k , num_samples , arma::fill::zeros );
  arma::vec num_fixed_shares_col( num_samples , arma::fill::zeros );
  for (int i = 0; i < num_samples; i++) {
    arma::vec shares_col = share_mat( arma::span::all , i );
    arma::vec indices_shares_col = indices_shares( arma::span::all , i );
    arma::uvec fixed_shares_col = find( indices_shares_col == 0 );
    double remaining_shares = 1 - accu( shares_col( fixed_shares_col ) );
    remaining_shares_mat( arma::span::all , i ).fill( remaining_shares );
    num_fixed_shares_col(i) = fixed_shares_col.n_elem;
  }

  arma::mat score_responses_slice( num_rows_response_mat , num_cols_response_mat , arma::fill::zeros );
  arma::mat score_trembles_slice( num_rows_response_mat , num_cols_response_mat , arma::fill::zeros );
  arma::mat score_contribution_responses_mat( num_ids , num_responses_to_est , arma::fill::zeros );
  arma::mat score_contribution_trembles_mat( num_ids , num_trembles_to_est , arma::fill::zeros );

  // initialize objects
  arma::mat SE_shares( k , num_samples , arma::fill::zeros );
  arma::mat shares_covar( k*num_samples , k*num_samples , arma::fill::zeros );
  arma::vec SE_responses( num_responses_to_est , arma::fill::zeros );
  arma::mat responses_covar( num_responses_to_est , num_responses_to_est , arma::fill::zeros );
  arma::vec SE_trembles( num_trembles_to_est , arma::fill::zeros );
  arma::mat trembles_covar( num_trembles_to_est  , num_trembles_to_est  , arma::fill::zeros );
  arma::vec SE_coefficients( num_coefficients - num_rows_coefficient_mat , arma::fill::zeros );
  arma::mat coefficients_covar( num_coefficients - num_rows_coefficient_mat  , num_coefficients - num_rows_coefficient_mat , arma::fill::zeros );

  arma::vec score_shares( k*num_samples , arma::fill::zeros );
  arma::vec score_responses( num_responses_to_est , arma::fill::zeros );
  arma::vec score_trembles( num_trembles_to_est , arma::fill::zeros );

  arma::mat fisher_info_shares( k*num_samples , k*num_samples , arma::fill::zeros );
  arma::mat fisher_info_responses( num_responses_to_est , num_responses_to_est , arma::fill::zeros );
  arma::mat fisher_info_trembles( num_trembles_to_est , num_trembles_to_est , arma::fill::zeros );

  // SE of mixed responses & trembles
  arma::mat r_weight_mat( num_rows_response_mat , num_cols_response_mat , arma::fill::zeros );
  arma::mat one_weight_mat( num_rows_response_mat , num_cols_response_mat , arma::fill::ones );
  r_weight_mat.fill( num_cols_response_mat - 1 );
  arma::mat tremble_weight_mat = ( one_weight_mat - response_mat)/r_weight_mat - response_mat;

  // create indices responses cube from mat
  int indices_responses_max = indices_responses.max();
  arma::cube indices_responses_cube( num_rows_response_mat , num_cols_response_mat , num_samples_responses , arma::fill::zeros );
  for (int i = 0; i < num_samples_responses; i++) {
    arma::mat indices_slice = indices_responses;
    indices_slice( find( indices_responses > 0 ) ) += i*(indices_responses_max);
    indices_responses_cube.slice(i) = indices_slice;
  }

  // create and fill response cube
  arma::cube response_cube( num_rows_response_mat , num_cols_response_mat , num_samples_responses , arma::fill::zeros );
  for (int i = 0; i < num_samples_responses; i++) {
    response_cube.slice(i) = response_mat;
  }
  for (int i = 0; i < num_responses_to_est ; i++) {
    response_cube( find( indices_responses_cube == i+1 ) ).fill( responses(i) );
  }

  // create indices trembles cube from mat
  int indices_trembles_max = indices_trembles.max();
  arma::cube indices_trembles_cube( num_rows_tremble_mat , num_cols_tremble_mat , num_samples_trembles , arma::fill::zeros );
  for (int i = 0; i < num_samples_trembles; i++) {
    arma::mat indices_slice = indices_trembles;
    indices_slice( find( indices_trembles > 0 ) ) += i*(indices_trembles_max);
    indices_trembles_cube.slice(i) = indices_slice;
  }

  // create and fill tremble cube
  arma::cube tremble_cube( tremble_mat.n_rows , tremble_mat.n_cols , num_samples_trembles , arma::fill::zeros );
  for (int i = 0; i < num_samples_trembles; i++) {
    tremble_cube.slice(i) = tremble_mat;

  }
  for (int i = 0; i < num_trembles_to_est ; i++) {
    tremble_cube( find( indices_trembles_cube == i+1 ) ).fill( trembles(i) );
  }

  // // correct i_shares_mat for remainder
  // arma::mat i_shares_mat_corr = i_shares_mat;
  // for (int i = 0; i < num_ids ; i++) {
  //   arma::rowvec individual_row = i_shares_mat.row( i );
  //   int share_index_id = ( sample_of_ids_shares( i ) - 1 );
  //   if( num_fixed_shares_col( share_index_id ) > 0 && num_fixed_shares_col( share_index_id ) != k ){
  //     arma::vec individual_row_corr = individual_row.t();
  //     arma::uvec shares_to_est_col = find( indices_shares( arma::span::all , share_index_id ) != 0 );
  //     arma::vec remaining_shares_col = remaining_shares_mat( arma::span::all , share_index_id );
  //     individual_row_corr( shares_to_est_col ) = remaining_shares_col( shares_to_est_col ) %  ( individual_row_corr( shares_to_est_col ) / accu( individual_row_corr( shares_to_est_col ) ) );
  //     individual_row = individual_row_corr.t();
  //   }
  //   i_shares_mat.row( i ) = individual_row;
  // }

  for (int i = 0; i < num_ids ; i++) {
    arma::vec score_contribution_responses( num_responses_to_est , arma::fill::zeros );
    arma::vec score_contribution_trembles( num_trembles_to_est , arma::fill::zeros );
    arma::mat lines_entity_k = repmat( i_shares_mat.row(i).t() , 1 , num_cols_response_mat );
    arma::mat entity_slice( num_rows_response_mat , num_cols_response_mat );
    for ( int l = 0; l < num_rows_response_mat; l++){
      entity_slice.row(l) = lines_entity_k.row( strat_id( l ) - 1 );
    }

    // fisher information of responses
    if( response == "mixed" && num_responses_to_est > 0 ){
      arma::mat response_mat_id = response_cube.slice( sample_of_ids_responses(i) - 1 );
      score_responses_slice = entity_slice % ( output_cube.slice(i) - ( response_mat_id % sum_outputs_cube.slice(i) ) );
      for (int j = 0; j < num_responses_to_est ; j++) {
        double contrib = accu( score_responses_slice( find( indices_responses_cube.slice( sample_of_ids_responses(i) - 1 ) == j+1 ) ) );
        if( arma::is_finite(contrib) ){
          score_contribution_responses(j) += contrib;
        }
      }
      score_contribution_responses_mat.row(i) = score_contribution_responses.t();
      score_responses += score_contribution_responses;
      fisher_info_responses += score_contribution_responses * score_contribution_responses.t();
    }

    // fisher information for trembles
    if( num_trembles_to_est > 0 ){
      arma::mat tremble_mat_id = tremble_cube.slice( sample_of_ids_trembles(i) - 1 );
      arma::mat pr_mat = response_mat % (1 - tremble_mat_id) + ( 1 - response_mat ) % ( tremble_mat_id / (num_outputs - 1) );   // probability for emission with trembles
      score_trembles_slice =  entity_slice % tremble_weight_mat % ( output_cube.slice(i)  /  pr_mat  ) ;
      for (int j = 0; j < num_trembles_to_est ; j++) {
        double contrib = accu( score_trembles_slice( find( indices_trembles_cube.slice( sample_of_ids_trembles(i) - 1 ) == j+1 ) ) );
        if( arma::is_finite(contrib) ){
          score_contribution_trembles(j) += contrib;
        }
      }
      score_contribution_trembles_mat.row(i) = score_contribution_trembles.t();
      score_trembles += score_contribution_trembles;
      fisher_info_trembles += score_contribution_trembles * score_contribution_trembles.t();
    }
  }

  // SEs of responses
  if( response == "mixed" && num_responses_to_est > 0 ){
    arma::mat inverse_fisher_info_responses = fisher_info_responses;
    if( false ){
      arma::vec unique_clusters = unique( cluster_id_vec );
      int num_clusters = unique_clusters.n_elem;
      arma::cube score_contributions_clusters( num_responses_to_est ,num_responses_to_est , num_clusters , arma::fill::zeros );
      for (int c = 0; c < num_clusters ; c++) {
        arma::mat score_contribution_cluster = sum( score_contribution_responses_mat.rows( find( cluster_id_vec == unique_clusters(c) ) ) , 0 );
        score_contributions_clusters.slice(c) = score_contribution_cluster.t() * score_contribution_cluster;
      }
      arma::mat meat_mat = sum( score_contributions_clusters , 2 );
      if( pinv( inverse_fisher_info_responses , fisher_info_responses ) ){
        responses_covar = inverse_fisher_info_responses * meat_mat * inverse_fisher_info_responses;
        SE_responses = sqrt( responses_covar.diag() );
      }
      else{
        responses_covar.fill(arma::datum::nan);
        SE_responses.fill(arma::datum::nan);
      }
    }
    else{
      if( pinv( inverse_fisher_info_responses , fisher_info_responses ) ){
        arma::mat identity_responses( num_responses_to_est , num_responses_to_est , arma::fill::eye );
        arma::mat rows_mat_responses = repmat( responses , 1 , num_responses_to_est );
        arma::mat cols_mat_responses = repmat( responses.t() , num_responses_to_est , 1 );
        arma::mat no_sum_constraint_responses = cols_mat_responses;
        for (int i = 0; i < num_responses_to_est ; i++) {
          for (int j = 0; j < num_responses_to_est ; j++) {
            arma::rowvec unique_indices_responses_row( num_cols_response_mat , arma::fill::zeros );
            for (int l = 0; l < num_rows_response_mat ; l++) {
              arma::rowvec indices_responses_row = indices_responses.row(l);
              if( any( indices_responses_row == i+1 ) ){
                indices_responses_row( find( indices_responses_row == i+1 ) );
                unique_indices_responses_row = indices_responses_row;
                l = num_rows_response_mat;
              }
            }
            if( any( unique_indices_responses_row == j+1 ) ){
              no_sum_constraint_responses(i,j) = 0;
            }
          }
        }
        identity_responses += no_sum_constraint_responses;
        arma::mat jacobian_responses = rows_mat_responses % ( identity_responses - cols_mat_responses );
        arma::mat covar_responses_mat = jacobian_responses * inverse_fisher_info_responses * jacobian_responses.t() ;
        responses_covar = covar_responses_mat;
        SE_responses = sqrt( covar_responses_mat.diag() );
      }
      else{
        SE_responses.fill(arma::datum::nan);
        responses_covar.fill(arma::datum::nan);
      }
    }

  }

  // SEs of trembles
  if( num_trembles_to_est > 0 ){
    arma::mat inverse_fisher_info_trembles = fisher_info_trembles;
    if( false ){
      arma::vec unique_clusters = unique( cluster_id_vec );
      int num_clusters = unique_clusters.n_elem;
      arma::cube score_contributions_clusters( num_trembles_to_est , num_trembles_to_est , num_clusters , arma::fill::zeros );
      for (int c = 0; c < num_clusters ; c++) {
        arma::mat score_contribution_cluster = sum( score_contribution_trembles_mat.rows( find( cluster_id_vec == unique_clusters(c) ) ) , 0 );
        score_contributions_clusters.slice(c) = score_contribution_cluster.t() * score_contribution_cluster;
      }
      arma::mat meat_mat = sum( score_contributions_clusters , 2 );
      if( pinv( inverse_fisher_info_trembles , fisher_info_trembles ) ){
        trembles_covar = inverse_fisher_info_trembles * meat_mat * inverse_fisher_info_trembles;
        SE_trembles = sqrt( trembles_covar.diag() );
      }
      else{
        trembles_covar.fill(arma::datum::nan);
        SE_trembles.fill(arma::datum::nan);
      }
    }
    else{
      arma::mat inverse_fisher_info_trembles = fisher_info_trembles;
      if( pinv( inverse_fisher_info_trembles , fisher_info_trembles ) ){
        trembles_covar = inverse_fisher_info_trembles;
        SE_trembles = sqrt( inverse_fisher_info_trembles.diag() );
      }
      else{
        SE_trembles.fill(arma::datum::nan);
      }
    }
  }

  if( LCR ){
    // SEs of coefficients
    arma::mat lower_hessian_mat = hessian_mat( arma::span( num_rows_coefficient_mat , num_coefficients-1 ) , arma::span( num_rows_coefficient_mat , num_coefficients-1 ) );
    arma::mat inverse_neg_lower_hessian_mat = -lower_hessian_mat;
    if( false ){
      arma::vec unique_clusters = unique( cluster_id_vec );
      int num_clusters = unique_clusters.n_elem;
      arma::cube score_contributions_clusters( num_coefficients , num_coefficients , num_clusters , arma::fill::zeros );
      for (int c = 0; c < num_clusters ; c++) {
        arma::mat score_contribution_cluster = sum( score_contribution_mat.rows( find( cluster_id_vec == unique_clusters(c) ) ) , 0 );
        score_contributions_clusters.slice(c) = score_contribution_cluster.t() * score_contribution_cluster;
      }
      arma::mat meat_mat = sum( score_contributions_clusters , 2 );
      arma::mat lower_meat_mat = meat_mat( arma::span( num_rows_coefficient_mat , num_coefficients-1 ) , arma::span( num_rows_coefficient_mat , num_coefficients-1 ) );
      if( pinv( inverse_neg_lower_hessian_mat , -lower_hessian_mat ) ){
        coefficients_covar = inverse_neg_lower_hessian_mat * lower_meat_mat * inverse_neg_lower_hessian_mat;
        SE_coefficients = sqrt( coefficients_covar.diag() );
      }
      else{
        coefficients_covar.fill(arma::datum::nan);
        SE_coefficients.fill(arma::datum::nan);
      }
    }
    else{
      arma::mat inverse_neg_lower_hessian_mat = lower_hessian_mat;
      if( pinv( inverse_neg_lower_hessian_mat , -lower_hessian_mat ) ){
        arma::mat covar_lower_mat_coefficients = inverse_neg_lower_hessian_mat;
        coefficients_covar = covar_lower_mat_coefficients;
        SE_coefficients = sqrt( covar_lower_mat_coefficients.diag() );
      }
      else{
        coefficients_covar.fill(arma::datum::nan);
        SE_coefficients.fill(arma::datum::nan);
      }
    }

    // SE of shares from covar mat coefficients via delta method
    arma::mat covar_mat_coefficients( num_coefficients , num_coefficients, arma::fill::zeros );
    arma::mat covar_lower_mat_coefficients = covar_mat_coefficients( arma::span( num_rows_coefficient_mat , num_coefficients-1 ) , arma::span( num_rows_coefficient_mat , num_coefficients-1 ) );
    arma::mat inverse_neg_hessian_mat = -hessian_mat;
    if( pinv( inverse_neg_hessian_mat , -hessian_mat ) ){
      covar_mat_coefficients = inverse_neg_hessian_mat;
    }
    else{
      covar_mat_coefficients.fill(arma::datum::nan);
    }
    int length_covariates = coefficient_mat.n_rows;
    arma::mat jacobian_mat_priors( k , num_coefficients , arma::fill::zeros );
    arma::mat zero_mat( k , num_coefficients , arma::fill::zeros );
    arma::rowvec ones_vec( length_covariates , arma::fill::ones );
    for (int j = 0; j < k ; j++) {
      zero_mat( j , arma::span( j*length_covariates , ( j*length_covariates + length_covariates - 1) ) ) = ones_vec;
    }
    for (int i = 0; i < num_ids ; i++) {
      arma::rowvec prior_i = individual_priors_mat.row(i);
      arma::mat prior_i_mat = repmat( prior_i.t() , 1 , num_coefficients );
      arma::mat prior_i_column_mat = repmat( prior_i.t() , 1 , length_covariates );
      arma::mat prior_i_column_vec( 1 , num_coefficients, arma::fill::zeros );
      arma::mat covariate_mat_i = repmat( covariate_mat.row(i) , k , k );
      for (int j = 0; j < k ; j++) {
        prior_i_column_vec( 0 , arma::span( (j*length_covariates) , (j*length_covariates + length_covariates - 1) ) ) = prior_i_column_mat.row(j);
      }
      arma::mat prior_i_column = repmat( prior_i_column_vec , k , 1 );
      jacobian_mat_priors += ( covariate_mat_i % (prior_i_mat % (zero_mat - prior_i_column) ) )/num_ids;
    }
    arma::mat covar_mat_shares = jacobian_mat_priors * covar_mat_coefficients * jacobian_mat_priors.t();
    SE_shares = sqrt( covar_mat_shares.diag() );
    shares_covar = covar_mat_shares;
  }
  else{  // LCR == false
    if( k > 1 ){
      for (int i = 0; i < num_ids ; i++) {
        int sample_of_ids_shares_i = sample_of_ids_shares(i);
        arma::rowvec individual_shares = share_mat.col(sample_of_ids_shares_i-1).t();
        arma::rowvec score_contributions_shares = i_shares_mat.row(i) - individual_shares;
        int start_index_id_fisher = k*(sample_of_ids_shares(i)-1);
        score_shares(arma::span(start_index_id_fisher,(start_index_id_fisher+k-1))) += score_contributions_shares.t();
        fisher_info_shares(arma::span(start_index_id_fisher,(start_index_id_fisher+k-1)),arma::span(start_index_id_fisher,(start_index_id_fisher+k-1))) += score_contributions_shares.t() * score_contributions_shares;
      }
      arma::mat inverse_fisher_info_shares = fisher_info_shares;
      if( pinv( inverse_fisher_info_shares , fisher_info_shares ) ){
        arma::mat identity_shares( k*num_samples , k*num_samples , arma::fill::eye );
        arma::mat rows_mat_shares( k*num_samples , k*num_samples , arma::fill::zeros );
        arma::mat cols_mat_shares( k*num_samples , k*num_samples , arma::fill::zeros );
        for (int sam = 0; sam < num_samples; sam++) {
          arma::vec shares_row( k*num_samples , arma::fill::zeros );
          shares_row( arma::span( k*sam , k*sam + k - 1 ) ) = shares_vec( arma::span( k*sam , k*sam + k - 1 ) );
          for (int line = (k*sam); line < (k*sam + k); line++) {
            rows_mat_shares.row(line) = shares_row.t();
            cols_mat_shares.col(line) = shares_row;
          }
        }
        arma::mat jacobian_shares = rows_mat_shares % ( identity_shares - cols_mat_shares ) ;
        arma::mat covar_shares_mat = jacobian_shares * inverse_fisher_info_shares * jacobian_shares.t() ;
        SE_shares = sqrt( covar_shares_mat.diag() );
        shares_covar = covar_shares_mat;
      }
      else{
        SE_shares.fill(arma::datum::nan);
        shares_covar.fill(arma::datum::nan);
      }
    }
  }

  // convergence checks based on first order condition
  if( num_shares_to_est  > 0 ){
    convergence(0,0) = max( abs( score_shares( shares_to_est ) ) );
  }
  if( num_responses_to_est > 0 ){
    convergence(0,1) = max(abs( score_responses ));
  }
  if( num_trembles_to_est > 0 ){
    convergence(0,2) = max(abs( score_trembles ));
  }
  if( LCR ){
    if( num_coefficients_to_est > 0 ){
      convergence(0,3) = max(abs( score_coefficients ));
    }
  }

  F(0,0) = SE_shares;
  F(1,0) = SE_responses;
  F(2,0) = SE_trembles;
  F(3,0) = SE_coefficients;
  F(4,0) = convergence;
  F(5,0) = shares_covar;
  F(6,0) = score_shares;
  F(7,0) = fisher_info_shares;
  F(8,0) = responses_covar;
  F(9,0) = score_responses;
  F(10,0) = fisher_info_responses;
  F(11,0) = trembles_covar;
  F(12,0) = score_trembles;
  F(13,0) = fisher_info_trembles;
  F(14,0) = coefficients_covar;
  F(15,0) = score_coefficients;
  F(16,0) = fisher_info_coefficients;

  return(F);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Main Function
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List stratEst_cpp(arma::mat data, arma::mat strategies, arma::vec sid , arma::mat shares , arma::vec trembles , arma::mat coefficient_mat, arma::mat covariates, bool LCR , arma::vec cluster, arma::vec quantile_vec , std::string response = "mixed", bool specific_shares = true, bool specific_responses = true, bool specific_trembles = true , bool specific_coefficients = true , std::string r_responses = "no", std::string r_trembles = "global", bool select_strategies = false, bool select_responses = false , bool select_trembles = false , int min_strategies = 1 , std::string crit = "bic", std::string SE = "analytic", int outer_runs = 10, double outer_tol_eval = 1e-10, int outer_max_eval = 1000, int inner_runs = 100, double inner_tol_eval = 1e-5, int inner_max_eval = 100, int LCR_runs = 100, double LCR_tol_eval = 1e-10, int LCR_max_eval = 1000, int BS_samples = 1000 , bool print_messages = true, bool integer_strategies = true, double newton_stepsize = 1 , bool penalty = false ) {

  if( print_messages == true ){
    Rcout<< "start estimation\n";
  }

  arma::field<arma::mat> R(26,1);
  int rows_data = data.n_rows;
  //int cols_data = data.n_cols;
  int rows_strategies = strategies.n_rows;
  int cols_strategies = strategies.n_cols;
  arma::vec state = strategies.col(0);
  int k = accu( state( find( state == 1 ) ) );
  int max_state = max( state );
  arma::vec n_ones( rows_data , 1 , arma::fill::ones );
  int LL_index = 6;
  int crit_index = 1;
  if( crit == "aic" ){
    crit_index = 0;
  }
  else if ( crit == "icl" ){
    crit_index = 2;
  }
  arma::mat complete_share_mat = shares;
  arma::mat complete_coefficients_mat = coefficient_mat;
  arma::mat failed_inner_runs( 1, 1 , arma::fill::zeros );
  arma::mat failed_outer_runs( 1, 1 , arma::fill::zeros );
  arma::mat failed_inner_LCR_runs( 1, 1 , arma::fill::zeros );
  arma::mat failed_outer_LCR_runs( 1, 1 , arma::fill::zeros );

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // check and transform function inputs
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // generate ids, supergames & rounds 1,2,...,max
  arma::vec id_vec = data.col(0);
  arma::vec match_vec = data.col(1);
  arma::vec round_vec = data.col(2);
  arma::uvec zeros_id = find( id_vec <= 0 );
  int num_zeros_id = zeros_id.n_elem;
  arma::uvec zeros_match = find( match_vec <= 0 );
  int num_zeros_match = zeros_match.n_elem;
  arma::uvec zeros_round = find( round_vec <= 0 );
  int num_zeros_round = zeros_round.n_elem;
  if ( num_zeros_id != 0 || num_zeros_match != 0 || num_zeros_round != 0 ){
    stop("The variables id, game (supergame), and period in of the data frame must contain values greater than zero.");
  }
  arma::vec unique_ids = unique( data.col(0) ) ;
  int num_ids = unique_ids.n_elem;
  for(int i = 0; i < num_ids; i++) {
    data.col(0).replace( unique_ids(i), i+1 );
  }
  id_vec = data.col(0);
  for(int i = 1; i <= num_ids; i++) {
    arma::vec matches_sbj = match_vec( find( data.col(0) == i ) );
    arma::vec unique_matches_sbj = unique( matches_sbj );
    int num_matches_sbj = unique_matches_sbj.n_elem;
    for(int j = 0; j < num_matches_sbj; j++) {
      matches_sbj.replace( unique_matches_sbj(j), j+1 );
      arma::vec rounds_sbj = round_vec( find( id_vec == i && match_vec == unique_matches_sbj(j) ) );
      arma::vec unique_rounds_sbj = unique( rounds_sbj );
      int num_unique_rounds_sbj = unique_rounds_sbj.n_elem;
      int num_rounds_sbj = rounds_sbj.n_elem;
      if( num_unique_rounds_sbj != num_rounds_sbj ){
        stop("The same period cannot occur several times in the same game for the same id.");
      }
      int num_rounds_sbj_match = unique_rounds_sbj.n_elem;
      for(int k = 0; k < num_rounds_sbj_match; k++) {
        rounds_sbj.replace( unique_rounds_sbj(k), k+1 );
      }
      round_vec( find( id_vec == i && match_vec == unique_matches_sbj(j)  ) ) = rounds_sbj;
    }
    match_vec( find( data.col(0) == i ) ) = matches_sbj;
  }
  data.col(1) = match_vec;
  data.col(2) = round_vec;

  // generate inputs 1,2,...,N
  arma::vec unique_inputs = unique( data.col(3) );
  arma::vec non_zero_unique_inputs = unique_inputs( find( unique_inputs != 0 ) );
  int num_non_zero_inputs = non_zero_unique_inputs.n_elem;
  for(int i = 0; i < num_non_zero_inputs; i++) {
    data.col(3).replace( non_zero_unique_inputs(i), i+1 );
  }

  // generate outputs 1,2,...,N
  arma::vec unique_outputs = unique( data.col(4) ) ;
  arma::vec non_zero_unique_outputs = unique_outputs( find( unique_outputs != 0 ) );
  int num_non_zero_outputs = non_zero_unique_outputs.n_elem;
  for(int i = 0; i < num_non_zero_outputs; i++) {
    data.col(4).replace( non_zero_unique_outputs(i) , i+1 );
  }

  //check fixed responses
  arma::vec output = data.col(4);
  arma::vec outputs = unique( output );
  int num_outputs = outputs.n_elem;
  arma::mat complete_response_mat( rows_strategies , num_outputs, arma::fill::zeros );
  arma::mat complete_tremble_mat( rows_strategies , num_outputs , arma::fill::zeros );
  arma::vec complete_tremble_vec( rows_strategies , arma::fill::zeros );
  complete_response_mat.fill(arma::datum::nan);
  complete_response_mat.tail_cols( num_non_zero_outputs ) = strategies.cols( 1 , num_non_zero_outputs );
  arma::uvec index_first_free( rows_strategies );
  for(int i = 0; i < rows_strategies; i++) {
    arma::rowvec responses_row = complete_response_mat.row( i );
    arma::uvec indices_free_responses_row = find_nonfinite( responses_row );
    if( indices_free_responses_row.n_elem > 0 ){
      index_first_free(i) = indices_free_responses_row(0);
      arma::vec fixed_responses_row = responses_row( find_finite( responses_row ) );
      arma::vec fixed_mixed_responses_row = fixed_responses_row( find( fixed_responses_row != 0 && fixed_responses_row != 1 ) );
      int num_fixed_responses_row =  fixed_responses_row.n_elem;
      if( num_fixed_responses_row > 0 ){
        if( arma::max( fixed_responses_row ) > 1+1e5  ){
          stop("Fixed responses cannot exceed one.");
        }
        else if( arma::max( fixed_responses_row ) > 1 && arma::max( fixed_responses_row ) <=  1+1e5 ){
          fixed_responses_row( find( fixed_responses_row == arma::max( fixed_responses_row ) ) ).fill(1);
        }
        if( arma::max( fixed_responses_row ) < 0  ){
          stop("Fixed responses cannot be negative.");
        }
        if( accu( fixed_responses_row ) > 1+1e5 ){
          stop("The sum of fixed shares cannot exceed one. stratEst cannot proceed with the current values.");
        }
        else if(  num_fixed_responses_row == num_outputs && accu( fixed_responses_row ) != 1 ){
          stop("The fixed shares must sum to one. stratEst cannot proceed with the current values.");
        }
        if( num_fixed_responses_row < (num_outputs-1)  && fixed_mixed_responses_row.n_elem != 0 && response == "pure" ){
          stop("It is not possible to estimate pure strategy parameters in a row where other parameters are fixed at mixed values. Estimate mixed parameters or change the fixed values to zero or one.");
        }
        if( num_fixed_responses_row < num_non_zero_outputs && ( r_responses != "no" || r_trembles != "no" ) ){
          stop("It is not possible to fix only a subset of response probabilities in a state of a strategy. Either fix all responses or no response.");
        }
      }
    }
  }

  // check samples
  arma::vec sample_id = data.col(5);
  arma::uvec zeros_samples = find( sample_id <= 0 );
  int num_zeros_samples = zeros_samples.n_elem;
  if ( num_zeros_samples != 0  ){
    stop("The variable sample in of the data frame must contain values greater than zero.");
  }
  for(int j = 0; j < num_ids; j++) {
    arma::vec unique_sample = unique( sample_id( find( id_vec == j+1 ) ) );
    int num_unique_sample = unique_sample.n_elem;
      if( num_unique_sample > 1 ){
        stop("The variable sample contains different values for the same subject.");
      }
  }

  // check fixed trembles
  arma::uvec indices_fixed_trembles = find_finite( trembles );
  int num_fixed_trembles = indices_fixed_trembles.n_elem;
  if( num_fixed_trembles >= 1 ){
    if( arma::max( trembles( indices_fixed_trembles ) ) > 1 ){
      stop("stratEst error: Fixed trembles cannot exceed one.");
    }
    if( arma::min( trembles( indices_fixed_trembles ) ) < 0  ){
      stop("stratEst error: Fixed trembles cannot be negative.");
    }
  }

  // generate samples 1,2,...,N
  arma::vec unique_samples = unique( sample_id ) ;
  int num_samples = unique_samples.n_elem ;
  for(int j = 0; j < num_samples; j++) {
    data.col(5).replace( unique_samples(j) , j+1 );
  }
  sample_id = data.col(5);

  //check fixed shares
  int num_rows_shares = complete_share_mat.n_rows;
  int num_cols_shares = complete_share_mat.n_cols;
  if( k !=  num_rows_shares ){
    stop("Matrix 'shares' must have as many rows as there are strategies.");
  }
  if( num_samples !=  num_cols_shares ){
    stop("Matrix 'shares' must have as many columns as there are samples.");
  }

  arma::vec fixed_shares = complete_share_mat( find_finite( complete_share_mat ) );
  int num_fixed_shares = fixed_shares.n_elem;
  if ( num_fixed_shares > 0 ){
    if( arma::max( fixed_shares ) > 1.001  ){
      stop("Fixed shares cannot exceed one.");
    }
    if( arma::min( fixed_shares ) < 0  ){
      stop("Fixed shares cannot be negative.");
    }
    arma::mat with_zero_share_mat = complete_share_mat;
    with_zero_share_mat.replace(arma::datum::nan, 0);
    arma::rowvec sums_of_shares = sum( with_zero_share_mat , 0 );
    if( any( sums_of_shares > 1.001 ) ){
      stop("The sum of all fixed values of a vector of 'shares' cannot exceed one. It is not possible to proceed with the current values.");
    }
    for(int j = 0; j < num_samples; j++) {
      arma::vec complete_share_col = complete_share_mat(arma::span::all,j);
      arma::vec fixed_shares_col = complete_share_col( find_finite( complete_share_col ) );
      int num_fixed_shares_col = fixed_shares_col.n_elem;
      if( num_fixed_shares_col == (k-1) ){
        complete_share_col.replace( arma::datum::nan , 1 - accu( fixed_shares_col ) );
      }
      complete_share_mat(arma::span::all,j) = complete_share_col;
    }
  }

  // check covariates
  int num_covariates = covariates.n_cols;
  arma::mat covariate_mat( num_ids , num_covariates , arma::fill::ones );
  if( LCR ){
    if( num_rows_shares == 1 ){
      stop("Latent class regression requires more than one strategy.");
    }
    if( num_samples > 1 ){
      stop("Latent class regression cannot be run with more than one sample.");
    }
    int cols_covariates = covariates.n_cols;
    int rows_covariates = covariates.n_rows;
    if( rows_covariates != rows_data  ){
      stop("Covariate matrix must have as many rows as data.");
    }
    if( num_fixed_shares > 0 ){
      stop("Shares cannot be fixed with covariates.");
    }
    arma::mat incomplete_covariate_mat( num_ids , covariates.n_cols , arma::fill::ones );
    for(int j = 0; j < num_ids; j++) {
      for(int c = 0; c < cols_covariates; c++) {
        arma::vec covariate = covariates.col(c);
        arma::vec unique_covariate = unique( covariate( find( id_vec == j+1 ) ) );
        int num_unique_covariate = unique_covariate.n_elem;
        if( num_unique_covariate > 1 ){
          stop("Covariate matrix contains different values of the same variable for the same id.");
        }
        else{
          incomplete_covariate_mat( j , c ) = unique_covariate(0);
        }
      }
    }
    // arma::mat intercept( num_ids , 1 , arma::fill::ones );
    // covariate_mat = join_rows( intercept , incomplete_covariate_mat );
    covariate_mat = incomplete_covariate_mat ;
  }

  // check cluster
  int num_clusters = 0;
  arma::vec cluster_id_vec( num_ids , arma::fill::zeros );
  bool CL = cluster.n_elem > 1;
  if( CL ){
    //SE = "bs";
    arma::uvec zeros_cluster = find( cluster <= 0 );
    int num_zeros_cluster = zeros_cluster.n_elem;
    if ( num_zeros_cluster != 0 ){
      stop("Cluster must contain values greater than zero.");
    }
    int elem_cluster = cluster.n_elem;
    if( elem_cluster != rows_data  ){
      stop("The cluster vector must have the same number of elements as there are rows in data.");
    }
    for(int j = 0; j < num_ids; j++) {
      arma::vec unique_cluster = unique( cluster( find( id_vec == j+1 ) ) );
      int num_unique_cluster = unique_cluster.n_elem;
      if( num_unique_cluster > 1 ){
        stop("The values of 'id' must be nested within clusters, i.e. the data of one individual appears only in one cluster.");
      }
      else{
        cluster_id_vec(j) = unique_cluster(0);
      }
    }
    // generate cluster 1,2,...,C
    arma::vec unique_cluster = unique( cluster );
    num_clusters = unique_cluster.n_elem;
    for(int j = 0; j < num_clusters; j++) {
      cluster.replace( unique_cluster(j), j+1 );
      cluster_id_vec.replace( unique_cluster(j), j+1 );
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // sort data matrix
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // // sort by id, supergame, round
  // arma::mat sorted_data( rows_data , cols_data , arma::fill::zeros );
  // int line = 0;
  // for( int i = 1; i <= num_ids; i++) {
  //   arma::vec matches_sbj = match_vec( find( data.col(0) == i ) );
  //   arma::vec unique_matches_sbj = unique( matches_sbj );
  //   int num_matches_sbj = unique_matches_sbj.n_elem;
  //   for(int j = 0; j < num_matches_sbj; j++) {
  //     matches_sbj.replace( unique_matches_sbj(j), j+1 );
  //     arma::vec rounds_sbj = round_vec( find( id_vec == i && match_vec == unique_matches_sbj(j) ) );
  //     arma::vec unique_rounds_sbj = unique( rounds_sbj );
  //     int num_rounds_sbj_match = unique_rounds_sbj.n_elem;
  //     for(int l = 0; l < num_rounds_sbj_match; l++) {
  //       sorted_data.row( line ) = data.rows( find( data.col(0) == i && data.col(1) == j+1 && data.col(2) == l+1 ) );
  //       line = line + 1;
  //     }
  //   }
  // }
  // data = sorted_data;


  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // initialize objects
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  arma::vec supergame = data.col(1);
  arma::vec period = data.col(2);
  arma::vec input = data.col(3);
  arma::vec inputs = unique( input );

  // calculate sample of ids for shares, responses, trembles
  arma::vec specific_sample_of_ids( num_ids , arma::fill::zeros );
  arma::vec sample_of_ids_shares( num_ids , arma::fill::zeros );
  arma::vec sample_of_ids_responses( num_ids , arma::fill::zeros );
  arma::vec sample_of_ids_trembles( num_ids , arma::fill::zeros );
  arma::vec num_ids_specific_sample( num_samples , arma::fill::zeros );
  arma::vec num_ids_sample_shares( num_samples , arma::fill::zeros );
  arma::vec num_ids_sample_responses( num_samples , arma::fill::zeros );
  arma::vec num_ids_sample_trembles( num_samples , arma::fill::zeros );
  int num_samples_responses = num_samples;
  int num_samples_trembles = num_samples;

  // sample specific objects
  if( specific_shares || specific_responses || specific_trembles || specific_coefficients ){
    for (int i = 0; i < num_ids; i++) {
      arma::vec unique_specific_sample_of_ids = unique( sample_id( find( id_vec == i+1 ) ) );
      specific_sample_of_ids(i) = unique_specific_sample_of_ids(0);
      num_ids_specific_sample( unique_specific_sample_of_ids(0) - 1 ) += 1;
    }
  }

  // sample specific shares
  if( specific_shares ){
    sample_of_ids_shares = specific_sample_of_ids;
    num_ids_sample_shares = num_ids_specific_sample;
  }else{
    sample_of_ids_shares.fill(1);
    num_ids_sample_shares.fill(num_ids);
  }

  // sample specific responses
  if( specific_responses ){
    sample_of_ids_responses = specific_sample_of_ids;
    num_ids_sample_responses = num_ids_specific_sample;
  }else{
    sample_of_ids_responses.fill(1);
    num_ids_sample_responses.fill(num_ids);
    num_samples_responses = 1;
  }

  // sample specific trembles
  if( specific_trembles ){
    sample_of_ids_trembles = specific_sample_of_ids;
    num_ids_sample_trembles = num_ids_specific_sample;
  }else{
    sample_of_ids_trembles.fill(1);
    num_ids_sample_trembles.fill(num_ids);
    num_samples_trembles = 1;
  }

  // calculate strategy ids
  arma::colvec complete_strat_id( rows_strategies , arma::fill::zeros );
  complete_strat_id( find( state == 1 ) ).ones();
  double current_id = 0;
  for (int i = 0; i < rows_strategies; i++) {
    current_id += complete_strat_id(i);
    complete_strat_id(i) = current_id;
  }

  // calculate state matrix for strategies
  arma::mat state_mat( rows_data , k , arma::fill::ones );
  arma::cube response_cube( max_state , num_non_zero_outputs , k , arma::fill::zeros );
  response_cube.fill(-1);
  for (int i = 0; i < k; i++) {
    arma::mat strategy = strategies.rows( find( complete_strat_id == (i+1) ) );
    response_cube( arma::span( 0 , ( strategy.n_rows - 1 ) ) , arma::span( 0 , (num_non_zero_outputs-1) ) , arma::span( i , i ) ) = strategy.cols( 1 , num_non_zero_outputs );
    for (int j = 0; j < rows_data; j++) {
      if( period(j) == 1){
        if( input(j) != 0 ){
          state_mat(j,i) = strategy( 0 , ( num_non_zero_outputs + input(j) ) );
        }
      }
      else{
        state_mat(j,i) = strategy( state_mat( j-1 , i ) - 1 , ( num_non_zero_outputs + input(j) ) ); //
      }
    }
  }

  // accumulate number of observed responses conditional on input (rows ids, cols states, slices strats)
  arma::cube complete_output_cube( rows_strategies , num_outputs , num_ids , arma::fill::zeros);
  arma::mat complete_sum_outputs_mat( rows_strategies , num_outputs , arma::fill::zeros);
  arma::cube complete_sum_outputs_cube( rows_strategies , num_outputs , num_ids , arma::fill::zeros);
  //arma::mat complete_output_ones( rows_data, output.n_cols , arma::fill::ones );
  // for( int l = 0; l < rows_strategies; l++) {
  //   for( int m = 0; m < num_outputs; m++) {
  //     for (int i = 0; i < num_ids; i++) {
  //       complete_output_cube(l,m,i) += accu( complete_output_ones.rows( find( id_vec == (i+1) && state_mat.col( complete_strat_id(l) - 1 ) == state(l) && output == unique_outputs(m) ) ) );
  //     }
  //   }
  // }
  // or
  // for (int i = 0; i < num_ids; i++) {
  //   arma::mat complete_output_slice(rows_strategies , num_outputs , arma::fill::zeros );
  //   for( int l = 0; l < rows_strategies; l++) {
  //     arma::rowvec complete_output_slice_row( num_outputs , arma::fill::zeros );
  //     arma::vec outputs_in_strategy_state = output( find( id_vec == (i+1) && state_mat.col( complete_strat_id(l) - 1 ) == state(l) ) );
  //     int num_outputs_in_strategy_state = outputs_in_strategy_state.n_elem;
  //     if( num_outputs_in_strategy_state > 0 ){
  //       for(int j = 0; j < num_outputs_in_strategy_state; j++){
  //         complete_output_slice_row( find( unique_outputs == outputs_in_strategy_state(j) ) ) += 1;
  //       }
  //       complete_output_slice.row(l) = complete_output_slice_row;
  //     }
  //   }
  //   complete_output_cube.slice(i) = complete_output_slice;
  // }

  //or
  int num_elem_outputs = output.n_elem;
  arma::mat output_as_mat( num_elem_outputs , num_outputs , arma::fill::zeros );
  for(int j = 0; j < num_elem_outputs; j++){
    arma::rowvec output_as_mat_row( num_outputs , arma::fill::zeros );
    output_as_mat_row( find( unique_outputs == output(j) ) ) += 1;
    output_as_mat.row(j) = output_as_mat_row;
  }
  for (int i = 0; i < num_ids; i++) {
    arma::mat complete_output_slice(rows_strategies , num_outputs , arma::fill::zeros );
    arma::uvec id_vec_i = id_vec == (i+1);
    for( int l = 0; l < rows_strategies; l++) {
      int complete_strat_id_l = complete_strat_id(l);
      arma::uvec state_mat_l = state_mat.col( complete_strat_id_l - 1 ) == state(l);
      arma::rowvec complete_output_slice_row( num_outputs , arma::fill::zeros );
      complete_output_slice_row = sum( output_as_mat.rows( find( state_mat_l && id_vec_i ) ) , 0 );
      complete_output_slice.row(l) = complete_output_slice_row;
    }
    complete_output_cube.slice(i) = complete_output_slice;
  }



  complete_sum_outputs_mat = sum( complete_output_cube , 1 );
  for (int i = 0; i < num_ids; i++) {
    complete_sum_outputs_cube.slice(i) = repmat( complete_sum_outputs_mat.col(i) , 1 , num_outputs );
  }

  // insert fixed trembles in complete tremble mat
  arma::mat trembles_mat = repmat( trembles , 1 , num_outputs );
  arma::umat fixed_trembles_mat = find_finite( trembles_mat );
  complete_tremble_mat( fixed_trembles_mat ) = trembles_mat( fixed_trembles_mat );

  // ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // indices and selection matrices
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // numbers of parameters to estimate
  int num_shares_to_est = 0;
  int num_responses_to_est = 0;
  int num_trembles_to_est = 0;
  int num_coefficients_to_est = 0;

  // generate complete indices shares
  arma::mat complete_indices_shares( k , num_samples , arma::fill::ones );
  if( specific_shares ){
    for (int i = 0; i < k; i++) {
      for (int j = 0; j < num_samples; j++) {
        complete_indices_shares(i,j) += i*num_samples + j;
      }
    }
  }
  else{
    for (int i = 0; i < k; i++) {
      for (int j = 0; j < num_samples; j++) {
        complete_indices_shares(i,j) += i;
      }
    }
  }

  // delete fixed shares from shares restriction mat
  arma::umat indices_fixed_shares = find_finite( complete_share_mat );
  complete_indices_shares( indices_fixed_shares ).fill(0);

  // normalize values in complete indices shares mat after deletion
  arma::vec unique_non_zero_restrictions_shares = unique( complete_indices_shares( find( complete_indices_shares != 0 ) ) );
  int num_unique_non_zero_restrictions_shares = unique_non_zero_restrictions_shares.n_elem;
  for (int i = 0; i < num_unique_non_zero_restrictions_shares; i++) {
    complete_indices_shares.replace( unique_non_zero_restrictions_shares(i) , i+1 );
  }

  arma::mat restriction_base_row( 1 , num_outputs , arma::fill::ones );
  for (int i = 0; i < num_outputs; i++) {
    restriction_base_row(0,i) += i;
  }
  arma::mat complete_indices_responses = repmat( restriction_base_row, rows_strategies , 1 );
  arma::mat complete_indices_trembles( rows_strategies , num_outputs , arma::fill::ones );
  if( r_responses == "no" ){
    for (int i = 0; i < rows_strategies; i++) {
      complete_indices_responses.row(i) += i*num_outputs;
    }
  }
  else if( r_responses == "strategies" ){
    for (int i = 0; i < rows_strategies; i++) {
      complete_indices_responses.row(i) += ( complete_strat_id(i) - 1 )*num_outputs;
    }
  }
  else if( r_responses == "states" ){
    for (int i = 0; i < rows_strategies; i++) {
      complete_indices_responses.row(i) += (state(i)-1)*num_outputs;
    }
  }
  if( r_trembles == "no" ){
    for (int i = 0; i < rows_strategies; i++) {
      complete_indices_trembles.row(i) += i;
    }
  }
  else if( r_trembles == "strategies" ){
    for (int i = 0; i < rows_strategies; i++) {
      complete_indices_trembles.row(i) += ( complete_strat_id(i) - 1 );
    }
  }
  else if( r_trembles == "states" ){
    for (int i = 0; i < rows_strategies; i++) {
      complete_indices_trembles.row(i) += (state(i)-1);
    }
  }
  arma::mat select_indices_responses = complete_indices_responses;
  arma::mat select_indices_trembles = complete_indices_trembles;

  // restriction matrices in case of selection
  if( select_responses ){
    complete_indices_responses = repmat( restriction_base_row, rows_strategies , 1 );
    for (int i = 0; i < rows_strategies; i++) {
      complete_indices_responses.row(i) += i*num_outputs;
    }
  }
  else if( select_trembles ){
    complete_indices_trembles.fill(1);
    for (int i = 0; i < rows_strategies; i++) {
      complete_indices_trembles.row(i) += i;
    }
  }

  // delete fixed responses and responses to sum from response restriction mat
  arma::umat indices_fixed_responses = find_finite( complete_response_mat );
  complete_indices_responses( indices_fixed_responses ).fill(0);
  arma::mat complete_responses_to_sum( complete_indices_responses.n_rows , complete_indices_responses.n_cols , arma::fill::zeros  );
  // //if the following is commented out response.par includes all parameters and not n-1 per row
  // for (int i = 0; i < rows_strategies; i++) {
  //   for (int j = 0; j < num_outputs; j++) {
  //     if( complete_indices_responses(i,j) > 0 ){
  //       complete_responses_to_sum(i,j) = 1;
  //       j = complete_indices_responses.n_cols;
  //     }
  //   }
  // }
  arma::umat complete_indices_responses_to_sum = find( complete_responses_to_sum == 1 );
  complete_indices_responses( complete_indices_responses_to_sum ).fill(0);
  arma::vec unique_non_zero_restrictions = unique( complete_indices_responses( find( complete_indices_responses != 0 ) ) );
  int num_unique_non_zero_restrictions = unique_non_zero_restrictions.n_elem;
  for (int i = 0; i < num_unique_non_zero_restrictions; i++) {
    complete_indices_responses.replace( unique_non_zero_restrictions(i) , i+1 );
  }

  // delete fixed trembles and mixed responses from tremble restriction vec
  arma::mat row_with_mixed_responses( rows_strategies , num_outputs , arma::fill::ones );
  arma::vec fixed_tremble( rows_strategies , arma::fill::ones );
  arma::mat row_with_fixed_tremble( rows_strategies , num_outputs , arma::fill::ones );
  fixed_tremble( find_nonfinite(trembles) ).fill(0);
  row_with_fixed_tremble = repmat( fixed_tremble , 1 , num_outputs );
  if( response == "mixed" ){
    for (int i = 0; i < rows_strategies; i++) {
      for (int j = 0; j < num_outputs; j++) {
        if( complete_response_mat(i,j) == 1 || complete_response_mat(i,j) == 0 ){
          row_with_mixed_responses.row(i).fill(0);
          j = complete_response_mat.n_cols;
        }
      }
    }
    complete_indices_trembles( find( row_with_mixed_responses == 1 ) ).fill(0);
  }
  else{
    arma::mat clean_response_mat = complete_response_mat;
    clean_response_mat.replace( arma::datum::nan, 0 );
    for (int i = 0; i < rows_strategies; i++) {
      for (int j = 0; j < num_outputs; j++) {
        if( clean_response_mat(i,j) == 1 || clean_response_mat(i,j) == 0 ){
          row_with_mixed_responses.row(i).fill(0);
          j = clean_response_mat.n_cols;
        }
      }
    }
    complete_indices_trembles( find( row_with_mixed_responses == 1 ) ).fill(0);
  }
  complete_indices_trembles( find( row_with_fixed_tremble == 1 ) ).fill(0);
  arma::vec unique_non_zero_trembles = unique( complete_indices_trembles( find( complete_indices_trembles != 0 ) ) );
  int num_unique_non_zero_trembles = unique_non_zero_trembles.n_elem;
  for (int i = 0; i < num_unique_non_zero_trembles; i++) {
    complete_indices_trembles.replace( unique_non_zero_trembles(i) , i+1 );
  }

  // check zero trembles
  for (int i = 0; i < rows_strategies; i++) {
    if( row_with_mixed_responses(i,0) == 1 && complete_tremble_mat(i,0) != 0 ){
      stop("stratEst error: Trembles of mixed response probabilities must be zero.");
    }
  }


  // generate complete indices coefficients
  arma::mat complete_indices_coefficients( num_covariates , k , arma::fill::zeros );
  arma::mat indices_coefficients( num_covariates , k , arma::fill::zeros );
  arma::umat indices_fixed_coefficients = find_finite( complete_coefficients_mat );
  if( LCR ){
    for (int i = 1; i < k; i++) {
      for (int j = 0; j < num_covariates; j++) {
        complete_indices_coefficients(j,i) += (i-1)*num_covariates + j + 1;
      }
    }

    // delete fixed coefficients from shares restriction mat and normalize
    complete_indices_coefficients( indices_fixed_coefficients ).fill(0);
    arma::vec unique_non_zero_restrictions_coefficients = unique( complete_indices_coefficients( find( complete_indices_coefficients != 0 ) ) );
    int num_unique_non_zero_restrictions_coefficients = unique_non_zero_restrictions_coefficients.n_elem;
    for (int i = 0; i < num_unique_non_zero_restrictions_coefficients; i++) {
      complete_indices_coefficients.replace( unique_non_zero_restrictions_coefficients(i) , i+1 );
    }
  }

  // preallocate result field
  int R_num_rows_strategies = 0;
  int R_num_strategies = 0;
  arma::cube R_output_cube( rows_strategies , num_outputs , num_ids , arma::fill::zeros );
  arma::cube R_sum_outputs_cube( rows_strategies , num_outputs , num_ids , arma::fill::zeros );
  arma::vec R_strat_id( rows_strategies , arma::fill::zeros );
  arma::mat R_shares_mat( k , num_samples , arma::fill::zeros );
  arma::mat R_indices_shares( rows_strategies , num_samples , arma::fill::zeros );
  arma::mat R_indices_responses( rows_strategies , num_outputs , arma::fill::zeros );
  arma::mat R_indices_trembles( rows_strategies , num_outputs , arma::fill::zeros );
  arma::mat R_responses_to_sum( rows_strategies , num_outputs , arma::fill::zeros );
  arma::mat R_indices_coefficients( num_covariates , k , arma::fill::zeros );
  arma::mat R_coefficient_mat( num_covariates , k , arma::fill::zeros );

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // strategy selection procedure
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  arma::mat complete_strat_id_mat = repmat( complete_strat_id , 1 , num_outputs );
  double crit_selected_min = arma::datum::inf;
  arma::vec survivors = unique( complete_strat_id );
  int kill_start = 0;

  // start strategies loop
  for (int K = k; K > 0; K--) {
    if( print_messages == true ){
      Rcout<< "model with " << K << " strategies\n";
    }
    int killed = 0;
    for (int kill = kill_start; kill <= K; kill++) {
      arma::vec survived = survivors( find( survivors != 0 ) );
      arma::vec zero_survivors( survived.n_elem + 1 , arma::fill::zeros );
      zero_survivors( arma::span( 1 , survived.n_elem ) ) = survived;
      arma::vec candidates = survived( find( survived != zero_survivors( kill ) ) );
      int num_candidates = candidates.n_elem;

      // create survivor objects
      arma::mat candidates_mat( rows_strategies , num_outputs , arma::fill::zeros );
      arma::cube candidates_cube( rows_strategies , num_outputs , num_ids , arma::fill::zeros );
      for (int j = 0; j < num_candidates; j++) {
        candidates_mat( find( complete_strat_id_mat == candidates(j) ) ).fill(1);
      }
      int num_rows_response_mat = accu( candidates_mat.col(0) );
      candidates_cube.each_slice() = candidates_mat;

      // create strat_id & share_mat & indices_shares
      arma::vec strat_id( num_rows_response_mat , arma::fill::zeros );
      arma::mat share_mat( num_candidates , num_samples, arma::fill::zeros );
      arma::mat indices_shares( num_candidates , num_samples , arma::fill::zeros );

      if( K == k && kill == 0 ){
        strat_id = complete_strat_id;
        share_mat = complete_share_mat;
        indices_shares = complete_indices_shares;
      }
      else{
        strat_id = complete_strat_id( find ( candidates_mat.col(0) == 1 ) );
        arma::vec unique_strat_id = unique( strat_id );
        int num_unique_strat_id = unique_strat_id.n_elem;
        for (int i = 0; i < num_unique_strat_id; i++) {
          strat_id.replace( unique_strat_id(i) , i+1 );
          share_mat( i , arma::span::all ) = complete_share_mat( unique_strat_id(i) - 1 , arma::span::all );
          indices_shares( i , arma::span::all ) = complete_indices_shares( unique_strat_id(i) - 1 , arma::span::all );
        }
      }

      arma::cube output_cube( num_rows_response_mat , num_outputs , num_ids , arma::fill::zeros );
      if( K == k && kill == 0 ){
        output_cube = complete_output_cube;
      }
      else{
        for (int i = 0; i < num_ids; i++) {
          arma::mat output_slice = complete_output_cube.slice(i);
          arma::vec output_cube_vec = output_slice( find ( candidates_mat == 1 ) );
          arma::mat output_mat = reshape( output_cube_vec , num_rows_response_mat , num_outputs );
          output_cube.slice(i) = output_mat;
        }
      }

      // create sum_outputs_cube
      arma::cube sum_outputs_cube( num_rows_response_mat , num_outputs , num_ids , arma::fill::zeros );
      if( K == k && kill == 0 ){
        sum_outputs_cube = complete_sum_outputs_cube;
      }
      else{
        for (int i = 0; i < num_ids; i++) {
          arma::mat sum_outputs_slice = complete_sum_outputs_cube.slice(i);
          arma::vec sum_outputs_cube_vec = sum_outputs_slice( find ( candidates_mat == 1 ) );
          arma::mat sum_outputs_mat = reshape( sum_outputs_cube_vec , num_rows_response_mat , num_outputs );
          sum_outputs_cube.slice(i) = sum_outputs_mat;
        }
      }

      // update indices_shares
      // update values in indices shares mat after deletion
      arma::vec unique_non_zero_restrictions_shares = unique( indices_shares( find( indices_shares != 0 ) ) );
      int num_unique_non_zero_restrictions_shares = unique_non_zero_restrictions_shares.n_elem;
      for (int i = 0; i < num_unique_non_zero_restrictions_shares; i++) {
        indices_shares.replace( unique_non_zero_restrictions_shares(i) , i+1 );
      }

      // create indices_responses
      arma::mat indices_responses( num_rows_response_mat , num_outputs , arma::fill::zeros );
      if( K == k && kill == 0 ){
        indices_responses = complete_indices_responses;
      }
      else{
        arma::vec indices_responses_vec = complete_indices_responses( find ( candidates_mat == 1 ) );
        indices_responses = reshape( indices_responses_vec , num_rows_response_mat , num_outputs );
        arma::vec unique_non_zero_indices_responses = unique( indices_responses( find( indices_responses != 0 ) ) );
        int num_unique_non_zero_indices_responses = unique_non_zero_indices_responses.n_elem;
        for (int i = 0; i < num_unique_non_zero_indices_responses; i++) {
          indices_responses.replace( unique_non_zero_indices_responses(i) , i+1 );
        }
      }

      // create indices_trembles
      arma::mat indices_trembles( num_rows_response_mat , num_outputs , arma::fill::zeros );
      if( K == k && kill == 0 ){
        indices_trembles = complete_indices_trembles;
      }
      else{
        arma::vec indices_trembles_vec = complete_indices_trembles( find ( candidates_mat == 1 ) );
        indices_trembles = reshape( indices_trembles_vec , num_rows_response_mat , num_outputs );
        arma::vec unique_non_zero_indices_trembles = unique( indices_trembles( find( indices_trembles != 0 ) ) );
        int num_unique_non_zero_indices_trembles = unique_non_zero_indices_trembles.n_elem;
        for (int i = 0; i < num_unique_non_zero_indices_trembles; i++) {
          indices_trembles.replace( unique_non_zero_indices_trembles(i) , i+1 );
        }
      }

      // create response_mat
      arma::mat response_mat( num_rows_response_mat , num_outputs , arma::fill::zeros );
      if( K == k && kill == 0 ){
        response_mat = complete_response_mat;
      }
      else{
        arma::vec response_mat_vec = complete_response_mat( find ( candidates_mat == 1 ) );
        response_mat = reshape( response_mat_vec , num_rows_response_mat , num_outputs );
      }

      // create tremble_mat
      arma::mat tremble_mat( num_rows_response_mat , num_outputs , arma::fill::zeros );
      if( K == k && kill == 0 ){
        tremble_mat = complete_tremble_mat;
      }
      else{
        arma::vec tremble_mat_vec = complete_tremble_mat( find ( candidates_mat == 1 ) );
        tremble_mat = reshape( tremble_mat_vec , num_rows_response_mat , num_outputs );
      }

      // create responses_to_sum
      arma::mat responses_to_sum( num_rows_response_mat , num_outputs , arma::fill::zeros );
      if( K == k && kill == 0 ){
        responses_to_sum = complete_responses_to_sum;
      }
      else{
        arma::vec responses_to_sum_vec = complete_responses_to_sum( find ( candidates_mat == 1 ) );
        responses_to_sum = reshape( responses_to_sum_vec , num_rows_response_mat , num_outputs );
      }
      arma::umat indices_responses_to_sum = find( responses_to_sum == 1 );

      // create coefficients_mat & indices_coefficients
      arma::mat coefficients_mat( num_covariates , num_candidates , arma::fill::zeros );
      arma::mat indices_coefficients_candidates( num_covariates , num_candidates , arma::fill::zeros );
      if( LCR ){
        if( K == k && kill == 0 ){
          coefficients_mat = complete_coefficients_mat;
          indices_coefficients = complete_indices_coefficients;
        }
        else{
          strat_id = complete_strat_id( find ( candidates_mat.col(0) == 1 ) );
          arma::vec unique_strat_id = unique( strat_id );
          int num_unique_strat_id = unique_strat_id.n_elem;
          for (int i = 0; i < num_unique_strat_id; i++) {
            strat_id.replace( unique_strat_id(i) , i+1 );
            coefficients_mat( arma::span::all , i ) = complete_coefficients_mat( arma::span::all , unique_strat_id(i) - 1  );
            indices_coefficients_candidates( arma::span::all , i ) = complete_indices_coefficients( arma::span::all , unique_strat_id(i) - 1 );
          }
          indices_coefficients_candidates.col(0).fill(0);
          arma::vec unique_non_zero_indices_coefficients = unique( indices_coefficients_candidates( find( indices_coefficients_candidates != 0 ) ) );
          int num_unique_non_zero_indices_coefficients = unique_non_zero_indices_coefficients.n_elem;
          for (int i = 0; i < num_unique_non_zero_indices_coefficients; i++) {
            indices_coefficients_candidates.replace( unique_non_zero_indices_coefficients(i) , i+1 );
          }
          indices_coefficients = indices_coefficients_candidates;
        }
      }

      // identify number of responses, trembles and shares to estimate
      num_shares_to_est = indices_shares.max();
      num_responses_to_est = indices_responses.max()*num_samples_responses;
      num_trembles_to_est = indices_trembles.max()*num_samples_trembles;
      num_coefficients_to_est = indices_coefficients.max();

      // generate indices
      arma::vec num_shares_to_est_col( num_samples , arma::fill::zeros );
      for (int j = 0; j < num_samples; j++) {
        arma::uvec shares_to_est_col = find_nonfinite( share_mat( arma::span::all , j ) );
        num_shares_to_est_col(j) = shares_to_est_col.n_elem;
      }

      // incomplete response mat where values of estimated parameters are added
      response_mat( find_nonfinite( response_mat ) ).fill(0);
      arma::mat incomplete_response_mat = response_mat;
      arma::vec remaining_response_vec = 1 - sum( incomplete_response_mat , 1 );
      arma::mat remaining_response_mat = repmat( remaining_response_vec , 1 , num_outputs );

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // start EM
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      // best result for k
      arma::field<arma::mat> R_k(26,1);

      // prepare outer runs (index or)
      arma::field<arma::mat> R_outer(26,1);
      double LL_outer_min = 0;

      // outer runs (index ro)
      for (int ro = 0; ro < outer_runs; ro++) {
        arma::field<arma::mat> R_inner(26,1);
        double LL_inner_min = 0;

        // inner runs (index ri)
        for (int ri = 0; ri < inner_runs; ri++) {
          if( print_messages == true ){
            if( kill == 0 ){
              if( select_strategies ){
                Rcout<< "estimate complete model (" << ro+1  << "/" << ri+1 << ")     \r";
              }
              else{
                Rcout<< "estimate model (" << ro+1  << "/" << ri+1 << ")     \r";
              }
            }
            else{
              if( integer_strategies == true ){
                Rcout<< "estimate model with " << K - 1 << " strategies " <<  " (" << ro+1  << "/" << ri+1 << ")     \r";
              }
              else{
                Rcout<< "estimate nested model " << kill <<  " (" << ro+1  << "/" << ri+1 << ")     \r";
              }
            }
          }
          //random starting points
          arma::mat start_share_mat = share_mat;
          arma::uvec shares_to_est = find( indices_shares > 0 );
          arma::vec zero_shares( num_shares_to_est , arma::fill::zeros );
          if( num_shares_to_est > 0 ){
            start_share_mat.elem( shares_to_est ).randu();
            for (int j = 0; j < num_samples; j++) {
              arma::vec start_shares_col = start_share_mat( arma::span::all , j );
              arma::vec indices_shares_col = indices_shares( arma::span::all , j );
              arma::uvec shares_to_est_col = find( indices_shares_col > 0 );
              if(  num_shares_to_est_col(j) == num_candidates ){
                start_share_mat( arma::span::all , j ) = start_shares_col / accu( start_shares_col );
              }
              else{
                arma::vec fixed_shares_col = start_shares_col( find( indices_shares_col == 0 ) );
                double sum_fixed_shares_col = accu( fixed_shares_col );
                double remaining_shares_col = 1 - sum_fixed_shares_col;
                arma::vec unfixed_shares_col = start_shares_col( shares_to_est_col );
                unfixed_shares_col = unfixed_shares_col / accu( unfixed_shares_col );
                start_shares_col( shares_to_est_col ) = unfixed_shares_col*remaining_shares_col;
                start_share_mat( arma::span::all, j ) = start_shares_col;
              }
            }
          }

          // random responses to start (normalized to remaining response)
          arma::vec start_responses = arma::randu( num_responses_to_est );

          // trembles to start
          arma::vec start_trembles = arma::randu( num_trembles_to_est );
          start_trembles /= 4;

          arma::field<arma::mat> R_temp = stratEst_EM( output_cube, sum_outputs_cube, strat_id, zero_shares, start_responses, start_trembles, start_share_mat, response_mat, tremble_mat, indices_shares, indices_responses, indices_trembles, responses_to_sum, response, 0 , inner_tol_eval, inner_max_eval, sample_of_ids_shares, num_ids_sample_shares, sample_of_ids_responses, num_ids_sample_responses, sample_of_ids_trembles, num_ids_sample_trembles );
          arma::mat LL_inner_temp = R_temp(LL_index,0);
          if( LL_inner_temp.has_nan() ){
            failed_inner_runs(0,0) += 1;
            LL_inner_temp.replace(arma::datum::nan, arma::datum::inf);
          }
          if ( ri == 0 || LL_inner_temp(0,0) < LL_inner_min ){
            R_inner = R_temp;
            LL_inner_min = LL_inner_temp(0,0);
          }
        }
        Rcpp::checkUserInterrupt();
        arma::mat pre_eval_vec = R_inner(10,0);
        R_outer = stratEst_EM( output_cube, sum_outputs_cube, strat_id, R_inner(0,0), R_inner(1,0), R_inner(2,0), R_inner(3,0), R_inner(4,0), R_inner(5,0), indices_shares, indices_responses, indices_trembles, responses_to_sum, response, pre_eval_vec(0,0) , outer_tol_eval, outer_max_eval, sample_of_ids_shares, num_ids_sample_shares, sample_of_ids_responses, num_ids_sample_responses, sample_of_ids_trembles, num_ids_sample_trembles );
        arma::mat LL_outer_mat = R_outer(LL_index,0);
        if( LL_outer_mat.has_nan() ){
          failed_outer_runs(0,0) += 1;
          LL_outer_mat.replace(arma::datum::nan, arma::datum::inf);
        }
        if ( ro == 0 || LL_outer_mat(0,0) < LL_outer_min ){
          R_k = R_outer;
          LL_outer_min = LL_outer_mat(0,0);
        }
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // response and tremble fusion procedures
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      arma::mat fused_indices_responses = indices_responses;
      arma::mat fused_indices_trembles = indices_trembles;
      arma::mat new_indices_responses = indices_responses;
      arma::mat new_indices_trembles = indices_trembles;

      arma::field<arma::mat> R_temp(26,1);
      arma::field<arma::mat> R_fused(26,1);

      int num_it_responses = 0;
      int num_it_trembles = 0;

      int trigger = 1;
      bool triggered = false;

      if( select_responses && num_responses_to_est > 1  ){

        if( print_messages == true ){
          Rcout << "\n";
          Rcout<< "select responses             \r";
        }

        R_fused = R_k;
        arma::vec crit_vals = R_k( 7 , 0 );
        double crit_min = crit_vals( crit_index );
        arma::mat pre_eval_vec = R_k(10,0);

        trigger = 1;
        int response_run = 1;

        // start response and tremble fusion loop
        while ( trigger == 1 ) {
          Rcpp::checkUserInterrupt();
          trigger = 0;
          triggered = false;

          arma::mat estimated_responses_mat = R_k(1,0);
          arma::vec estimated_responses = estimated_responses_mat.col(0);
          int num_estimated_responses = estimated_responses.n_elem;
          int num_col_indices = indices_responses.n_cols;
          int max_indices_responses = indices_responses.max();
          num_it_responses = max_indices_responses/num_col_indices;
          arma::vec responses_indices_vec = arma::linspace(1,num_estimated_responses,num_estimated_responses);
          if( num_it_responses > 1 ){
            int counter = 0;
             for (int r1 = 0; r1 < (num_it_responses-1); r1++) {
              for (int r2 = r1+1; r2 < num_it_responses; r2++) {
                counter = counter + 1;
                if( print_messages == true ){
                  Rcout<< "select responses (" << response_run << "/" << counter << ")     \r";
                }
                fused_indices_responses = indices_responses;
                arma::vec indices_ones( num_estimated_responses , arma::fill::zeros );
                arma::vec c1 = select_indices_responses( find( indices_responses  == r1*num_col_indices ) );
                arma::vec c2 = select_indices_responses( find( indices_responses  == r2*num_col_indices ) );
                if( approx_equal(c1, c2, "absdiff", 0.1 ) ){
                  arma::vec fused_responses = estimated_responses;
                  for(int r3 = 1; r3 <= num_col_indices; r3++) {
                    for(int sam = 0; sam < num_samples_responses; sam++) {
                      arma::vec estimated_r1 = estimated_responses( find( responses_indices_vec  == sam*max_indices_responses + r1*num_col_indices + r3 ) );
                      arma::vec estimated_r2 = estimated_responses( find( responses_indices_vec  == sam*max_indices_responses + r2*num_col_indices + r3 ) );
                      fused_responses( find( responses_indices_vec  == sam*max_indices_responses + r1*num_col_indices + r3 ) ) = (estimated_r1 + estimated_r2)/2;
                      indices_ones( find( responses_indices_vec  == sam*max_indices_responses + r2*num_col_indices + r3 ) ).fill(1);
                    }
                      fused_indices_responses.replace( r2*num_col_indices + r3 , r1*num_col_indices + r3 );
                  }
                  fused_responses.shed_rows( find( indices_ones == 1 ) );
                  arma::vec unique_non_zero_restrictions = unique( fused_indices_responses( find( fused_indices_responses != 0 ) ) );
                  int num_unique_non_zero_restrictions = unique_non_zero_restrictions.n_elem;
                  for (int i = 0; i < num_unique_non_zero_restrictions; i++) {
                    fused_indices_responses.replace( unique_non_zero_restrictions(i) , i+1 );
                  }
                  arma::mat pre_eval_mat = R_k(10,0);
                  R_temp = stratEst_EM(output_cube, sum_outputs_cube, strat_id, R_k(0,0), fused_responses, R_k(2,0), R_k(3,0), R_k(4,0), R_k(5,0), indices_shares, fused_indices_responses, indices_trembles, responses_to_sum, response, pre_eval_mat(0,0), outer_tol_eval, outer_max_eval, sample_of_ids_shares, num_ids_sample_shares, sample_of_ids_responses, num_ids_sample_responses, sample_of_ids_trembles, num_ids_sample_trembles );
                  arma::vec crit_vals = R_temp( 7 , 0 );
                  double crit_fused = crit_vals( crit_index );
                  if( std::isnan( crit_fused ) ){
                    crit_fused = arma::datum::inf;
                  }
                  if ( crit_fused < crit_min  ){
                    R_fused = R_temp;
                    crit_min = crit_fused;
                    new_indices_responses = fused_indices_responses;
                    trigger = 1;
                    triggered = true;
                  }
                }
              }
            }
          }
          if( triggered ){
            response_run = response_run + 1;
          }
          R_k = R_fused;
          indices_responses = new_indices_responses;
        }
      }  // end response fusion loop

      if( select_trembles && num_trembles_to_est > 1  ){

        if( print_messages == true ){
          Rcout << "\n";
          Rcout<< "select trembles             \r";
        }

        R_fused = R_k;
        arma::vec crit_vals = R_k( 7 , 0 );
        double crit_min = crit_vals( crit_index );
        arma::mat pre_eval_vec = R_k(10,0);
        num_it_trembles = R_k(2,0).n_elem;

        trigger = 1;
        int tremble_run = 1;

        // start response and tremble fusion loop
        while ( trigger == 1 ) {
          Rcpp::checkUserInterrupt();
          trigger = 0;
          triggered = false;

          // check every combination of trembles
          int max_indices_trembles = indices_trembles.max();
          arma::vec unique_non_zero_trembles = unique( indices_trembles( find( indices_trembles != 0 ) ) );
          arma::vec estimated_trembles = R_k(2,0);
          int num_estimated_trembles = estimated_trembles.n_elem;
          num_it_trembles = unique_non_zero_trembles.n_elem;
          arma::vec trembles_indices_vec = arma::linspace(1,num_estimated_trembles,num_estimated_trembles);
          int counter = 0;
          if( num_it_trembles >= 1 ){
            for (int t1 = 0; t1 < (num_it_trembles-1); t1++) {
              for (int t2 = t1+1; t2 < num_it_trembles; t2++) {
                counter = counter + 1;
                if( print_messages == true ){
                  Rcout<< "select trembles (" << tremble_run << "/" << counter << ")     \r";
                }
                fused_indices_trembles = indices_trembles;
                arma::vec indices_ones( num_estimated_trembles , arma::fill::zeros );
                arma::vec c1 = select_indices_trembles( find( indices_trembles == t1+1 ) );
                arma::vec c2 = select_indices_trembles( find( indices_trembles == t2+1 ) );
                if(  c1(0) == c2(0) ){
                  arma::vec fused_trembles = estimated_trembles;
                  for(int sam = 0; sam < num_samples_trembles; sam++) {
                    //fused_trembles(t1) = (estimated_trembles(sam*max_indices_trembles + t1) + estimated_trembles(sam*max_indices_trembles + t2))/2;
                    indices_ones( find( trembles_indices_vec == ( sam*max_indices_trembles + t2 + 1 ) ) ).fill(1);
                  }
                  fused_indices_trembles.replace( (t2+1) , (t1+1) );
                  fused_trembles.shed_rows( find( indices_ones == 1 ) );
                  fused_trembles.fill( arma::fill::randu );
                  fused_trembles /= 4;
                  arma::vec unique_non_zero_restrictions = unique( fused_indices_trembles( find( fused_indices_trembles != 0 ) ) );
                  int num_unique_non_zero_restrictions = unique_non_zero_restrictions.n_elem;
                  for (int i = 0; i < num_unique_non_zero_restrictions; i++) {
                    fused_indices_trembles.replace( unique_non_zero_restrictions(i) , i+1 );
                  }
                  arma::mat pre_eval_mat = R_k(12,0);
                  R_temp = stratEst_EM(output_cube, sum_outputs_cube, strat_id, R_k(0,0), R_k(1,0), fused_trembles, R_k(3,0), R_k(4,0), R_k(5,0), indices_shares, indices_responses, fused_indices_trembles, responses_to_sum, response, pre_eval_mat(0,0), outer_tol_eval, outer_max_eval, sample_of_ids_shares, num_ids_sample_shares, sample_of_ids_responses, num_ids_sample_responses, sample_of_ids_trembles, num_ids_sample_trembles );
                  arma::vec crit_vals = R_temp( 7 , 0 );
                  double crit_fused = crit_vals( crit_index );
                  if( std::isnan( crit_fused ) ){
                    crit_fused = arma::datum::inf;
                  }
                  if ( crit_fused < crit_min  ){
                    R_fused = R_temp;
                    crit_min = crit_fused;
                    new_indices_responses = indices_responses;
                    new_indices_trembles = fused_indices_trembles;
                    trigger = 1;
                    triggered = true;
                  }
                }
              }
            }
          }
          if( triggered ){
            tremble_run = tremble_run + 1;
          }
          R_k = R_fused;
          indices_trembles = new_indices_trembles;
        }
      }  // end tremble fusion loop

      //check if R_k is better than current best R
      arma::vec crit_vals = R_k( 7 , 0 );
      double crit_selected = crit_vals( crit_index );
      if( std::isnan( crit_selected ) ){
        crit_selected = arma::datum::inf;
      }
      if( print_messages == true && select_strategies ){
        Rcout<< "\n";
        Rcout<< crit << ": " << crit_selected;
        Rcout<< "\n";
      }
      if ( crit_selected < crit_selected_min || ( K == k && kill == 0 ) ){
        R = R_k;
        killed = kill;
        crit_selected_min = crit_selected;
        R_num_rows_strategies = num_rows_response_mat;
        R_num_strategies = num_candidates;
        R_output_cube( arma::span( 0 , num_rows_response_mat - 1 ) , arma::span::all , arma::span::all ) = output_cube;
        R_sum_outputs_cube( arma::span( 0 , num_rows_response_mat - 1 ) , arma::span::all , arma::span::all ) = sum_outputs_cube;
        R_strat_id( arma::span( 0 , num_rows_response_mat - 1 ) ) = strat_id;
        R_shares_mat( arma::span( 0 , num_candidates - 1 ) , arma::span::all ) = share_mat;
        R_indices_shares( arma::span( 0 , num_candidates - 1 ) , arma::span::all ) = indices_shares;
        R_indices_responses( arma::span( 0 , num_rows_response_mat - 1 ) , arma::span::all ) = indices_responses;
        R_indices_trembles( arma::span( 0 , num_rows_response_mat - 1 ) , arma::span::all ) = indices_trembles;
        R_responses_to_sum( arma::span( 0 , num_rows_response_mat - 1 ) , arma::span::all ) = responses_to_sum;
        if( LCR ){
          R_indices_coefficients( arma::span::all , arma::span( 0 , num_candidates - 1 ) ) = indices_coefficients;
          R_coefficient_mat( arma::span::all , arma::span( 0 , num_candidates - 1 ) ) = coefficients_mat;
        }
      }
      if( select_strategies == false || K == 1 ){
        kill = K+1;
      }
      if( integer_strategies == true && kill == 1  ){
        kill = K;
      }

    } // end kill loop


    if( killed == 0 || select_strategies == false ){
      K = 0;
      if( print_messages == true && select_strategies ){
        Rcout<< "\r";
        Rcout<< "end of strategy selection: no strategy can be eliminated";
      }
    }
    else if ( K <= 2 ){
      arma::vec survived = survivors( find( survivors != 0 ) );
      survivors( find( survivors == survived( killed - 1 ) ) ).fill(0);
      K = 0;
    }
    else{
      arma::vec survived = survivors( find( survivors != 0 ) );
      survivors( find( survivors == survived( killed - 1 ) ) ).fill(0);
      kill_start = 1;
      if( print_messages == true ){
        Rcout<< "\r";
        if( integer_strategies == true ){
          Rcout<< "eliminate one strategy (information criterion)";
        }
        else{
          Rcout<< "eliminate strategy " << survived( killed - 1 ) << " (information criterion)";
        }
      }
    }
    if( print_messages == true ){
      Rcout<< "\n";
    }
  } // end strategy selection loop

  // document which strategy survived the selection
  int num_elem_sid = sid.n_elem;
  int num_survivors = survivors.n_elem;
  arma::vec sid_zeros( num_elem_sid , arma::fill::zeros );
  for (int suv = 0; suv < num_survivors; suv++) {
    sid_zeros( find ( sid == survivors(suv) ) ).fill(1);
  }
  sid = sid( find( sid_zeros == 1 ) );

  //store fused objects
  arma::mat shares_new = R(3,0);
  k = shares_new.n_rows;
  arma::cube output_cube = R_output_cube( arma::span( 0 , R_num_rows_strategies - 1 ) , arma::span::all , arma::span::all );
  arma::cube sum_outputs_cube = R_sum_outputs_cube( arma::span( 0 , R_num_rows_strategies - 1 ) , arma::span::all , arma::span::all );
  arma::vec strat_id = R_strat_id( arma::span( 0 , R_num_rows_strategies - 1 ) );
  arma::umat shares_to_est = find_nonfinite( R_shares_mat( arma::span( 0 , R_num_strategies - 1 ) , arma::span::all ) );
  num_shares_to_est = shares_to_est.n_elem;
  arma::mat indices_shares = R_indices_shares( arma::span( 0 , k - 1 ) , arma::span::all );
  arma::mat indices_responses = R_indices_responses( arma::span( 0 , R_num_rows_strategies - 1 ) , arma::span::all );
  arma::mat indices_trembles = R_indices_trembles( arma::span( 0 , R_num_rows_strategies - 1 ) , arma::span::all );
  arma::mat responses_to_sum = R_responses_to_sum( arma::span( 0 , R_num_rows_strategies - 1 ) , arma::span::all );
  if( LCR ){
    indices_coefficients = R_indices_coefficients( arma::span::all , arma::span( 0 , k - 1 ) );
    coefficient_mat = R_coefficient_mat( arma::span::all , arma::span( 0 , k - 1 ) );
    num_coefficients_to_est = indices_coefficients.max();
  }
  arma::mat responses_new = R(1,0);
  arma::mat trembles_new = R(2,0);

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // latent class regression
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if( LCR ){

    arma::field<arma::mat> R_LCR(26,1);

    if( k == 1 ){
      Rcout<<"warning: only one strategy remains after selection. Cannot proceed with latent class regeression.\n";
      LCR = false;
    }
    else{
      bool replaced = false;
      if( print_messages == true ){
        Rcout<<"estimate lcr model\n";
      }
      //initialize empty coefficients
      arma::mat shares_no_lcr = R(3,0);
      arma::mat start_coefficients_mat  = reshape( indices_coefficients , num_coefficients_to_est , 1 );
      arma::uvec coefficients_to_est = find( start_coefficients_mat.col(0) > 0 );
      //arma::uvec coefficients_to_est = find_finite( start_coefficients );

      //arma::uvec coefficients_to_est = find( indices_coefficients > 0 );
      arma::mat R_shares = R(3,0);
      arma::vec start_intercepts = log( R_shares( arma::span( 1 , ( k-1 ) ) , 0 ) / R_shares(0,0) );
      arma::uvec shares_to_est_lcr = find( indices_shares > 0 );
      arma::mat pre_eval_vec = R(10,0);
      arma::field<arma::mat> R_LCR_inner(26,1);
      arma::vec LL_inner_min( 1 , 1 , arma::fill::zeros );
      LL_inner_min.fill( arma::datum::inf );
      // inner runs (index ri)
      for (int ri = 0; ri < LCR_runs; ri++) {
        Rcpp::checkUserInterrupt();
        if( print_messages == true ){
          Rcout<< "lcr run " << ri+1  << " (of " << LCR_runs << ")     \r";
        }

        //initialize coefficients to start
        arma::vec start_coefficients( num_coefficients_to_est , arma::fill::randn );
        start_coefficients = start_coefficients/10;

        // generate start intercepts
        arma::rowvec indices_coefficients_first_row = indices_coefficients.row(0);
        for (int i = 1; i < k; i++) {
          int target_index = indices_coefficients_first_row(i);
          if( target_index >= 1 ){
            start_coefficients( target_index - 1 ) = start_intercepts( i-1 );
          }
        }

        arma::mat pre_eval_vec = R(10,0);
        arma::field<arma::mat> R_temp = stratEst_LCR_EM( output_cube, sum_outputs_cube, strat_id, covariate_mat, R(0,0), R(1,0), R(2,0), start_coefficients, R(3,0), R(4,0), R(5,0), coefficient_mat, shares_to_est_lcr, indices_responses, indices_trembles, indices_coefficients, true, coefficients_to_est, responses_to_sum, response, pre_eval_vec(0,0), LCR_tol_eval, LCR_max_eval, newton_stepsize, penalty, sample_of_ids_shares, num_ids_sample_shares, sample_of_ids_responses, num_ids_sample_responses, sample_of_ids_trembles, num_ids_sample_trembles  );
        arma::mat LL_inner_temp = R_temp(LL_index,0);
        arma::mat R_temp_shares = R_temp(3,0);
        arma::mat R_temp_coefficients = R_temp(14,0);
        arma::mat R_temp_responses = R_temp(1,0);
        arma::mat R_temp_trembles = R_temp(2,0);
        if( LL_inner_temp.has_nan() || R_temp_shares.has_nan() || R_temp_coefficients.has_nan() || R_temp_responses.has_nan() || R_temp_trembles.has_nan() ){
          failed_inner_LCR_runs(0,0) += 1;
          LL_inner_temp.fill( arma::datum::inf );
        }
        if (  LL_inner_temp(0,0) > 0 && LL_inner_temp(0,0) < LL_inner_min(0,0) ){
          R_LCR = R_temp;
          LL_inner_min(0,0) = LL_inner_temp(0,0);
          replaced = true;
        }
      }
      if( replaced == false ){
        warning("Latent class regression failed. Likelihood values are infinite.");
        LCR = false;
      }
      else{
        R = R_LCR;
      }
      if( print_messages == true ){
        Rcout<< "\n";
      }
      // if( arma::max( vectorise(shares_no_lcr) - vectorise(R(3,0)) ) > 0.05 ){
      //   Rcout<<"warning: shares of latent class model deviate substantially from the model without covariates.\n";
      // }
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // analytical standard erros & convergence
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  arma::field<arma::mat> R_SE(17,1);
  R_SE = stratEst_SE( output_cube, sum_outputs_cube, strat_id, covariate_mat, R(3,0), R(1,0), R(2,0), R(14,0), R(4,0), R(5,0), R(15,0), R(13,0), R(19,0), R(17,0), R(18,0), R(16,0), indices_shares, indices_responses, indices_trembles, response, R(20,0), CL, cluster_id_vec, LCR, sample_of_ids_shares, num_ids_sample_shares, sample_of_ids_responses, num_ids_sample_responses, sample_of_ids_trembles, num_ids_sample_trembles );

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // bootstrapped standard errors
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // preallocate result matrices
  shares_to_est = find( indices_shares > 0 );
  num_shares_to_est = indices_shares.max();
  num_responses_to_est = R(1,0).n_elem;
  num_trembles_to_est = R(2,0).n_elem;
  num_coefficients_to_est = R(14,0).n_elem;
  arma::mat BS_shares_SE( num_shares_to_est , 1 , arma::fill::zeros );
  arma::mat BS_coefficients_SE( num_coefficients_to_est , 1 , arma::fill::zeros  );
  arma::mat BS_responses_SE( num_responses_to_est , 1 , arma::fill::zeros );
  arma::mat BS_trembles_SE( num_trembles_to_est , 1 , arma::fill::zeros  );

  arma::mat BS_shares_SE_mat( num_shares_to_est , BS_samples , arma::fill::zeros );
  arma::mat BS_coefficients_SE_mat( num_coefficients_to_est , BS_samples , arma::fill::zeros  );
  arma::mat BS_responses_SE_mat( num_responses_to_est , BS_samples , arma::fill::zeros );
  arma::mat BS_trembles_SE_mat( num_trembles_to_est , BS_samples , arma::fill::zeros  );
  arma::mat BS_check_mat( 1 , BS_samples , arma::fill::zeros );

  int num_quantiles = quantile_vec.n_elem;

  arma::mat BS_shares_quantiles( num_shares_to_est , num_quantiles , arma::fill::zeros );
  arma::mat BS_coefficients_quantiles( num_coefficients_to_est , num_quantiles , arma::fill::zeros  );
  arma::mat BS_responses_quantiles( num_responses_to_est , num_quantiles , arma::fill::zeros );
  arma::mat BS_trembles_quantiles( num_trembles_to_est , num_quantiles , arma::fill::zeros  );

  if( SE == "bootstrap" ){
    arma::mat estimated_shares = R(0,0);
    arma::mat estimated_coefficients = R(14,0);
    arma::mat estimated_responses = R(1,0);
    arma::mat estimated_trembles = R(2,0);
    if( response == "pure" ){
      num_responses_to_est = 0;
    }

    int BS_samples_shares = BS_samples;
    int BS_samples_responses = BS_samples;
    int BS_samples_trembles = BS_samples;
    int BS_samples_coefficients = BS_samples;
    if( print_messages == true ){
      Rcout<<"start bootstrap\n";
    }
    for (int i = 0; i < BS_samples; i++) {
      Rcpp::checkUserInterrupt();
      if( print_messages == true ){
        Rcout<< "sample " << i+1 << " (of " << BS_samples <<")\r";
      }
      arma::field<arma::mat> R_LCR_BS(26,1);
      arma::field<arma::mat> R_NO_LCR_BS(26,1);
      arma::field<arma::mat> R_CHECK_BS(26,1);

      int num_ids_to_sample = num_ids;

      // sample clusters if CL true
      arma::vec sampled_clusters( num_clusters , arma::fill::randu );
      if( CL ){
        sampled_clusters *= num_clusters;
        sampled_clusters = ceil( sampled_clusters );
        int ids_in_sampled_clusters = 0;
        for (int j = 0; j < num_clusters; j++) {
          arma::vec unique_ids_in_sampled_cluster = unique( id_vec( find( cluster == sampled_clusters(j) ) ) );
          ids_in_sampled_clusters = ids_in_sampled_clusters + unique_ids_in_sampled_cluster.n_elem;
        }
        num_ids_to_sample = ids_in_sampled_clusters;
      }
      arma::vec sampled_ids( num_ids_to_sample , arma::fill::zeros );

      // sample ids
      if( CL ){
        arma::mat sampled_ids_mat = unique( id_vec( find( cluster == sampled_clusters(0) ) ) );
        for (int j = 1; j < num_clusters; j++) {
          arma::mat sampled_ids_mat_old = sampled_ids_mat;
          arma::mat sampled_ids_mat_new = unique( id_vec( find( cluster == sampled_clusters(j) ) ) );
          sampled_ids_mat = join_cols(sampled_ids_mat_old,sampled_ids_mat_new);
        }
        sampled_ids = sampled_ids_mat.col(0);
      }
      else{
        sampled_ids.randu();
        sampled_ids *= num_ids_to_sample;
        sampled_ids = ceil( sampled_ids );
      }

      int num_sampled_ids = sampled_ids.n_elem;
      arma::cube output_cube_BS_sample( output_cube.n_rows , output_cube.n_cols , num_sampled_ids , arma::fill::zeros );
      arma::cube sum_outputs_cube_BS_sample( sum_outputs_cube.n_rows , sum_outputs_cube.n_cols , num_sampled_ids , arma::fill::zeros );
      arma::mat covariate_mat_BS_sample( num_sampled_ids , covariate_mat.n_cols , arma::fill::zeros );
      arma::vec sample_of_ids_shares_BS( num_sampled_ids , arma::fill::zeros );
      arma::vec num_ids_sample_shares_BS( num_samples , arma::fill::zeros );
      arma::vec sample_of_ids_responses_BS( num_sampled_ids , arma::fill::zeros );
      arma::vec num_ids_sample_responses_BS( num_samples , arma::fill::zeros );
      arma::vec sample_of_ids_trembles_BS( num_sampled_ids , arma::fill::zeros );
      arma::vec num_ids_sample_trembles_BS( num_samples , arma::fill::zeros );
      for (int j = 0; j < num_sampled_ids; j++) {
        int sampled_ids_j = sampled_ids(j);
        output_cube_BS_sample.slice(j) = output_cube.slice( sampled_ids_j - 1 );
        sum_outputs_cube_BS_sample.slice(j) = sum_outputs_cube.slice( sampled_ids_j - 1 );
        arma::vec unique_sample_of_ids_shares_BS = unique( sample_id( find( id_vec == sampled_ids_j ) ) );

        int sample_of_ids_shares_j = sample_of_ids_shares( sampled_ids_j - 1 );
        sample_of_ids_shares_BS(j) =  sample_of_ids_shares_j;
        num_ids_sample_shares_BS( sample_of_ids_shares_j - 1 ) += 1;

        int sample_of_ids_responses_j = sample_of_ids_responses( sampled_ids_j - 1 );
        sample_of_ids_responses_BS(j) = sample_of_ids_responses_j;
        num_ids_sample_responses_BS( sample_of_ids_responses_j - 1 ) += 1;

        int sample_of_ids_trembles_j = sample_of_ids_trembles( sampled_ids_j - 1 );
        sample_of_ids_trembles_BS(j) = sample_of_ids_trembles_j;
        num_ids_sample_trembles_BS( sample_of_ids_trembles_j - 1 ) += 1;

        if( LCR ){
          covariate_mat_BS_sample.row(j) = covariate_mat.row( sampled_ids_j - 1 );
        }
      }

        arma::mat pre_eval_mat( 1 , 1 , arma::fill::zeros );
        arma::mat indices_shares_zero( indices_shares.n_rows , indices_shares.n_cols , arma::fill::zeros );
        arma::mat indices_responses_zero( indices_responses.n_rows , indices_responses.n_cols , arma::fill::zeros );
        arma::mat indices_trembles_zero( indices_trembles.n_rows , indices_trembles.n_cols , arma::fill::zeros );
        arma::mat indices_coefficients_zero( indices_coefficients.n_rows , indices_coefficients.n_cols , arma::fill::zeros );

        arma::vec empty_shares( 0 , arma::fill::none );
        arma::vec empty_responses( 0 , arma::fill::none );
        arma::vec empty_trembles( 0 , arma::fill::none );

        // bootstrap shares or coefficients
        if( num_shares_to_est > 0 ){
          arma::mat estimated_shares_BS( k , num_samples , arma::fill::zeros );
          if( LCR == false ){
            R_NO_LCR_BS = stratEst_EM( output_cube_BS_sample, sum_outputs_cube_BS_sample, strat_id, R(0,0), empty_responses, empty_trembles, R(3,0), R(4,0), R(5,0), indices_shares, indices_responses_zero, indices_trembles_zero, responses_to_sum, response, 0, outer_tol_eval, outer_max_eval, sample_of_ids_shares_BS, num_ids_sample_shares_BS, sample_of_ids_responses_BS, num_ids_sample_responses_BS, sample_of_ids_trembles_BS, num_ids_sample_trembles_BS );
            estimated_shares_BS = R_NO_LCR_BS(0,0);
            if( estimated_shares_BS.is_finite() && estimated_shares_BS.n_elem > 0 ){
              BS_shares_SE_mat( arma::span::all , i ) = estimated_shares_BS;
              BS_shares_SE += ( ( estimated_shares_BS - estimated_shares) % (estimated_shares_BS - estimated_shares ) );
            }
            else{
              BS_samples_shares = BS_samples_shares - 1;
            }
          }
          else{
            arma::mat estimated_coefficients_BS( covariate_mat.n_cols * (k-1) , 1 , arma::fill::zeros );
            arma::uvec coefficients_to_est_bs = find_finite( R(14,0) );
            R_LCR_BS = stratEst_LCR_EM( output_cube_BS_sample, sum_outputs_cube_BS_sample, strat_id, covariate_mat_BS_sample, R(0,0), empty_responses, empty_trembles, R(14,0), R(3,0), R(4,0), R(5,0), R(15,0), shares_to_est, indices_responses_zero, indices_trembles_zero, indices_coefficients, true, coefficients_to_est_bs, responses_to_sum, response, 0 , outer_tol_eval, outer_max_eval, newton_stepsize, true, sample_of_ids_shares_BS, num_ids_sample_shares_BS, sample_of_ids_responses_BS, num_ids_sample_responses_BS, sample_of_ids_trembles_BS, num_ids_sample_trembles_BS );
            estimated_coefficients_BS = R_LCR_BS(14,0);
            estimated_shares_BS = R_LCR_BS(3,0);
            if( estimated_shares_BS.is_finite() && estimated_coefficients_BS.is_finite() ){
              BS_shares_SE_mat( arma::span::all , i ) = estimated_shares_BS( shares_to_est );
              BS_coefficients_SE_mat( arma::span::all , i ) = estimated_coefficients_BS( coefficients_to_est_bs );
              BS_shares_SE += ( ( estimated_shares_BS( shares_to_est ) - estimated_shares( shares_to_est )) % (estimated_shares_BS( shares_to_est ) - estimated_shares( shares_to_est ) ) );
              BS_coefficients_SE += ( ( estimated_coefficients_BS( coefficients_to_est_bs ) - estimated_coefficients( coefficients_to_est_bs )) % (estimated_coefficients_BS( coefficients_to_est_bs ) - estimated_coefficients( coefficients_to_est_bs ) ) );
            }
            else{
              BS_samples_coefficients = BS_samples_coefficients - 1;
            }
          }

        }

        // bootstrap responses & trembles
        if( response != "pure" && num_responses_to_est > 0 ){
          arma::vec estimated_responses = R(1,0);
          arma::vec estimated_responses_BS( estimated_responses.n_elem , arma::fill::zeros );
          arma::mat indices_responses_BS = indices_responses;
          indices_responses_BS.fill(0);
          int num_estimated_responses = estimated_responses.n_elem;
          int num_col_indices = indices_responses.n_cols;
          int max_indices_responses = indices_responses.max();
          int num_it_responses = max_indices_responses/num_col_indices;
          arma::vec responses_indices_vec = arma::linspace(1,num_estimated_responses,num_estimated_responses);
          arma::vec responses_BS( num_col_indices*num_samples_responses , arma::fill::zeros );
          for (int r1 = 0; r1 < num_it_responses; r1++) {
            arma::vec indices_ones( num_estimated_responses , arma::fill::zeros );
            for(int r3 = 1; r3 <= num_col_indices; r3++) {
              for(int sam = 0; sam < num_samples_responses; sam++) {
                indices_ones( find( responses_indices_vec  == sam*max_indices_responses + r1*num_col_indices + r3 ) ).fill(1);
              }
              indices_responses_BS( find( indices_responses == r1*num_col_indices + r3 ) ).fill(r3);

            }
            responses_BS = estimated_responses( find( indices_ones == 1 ) );
            if( LCR ){
              arma::uvec coefficients_to_est_bs = find_finite( R(14,0) );
              R_LCR_BS = stratEst_LCR_EM( output_cube_BS_sample, sum_outputs_cube_BS_sample, strat_id, covariate_mat_BS_sample, empty_shares, responses_BS, empty_trembles, R(14,0), R(3,0), R(4,0), R(5,0), R(15,0), shares_to_est, indices_responses_BS, indices_trembles_zero, indices_coefficients_zero, false, coefficients_to_est_bs, responses_to_sum, response, 0 , outer_tol_eval, outer_max_eval, newton_stepsize, penalty, sample_of_ids_shares_BS, num_ids_sample_shares_BS, sample_of_ids_responses_BS, num_ids_sample_responses_BS, sample_of_ids_trembles_BS, num_ids_sample_trembles_BS );
              arma::mat response_estimate_BS = R_LCR_BS(1,0);
              estimated_responses_BS( find( indices_ones == 1 ) ) = response_estimate_BS;
            }else{
              R_NO_LCR_BS = stratEst_EM( output_cube_BS_sample, sum_outputs_cube_BS_sample, strat_id, empty_shares, responses_BS, empty_trembles, R(3,0), R(4,0), R(5,0), indices_shares_zero, indices_responses_BS, indices_trembles_zero, responses_to_sum, response, 0, outer_tol_eval, outer_max_eval, sample_of_ids_shares_BS, num_ids_sample_shares_BS, sample_of_ids_responses_BS, num_ids_sample_responses_BS, sample_of_ids_trembles_BS, num_ids_sample_trembles_BS );
              arma::mat response_estimate_BS = R_NO_LCR_BS(1,0);
              estimated_responses_BS( find( indices_ones == 1 ) ) = response_estimate_BS;
            }
          }
          if( estimated_responses_BS.is_finite() ){
            BS_responses_SE_mat( arma::span::all , i ) = estimated_responses_BS;
            BS_responses_SE += ( (estimated_responses_BS - estimated_responses) % (estimated_responses_BS - estimated_responses ) );
          }
          else{
            BS_samples_responses = BS_samples_responses - 1;
          }
        }

        int num_trembles_per_sample = indices_trembles.max();
        arma::mat estimated_trembles_BS( estimated_trembles.n_elem , 1 , arma::fill::zeros );
        arma::mat trembles_BS( num_samples_trembles , 1 );
        for (int j = 0; j < num_trembles_per_sample; j++) {
          for( int k = 0; k < num_samples_trembles; k++){
            trembles_BS(k,0) = ( estimated_trembles(j+k*num_trembles_per_sample,0) );
          }
          arma::mat indices_trembles_BS = indices_trembles_zero;
          indices_trembles_BS( find( indices_trembles == j+1 ) ).fill(1);
          if( LCR ){
            arma::uvec coefficients_to_est = find_finite( R(14,0) );
            R_LCR_BS = stratEst_LCR_EM( output_cube_BS_sample, sum_outputs_cube_BS_sample, strat_id, covariate_mat_BS_sample, empty_shares, empty_responses, trembles_BS, R(14,0), R(3,0), R(4,0), R(5,0), R(15,0), shares_to_est, indices_responses_zero, indices_trembles_BS, indices_coefficients_zero, false, coefficients_to_est, responses_to_sum, response, 0, outer_tol_eval, outer_max_eval, newton_stepsize, penalty, sample_of_ids_shares_BS, num_ids_sample_shares_BS, sample_of_ids_responses_BS, num_ids_sample_responses_BS, sample_of_ids_trembles_BS, num_ids_sample_trembles_BS );
            arma::mat trembles_estimate_BS = R_LCR_BS(2,0);
            for( int k = 0; k < num_samples_trembles; k++){
              estimated_trembles_BS(j+k*num_trembles_per_sample,0) = trembles_estimate_BS(k,0);
            }
          }
          else{
            R_NO_LCR_BS = stratEst_EM( output_cube_BS_sample, sum_outputs_cube_BS_sample, strat_id, empty_shares, empty_responses, trembles_BS, R(3,0), R(4,0), R(5,0), indices_shares_zero, indices_responses_zero, indices_trembles_BS, responses_to_sum, response, 0 , inner_tol_eval, inner_max_eval, sample_of_ids_shares_BS, num_ids_sample_shares_BS, sample_of_ids_responses_BS, num_ids_sample_responses_BS, sample_of_ids_trembles_BS, num_ids_sample_trembles_BS );
            arma::mat trembles_estimate_BS = R_NO_LCR_BS(2,0);
            for( int k = 0; k < num_samples_trembles; k++){
              estimated_trembles_BS(j+k*num_trembles_per_sample,0) = trembles_estimate_BS(k,0);
            }
          }
        }
        if( estimated_trembles_BS.is_finite() ){
          BS_trembles_SE_mat( arma::span::all , i ) = estimated_trembles_BS;
          BS_trembles_SE += ( (estimated_trembles_BS - estimated_trembles) % (estimated_trembles_BS - estimated_trembles ) );
        }
        else{
          BS_samples_trembles = BS_samples_trembles - 1;
        }

    }
    // store bootstrap results
    if( num_shares_to_est > 0 ){
      if( LCR ){
        BS_shares_SE = sqrt( BS_shares_SE/BS_samples_coefficients );
        BS_coefficients_SE = sqrt( BS_coefficients_SE/BS_samples_coefficients );
        BS_shares_quantiles = arma::quantile( BS_shares_SE_mat , quantile_vec , 1  );
        BS_coefficients_quantiles = arma::quantile( BS_coefficients_SE_mat , quantile_vec , 1  );
        if( BS_samples_coefficients/BS_samples < 0.9 ){
          BS_shares_SE.fill(-1);
          BS_coefficients_SE.fill(-1);
          BS_shares_quantiles.fill(-1);
          BS_coefficients_quantiles.fill(-1);
          Rcout<<"stratEst warning: Bootstrap for coefficients failed. More than 10 percent of the bootstrap samples did not produce estimates.\n";
        }
      }
      else{
        BS_shares_SE = sqrt( BS_shares_SE/BS_samples_shares );
        BS_shares_quantiles = arma::quantile( BS_shares_SE_mat , quantile_vec , 1  );
        if( BS_samples_shares/BS_samples < 0.9 ){
          BS_shares_SE.fill(-1);
          BS_shares_quantiles.fill(-1);
          Rcout<<"stratEst warning: Bootstrap for shares failed. More than 10 percent of the bootstrap samples did not produce estimates.\n";
        }
      }
    }
    if( num_responses_to_est > 0 ){
      BS_responses_SE = sqrt( BS_responses_SE/BS_samples_responses );
      BS_responses_quantiles = arma::quantile( BS_responses_SE_mat , quantile_vec , 1  );
      BS_trembles_SE = sqrt( BS_trembles_SE/BS_samples_trembles );
      BS_trembles_quantiles = arma::quantile( BS_trembles_SE_mat , quantile_vec , 1  );
      if( BS_samples_responses/BS_samples < 0.9 ){
        BS_responses_SE.fill(-1);
        BS_responses_quantiles.fill(-1);
        Rcout<<"stratEst warning: Bootstrap for responses failed. More than 10 percent of the bootstrap samples did not produce estimates.\n";
      }
    }
    if( num_trembles_to_est > 0 ){
      if( BS_samples_trembles/BS_samples < 0.9 ){
        BS_trembles_SE.fill(-1);
        BS_trembles_quantiles.fill(-1);
        Rcout<<"stratEst warning: Bootstrap for trembles failed. More than 10 percent of the bootstrap samples did not produce estimates.\n";
      }
    }
    if( print_messages == true ){
      Rcout<< "\n";
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // output
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  num_responses_to_est = R(1,0).n_elem;

  // rbind trembles_mat
  arma::mat new_trembles = R(2,0);
  arma::mat tremble_mat = R(5,0);
  int indices_trembles_max = indices_trembles.max();
  int it_rows = indices_trembles.n_rows;
  int it_cols = indices_trembles.n_cols;
  arma::mat indices_trembles_long( it_rows*num_samples_trembles , it_cols  , arma::fill::zeros );
  arma::mat trembles_long( it_rows*num_samples_trembles , it_cols  , arma::fill::zeros );
  for (int i = 0; i < num_samples_trembles; i++) {
    arma::mat indices_trembles_mat = indices_trembles;//
    indices_trembles_mat += i*(indices_trembles_max);
    arma::mat trembles_mat = tremble_mat;
    for( int j = 0; j < num_trembles_to_est; j++) {
      trembles_mat( find( indices_trembles_mat == j+1 ) ).fill( new_trembles(j) );
    }
    trembles_long( arma::span( it_rows*i  , it_rows*i + it_rows - 1 ) , arma::span::all  ) = trembles_mat;
    indices_trembles_long( arma::span( it_rows*i , it_rows*i + it_rows - 1 ) , arma::span::all  ) = indices_trembles_mat;
  }

  // update strategy parameters
  arma::mat strategy_responses = R(4,0);
  arma::vec finally_survived = survivors( find( survivors != 0 ) );
  int num_finally_survived = finally_survived.n_elem;
  arma::vec strat_id_survived( rows_strategies , arma::fill::zeros );
  for (int i = 0; i < num_finally_survived; i++) {
    strat_id_survived( find( complete_strat_id == finally_survived(i) ) ).fill(1);
  }
  arma::mat final_strategies( strategy_responses.n_rows , cols_strategies , arma::fill::zeros  );
  arma::mat strat_id_survived_mat = repmat( strat_id_survived , 1 , cols_strategies );
  arma::vec strategies_alive = strategies( find( strat_id_survived_mat == 1 ) );
  final_strategies = reshape( strategies_alive , strategy_responses.n_rows , cols_strategies );
  final_strategies.cols( 1 , num_non_zero_outputs ) = strategy_responses.cols( num_outputs - num_non_zero_outputs , num_outputs - 1 );

  // create objects for samples
  arma::mat response_long_mat = R(4,0);
  arma::mat tremble_long_mat = R(5,0);
  arma::mat indices_responses_long_mat = indices_responses;
  arma::mat indices_trembles_long_mat = indices_trembles;

  if( ( num_samples_trembles > 1 && num_trembles_to_est > 0 ) || ( num_samples_responses > 1 && num_responses_to_est > 0 ) ){

    final_strategies = repmat( final_strategies , num_samples , 1 );
    response_long_mat = repmat( response_long_mat , num_samples , 1 );
    tremble_long_mat = repmat( tremble_long_mat , num_samples , 1 );

    if( num_samples_responses > 1 ){
      // create indices responses long mat from mat
      int indices_responses_max = indices_responses.max();
      for (int i = 1; i < num_samples_responses; i++) {
        arma::mat next_indices = indices_responses;
        next_indices(find(next_indices > 0)) += i*(indices_responses_max);
        indices_responses_long_mat = join_cols( indices_responses_long_mat , next_indices );
      }
    }
    else{
      indices_responses_long_mat = repmat( indices_responses , num_samples , 1 );
    }

    if( num_samples_trembles > 1 ){
      // create indices trembles long mat from mat
      int indices_trembles_max = indices_trembles.max();
      for (int i = 1; i < num_samples_trembles; i++) {
        arma::mat next_indices = indices_trembles;
        next_indices(find(next_indices > 0)) += i*(indices_trembles_max);
        indices_trembles_long_mat = join_cols( indices_trembles_long_mat , next_indices );
      }
    }
    else{
      indices_trembles_long_mat = repmat( indices_trembles , num_samples , 1 );
    }

    // fill long matrices
    arma::mat responses_par = R(1,0);
    arma::mat trembles_par = R(2,0);
    arma::vec responses_par_vec = responses_par.col(0);
    arma::vec trembles_par_vec = trembles_par.col(0);

    if( num_samples_responses > 1 ){
      for (int i = 0; i < num_responses_to_est ; i++) {
        response_long_mat( find( indices_responses_long_mat == i+1 ) ).fill( responses_par_vec(i) );
      }
    }
    if( num_samples_trembles > 1 ){
      for (int i = 0; i < num_trembles_to_est ; i++) {
        tremble_long_mat( find( indices_trembles_long_mat == i+1 ) ).fill( trembles_par_vec(i) );
      }
    }

    // update strategies
    final_strategies.cols( 1 , num_non_zero_outputs ) = response_long_mat.cols( num_outputs - num_non_zero_outputs , num_outputs - 1 );

  }

  // global chi square
  arma::mat chi_mat( 1 , 2 , arma::fill::zeros );
  arma::mat observed = R(21,0);
  arma::mat probabilities_mat = response_long_mat % (1 - tremble_long_mat) + ( 1 - response_long_mat ) % ( tremble_long_mat / (num_outputs - 1) );
  arma::mat sum_observed_mat = sum( observed , 1 );
  arma::mat expected = probabilities_mat;
  expected.each_col() %= sum_observed_mat;
  arma::umat non_zero = find( expected != 0 );
  int num_non_zero = non_zero.n_elem;
  double chi = 0;
  if( num_non_zero > 0 ){
    chi = accu(((observed(non_zero) - expected(non_zero)) % (observed(non_zero) - expected(non_zero)))/expected(non_zero) );
  }
  if( arma::is_finite(chi) ){
    chi_mat(0,0) = chi;
  }

  //local chi square
  arma::mat chi_local( 1 , k , arma::fill::zeros );
  arma::mat strat_id_mat = repmat( strat_id , 1 , k );
  arma::mat assignments = R(13.0);
  for (int i = 0; i < num_ids; i++) {
    arma::rowvec target_row = assignments.row(i);
    arma::uword index_id = index_max( target_row );
    arma::vec v = arma::linspace(1,k,k);
    int k_id = v(index_id);
    arma::mat observed_k = output_cube.slice(i);
    arma::mat sum_observed_k = sum_outputs_cube.slice(i);
    arma::mat observed_id = observed_k.rows( find( strat_id == k_id ) );
    arma::mat sum_observed_id = sum_observed_k.rows( find( strat_id == k_id ) );
    arma::mat expected_id = probabilities_mat.rows( find( strat_id == k_id ) );
    expected_id %= sum_observed_id;
    arma::umat non_zero_id = find( expected_id != 0 );
    int num_non_zero_id = non_zero_id.n_elem;
    double chi_id = 0;
    if( num_non_zero_id > 0 ){
      chi_id = accu(((observed_id(non_zero_id) - expected_id(non_zero_id)) % (observed_id(non_zero_id) - expected_id(non_zero_id)))/expected_id(non_zero_id) );
    }
    if( k_id > 0 && arma::is_finite(chi_id) ){
      chi_local(k_id-1) += chi_id;
    }
  }

  arma::mat SE_shares = R_SE(0,0);
  arma::mat covar_shares = R_SE(5,0);
  arma::mat score_shares = R_SE(6,0);
  arma::mat fisher_shares = R_SE(7,0);

  //prepare output
  arma::mat final_shares_vec = R(3,0);
  arma::mat final_shares = R(0,0);
  arma::mat final_responses = R(1,0);
  arma::mat final_trembles = R(2,0);
  arma::mat final_indices_shares = indices_shares;
  arma::mat final_indices_responses = indices_responses_long_mat;
  arma::mat final_indices_trembles = indices_trembles_long_mat.col(0);
  arma::mat final_tremble_vec = tremble_long_mat.col(0);
  arma::mat final_response_mat = response_long_mat;
  arma::mat final_LL =  R(6,0);
  arma::mat final_crit = R(7,0);
  arma::mat final_eval = R(10,0);
  arma::mat final_eps = R(11,0);
  arma::mat final_entropy = R(12,0);
  arma::mat final_i_class = R(13,0);
  arma::mat final_coefficients = R(14,0);
  arma::mat final_coefficient_mat = R(15,0);
  arma::mat final_individual_priors = R(16,0);
  arma::mat final_state_obs = R(21,0);
  arma::mat final_free_params = R(25,0);

  //final stats
  arma::mat final_SE_responses = R_SE(1,0);
  arma::mat final_SE_trembles = R_SE(2,0);
  arma::mat SE_coefficients = R_SE(3,0);
  arma::mat final_convergence = R_SE(4,0);
  arma::mat final_covar_responses = R_SE(8,0);
  arma::mat final_score_responses = R_SE(9,0);
  arma::mat final_fisher_responses = R_SE(10,0);
  arma::mat final_covar_trembles = R_SE(11,0);
  arma::mat final_score_trembles = R_SE(12,0);
  arma::mat final_fisher_trembles = R_SE(13,0);
  arma::mat covar_coefficients = R_SE(14,0);
  arma::mat score_coefficients = R_SE(15,0);
  arma::mat fisher_coefficients = R_SE(16,0);

  arma::mat final_quantiles_shares = BS_shares_quantiles;
  arma::mat final_quantiles_coefficients = BS_coefficients_quantiles;
  arma::mat final_quantiles_responses = BS_responses_quantiles;
  arma::mat final_quantiles_trembles = BS_trembles_quantiles;

  // eliminate fixed shares.par elements
  if( specific_shares == 0 && num_shares_to_est > 0 ){
    arma::mat zero_indices_shares = final_indices_shares;
    zero_indices_shares.fill(0);
    zero_indices_shares.col(0) = final_indices_shares.col(0);
    final_indices_shares = zero_indices_shares;
  }

  final_shares = final_shares_vec( find( final_indices_shares > 0 ) );
  int num_final_shares = final_shares.n_elem;
  arma::mat final_score_shares = score_shares( find( final_indices_shares > 0 ) );
  arma::vec indices_shares_vec = vectorise( final_indices_shares );
  arma::mat final_SE_shares = SE_shares( find( indices_shares_vec > 0 ) );
  arma::mat shares_indices_mat = indices_shares_vec * indices_shares_vec.t() ;
  arma::mat final_fisher_shares = fisher_shares( find( shares_indices_mat > 0 ) );
  arma::mat final_covar_shares = covar_shares( find( shares_indices_mat > 0 ) ) ;
  final_fisher_shares.reshape( num_final_shares , num_final_shares );
  final_covar_shares.reshape( num_final_shares , num_final_shares );
  if( num_final_shares > 0 ){
    final_quantiles_shares = final_quantiles_shares( arma::span(0,num_final_shares-1) , arma::span::all );
  }

  if( specific_shares == 0 && num_shares_to_est > 0 ){
    final_indices_shares = final_indices_shares.col(0);
  }

  arma::mat final_SE_coefficients = SE_coefficients;
  arma::mat final_covar_coefficients = covar_coefficients;
  arma::mat final_score_coefficients = score_coefficients;
  arma::mat final_fisher_coefficients = fisher_coefficients;

  if( LCR ){
    // eliminate fixed coefficients.par elements
    arma::vec indices_coefficients_short_vec = vectorise( indices_coefficients( arma::span::all , arma::span( 0 , indices_coefficients.n_cols - 1 ) ) );
    arma::uvec positive_coefficients = find( indices_coefficients_short_vec > 0 );
    int num_positive_coefficients = positive_coefficients.n_elem;
    indices_coefficients_short_vec = indices_coefficients_short_vec( find( indices_coefficients_short_vec > 0 ) );
    arma::mat indices_coefficients_short_vec_mat = indices_coefficients_short_vec * indices_coefficients_short_vec.t();
    final_SE_coefficients = SE_coefficients( find( indices_coefficients_short_vec > 0 ) );
    final_covar_coefficients = covar_coefficients( find( indices_coefficients_short_vec_mat > 0 ) );
    final_covar_coefficients.reshape( num_positive_coefficients , num_positive_coefficients );
    final_score_coefficients = score_coefficients( find( indices_coefficients_short_vec > 0 ) );
    final_fisher_coefficients = fisher_coefficients( find( indices_coefficients_short_vec_mat > 0 ) ) ;
    final_fisher_coefficients.reshape( num_positive_coefficients , num_positive_coefficients );
  }

  //final fit variables
  arma::mat final_fit( 1 , 7+k , arma::fill::zeros );
  final_fit(0,0) = -final_LL(0,0);
  final_fit(0,1) = final_crit( crit_index );
  final_fit(0,2) = final_entropy(0,0);
  final_fit(0,arma::span(3,5) ) = final_crit.t();
  final_fit(0,6) = chi_mat(0,0);
  final_fit(0,arma::span(7,7+k-1)) = chi_local;
  arma::mat final_solver( 1 , 3 , arma::fill::zeros );
  final_solver(0,0) = final_eval(0,0);
  final_solver(0,1) = final_eps(0,0);
  if( SE == "bootstrap" ){
    if( LCR ){
      final_SE_shares = BS_shares_SE;
      final_SE_coefficients = BS_coefficients_SE;
    }
    else{
      if( num_shares_to_est > 0 ){
        final_SE_shares = BS_shares_SE; //BS_shares_SE( find( indices_shares_vec > 0 ) );
      }
    }
    if( num_responses_to_est > 0 && response == "mixed"){
      final_SE_responses = BS_responses_SE;
    }
    if( num_trembles_to_est > 0 ){
      final_SE_trembles = BS_trembles_SE;
    }
  }

  List shares_list = List::create( Named("shares.par") = final_shares, Named("shares.indices") = final_indices_shares, Named("shares.se") = final_SE_shares, Named("shares.score") =  final_score_shares, Named("shares.covar") = final_covar_shares, Named("shares.fisher") = final_fisher_shares , Named("shares.quantiles") = final_quantiles_shares );
  List responses_list = List::create( Named("responses.par") = final_responses, Named("responses.indices") = final_indices_responses, Named("responses.se") = final_SE_responses, Named("responses.covar") = final_covar_responses, Named("responses.score") = final_score_responses, Named("responses.fisher") = final_fisher_responses, Named("responses.quantiles") = final_quantiles_responses  );
  List trembles_list = List::create( Named("trembles.par") = final_trembles, Named("trembles.indices") = final_indices_trembles, Named("trembles.se") = final_SE_trembles, Named("trembles.covar") = final_covar_trembles, Named("trembles.score") = final_score_trembles, Named("trembles.fisher") = final_fisher_trembles, Named("trembles.quantiles") = final_quantiles_trembles );
  List coefficients_list = List::create( Named("coefficients.par") = final_coefficients,  Named("coefficients.se") = final_SE_coefficients, Named("coefficients.covar") = final_covar_coefficients, Named("coefficients.score") = final_score_coefficients, Named("coefficients.fisher") = final_fisher_coefficients, Named("coefficients.quantiles") = final_quantiles_coefficients, Named("covariate.mat") = covariate_mat  );

  List S = List::create( Named("shares") = final_shares_vec, Named("strategies") = final_strategies,  Named("responses") = final_response_mat, Named("trembles") = final_tremble_vec,  Named("coefficients") = final_coefficient_mat, Named("fit") = final_fit, Named("solver") = final_solver, Named("state.obs") = final_state_obs, Named("assignments") = final_i_class, Named("priors") = final_individual_priors, Named("convergence") = final_convergence, Named("sid") = sid, Named("selected.strats") = survivors, Named("lcr") = LCR , Named("n.par") = final_free_params , Named("shares.list") = shares_list, Named("responses.list") = responses_list, Named("trembles.list") = trembles_list, Named("coefficients.list") = coefficients_list );

  if( print_messages == true ){
    Rcout<< "---\n";
  }

  return(S);
}






