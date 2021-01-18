/*
 * Author: S. Weber
 *
 * Various utlities
 */

// turn a slicing variable for a ragged array
// S = {5, 6, 11}
// into
// Si = {0, 5, 5 + 6, 5+6+11} + 1
// such that we can index the ragged array A as
// A[Si[u] : Si[u+1]-1]
// for the uth unit
int[] make_slice_index(int[] S) {
  int Si[size(S)+1];
  int cv = 1;
  Si[1] = cv;
  for(i in 1:size(S)) {
    cv = cv + S[i];
    Si[i+1] = cv;
  }
  return(Si);
}

// helper function for rle_int
int rle_elem_count(int[] set) {
  int U = 1;
  for(i in 2:num_elements(set)) {
    if(set[i-1] != set[i])
      U = U + 1;
  }
  return(U);
}

// helper function for rle_int
int rle_elem_count_vector(vector set) {
  int U = 1;
  for(i in 2:num_elements(set)) {
    if(set[i-1] != set[i])
      U = U + 1;
  }
  return(U);
}

// repeated length encoding, see rle in R
// this function returns the # of times elements are repeated ...
int[] rle_int(int[] set) {
  int res[rle_elem_count(set)];
  int c = 1;
  res[1] = 1;
  for(i in 2:num_elements(set)) {
    if(set[i-1] == set[i]) {
      res[c] = res[c] + 1;
    } else {
      c = c + 1;
      res[c] = 1;
    }
  }
  return(res);
}

// ... and this function returns which elements are repeated for ints
int[] rle_elem_int(int[] set) {
  int N = rle_elem_count(set);
  int first_ind[N] = make_slice_index(rle_int(set))[1:N];

  return set[first_ind];
}

/* checks if in an id vector a given id is used in multiple blocks. If
   so, a reject is fired to ensure that the id vector is blocked into
   sorted groups.
 */
void check_duplicate_ids(int[] id) {
  int N = rle_elem_count(id);
  int sorted_ids[N] = sort_asc(rle_elem_int(id));
  int cid = sorted_ids[1];
  for(i in 1:N-1) {
    if(sorted_ids[i] == sorted_ids[i+1])
      reject("ID ", sorted_ids[i], " occurs multiple times within id vector.");
  }
}

int[] decimal2base(int decimal, int digits, int base) {
  int base_rep[digits];
  int current = decimal;
  for(i in 1:digits) {
    base_rep[i] = current % base;
    current = current / base;
  }
  return base_rep;
}

int power_int(int number, int power);
  
int power_int(int number, int power) {
  if(power < 0)
    reject("Cannot raise an integer to a negative power and expect an integer result.");
  if (power == 0)
    return 1;
  else
    return number * power_int(number, power - 1);
}

// count the number of unique elements
int cardinality_int(int[] elems) {
  return rle_elem_count(sort_asc(elems));
}

// count the number of unique elements
int cardinality_vector(vector elems) {
  return rle_elem_count_vector(sort_asc(elems));
}

// create an integer sequence
int[] seq_int(int start, int end) {
  int N = end - start + 1;
  int seq[N];
  for(i in 1:N) seq[i] = i + start - 1;
  return(seq);
}

// return an int vector with each element in the input repeated as
// many times as given in the each argument
int[] rep_each(int[] set, int[] each) {
  int N = sum(each);
  int replicated[N];
  int p = 1;

  for(i in 1:size(set)) {
    replicated[p:p+each[i]-1] = rep_array(set[i], each[i]);
    p = p + each[i];
  }

  return(replicated);
}

/* calculate the absolute value of a - b in log-space with log(a)
   and log(b) given. Does so by squaring and taking the root, i.e.
   
   la = log(a)
   lb = log(b)
   
   sqrt( (a - b)^2 ) = sqrt( a^2 - 2 * a * b + b^2 )
   
   <=> 0.5 * log_diff_exp(log_sum_exp(2*la, 2*lb), log(2) + la + lb)
*/
real log_diff_exp_abs(real la, real lb) {
  return(0.5 * log_diff_exp(log_sum_exp(2*la, 2*lb), log(2) + la + lb));
}


/* find_interval, see findInterval from R
 * i.e. find the ranks of x in sorted; sorted is assumed to be weakly-sorted.
 */
int[] find_interval_slow(vector x, vector sorted) {
  int res[num_elements(x)];
  // very brutal and ineffcient way of doing this, but well, it's
  // C++ speed...
  for(i in 1:num_elements(x)) {
    res[i] = rank(append_row(rep_vector(x[i], 1), sorted), 1);
  }
  return(res);
}

/* faster version which uses bisectioning search
 */
int find_interval_elem(real x, vector sorted, int start_ind) {
  int res;
  int N;
  int max_iter;
  real left;
  real right;
  int left_ind;
  int right_ind;
  int iter;
    
  N = num_elements(sorted);
  
  if(N == 0) return(0);
  
  left_ind  = start_ind;
  right_ind = N;
  
  max_iter = 100 * N;
  left  = sorted[left_ind ] - x;
  right = sorted[right_ind] - x;
  
  if(0 <= left)  return(left_ind-1);
  if(0 == right) return(N-1);
  if(0 >  right) return(N);
  
  iter = 1;
  while((right_ind - left_ind) > 1  && iter != max_iter) {
    int mid_ind;
    real mid;
    // is there a controlled way without being yelled at with a
    // warning?
    mid_ind = (left_ind + right_ind) / 2;
    mid = sorted[mid_ind] - x;
    if (mid == 0) return(mid_ind-1);
    if (left  * mid < 0) { right = mid; right_ind = mid_ind; }
    if (right * mid < 0) { left  = mid; left_ind  = mid_ind; }
    iter = iter + 1;
  }
  if(iter == max_iter)
    print("Maximum number of iterations reached.");
  return(left_ind);
}

int[] find_interval(vector x, vector sorted) {
  int res[num_elements(x)];
  for(i in 1:num_elements(x)) {
    res[i] = find_interval_elem(x[i], sorted, 1);
  }
  return(res);
}

// takes as input x an ascending sorted vector x which allows to
// move the left starting index to be moved
int[] find_interval_asc(vector x, vector sorted) {
  int res[num_elements(x)];
  int last;
  last = 1;
  for(i in 1:num_elements(x)) {
    res[i] = find_interval_elem(x[i], sorted, last);
    if(res[i] > 0) last = res[i];
  }
  return(res);
}

int[] find_interval_blocked(int[] vals_M, vector vals, int[] sorted_M, vector sorted) {
  int res[num_elements(vals)];
  int M;
  int v;
  int s;
  M = num_elements(vals_M);
  v = 1;
  s = 1;
  for(m in 1:M) {
    int temp[vals_M[m]];
    temp = find_interval(segment(vals, v, vals_M[m]), segment(sorted, s, sorted_M[m]));
    for(n in 1:vals_M[m])
      res[v + n - 1] = temp[n];
    v = v + vals_M[m];
    s = s + sorted_M[m];
  }
  return(res);
}

// count number times elem appears in test set
int count_elem(int[] test, int elem) {
  int count;
  count = 0;
  for(i in 1:num_elements(test))
    if(test[i] == elem)
      count = count + 1;
  return(count);
}

// count number times elems appears in test set
int[] count_elems(int[] test, int[] elems) {
  int counts[num_elements(elems)];
  for(i in 1:num_elements(elems))
    counts[i] = count_elem(test, elems[i]);
  return(counts);
}

// find elements in test which are equal to elem
int[] which_elem(int[] test, int elem) {
  int res[count_elem(test, elem)];
  int ci;
  ci = 1;
  for(i in 1:num_elements(test))
    if(test[i] == elem) {
      res[ci] = i;
      ci = ci + 1;
    }
  return(res);
}

// divide fac by div and return the rounded down integer
int floor_div_int(real fac, real div) {
  int count;
  if(fac < 0)
    reject("floor_div_int only works for positive values.");
  count = 1;
  while(count * div <= fac) { count = count + 1; }
  count = count - 1;
  return count;
}

int[] count_obs_event_free(int[] obs_timeRank, int ndose) {
  int dose_next_obs[ndose];
  int o;
  int O;
  dose_next_obs = rep_array(0, ndose);
  o = 0;
  O = size(obs_timeRank);
  while (o < O && obs_timeRank[o+1] == 0) { o = o + 1; }
  for (i in 1:ndose) {
    int count;
    count = 0;
    while(o < O && obs_timeRank[o+1] == i) {
      o = o + 1;
      count = count + 1;
    }
    dose_next_obs[i] = count;
  }
  return(dose_next_obs);
}

int[] count_obs_event_free_blocked(int[] M, int[] obs_timeRank, int[] ndose) {
  int dose_next_obs[sum(ndose)];
  int l;
  int ld;
  dose_next_obs = rep_array(0, sum(ndose));
  l = 1;
  ld = 1;
  for (i in 1:size(M)) {
    int u;
    int ud;
    u = l + M[i] - 1;
    ud = ld + ndose[i] - 1;
    dose_next_obs[ld:ud] = count_obs_event_free(obs_timeRank[l:u], ndose[i]);
    l = u + 1;
    ld = ud + 1;
  }
  return(dose_next_obs);
}

/*
 * subsets the input data structure to the indices given as second
 * argument
 */
int[] subset_int(int[] cand, int[] ind_set) {
  int out[size(ind_set)];
  for(i in 1:size(ind_set))
    out[i] = cand[ind_set[i]];
  return out;
}
  
vector subset_vec(vector cand, int[] ind_set) {
  vector[size(ind_set)] out;
  for(i in 1:size(ind_set))
    out[i] = cand[ind_set[i]];
  return out;
}
  
matrix subset_matrix(matrix cand, int[] ind_set) {
  matrix[size(ind_set),cols(cand)] out;
  for(i in 1:size(ind_set))
    out[i] = cand[ind_set[i]];
  return out;
}
  
// check that assumption of 1...J labeling of rows hold and warn if
// not
void check_ids(int[] id) {
  int cid = 0;
  int warned = 0;
  cid = 0;
  for(n in 1:num_elements(id)) {
    if(id[n] != cid) {
      if(id[n] != cid + 1) {
        if(!warned)
          print("WARNING: id vector not correctly sorted, i.e. not in range 1..J. Consider using the cid vector internally.");
        warned = 1;
      } else {
        cid = cid + 1;
      }
    }
  }
  if(max(id) != cid)
    print("WARNING: Last patient's id not equal to max(id).");
}

// check that addl dose coding is correct. That is, we can only have a
// single active addl dose such that addl dosing must be off whenever
// the next dose occurs.
void check_addl_dosing(vector dose_time, vector dose_tau, int[] dose_addl) {
  int D = num_elements(dose_time);

  for(d in 2:D) {
    if(dose_time[d] < (dose_time[d-1] + dose_tau[d-1] * dose_addl[d-1]))
      reject("Forbidden overlapping dosing records found.");
  }
}

void check_addl_dosing_blocked(int[] dose_M, vector dose_time, vector dose_tau, int[] dose_addl) {
  int M = num_elements(dose_M);
  int b = 1;
  for(m in 1:M) {
    check_addl_dosing(segment(dose_time, b, dose_M[m]),
                      segment(dose_tau, b, dose_M[m]),
                      segment(dose_addl, b, dose_M[m]));
    b = b + dose_M[m];
  }
}

