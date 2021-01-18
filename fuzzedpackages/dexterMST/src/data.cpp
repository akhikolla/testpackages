#include <RcppArmadillo.h>
#include <stack>

#define ALL 0
#define LAST 1

using namespace Rcpp;


// faster factor creation
// http://gallery.rcpp.org/articles/fast-factor-generation
// to do: levels unsorted, does it work?
template <int RTYPE>
IntegerVector fast_factor_template( const Vector<RTYPE>& x, bool as_int ) 
{
    Vector<RTYPE> levs = sort_unique(x);
    IntegerVector out = match(x, levs);
	if(!as_int)
	{
		out.attr("levels") = as<CharacterVector>(levs);
		out.attr("class") = "factor";
	}
    return out;
}

template <int RTYPE>
IntegerVector fast_factor_template( const Vector<RTYPE>& x, const Vector<RTYPE>& levs, bool as_int ) 
{
    IntegerVector out = match(x, levs);
	if(!as_int)
	{
		out.attr("levels") = as<CharacterVector>(levs);
		out.attr("class") = "factor";
	}
    return out;
}

// [[Rcpp::export]]
SEXP fast_factor( SEXP x, bool as_int) 
{
    switch( TYPEOF(x) ) {
    case INTSXP: return fast_factor_template<INTSXP>(x, as_int);
    case REALSXP: return fast_factor_template<REALSXP>(x, as_int);
    case STRSXP: return fast_factor_template<STRSXP>(x, as_int);
    }
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP fast_factor_lev( SEXP x, SEXP levs, bool as_int) 
{
    switch( TYPEOF(x) ) {
    case INTSXP: return fast_factor_template<INTSXP>(x, levs, as_int);
    case REALSXP: return fast_factor_template<REALSXP>(x, levs, as_int);
    case STRSXP: return fast_factor_template<STRSXP>(x, levs, as_int);
    }
    return R_NilValue;
}


//[[Rcpp::export]]
IntegerVector bid_c(const std::vector<std::string>& test_id, const std::vector<std::string>& booklet_id,
					const std::vector<std::string>& test_lev, const std::vector<std::string>& booklet_lev)
{
	typedef std::pair<std::string, std::string> key;
	
	struct key_hash
	{
		std::size_t operator()(const key& k) const
		{
			return std::hash<std::string>{}(k.first) ^ std::hash<std::string>{}(k.second);
		}	
	};
	
	std::unordered_map<key, int, key_hash> bid;
	
	const int ll = test_lev.size();
	const int n = test_id.size();
	
	IntegerVector out(n);
	
	// train
	for(int i=0; i<ll;i++)
		bid[std::make_pair(test_lev[i], booklet_lev[i])] = i+1;
	
# pragma omp parallel for
	for(int i=0; i<n; i++)
	{
		out[i] = bid[std::make_pair(test_id[i], booklet_id[i])];	
	}

	return out;
}

// [[Rcpp::export]]
bool is_person_booklet_sorted(const IntegerVector& booklet_id, const IntegerVector& person_id)
{
	const int n = booklet_id.length();
	std::atomic<bool> sorted(true);
	// prefer parallel time worst case over fast return if not sorted
#pragma omp parallel for
	for(int i=1; i<n; i++)
	{
		if((booklet_id[i] < booklet_id[i-1] && person_id[i] == person_id[i-1]) ||
		   (person_id[i] < person_id[i-1]))
		   sorted = false; 
	}
	return sorted;
}



// [[Rcpp::export]]
IntegerVector mutate_booklet_score(const IntegerVector& person_id, const IntegerVector& booklet_id, const IntegerVector& item_score)
{
	const int nr = person_id.length();
	int booklet = booklet_id[0];
	int person = person_id[0];
	int start = 0;
	int ss=0;

	IntegerVector booklet_score(nr);
	
	for(int r=0;r<nr;r++)
	{
		if(person != person_id[r] || booklet != booklet_id[r])
		{
			std::fill(booklet_score.begin() + start, booklet_score.begin() + r, ss);
			start = r;
			ss = 0;
			person = person_id[r];
			booklet = booklet_id[r];
		}
		ss += item_score[r];
	}
	//finally
	std::fill(booklet_score.begin() + start, booklet_score.end(), ss);
	
	return booklet_score;
}


/*

1) bestaand boekjes design +
	a) nieuwe itemsamenstelling in boekje (in modules maar geen overlap tussen modules dus kan per boekje)
	b) behaalde scoreputen die niet mee worden gerekend per module (dus patroon) m.u.v. de laatste module (die maken niet uit)
	c) leerlingen die geen enkele geincludeerde respons hebben worden verwijderd
	= nieuwe boekjes

2) per nieuw boekje willen we
	a) design (items) per module
	b) nieuwe routingregels
	c) identificatie van oorspronkelijk boekje waaraan ontleend
	d) aantal ll

3) per overgebleven leerling willen we
	a) welk nieuw boekje
	b) somscore geincludeerde responses
	c) afzonderlijke item scores

dus ~~~~
unsorted map van vector lengte max modules, bwaag
is er een minimal perfect hash?

key:
bool vector lengte items*max_score+1?
boekjes kunnen zelfde items hebben, samenvoegen is niet de bedoeling dus + 2 integers
maar 3 opties: score, geen score, niet aangeboden, aangeboden geincludeerd, niet geincludeerd , aaaaah

meest haalbaar:
bool vector items + int test + int booklet + reversible hash combo som geexcludeerd (value kan reverse zijn) 

*/


// [[Rcpp::export]]
List make_booklets_unsafe(
				const IntegerVector& person_id, IntegerVector& booklet_id, 
				const IntegerVector& module_nbr, const IntegerVector& item_id, 
				const IntegerVector& item_score, 
				IntegerVector& booklet_score, IntegerVector& include_rsp, 				
				const IntegerVector& bk_nmod ) 
{
	const int M = person_id.length();
	const int nit = as<CharacterVector>(item_id.attr("levels")).length();
	
	static int nmod = max(bk_nmod);
	
	std::vector<int> mod_omit(nmod+1,0);
	std::vector<bool> bk(nit+1, false);
	
	typedef std::tuple<std::vector<bool>, int, std::vector<int>> key;
	
	struct key_hash
	{
		const int sz = nmod;
		std::size_t operator()(const key& k) const
		{
			// hash bool vector xor int
			auto hash1 = std::hash<std::vector<bool>>{}(std::get<0>(k)); 
			auto hash2 = std::hash<int>{}(std::get<1>(k)); 
			hash1 ^= hash2 + 0x9e3779b9 + (hash1 << 6) + (hash1 >> 2);
			// first and last last one is guaranteed zero since min(modnbr)=1 and max score is set to zero 
			for(int m=1;m<sz; m++)
				hash1 ^= std::get<2>(k)[m]+ 0x9e3779b9 + (hash1 << 6) + (hash1 >> 2);			

			return hash1;
		}	
	};


	std::unordered_map<key, int, key_hash> booklets;
	
	
	int ss = 0;
	int bknum = 1;
	int person_start=0;
	
	if(include_rsp[0] == 1)
	{
		ss = item_score[0];
		bk[item_id[0]] = true;
	}
	else
	{
		mod_omit[module_nbr[0]] = item_score[0];	
	}
	
	for(int r=1; r < M; r++)
	{
		if(person_id[r] != person_id[r-1] || booklet_id[r] != booklet_id[r-1] )
		{
			
			if(std::any_of(bk.begin(), bk.end(), [](bool b){return b;}))
			{
				std::fill(booklet_score.begin()+person_start, booklet_score.begin() + r, ss);

				for(int m = bk_nmod[booklet_id[r-1]]; m <= nmod; m++)
					mod_omit[m] = 0;
				
				const auto& ret = booklets.insert( std::make_pair(std::make_tuple(bk, booklet_id[r-1], mod_omit), bknum));
				
				std::fill(booklet_id.begin()+person_start, booklet_id.begin() + r, ret.first->second); // moving to unique booklets
				
				if (ret.second) // did not already exist
					bknum++;
	
			}
			else
			{
				std::fill(include_rsp.begin()+person_start, include_rsp.begin() + r, 0);			
			}
			person_start = r;			
			std::fill(mod_omit.begin(), mod_omit.end(), 0);
			std::fill(bk.begin(), bk.end(), false);
			ss = 0;			
		}
		if(include_rsp[r] == 1)
		{
			bk[item_id[r]] = true;
			ss += item_score[r];
		}
		else
		{
			mod_omit[module_nbr[r]] += item_score[r];		
		}		
		
	}
	//finally
	if(std::any_of(bk.begin(), bk.end(), [](bool b){return b;}))
	{
		std::fill(booklet_score.begin()+person_start, booklet_score.end(), ss);
		
		for(int m = bk_nmod[booklet_id[M-1]]; m <= nmod; m++)
			mod_omit[m] = 0;

		const auto& ret = booklets.insert( std::make_pair(std::make_tuple(bk,  booklet_id[M-1], mod_omit), bknum));
				
		std::fill(booklet_id.begin()+person_start, booklet_id.end(), ret.first->second); 
	}
	else
	{
		std::fill(include_rsp.begin()+person_start, include_rsp.end() , 0);			
	}

	// booklet_id.bid is no longer a factor but a numeric id, remove levels
	booklet_id.attr("levels") = R_NilValue;
	booklet_id.attr("class") = "integer";	

	// transform hashtable to two dataframes

	const int rv_des = nit * booklets.size();
	const int rv_mod = nmod * booklets.size();
	
	std::vector<int> des_booklet_id, des_item_id, des_old_booklet_id; // test and booklet are necessary to couple items to modules
	std::vector<int>  mod_booklet_id, mod_module_nbr, mod_excl, mod_old_booklet_id;
	
	des_booklet_id.reserve(rv_des);
	des_item_id.reserve(rv_des);
	des_old_booklet_id.reserve(rv_des);
	
	mod_booklet_id.reserve(rv_mod);
	mod_module_nbr.reserve(rv_mod);
	mod_excl.reserve(rv_mod);
	mod_old_booklet_id.reserve(rv_mod);
	
	for(auto& iter: booklets )
	{
		int bk = iter.second;
		const auto& bki =  std::get<0>(iter.first);
		
		int old_bk = std::get<1>(iter.first);
		
		for(int i=1;i<=nit;i++) if(bki[i])
		{
			des_booklet_id.push_back(bk);
			des_item_id.push_back(i);
			des_old_booklet_id.push_back(old_bk);
		}
		
		int nm =  bk_nmod[old_bk];
		
		const auto& excl = std::get<2>(iter.first);
		
		for(int m = 1; m <= nm; m++)
		{
			mod_booklet_id.push_back(bk);
			mod_module_nbr.push_back(m);
			mod_excl.push_back(excl[m]);
			mod_old_booklet_id.push_back(old_bk);
		}		
	}	
	
	des_booklet_id.shrink_to_fit();
	des_item_id.shrink_to_fit();
	des_old_booklet_id.shrink_to_fit();
	
	mod_booklet_id.shrink_to_fit();
	mod_module_nbr.shrink_to_fit();
	mod_excl.shrink_to_fit();
	mod_old_booklet_id.shrink_to_fit();
	
	// to do: error when completely excluding modules? possibly best to infer from design in r
	
	return List::create(
		Named("booklets_items") = DataFrame::create(Named("bid") = des_booklet_id, Named("item_id") = des_item_id,
													Named("old_bid") = des_old_booklet_id),
		Named("booklets_modules") = DataFrame::create(Named("bid") = mod_booklet_id, Named("module_nbr") = mod_module_nbr, Named("mod_excluded") = mod_excl,
													  Named("old_bid") = mod_old_booklet_id));	

}


// [[Rcpp::export]]
List suf_stats_nrm_c(const IntegerVector& booklet_id, const IntegerVector& booklet_score, const IntegerVector& item_id, const IntegerVector& item_score, const int nit, const int max_score)
{
	typedef std::tuple<int, int, int> key;
	typedef std::pair<int, int> val;
	
	const int n = booklet_id.length();	
	
	// prepare map for plt
	struct int3_hash
	{
		std::size_t operator()(const key& k) const
		{
			// booklet, booklet_score, item
			return std::hash<int>()(std::get<0>(k) * 128 + std::get<1>(k)  + std::get<2>(k) * 8192);
		}	
	};

	std::unordered_map<key, val, int3_hash> plt;

	// reserve space for ssIS
	const int nscore = max_score + 1;	
	const int ssIS_len = nit * nscore;
	std::vector<int> ssI(ssIS_len, 0), ssi_item, ssi_item_score;
	ssi_item.reserve(ssIS_len);
	ssi_item_score.reserve(ssIS_len);
	
	// main loop
	for(int i=0; i<n; i++)
	{
		// plt
		auto& s = plt[std::forward_as_tuple(booklet_id[i], booklet_score[i], item_id[i])]; 
		std::get<0>(s) += item_score[i];
		std::get<1>(s)++;
		// ssI
		ssI[(item_id[i]-1) * nscore + item_score[i]]++;
	}
	
	// convert map to data.frame for plot
	const int sz = plt.size();
	IntegerVector plt_booklet_id(sz), plt_item_id(sz), plt_booklet_score(sz), plt_n(sz);
	NumericVector plt_mean(sz);
	
	int i=0;
	for(std::unordered_map<key, val, int3_hash>::iterator iter = plt.begin(); iter != plt.end(); ++iter)
	{
		auto& k = iter->first;
		auto& s = iter->second;
		
		plt_booklet_id[i] = std::get<0>(k);
		plt_booklet_score[i] = std::get<1>(k);
		plt_item_id[i] = std::get<2>(k);
		plt_n[i] = std::get<1>(s);
		plt_mean[i++] = ((double)(std::get<0>(s)))/std::get<1>(s);		
	}
	
	if(Rf_isFactor(booklet_id))
	{
		plt_booklet_id.attr("levels") = booklet_id.attr("levels");	
		plt_booklet_id.attr("class") = "factor";
	}
	plt_item_id.attr("levels") = item_id.attr("levels");	
	plt_item_id.attr("class") = "factor";
	
	// shrink ssIS to remove non-existing scores
	int indx = 0;
	for(int itm=0; itm<nit; itm++)
	{
		for(int s=0; s<=max_score; s++)
		{
			if(ssI[itm * nscore + s] > 0)
			{
				ssI[indx++] = ssI[itm * nscore + s];
				ssi_item.push_back(itm+1); // 1-indexed
				ssi_item_score.push_back(s);
			}
		}
	}
	
	ssi_item.shrink_to_fit();
	ssi_item_score.shrink_to_fit();
	ssI.resize(indx);
	
	return List::create(
		Named("plt") = DataFrame::create(Named("bid") = plt_booklet_id, Named("booklet_score") = plt_booklet_score, 
										 Named("item_id") = plt_item_id, Named("meanScore") = plt_mean, Named("N") = plt_n),
		Named("ssIS") = DataFrame::create(Named("item_id") = ssi_item, Named("item_score") = ssi_item_score, Named("sufI") = ssI));	

}



// sumscore for single booklet
// [[Rcpp::export]]
std::vector<int> im_booklet_score(const IntegerVector& person_id, const IntegerVector& item_score)
{
	const int n = person_id.length();
	std::vector<int> booklet_score(n);
	
	int s = item_score[0];	
	int start = 0;
	
	for(int i=1; i<n; i++)
	{
		if(person_id[i]==person_id[i-1])
		{
			s+=item_score[i];
		}
		else
		{
			std::fill(booklet_score.begin()+start, booklet_score.begin()+i, s);
			s=item_score[i];
			start=i;
		}
	}
	std::fill(booklet_score.begin()+start, booklet_score.end(), s);
	return booklet_score;
}


// [[Rcpp::export]]
List suf_stats_im_c(const IntegerVector& booklet_score, const IntegerVector& item_id, const IntegerVector& item_score, const int nit, const int max_score)
{
	const int szi = 1+max_score;

	const int n = booklet_score.length();
	
	std::vector<int> ss_is_c(szi+nit*szi,0), ss_is_i(szi+nit*szi,0);
	
	const int bmax = nit*max_score;
	
	std::vector<int> scoretab(bmax+1,0), plt_sum((bmax+2)*nit,0);
	
	for(int r=0;r<n;r++)
	{
		int ii = item_id[r] * szi + item_score[r];
		if(item_id[r] == 1)
			scoretab[booklet_score[r]]++; 
		plt_sum[booklet_score[r] * nit + item_id[r]] += item_score[r]; 
		ss_is_c[ii] += item_score[r] * booklet_score[r];
		ss_is_i[ii]++;	
	}
	
	std::vector<int> ssi_item_id, ssi_item_score, plt_booklet_score, plt_item_id,plt_n;
	std::vector<double> plt_ms;
	
	plt_ms.reserve(nit*(1+bmax));
	plt_n.reserve(nit*(1+bmax));
	plt_item_id.reserve(nit*(1+bmax));
	plt_booklet_score.reserve(nit*(1+bmax));
	
	ssi_item_id.reserve(nit*max_score);
	ssi_item_score.reserve(nit*max_score);

	int sindx=0,pindx=0;
	for(int i=1; i<=nit; i++)
	{
		for(int s=1; s<=max_score; s++)
		{
			int ii = i * szi + s;
			if(ss_is_i[ii]>0)
			{
				ssi_item_score.push_back(s);
				ssi_item_id.push_back(i);
				ss_is_i[sindx] = ss_is_i[ii];
				ss_is_c[sindx++] = ss_is_c[ii];
			}
		}
		for(int bs=0; bs<bmax; bs++)
		{
			int bi = bs*nit+i;
			if(scoretab[bs]>0)
			{
				plt_n.push_back(scoretab[bs]);
				plt_booklet_score.push_back(bs);
				plt_item_id.push_back(i);
				plt_ms.push_back(plt_sum[bi]/((double)scoretab[bs]));
				pindx++;
			}
		}
	}		
	
	plt_booklet_score.shrink_to_fit();
	plt_ms.shrink_to_fit();
	plt_item_id.shrink_to_fit();
	plt_n.shrink_to_fit();
	// to do: check if resize is guaranteed (so we may read size)
	ssi_item_id.shrink_to_fit();
	ssi_item_score.shrink_to_fit();
	ss_is_i.resize(ssi_item_id.size());
	ss_is_c.resize(ssi_item_id.size());
	
	
	
	
	return List::create(
		Named("plt") = DataFrame::create(Named("booklet_score") = plt_booklet_score, 
										 Named("item_id") = plt_item_id, Named("mean_item_score") = plt_ms, Named("n") = plt_n),
		Named("ssIS") = DataFrame::create(Named("item_id") = ssi_item_id, Named("item_score") = ssi_item_score, 
											Named("sufI") = ss_is_i, Named("sufC") = ss_is_c),
		Named("scoretab") = scoretab);			

}


// [[Rcpp::export]]
bool is_connected_C( const IntegerMatrix& A)
{
	const int n = A.ncol();		
	std::vector<bool> visited(n, false);
	std::stack<int> st;

	st.push(0);
	
	// non recursive DFS
	while(!st.empty())
	{
		int s = st.top(); 
        st.pop();
		visited[s] = true;
		
		for(int i=0; i<n; i++) 
			if(A(i,s)>0 && !visited[i])
				st.push(i);
	}
	for(int i=0; i<n; i++)
		if(!visited[i])
			return false;
	
	return(true);
}

