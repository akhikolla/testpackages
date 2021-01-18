#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>

// #include <R.h>
// #include <Rinternals.h>
// #include <R_ext/Rdynload.h>
#include "kmcuda.h"

// need this to avoid assert conflicts with R
#define NDEBUG

#include <RcppEigen.h>

namespace {
std::unordered_map<std::string, SEXP> parse_args(
		const std::unordered_set<std::string> &allowed,
		std::initializer_list<std::string> required, SEXP args) {
	std::unordered_map<std::string, SEXP> result;
	args = CDR(args);
	bool pure = true;
	for (unsigned i = 0; args != R_NilValue; i++, args = CDR(args)) {
		SEXP arg = CAR(args);
		if (Rf_isNull(TAG(args))) {
			if (pure && i == result.size()) {
				result.emplace(required.begin()[i], arg);
			} else {
				Rcpp::stop("positional argument follows keyword argument");
			}
		} else {
			pure = false;
			const char *name = CHAR(PRINTNAME(TAG(args)));
			if (allowed.find(name) == allowed.end()) {
				Rcpp::stop("got an unexpected keyword argument \"%s\"", name);
			}
			result.emplace(name, arg);
		}
	}
	return result;
}
	
	template<typename R>
	R parse_dict(const std::unordered_map<std::string, R>& dict, const char *arg,
              SEXP name) {
		if (!Rf_isString(name)) {
			Rcpp::stop("\"%s\" name must be a string", arg);
		}
		const char *init_name = CHAR(Rf_asChar(name));
		auto init_iter = dict.find(init_name);
		if (init_iter == dict.end()) {
			Rcpp::stop("\"%s\" = \"%s\" is not supported", arg, init_name);
		}
		return init_iter->second;
	}
	
	int parse_int(SEXP value) {
		if (Rf_isInteger(value)) {
			return INTEGER(value)[0];
		}
		return REAL(value)[0];
	}
	
	int parse_int(const std::unordered_map<std::string, SEXP> &kwargs,
               const std::string &name, int def) {
		auto iter = kwargs.find(name);
		if (iter == kwargs.end()) {
			return def;
		}
		if (!Rf_isNumeric(iter->second)) {
			Rcpp::stop("\"%s\" must be an integer", name.c_str());
		}
		return parse_int(iter->second);
	}
	
// 	void parse_samples(
// 			const std::unordered_map<std::string, SEXP> &kwargs,
// 			std::unique_ptr<float[]> *samples, uint32_t *samples_size,
// 			uint16_t *features_size) {
// 		std::unique_ptr<SEXP[]> samples_chunks;
// 		int chunks_size = 0;
// 		{
// 			auto samples_iter = kwargs.find("samples");
// 			if (samples_iter == kwargs.end()) {
// 				Rcpp::stop("\"samples\" must be defined");
// 			}
// 			SEXP samples_obj = samples_iter->second;
// 			if (Rf_isReal(samples_obj)) {
// 				chunks_size = 1;
// 				samples_chunks.reset(new SEXP[1]);
// 				samples_chunks[0] = samples_obj;
// 			} else if (TYPEOF(samples_obj) == VECSXP) {
// 				chunks_size = Rf_length(samples_obj);
// 				samples_chunks.reset(new SEXP[chunks_size]);
// 				for (int i = 0; i < chunks_size; i++) {
// 					samples_chunks[i] = VECTOR_ELT(samples_obj, i);
// 				}
// 			} else {
// 				Rcpp::stop("\"samples\" must be a 2D real matrix or a vector of 2D real matrices");
// 			}
// 		}
// 		*samples_size = 0;
// 		for (int i = 0; i < chunks_size; i++) {
// 			if (!Rf_isReal(samples_chunks[i])) {
// 				Rcpp::stop("\"samples\" must be a 2D real matrix or a vector of 2D real matrices");
// 			}
// 			SEXP dims = Rf_getAttrib(samples_chunks[i], R_DimSymbol);
// 			if (Rf_length(dims) != 2) {
// 				Rcpp::stop("\"samples\" must be a 2D real matrix or a vector of 2D real matrices");
// 			}
// 			int samples_size_i = INTEGER(dims)[0];
// 			if (static_cast<int64_t>(*samples_size) + samples_size_i > INT32_MAX) {
// 				Rcpp::stop("too many samples (>INT32_MAX)");
// 			}
// 			*samples_size += samples_size_i;
// 			int features_size_i = INTEGER(dims)[1];
// 			if (features_size_i > UINT16_MAX) {
// 				Rcpp::stop("too many features (>UINT16_MAX)");
// 			}
// 			if (i == 0) {
// 				*features_size = features_size_i;
// 			} else if (*features_size != features_size_i) {
// 				Rcpp::stop("\"samples\" vector contains matrices with different number of columns");
// 			}
// 		}
// 		samples->reset(new float[
//                    static_cast<uint64_t>(*samples_size) * *features_size]);
// 		float *samples_float = samples->get();
// 		{
// 			int offset = 0;
// 			for (int i = 0; i < chunks_size; i++) {
// 				double *samples_double = REAL(samples_chunks[i]);
// 				SEXP dims = Rf_getAttrib(samples_chunks[i], R_DimSymbol);
// 				uint32_t fsize = *features_size;
// 				uint32_t ssize = INTEGER(dims)[0];
// #pragma omp parallel for
// 				for (uint64_t f = 0; f < fsize; f++) {
// 					for (uint64_t s = 0; s < ssize; s++) {
// 						samples_float[offset + s * fsize + f] = samples_double[f * ssize + s];
// 					}
// 				}
// 				offset += ssize * fsize;
// 			}
// 		}
// 	}
	
	KMCUDADistanceMetric parse_metric(
			const std::unordered_map<std::string, SEXP> &kwargs) {
		KMCUDADistanceMetric metric = kmcudaDistanceMetricL2;
		auto metric_iter = kwargs.find("metric");
		if (metric_iter != kwargs.end()) {
			metric = parse_dict(kmcuda::metrics, "metric", metric_iter->second);
		}
		return metric;
	}
	
	KMCUDADistanceMetric parse_metric(SEXP metric_sexp){
		KMCUDADistanceMetric metric;
		if(Rf_isNull(metric_sexp)){
			metric = kmcudaDistanceMetricL2;	
		}else{
			metric = parse_dict(kmcuda::metrics, "metric", metric_sexp);	
		}
		return metric;
	}
	
	// int parse_device(
	// 		const std::unordered_map<std::string, SEXP> &kwargs) {
	// 	int device = parse_int(kwargs, "device", 0);
	// 	if (device < 0) {
	// 		Rcpp::stop("\"device\" may not be negative");
	// 	}
	// 	return device;
	// }
	
	// static const std::unordered_set<std::string> kmeans_kwargs {
	// 	"samples", "clusters", "tolerance", "init", "yinyang_t", "metric",
	// 	"average_distance", "seed", "device", "verbosity"
	// };
	
	// static const std::unordered_set<std::string> knn_kwargs {
	// 	"k", "samples", "centroids", "assignments", "metric", "device",
	// 	"verbosity"
	// };
}

// extern "C" {
	
// static SEXP r_kmeans_cuda(SEXP args);
// static SEXP r_knn_cuda(SEXP args);

// static R_ExternalMethodDef externalMethods[] = {
// 	{"kmeans_cuda", (DL_FUNC) &r_kmeans_cuda, -1},
// 	{"knn_cuda", (DL_FUNC) &r_knn_cuda, -1},
// 	{NULL, NULL, 0}
// };

// void R_init_libKMCUDA(DllInfo *info) {
// 	R_registerRoutines(info, NULL, NULL, NULL, externalMethods);
// }

// [[Rcpp::export]]
SEXP r_kmeans_cuda(
		Rcpp::NumericMatrix samples,
		int clusters,
		SEXP tolerance_sexp,
		SEXP init_sexp,
		SEXP yinyang_t_sexp,
		SEXP metric_sexp,
		bool average_distance_bool,
		SEXP seed_sexp,
		int device,
		int verbosity) 
{
	// auto kwargs = parse_args(kmeans_kwargs, {"samples", "clusters"}, args);
	// std::unique_ptr<float[]> samples;
	uint32_t samples_size;
	uint16_t features_size;
	// parse_samples(kwargs, &samples, &samples_size, &features_size);
	samples_size = samples.nrow();
	features_size = samples.ncol();
	
	// SEXP clusters_obj = kwargs["clusters"];
	// if (!Rf_isNumeric(clusters_obj)) {
	// 	Rcpp::stop("\"clusters\" must be a positive integer");
	// }
	// int clusters_size = parse_int(clusters_obj);
	// if (clusters_size <= 0) {
	// 	Rcpp::stop("\"clusters\" must be a positive integer");
	// }
	
	int clusters_size = clusters;
	
	if (static_cast<uint64_t>(clusters_size) * features_size > INT32_MAX
      || static_cast<uint32_t>(clusters_size) >= samples_size) {
		Rcpp::stop("\"clusters\" is too big");
	}
	
	
	auto centroids = std::unique_ptr<float[]>(
		new float[clusters_size * features_size]);
	
	KMCUDAInitMethod init = kmcudaInitMethodPlusPlus;
	int afkmc2_m = 0;
	
	if(Rf_isString(init_sexp)){
		init = parse_dict(kmcuda::init_methods, "init", init_sexp);
	}else if (Rf_isReal(init_sexp)){
		init = kmcudaInitMethodImport;
		SEXP dims = Rf_getAttrib(init_sexp, R_DimSymbol);
		if (Rf_length(dims) != 2
        || INTEGER(dims)[0] != clusters_size
        || INTEGER(dims)[1] != features_size) {
        Rcpp::stop("invalid centroids dimensions in \"init\"");
		}
		
		double *centroids_double = REAL(init_sexp);
		#pragma omp parallel for
		for (uint64_t f = 0; f < features_size; f++) {
			for (int64_t c = 0; c < clusters_size; c++) {
				centroids[c * features_size + f] = centroids_double[f * clusters_size + c];
			}
		}
	}else{
		Rcpp::stop("\"init\" must be either a string or a list or a 2D matrix");
	}
	
// 	auto init_iter = kwargs.find("init");
// 	if (init_iter != kwargs.end()) {
// 		if (Rf_isString(init_iter->second)) {
// 			init = parse_dict(kmcuda::init_methods, "init", init_iter->second);
// 		} else if (TYPEOF(init_iter->second) == VECSXP) {
// 			if (Rf_length(init_iter->second) == 0) {
// 				Rcpp::stop("\"init\" may not be an empty list");
// 			}
// 			init = parse_dict(kmcuda::init_methods, "init", CAR(init_iter->second));
// 			if (init == kmcudaInitMethodAFKMC2 && Rf_length(init_iter->second) > 1) {
// 				SEXP afkmc2_m_obj = CAAR(init_iter->second);
// 				if (!Rf_isNumeric(afkmc2_m_obj)) {
// 					Rcpp::stop("\"init\" = %s: parameter must be a positive integer",
//           CHAR(Rf_asChar(CAR(init_iter->second))));
// 				}
// 				afkmc2_m = parse_int(afkmc2_m_obj);
// 				if (afkmc2_m <= 0) {
// 					Rcpp::stop("\"init\" = %s: parameter must be a positive integer",
//           CHAR(Rf_asChar(CAR(init_iter->second))));
// 				}
// 			} else if (Rf_length(init_iter->second) != 1) {
// 				Rcpp::stop("\"init\" has wrong number of parameters");
// 			}
// 		} else if (Rf_isReal(init_iter->second)) {
// 			init = kmcudaInitMethodImport;
// 			SEXP dims = Rf_getAttrib(init_iter->second, R_DimSymbol);
// 			if (Rf_length(dims) != 2
//         || INTEGER(dims)[0] != clusters_size
//         || INTEGER(dims)[1] != features_size) {
//         Rcpp::stop("invalid centroids dimensions in \"init\"");
// 			}
// 			double *centroids_double = REAL(init_iter->second);
// #pragma omp parallel for
// 			for (uint64_t f = 0; f < features_size; f++) {
// 				for (int64_t c = 0; c < clusters_size; c++) {
// 					centroids[c * features_size + f] = centroids_double[f * clusters_size + c];
// 				}
// 			}
// 		} else {
// 			Rcpp::stop("\"init\" must be either a string or a list or a 2D matrix");
// 		}
// 	}
	
	float tolerance = Rcpp::as<float>(tolerance_sexp);
	// auto tolerance_iter = kwargs.find("tolerance");
	// if (tolerance_iter != kwargs.end()) {
	// 	if (!Rf_isReal(tolerance_iter->second)) {
	// 		Rcpp::stop("\"tolerance\" must be a real value");
	// 	}
	// 	tolerance = REAL(tolerance_iter->second)[0];
	// 	if (tolerance < 0 || tolerance > 1) {
	// 		Rcpp::stop("\"tolerance\" must be in [0, 1]");
	// 	}
	// }
	
	float yinyang_t = Rcpp::as<float>(yinyang_t_sexp);
	
	// auto yinyang_t_iter = kwargs.find("yinyang_t");
	// if (yinyang_t_iter != kwargs.end()) {
	// 	if (!Rf_isReal(yinyang_t_iter->second)) {
	// 		Rcpp::stop("\"yinyang_t\" must be a real value");
	// 	}
	// 	yinyang_t = REAL(yinyang_t_iter->second)[0];
	// 	if (yinyang_t < 0 || yinyang_t > 0.5) {
	// 		Rcpp::stop("\"tolerance\" must be in [0, 0.5]");
	// 	}
	// }
	
	// KMCUDADistanceMetric metric = parse_metric(kwargs);
	KMCUDADistanceMetric metric = parse_metric(metric_sexp);
	
	// uint32_t seed = parse_int(kwargs, "seed", time(NULL));
	uint32_t seed = Rcpp::as<uint32_t>(seed_sexp);
	
	// int device = parse_device(kwargs);
	// int verbosity = parse_int(kwargs, "verbosity", 0);
	
	float average_distance, *average_distance_ptr = nullptr;
	if(average_distance_bool){
		average_distance_ptr = &average_distance;
	}
	// auto average_distance_iter = kwargs.find("average_distance");
	// if (average_distance_iter != kwargs.end()) {
	// 	if (LOGICAL(average_distance_iter->second)[0]) {
	// 		average_distance_ptr = &average_distance;
	// 	}
	// }
	
	auto assignments = std::unique_ptr<uint32_t[]>(new uint32_t[samples_size]);
	
	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> samples_eigen =
		Rcpp::as<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> >(samples);
	
	auto result = kmeans_cuda(
		init, &afkmc2_m, tolerance, yinyang_t, metric, samples_size, features_size,
		clusters_size, seed, device, -1, 0, verbosity, &samples_eigen(0),
		centroids.get(), assignments.get(), average_distance_ptr);
	
	if (result != kmcudaSuccess) {
		Rcpp::stop("kmeans_cuda error %d %s%s", result,
       kmcuda::statuses.find(result)->second, (verbosity > 0)? "" : "; "
         "\"verbosity\" = 2 would reveal the details");
	}
	// SEXP centroids2 = PROTECT(Rf_allocMatrix(REALSXP, clusters_size, features_size));
	Rcpp::NumericMatrix centroids2(clusters_size, features_size);
	
	// double *centroids_double = REAL(centroids2);
	float *centroids_float = centroids.get();
#pragma omp parallel for
	for (uint64_t f = 0; f < features_size; f++) {
		for (int64_t c = 0; c < clusters_size; c++) {
			// centroids_double[f * clusters_size + c] = centroids_float[c * features_size + f];
			// centroids2[c, f] = centroids_float[c * features_size + f];
			centroids2[f * clusters_size + c] = centroids_float[c * features_size + f];
		}
	}
	// SEXP assignments2 = PROTECT(Rf_allocVector(INTSXP, samples_size));
	Rcpp::IntegerVector assignments2(samples_size);
	// assignments2.
	
	// uint32_t *assignments_ptr = assignments.get();
	int *assignments2_ptr = INTEGER(assignments2);
#ifndef __APPLE__
#pragma omp parallel for simd
	for (uint32_t i = 0; i < samples_size; i++) {
		assignments2_ptr[i] = assignments[i] + 1;
	}
#else
#pragma omp simd
	for (uint32_t i = 0; i < samples_size; i++) {
		assignments2_ptr[i] = assignments[i] + 1;
	}
#endif
	// SEXP tuple = PROTECT(Rf_allocVector(VECSXP, 2 + (average_distance_ptr != nullptr)));
	// SET_VECTOR_ELT(tuple, 0, centroids2);
	// SET_VECTOR_ELT(tuple, 1, assignments2);
	// SEXP names = PROTECT(Rf_allocVector(
	// 	STRSXP, 2 + (average_distance_ptr != nullptr)));
	// SET_STRING_ELT(names, 0, Rf_mkChar("centroids"));
	// SET_STRING_ELT(names, 1, Rf_mkChar("assignments"));
	// if (average_distance_ptr != nullptr) {
	// 	SEXP average_distance2 = PROTECT(Rf_allocVector(REALSXP, 1));
	// 	REAL(average_distance2)[0] = average_distance;
	// 	SET_VECTOR_ELT(tuple, 2, average_distance2);
	// 	SET_STRING_ELT(names, 2, Rf_mkChar("average_distance"));
	// }
	// Rf_setAttrib(tuple, R_NamesSymbol, names);
	// UNPROTECT(4 + (average_distance_ptr != nullptr));
	// return tuple;
	
	if(average_distance_ptr != nullptr){
		return Rcpp::List::create(Rcpp::Named("centroids") = Rcpp::wrap(centroids2),
                            Rcpp::Named("assignments") = Rcpp::wrap(assignments2),
                            Rcpp::Named("average_distance") = Rcpp::wrap(average_distance)
		);
	}else{
		return Rcpp::List::create(Rcpp::Named("centroids") = Rcpp::wrap(centroids2),
                            Rcpp::Named("assignments") = Rcpp::wrap(assignments2)
		);
	}
	
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix r_knn_cuda(
		int k,
		Rcpp::NumericMatrix samples,
		Rcpp::NumericMatrix centroids_sexp,
		Rcpp::IntegerVector assignments_sexp,
		SEXP metric_sexp,
		int device,
		int verbosity
		) {
	// auto kwargs = parse_args(
	// 	knn_kwargs, {"k", "samples", "centroids", "assignments"}, args);
	// int k = parse_int(kwargs, "k", 0);
	// if (k <= 0) {
	// 	Rcpp::stop("\"k\" must be positive");
	// }
	// std::unique_ptr<float[]> samples;
	uint32_t samples_size;
	uint16_t features_size;
	// parse_samples(kwargs, &samples, &samples_size, &features_size);
	samples_size = samples.nrow();
	features_size = samples.ncol();
	
	if (static_cast<uint64_t>(samples_size) * k > INT32_MAX) {
		Rcpp::stop("too big \"k\": dim(samples)[0] * k > INT32_MAX");
	}
	// auto centroids_iter = kwargs.find("centroids");
	// if (centroids_iter == kwargs.end()) {
	// 	Rcpp::stop("\"centroids\" must be specified");
	// }
	// if (!Rf_isReal(centroids_iter->second)) {
	// 	Rcpp::stop("\"centroids\" must be a 2D real matrix");
	// }
	// SEXP dims = Rf_getAttrib(centroids_iter->second, R_DimSymbol);
	// if (Rf_length(dims) != 2 || INTEGER(dims)[1] != features_size) {
	// 	Rcpp::stop("invalid \"centroids\"'s dimensions");
	// }
	// int clusters_size = INTEGER(dims)[0];
	int clusters_size = centroids_sexp.nrow();
	std::unique_ptr<float[]> centroids(new float[clusters_size * features_size]);
	// double *centroids_double = REAL(centroids_iter->second);
	float *centroids_float = centroids.get();
#pragma omp parallel for
	for (uint64_t f = 0; f < features_size; f++) {
		for (int64_t c = 0; c < clusters_size; c++) {
			centroids_float[c * features_size + f] = centroids_sexp[f * clusters_size + c];
		}
	}
	// auto assignments_iter = kwargs.find("assignments");
	// if (assignments_iter == kwargs.end()) {
	// 	Rcpp::stop("\"assignments\" must be specified");
	// }
	// if (!Rf_isInteger(assignments_iter->second)) {
	// 	Rcpp::stop("\"assignments\" must be an integer vector");
	// }
	// if (static_cast<uint32_t>(Rf_length(assignments_iter->second)) != samples_size) {
	// 	Rcpp::stop("invalid \"assignments\"'s length");
	// }
	
	std::unique_ptr<uint32_t[]> assignments(new uint32_t[samples_size]);
	// int *assignments_obj_ptr = INTEGER(assignments_iter->second);
	uint32_t *assignments_ptr = assignments.get();
#ifndef __APPLE__
#pragma omp parallel for simd
	for (uint32_t i = 0; i < samples_size; i++) {
		assignments_ptr[i] = assignments_sexp[i] - 1;
	}
#else
#pragma omp simd
	for (uint32_t i = 0; i < samples_size; i++) {
		assignments_ptr[i] = assignments_sexp[i] - 1;
	}
#endif
	
	// KMCUDADistanceMetric metric = parse_metric(kwargs);
	KMCUDADistanceMetric metric = parse_metric(metric_sexp);
	
	// int device = parse_device(kwargs);
	// int verbosity = parse_int(kwargs, "verbosity", 0);
	
	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> samples_eigen =
		Rcpp::as<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> >(samples);
	
	std::unique_ptr<uint32_t[]> neighbors(new uint32_t[samples_size * k]);
	auto result = knn_cuda(
		k, metric, samples_size, features_size, clusters_size, device, -1, 0,
		verbosity, &samples_eigen(0), centroids.get(), assignments_ptr, neighbors.get());
	if (result != kmcudaSuccess) {
		Rcpp::stop("knn_cuda error %d %s%s", result,
       kmcuda::statuses.find(result)->second, (verbosity > 0)? "" : "; "
         "\"verbosity\" = 2 would reveal the details");
	}
	
	// SEXP neighbors_obj = PROTECT(Rf_allocMatrix(INTSXP, samples_size, k));
	Rcpp::IntegerMatrix neighbors_obj(samples_size, k);
	
	const uint32_t *neighbors_ptr = neighbors.get();
	// int *neighbors_obj_ptr = INTEGER(neighbors_obj);
#pragma omp parallel for
	for (int i = 0; i < k; i++) {
		for (uint32_t s = 0; s < samples_size; s++) {
			neighbors_obj[i * samples_size + s] = neighbors_ptr[s * k + i] + 1;
		}
	}
	// UNPROTECT(1);
	return neighbors_obj;
}
	
// }  // extern "C"
