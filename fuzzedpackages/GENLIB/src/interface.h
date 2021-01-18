////////////////////////////////////////////////////////////////////////////
// GenLib.h : main header file for GENLIB.DLL
// See GenLib.cxx for implementation.
// See GenLib.ssc for implementation of the corresponding S function.
//
// GENLIBDLL_DLL_API assists exporting/importing functions defined in GENLIB.DLL.
// It is for exporting if this header file is included in the source file (.c/.cxx) 
// that defines these functions, else importing.
//
// Symbols are exported if SP_CHAPTER_GENLIB is defined, else imported.
// SP_CHAPTER_GENLIB is automatically defined by the SPLUS Chapter Wizard. 
// Look at the setting..., C/C++ tap. 
////////////////////////////////////////////////////////////////////////////

#if !defined(S_GENLIBDLL_INCLUDED_)
#define S_GENLIBDLL_INCLUDED_

#undef GENLIBDLL_DLL_API
/*
#if defined(WIN32) && defined(SP_CHAPTER_GENLIBDLL)
	#define GENLIBDLL_DLL_API(returnType) __declspec(dllexport) returnType __stdcall
#elif defined(WIN32)
	#define GENLIBDLL_DLL_API(returnType) __declspec(dllimport) returnType __stdcall
#else
	#define GENLIBDLL_DLL_API(returnType) extern returnType
#endif
*/
#define GENLIBDLL_DLL_API(returnType) extern returnType
//#ifdef __cplusplus
extern "C" {
//#endif

GENLIBDLL_DLL_API(SEXP) SPLUSFlushCacheGenealogie();
GENLIBDLL_DLL_API(SEXP) SPLUSGetTimer(SEXP TimeInSec);
GENLIBDLL_DLL_API(SEXP) SPLUSValidateGenealogie(SEXP Genealogie, SEXP isValid);
GENLIBDLL_DLL_API(SEXP) SPLUSChangeMaxProcessingTime(SEXP newMaximum, SEXP oldMaximum);

GENLIBDLL_DLL_API(SEXP) SPLUSPhiMatrix(SEXP Genealogie, SEXP proposant, SEXP NProposant,SEXP Niveau, SEXP pdRetour, SEXP printit);

GENLIBDLL_DLL_API(SEXP) SPLUSPhiMatrixMT(SEXP Genealogie, SEXP proposant, SEXP NProposant,SEXP Niveau, SEXP pdRetour, SEXP printit);

GENLIBDLL_DLL_API(SEXP) SPLUSPhis( SEXP Genealogie, SEXP proposant, SEXP NProposant, SEXP NiveauMin, SEXP NiveauMax, SEXP pdRetour, 
							SEXP MatrixArray, SEXP printit);

GENLIBDLL_DLL_API(SEXP) SPLUSPhisMT(SEXP Genealogie, SEXP proposant, SEXP NProposant,SEXP NiveauMin,SEXP NiveauMax, SEXP pdRetour,
							 SEXP MatrixArray, SEXP printit);

GENLIBDLL_DLL_API(SEXP) SPLUSF(SEXP Genealogie, SEXP proposant, SEXP NProposant,SEXP Niveau, SEXP pdRetour, SEXP printit);

GENLIBDLL_DLL_API(SEXP) SPLUSFS(SEXP Genealogie, SEXP proposant, SEXP NProposant,SEXP NiveauMin,SEXP NiveauMax, SEXP pdRetour, SEXP printit);

GENLIBDLL_DLL_API(SEXP) SPLUSChild(SEXP Genealogie, SEXP plProposant,SEXP lNProposant, SEXP retour);

//GENLIBDLL_DLL_API(SEXP) SPLUSTestEbranche(SEXP sGenealogie, SEXP sProposant, SEXP sNProposant, SEXP sAncetre, SEXP sNAncetre, SEXP sRetour, SEXP sTaille);
GENLIBDLL_DLL_API(SEXP) SPLUSebranche(SEXP sGenealogie, SEXP sProposant, SEXP sNProposant, SEXP sAncetre, SEXP sNAncetre, SEXP sRetour, SEXP sTaille);  

GENLIBDLL_DLL_API(SEXP) SPLUSnumeroGen(SEXP Genealogie, SEXP plProposant,SEXP NProposant, SEXP retour);
GENLIBDLL_DLL_API(SEXP) SPLUSnumGenMin(SEXP Genealogie, SEXP plProposant,SEXP NProposant, SEXP retour);
GENLIBDLL_DLL_API(SEXP) SPLUSnumGenMoy(SEXP Genealogie, SEXP plProposant,SEXP NProposant, SEXP retour);

GENLIBDLL_DLL_API(SEXP) SPLUSConGen(SEXP slGenealogie, SEXP slProposant, SEXP sNProposant, SEXP slAncetre, SEXP sNAncetre, SEXP sdRetour, SEXP sprintit);

GENLIBDLL_DLL_API(SEXP) SPLUSConGenPLUS(SEXP plGenealogie, SEXP plProposant,SEXP lNProposant, SEXP plAncetre, SEXP lNAncetre, SEXP pdSexe,
								SEXP pdRetour, SEXP printit);

GENLIBDLL_DLL_API(SEXP) SPLUSCGCumul(SEXP Genealogie, SEXP plProposant,SEXP lNProposant, SEXP plAncetre, SEXP lNAncetre, SEXP AncRet,
							  SEXP pdRetour  ,SEXP pdRetourCumul, SEXP printit);

GENLIBDLL_DLL_API(SEXP) SPLUSCGCumuldirect(SEXP matriceCG, SEXP lNProposant, SEXP plAncetre, SEXP lNAncetre, SEXP AncRet, SEXP  pdSomAnc, SEXP pdSomCumul);

//GENLIBDLL_DLL_API(SEXP) SPLUSSimul(SEXP Genealogie, SEXP proposant, SEXP etatproposant,SEXP nproposant, SEXP ancetre, SEXP etatancetre, SEXP nancetre,
GENLIBDLL_DLL_API(SEXP) SPLUSSimul(SEXP Genealogie, SEXP proposant, SEXP etatproposant,SEXP nproposant, SEXP ancetre, SEXP etatancetre, SEXP nancetre,
							SEXP nSimul, SEXP pdRetConj, SEXP pdRetSimul, SEXP pdRetProp, SEXP sprobRecomb, SEXP sprobSurvieHomo,
							SEXP PrintProgress);

GENLIBDLL_DLL_API(SEXP) SPLUSSimulSingle(SEXP Genealogie, SEXP proposant, SEXP nproposant, SEXP ancetre, SEXP etatancetre, SEXP nancetre,
								 SEXP NSimul, SEXP pdRetour,SEXP PrintProgress);

GENLIBDLL_DLL_API(SEXP) SPLUSSimulSingleFreq(SEXP Genealogie, SEXP proposant, SEXP nproposant, SEXP ancetre, SEXP etatancetre, SEXP nancetre,
									SEXP NSimul, SEXP pdRetour,SEXP PrintProgress);

GENLIBDLL_DLL_API(SEXP) SPLUSSimulSingleProb(SEXP SGenealogie, SEXP SplProposant, SEXP SlNProposant, SEXP SplAncetre,SEXP SlNAncetres,
									SEXP SplAncEtat,SEXP mtProb, SEXP SlSimul,SEXP SPrintProgress);

GENLIBDLL_DLL_API(SEXP) SPLUSSimulSingleFct(SEXP SGenealogie, SEXP Sindividus, SEXP SplAncetre, SEXP SplAncEtatAll1, SEXP SplAncEtatAll2, 
								    SEXP SlNAncetre, SEXP SlSimul,SEXP SfctSousGrp,SEXP Sprintprogress);

GENLIBDLL_DLL_API(SEXP) SPLUSProb(SEXP Genealogie, SEXP proposant, SEXP etatproposant,SEXP nproposant, SEXP ancetre, SEXP etatancetre, SEXP nancetre,
						    SEXP pdRetConj, SEXP pdRetSimul,SEXP PrintProgress,SEXP onlyConj);

GENLIBDLL_DLL_API(SEXP) SPLUSCoeffApparentement(SEXP Genealogie, SEXP proposant, SEXP Nproposant, SEXP ancetre, SEXP Retour,
									   SEXP DuppDetection, SEXP printprogress);

GENLIBDLL_DLL_API(SEXP) SPLUSCALLCreerObjetGenealogie(SEXP SIndividu, SEXP SPere, SEXP SMere, SEXP SSexe);
//GENLIBDLL_DLL_API(SEXP) SPLUSCALLCreerObjetGenealogie(std::vector<double> SIndividu,std::vector<double> SPere,
//									  std::vector<double> SMere, std::vector<double> SSexe);
GENLIBDLL_DLL_API(SEXP) SPLUSOutgen(SEXP genealogie, SEXP plRetIndividu,SEXP plRetPere,SEXP plRetMere, SEXP plRetSexe, SEXP mustsort);

GENLIBDLL_DLL_API(SEXP) SPLUSOutIndice(SEXP genealogie, SEXP plRetIndividu,SEXP plRetPere,SEXP plRetMere, SEXP plRetSexe, SEXP mustsort);

GENLIBDLL_DLL_API(SEXP) SPLUSFondParGen(SEXP Genealogie,SEXP prop,SEXP nbProp, SEXP retour);

//#ifdef __cplusplus
}
//#endif

#endif //S_GENLIBDLL_INCLUDED_
