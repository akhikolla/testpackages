///////////////////////////////////////////////////////////////////////////
// genlibDLL.h : main header file for GENLIBDLL.DLL
// See genlibDLL.cxx for implementation.
// See genlibDLL.ssc for implementation of the corresponding S function.
//
// GENLIBDLL_DLL_API assists exporting/importing functions defined in GENLIBDLL.DLL.
// It is for exporting if this header file is included in the source file (.c/.cxx) 
// that defines these functions, else importing.
//
// Symbols are exported if SP_CHAPTER_GENLIBDLL is defined, else imported.
// SP_CHAPTER_GENLIBDLL is automatically defined by the SPLUS Chapter Wizard. 
// Look at the setting..., C/C++ tap. 
////////////////////////////////////////////////////////////////////////////

#ifndef S_GENLIBDLL_INCLUDED_
//#if !defined(S_GENLIBDLL_INCLUDED_)
#define S_GENLIBDLL_INCLUDED_

#undef GENLIBDLL_DLL_API
#if defined(WIN32) && defined(SP_CHAPTER_GENLIBDLL)
	#define GENLIBDLL_DLL_API(returnType) __declspec(dllexport) returnType __stdcall
#elif defined(WIN32)
	#define GENLIBDLL_DLL_API(returnType) __declspec(dllimport) returnType __stdcall
#else
	#define GENLIBDLL_DLL_API(returnType) extern returnType
#endif

#ifdef __cplusplus
extern "C" {
#endif

__declspec (.dll) void genlibDLLC(double*, double*, int* ); 
__declspec (S.dll) (SEXP*) genlibDLLCall(SEXP* );  //__declspec (S.dll) (s_object*) genlibDLLCall(s_object* ); 
__declspec (S.dll) (SEXP*) genlibDLLCall2(SEXP* );  //__declspec (S.dll) (s_object*) genlibDLLCall2(s_object* ); 
GENLIBDLL_DLL_API(void) genlibDLLC(double*, double*, int* ); 
GENLIBDLL_DLL_API(s_object*) genlibDLLCall(SEXP* );  // genlibDLLCall(s_object* ); 
GENLIBDLL_DLL_API(s_object*) genlibDLLCall2(SEXP* ); //genlibDLLCall2(s_object* ); 

#ifdef __cplusplus
}
#endif

#endif //S_GENLIBDLL_INCLUDED_
