/*--------------------------------------------------------------------------
 * jdsystem.h - definitions/declarations for system dependent functions.
 *
 *       (C) Jurgen Doornik 1994-2003
 *
 *--------------------------------------------------------------------------*/

#ifndef INC_JDSYSTEM
#define INC_JDSYSTEM

#ifdef __cplusplus
extern "C" {
#endif

enum DllErrors
{
	DLL_OK = 0, DLL_NOLOADLIB, DLL_NODYNFUNC, DLL_NODYNAMIC
};	
	
/*======================== Unicode related settings ========================*/
/* define JDUNICODE if UNICODE or _UNICODE is defined, and vice versa */
#if defined(_UNICODE) || defined(UNICODE)
  #ifndef JDUNICODE
    #define JDUNICODE
  #endif
#else /* !UNICODE */
  #undef JDUNICODE
#endif /* UNICODE/!UNICODE */

#if defined(JDUNICODE)
  #ifndef _UNICODE
    #define _UNICODE
  #endif
  #ifndef UNICODE
    #define UNICODE
  #endif
#endif /* JDUNICODE */

/*======================= platform specific settings =======================*/
#if defined(_WINDOWS) || defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(WIN64) || defined(_WIN64) || defined(__WIN64__)

  /* _M_WIN32 for any 32 AND 64 bit Windows, _M_WIN64 only for 64 bit */
  #define _M_WIN32
  #if defined(WIN64) || defined(_WIN64) || defined(__WIN64__)
    #define _M_WIN64
  #endif

  /* Windows compiler specific */
  #if defined(_MSC_VER)

    #define JDINLINE 	__inline
    #define JDCALL 		__stdcall
    #define JDCALLC 	__cdecl

    /* for now, no warnings for `deprecated' C standard library functions */
    #if _MSC_VER >= 1300
      #pragma warning(disable : 4996)
    #endif

  #elif defined(__MINGW32__)

    #define JDINLINE 	inline
    #define JDCALL		__stdcall
    #define JDCALLC 	__cdecl

  #elif defined (__BORLANDC__)

    #define JDINLINE
    #define JDCALL 		__stdcall
    #define JDCALLC 	__cdecl
	#define TCHAR _TCHAR

  #else

    #define JDINLINE 	
    #define JDCALL		__stdcall
    #define JDCALLC		__cdecl

  #endif /* endif Windows compiler specific */

  #define JDUSE_TCHAR_H
  #define EOL_STR "\r\n"                             /* end of line marking */
  #define JDPATH_SEP '\\'
  #define JDDRIVE_SEP ':'
  #define JDENV_SEP ';'
  #define JDDLL_EXT ".dll"
  #define SYSTEM_ID "Windows"				    /* NB: no _T() in this file */
  #if defined(_M_WIN64)
	#define JD64BIT
  #endif

#else

  /* defaults for unix-based systems */
  #define JDIS_UNIX
  #if defined(__x86_64__) || defined(_LP64)
    #define JD64BIT
  #endif
  #define JDINLINE
  #define JDCALL
  #define JDCALLC
  #define JDINLINE
  #define stricmp strcasecmp
  #define strnicmp strncasecmp
  #define JDUSE_WCHAR_H								/* no tchar.h available */
  #define EOL_STR "\n"                               /* end of line marking */
  #define JDPATH_SEP '/'
  #define JDDRIVE_SEP ':'             /* no drive separator, use : as dummy */
  #define JDENV_SEP ':'
  #define JDGETTIMEOFDAY	     /* use gettimeofday for more accurate time */
  
  /* specific unix-like systems */
  #if defined(__linux__)
  
    /* GCC for Linux on i386 or X86-64 */
    #define _M_LINUX
  
    #define SYSTEM_ID "Linux"
    #define JDDLL_EXT ".so"
    #define JD_DLOPEN
	#define JDUSE_GLIBC2_13
	
	#if defined(JDUSE_GLIBC2_13)
      #define JDCLOCK_GETTIME_SYSCALL			/* avoid GLIBC 2.17 for now */
	#else
      #define JDCLOCK_GETTIME			    		   /* use clock_gettime */
	#endif

  #elif defined(_M_OSX)                       				    
  
    /* Mac OS-X on Intel, gcc or clang */
    #define SYSTEM_ID "OS_X"
    #define JDDLL_EXT ".so"	   /* use .so for Ox plugins, .dylib for ox DLL */
    #define JDDLL_SUFFIX "_osx"	           /* default suffix for Ox plugins */
    #define JD_DLOPEN
  
    #define wcscasecmp  WcsCaseCmp
    #define wcsncasecmp WcsNCaseCmp

	#ifdef _NONSTD_SOURCE
      /* BSD/OSX define _T() as something different in ctype.h, override it */
      #include <ctype.h>
      #undef _T
  	#endif
	
  #elif defined(__asmjs__)
    /* Emscripten targeting asm.js */
    #define _M_ASMJS
  
    #define SYSTEM_ID "asm.js"
    #define JDDLL_EXT ".js"
    #define JDDLL_UNSUPPORTED
    #define JD_SANDBOX
    #define JDCLOCK_GETTIME							     /* does this work? */
  
    #define wcscasecmp  wcscmp
    #define wcsncasecmp wcsncmp
  
  #elif defined(__pnacl__)
    /* portable native client */
    #define _M_PNACL
  
    #define SYSTEM_ID "PNaCl"
    #define JDDLL_EXT ".so"
    #define JDDLL_UNSUPPORTED
    #define JD_SANDBOX
    #define JD_NODIRENT_DTYPE
  
  #elif defined(__sun)
    /* Solaris on Sparc or x86, neither currently supported */
    #define _M_SUN
  
    #include <sys/types.h>							 /* for _LP64 or _ILP32 */
    #if defined(__sparc)
      #define JDBIG_ENDIAN                  /* this is a big endian machine */
    #endif
  
    #define SYSTEM_ID "Sun"
    #define JDDLL_EXT ".so"
    #if defined(__sparc)
      #define JDDLL_SUFFIX "_sparc"	       /* default suffix for Ox plugins */
    #else
      #define JDDLL_SUFFIX "_sunx86"       /* default suffix for Ox plugins */
    #endif
    #define JD_DLOPEN
  
    #define wcscasecmp  wcscmp
    #define wcsncasecmp wcsncmp
    #define JD_NODIRENT_DTYPE
  
  #else										/* generic ANSI-C Unix platform */
  
    #define JDDLL_EXT ".so"
    #define JDDLL_SUFFIX "_gen"	           /* default suffix for Ox plugins */
    #define JD_DLOPEN
    #define SYSTEM_ID "Generic"
  
  #endif  /* specific unix-like systems */
#endif


#ifndef max
#define max(a,b)  (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
#define min(a,b)  (((a) < (b)) ? (a) : (b))
#endif

/*===================== platform specific 64-bit integer ===================*/
/* 64-bit Linux uses LP64:  int=4 long=pointer=8 long long=8 bytes
   64-bit Windows uses P64:	int=long=4 pointer=8 long long=8
   64-bit uses ILP32:		int=long=pointer=4 long long=8
*/
#ifdef __LP64__
	typedef unsigned long UINT64;
	typedef          long INT64;
	#define LIT_UINT64(c) (c##ul)
	#define LIT_INT64(c)  (c##l)
#elif defined(_MSC_VER) || defined (__BORLANDC__)
	typedef unsigned __int64 UINT64;
	typedef          __int64 INT64;
	#define LIT_UINT64(c) (c##ui64)
	#define LIT_INT64(c)  (c##i64)
#else 
	typedef unsigned long long UINT64;
	typedef          long long INT64;
	#define LIT_UINT64(c) (c##ull)
	#define LIT_INT64(c)  (c##ll)
#endif

#ifdef __cplusplus
}
#endif

#endif  /* INC_JDSYSTEM */
                             
