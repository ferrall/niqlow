/*--------------------------------------------------------------------------
 * jdtypes.h - general typedefs and constants
 *
 *       (C) Jurgen Doornik 1990-1997
 *
 *--------------------------------------------------------------------------*/

/*------------------ revision history: -------------------------------------
 * 17-15-95: Only defining M_ if not already existing
 * 20-10-94: Changed definition of SQR, as x*x instead of pow(x,2)
 *  5-08-93: Changed definition of M_PI, added some functions of pi.
 *--------------------------------------------------------------------------*/

#include <stdlib.h>                                     /* some need size_t */
#ifndef JDCALL
#include "jdsystem.h"
#endif


#ifndef INC_JDTYPES
#define INC_JDTYPES

#ifdef __cplusplus
extern "C" {
#endif

/*=========================== typedefs =====================================*/
/* general typedefs                  prefix Hungarian Notation*/
typedef unsigned char   	BYTE;                                      /* b */
#define BOOL jdbool
typedef int    		  		jdbool;                                  /* f,b */
#ifdef JD_LDOUBLE              /* on some machines long double is very slow */
  typedef long double     	LDOUBLE;                                 /* ld  */  
#else										  /* LDOUBLE defaults to double */
  typedef double          	LDOUBLE;  								 /* ld  */ 
#endif
typedef unsigned int    	UINT;     								 /* ui  */
typedef int  (JDCALL * PFI)(void);	 /* ptr to function returning int, pfi  */
typedef double **			MATRIX; 							   /* m,mat */
typedef double *			VECTOR;						           /* v,vec */
typedef int    **			INTMAT;							          /* im */
typedef int    *			INTVEC;							          /* iv */
#if defined(JD64BIT)
  typedef INT64				VECIDX;			 /* index for vectorized MATRIX */
#else
  typedef int				VECIDX;			 /* index for vectorized MATRIX */
#endif

typedef struct
{
	int  year1, period1, year2, period2, freq;
} SAMPLE;

typedef struct
{
	BOOL	isdated;
	union
	{	SAMPLE 	sam;
		double  dates[3];						    /* min, max, obsperweek */
	} db;
	union
	{	SAMPLE 	sam;
		double  dates[2];
	} sel;
} SUBSAMPLE;


#if defined(JDUSE_TCHAR_H)
  #include <tchar.h>
#elif defined(JDUSE_WCHAR_H)
  #if defined(JDUNICODE)
    #include <wchar.h>
    typedef wchar_t		TCHAR;
  #else
    typedef char 		TCHAR;
  #endif
#elif !defined(JDUNICODE)
    typedef char 		TCHAR;
#endif


#ifndef FALSE
#define  FALSE  0
#define  TRUE   1
#endif

#define FREQ_ANNUAL     1
#define FREQ_QUARTERLY  4
#define FREQ_MONTHLY   12
#define FREQ_WEEKLY52  52

                     /* Bit flags: FLAGi is that value with bit number i on */
#ifndef FLAG0
#define  FLAG0      0x0001
#define  FLAG1      0x0002
#define  FLAG2      0x0004
#define  FLAG3      0x0008
#define  FLAG4      0x0010
#define  FLAG5      0x0020
#define  FLAG6      0x0040
#define  FLAG7      0x0080
#define  FLAG8      0x0100
#define  FLAG9      0x0200
#define  FLAG10     0x0400
#define  FLAG11     0x0800
#define  FLAG12     0x1000
#define  FLAG13     0x2000
#define  FLAG14     0x4000
#define  FLAG15     0x8000
#define  FLAG16   0x010000
#define  FLAG17   0x020000
#define  FLAG18   0x040000
#define  FLAG19   0x080000
#define  FLAG20   0x100000
#define  FLAG21   0x200000
#define  FLAG22   0x400000
#define  FLAG23   0x800000
#define  FLAG24  0x1000000
#define  FLAG25  0x2000000
#define  FLAG26  0x4000000
#define  FLAG27  0x8000000
#define  FLAG28 0x10000000
#define  FLAG29 0x20000000
#define  FLAG30 0x40000000
#define  FLAG31 0x80000000
#endif

#define FLT_MIN_E_EXP  (-82)                              /* min e exponent */
#define FLT_MAX_E_EXP    88                               /* max e exponent */
#define DBL_MIN_E_EXP (-706)                              /* min e exponent */
#define DBL_MAX_E_EXP   709                               /* max e exponent */

#ifndef M_PI
#define M_PI        3.1415926535897932384626433832795028841972        /* pi */
#endif
#ifndef M_2PI
#define M_2PI       6.2831853071795864769252867665590057683944      /* 2*pi */
#endif
#ifndef M_SQRT2PI
#define M_SQRT2PI   2.50662827463100050241576528481104525300700/* sqrt(2*pi)*/
#endif
#ifndef M_PI_2
#define M_PI_2      1.5707963267948966192313216916397514420986      /* pi/2 */
#endif
#ifndef M_1_PI
#define M_1_PI      0.3183098861837906715377675267450287240689      /* 1/pi */
#endif
#ifndef M_E
#define M_E         2.71828182845904523543                             /* e */
#endif
#ifndef M_EULER
#define M_EULER     0.57721566490153286061       /* Euler's constant, gamma */
#endif
#ifndef M_SQRT2INV
#define M_SQRT2INV  0.70710678118654752440						/* 1/sqrt(2)*/
#endif

#ifndef SQR
#define  SQR(x) ((x) * (x))
#endif
#ifndef ROUND
#define  ROUND(x) floor(x + 0.5)
#endif

#define M_TICKS_PER_DAY (24 * 60 * 60 * 100.0)

enum JdShowModes
{	SM_ACTIVATE_VIEW 	= FLAG0,
	SM_ACTIVATE_SERVER 	= FLAG1,
	SM_UPDATE 			= FLAG2,
	SM_ACTIVATE_LINK 	= FLAG8
};
enum JdTextModes
{	TEXT_ASCII,
	TEXT_UTF8,
	TEXT_UTF16LE,
	TEXT_UTF16BE,
	TEXT_UTF32LE,
	TEXT_UTF32BE
};

/* jdsystem.c */
double JDCALL DInf(void);
double JDCALL DNaN(void);
BOOL   JDCALL FIsInf(double d);
BOOL   JDCALL FIsInfOrNaN(double d);
BOOL   JDCALL FIsNaN(double d);
void   JDCALL SetInf(double *pd);
void   JDCALL SetNaN(double *pd);

#ifdef __cplusplus
}
#endif

#endif  /* INC_JDTYPES */

