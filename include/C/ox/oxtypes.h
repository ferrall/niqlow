/*--------------------------------------------------------------------------
 * oxtypes.h - definitions/declarations for Ox types, requires jdtypes.h
 *
 *       (C) Jurgen Doornik 1992-1995
 *
 *--------------------------------------------------------------------------*/

#ifndef FILE
#include <stdio.h>
#endif
#ifndef JDCALL
#include "jdsystem.h"
#endif
#ifndef FLAG0
#include "jdtypes.h"
#endif

#ifndef INC_OXTYPES
#define INC_OXTYPES

#ifdef __cplusplus
extern "C" {
#endif

/*========================== calling convention ============================*/
#define OXCALL  JDCALL
#define OXCALLC JDCALLC

/*========================== enums/constants ===============================*/
#define OXTYPEBITS 0x00ff                            /* righmost eight bits */
#define OXFLAGBITS (~OXTYPEBITS)
#define OXTYPEBITCOUNT  8                             		  /* eight bits */

enum OxTypes                                              /* variable types */
{   OX_INT       = 1,
    OX_DOUBLE    = 2,
    OX_MATRIX    = 3,
    OX_STRING    = 4,
    OX_ARRAY     = 5,
    OX_FUNCTION  = 6,
    OX_CLASS     = 7,
    OX_VECTOR    = 8,
    OX_INTFUNC   = 9,
    OX_RANGE     = 10,
    OX_FILE      = 11,
    OX_IMPORT    = 12,								   /* for link purposes */
	OX_LAMBDA	 = 13,
	OX_BLOB		 = 14,
	OX_RETURN	 = 64,

	OX_NULL      = FLAG8,                       /* not yet assigned a value */
    OX_VALUE     = FLAG9,			 /* run-time only: owns memory of value */
    OX_CONST     = FLAG10,                                   /* is constant */
	OX_RESERVED  = FLAG11,
	OX_EXTERN    = FLAG12,  	/* is definition only, is error at run-time */
    OX_GLOBAL    = FLAG13,                          /* is a global variable */
    OX_STATDECL  = FLAG14,                         /* is static declaration */
    OX_INLINE    = FLAG15,								 /* inline function */
	  OX_KEYWORD = FLAG15,				   /* Gauss uses inline for keyword */
    OX_MEMBER    = FLAG16,								 /* is class member */
	OX_STATIC    = FLAG17,                           	/* is static member */
    OX_VIRTUAL   = FLAG18,						   /* is a virtual function */
    OX_PUBLIC    = FLAG19, 			             /* is a public data member */
	OX_INDEX     = FLAG20,                                   /* is an index */
    OX_ADDRESS   = FLAG21,								   /* is an address */
    OX_ARGUMENT  = FLAG22,                   /* parser only: is an argument */
    OX_VARARGS   = FLAG23,   /* function has variable no of arguments (...) */
    OX_SERIAL    = FLAG24,/* is a serial variable (parser only) or function */
    OX_VECMAT    = FLAG25,		     /* run-time only: left index is matrix */
    OX_VECRANGE  = FLAG26,		      /* run-time only: left index is range */
    OX_IDXSCALAR = FLAG27, 		      /* run-time only: a scalar is indexed */
	OX_INTERNAL  = FLAG28	 /* parser only: Ox stdlib function declaration */
};

/*============================== typedefs ==================================*/
typedef struct oxarray
{   int    size;
    struct oxvalue *data;
} OxARRAY;

typedef struct oxstring
{   int    size;
    char  *data;
} OxSTRING;

typedef struct oxclass
{   int    size;                                  /* number of data members */
    int    base;        /* -1 or index of base class in global symbol table */
    int    self;			      /* index of itself in global symbol table */
    int    cfunc;								     /* number of functions */	
    struct oxsymbol *data;                              /* data member list */
    struct oxsymbol *func;          /* function list, ends in ->symbol NULL */
} OxCLASS;

typedef struct oxrange
{   int    iMin, iMax;                                             /* range */
} OxRANGE;

typedef struct oxindex			   /* used for storing left index in matrix */
{	union      /* type will be OX_VECTOR, with flags OX_VECMAT, OX_VECRANGE */
    {	int    index;                           /* integer index is default */
        struct oxrange   rval;          /* range index, flag is OX_VECRANGE */
		struct oxvalue	*pval;           /* matrix index, flag is OX_VECMAT */
    } t;
} OxINDEX;

typedef struct oxmatrix
{   int    r, c;                                    /* dimensions of matrix */
    MATRIX data;
	struct oxindex  idx;
} OxMAT;

typedef struct oxfunctioninfo	   /* allocated separately, because too big */
{
    char  *desc;                 /* arg_types$function_name@external_link\0 */
           /* arg_types is 'c' or ' ', one entry for each declared argument */
          /* $ indicates end of argument string, @external_link is optional */
	int	   iFileNameLit;		  /* location of file name in literal stack */
    struct oxsymbol *args;      	               /* NULL or argument list */
	
} OxFUNCTIONINFO;

typedef struct oxfunction
{   int    cargs;/* declared argument count; negative indicates last was ...*/
    union                                                          /* value */
    {
	    int    data;                             /* entry into pcode buffer */
	    void (OXCALL * pFunc)(struct oxvalue *,/* internal function:address */
			struct oxvalue *, int);
	} code;
	struct oxclass 	*clthis;/* runtime, if !NULL:function is member of this */
	struct oxfunctioninfo *pinfo;		   /* ptr to relevant function info */
} OxFUNCTION;

typedef struct oxlambda
{   int    cargs;/* declared argument count; negative indicates last was ...*/
    int    data;	                             /* entry into pcode buffer */
	int	   bpstart;
	int    ap, bp;							 /* context of holding function */
} OxLAMBDA;

typedef struct oxfile
{   FILE  *fp;                                              /* file pointer */
    int    fmode;                             /* internal status/type flags */
    int    *prc;     /* ptr to matrix dimension (known for some file types) */
	void   *zip;					 /* opaque ptr to zip related structure */
} OxFILE;

typedef struct oxreturn			         /* return info after function call */
{   int    ap, bp;										/* only runtime use */
	int	   cpcode;
} OxRETURN;

typedef struct oximport
{   char  *name;							                 /* import name */
	int	   status;										   /* import status */
	int	   type;											 /* import type */
    char  *path;							                 /* import path */
} OxIMPORT;

typedef struct oxblob			                       /* storage container */
{   int i1;
	int i2;
	union
    {	int    ival;
        double dval;
		INT64  i64val;
        void   *pval;
		struct oxvalue *oxval;
    } t1;
    union
    {	int    ival;
        double dval;
		INT64  i64val;
        void   *pval;
		struct oxvalue *oxval;
    } t2;
    union
    {	int    ival;
        double dval;
		INT64  i64val;
        void   *pval;
		struct oxvalue *oxval;
    } t3;
} OxBLOB;

typedef struct oxvalue                        /* value descriptor structure */
{   int type;                                              /* variable type */
    union                                                          /* value */
    {
        double dval;                                              /* double */
        int    ival;                                             /* integer */
        struct oxstring    sval;                                  /* string */
        struct oxmatrix    mval;                                  /* matrix */
        struct oxrange     rval;                                   /* range */
        struct oxarray     aval;                                   /* array */
        struct oxfunction  fval;                                /* function */
        struct oxclass     clval;                                  /* class */
        struct oxfile      ioval;                              /* i/o files */
		struct oxreturn	   retval;			        /* function call return */
		struct oximport	   impval;						     /* import item */
		struct oxlambda    lambda;						 /* lambda function */
		struct oxblob      blob;								    /* blob */
    } t;
} OxVALUE;

typedef struct oxsymbol                                     /* symbol table */
{   char  *symbol;                                 /* points to symbol name */
    struct oxvalue val;                                  /* value of result */
} OxSYMBOL;

typedef void (OXCALL * OxFUNCP)(OxVALUE *, OxVALUE *, int);


/*========================== macro's =======================================*/
#define HASFLAG(i,fl)      ( (i) & (fl) )
#define ADDFLAG(i,fl)      ( (i) |= ((fl) & OXFLAGBITS) )
#define CLEARFLAG(i,fl)    ( (i) &= (OXTYPEBITS | ~(fl)) )
#define COPYFLAGS(i,i0,fl) ( (i) = (((i) & (OXTYPEBITS | ~(fl))) | ((i0) & (fl))))
#define HASTYPE(i,tp)      ( ((i) & OXTYPEBITS) == (tp) )
#define GETTYPE(i)         ( (i) & OXTYPEBITS)
#define SETTYPE(i,tp)      ( (i) = ((i) & OXFLAGBITS) | (tp) )

#define GETPVTYPE(pv)     ( GETTYPE((pv)->type) )
#define SETPVTYPEFLAGS(pv,tp,fl)   	( (pv)->type = (tp) | (fl) )

#define ISINT(pv)         ( GETTYPE((pv)->type) == OX_INT )
#define ISDOUBLE(pv)      ( GETTYPE((pv)->type) == OX_DOUBLE )
#define ISMATRIX(pv)      ( GETTYPE((pv)->type) == OX_MATRIX )
#define ISSTRING(pv)      ( GETTYPE((pv)->type) == OX_STRING )
#define ISARRAY(pv)       ( GETTYPE((pv)->type) == OX_ARRAY )
#define ISFUNCTION(pv)    ( GETTYPE((pv)->type) == OX_FUNCTION )
#define ISCLASS(pv)       ( GETTYPE((pv)->type) == OX_CLASS )
#define ISVECTOR(pv)      ( GETTYPE((pv)->type) == OX_VECTOR )
#define ISINTFUNC(pv)     ( GETTYPE((pv)->type) == OX_INTFUNC )
#define ISRANGE(pv)       ( GETTYPE((pv)->type) == OX_RANGE )
#define ISFILE(pv)        ( GETTYPE((pv)->type) == OX_FILE )
#define ISIMPORT(pv)      ( GETTYPE((pv)->type) == OX_IMPORT )
#define ISLAMBDA(pv)      ( GETTYPE((pv)->type) == OX_LAMBDA )
#define ISBLOB(pv)        ( GETTYPE((pv)->type) == OX_BLOB )

#define ISNULL(pv)        ( HASFLAG((pv)->type, OX_NULL) )
#define ISVALUE(pv)       ( HASFLAG((pv)->type, OX_VALUE) )
#define ISCONST(pv)       ( HASFLAG((pv)->type, OX_CONST) )
#define ISEXTERN(pv)      ( HASFLAG((pv)->type, OX_EXTERN) )
#define ISGLOBAL(pv)      ( HASFLAG((pv)->type, OX_GLOBAL) )
#define ISSTATDECL(pv)    ( HASFLAG((pv)->type, OX_STATDECL) )
#define ISINLINE(pv)      ( HASFLAG((pv)->type, OX_INLINE) )
#define ISMEMBER(pv)      ( HASFLAG((pv)->type, OX_MEMBER) )
#define ISSTATIC(pv)      ( HASFLAG((pv)->type, OX_STATIC) )
#define ISVIRTUAL(pv)     ( HASFLAG((pv)->type, OX_VIRTUAL) )
#define ISPUBLIC(pv)      ( HASFLAG((pv)->type, OX_PUBLIC) )
#define ISINDEX(pv)       ( HASFLAG((pv)->type, OX_INDEX) )
#define ISADDRESS(pv)     ( HASFLAG((pv)->type, OX_ADDRESS) )
#define ISVARARGS(pv)     ( HASFLAG((pv)->type, OX_VARARGS) )
#define ISSERIAL(pv)      ( HASFLAG((pv)->type, OX_SERIAL) )
#define ISVECMAT(pv)      ( HASFLAG((pv)->type, OX_VECMAT) )
#define ISVECRANGE(pv)    ( HASFLAG((pv)->type, OX_VECRANGE) )
#define ISIDXSCALAR(pv)   ( HASFLAG((pv)->type, OX_IDXSCALAR) )
#define ISINTERNAL(pv)    ( HASFLAG((pv)->type, OX_INTERNAL) )

#define ISSCALAR(pv)      (ISINT(pv) || ISDOUBLE(pv))
#define ISSCALMAT(pv)     (ISMATRIX(pv) && (pv)->t.mval.r == 1 && (pv)->t.mval.c == 1)
#define ISNONSCALMAT(pv)  (ISMATRIX(pv) && ((pv)->t.mval.r > 1 || (pv)->t.mval.c > 1))
#define MATISSCALAR(pv)   ((pv)->t.mval.r == 1 && (pv)->t.mval.c == 1)
#define MATISSCALEMPTY(pv)((pv)->t.mval.r <= 1 && (pv)->t.mval.c <= 1)
#define MATISEMPTY(pv)    ((pv)->t.mval.r == 0 || (pv)->t.mval.c == 0)
#define ISARITHMETIC(pv)  (ISINT(pv) || ISDOUBLE(pv) || ISMATRIX(pv))


#define OxInt(pv,i)       (pv)[i].t.ival
#define OxDbl(pv,i)       (pv)[i].t.dval
#define OxDouble(pv,i)    (pv)[i].t.dval
#define OxMat(pv,i)       (pv)[i].t.mval.data
#define OxMatc(pv,i)      (pv)[i].t.mval.c
#define OxMatr(pv,i)      (pv)[i].t.mval.r
#define OxMatrc(pv,i)     ((VECIDX)((pv)[i].t.mval.r) * (VECIDX)((pv)[i].t.mval.c))
#define OxStr(pv,i)       (pv)[i].t.sval.data
#define OxStrLen(pv,i)    (pv)[i].t.sval.size
#define OxArray(pv,i)     (pv)[i].t.aval.data
#define OxArrayLen(pv,i)  (pv)[i].t.aval.size
#define OxZero(pv,i)      (pv)[i].type = OX_INT, (pv)[i].t.ival = 0
#define OxSetInt(pv,i,j)  (pv)[i].type = OX_INT, (pv)[i].t.ival = j
#define OxSetDbl(pv,i,d)  (pv)[i].type = OX_DOUBLE, (pv)[i].t.dval = d
#define OxSetDouble(pv,i,d) (pv)[i].type = OX_DOUBLE, (pv)[i].t.dval = d
#define OxSetMatPtr(pv,i,m,cr,cc) (pv)[i].type = OX_MATRIX, (pv)[i].t.mval.r = cr, (pv)[i].t.mval.c = cc, (pv)[i].t.mval.data = m
#define OxSetStrPtr(pv,i,s,cs) (pv)[i].type = OX_STRING, (pv)[i].t.sval.data = s, (pv)[i].t.mval.size = cs
#define OxSetAddress(pv,i,pva) (pv)[i].type = OX_ARRAY | OX_ADDRESS, OxArraySize(pv+i) = 1, OxArrayData(pv+i) = pva

#define OxArraySize(pv)   ((pv)->t.aval.size)
#define OxArrayData(pv)   ((pv)->t.aval.data)
#define OxStringSize(pv)  ((pv)->t.sval.size)
#define OxStringData(pv)  ((pv)->t.sval.data)
#define OxClassSize(pv)   ((pv)->t.clval.size)
#define OxFuncReqArgs(pv) ((pv)->t.fval.cargs)
#define OxFuncVa_args(pv) (ISVARARGS(pv))
#define OxFuncInfo(pv)	  ((pv)->t.fval.pinfo)
#define OxFuncDesc(pv)	  ((pv)->t.fval.pinfo->desc)
#define OxFuncName(pv)	  ((pv)->t.fval.pinfo->desc + OxFuncReqArgs(pv) + 1)
#define OxFuncAddress(pv) ((pv)->t.fval.code.pFunc)
#define OxFuncPcode(pv)   ((pv)->t.fval.code.data)
#define OxFuncFileNameLit(pv) ((pv)->t.fval.pinfo->iFileNameLit)
#define OxFuncArgList(pv) ((pv)->t.fval.pinfo->args)
#define OxImportName(pv)  ((pv)->t.impval.name)


#ifdef __cplusplus
}
#endif

#endif  /* INC_OXTYPES */


