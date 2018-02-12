/*--------------------------------------------------------------------------
 * jdmatrix.h - definitions/declarations for matrix allocation, and basic
 *              initialization.
 *
 *       Jurgen Doornik 1990-1998
 *
 *--------------------------------------------------------------------------*/

#ifndef INC_JDMATRIX
#define INC_JDMATRIX

#ifdef __cplusplus
extern "C" {
#endif

#define MatAlloc        MatAllocBlock
#define CMatAlloc       CMatAllocBlock
#define MatFree(m,r,c)  MatFreeBlock(m)
#define VecAlloc(dim)   ((VECTOR)malloc( (dim) * sizeof(double)))
#define VecFree(vX)     free(vX)
#define VecDup(vSrc,cX) VecCopy(VecAlloc(cX), vSrc, cX)

#ifdef JDDEBUG_
#include "jddebug.h"		/* replaces matrix allocators by debug versions */
#endif

typedef struct MatrixCache
{
	int 		 cCacheSize;
	int 		 cCacheMatrixSize;
	MATRIX 		 *cache;	 	                         /* cached matrices */
	int    		 *acRowCache;                       /* no of rows in matrix */
	int    		 *acColCache;                    /* no of columns in matrix */
	unsigned int *auiCacheCtr;                   /* counter for cached item */
	int 		 iCacheFree;             /* entry of known free slot, or -1 */
	unsigned int uiCacheCtr;					      /* cache item counter */
	int  		 cCacheSizeUsed, cCacheMatrixSizeUsed; /* actual dimensions */
	int  		 cCacheHit, cCacheMiss, cCacheFlush;		  /* usage info */
} JDMATCACHE;


/*================== function prototypes & macro's =========================*/
#define IntVecAlloc(dim)      ((INTVEC)malloc( (dim) * sizeof(int)))
#define IntVecFree(aiX)       free(aiX)
#define IntMatAlloc(r, c)     IntMatAllocBlock(r, c)
#define IntMatFree(mat, r, c) IntMatFreeBlock(mat)

/* jdstdlib.c */
MATRIX JDCALL MatCpy(MATRIX mDest, MATRIX mSrc, int cR, int cC);
MATRIX JDCALL MatCopy(MATRIX mDest, MATRIX mSrc, int cR, int cC);
void   JDCALL MatCopyVecr(MATRIX mX, VECTOR vY, int cR, int cC);
void   JDCALL MatCopyVecc(MATRIX mX, VECTOR vY, int cR, int cC);
void   JDCALL VecrCopyMat(VECTOR vDest, MATRIX mSrc, int cR, int cC);
void   JDCALL VeccCopyMat(VECTOR vDest, MATRIX mSrc, int cR, int cC);
void   JDCALL MatSetAt(MATRIX mX, double d, int i, int j);
double JDCALL MatGetAt(MATRIX mX, int i, int j);
MATRIX JDCALL MatDup(MATRIX mSrc, int cR, int cC);
MATRIX JDCALL MatI(MATRIX mX, int cX);
MATRIX JDCALL MatZero(MATRIX mX, int cR, int cC);
MATRIX JDCALL MatNaN(MATRIX mX, int cR, int cC);
VECTOR JDCALL VecAllocBlock(int cX);
void   JDCALL VecFreeBlock(VECTOR vX);
VECTOR JDCALL VecDupBlock(VECTOR vSrc, int cX);
VECTOR JDCALL VecCopy(VECTOR vDest, VECTOR vSrc, int cX);
MATRIX JDCALL MatAllocBlock(int cR, int cC);
void   JDCALL MatFreeBlock(MATRIX mX);
INTMAT JDCALL IntMatAllocBlock(int r, int c);
void   JDCALL IntMatFreeBlock(INTMAT mat);

/* jdstdlib.c - matrix cache */
void   JDCALL MatCacheFlush(JDMATCACHE *pcache);
void   JDCALL MatCacheExit(JDMATCACHE *pcache);
void   JDCALL MatCacheStatus(JDMATCACHE *pcache);
void   JDCALL SetMatCacheSize(int cEntries, int cSize);
void   JDCALL MatCacheInit(JDMATCACHE *pcache);
MATRIX JDCALL CMatAllocBlock(int cR, int cC);
MATRIX JDCALL CMatDup(MATRIX mSrc, int cR, int cC);
void   JDCALL CMatFree(MATRIX mX, int cR, int cC);
void   JDCALL SetMatCacheLocker(JDMATCACHE * (JDCALL * pfLocker)(int lockit));

/* jdsort.c */
int    JDCALL SortVec(VECTOR vX, int cX);
int    JDCALL SortVec2(VECTOR vX, int cX, VECTOR vY);
int    JDCALL SortVec2DropNaN(VECTOR vX, int cX, VECTOR vY);
int    JDCALL SortMatCol(MATRIX mX, int iCol, int cX);
int    JDCALL SortmXtByVec(int cT, VECTOR vBy, MATRIX mXt, int cX);
int    JDCALL SortmXByCol(int iCol, MATRIX mX, int cT, int cX);
void   JDCALL SortmXByColumns(int *piIdx, int cIdx, MATRIX mX, int cT, int cX);
void   JDCALL SortVecByVec(int cT, VECTOR vBy, VECTOR vAc);

#ifdef __cplusplus
}
#endif

#endif  /* INC_JDMATRIX */


