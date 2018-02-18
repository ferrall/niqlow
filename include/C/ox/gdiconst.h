/*--------------------------------------------------------------------------
 * gdiconst.h - constants for the GDI plot library
 *
 *       Jurgen Doornik 1994-98
 *--------------------------------------------------------------------------*/

#ifndef INC_GDICONST
#define INC_GDICONST

#ifdef __cplusplus
extern "C" {
#endif

/*=========================== gdi constants ================================*/
#define GDI_CXPIXEL       15000                       /* bottom left is 0,0 */
#define GDI_CYPIXEL       10000
#define GDI_MX_AREA          36
#define GDI_MN_WIN            3
#define GDI_MX_WIN          256

#define GDI_BLOCKSIZE       256
#define GDI_MXAUTOLINES      16
#define GDI_MXAUTOMODES		  4
#define GDI_MXFONTS          10
#define GDI_MXPALETTE	      8

enum GDImodes									  /* display (screen) modes */
{
	GDI_MODE_COLOR,							   /* colour always the default */
	GDI_MODE_BWG
};
enum GDIprintmodes						 /* print modes: four possibilities */
{
	GDI_PRINT_BW,
	GDI_PRINT_BWG,
	GDI_PRINT_GRAY,
	GDI_PRINT_COLOR								  /* colour now the default */
};
enum GDIcodes
{
    GDI_NOP,                                 /* value 0: no operation: skip */
    GDI_MOVETO,                                   /* x: 0-10000, y: 0-10000 */
	GDI_RMOVETO,
    GDI_LINETO,                                   /* x: 0-10000, y: 0-10000 */
	GDI_RLINETO,
	GDI_BOXTO,
    GDI_LINETYPE,                                    /* type: 0-10, 0=solid */
    GDI_LINEWIDTH,                                        /* width: 0-10000 */
    GDI_COLOR,                                             /* colour: 0-256 */
    GDI_COLORFILL,                   /* fill colour in RGB format, -1: none */
    GDI_TEXTDIRECTION,                                  /* direction: 0-360 */
    GDI_TEXTFONT,                                 /* text font number, size */
    GDI_TEXT,                                            /* text: ends in 0 */
	GDI_TEXTTO,                     /* text: ends in 0, update the position */
	GDI_TEXTTOSYM,               /* character from symbol font, index, size */
	GDI_TEXTTOCH,	     	          /* character from current font, index */
	GDI_TEXTTABLE,       /* table mode, 0=start, 1= new line, 1=next column */
    GDI_SYMBOL,            /* symbol type (GDI_SYMBOL is not emitted):      */
        GDI_FILLBOX,                                          /* filled box */
        GDI_BOX,                                                /* open box */
        GDI_PLUS,                                                   /* plus */
        GDI_DASH,                                                   /* dash */
        GDI_CIRCLE,                                               /* circle */
        GDI_NOSYMBOL,
        GDI_FILLCIRCLE,
        GDI_TRIANGLE,
        GDI_FILLTRIANGLE,
        GDI_DIAMOND,
        GDI_FILLDIAMOND,
        GDI_CROSS,
        GDI_TRIANGLEUP,
        GDI_FILLTRIANGLEUP,
        GDI_CNTSYMBOL,						    /* upper bound, not emitted */
    GDI_BEGINPATH,						                      /* begin path */
    GDI_ENDPATH,               /* end path, stroke, fill or stroke and fill */
	GDI_FILLBOXTO,
    GDI_BEGINTEXT,						                /* begin text block */
    GDI_ENDTEXT,							              /* end text block */
    GDI_ENDS                                                 /* end of code */
};
enum GDIsymbolCodes
{
    PL_FILLBOX,                                               /* filled box */
    PL_BOX,                                                     /* open box */
    PL_PLUS,                                                        /* plus */
    PL_DASH,                                                        /* dash */
    PL_CIRCLE,                                                      /* dash */
    PL_NOSYMBOL,			/* leave at value 5 for backwards compatibility */
	PL_LINE = PL_NOSYMBOL,
    PL_FILLCIRCLE,
    PL_TRIANGLE,
    PL_FILLTRIANGLE,
    PL_DIAMOND,
    PL_FILLDIAMOND,
    PL_CROSS,
    PL_TRIANGLEUP,
    PL_FILLTRIANGLEUP,
	PL_CNTSYMBOL
};
enum GDItypeCodes
{
	TP_SOLID,
	TP_DOTTED,
	TP_DASHED,
	TP_LDASHED,
	TP_USER
};
enum GDIlinetypeCodes
{
	LTP_LINE,
	LTP_SYMBOLS,
	LTP_LINESYMBOLS,
	LTP_INDEX,
	LTP_INDEXSYMBOLS,
	LTP_BARS,
	LTP_SHADING,
	LTP_MAX = LTP_SHADING
};
enum GDIpathCodes
{
	PATH_STROKE,
	PATH_FILL,
	PATH_STROKE_FILL
};	
enum GDItextTableCodes
{
	COL_INIT = 0, COL_TEST = 1, COL_EXIT = 2, COL_NEXT = 3
};
/*=========================== gdidraw constants ============================*/
enum GDIaxisScaleCodes
{
	AXIS_LINEAR,
	AXIS_LOG,
	AXIS_LOG10,
	AXIS_SCALED,
	AXIS_DATE
};
enum TextTypes
{
	TEXT_TEXT,
	TEXT_TITLE,
	TEXT_XLABEL,
	TEXT_YLABEL,
	TEXT_ZLABEL
};

enum DrawTypes                                              /* draw objects */
{
    DRAW_NONE     = 0,
    DRAW_XVEC     = 1,
    DRAW_AXIS     = 2,
	DRAW_LINE     = 4,
    DRAW_HIST     = 5,
    DRAW_AREA     = 6,
    DRAW_TEXT     = 7,
    DRAW_LEGEND   = 8
};
enum AnchorTypes
{	ANCHOR_MIN, ANCHOR_MAX, ANCHOR_USER
};
enum ZModeTypes
{	ZMODE_NONE = 0, ZMODE_SYMBOL, ZMODE_VALUE, ZMODE_BAR, ZMODE_BAND,
	ZMODE_HILO, ZMODE_3D, ZMODE_FAN,
	ZMODE_MESH = 256, ZMODE_TRIANGLE
};
enum ContourTypes
{	CONTOUR_NONE = 0, CONTOUR_ADD, CONTOUR_2D
};

enum DrawSetOptions
{
	SET_DEFAULT,		SET_MARGIN,			SET_BOX,
	SET_FONT,			SET_AXIS,			SET_AXISLINE,
	SET_GRID,			SET_LEGEND,			SET_LEGENDHIDE,
	SET_HISTOGRAM,		SET_COLOR,			SET_SYMBOL,
	SET_LINE,			SET_COLORMODEL,		SET_PRINTPAGE,
	SET_PALETTE_MAX,	SET_PALETTE_MIN,	SET_XYSTYLE,
	SET_PAPERCOLOR,		SET_BWG,			SET_LINEBWG,
	SET_AXISZEROS,		SET_AXISFORMAT,		SET_LEGENDRESIZE,
	SET_LEGENDFONTSIZE, SET_LEGENDCOLOR
};

enum DrawSetAdjustment
{
	ADJ_COLOR,		ADJ_SYMBOL, 	ADJ_INDEX,
	ADJ_SCALE,		ADJ_MINMAX, 	ADJ_AXISHIDE,
	ADJ_AXISLINE,	ADJ_AXISGRID,	ADJ_SYMBOLUSE,
	ADJ_AREA_X,     ADJ_AREA_Y,     ADJ_AREA_P,
	ADJ_AREAMATRIX, ADJ_TEXTANGLE,	ADJ_AXISLABEL,
	ADJ_AREA_3D,	ADJ_3D_ROTATE,	ADJ_3D_RESET,
	ADJ_AREA_Z,     ADJ_AXISSCALE,	ADJ_AXISCENTRE,
	ADJ_AREA_3DSET, ADJ_PAPERSCALE, ADJ_LEGEND,
	ADJ_COLORMODEL, ADJ_PAPERCOLOR, ADJ_AREASCOLOR,
	ADJ_REGRESSION
};


#ifdef __cplusplus
}
#endif

#endif  /* INC_GDICONST */

