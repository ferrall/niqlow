/* This file is part of niqlow. Copyright (C) 2012 Christopher Ferrall */
#include <oxstd.h>

struct LinePoint
	{
	decl step;
	decl f;
	}

struct  LineMax
   		{
		/** . @internal **/			static 	const 	decl 	tiny             = 1.0e-14;
		/** . @internal **/        	static 	const 	decl 	tolerance        = 1E-3;
		/** . @internal **/        	static 	const 	decl 	glimit           = 10.0;
		/** . @internal **/        	static 	const 	decl 	gold             = 1.61803399;
		/** . @internal **/        	static 	const 	decl 	rgold            = .61803399;
		/** . @internal **/        	static 	const 	decl 	cgold            = 1-.61803399;
		/** . @internal **/        	static 	const 	decl 	maxstp           = 5.0;

		protected:
				const	decl	p1,p2,p3,p4,p5,p6;
						decl 	O;
						decl	maxiter,improved;
						decl 	Delta;
						decl    q,a,b;
		
		LineMax();
		~LineMax();
		Iterate(O,Delta,maxiter);
		Try(pt,step);
		Bracket();
		Golden();
		}
