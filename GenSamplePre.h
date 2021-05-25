#pragma once
#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>

#include <NTL/mat_ZZ.h>
#include <time.h>
#include <NTL/LLL.h>

#include <NTL/vec_RR.h>
#include <NTL/RR.h> 
#include <random>
#include "SamplePre.h"
#include "SampleD.h"

#define PI 3.14159
NTL_CLIENT

void GenSamplePre(vec_ZZ& e, const mat_ZZ& A, const mat_ZZ& BS, const ZZ& s, const vec_ZZ& y, const RR& r)
{
	mat_ZZ B, A1,A0;
	vec_ZZ z, e1, x, e2;
	long m1 = A.NumCols();
	long n1 = A.NumRows();
	long m2 = BS.NumCols();
	long k = m1 / m2;
// cout << "test1" << endl;
	ident(B, (k - to_long(s))*m2);
	
	vec_RR c;
	c.SetLength(B.NumCols());
	SampleD(e1, B, r, c);
	
	A1.SetDims(n1, m1-m2);

	long i, j;
	for (i = 1; i <= n1; i++)
	{
		for (j = 1; j <= m1 - m2; j++)
		{

			A1(i, j) = A(i, j + m2);
		}
	}
// cout << "test2" << endl;
	//z = y - A1*e1;

	mul(x, A1, e1);

	
	sub(z, y, x);
	//ZZ q;
	//q = 10;
	//ZZ_p::init(q);
	A0.SetDims(n1, m2);
	for (i = 1; i <= n1; i++)
		for (j = 1; j <= m2; j++) {

			A0(i, j) = A(i, j);
		}
	SamplePre(e2, A0, BS, r, z);
	e.SetLength(m1);
	for (i = 1; i <= m2; i++)
	{
		//e(i) = to_ZZ_p(e2(i));
		e(i) = e2(i);
	}
// cout << "test3" << endl;

	for (i = m2 + 1; i <= m1; i++)
	{
		//e(i) = to_ZZ_p(e1(i - m2));
		e(i) = e1(i-m2);
	}
// cout << "test4" << endl;
}