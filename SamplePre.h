#pragma once
#pragma once
#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <time.h>
#include <NTL/LLL.h>
#include <NTL/vec_RR.h>
#include <NTL/RR.h> 
#include <random>
#include "SampleD.h"

#define PI 3.14159
NTL_CLIENT



//输出采样格点向量 v
void SamplePre(vec_ZZ& v, const mat_ZZ& A, const mat_ZZ& B, const RR& s, const vec_ZZ& u)
{
	long n = A.NumRows();
	long m = A.NumCols();
	vec_ZZ t, x;
	vec_RR t1, t2;
	mat_ZZ C;
	t.SetLength(m);
	t1.SetLength(m);
	t2.SetLength(m);
	transpose(C, A);
	long reduce = 0;

	LatticeSolve(t, C, u, reduce);
	if (IsZero(t))
		Error("LatticeSolve: bad reduce parameter");

	long i;
	for (i = 1; i <= m; i++)
	{
		t1(i) = to_RR(t(i));
	}
	sub(t1, t2, t1);
	//	negate(t1, t1);
	SampleD(x, B, s, t1);       //得到 x
	add(v, x, t);



}