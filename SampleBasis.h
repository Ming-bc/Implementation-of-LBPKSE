#pragma once
#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>

#include <NTL/mat_ZZ.h>
#include <time.h>
#include <NTL/LLL.h>

#include <NTL/vec_RR.h>
#include <NTL/RR.h> 
#include <random>
#include "GenSamplePre.h"


#define PI 3.14159
NTL_CLIENT




void SampleBasis(mat_ZZ& B, const mat_ZZ& A, const mat_ZZ& BS, const ZZ& s, const RR& L)
{
	vec_ZZ y;
	RR z, r1;
	long m1 = A.NumCols();
	long n1 = A.NumRows();
	y.SetLength(n1);
// cout << "test" << endl;
	for (long i = 1; i <= m1; i++)
	{
// cout << i << endl;
		GenSamplePre(B(i), A, BS, s, y, L);
	}
// cout << "test" << endl;
	transpose(B,B);
}