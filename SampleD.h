#pragma once
#ifndef SAMPLED_H_
#define SAMPLED_H_

#include "SampleZ.h"

#define PI 3.14159
NTL_CLIENT


//输入： n 维格基B, 采样参数 s， 中心点 c（是个 n 维的向量）
//输出：采样的格点 v
void SampleD(vec_ZZ& v, const mat_ZZ B,  const RR& s, const vec_RR& c)
{
	long n1 = B.NumRows();
	long m1 = B.NumCols();
// cout << n1 << " " << m1 << endl;
	mat_RR B1, C, B2;
	B1.SetDims(n1, m1);    //定义矩阵的维数
	B2.SetDims(n1, m1);
	C.SetDims(n1, m1);
	C(n1) = c;
//	cout << "C =\n" << C << endl;
	long i, j;
// cout << "test1" << endl;
	for (i = 1; i <= n1; i++) {
		for (j = 1; j <= m1; j++) {
			B1(i, j) = to_RR(B(i, j));
		}
	}
	B2 = B1;
	/*
	for (i = 1; i <= n1; i++)
	{
		for (j = 1; j <= m1; j++) {
			B2(i, j) = to_RR(B1(i, j));
		}
	}
	*/
// cout << "test2" << endl;
	mat_ZZ V, V1;
	V.SetDims(n1, m1);
	V1.SetDims(m1, n1);
	mat_RR mu;
	vec_RR c1, y1, g4;
	vec_ZZ z, g3;
	mu.SetDims(n1, m1);
	c1.SetLength(m1);     //定义向量的长度（范数）
	y1.SetLength(m1);
	z.SetLength(m1);
	g3.SetLength(m1);
	g4.SetLength(m1);
	RR g1, g2,  c2, s1, s2;
// cout << "test3" << endl;
	transpose(V1, B);
	ComputeGS(V1, mu, c1);    //求矩阵B 的正交化，mu 是正交化元素中减数向量的系数组成的矩阵， c1 是贝塔的平方组成的向量	
	transpose(B1, B1);
	for (i = 1; i <= n1; i++) {
		for (j = 1; j <= i - 1; j++) {
			mul(y1, mu(i, j), B1(j));      //y1 = mu(i, j) * B1(j), 正交化元素中的减数，参照施密特正交化公式
			sub(B1(i), B1(i), y1);       //循环之后，得到一个累加的结果
		}
	}
	//经过以上步骤，得到的 B1 是格基 B 经过施密特正交化后的矩阵
// cout << "test4" << endl;

	//第一步
	for (i = n1; i >= 1; i--)
	{
// cout << i << endl;
		InnerProduct(g1, C(i), B1(i));   //内积，得到一个数g1
		InnerProduct(g2, B1(i), B1(i));   //
		div(c2, g1, g2);     // c2 = c'
ZZ qq;
qq = to_ZZ(17); // 5,31
ZZ_p::init(qq);
c2 = to_RR(rem(to_ZZ(c2), to_long(qq)));
		div(s1, s, SqrRoot(g2));    //s1 =s'
// cout << "test41" << endl;
		SampleZ(z(i),  s1,  c2, n1);   //z(i)
// cout << "test42" << endl;
		mul(g3, z(i), V1(i));     //z(i) b(i）
		mul(g4, to_RR( z(i)), B2(i));  //
		sub(C(i - 1), C(i), g4);
		add(V(i - 1), V(i), g3);
	}
//	第二步：输出v0
	v = V(0);
// cout << "test5" << endl;
}
#endif