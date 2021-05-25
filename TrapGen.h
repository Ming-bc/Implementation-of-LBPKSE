#pragma once
#ifndef TRAPGEN_H_
#define TRAPGEN_H_

#include <NTL/ZZ.h>        //��������������
#include <NTL/vec_ZZ.h>     //��������
#include <NTL/mat_poly_ZZ.h>     //����ʽϵ������
#include <NTL/mat_ZZ.h>      //��������
#include <time.h>      
#include <NTL/LLL.h>
#include <NTL/matrix.h>  //����SetDims,�������ά��
#include <NTL/new.h>  
#include <NTL/lzz_p.h>    // ����ģ p
#include <NTL/lip.h>      //
#include <NTL/mat_lzz_p.h>    //����Ԫ��ģp
#include <NTL/vec_RR.h>  //ʵ������
#include <NTL/RR.h>   //ʵ����
#include <NTL/mat_ZZ_p.h>
#include<NTL/HNF.h>
//#include"Negate.h"
#include <random>     //�����
#include <fstream>     
#include <NTL/tools.h>
#define PI 3.14159
NTL_CLIENT


// 用于生成矩阵A以及相应格的短基
//输入n,q，输出矩阵A0,短基B0
void TrapGen(mat_ZZ& A0, mat_ZZ& B0, long n, const ZZ& q)
{
	RR z;
	div(z, log(to_RR(q)), log(to_RR(2)));      //z = log(to_RR(q)) / log(to_RR(2)) = log2(q) 

	ZZ r, r1, L, b;     //
	RR l, r2, r3, a;   //
	mat_ZZ A1, A2, H, H1, I, G, G1, P, B, R, e, T,Y;
	mat_ZZ_p X;
	ZZ_p::init(q);     //q赋值给模数 p
	vec_ZZ  c;  

	r = 2;
	div(r2, log(to_RR(r)), log(to_RR(2)));      // r2 = log(to_RR(r)) /log(to_RR(2)) = log2(r)
	div(r3, log(to_RR(q)), log(to_RR(2)));     // r3  = log(to_RR(q)) / log(to_RR(2)) =  log2(q)
	div(l, r3, r2);
	CeilToZZ(L, l);    //向上取整， L = l, truncated to +infinity 

	long m1 = to_long(2)* n*to_long(z);      //m1 =2 n log2(q)
	// long m2 = to_long(4)* n*to_long(z);    //m2 >m1 *L
	long m2 = m1 * to_long(L);    //m2 >m1 *L
	long m;
	long i, j, k;
	A1.SetDims(n, m1);    //
// cout << "test1" << endl;
	
	label:
	for (i = 1; i <= n; i++)
		for (j = 1; j <= m1; j++) 
		{
			RandomBnd(A1(i, j), q);  //生成A1, A1里面的元素在 0~q-1 之间
		}
	image(X, to_mat_ZZ_p(A1));     //X是A1的行阶梯形矩阵，也是上三角矩阵，即A'
//	cout << "\nX= \n" << X << endl;
	long rX;
	rX = gauss(X);

	mat_ZZ_p M;
	ZZ_p d;
	
	M.SetDims(n, n);
	transpose(X, X);              //转置X
	for (i = 1; i <= n; i++)
	{
		M(i) = X(i);           //取转置后X的前 n 行
	}
	transpose(M, M);     //再将M转置回来
	long rM;
	rM = gauss(M);      //M的秩rM
	if (rM != n)            //如果M的秩不等于n，返回label，重新生成X
	{
		goto label;
	}
//	cout << "\nX= \n" << X << endl;
//	cout <<"\nM= \n" << M << endl;             

	mat_ZZ_p Q;
	mat_ZZ_p x4;        //
	x4.SetDims(m1 - n, n);
	Q.SetDims(1, n);              //生成1 *n的零矩阵
//	cout << "Q = " << Q << endl;
	for (j = n + 1; j <= m1; j++)
	{
		//	negate(X(j), X(j));
		sub(X(j), Q(1), X(j));                //X(j) = - X(j)
		solve(d, M, x4(j - n), X(j));   //    如果d!=0; M*x4 = X(j) , d = determinant(M)   // 解线性方程组Ax = b
	}
	transpose(x4, x4);  
//	cout << "x4= \n" << x4 << endl;

	mat_ZZ_p x5;
	ident(x5, m1 - n);        //生成一个 （m1-n）*(m1-n)的单位矩阵
//	cout << "\nx5= \n" << x5 << endl;
	mat_ZZ_p x6;        //把x4和x5组合成x6， x5放到x4的下方
	x6.SetDims(m1, m1 - n);          //x6是一个 m1 *(m1-n)的矩阵
	for (i = 1; i <= n; i++) {
		x6(i) = x4(i);
	}
	for (i = n + 1; i <= m1; i++)
	{
		x6(i) = x5(i - n);
	}
//	cout <<  "\nx6 = \n" << x6 << endl;

	mat_ZZ x7;
	diag(x7, m1, q);        //生成一个对角阵，对角线上的元素是q
//	cout << "\nx7 = \n" << x7 << endl;
	transpose(x6, x6);            //转置x6, (m1-n)* m1 
	for (i = 1; i <= m1 - n; i++)
	{
		x7(i) = to_vec_ZZ(x6(i));     //将转置后x6的 (m1-n)行，即整个x6, 放到 x7 的前 (m1-n)行  ，也就是将x7 的前 (m1-n)行替换为转置后的x6
	}
	transpose(x7, x7);    //转置x7，还回原来的样子
//	cout << "\nx7=\n"<<x7 << endl;
//	cout << "A1=" << A1 << endl;                


	ZZ D;
	determinant(D, x7);
	if (D == 0)
	{
		goto label;       //如果D = 0, 跳回 label 重新执行（防止计算Hermite矩阵时 D=0）
	}


// cout << "test2" << endl;

	//生成矩阵A1的Hermite normal form matrix	
	transpose(x7, x7);
	HNF(H, x7, D);
	transpose(H, H);
//	cout << H << endl;

	mat_ZZ AH;
	mul(AH, A1, H);
//	cout << "AH= \n" << AH << endl;


	ident(I, m1);      //I为 m1*m1 阶的单位矩阵
	sub(H1, H, I);    //H1 = H - I
	mat_ZZ_p H4;
	H4 = to_mat_ZZ_p(H1);   //相
	for (i = 1; i <= m1; i++)
	{
		H1(i) = to_vec_ZZ(H4(i));
	}
//	cout << "H1= \n" << H1 << endl;
// cout << "test3" << endl;
	//定义G
	G.SetDims(m1, m2);      //定义G矩阵的维数是m1 *m2
	mat_ZZ W;
	W.SetDims(m1, m1);
	transpose(G, G);
// cout << m1 << " " << m2 << endl;
	transpose(H1, H1);
	for (i = 1; i <= m1; i++)
	{
// cout << i << endl;
		G(i*to_long(l)) = H1(i);
	}
//	cout << "\nl = \n" << to_long(l) << endl;
// cout << "test31" << endl;
	for (i = 1; i <= m1; i++)
	{
	//	for (j = 1; j <= m1; j++)
		for (j = to_long(l) - 1; j >= 1; j--)
		{
			power(r1, r, to_long(l) - j); // x = a^e (e >= 0)
			div(a, to_RR(H1(i, j)), to_RR(r1));
			FloorToZZ(W(i,j), a);  // b = a, truncated to -infinity 
		}
		for (j = i * to_long(l) - 1; j >= (i - 1)*to_long(l) + 1; j--)
		{
			G(j) = W(i);
		}
	}
	transpose(G, G);

//	cout << "\nG=\n" << G << endl;
// cout << "test4" << endl;
	//定义P
	P.SetDims(m2, m1);
	ident(e, m2);      //e为m2阶的单位矩阵
	transpose(P, P);
	for (j = 1; j <= m1; j++) {
		P(j) = e(j*to_long(l));
	}
	transpose(P, P);

//	cout << "\nP = \n" << P << endl;
	mat_ZZ GP;
	mul(GP, G, P);
//	cout << "\nGP = \n" << GP << endl;


   //定义B
	B.SetDims(m2, m2);
	ident(B, m2);
	T.SetDims(to_long(l), to_long(l));     //生成维度为 l *l 的矩阵T
	ident(T, to_long(l));     //阶为l的单位矩阵T
	for (i = 1; i <= to_long(l); i++) {
		for (j = i + 1; j <= to_long(l); j++) {
			T(i, j) = - r;
		}
	}

	for (k = 1; k <= m1; k++) {
		for (i = 1; i <= to_long(l); i++) {
			for (j = i; j <= to_long(l); j++) {
				B((k - 1)*to_long(l) + i, (k - 1)*to_long(l) + j) = T(i, j);
			}
		}
	}

//	cout << "\nB = \n" << B<< endl;
// cout << "test5" << endl;

	//定义R
	R.SetDims(m1, m2);
	c.SetLength(4);       //生成一个长度为4的向量
	c(0) = -1;
	c(1) = c(2) = 0;
	c(3) = 1;
	for (i = 1; i <= m1; i++) {
		for (j = 1; j <= m2; j++) {
			R(i, j) = c(RandomBnd(4));  //   生成的元素在 0~3 之间
		}
	}



	mat_ZZ X1, X2, X3, A, S, S1, S2;
	ZZ_p::init(q);
	mat_ZZ_p A22;
	add(X1, R, G);     //X1= R+G
	mul(X2, A1, X1);   //X2 = A1 *X1 =A1*(R+G)

	mat_ZZ Q1;
	Q1.SetDims(n, m2);     //生成一个n *m1的零矩阵
	sub(A2, Q1, X2);    //A2 = -X2 =  - A1*(R+G)
//	negate(A2, X2);     //A2 = -X2 =  - A1*(R+G)

	A22 = to_mat_ZZ_p(A2); 
	m = m1 + m2;
	A.SetDims(n, m);
	transpose(A, A);    //转置
	transpose(A1, A1);
	transpose(A22, A22);
	for (j = 1; j <= m1; j++)
	{
		A(j) = A1(j);
	}
	for (j = m1 + 1; j <= m; j++) 
	{
		A(j) = to_vec_ZZ(A22(j - m1));
	}
	transpose(A, A);
	A0 = A;

	S.SetDims(m, m);
	mul(S1, X1, B);        //D = (R+G)*B = X1 *B 
	mul(X3, R, P);
	sub(S2, X3, I);         //V = RP -I

	for (i = 1; i <= m1; i++) {
		for (j = 1; j <= m2; j++) {
			S(i, j) = S1(i, j);     
		}
	}

	for (i = 1; i <= m1; i++) {
		for (j = 1; j <= m1; j++) {
			S(i, j + m2) = S2(i, j);
		}
	}

	for (i = 1; i <= m2; i++) {
		for (j = 1; j <= m2; j++) {
			S(i + m1, j) = B(i, j);
		}
	}

	for (i = 1; i <= m2; i++) {
		for (j = 1; j <= m1; j++) {
			S(i + m1, j + m2) = P(i, j);
		}
	}
	
	B0 = S;   //短基
}

#endif