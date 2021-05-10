#pragma once
#pragma once
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

// �������ɾ���A�Լ���Ӧ��Ķ̻�
//����n,q���������A0,�̻�B0
void TrapGen(mat_ZZ& A0, mat_ZZ& B0, long n, const ZZ& q)
{
	RR z;
	div(z, log(to_RR(q)), log(to_RR(2)));      //z = log(to_RR(q)) / log(to_RR(2)) = log2(q) 

	ZZ r, r1, L, b;     //
	RR l, r2, r3, a;   //
	mat_ZZ A1, A2, H, H1, I, G, G1, P, B, R, e, T, Y;
	mat_ZZ_p X;
	ZZ_p::init(q);     //q��ֵ��ģ�� p
	vec_ZZ  c;

	r = 2;
	div(r2, log(to_RR(r)), log(to_RR(2)));      // r2 = log(to_RR(r)) /log(to_RR(2)) = log2(r)
	div(r3, log(to_RR(q)), log(to_RR(2)));     // r3  = log(to_RR(q)) / log(to_RR(2)) =  log2(q)
	div(l, r3, r2);
	CeilToZZ(L, l);    //����ȡ���� L = l, truncated to +infinity 

	long m1 = to_long(2)* n*to_long(z);      //m1 =2 n log2(q)
	long m2 = m1 * to_long(L);    //m2 >m1 *L
	long m;
	long i, j, k;
	A1.SetDims(n, m1);    //


label:
	for (i = 1; i <= n; i++)
		for (j = 1; j <= m1; j++)
		{
			RandomBnd(A1(i, j), q);  //����A1, A1�����Ԫ���� 0~q-1 ֮��
		}
	image(X, to_mat_ZZ_p(A1));     //X��A1���н����ξ���Ҳ�������Ǿ��󣬼�A'
//	cout << "\nX= \n" << X << endl;
	long rX;
	rX = gauss(X);

	mat_ZZ_p M;
	ZZ_p d;

	M.SetDims(n, n);
	transpose(X, X);              //ת��X
	for (i = 1; i <= n; i++)
	{
		M(i) = X(i);           //ȡת�ú�X��ǰ n ��
	}
	transpose(M, M);     //�ٽ�Mת�û���
	long rM;
	rM = gauss(M);      //M����rM
	if (rM != n)            //���M���Ȳ�����n������label����������X
	{
		goto label;
	}
	//	cout << "\nX= \n" << X << endl;
	//	cout <<"\nM= \n" << M << endl;             

	mat_ZZ_p Q;
	mat_ZZ_p x4;        //
	x4.SetDims(m1 - n, n);
	Q.SetDims(1, n);              //����1 *n�������
//	cout << "Q = " << Q << endl;
	for (j = n + 1; j <= m1; j++)
	{
		//	negate(X(j), X(j));
		sub(X(j), Q(1), X(j));                //X(j) = - X(j)
		solve(d, M, x4(j - n), X(j));   //    ���d!=0; M*x4 = X(j) , d = determinant(M)   // �����Է�����Ax = b
	}
	transpose(x4, x4);
	//	cout << "x4= \n" << x4 << endl;

	mat_ZZ_p x5;
	ident(x5, m1 - n);        //����һ�� ��m1-n��*(m1-n)�ĵ�λ����
//	cout << "\nx5= \n" << x5 << endl;
	mat_ZZ_p x6;        //��x4��x5��ϳ�x6�� x5�ŵ�x4���·�
	x6.SetDims(m1, m1 - n);          //x6��һ�� m1 *(m1-n)�ľ���
	for (i = 1; i <= n; i++) {
		x6(i) = x4(i);
	}
	for (i = n + 1; i <= m1; i++)
	{
		x6(i) = x5(i - n);
	}
	//	cout <<  "\nx6 = \n" << x6 << endl;

	mat_ZZ x7;
	diag(x7, m1, q);        //����һ���Խ��󣬶Խ����ϵ�Ԫ����q
//	cout << "\nx7 = \n" << x7 << endl;
	transpose(x6, x6);            //ת��x6, (m1-n)* m1 
	for (i = 1; i <= m1 - n; i++)
	{
		x7(i) = to_vec_ZZ(x6(i));     //��ת�ú�x6�� (m1-n)�У�������x6, �ŵ� x7 ��ǰ (m1-n)��  ��Ҳ���ǽ�x7 ��ǰ (m1-n)���滻Ϊת�ú��x6
	}
	transpose(x7, x7);    //ת��x7������ԭ��������
//	cout << "\nx7=\n"<<x7 << endl;
//	cout << "A1=" << A1 << endl;                


	ZZ D;
	determinant(D, x7);
	if (D == 0)
	{
		goto label;       //���D = 0, ���� label ����ִ�У���ֹ����Hermite����ʱ D=0��
	}



	//���ɾ���A1��Hermite normal form matrix	
	transpose(x7, x7);
	HNF(H, x7, D);
	transpose(H, H);
	//	cout << H << endl;

	mat_ZZ AH;
	mul(AH, A1, H);
	//	cout << "AH= \n" << AH << endl;


	ident(I, m1);      //IΪ m1*m1 �׵ĵ�λ����
	sub(H1, H, I);    //H1 = H - I
	mat_ZZ_p H4;
	H4 = to_mat_ZZ_p(H1);   //��
	for (i = 1; i <= m1; i++)
	{
		H1(i) = to_vec_ZZ(H4(i));
	}
	//	cout << "H1= \n" << H1 << endl;

		//����G
	G.SetDims(m1, m2);      //����G�����ά����m1 *m2
	mat_ZZ W;
	W.SetDims(m1, m1);
	transpose(G, G);
	transpose(H1, H1);
	for (i = 1; i <= m1; i++)
	{
		G(i*to_long(l)) = H1(i);
	}
	//	cout << "\nl = \n" << to_long(l) << endl;

	for (i = 1; i <= m1; i++)
	{
		//	for (j = 1; j <= m1; j++)
		for (j = to_long(l) - 1; j >= 1; j--)
		{
			power(r1, r, to_long(l) - j); // x = a^e (e >= 0)
			div(a, to_RR(H1(i, j)), to_RR(r1));
			FloorToZZ(W(i, j), a);  // b = a, truncated to -infinity 
		}
		for (j = i * to_long(l) - 1; j >= (i - 1)*to_long(l) + 1; j--)
		{
			G(j) = W(i);
		}
	}
	transpose(G, G);

	//	cout << "\nG=\n" << G << endl;

		//����P
	P.SetDims(m2, m1);
	ident(e, m2);      //eΪm2�׵ĵ�λ����
	transpose(P, P);
	for (j = 1; j <= m1; j++) {

		P(j) = e(j*to_long(l));
	}
	transpose(P, P);

	//	cout << "\nP = \n" << P << endl;
	mat_ZZ GP;
	mul(GP, G, P);
	//	cout << "\nGP = \n" << GP << endl;


	   //����B
	B.SetDims(m2, m2);
	ident(B, m2);
	T.SetDims(to_long(l), to_long(l));     //����ά��Ϊ l *l �ľ���T
	ident(T, to_long(l));     //��Ϊl�ĵ�λ����T
	for (i = 1; i <= to_long(l); i++)
		for (j = i + 1; j <= to_long(l); j++)
		{
			T(i, j) = -r;
		}

	for (k = 1; k <= m1; k++)
	{
		for (i = 1; i <= to_long(l); i++)
			for (j = i; j <= to_long(l); j++) {

				B((k - 1)*to_long(l) + i, (k - 1)*to_long(l) + j) = T(i, j);
			}
	}

	//	cout << "\nB = \n" << B<< endl;


		//����R
	R.SetDims(m1, m2);
	c.SetLength(4);       //����һ������Ϊ4������
	c(0) = -1;
	c(1) = c(2) = 0;
	c(3) = 1;
	for (i = 1; i <= m1; i++)
		for (j = 1; j <= m2; j++)
		{
			R(i, j) = c(RandomBnd(4));  //   ���ɵ�Ԫ���� 0~3 ֮��
		}



	mat_ZZ X1, X2, X3, A, S, S1, S2;
	ZZ_p::init(q);
	mat_ZZ_p A22;
	add(X1, R, G);     //X1= R+G
	mul(X2, A1, X1);   //X2 = A1 *X1 =A1*(R+G)

	mat_ZZ Q1;
	Q1.SetDims(n, m2);     //����һ��n *m1�������
	sub(A2, Q1, X2);    //A2 = -X2 =  - A1*(R+G)
//	negate(A2, X2);     //A2 = -X2 =  - A1*(R+G)

	A22 = to_mat_ZZ_p(A2);
	m = m1 + m2;
	A.SetDims(n, m);
	transpose(A, A);    //ת��
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
	A0 = A;                      //�õ�����A0

	S.SetDims(m, m);
	mul(S1, X1, B);        //D = (R+G)*B = X1 *B 
	mul(X3, R, P);
	sub(S2, X3, I);         //V = RP -I
	for (i = 1; i <= m1; i++)
		for (j = 1; j <= m2; j++)
		{
			S(i, j) = S1(i, j);
		}

	for (i = 1; i <= m1; i++)
		for (j = 1; j <= m1; j++)
		{
			S(i, j + m2) = S2(i, j);
		}

	for (i = 1; i <= m2; i++)
		for (j = 1; j <= m2; j++)
		{
			S(i + m1, j) = B(i, j);
		}

	for (i = 1; i <= m2; i++)
		for (j = 1; j <= m1; j++)
		{
			S(i + m1, j + m2) = P(i, j);
		}

	B0 = S;   //�̻�



}