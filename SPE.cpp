#include<iostream>
#include <NTL/ZZ.h>        //��������������
#include <NTL/vec_ZZ.h>     //��������
#include <NTL/mat_poly_ZZ.h>     //����ʽϵ������
#include <NTL/mat_ZZ.h>      //��������
#include <NTL/mat_ZZ_p.h>
#include <time.h>      
#include<NTL/ZZ_p.h>
#include <NTL/LLL.h>
#include <NTL/matrix.h>  //����SetDims,�������ά��
#include <NTL/new.h>  
#include <NTL/lzz_p.h>    // ����ģ p
#include <NTL/lip.h>      //
#include <NTL/mat_lzz_p.h>    //����Ԫ��ģp
#include <NTL/vec_RR.h>  //ʵ������
#include <NTL/RR.h>   //ʵ����
#include<NTL/HNF.h>
#include <random>     //�����
#include <fstream>     
#include<math.h>
#include<vector>
#include<NTL/tools.h>
#include"TrapGen.h"
#include"SampleLeft.h"
#include <stdlib.h>
#include"errordistribution.h"
NTL_CLIENT


int main()
{
//	double time1, time2;
	RR s = to_RR(8);
	mat_ZZ A0, B0;                          //����A0�Ͷ̻�B0
	long n = 10;
	ZZ q;
	q = to_ZZ(127);
//	time1 = GetTime();

	//Setup
	//1.
	cout << "Setup begin." << endl;
	TrapGen(A0, B0, n, q);
	ZZ_p::init(q);

	// cout << "pk=\n" << endl;
	// cout << "A0 = \n" << A0 << endl;
//	cout << "B0 = \n" << B0 << endl;

	//2.
	long L = 5; 
	long n1 = A0.NumRows();        //A0���� = n1 = n
	long m1 = A0.NumCols();
//	cout << "n1 = \n" << n1 << endl;
//	cout << "m1 = \n" << m1 << endl;

	vector<Mat<ZZ_p> > matrixVector;       	//matrixVector
	for (long i = 0; i < L + 1; i++)
	{
		mat_ZZ_p D;
		D.SetDims(n1, m1);
		random(D, n1, m1);       //�������һ������
		matrixVector.push_back(D);      //��D���ӵ�matrixVector�ṹ�ĺ��棬������������Ϊ���������һ��Ԫ��
	//	cout << "D" << i << "=\n" << D<< endl;
		// cout << "D" << i+1 << "=\n" << matrixVector[i] << endl;

	}

//3.
	vec_ZZ u;
	u.SetLength(n1);
	for (long j = 1; j <= n1; j++)
	{
		RandomBnd(u(j), q);
	//	cout << "u=\n" << u(j) << endl;

	}
	// cout << "u=\n" << u << endl;

	//4.
	// cout << "msk = \n" << B0 << endl;


    // KeyGen
	// 1.
	cout << "SKGen begin:" << endl;
	mat_ZZ_p sumA, sumAB;   //����0����
	sumA.SetDims(n1, m1);
	sumAB.SetDims(n1, m1);
	string string1 = "10111";
    // cout << "���볤��Ϊ "<< L<< " ��01�ַ���:" << endl;
	// cin >> string1;
	for (long i = 0; i <L; i++) 
	{
		if (string1[i] == '1')
		{
			add(sumAB, matrixVector[L], matrixVector[i]);
			add(sumA, sumA , sumAB);
		}
	}
	// cout << "As=\n" << sumA << endl;
	

	//2.
	vec_ZZ e;
	mat_ZZ As;
	As.SetDims(n1, m1);
	for (long i = 1; i <=n1; i++)
	{
		As(i) =to_vec_ZZ( sumA(i));      //ת�����������
	}
//	cout << "As==\n" << As << endl;
	SampleLeft(e, A0, As, B0, s, u);
	//3.
	// cout << "sks=\n" << e << endl;

	/*
	���� ��A0|As)e = u
	vec_ZZ e1, e2;
	e1.SetLength(m1);
	e2.SetLength(m1);
	for (long i = 1; i <= m1; i++)
	{
		e1(i) = e(i);
	}
	for (long i = m1 + 1; i <= m1 + m1; i++)
	{
		e2(i - m1) = e(i);
	}
	
	vec_ZZ X1, X2,X3;
	vec_ZZ_p X4;
	X1.SetLength(n1);
	X2.SetLength(n1);
	mul(X1, A0, e1);
	mul(X2, As, e2);
	add(X3, X1, X2);
	X4 = to_vec_ZZ_p(X3);
	cout << "X3=\n" << X4 << endl;
	cout << "u=\n" << u << endl;
	*/


	//Encrypt
	cout << "Encrypt begin:" << endl;
	vec_ZZ v;
	ZZ M;
	string string2 = "10111";
	// cout << "���볤��Ϊ " << L << " ��01�ַ���:" << endl;
	// cin >> string2;
	v.SetLength(n1);
	for (long i = 1; i <= n1; i++) {
		RandomBnd(v(i), q);
	}
	RandomBnd(M, to_ZZ(2));     //Mȡ0����1

	
	mat_ZZ R;
	vec_ZZ F;
	vec_ZZ_p y;
	y.SetLength(m1);
	R.SetDims(m1, m1);
	F.SetLength(2);       //����һ������Ϊ4������
	F(0) = -1;
	F(1) =  1;
	for (long i = 1; i <= m1; i++)
	{
		errordistribution(y(i), q);                     //�õ� y

	}
	mat_ZZ_p Z1;
	Z1.SetDims(L, m1);
	for (long k = 1; k <= L; k++)
	{
		for (long i = 1; i <= m1; i++)
		{
			for (long j = 1; j <= m1; j++)
			{
				R(i, j) = F(RandomBnd(2));  //   ���ɵ�Ԫ���� 0~1֮��
			}
		}
		transpose(R, R);
		mul(Z1(k), to_mat_ZZ_p(R), y);   //�õ�Z(k) =Z(i)

	}

	vec_ZZ_p x1;      //xһ��Ԫ�ص�����
	x1.SetLength(1);
	errordistribution(x1(1), q);    //�õ�x


	//4.
	ZZ Q1,Q3,Q4,Q5;
	ZZ P;
	RR Q2;
	InnerProduct(Q1,u,v);
	div(Q2, to_RR(q), 2);
	FloorToZZ(Q3, Q2);        //
	mul(Q4, M, Q3);    //Q4
	vec_ZZ Q6;
	Q6.SetLength(1);
	Q6 = to_vec_ZZ(x1);
	add(Q5, Q1, Q6(1));
	add(P, Q5, Q4);

	ZZ p;
	p = rem(P, to_long(q));    //�õ���Χ�ڵ�p

	vec_ZZ  C0, E1;
	C0.SetLength(m1);
	E1.SetLength(n1);
	mat_ZZ E2;
	E2.SetDims(m1, n1);
	transpose(E2, A0);
	mul(E1, E2, v);
	add(C0, E1, to_vec_ZZ(y));
	vec_ZZ_p c0;
	c0 = to_vec_ZZ_p(C0);      //c0���Ƿ��Ϸ�Χ��c0


	//5.
	mat_ZZ_p G1, G2;
	mat_ZZ_p G3;
	G1.SetDims(n1, m1);
	G2.SetDims(m1, n1);
	G3.SetDims(L, m1);
	for (long i = 0; i < L; i++)
	{
		if (string2[i] == '1')
		{
			add(G1, matrixVector[L], matrixVector[i]);
			transpose(G2, G1);
			mul(G3(i+1), G2, to_vec_ZZ_p(v));
			add(G3(i+1), G3(i+1), Z1(i+1));
		}
	}
	// cout << "CT=\n" << endl;
	// cout << "p= " << p << endl;
	// cout << "co= " << c0 << endl;
	// cout << "G3=\n" << G3 << endl;


	//Decrypt
	//1.
	cout << "Decrypt  begin:" << endl;
	vec_ZZ_p G4;
	vec_ZZ_p G5;
	ZZ G6, W;

	G5.SetLength(2 * m1);
	G4.SetLength(m1);
	for (long i = 0; i < L; i++)
	{
		if (string1[i] == '1')
		{
			add(G4, G4, G3(i+1));             //�õ�һ���ۼӺ�
		}
	}
	// cout << "G4=\n" << G4 << endl;


	for (long i = 1; i <= m1; i++)
	{
		G5(i) = c0(i);
	}
	for (long i = m1 + 1; i <= m1 + m1; i++)
	{
		G5(i) = G4(i-m1);
	}

	// cout << "G5=\n" << G5 << endl;
	InnerProduct(G6, e, to_vec_ZZ(G5));
	sub(W, p, G6);

	ZZ w1;
	w1 = rem(W, to_long(q));      //w1 =w �ڷ�Χ֮��

	//2.
	ZZ G7, G8;
	sub(G7, w1, Q3);
	abs(G8, G7);         //ȡ����ֵ

	// cout << "G8=\n" << G8 << endl;
	// cout << "w1=\n" << w1 << endl;
	// cout << "M=\n" << M << endl;


	RR T1;
	ZZ T2;
	div(T1, to_RR(q), to_RR(4));
	FloorToZZ(T2, T1);             //q/4 

	ZZ M5;
	if (G8 < T2)
	{
		M5 = to_ZZ(1);
	}
	else
	{
		M5 = to_ZZ(0);
	}

	if (M5 == M)
	{
		cout << "Correct." << endl;
	}
	else
	{
		cout << " Error." << endl;

	}

	return 0;


}