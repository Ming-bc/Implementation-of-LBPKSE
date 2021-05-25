#include<iostream>
#include <NTL/ZZ.h>        //��������������
#include <NTL/vec_ZZ.h>     //��������
#include <NTL/mat_poly_ZZ.h>     //����ʽϵ������
#include <NTL/mat_ZZ.h>      //��������
#include <NTL/mat_ZZ_p.h>
#include <time.h>      
#include <NTL/ZZ_p.h>
#include <NTL/LLL.h>
#include <NTL/matrix.h>  //����SetDims,�������ά��
#include <NTL/new.h>  
#include <NTL/lzz_p.h>    // ����ģ p
#include <NTL/lip.h>      //
#include <NTL/mat_lzz_p.h>    //����Ԫ��ģp
#include <NTL/vec_RR.h>  //ʵ������
#include <NTL/RR.h>   //ʵ����
#include <NTL/HNF.h>
#include <random>     //�����
#include <fstream>     
#include <math.h>
#include <vector>
  
#include <stdio.h>
#include <fstream>
#include <string>
#include <sstream>
#include <bitset>
#include <sys/time.h>


#include<NTL/tools.h>
#include"TrapGen.h"
#include"SampleLeft.h"
#include"SampleBasis.h"

#include <stdlib.h>
#include"errordistribution.h"

NTL_CLIENT






int main()
{
int loop = 10;
for(int i = 0; i < loop; i ++){

	RR s = to_RR(8); // 8
	mat_ZZ A0, B0;                          //����A0�Ͷ̻�B0
	long n = 5; // 5
	ZZ q;
	q = to_ZZ(17); // 17
	timeval time1, time2, time3, time4, time5, time6, time7, time8, time9, time10, time11, time12;
	timeval timeInit1, timeInit2;
	//Setup
// cout << "Setup begin." << endl;
	//1.
	gettimeofday(&time1, NULL);
	TrapGen(A0, B0, n, q);
	ZZ_p::init(q);


	// cout << "pk=\n" << endl;
	// cout << "A0 = \n" << A0 << endl;
	// cout << "B0 = \n" << B0 << endl;

	// 2.
	long L = 15;
	string string1 = "101111111111111";
	string string2 = "101101111111111";
	long n1 = A0.NumRows();        //A0���� = n1 = n
	long m1 = A0.NumCols();
	// cout << "L = \n" << L << endl;
	// cout << "n1 = " << n1 << endl << "m1 = " << m1 << endl;

	vector<Mat<ZZ_p> > matrixVector;       	//matrixVector
	for (long i = 0; i < L + 1; i++)
	{
		mat_ZZ_p D;
		D.SetDims(n1, m1);
		random(D, n1, m1);       //�������һ������
		matrixVector.push_back(D);      //��D���ӵ�matrixVector�ṹ�ĺ��棬������������Ϊ���������һ��Ԫ��
		// cout << "D" << i << "=\n" << D<< endl;
		// cout << "D" << i+1 << "=\n" << matrixVector[i] << endl;

	}

// 3.
	vec_ZZ u;
	u.SetLength(n1);
	for (long j = 1; j <= n1; j++)
	{
		RandomBnd(u(j), q);
		// cout << "u=\n" << u(j) << endl;

	}
	// cout << "u=\n" << u << endl;

	//4.
	// cout << "msk = \n" << B0 << endl;
	gettimeofday(&time2, NULL);

    // KeyGen
	// 1.
// cout << "SKGen begin." << endl;
	gettimeofday(&time3, NULL);
	mat_ZZ_p sumA, sumAB;   //����0����
	sumA.SetDims(n1, m1);
	sumAB.SetDims(n1, m1);
	// string string1;
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
	mat_ZZ As,B,AAs;
	As.SetDims(n1, m1);
	AAs.SetDims(n1,2*m1);
	B.SetDims(2*m1,2*m1);
	for (long i = 1; i <=n1; i++)
	{
		As(i) =to_vec_ZZ( sumA(i));      //ת�����������
	}
	// cout << "As==\n" << As << endl;
	SampleLeft(e, A0, As, B0, s, u);

	for (long i = 1; i <= n1; i++)
	{
		for (long j = 1; j <= m1; j++)
		{
			AAs(i, j) = A0(i, j);
		}
	}
	for (long i = 1; i <= n1; i++)
	{
		for (long j = m1+1; j <=2* m1; j++)
		{
			AAs(i, j) = As(i, j-m1);
		}
	}
	SampleBasis(B,AAs, B0, to_ZZ(1), s);
	

	//3.
	// cout << "sks=\n" << e << endl;
	// cout << "e=\n" << e << endl;
	// cout << "B=\n" << B << endl;
	gettimeofday(&time4, NULL);

	//Encrypt
	gettimeofday(&time5, NULL);
	vec_ZZ v;
	ZZ M;
	// string string2;
// cout << "Encrypt  begin." << endl;
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
	x1.SetLength(2);
	errordistribution(x1(1), q);
	errordistribution(x1(2), q);    //�õ�x


	//4.
	ZZ Q1,Q11,Q3,Q4,Q5;
	ZZ P1,P2;
	RR Q2;
	vec_ZZ kw;//keyword
	kw.SetLength(n1);
	for (long i = 1; i <= n1; i++) {
		RandomBnd(kw(i), q);
	}

	InnerProduct(Q1,u,v);
	InnerProduct(Q11, kw, v);
	div(Q2, to_RR(q), 2);
	FloorToZZ(Q3, Q2);        //
	mul(Q4, M, Q3);    //Q4
	vec_ZZ Q6;
	Q6.SetLength(2);
	Q6 = to_vec_ZZ(x1);
	add(Q5, Q1, Q6(1));
	add(P1, Q5, Q4);
	add(P2, Q11, Q6(2));
	// cout << "P1= " << P1 << endl;
	// cout << "P2= " << P2 << endl;
	ZZ p1,p2;
	p1 = rem(P1, to_long(q));    //�õ���Χ�ڵ�p
	p2 = rem(P2, to_long(q));

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
	// cout << "p= " << p1 << endl;
	// cout << "p1= " << p2 << endl;
	// cout << "co= " << c0 << endl;
	// cout << "G3=\n" << G3 << endl;
	gettimeofday(&time6, NULL);

	//Trapdoor
// cout << "Trapdoor begin." << endl;
	gettimeofday(&time7, NULL);
	vec_ZZ te;
	vec_ZZ kw1;//keyword
	te.SetLength(2*m1);
	kw1.SetLength(n1);
	for (long i = 1; i <= n1; i++) {
		RandomBnd(kw1(i), q);
	}
	SamplePre(te, AAs, B, s, kw1);
	// cout << "te= " << te << endl;
	gettimeofday(&time8, NULL);

	//Search
// cout << "Search begin." << endl;
	//gettimeofday(&time9, NULL);
	vec_ZZ_p G4, G5;
	ZZ G62, W2;
gettimeofday(&timeInit1, NULL);
gettimeofday(&time9, NULL);
	G5.SetLength(2 * m1);
	G4.SetLength(m1);
gettimeofday(&timeInit2, NULL);
	for (long i = 0; i < L; i++)
	{
		if (string1[i] == '1')
		{
			add(G4, G4, G3(i + 1));             //�õ�һ���ۼӺ�
		}
	}
	// cout << "G4=\n" << G4 << endl;

	for (long i = 1; i <= m1; i++)
	{
		G5(i) = c0(i);
	}
	for (long i = m1 + 1; i <= m1 + m1; i++)
	{
		G5(i) = G4(i - m1);
	}

	// cout << "G5=\n" << G5 << endl;
// timeval timeTemp;
// gettimeofday(&timeTemp, NULL);
// cout << timeTemp.tv_usec << endl;
	InnerProduct(G62, te, to_vec_ZZ(G5));
	sub(W2, p2, G62);
// gettimeofday(&timeTemp, NULL);
// cout << timeTemp.tv_usec << endl;
gettimeofday(&time10, NULL);

	ZZ w2, T2;
	RR T1;
	w2 = rem(W2, to_long(q));
	div(T1, to_RR(q), to_RR(4));
	FloorToZZ(T2, T1);

	ZZ M6;
	if (w2 < T2)
	{
		M6 = to_ZZ(1);
	}
	else
	{
		M6 = to_ZZ(0);
	}

	// cout << "M6=\n" << M6 << endl;
	// gettimeofday(&time10, NULL);

	//Decrypt
// cout << "Decrypt begin." << endl;
	//1.
	vec_ZZ_p L4, L5;
	// gettimeofday(&time11, NULL);
	ZZ L6, W1;

//gettimeofday(&time11, NULL);
	L5.SetLength(2 * m1);
	L4.SetLength(m1);
 gettimeofday(&time11, NULL);
	for (long i = 0; i < L; i++)
	{
		if (string1[i] == '1')
		{
			add(L4, L4, G3(i+1));             //�õ�һ���ۼӺ�
		}
	}
	// cout << "L4=\n" << L4 << endl;

	for (long i = 1; i <= m1; i++)
	{
		L5(i) = c0(i);
	}
	for (long i = m1 + 1; i <= m1 + m1; i++)
	{
		L5(i) = L4(i-m1);
	}

	// cout << "G5=\n" << G5 << endl;
// gettimeofday(&timeTemp, NULL);
// cout << timeTemp.tv_usec << endl;
	InnerProduct(L6, e, to_vec_ZZ(L5));
	sub(W1, p1, L6);
// gettimeofday(&timeTemp, NULL);
// cout << timeTemp.tv_usec << endl;
gettimeofday(&time12, NULL);
	ZZ w1;
	w1 = rem(W1, to_long(q));      //w1 =w �ڷ�Χ֮��

	//2.
	ZZ G7, G8;
	sub(G7, w1, Q3);
	abs(G8, G7);         //ȡ����ֵ

	// cout << "G8=\n" << G8 << endl;
	// cout << "w1=\n" << w1 << endl;
	// cout << "M=\n" << M << endl;
	// gettimeofday(&time12, NULL);

	//RR T1;
	//ZZ T2;
	//div(T1, to_RR(q), to_RR(4));
	//FloorToZZ(T2, T1);          //q/4 
   
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
		// cout << "Decrypt correct." << endl;
	}
	else
	{
		// cout << "Decrypt error." << endl;
	}

	cout << endl << "Excution time-" << i << endl;
	cout << "Setup: " << (1000000 * (time2.tv_sec - time1.tv_sec) + (time2.tv_usec - time1.tv_usec)) << " us" << endl;
    cout << "KeyGen: " << (1000000 * (time4.tv_sec - time3.tv_sec) + (time4.tv_usec - time3.tv_usec)) << " us" << endl;
    cout << "Encrypt: " << (1000000 * (time6.tv_sec - time5.tv_sec) + (time6.tv_usec - time5.tv_usec)) << " us" << endl;
	cout << "Trapdoor: " << (1000000 * (time8.tv_sec - time7.tv_sec) + (time8.tv_usec - time7.tv_usec)) << " us" << endl;
	cout << "Search: " << (1000000 * (time10.tv_sec - time9.tv_sec) + (time10.tv_usec - time9.tv_usec)) << " us" << endl;
	cout << "Decrypt: " << (1000000 * (time12.tv_sec - time11.tv_sec) + (time12.tv_usec - time11.tv_usec))
							+ (1000000 * (timeInit2.tv_sec - timeInit1.tv_sec) + (timeInit2.tv_usec - timeInit1.tv_usec)) << " us" << endl;
	//return 0;
}
return 0;
}