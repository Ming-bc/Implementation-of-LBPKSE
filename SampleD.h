#pragma once
#pragma once
#include "SampleZ.h"
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
#include <random>     //�����
#include <fstream>     
#include<stdlib.h>
#include <NTL/tools.h>
#define PI 3.14159
NTL_CLIENT


//���룺 n ά���B, �������� s�� ���ĵ� c���Ǹ� n ά��������
//����������ĸ�� v
void SampleD(vec_ZZ& v, const mat_ZZ B, const RR& s, const vec_RR& c)
{
	long n1 = B.NumRows();
	long m1 = B.NumCols();
	mat_RR B1, C, B2;
	B1.SetDims(n1, m1);    //��������ά��
	B2.SetDims(n1, m1);
	C.SetDims(n1, m1);
	C(n1) = c;
	//	cout << "C =\n" << C << endl;
	long i, j;

	for (i = 1; i <= n1; i++)
	{
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
	mat_ZZ V, V1;
	V.SetDims(n1, m1);
	V1.SetDims(n1, m1);
	mat_RR mu;
	vec_RR c1;
	vec_RR y1, g4;
	vec_ZZ z, g3;
	c1.SetLength(m1);     //���������ĳ��ȣ�������
	y1.SetLength(m1);
	z.SetLength(m1);
	g3.SetLength(m1);
	g4.SetLength(m1);
	RR g1, g2, c2, s1, s2;

	transpose(V1, B);
	ComputeGS(V1, mu, c1);    //�����B ����������mu ��������Ԫ���м���������ϵ����ɵľ��� c1 �Ǳ�����ƽ����ɵ�����
	transpose(B1, B1);

	for (i = 1; i <= n1; i++)
	{
		for (j = 1; j <= i - 1; j++)
		{
			mul(y1, mu(i, j), B1(j));      //y1 = mu(i, j) * B1(j), ������Ԫ���еļ���������ʩ������������ʽ
			sub(B1(i), B1(i), y1);       //ѭ��֮�󣬵õ�һ���ۼӵĽ��
		}
	}
	//�������ϲ��裬�õ��� B1 �Ǹ�� B ����ʩ������������ľ���


	//��һ��
	for (i = n1; i >= 1; i--)
	{
		InnerProduct(g1, C(i), B1(i));   //�ڻ����õ�һ����g1
		InnerProduct(g2, B1(i), B1(i));   //
		div(c2, g1, g2);     // c2 = c'
		div(s1, s, SqrRoot(g2));    //s1 =s'
		SampleZ(z(i), s1, c2, n1);   //z(i)

		mul(g3, z(i), V1(i));     //z(i) b(i��
		mul(g4, to_RR(z(i)), B2(i));  //
		sub(C(i - 1), C(i), g4);
		add(V(i - 1), V(i), g3);
	}
	//	�ڶ��������v0
	v = V(0);

}
