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
#include <random>     //�����
#include <fstream>     
#include<stdlib.h>
#include <NTL/tools.h>
#define PI 3.14159
NTL_CLIENT

// ���룺�������� s, �������� c, ��ȫ���� n
//����� ����ɢ��˹�ֲ��еó�����ֵ x
void SampleZ(ZZ& x, const RR& s, const RR& c, long n)
{
	RR res, f, t, z, g, h, b;
	ZZ st, f1, y;
	vec_ZZ  a;
	long  d;
	div(res, log(to_RR(n)), log(to_RR(2)));   //res = log2(n)
	SqrRoot(t, res);    //
//	cout << "\nt=\n" << t << endl;
	mul(z, s, t);

	sub(f, c, z);            //f = c - t(n) 
	CeilToZZ(f1, f);   // ȡ������
	add(b, c, z);           //c + t(n)  
	FloorToZZ(st, b);    //c + t(n)  ȡ������

	d = to_long(st - f1 + 1);  //   Ԫ�ظ���
	a.SetLength(d);      //������ a ����Ϊһ������Ϊ d ������

	long i;
	for (i = 1; i <= d; i++)
	{
		a(i) = f1;
		f1++;
	}


loop:
	y = a(to_long(1) + RandomBnd(d));    //y���ǲ���ֵx
	div(g, to_RR(PI*sqr(abs(to_RR(y) - c))), sqr(s));
	//	cout << "\ng=\n" << g << endl;
	sub(g, 0, g);         //g = -g
	exp(h, g);
	//	cout << "\nh=\n" << h << endl;
	if (RandomBnd(to_ZZ(inv(h))) == 0)    //�� һ�� �������
		x = y;
	else
	{
		goto loop;
	}
}