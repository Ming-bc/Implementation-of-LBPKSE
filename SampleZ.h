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

// 输入：采样参数 s, 采样中心 c, 安全参数 n
//输出： 从离散高斯分布中得出采样值 x
void SampleZ(ZZ& x, const RR& s, const RR& c, long n)
{
	RR res, f, t, z, g, h, b;
	ZZ st, f1, y;
	vec_ZZ  a;
	long  d;
// cout << "testZ1" << endl;
	div(res, log(to_RR(n)), log(to_RR(2)));   //res = log2(n)
	SqrRoot(t, res);    //
	mul(z, s, t);
	sub(f, c, z);            //f = c - t(n) 
	CeilToZZ(f1, f);   // 取上整数
	add(b, c, z);           //c + t(n)  
	FloorToZZ(st, b);    //c + t(n)  取下整数
// cout << "st: " << st << " f1: " << f1 << endl;
// cout << "d: " << st - f1 + 1 << endl;
	d = to_long(st - f1 + 1);  //   元素个数
	a.SetLength(d);      //将向量 a 设置为一个长度为 d 的向量
	long i;
	for (i = 1; i <= d; i++)
	{
		a(i) = f1;
		f1++;
	}
// cout << "testZ4" << endl;

while(true){
	y = a(to_long(1)+RandomBnd(d));    //y就是采样值x
	div(g, to_RR(PI*sqr(abs(to_RR(y) - c))), sqr(s));
//	cout << "\ng=\n" << g << endl;
	sub(g, 0, g);         //g = -g
	exp(h, g);   
//	cout << "\nh=\n" << h << endl;
	if (RandomBnd(to_ZZ(inv(h))) == 0){    //以 一定 概率输出
		x = y;
		break;
	}
	else {
		continue;
	}
}
// cout << "testZ5" << endl;
}
