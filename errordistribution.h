#pragma once
#pragma once
#pragma once
#include <iostream>  
#include <time.h>  
#include <iomanip>  
#include <math.h>  
#include <fstream>  
#include<NTL/ZZ_p.h>
#include <NTL/ZZ.h>
#define PI 3.14159  
NTL_CLIENT

using namespace std;
void UNIFORM(double *);  //UINFORM��������  

int x = 0;   //���ﶨ��xһ��ȫ�ֱ������ҳ�ʼ��ֵ0������Ĺ��ý������Ӻ���UNIFORM�е�������  
int round_double(double number)
{
	return (number > 0.0) ? floor(number + 0.5) : ceil(number - 0.5);
}

void errordistribution(ZZ_p& r2, const ZZ& q)
{
	//	int i, j;
		//ZZ r2;
	double a, A, B, C, E, D, r1;//
	double uni[2];
	double *p;
	srand((unsigned)time(NULL));  //��������Ӳ���ϵͳʱ��  
	a = 0.05;
	E = 0;
	D = (a / sqrt(2 * PI))*(a / sqrt(2 * PI));           //theta ��ƽ��
	ZZ_p::init(q);

	UNIFORM(&uni[0]);  //����UNIFORM��������2�����ȷֲ������������������nui[2]  
	A = sqrt((-2)*log(uni[0]));
	B = 2 * PI*uni[1];
	C = A * cos(B);
	r1 = E + C * D;    //E,D�ֱ��������ͷ���  
	if (isinf(r1) == 0)                            //inf һ��������С���½�
	{
		r2 = to_ZZ_p(round(to_double(q)*r1));
	}
	//	r2 = r1;
		//round
		//rem( r,  r2,q);
}
void UNIFORM(double *p)
{
	int i, a;
	double f;
	for (i = 0; i < 2; i++, x = x + 689)
	{
		a = rand() + x;  //����689����Ϊϵͳ����������ĸ���Ƶ��ԶԶ����������ú�����ʱ��  
		a = a % 1000;
		f = (double)a;
		f = f / 1000.0;
		*p = f;
		p++;
	}
}

