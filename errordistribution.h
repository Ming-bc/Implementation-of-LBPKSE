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
void UNIFORM(double *);  //UINFORM函数声明  

int x = 0;   //这里定义x一个全局变量并且初始付值0，这个的功用将会在子函数UNIFORM中得以体现  
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
	srand((unsigned)time(NULL));  //随机数种子采用系统时钟  
	a = 0.05;
	E = 0;
	D = (a / sqrt(2 * PI))*(a / sqrt(2 * PI));           //theta 的平方
	ZZ_p::init(q);

	UNIFORM(&uni[0]);  //调用UNIFORM函数产生2个均匀分布的随机数并存入数组nui[2]  
	A = sqrt((-2)*log(uni[0]));
	B = 2 * PI*uni[1];
	C = A * cos(B);
	r1 = E + C * D;    //E,D分别是期望和方差  
	if (isinf(r1) == 0)                            //inf 一个集合最小的下界
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
		a = rand() + x;  //加上689是因为系统产生随机数的更换频率远远不及程序调用函数的时间  
		a = a % 1000;
		f = (double)a;
		f = f / 1000.0;
		*p = f;
		p++;
	}
}

