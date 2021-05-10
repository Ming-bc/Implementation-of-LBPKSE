#pragma once
#include<stdlib.h>
#include<stdio.h>
#include<iostream>

struct emxArray_cell_wrap_0
{
	cell_wrap_0 *data;
	int *size;
	int allocatedSize;
	int numDimensions;
	boolean_T canFreeData;
};

typedef struct {
	emxArray_real_T *f1;
} cell_wrap_0;

struct emxArray_real_T
{
	double *data;
	int *size;
	int allocatedSize;
	int numDimensions;
	boolean_T canFreeData;
};
