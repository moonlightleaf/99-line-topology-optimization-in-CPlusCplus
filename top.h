//拓扑优化迭代求解的函数主体
#pragma once
#include<iostream>
#include<opencv.hpp>
#include"OC.h"
#include"check.h"
#include"FE.h"
void top(const int nelx, const int nely, const float volfrac, const float penal, const float rmin);
