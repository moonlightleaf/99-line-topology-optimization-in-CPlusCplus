//求解单元刚度矩阵子程序
#pragma once
#include<opencv.hpp>
#include<vector>
cv::Mat lk(const float E = 1.0f, const float nu = 0.3f);