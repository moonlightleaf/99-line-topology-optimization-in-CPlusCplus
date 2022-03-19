//无关网格敏度过滤子程序
#pragma once
#include<opencv.hpp>
cv::Mat check(const int nelx, const int nely, const float rmin, const cv::Mat& x, const cv::Mat& dc);