//采用优化准则法(OC)迭代子程序
#pragma once
#include<opencv.hpp>
cv::Mat OC(const int nelx, const int nely, const cv::Mat& x, const float volfrac, const cv::Mat& dc);
