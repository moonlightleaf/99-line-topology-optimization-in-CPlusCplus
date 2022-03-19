#include"check.h"
double max(const double a, const double b) {
	return a > b ? a : b;
}
double min(const double a, const double b) {
	return a < b ? a : b;
}
cv::Mat check(const int nelx, const int nely, const float rmin, const cv::Mat& x, const cv::Mat& dc) {
	//dcn用来保存更新的目标函数的灵敏度
	cv::Mat dcn = cv::Mat::zeros(cv::Size(nelx, nely), CV_32FC1);
	for (int i = 0; i < nelx; ++i) {
		for (int j = 0; j < nely; ++j) {
			float sum = 0.f;
			//在过滤半径定义的范围内遍历
			int rminfloor = rmin;
			int kLow;
			for (int k = max(i - rminfloor, 0); k <= min(i + rminfloor, nelx-1); ++k) {
				for (int l = max(j - rminfloor, 0); l <= min(j + rminfloor, nely - 1); ++l) {
					//fac即公式中的卷积算子Hf
					float fac = rmin - sqrt((i - k)*(i - k) + (j - l)*(j - l));
					sum += max(fac, 0);
					dcn.at<float>(j, i) = dcn.at<float>(j, i) + max(0, fac)*x.at<float>(l, k)*dc.at<float>(l, k);
				}
			}
			dcn.at<float>(j, i) = dcn.at<float>(j, i) / (x.at<float>(j, i)*sum);
		}
	}
	return dcn;
}