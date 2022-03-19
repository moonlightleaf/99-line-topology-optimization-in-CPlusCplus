//求解单元刚度矩阵子程序
#include"lk.h"
cv::Mat lk(const float E, const float nu) {
	const std::vector<float> k = { 1.0f / 2 - nu / 6, 1.0f / 8 + nu / 8, -1.0f / 4 - nu / 12, -1.0f / 8 + 3.0f * nu / 8,
		-1.0f / 4 + nu / 12, -1.0f / 8 - nu / 8, nu / 6, 1.0f / 8 - 3.0f * nu / 8 };
	float arrKE[8][8] = {
		k[0],k[1],k[2],k[3],k[4],k[5],k[6],k[7],
		k[1],k[0],k[7],k[6],k[5],k[4],k[3],k[2],
		k[2],k[7],k[0],k[5],k[6],k[3],k[4],k[1],
		k[3],k[6],k[5],k[0],k[7],k[2],k[1],k[4],
		k[4],k[5],k[6],k[7],k[0],k[1],k[2],k[3],
		k[5],k[4],k[3],k[2],k[1],k[0],k[7],k[6],
		k[6],k[3],k[4],k[1],k[2],k[7],k[0],k[5],
		k[7],k[2],k[1],k[4],k[3],k[6],k[5],k[0]
	};
	
	cv::Mat KE(cv::Size(8, 8), CV_32FC1);
	for (int i = 0; i < 8; ++i) {
		for (int j = 0; j < 8; ++j) {
			KE.at<float>(i, j) = E / (1 - nu * nu) * arrKE[i][j];
		}
	}
	return KE;
}