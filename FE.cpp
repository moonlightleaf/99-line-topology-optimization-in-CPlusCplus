#include"FE.h"
cv::Mat* FE(const int nelx, const int nely, const cv::Mat& x, const float penal, cv::Mat* ptrU) {
	//计算单元刚度矩阵
	cv::Mat *ptrKE = new cv::Mat(lk());
	//总体刚度矩阵
	cv::Mat *ptrK = new cv::Mat(cv::Mat::zeros(cv::Size(2 * (nelx + 1)*(nely + 1), 2 * (nelx + 1)*(nely + 1)), CV_32FC1));
	//力矩阵
	cv::Mat *ptrF = new cv::Mat(cv::Mat::zeros(cv::Size(1, 2 * (nelx + 1)*(nely + 1)), CV_32FC1));
	//全局节点位移矩阵
	*ptrU = cv::Mat::zeros(cv::Size(1, 2 * (nelx + 1)*(nely + 1)), CV_32FC1);
	for (int elx = 0; elx < nelx; ++elx) {
		for (int ely = 0; ely < nely; ++ely) {
			//计算左上角n1、右上角n2单元节点编号,单元的编号规则是从0开始纵向逐列递增，符合cpp编码习惯，而非matlab从1编码
			int n1 = (nely + 1)*elx + ely;
			int n2 = (nely + 1)*(elx + 1) + ely;
			//组装总体刚度矩阵时需要插入的位置索引
			std::vector<int> edof = { 2 * n1, 2 * n1 + 1, 2 * n2,2 * n2 + 1, 2 * n2 + 2, 2 * n2 + 3, 2 * n1 + 2, 2 * n1 + 3 };
			for (int i = 0; i < 8; ++i) {
				for (int j = 0; j < 8; ++j) {
					(*ptrK).at<float>(edof[i], edof[j]) += pow(x.at<float>(ely, elx), penal)*(*ptrKE).at<float>(i, j);
				}
			}
		}
	}
	//施加载荷,本例应用了一个左上角的垂直单元力
	(*ptrF).at<float>(1, 0) = -1;
	//施加约束,左边第一列和右下角约束,已经固定的自由度在全局位移矩阵中的索引下标，约束左侧x和右下角的y方向自由度
	std::vector<int> fixeddofs;
	for (int i = 0; i < nely; ++i) {
		fixeddofs.push_back(i * 2);
	}
	fixeddofs.push_back(2 * (nelx + 1)*(nely + 1) - 1);
	//剩下的无约束自由度在全局位移矩阵中的下标索引
	std::list<int> prefreedofs;
	for (int i = 0; i < 2 * (nelx + 1)*(nely + 1); ++i) {
		prefreedofs.push_back(i);
	}
	for (const auto& i : fixeddofs) {
		prefreedofs.remove(i);
	}
	std::vector<int> freedofs(prefreedofs.begin(), prefreedofs.end());//存成vector方便接下来快速随机访问
	//求解线性方程组，得到各节点x、y方向的位移值储存在U中
	cv::Mat *ptrcalculatingK = new cv::Mat(cv::Mat::zeros(cv::Size(freedofs.size(), freedofs.size()), CV_32FC1));
	//cv::Mat calculatingK = cv::Mat::zeros(cv::Size(freedofs.size(), freedofs.size()), CV_32FC1);
	for (int i = 0; i < freedofs.size(); ++i) {
		for (int j = 0; j < freedofs.size(); ++j) {
			(*ptrcalculatingK).at<float>(i, j) = (*ptrK).at<float>(freedofs[i], freedofs[j]);
		}
	}
	cv::Mat *ptrcalculatingF = new cv::Mat(cv::Mat::zeros(cv::Size(1, freedofs.size()), CV_32FC1));
	//cv::Mat calculatingF = cv::Mat::zeros(cv::Size(1, freedofs.size()), CV_32FC1);
	for (int i = 0; i < freedofs.size(); ++i) {
		(*ptrcalculatingF).at<float>(i, 0) = (*ptrF).at<float>(freedofs[i], 0);
	}
	cv::Mat *ptrcalculatingU = new cv::Mat(cv::Mat::zeros(cv::Size(1, freedofs.size()), CV_32FC1));
	//cv::Mat calculatingU = cv::Mat::zeros(cv::Size(1, freedofs.size()), CV_32FC1);
	*ptrcalculatingU = (*ptrcalculatingK).inv(cv::DECOMP_SVD)*(*ptrcalculatingF);
	
	//得到用于计算的calculatingU后，将其按位置填回原来的全局U
	for (int i = 0; i < freedofs.size(); ++i) {
		(*ptrU).at<float>(freedofs[i], 0) = (*ptrcalculatingU).at<float>(i, 0);
	}
	delete ptrKE;
	delete ptrK;
	delete ptrF;
	
	delete ptrcalculatingF;
	delete ptrcalculatingK;
	delete ptrcalculatingU;
	/*std::cout << "U:" << std::endl;
	std::cout << (*ptrU).t() << std::endl;*/
	return ptrU;
}
