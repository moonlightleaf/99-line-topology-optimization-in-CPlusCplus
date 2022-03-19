#include"top.h"
void top(const int nelx, const int nely, const float volfrac, const float penal, const float rmin) {
	//初始化设计变量x为体积分数
	cv::Mat x = cv::Mat::ones(cv::Size(nelx, nely), CV_32FC1);
	x = x * volfrac;
	//储存循环次数
	int loop = 0;
	//储存每次迭代后设计变量的变化量，用以判定是否收敛
	float change = 1.f;
	//用于保存上一次迭代后的设计变量
	cv::Mat xold(nely, nelx, CV_32FC1);
	//单元刚度矩阵
	const cv::Mat KE = lk().clone();
	/*****************************************************************************************************************************/
	//一些变量提前申请一块内存，防止内存不够用导致计算结果错误
	
	//设计变量，x和xold   还没用上
	cv::Mat *ptrx = new cv::Mat(cv::Mat::ones(cv::Size(nelx, nely), CV_32FC1));
	*ptrx = *ptrx * volfrac;
	cv::Mat *ptrxold = new cv::Mat(*ptrx);
	//全局节点位移矩阵
	cv::Mat *ptrU = new cv::Mat(cv::Mat::zeros(cv::Size(1, 2 * (nelx + 1)*(nely + 1)), CV_32FC1));
	/*****************************************************************************************************************************/

	//循环主体
	while (change > 0.01) {
		++loop;

		//先保存上一次的设计变量
		xold = x.clone();
		//每次迭代都进行一次有限元分析，计算结点位移，并储存在全局位移向量U中
		ptrU = FE(nelx, nely, x, penal, ptrU);
		//用于保存目标函数的变量
		float c = 0.f;
		//用于保存敏度
		cv::Mat dc(nely, nelx, CV_32FC1);
		dc = 0.f*dc;

		//遍历设计域矩阵元素，从左上角开始，一列一列
		for (int ely = 0; ely < nely; ++ely) {
			for (int elx = 0; elx < nelx; ++elx) {
				//计算左上角n1、右上角n2单元节点编号,单元的编号规则是从0开始纵向逐列递增，符合cpp编码习惯，而非matlab从1编码
				int n1 = (nely + 1)*elx + ely;
				int n2 = (nely + 1)*(elx + 1) + ely;
				//拿到当前单元的位移向量，edof是当前单元自由度在全局位移矩阵中的下标索引
				std::vector<int> edof = { 2 * n1, 2 * n1 + 1, 2 * n2,2 * n2 + 1, 2 * n2 + 2, 2 * n2 + 3, 2 * n1 + 2, 2 * n1 + 3 };
				cv::Mat Ue(8, 1, CV_32FC1);//用来存储单元位移向量
				for (int i = 0; i < 8; ++i) {
					Ue.at<float>(i, 0) = (*ptrU).at<float>(edof[i], 0);
				}
				
				//计算总体结构柔度值
				cv::Mat temp = Ue.t()*KE*Ue;
				c += pow(x.at<float>(ely, elx), penal)*temp.at<float>(0, 0);
				//计算每个敏度值填入dc对应位置
				dc.at<float>(ely, elx) = -1 * penal*pow(x.at<float>(ely, elx), penal - 1)*temp.at<float>(0, 0);
			}
		}
		//无关网格敏度过滤
		dc = check(nelx, nely, rmin, x, dc).clone();
		//采用OC求解设计变量
		x = OC(nelx, nely, x, volfrac, dc);
		//变化量最大的设计变量的change值；有正有负，分别找最大最小，取二者中绝对值最大的作为change值
		cv::Mat xchange = x - xold;
		double changeMin, changeMax;
		cv::minMaxLoc(xchange, &changeMin, &changeMax);
		change = abs(changeMax) > abs(changeMin) ? abs(changeMax) : abs(changeMin);

		//打印信息
		std::cout << "循环次数： " << loop << " 目标函数： " << c 
			<< " 体积分数： " << cv::sum(x)[0] / (nelx*nely) << " 设计变量变化： " << change << std::endl;
	}
	std::cout << x << std::endl;
}