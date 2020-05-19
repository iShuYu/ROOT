#include <iostream>
#include <TCanvas.h>
#include <TRandom.h>
#include <TMinuit.h>

using namespace std;
using namespace TMath;
const Int_t nentries = 100;
const Int_t N = 5;	//数据个数
const Int_t npar = 2;	//拟合参数的个数

Double_t x[N] = {1000., 828., 800., 600., 300.};
Double_t y[N] = {1500., 1340., 1328., 1172., 800.};
Double_t sigma = 15.;

void fcn(Int_t &npar, Double_t *deriv, Double_t &chi2, Double_t *par, Int_t iflag)
//这里npar是参数个数，*par是指针数组，这是固定格式哈，虽然可能有的没有用或者暂时不知道怎么用，比如deriv和iflag
{
	Double_t logL = 0;
	Double_t alpha, beta;
	alpha = par[0];
	beta = par[1];
	Double_t logf;
	for (Int_t i=0; i<N; i++)
	{
		//这题是高斯分布，所以logf是这种形式，指数分布类推，然后三种情况分别对应了这里的一二三题。
		logf = pow((y[i]-pow(x[i],beta)*alpha),2)/sigma/sigma;
		//logf = pow((y[i]-(alpha*x[i]+beta*pow(x[i],2))),2)/sigma/sigma;
		//logf = pow((y[i]-alpha*x[i]),2)/sigma/sigma;
		logL = logL + logf;
	}
	chi2 = logL;
	//无论如何，按照这个格式，把logL的和传输到第三个参数的位置即可
}
 
void leastsquare()
{
	//拟合部分
	TMinuit	*minuit = new TMinuit(npar);	//这里用指针或者不用都行
	minuit -> SetFCN(fcn);	//定义目标函数
	Double_t arglist[10];
	Int_t ierflg = 0;
	arglist[0] = 1;
	minuit -> mnexcm("SET ERR", arglist, 1, ierflg);
	//设置初始值和步长(可以根据数据自己估计个起始值)
	//这里的mnparm函数第一个是参数序数，第二个是起始值，第三个是步长，后面几个不是很懂是什么
	static Double_t vstart[npar] = {1., 1.};
	static Double_t step[npar] = {0.1, 0.1};
	minuit -> mnparm(0, "alpha", vstart[0], step[0], 0, 0, ierflg);
	minuit -> mnparm(1, "beta", vstart[1], step[1], 0, 0, ierflg);
	//现在进行最小化
	arglist[0] = 500;
	arglist[1] = 1.;
//	minuit -> mnexcm("SIMPLEX", arglist, 0, ierflg);
	minuit -> mnexcm("MIGRAD", arglist, 0, ierflg);
	//这里课件给的是minuit->Migrad()，效果可能是一样的
//	minuit -> mnexcm("HESSE", arglist, 0, ierflg);
	//输出结果
	Double_t fmin, fedm, errdef, cov_matrix[npar][npar];
	Double_t alpha, beta, alpha_err, beta_err;
	Int_t nvpar, nparx, icstat;
	minuit -> mnstat(fmin, fedm, errdef, nvpar, nparx, icstat);
	minuit -> GetParameter(0, alpha, alpha_err);
	minuit -> GetParameter(1, beta, beta_err);
	//注意这里的几个函数mnstat，Getparameter都是赋值函数，把计算结果给这些变量
	minuit -> mnemat(&cov_matrix[0][0], npar);	//协方差矩阵，“&”意思是说从这个元素开始依次赋值的意思
	Double_t rho = cov_matrix[0][1] / (alpha_err * beta_err);	//alpha、beta相关系数
	cout << "alpha =" << alpha << "+/-" << alpha_err << endl;
	cout << "beta =" << beta << "+/-" << beta_err << endl;
	cout << "cov[alpha][beta] =" << cov_matrix[0][1] << endl;
	cout << "rho[alpha][beta] =" << rho << endl;
	cout << "Log L =" << fmin << endl;
}


