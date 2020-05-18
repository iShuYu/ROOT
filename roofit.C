#include <iostream>
#include <fstream>
#include <RooAbsReal.h>
#include <TTree.h>
#include <math.h>
#include <string>
#include <RooFit.h>
#include <RooRealVar.h>
#include <TROOT.h>
#include <RooGaussian.h>
#include <TFile.h>
#include <RooDataSet.h>
using namespace std;
using namespace RooFit;
int main(int argc, char* argv[])
{
	//INPUT:(use flag as a switch)
	ifstream config(argv[1]);
	int num_experiments, flag, N_entries;
	string output;
	config >> N_entries >> num_experiments >> output >> flag;

	//OUTPUT:
	TFile *FO = new TFile(output.c_str(),"recreate");
	TTree TO("tree","");
	double sigma, mean, likelihood, lambda;
	TO.Branch("sigma", &sigma, "sigma/D");
	TO.Branch("mean", &mean, "mean/D");
	TO.Branch("likelihood", &likelihood, "likelihood/D");
	TO.Branch("lambda",&lambda,"lambda/D");
	
	//FIT MODEL and DATA MODEL:
	RooRealVar x("x", "variable", -3., 3);
	RooRealVar Mean("Mean", "Mean Value", 0., -1., 1.);
	RooRealVar Sigma("Sigma", "Variance", 1., 0., 3.);
	RooRealVar DM("DM", "DATA MEAN", 0.);
	RooRealVar DS("DS", "DATA SIGMA", 1.);
	RooGaussian Model("Model", "gauss pdf model", x, Mean, Sigma);
	RooGaussian Model_D("MD", "", x, DM, DS);
	
	//SWITCH
	if (flag == 1) 
	{
		Mean.setVal(0.05);
		Mean.setConstant(1);
	}
	if (flag == 2) 
        {
                Mean.setVal(0.2);
                Mean.setConstant(1);
        }
	if (flag == 3) 
        {
                Mean.setVal(0.5);
                Mean.setConstant(1);
        }
	//LOOP:
	for (int i=0; i<num_experiments; i++)
	{
		//DATA SET
		RooDataSet *data = Model_D.generate(x,N_entries);

		//FIT
		Model.fitTo(*data);
		
		//SAVE THE RESULT
		TO.GetEntry(i);
		mean = Mean.getVal();
		sigma = Sigma.getVal();
		RooAbsReal *nll = Model.createNLL(*data);
		likelihood = nll -> getVal();
		RooAbsReal *pll = nll -> createProfile(Mean);
		lambda = (pll -> getVal()) / (nll -> getVal());
		TO.Fill();
	}

	TO.Write();
	FO -> Close();
}




