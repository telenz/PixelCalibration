#include <iostream>
#include  <cstdio>
#include "TProfile.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TObjString.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TROOT.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TSystemDirectory.h"
#include "TList.h"
#include "langaus.C"


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>

const Int_t nOfPixelModules = 1440;

typedef struct { 

  unsigned int nEntries;
  unsigned int firstRun;
  unsigned int lastRun;
  UInt_t   DetId[nOfPixelModules];
  UChar_t  SubDet[nOfPixelModules];
  Double_t Gain[nOfPixelModules]; 
  Double_t FitWidth[nOfPixelModules]; 
  Double_t FitChi2NDF[nOfPixelModules]; 

  
} CalibrationStep_t;


int makeFitsForCalibration(TString calibStep="mc")
{
 
  clock_t start, stop;
  double t = 0.0;

  UInt_t   detid;
  UChar_t  subdet;
  Int_t index;

  std::map<int,TH1F*> proj;
  std::map<int,UChar_t> subDetector;

  Double_t sv[4];
  TH2F *Charge_Vs_Index;


  TFile *_file0 = TFile::Open("DQM_" + calibStep + ".root");

  _file0->GetObject("Charge_Vs_Index",Charge_Vs_Index);
  
  
  TChain* chain = new TChain("alcaSiStripGainsHarvester/APVGain");
  // arbitrary just for detid and subdet and index
  chain->Add("Gains_Tree_C1_CL2_New.root");
  chain->SetBranchAddress("Index",&index);
  chain->SetBranchAddress("DetId",&detid);
  chain->SetBranchAddress("SubDet",&subdet);
  
  TFile fTree("PixelGain_" + calibStep + ".root","recreate");
  TTree t2("PixelGain","a Tree with Pixel Calibration Factors");
  CalibrationStep_t step;
  t2.Branch("firstRun",&step.firstRun,"firstRun/I");
  t2.Branch("lastRun",&step.lastRun,"lastRun/I");
  //t2.Branch("nEntries",&step.nEntries,"nEntries/I");
  TString str = (TString) "DetId[" + (long) nOfPixelModules + (TString) "]/i";
  t2.Branch("DetId",step.DetId,str);
  str = (TString) "SubDet[" + (long) nOfPixelModules + (TString) "]/b";
  t2.Branch("SubDet",step.SubDet,str);
  str = (TString) "Gain[" + (long) nOfPixelModules + (TString) "]/D";
  t2.Branch("Gain",step.Gain,str);
  str = (TString) "FitWidth[" + (long) nOfPixelModules + (TString) "]/D";
  t2.Branch("FitWidth",step.FitWidth,str);
  str = (TString) "FitChi2NDF[" + (long) nOfPixelModules + (TString) "]/D";
  t2.Branch("FitChi2NDF",step.FitChi2NDF,str);
 
  bool iniAllTogether=true;
  TH1F* aux=0;  
  for(int n=0; n<chain->GetEntries(); n++)
    {
      chain->GetEntry(n);
      
      if(subdet>2) continue;
      
      if(proj[detid] == 0){

	aux = (TH1F*) Charge_Vs_Index->ProjectionY("proj1",Charge_Vs_Index->GetXaxis()->FindBin(index),Charge_Vs_Index->GetXaxis()->FindBin(index),"e");

	proj[detid] = (TH1F*)aux->Clone("proj1");
	subDetector[detid] = subdet;
	

	if(iniAllTogether){
	  proj[1]=(TH1F*)aux->Clone("proj1");
	  iniAllTogether=false;
	}
	

      }
      else {

	aux=(TH1F*) Charge_Vs_Index->ProjectionY("proj1",Charge_Vs_Index->GetXaxis()->FindBin(index),Charge_Vs_Index->GetXaxis()->FindBin(index),"e");
	proj[detid]->Add(aux,1);
      }
      
      proj[1]->Add(aux,1);
      
    }

  // Fitting procedure
  TF1 *ffit = new TF1("ffit",langaufun,50,4000,4);
  ffit->SetLineColor(2);
  ffit->SetParLimits(1, 0, 1000); // Landau MPV
  ffit->SetParLimits(3, 0, 1000); // Gaussian Sigma

  sv[0]=25.; sv[1]=300; sv[2]=10000; sv[3]=35.;
  ffit->SetParameters(sv);  
    
  assert((start = clock())!=-1);
 
  std::map<int, TH1F*>::iterator projIter;
  int idx = 0;
  for(projIter = proj.begin(); projIter != proj.end(); projIter++) {

    TH1F* ex = projIter->second;
    
    if(ex->GetEntries()<25){

      step.DetId[idx]        = projIter->first;
      step.Gain[idx]         = 1.;
      step.FitWidth[idx]     = -1.;
      step.FitChi2NDF[idx]   = -1.;
      step.SubDet[idx]       = subDetector[projIter->first];
      idx+=1;
      continue;
      
    }

    if(ex->GetEntries()<2000.){
      cout<<"ex->GetEntries() = "<<ex->GetEntries()<<endl;
      ex->Rebin(2);
    }


    sv[0]=25.; sv[1]=ex->GetBinCenter(ex->GetMaximumBin()); sv[2]=ex->Integral(); sv[3]=35.;
    ffit->SetParameters(sv);  
    

    //####################### Fit ##########################################
    bool firstFit=true;
    int nFit=0;
    double lowRange=0;
    double highRange=0;
    ex->Fit("ffit","QB","",50,4000); 
    
    while((ffit->GetChisquare()/ffit->GetNDF() > 5 || firstFit || lowRange<0) && nFit<5){
      lowRange  = ffit->GetMaximumX()-1.5*(ffit->GetParameter(0)+ffit->GetParameter(3));
      highRange = ffit->GetMaximumX()+1.5*(ffit->GetParameter(0)+ffit->GetParameter(3));
      if(lowRange<=0)                      lowRange=100.;
      if(highRange>100000 || highRange<=0) highRange=500.;
      if(lowRange>ex->GetBinCenter(ex->GetMaximumBin())){
	lowRange  = ex->GetBinCenter(ex->GetMaximumBin())-100;
	cout<<"low range wrong!"<<endl;
	cout<<"detid = "<<projIter->first<<endl;
      }
      if(highRange<ex->GetBinCenter(ex->GetMaximumBin())){
	highRange = ex->GetBinCenter(ex->GetMaximumBin())+100;
	cout<<"high range wrong!"<<endl;
	cout<<"detid = "<<projIter->first<<endl;
      }
      ex->Fit("ffit","QB","",lowRange,highRange);  
      firstFit=false;
      nFit +=1;
    }
    cout<<"detid = "<<projIter->first<<endl;
    cout<<"nFit = "<<nFit<<endl;
    
    if(ffit->GetChisquare()/ffit->GetNDF() > 3 ){
 
      cout<<"Set new parameters"<<endl;
      sv[0]=25; sv[1]=ex->GetBinCenter(ex->GetMaximumBin()); sv[2]=ex->Integral(); sv[3]=35.;
      ffit->SetParameters(sv);  

      firstFit=true;
      nFit=0;
      ex->Fit("ffit","QB","",50,4000); 

      while((ffit->GetChisquare()/ffit->GetNDF() > 3 || firstFit || lowRange<0) && nFit<5){
	lowRange  = ffit->GetMaximumX()-0.8*(ffit->GetParameter(0)+ffit->GetParameter(3));
	highRange = ffit->GetMaximumX()+0.8*(ffit->GetParameter(0)+ffit->GetParameter(3));
	if(lowRange<=0)                      lowRange=100.;
	if(highRange>100000 || highRange<=0) highRange=500.;
	if(lowRange>ex->GetBinCenter(ex->GetMaximumBin())){
	  lowRange  = ex->GetBinCenter(ex->GetMaximumBin())-100;
	  cout<<"low range wrong!"<<endl;
	  cout<<"detid = "<<projIter->first<<endl;
	}
	if(highRange<ex->GetBinCenter(ex->GetMaximumBin())){
	  highRange = ex->GetBinCenter(ex->GetMaximumBin())+100;
	  cout<<"high range wrong!"<<endl;
	  cout<<"detid = "<<projIter->first<<endl;
	}
	ex->Fit("ffit","QB","",lowRange,highRange);  
	firstFit=false;
	nFit +=1;
      }

 
    }

    if(lowRange>ex->GetBinCenter(ex->GetMaximumBin())){
      lowRange  = ex->GetBinCenter(ex->GetMaximumBin())-100;
      cout<<"Fit did not work!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
    }
    if(highRange<ex->GetBinCenter(ex->GetMaximumBin())){
      highRange = ex->GetBinCenter(ex->GetMaximumBin())+100;
      cout<<"Fit did not work!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
    }
    

    /*
    cout<<"Landau Width      = "<<ffit->GetParameter(0)<<endl;
    cout<<"MPV               = "<<ffit->GetParameter(1)<<endl;
    cout<<"Maximum           = "<<ffit->GetMaximumX()<<endl;
    cout<<"Area              = "<<ffit->GetParameter(2)<<endl;
    cout<<"Gaussian sigma    = "<<ffit->GetParameter(3)<<endl;
    cout<<"MPV/300           = "<<ffit->GetParameter(1)/300<<endl;
    cout<<"Max/300           = "<<ffit->GetMaximumX()/300<<endl;
    cout<<"max bin           = "<<ex->GetBinCenter(ex->GetMaximumBin())<<endl;
    cout<<"Chi2/ndof         = "<<ffit->GetChisquare()/ffit->GetNDF()<<endl<<endl;
    */


    if(projIter->first!=1){
      
      step.DetId[idx]      = projIter->first;
      step.Gain[idx]       = ffit->GetMaximumX()/300.;
      step.FitWidth[idx]   = ffit->GetParameter(0)/300.;
      step.FitChi2NDF[idx] = ffit->GetChisquare()/ffit->GetNDF();
      step.SubDet[idx]     = subDetector[projIter->first];
      idx+=1;
    }
   
    if(ffit->GetChisquare()/ffit->GetNDF()>5){
      
      cout<<"Not a good fit!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
      cout<<projIter->first<<endl;
      cout<<"Landau Width      = "<<ffit->GetParameter(0)<<endl;
      cout<<"MPV               = "<<ffit->GetParameter(1)<<endl;
      cout<<"Maximum           = "<<ffit->GetMaximumX()<<endl;
      cout<<"Area              = "<<ffit->GetParameter(2)<<endl;
      cout<<"Gaussian sigma    = "<<ffit->GetParameter(3)<<endl;
      cout<<"MPV/300           = "<<ffit->GetParameter(1)/300<<endl;
      cout<<"Max/300           = "<<ffit->GetMaximumX()/300<<endl;
      cout<<"Chi2/ndof         = "<<ffit->GetChisquare()/ffit->GetNDF()<<endl<<endl;
    }
    //##################################################################################
  
    TFile f(Form("fits/" + calibStep + "/%i.root",projIter->first),"recreate");
    ex->Write();
    f.Close();
    
  }

  t2.Fill();

  stop = clock();
  t = (double) (stop-start)/CLOCKS_PER_SEC;
  printf("Run time: %f\n", t);

  fTree.cd();
  t2.Write();

  return 0;
}
