#include <iostream>
#include  <cstdio>
#include "TFile.h"
#include "TString.h"
#include "TROOT.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TH1D.h"
#include "TF1.h"


//const Int_t nOfStripModules = 15148;
const Int_t nOfPixelModules   = 1440;
const Int_t nOfModules        = 16588;
const double globalPixelMC    = 3.50843; // <- modulewise, 1 factor:3.48222;
const double globalStrip      = 3.30271;

typedef struct { 

  unsigned int nEntries;
  unsigned int firstRun;
  unsigned int lastRun;
  unsigned int Index[nOfModules];
  UInt_t       DetId[nOfModules];
  UChar_t      SubDet[nOfModules];
  UChar_t      APVId[nOfModules];
  Double_t     Gain[nOfModules]; 
  Double_t     GainWoGlobalScaling[nOfModules];
  Double_t     FitWidth[nOfModules]; 
  
} CalibrationStep_t;


int makeCalibrationTree()
{

  TString Cstep[1];

  int runRangeStart = 0;
  int runRangeEnd   = 30000000;

  TFile f("MC8TeVGains.root","recreate");
  TDirectory *dir = f.mkdir("SiStripCalib");
  dir->cd();
  TTree t2("APVGain","a Tree with Pixel Calibration Factors");

  //  TFile f("MC8TeVGains.root","recreate");
  //TTree t2("PixelGain","a Tree with Pixel Calibration Factors");
  CalibrationStep_t step;

  t2.Branch("firstRun",&step.firstRun,"firstRun/I");
  t2.Branch("lastRun",&step.lastRun,"lastRun/I");
  TString str = (TString) "DetId[" + (long) nOfModules + (TString) "]/i";
  t2.Branch("DetId",step.DetId,str);
  str = (TString) "Index[" + (long) nOfModules + (TString) "]/i";
  t2.Branch("Index",step.Index,str);
  str = (TString) "SubDet[" + (long) nOfModules + (TString) "]/b";
  t2.Branch("SubDet",step.SubDet,str);
  str = (TString) "APVId[" + (long) nOfModules + (TString) "]/b";
  t2.Branch("APVId",step.APVId,str);
  str = (TString) "Gain[" + (long) nOfModules + (TString) "]/D";
  t2.Branch("Gain",step.Gain,str);
  str = (TString) "GainWoGlobalScaling[" + (long) nOfModules + (TString) "]/D";
  t2.Branch("GainWoGlobalScaling",step.GainWoGlobalScaling,str);
  str = (TString) "FitWidth[" + (long) nOfModules + (TString) "]/D";
  t2.Branch("FitWidth",step.FitWidth,str);

  std::map<int,double> stripModules;  

  // Read input tree
  TChain* chain = new TChain("PixelGain");
  chain->Add("PixelGain_MC.root");

  UInt_t   detid[nOfPixelModules];
  UChar_t  subdet[nOfPixelModules];
  Double_t gain[nOfPixelModules];
  Double_t width[nOfPixelModules];

  chain->SetBranchAddress("DetId",detid);
  chain->SetBranchAddress("SubDet",subdet);
  chain->SetBranchAddress("Gain",gain);
  chain->SetBranchAddress("FitWidth",width);
    
  step.firstRun = runRangeStart;
  step.lastRun  = runRangeEnd;

  int idx=0;
  for(int n=0; n<chain->GetEntries(); n++)
      {
	chain->GetEntry(n);
	
	for(int j=0;j<nOfPixelModules;j++){
	  
	  step.DetId[idx]                 = detid[j];
	  step.SubDet[idx]                = subdet[j];
	  step.APVId[idx]                 = 0;
	  step.Gain[idx]                  = gain[j]*globalPixelMC/globalStrip;
	  step.GainWoGlobalScaling[idx]   = gain[j];
	  step.FitWidth[idx]              = width[j];
	  step.Index[idx]                 = idx;
	  idx += 1;
	}
      }


  // Read input tree
  TChain* chainStrip = new TChain("SiStripCalib/APVGain;2");
  chainStrip->Add("MC7TeVGains.root");
  
  UInt_t   detidStrip;
  UChar_t  apvidStrip;
  UInt_t   indexStrip;
  UChar_t  subdetStrip;
  Double_t gainStrip;
  Float_t widthStrip;

  chainStrip->SetBranchAddress("DetId"    , &detidStrip);
  chainStrip->SetBranchAddress("Index"    , &indexStrip);
  chainStrip->SetBranchAddress("SubDet"   , &subdetStrip);
  chainStrip->SetBranchAddress("Gain"     , &gainStrip);
  chainStrip->SetBranchAddress("FitWidth" , &widthStrip);
  chainStrip->SetBranchAddress("APVId"    , &apvidStrip);
     
  for(int n=0; n<chainStrip->GetEntries(); n++)
    {
      chainStrip->GetEntry(n);
      
      
      if((int)subdetStrip<3) continue;
      
      stripModules[detidStrip] = gainStrip;
    }
  
  
  for(int n=0; n<chainStrip->GetEntries(); n++)
    {
      chainStrip->GetEntry(n);
      
      if((int)subdetStrip<3 || apvidStrip!=0) continue;

      step.DetId[idx]                 = detidStrip;
      step.SubDet[idx]                = subdetStrip;
      step.APVId[idx]                 = apvidStrip;
      step.Gain[idx]                  = stripModules[detidStrip];
      step.GainWoGlobalScaling[idx]   = stripModules[detidStrip];
      step.FitWidth[idx]              = widthStrip/300.;
      step.Index[idx]                 = idx;
      idx += 1;

    }

  // fill the Tree with current step parameters
  t2.Fill();

  

 
  //save the Tree header. The file will be automatically closed
  //when going out of the function scope
  t2.Write();
  
  return 0;


}
