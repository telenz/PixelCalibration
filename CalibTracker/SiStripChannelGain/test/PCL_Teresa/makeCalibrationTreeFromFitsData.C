#include <iostream>
#include  <cstdio>
#include "TFile.h"
#include "TString.h"
#include "TROOT.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"




//const Int_t nOfStripModules = 15148;
const Int_t nOfPixelModules   = 1440;
const Int_t nOfModules        = 16588;
const double globalPixelData  = 3.5090;
const double globalStrip      = 3.3027;

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

  TString Cstep[5];
  Cstep[0]="AB_CL2";
  Cstep[1]="C1_CL2";
  Cstep[2]="C2_CL2";
  Cstep[3]="D1_CL2";
  Cstep[4]="D2_CL2";

  Cstep[0]="AB";
  Cstep[1]="C1";
  Cstep[2]="C2";
  Cstep[3]="D1";
  Cstep[4]="D2";

  int runRangeStart[5];
  runRangeStart[0] = 190456; 
  runRangeStart[1] = 198022; 
  runRangeStart[2] = 199356;
  runRangeStart[3] = 203777; 
  runRangeStart[4] = 205826; 
  int runRangeEnd[5];
  runRangeEnd[0] = 196531;
  runRangeEnd[1] = 199336;
  runRangeEnd[2] = 203742;
  runRangeEnd[3] = 205781;
  runRangeEnd[4] = 208686;

  TFile f("Data8TeVGains.root","recreate");
  TDirectory *dir = f.mkdir("SiStripCalib");
  dir->cd();
  TTree t2("APVGain","a Tree with Pixel Calibration Factors");
  //  TTree t2("SiStripCalib/APVGain","a Tree with Pixel Calibration Factors");
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

  for(int i=0; i<5; i++){

    // Read input tree
    TChain* chain = new TChain("PixelGain");
    chain->Add("PixelGain_" + Cstep[i] + ".root");

    UInt_t   detid[nOfPixelModules];
    UChar_t  subdet[nOfPixelModules];
    Double_t gain[nOfPixelModules];
    Double_t width[nOfPixelModules];

    chain->SetBranchAddress("DetId",detid);
    chain->SetBranchAddress("SubDet",subdet);
    chain->SetBranchAddress("Gain",gain);
    chain->SetBranchAddress("FitWidth",width);

    // Read input tree
    TChain* chainStrip = new TChain("alcaSiStripGainsHarvester/APVGain");
    chainStrip->Add("Gains_Tree_AB_CL2_New.root");

    UInt_t   detidStrip;
    UChar_t  apvidStrip;
    UInt_t   indexStrip;
    UChar_t  subdetStrip;
    Double_t gainStrip;
    Double_t widthStrip;

    chainStrip->SetBranchAddress("DetId",&detidStrip);
    chainStrip->SetBranchAddress("Index",&indexStrip);
    chainStrip->SetBranchAddress("SubDet",&subdetStrip);
    chainStrip->SetBranchAddress("Gain",&gainStrip);
    chainStrip->SetBranchAddress("FitWidth",&widthStrip);
    chainStrip->SetBranchAddress("APVId",&apvidStrip);
  
  
    
    step.firstRun = runRangeStart[i];
    step.lastRun  = runRangeEnd[i];

    int countPixelModules = 0;

    int nMax = chain->GetEntries();

    int idx=0;
    for(int n=0; n<nMax; n++)
      {
	chain->GetEntry(n);

	for(int j=0;j<nOfPixelModules;j++){
	  
	  step.DetId[idx]                 = detid[j];
	  step.SubDet[idx]                = subdet[j];
	  step.APVId[idx]                 = 0;
	  step.Gain[idx]                  = gain[j]*globalPixelData/globalStrip;
	  step.GainWoGlobalScaling[idx]   = gain[j];
	  step.FitWidth[idx]              = width[j];
	  step.Index[idx]                 = idx;
	  idx += 1;
	  countPixelModules += 1;
	}

      }



    for(int n=0; n<chainStrip->GetEntries(); n++)
      {
	chainStrip->GetEntry(n);


	  
	if((int)subdetStrip<3 || (int)apvidStrip!=0) continue;

	  step.DetId[idx]                 = detidStrip;
	  step.SubDet[idx]                = subdetStrip;
	  step.APVId[idx]                 = 0;
	  step.Gain[idx]                  = 1.;
	  step.GainWoGlobalScaling[idx]   = 1.;
	  step.FitWidth[idx]              = 1;
	  step.Index[idx]                 = idx;
	  idx += 1;
	  countPixelModules += 1;
	  //cout<<"idx = "<<idx<<endl;
	  //cout<<"gain["<<j<<"] = "<<gain[j]<<endl;


      }


    // fill the Tree with current step parameters
    
    t2.Fill();
  
  }
 
  //save the Tree header. The file will be automatically closed
  //when going out of the function scope
  t2.Write();
  
  return 0;


}
