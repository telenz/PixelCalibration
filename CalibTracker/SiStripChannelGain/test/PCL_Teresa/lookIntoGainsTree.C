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
#include "TMath.h"


int lookIntoGainsTree(const TString filename="Gains_Tree_AB_CL0.root")
{

  //TCanvas* canvas;
  TH1F *NEntries;
  TTree* tree;
  TFile *_file = TFile::Open(filename,"READ");
  if(!_file) cout<<"No file available!"<<endl; 
  _file->GetObject("alcaSiStripGainsHarvester/APVGain",tree);
  if(!tree) cout<<"No tree available!"<<endl; 
  tree->Draw("NEntries");
  double mean = TMath::Mean(tree->GetSelectedRows(),tree->GetV1()); 
  cout<<"Mean of "<<filename<<" = "<<mean<<endl<<endl;
  
 
  return 0;
}
