#include <iostream>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include "TFile.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"


int addDQMFiles(TString step = "AB"){

  TH2F* Charge_Vs_Index=0;
  TH2F* together=0;
  int i=0;
  int run;
  TFile *f;
  TString strFile;
  ifstream inputFile("runlist" + step + ".txt");
  while(inputFile>>run){
    cout<<"run = "<<run<<endl;
    strFile.Form("DQM_V0001_R000%i__Express__PCLTest__ALCAPROMPT.root",run);
    cout<<"Open "<<strFile<<endl;
    f = TFile::Open(strFile);
    if(!f){
      cout<<"continue"<<endl;
      continue;
    }
    TString strHisto;
    strHisto.Form("DQMData/Run %i/AlCaReco/Run summary/SiStripGains/Charge_Vs_Index;1",run);
    cout<<"Get "<<strHisto<<endl;
    if(i==0){
      f->GetObject(strHisto,together);
      together->SetDirectory(0);
      if(!together){
	cout<<"No histogram available"<<endl;
	continue;
      }
    }
    else{
      f->GetObject(strHisto,Charge_Vs_Index);
      together->Add(Charge_Vs_Index,1);
    }
    delete f;
    i++;
  }

  TFile fSave("DQM_" + step + ".root","recreate");
  together->Write();
  fSave.Close();

  return 0;
}
