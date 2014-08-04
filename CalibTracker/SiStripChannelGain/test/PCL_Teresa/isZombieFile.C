#include "TFile.h"

int isZombieFile(TString path){
  
  TFile file(path);
  if (file.IsZombie()) {
    return 1;
  }
  
  return 0;
  
}
