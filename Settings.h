#ifndef SETTINGS
#define SETTINGS

#include <string>
#include <iostream>
#include "TH1F.h"
#include "TMutex.h"

class Settings{
  public:

  std::string inputFile = "/mnt/hadoop/cms/store/user/luck/2015-Data-promptRECO-photonSkims/pp-photonHLTFilter-v0-HiForest/0.root";

  int jetRadiusTimes10 = 3;

  float phoPtCut = 60;
  float phoEtaCut = 1.44;
  float jetPtCut = 30;
  float jetEtaCut = 2;


  //for controlling threads
  
  int nextEvt[2] = {-1,-1};
  bool isQueueFull = false;
  bool isEvtThreadDone = false;
  
  Settings();

  private:

};

Settings::Settings(){
  std::cout << "Initializing Settings!" << std::endl; 
}
#endif
