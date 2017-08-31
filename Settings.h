#ifndef SETTINGS
#define SETTINGS

#include <string>
#include <iostream>
#include "TH1F.h"
#include "TMutex.h"

class Settings{
  public:

  std::string ppInputFile = "/mnt/hadoop/cms/store/user/luck/2015-Data-promptRECO-photonSkims/pp-photonHLTFilter-v0-HiForest/0.root";
  //std::string PbPbInputFile = "/mnt/hadoop/cms/store/user/richard/2015-Data-promptRECO-photonSkims/HIPhoton40AndZ/PbPb-photonHLTFilter-v3/160202_145715/000*/*.root";

  bool isPP = true;
  int nEvtToProcess = -1;//set to -1 to get full data sample

  int jetRadiusTimes10 = 3;
  static const int nSoftDrops = 5;
  float subjetDRCut = 0.1;
  float softDropZ[nSoftDrops] = {0.1,0.2,0.5,0.1,0.2};
  float softDropBeta[nSoftDrops] = {0,0,1.5,1,1};

  std::string phoTrigger = "HLT_HISinglePhoton40_Eta1p5_v1";
  float phoPtCut = 60;
  float phoEtaCut = 1.44;
  float jetPtCut = 30;
  float jetEtaCut = 1.6;
  float dPhiCut = 2.74889357; //7pi/8


  //for controlling threads
  std::pair< float, float > photon[2];//stores phi,pt for passing between threads
  std::vector< std::vector< float > > b2bJets[2];//stores eta,phi,pt of back to back jets for passing between threads
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
