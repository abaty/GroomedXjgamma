#include "Settings.h"
#include "include/Tools.h"
#include "fastjet/ClusterSequence.hh"
#include "include/softDropGroomer.hh"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TMath.h"
#include "TThread.h"
//#include "TMutex.h"
#include "TCondition.h"
#include <vector>
#include <string>
#include <iostream>
#include <stdlib.h> 

using namespace fastjet;

TMutex mtx1;

void *findEvts(void *ptr){
  Settings * s = (Settings*) ptr;
  TFile * f = TFile::Open(s->inputFile.c_str(),"read");
  TTree * skim = (TTree*)f->Get("skimanalysis/HltTree");
  TTree * photons = (TTree*)f->Get("ggHiNtuplizerGED/EventTree");  
  TTree * evt = (TTree*)f->Get("hiEvtAnalyzer/HiTree");
  TTree * hlt = (TTree*)f->Get("hltanalysis/HltTree");
  std::cout << Form("ak%dPFJetAnalyzer/t",s->jetRadiusTimes10) << std::endl;
  TTree * jet = (TTree*)f->Get(Form("ak%dPFJetAnalyzer/t",s->jetRadiusTimes10));
  skim->AddFriend(evt); 

  int pVtx = 0;
  int beam = 0;
  int HBHE = 0;
  int trig = 0; 
  float vz = -99;
 
  std::vector< float > * seedTime = 0;
  std::vector< float > * swissCrx = 0;
  std::vector< float > * sigmaIEtaIEta = 0;
  std::vector< float > * ecalR4 = 0;
  std::vector< float > * hcalR4 = 0;
  std::vector< float > * trkR4 = 0;
  std::vector< float > * hOverE = 0;
  std::vector< float > * phoEt = 0;
  std::vector< float > * phoPhi = 0;
  std::vector< float > * phoEta = 0;
  float jtpt[1000];
  float jtphi[1000];
  float jteta[1000];
  int nref;

  skim->SetBranchAddress("pPAprimaryVertexFilter",&pVtx);
  skim->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&HBHE);
  skim->SetBranchAddress("pBeamScrapingFilter",&beam);
  photons->SetBranchAddress("pho_seedTime",&seedTime);
  photons->SetBranchAddress("pho_swissCrx",&swissCrx);
  photons->SetBranchAddress("phoSigmaIEtaIEta",&sigmaIEtaIEta);
  photons->SetBranchAddress("pho_ecalClusterIsoR4",&ecalR4);
  photons->SetBranchAddress("pho_hcalRechitIsoR4",&hcalR4);
  photons->SetBranchAddress("pho_trackIsoR4PtCut20",&trkR4);
  photons->SetBranchAddress("phoHoverE",&hOverE);
  photons->SetBranchAddress("phoEt",&phoEt);
  photons->SetBranchAddress("phoPhi",&phoPhi);
  photons->SetBranchAddress("phoEta",&phoEta);
  evt->SetBranchAddress("vz",&vz);
  hlt->SetBranchAddress("HLT_HISinglePhoton50_Eta3p1_v1",&trig);
  jet->SetBranchAddress("nref",&nref);
  jet->SetBranchAddress("jtpt",&jtpt);
  jet->SetBranchAddress("jtphi",&jtphi);
  jet->SetBranchAddress("jteta",&jteta);


  for(int i = 0; i<100000; i++){//i<skim->GetEntries(); i++){
    if(i%10000 == 0) std::cout << i << "/" << skim->GetEntries() << std::endl;
    hlt->GetEntry(i);
    if(!trig) continue;
    skim->GetEntry(i);
    if(!(pVtx && beam && HBHE)) continue;
    if(TMath::Abs(vz)>15) continue;

    photons->GetEntry(i);
    std::vector< float > goodPhotonPhi;
    for(unsigned int j = 0; j<phoEt->size(); j++){
      if(phoEt->at(j) < s->phoPtCut) continue;
      if(TMath::Abs(phoEta->at(j)) > s->phoEtaCut) continue;
      if(!(TMath::Abs(seedTime->at(j))<3 && swissCrx->at(j)<0.9 && sigmaIEtaIEta->at(j)>0.002)) continue;//spike rejection
      if(!((ecalR4->at(j)+hcalR4->at(j)+trkR4->at(j)) < 1 && (hOverE->at(j)<0.1))) continue;//isolation
      if(!(sigmaIEtaIEta->at(j)<0.01)) continue;//sigmaietaieta
      goodPhotonPhi.push_back(phoPhi->at(j));
    }
    if(goodPhotonPhi.size()==0) continue;

    jet->GetEntry(i);
    bool hasBackToBackJet = false;
    //caution, assumes jet tree is ordered by pt for speed!
    for(unsigned int j = 0; j<nref; j++){
      if(jtpt[j] < s->jetPtCut) break;
      if(TMath::Abs(jteta[j]) > s->jetEtaCut) continue;
      for(unsigned int p = 0; p<goodPhotonPhi.size(); p++){
        if(TMath::ACos(TMath::Cos(goodPhotonPhi.at(p)-jtphi[j])) > 5*TMath::Pi()/6.0) hasBackToBackJet = true; 
      } 
    }
    if(!hasBackToBackJet) continue;

    //put cuts above here
    //controls placing event index into a buffer accesible by the other thread
    mtx1.Lock();
      while(s->isQueueFull){
        mtx1.UnLock();
        sleep(0.1);
        //std::cout << "Thread 2 waiting" << std::endl;
        mtx1.Lock();
      }
      //std::cout << "placing evt " << i << std::endl;
      if(s->nextEvt[0] == -1) s->nextEvt[0] = i;
      else if(s->nextEvt[1] == -1){
        s->nextEvt[1] = i;
        s->isQueueFull = true;
      } 
    mtx1.UnLock();    
  }
  std::cout << "Thread 1 Done!" << std::endl;
  mtx1.Lock();
    s->isEvtThreadDone = true;
  mtx1.UnLock();
}


void runGrooming(){
  Settings s = Settings();
  std::cout << "Initializing evt search thread" << std::endl;
 
  //setup



  //thread setup and control
  TThread * thread1 = new TThread("thread1",findEvts, (void*)(&s));
  thread1->Run();
  bool isDone = false;

  TFile * f = TFile::Open(s.inputFile.c_str(),"read");
  TTree * pf = (TTree*)f->Get("pfcandAnalyzer/pfTree");

  std::vector< int > * id = 0;
  std::vector< float > * pt = 0;
  std::vector< float > * eta = 0;
  std::vector< float > * phi = 0;
  std::vector < float > * e = 0;

  pf->SetBranchAddress("pfId",&id);
  pf->SetBranchAddress("pfPt",&pt);
  pf->SetBranchAddress("pfEta",&eta);
  pf->SetBranchAddress("pfPhi",&phi);
  pf->SetBranchAddress("pfEnergy",&e);

  while(true){
    sleep(0.01);
    //std::cout << "Thread 1 waiting" << std::endl;
    mtx1.Lock();
      if(s.nextEvt[0] != -1){
        int thisEvt = s.nextEvt[0];
        s.nextEvt[0] = s.nextEvt[1];
        s.nextEvt[1] = -1;
        s.isQueueFull = false;
        if(s.isEvtThreadDone && s.nextEvt[0]==-1) isDone = true;
        mtx1.UnLock();

        //analysis code here
        pf->GetEntry(thisEvt);
        std::vector<PseudoJet> particles;
        for(unsigned int i = 0; i<id->size(); i++){
          if(TMath::Abs(eta->at(i))>2.9) continue;
          particles.push_back( PseudoJet(pt->at(i)*TMath::Cos(phi->at(i)),pt->at(i)*TMath::Sin(phi->at(i)),pt->at(i)/(TMath::Tan(2*TMath::ATan(TMath::Exp(-eta->at(i))))),e->at(i)));//FIXME change to TLorentzVector constructor
        }

        JetDefinition jet_def(antikt_algorithm, s.jetRadiusTimes10/10.0);
        ClusterSequence cs(particles, jet_def);
        std::vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(30.0));
        std::cout <<   "        pt y phi" << std::endl;
        for (unsigned i = 0; i < jets.size(); i++) {
          std::cout << "jet " << i << ": "<< jets[i].pt() << " " 
                   << jets[i].rap() << " " << jets[i].phi() << std::endl;
        }

        //std::cout << "done with " << thisEvt << std::endl;

        //end analysis code
        if(isDone) break;
      }else{
        if(s.isEvtThreadDone && s.nextEvt[0]==-1) isDone = true;
        mtx1.UnLock();
        if(isDone) break;
      }

  }
  thread1->Join();
}

int main(){
  runGrooming();
  return 1;
}
