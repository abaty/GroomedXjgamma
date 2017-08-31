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
#include "TLorentzVector.h"
#include <vector>
#include <string>
#include <iostream>
#include <stdlib.h> 

using namespace fastjet;

TMutex mtx1;

void *findEvts(void *ptr){
  Settings * s = (Settings*) ptr;
  TFile * f;
  if(s->isPP) f = TFile::Open(s->ppInputFile.c_str(),"read");
  else       f = TFile::Open(s->ppInputFile.c_str(),"read");
  TTree * skim = (TTree*)f->Get("skimanalysis/HltTree");
  TTree * photons; 
  if(s->isPP) photons = (TTree*)f->Get("ggHiNtuplizerGED/EventTree");  
  else        photons = (TTree*)f->Get("ggHiNtuplizer/EventTree");  
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
  hlt->SetBranchAddress(s->phoTrigger.c_str(),&trig);
  jet->SetBranchAddress("nref",&nref);
  jet->SetBranchAddress("jtpt",&jtpt);
  jet->SetBranchAddress("jtphi",&jtphi);
  jet->SetBranchAddress("jteta",&jteta);


  for(int i = 0; i< ((s->nEvtToProcess<0)?skim->GetEntries():s->nEvtToProcess); i++){
    if(i%10000 == 0) std::cout << i << "/" << skim->GetEntries() << std::endl;
    hlt->GetEntry(i);
    if(!trig) continue;
    skim->GetEntry(i);
    if(!(pVtx && beam && HBHE)) continue;
    if(TMath::Abs(vz)>15) continue;

    photons->GetEntry(i);
    std::vector< std::pair< float, float > > goodPhotons;
    for(unsigned int j = 0; j<phoEt->size(); j++){
      if(phoEt->at(j) < s->phoPtCut) continue;
      if(TMath::Abs(phoEta->at(j)) > s->phoEtaCut) continue;
      if(!(TMath::Abs(seedTime->at(j))<3 && swissCrx->at(j)<0.9 && sigmaIEtaIEta->at(j)>0.002)) continue;//spike rejection
      if(!((ecalR4->at(j)+hcalR4->at(j)+trkR4->at(j)) < 1 && (hOverE->at(j)<0.1))) continue;//isolation
      if(!(sigmaIEtaIEta->at(j)<0.01)) continue;//sigmaietaieta
      goodPhotons.push_back(std::pair<float, float>(phoPhi->at(j), phoEt->at(j)));
    }
    if(goodPhotons.size()==0) continue;

    jet->GetEntry(i);
    bool hasBackToBackJet = false;
    int photonIndx = -1;
    std::vector< std::vector< float > > jetVars;
    //caution, assumes jet tree is ordered by pt for speed!
    for(unsigned int j = 0; j<nref; j++){
      if(jtpt[j] < s->jetPtCut) break;
      if(TMath::Abs(jteta[j]) > s->jetEtaCut) continue;
      for(unsigned int p = 0; p<goodPhotons.size(); p++){
        if(TMath::ACos(TMath::Cos(goodPhotons.at(p).first-jtphi[j])) > s->dPhiCut){
          hasBackToBackJet = true; 
          photonIndx = p;
          std::vector< float > tempJet;
          tempJet.push_back(jteta[j]);
          tempJet.push_back(jtphi[j]);
          tempJet.push_back(jtpt[j]);
          jetVars.push_back(tempJet);
          continue;
        }
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
    if(s->nextEvt[0] == -1){
      s->nextEvt[0] = i;
      (s->photon)[0] = goodPhotons.at(photonIndx);
      (s->b2bJets)[0] = jetVars;
    } else if(s->nextEvt[1] == -1){
      s->nextEvt[1] = i;
      (s->photon)[1] = goodPhotons.at(photonIndx);
      (s->b2bJets)[1] = jetVars;
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
  TH1::SetDefaultSumw2();
  Settings s = Settings();
  std::cout << "Initializing evt search thread" << std::endl;
 
  //thread setup and control
  TThread * thread1 = new TThread("thread1",findEvts, (void*)(&s));
  thread1->Run();
  bool isDone = false;

  TFile * out = TFile::Open("output.root","recreate");
  TH1D * xjg_LT[s.nSoftDrops];
  TH1D * xjg_GT[s.nSoftDrops];
  for(int i = 0; i<s.nSoftDrops; i++){
    xjg_LT[i] = new TH1D(Form("xjg_lt_%d",i),Form("xjg_lt_%d",i),20,0,2);
    xjg_GT[i] = new TH1D(Form("xjg_gt_%d",i),Form("xjg_gt_%d",i),20,0,2);
  } 

  //loading pfcands
  TFile * f;
  if(s.isPP) f = TFile::Open(s->ppInputFile.c_str(),"read");
  else       f = TFile::Open(s->ppInputFile.c_str(),"read");
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

  int nPho = 0;
  float photonPhi = -1;
  float photonPt = -1;
  std::vector<  std::vector< float > > jetVars;//contains a list of jets w/ eta/phi/pt


  while(true){
    sleep(0.01);
    //std::cout << "Thread 1 waiting" << std::endl;
    mtx1.Lock();
      if(s.nextEvt[0] != -1){
        int thisEvt = s.nextEvt[0];
        photonPhi = (s.photon)[0].first;
        photonPt = (s.photon)[0].second;
        jetVars = (s.b2bJets)[0];
        s.nextEvt[0] = s.nextEvt[1];
        (s.photon)[0] = (s.photon)[1];
        (s.b2bJets)[0] = (s.b2bJets)[1];
        s.nextEvt[1] = -1; 
        s.isQueueFull = false;
        if(s.isEvtThreadDone && s.nextEvt[0]==-1) isDone = true;
        mtx1.UnLock();

        //analysis code here
        std::cout << "\nPhoton Pt: " << photonPt << " Photon Phi: " << photonPhi << std::endl;

        //pfcandSetup for fastJet
        pf->GetEntry(thisEvt);
        nPho++;
        std::vector<PseudoJet> particles;
        for(unsigned int i = 0; i<id->size(); i++){
          TLorentzVector vec = TLorentzVector();
          vec.SetPtEtaPhiE(pt->at(i),eta->at(i), phi->at(i), e->at(i));
          particles.push_back( PseudoJet(vec.Px(),vec.Py(),vec.Pz(), vec.E()) );
        }

        //clustering algorithm
        JetDefinition jet_def(antikt_algorithm, s.jetRadiusTimes10/10.0);
        ClusterSequence cs(particles, jet_def);
        std::vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

        //matching and SoftDrop
        std::vector<PseudoJet> matchedJets;
        softDropGroomer sd[s.nSoftDrops];
        for(int k = 0; k<s.nSoftDrops; k++) sd[k] = softDropGroomer(s.softDropZ[k], s.softDropBeta[k], s.jetRadiusTimes10/10.0);

        for(unsigned int j = 0; j < jetVars.size(); j++){
          bool isMatched = false;
          std::cout << "Forest Jet " << j << " Pt: " << jetVars.at(j).at(2) << " Eta: " << jetVars.at(j).at(0) 
                  << " Phi: " << jetVars.at(j).at(1) << std::endl;

          for (unsigned int i = 0; i < jets.size(); i++) {
            if(TMath::Abs(jets[i].eta()-jetVars.at(j).at(0))>0.05) continue;
            if(TMath::Abs(jets[i].phi_std()-jetVars.at(j).at(1))>0.05) continue;
            isMatched = true;
            std::cout << "FJ jet: " << i << ": "<< jets[i].pt() << " " 
                   << jets[i].eta() << " " << jets[i].phi_std() << std::endl;
            for(int k = 0; k<s.nSoftDrops; k++){
              std::vector<PseudoJet> groomedJets = sd[k].doGrooming(jets);

              if(sd[k].getDR12().at(i) <= s.subjetDRCut) xjg_LT[k]->Fill(jetVars.at(j).at(2)/photonPt);
              else xjg_GT[k]->Fill(jetVars.at(j).at(2)/photonPt);
            }
            continue; 
          }
          /*if(!isMatched){
            std::cout << "No Match!!!" << std::endl;
            std::cout << "Forest jet: " << j << ": "<< jetVars.at(j).at(2) << " "
                   << jetVars.at(j).at(0) << " " << jetVars.at(j).at(1) << std::endl;
            for (unsigned int i = 0; i < jets.size(); i++) {
              std::cout << "FJ jet: " << i << ": "<< jets[i].pt() << " "
                   << jets[i].eta() << " " << jets[i].phi() << std::endl;
            }
          }*/
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
  for(int i = 0; i<s.nSoftDrops; i++){
    xjg_LT[i]->Scale(1/(float)nPho);
    xjg_GT[i]->Scale(1/(float)nPho);
    xjg_GT[i]->SetLineColor(kRed);
  } 
  out->Write();
}

int main(){
  runGrooming();
  return 1;
}
