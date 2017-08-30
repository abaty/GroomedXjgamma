#include "fastjet/ClusterSequence.hh"
#include "include/softDropGroomer.hh"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TMath.h"
#include <iostream>
#include <vector>
using namespace fastjet;

void runFJ(){
//  TFile * f = TFile::Open("/data/cmcginn/HighPt10Skims/outSkimFile_99p2Percent.root","read");//dijet tree
  TFile * f = TFile::Open("/mnt/hadoop/cms/store/user/luck/2015-Data-promptRECO-photonSkims/pp-photonHLTFilter-v0-HiForest/0.root","read");
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

  pf->GetEntry(1);
  std::cout << id->at(0) << " " << pt->at(0) << " " << eta->at(0) << " " << phi->at(0) << " " << e->at(0) << std::endl;

  std::vector<PseudoJet> particles;
  for(unsigned int i = 0; i<id->size(); i++){
    if(TMath::Abs(eta->at(i))>2.6) continue;
    particles.push_back( PseudoJet(pt->at(i)*TMath::Cos(phi->at(i)),pt->at(i)*TMath::Sin(phi->at(i)),pt->at(i)/(TMath::Tan(2*TMath::ATan(TMath::Exp(-eta->at(i))))),e->at(i)));
  }

  JetDefinition jet_def(antikt_algorithm, 0.4);
  ClusterSequence cs(particles, jet_def);
  std::vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(20.0));
  std::cout << "Clustering with " << jet_def.description() << std::endl;

  std::cout <<   "        pt y phi" << std::endl;
  for (unsigned i = 0; i < jets.size(); i++) {
    std::cout << "jet " << i << ": "<< jets[i].pt() << " " 
                   << jets[i].rap() << " " << jets[i].phi() << std::endl;
    std::vector<PseudoJet> constituents = jets[i].constituents();
    for (unsigned j = 0; j < constituents.size(); j++) {
      std::cout << "    constituent " << j << "'s pt: " << constituents[j].pt()
           << std::endl;
    }
  }

  softDropGroomer sd = softDropGroomer(0.1,0,0.4);
  std::vector<PseudoJet> groomedJets = sd.doGrooming(jets);
  std::cout << jets.at(1).pt() << " " << sd.getDR12().at(1) << std::endl;   

  /* 
  if( groomedJets.at(1).has_pieces() ) {
        std::vector<PseudoJet> subjets = groomedJets.at(1).pieces();
        fastjet::PseudoJet subjet1 = subjets[0];
        fastjet::PseudoJet subjet2 = subjets[1];
        std::cout << subjet1.pt() << " " << subjet1.eta() << " " << subjet1.phi() << std::endl;
        std::cout << subjet2.pt() << " " << subjet2.eta() << " " << subjet2.phi() << std::endl;
  }*/

  /*
  std::vector<PseudoJet> jetConstituents = jets[0].constituents();
  JetDefinition CAjet_def(fastjet::cambridge_algorithm, 999);
  fastjet::ClusterSequence CAcs(jetConstituents, CAjet_def);
  std::vector<fastjet::PseudoJet> tempJets = fastjet::sorted_by_pt(CAcs.inclusive_jets());
  if(tempJets.size()<1) return;
  
  fastjet::contrib::SoftDrop * sd = new fastjet::contrib::SoftDrop(0.1, 0, 0.4 );
  fastjet::PseudoJet transformedJet = tempJets[0];
  if ( transformedJet == 0 ) return;
  transformedJet = (*sd)(transformedJet);
  std::vector<fastjet::PseudoJet> subjets;
  std::cout << "here" << std::endl;*/
  /*if ( transformedJet.has_pieces() ) {
        subjets = transformedJet.pieces();
        fastjet::PseudoJet subjet1 = subjets[0];
        fastjet::PseudoJet subjet2 = subjets[1];
        std::cout << subjet1.pt() << " " << subjet1.eta() << " " << subjet1.phi() << std::endl;
        std::cout << subjet2.pt() << " " << subjet2.eta() << " " << subjet2.phi() << std::endl;
  }*/
}

int main(){
  
  runFJ();
  return 1;
}
