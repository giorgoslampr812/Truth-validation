#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

void analyse_1lep() {

    TFile *f = TFile::Open("Run2_WW_364848bv_100k.root");
    TTree *t = (TTree*)f->Get("events");

    // Leptons
    std::vector<double> *lep_pt=0,*lep_eta=0,*lep_phi=0,*lep_m=0;
    std::vector<int> *lep_pid=0;

    t->SetBranchAddress("lep_pt",&lep_pt);
    t->SetBranchAddress("lep_eta",&lep_eta);
    t->SetBranchAddress("lep_phi",&lep_phi);
    t->SetBranchAddress("lep_m",&lep_m);
    t->SetBranchAddress("lep_pid",&lep_pid);

    // Small-R jets
    std::vector<double> *jet_pt=0,*jet_eta=0,*jet_phi=0,*jet_m=0;

    t->SetBranchAddress("jet_pt",&jet_pt);
    t->SetBranchAddress("jet_eta",&jet_eta);
    t->SetBranchAddress("jet_phi",&jet_phi);
    t->SetBranchAddress("jet_m",&jet_m);

    // Fat jets
    std::vector<double> *jet_fat_pt=0,*jet_fat_eta=0,*jet_fat_phi=0,*jet_fat_m=0;

    t->SetBranchAddress("jet_fat_pt",&jet_fat_pt);
    t->SetBranchAddress("jet_fat_eta",&jet_fat_eta);
    t->SetBranchAddress("jet_fat_phi",&jet_fat_phi);
    t->SetBranchAddress("jet_fat_m",&jet_fat_m);

    //Missing Et
    double ETmiss;
    double Neutrino_pt;
    t->SetBranchAddress("Miss_ET",&ETmiss);
    t->SetBranchAddress("Neutrino_pt",&Neutrino_pt);
    double Btags;
    t->SetBranchAddress("Btags",&Btags);
    Long64_t nentries = t->GetEntries();

    std::cout<<"Total events: "<<nentries<<std::endl;
    
    // Generator weight
    std::vector<float> *genWeights=0;
    t->SetBranchAddress("genWeights",&genWeights);
    double Gen_sum[146];

    // Cutflow histogram
    TH1F *cutflow = new TH1F("cutflowr","Cutflowr",6,0,6);
    cutflow->GetXaxis()->SetBinLabel(1,"All events");
    cutflow->GetXaxis()->SetBinLabel(2,"Lepton fiducial");
    cutflow->GetXaxis()->SetBinLabel(3,"4 jets pass");
    cutflow->GetXaxis()->SetBinLabel(4,"Tag jets pass");
    cutflow->GetXaxis()->SetBinLabel(5,"W boson pass");
    cutflow->GetXaxis()->SetBinLabel(6,"mjjj > 220");

    TH1F *cutflow2 = new TH1F("cutflowm","Cutflowm",6,0,6);
    cutflow2->GetXaxis()->SetBinLabel(1,"All events");
    cutflow2->GetXaxis()->SetBinLabel(2,"Lepton fiducial");
    cutflow2->GetXaxis()->SetBinLabel(3,"4 jets pass");
    cutflow2->GetXaxis()->SetBinLabel(4,"Tag jets pass");
    cutflow2->GetXaxis()->SetBinLabel(5,"W boson pass_fat");

    for(Long64_t i=0;i<nentries;i++){

        t->GetEntry(i);

        //--------------------------------
        // Lepton fiducial cut
        //--------------------------------
        cutflow->Fill(0); // all events
        cutflow2->Fill(0);
        for(size_t k=0;k<genWeights->size();k++){
        Gen_sum[k] = Gen_sum[k] + genWeights->at(k);
        }	
        int good_leptons=0;
	bool flag=false;
        if (lep_pt->size()<1) continue;
     //   if((lep_pt->at(0)>27) && (fabs(lep_eta->at(0))<2.5) && (Neutrino_pt > 80)){ 
        for(size_t j=0;j<jet_pt->size();j++){
	 double eta1=jet_eta->at(j);
         double phi1=jet_phi->at(j);
         double deta4 = eta1 - lep_eta->at(0);
         double dphi4 = phi1 - lep_phi->at(0);
         double DR4 = sqrt(deta4*deta4 + dphi4*dphi4);
	 if((lep_pt->at(0)>27) && (fabs(lep_eta->at(0))<2.5) && (Neutrino_pt > 80) && (Neutrino_pt <200) && (DR4>0.4)) flag = true;
	}
	if (flag){
        cutflow->Fill(1); // lepton
        cutflow2->Fill(1);
	}else {
        continue;
	}
        //--------------------------------
        // Jet fiducial selection
        //--------------------------------
        
        std::vector<int> good_jets;

        for(size_t j=0;j<jet_pt->size();j++){

            double pt=jet_pt->at(j);
	    double eta = abs(jet_eta->at(j));
            double eta1=jet_eta->at(j);
	    double phi1=jet_phi->at(j);
	    double deta1 = eta1 - lep_eta->at(0);
            double dphi1 = phi1 - lep_phi->at(0);
            double DR1 = sqrt(deta1*deta1 + dphi1*dphi1); 
            if(eta<2.5 && pt>20 && DR1>0.2) good_jets.push_back(j); //Jet overlap may help
            if(eta>2.5 && eta<4.5 && pt>30 && DR1>0.2) good_jets.push_back(j); // Jet overlap may help
	    //if(eta<2.5 && pt>20) good_jets.push_back(j); //Jet overlap may help
            //if(eta>2.5 && eta<4.5 && pt>30) good_jets.push_back(j); // Jet overlap may help
        }

        if(good_jets.size()<4) continue;
        cutflow->Fill(2); // four jets
        cutflow2->Fill(2);
        //--------------------------------
        // Sort jets by pT
        //--------------------------------
        
        std::sort(good_jets.begin(),good_jets.end(),
        [&](int a,int b){return jet_pt->at(a)>jet_pt->at(b);});

        int j1=good_jets[0];
        int j2=good_jets[1];
        int j3=good_jets[2];
        int j4=good_jets[3];

        //--------------------------------
        // Build jet 4-vectors
        //--------------------------------

        TLorentzVector J[4];
        int idx[4]={j1,j2,j3,j4};

        for(int k=0;k<4;k++){

            J[k].SetPtEtaPhiM(
                jet_pt->at(idx[k]),
                jet_eta->at(idx[k]),
                jet_phi->at(idx[k]),
                jet_m->at(idx[k])
            );
        }

        //--------------------------------
        // Find highest dijet mass
        //--------------------------------

        double max_mjj=-1;
        int tag1=-1,tag2=-1;

        for(int a=0;a<4;a++){
        for(int b=a+1;b<4;b++){

            double mjj=(J[a]+J[b]).M();
            if(J[a].Eta()*J[b].Eta()<0){
             if(mjj>max_mjj){

                max_mjj=mjj;
                tag1=a;
                tag2=b;
            }
	  }
         }
	}

        TLorentzVector jet_tag1=J[tag1];
        TLorentzVector jet_tag2=J[tag2];
        
        //--------------------------------
        // VBS tagging jet cuts
        //--------------------------------

        if(max_mjj<400) continue;
        if(jet_tag1.Pt()<30) continue;
        if(jet_tag2.Pt()<30) continue;
        cutflow->Fill(3); // jet tag
        cutflow2->Fill(3); 
        
        //--------------------------------
        // Remaining jets = W jets
        //--------------------------------

        std::vector<TLorentzVector> wjets;

        for(int k=0;k<4;k++){

            if(k==tag1 || k==tag2) continue;

            wjets.push_back(J[k]);
        }
   
        //--------------------------------
        // Sort W jets by pt
        //--------------------------------

        std::sort(wjets.begin(),wjets.end(),
        [](TLorentzVector a, TLorentzVector b){return a.Pt()>b.Pt();});

        TLorentzVector wj1=wjets[0];
        TLorentzVector wj2=wjets[1];
  
        double fat_pt=-1;
        double fat_eta=0;
        double fat_m=0;
	double fat_phi = 0;

        for(size_t j=0;j<jet_fat_pt->size();j++){

            double pt=jet_fat_pt->at(j);
            double eta=fabs(jet_fat_eta->at(j));
	    double ma = jet_fat_m->at(j);
	    
            double eta2 = jet_fat_eta->at(j);
            double phi2 = jet_fat_phi->at(j);
            double deta1 = eta2 - jet_tag1.Eta();
            double dphi1 = phi2 - jet_tag1.Phi();
            double DR1 = sqrt(deta1*deta1 + dphi1*dphi1);
            double deta2 = eta2 - jet_tag2.Eta();
            double dphi2 = phi2 - jet_tag2.Phi();
            double DR2 = sqrt(deta2*deta2 + dphi2*dphi2);
	    double deta3 = eta2 - lep_eta->at(0);
            double dphi3 = phi2 - lep_phi->at(0);
            double DR3 = sqrt(deta3*deta3 + dphi3*dphi3);
            if((pt>200) && (eta<2.0) && (DR3>1.0)){
                if(pt>fat_pt){

                    fat_pt=pt;
                    fat_eta=jet_fat_eta->at(j);
                    fat_m=jet_fat_m->at(j);
		    fat_phi = jet_fat_phi->at(j);
                }
            }
        }
	double deta1 = fat_eta - jet_tag1.Eta();
        double dphi1 = fat_phi - jet_tag1.Phi();
        double DR1 = sqrt(deta1*deta1 + dphi1*dphi1);
        double deta2 = fat_eta - jet_tag2.Eta();
        double dphi2 = fat_phi - jet_tag2.Phi();
        double DR2 = sqrt(deta2*deta2 + dphi2*dphi2);
	double DR_av = (DR1+DR2)/2;
	double deta3 = fat_eta - lep_eta->at(0);
        double dphi3 = fat_phi - lep_phi->at(0);
        double DR3 = sqrt(deta3*deta3 + dphi3*dphi3);
	if ((fat_m > 64) && (fat_m < 106) && (Btags==0) && (DR1>1.4) && (DR2>1.4)){
	  cutflow2->Fill(4); // Fat jet mass 
	  continue;
	}else {
		
        //--------------------------------
        // W jet cuts
        //--------------------------------

        if(wj1.Pt() < 40) continue;
        if(wj2.Pt() < 20) continue;

        double mwjj = (wj1 + wj2).M();

        if(mwjj < 64 || mwjj > 106) continue;
        cutflow->Fill(4); // W boson cut
        // compute all 3-jet invariant masses
        double m123 = (J[0] + J[1] + J[2]).M();
        double m124 = (J[0] + J[1] + J[3]).M();
        double m134 = (J[0] + J[2] + J[3]).M();
        double m234 = (J[1] + J[2] + J[3]).M();
        
        // apply mjjj cut
	if(Btags>0) continue;
        if(m123 <= 220) continue;
        if(m124 <= 220) continue;
        if(m134 <= 220) continue;
        if(m234 <= 220) continue;
        cutflow->Fill(5); // Mjjj cut
	
      }
    }

    for(size_t k=0;k<146;k++){
            std::cout << Gen_sum[k] <<",";
    }

    TCanvas *c2 = new TCanvas("c2","Cutflowr",800,600);
    cutflow->Draw("hist text");

    TCanvas *c1 = new TCanvas("c1","Cutflowm",800,600);
    cutflow2->Draw("hist text");  
}
