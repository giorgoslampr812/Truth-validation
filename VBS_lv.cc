// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DirectFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/InvisibleFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/HeavyHadrons.hh"
//#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/TauFinder.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

#include "Rivet/Tools/Cutflow.hh"
#include "TFile.h"
#include <TTree.h>
#include <fstream>
#include <iostream>
#include "TH1.h"
#include "TParameter.h"

namespace Rivet {


  /// @brief Add a short analysis description here
  class VBS_lv : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(VBS_lv);


    /// @name Analysis methods
    /// @{
    bool is_w_jet(const Jet& j, const std::vector<Jet>& w_jets) {
    for (const auto& wj : w_jets) {
        // Compare constituent counts first for quick rejection
        if (j.particles().size() != wj.particles().size()) continue;
        
        // Full constituent comparison
        bool all_match = true;
        for (size_t i = 0; i < j.particles().size(); ++i) {
            if (j.particles()[i].genParticle() != wj.particles()[i].genParticle()) {
                all_match = false;
                break;
            }
          }
        if (all_match) return true;
        }
        return false;
      }
    
    bool hasAncestor(const std::shared_ptr<const HepMC3::GenParticle>& gp, int pdgid, int maxDepth = 10) {
  if (!gp || maxDepth <= 0) return false;

  auto vtx = gp->production_vertex();
  if (!vtx) return false;

    for (const auto& parent : vtx->particles_in()) {
      if (!parent) continue;
      if (parent->pdg_id() == pdgid) return true;
      if (hasAncestor(parent, pdgid, maxDepth - 1)) return true;
      }
    return false;
  }
    
  double jetWidth(const Jet& jet) {
    double sumPt = 0.0, sumPtDR = 0.0;
    const FourMomentum& jMom = jet.momentum();
    for (const Particle& p : jet.constituents()) {
      double dR = deltaR(jMom, p.momentum());
      double pt = p.momentum().pT();
      sumPt += pt;
      sumPtDR += pt * dR;
    }
  return (sumPt > 0) ? sumPtDR / sumPt : 0.0;
  }
  
    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs(Cuts::abseta < 5.0);
      declare(fs, "FS");
      
      PromptFinalState photons(Cuts::abspid == PID::PHOTON); 
      // Electrons and muons in Fiducial PS
    PromptFinalState bare_lep(Cuts::abspid == PID::ELECTRON ||  Cuts::abspid == PID::MUON, true);
    const Cut cuts_lep = Cuts::pT > 0 *GeV && Cuts::abseta < 5.0 ; 
    //LeptonFinder all_dressed_lep(bare_lep,0.1, cuts_lep );
    DressedLeptons all_dressed_lep(photons, bare_lep,0.1, cuts_lep );
    declare(all_dressed_lep, "DressedLeptons" ); 
    
    // Missing momentum
    declare(MissingMomentum(fs), "MET");

    VetoedFinalState ivfs(fs);
    ivfs.addVetoOnThisFinalState(VisibleFinalState(fs));
    declare(ivfs, "InvisibleFS");

    IdentifiedFinalState nu_id;
    nu_id.acceptNeutrinos();
    PromptFinalState neutrinos(nu_id);
    neutrinos.acceptTauDecays(false);
    declare(neutrinos, "Neutrinos");

       //Jet forming
       VetoedFinalState vfs(FinalState(Cuts::abseta < 5.0));
       vfs.addVetoOnThisFinalState(all_dressed_lep);
       vfs.addVetoOnThisFinalState(neutrinos);
 
       FastJets jets(vfs, FastJets::ANTIKT, 0.4);
       declare(jets, "Jets");
       FastJets jets_fat(vfs, FastJets::ANTIKT, 1.0);
       declare(jets_fat, "Jets_fat");
       

      // Book histograms
      book(cutflow_resolved, "cutflow_resolved", {0, 1, 2, 3, 4});//{"AllEvents" "MinLeptoncount" "Min Jet count" "VBS jet Pair"  "VBS jet pair mass selection" "resolved signal jets" "W boson selection" "same lepton flavour" "Z boson selection"});
      // Specify custom binnings
//--------tagging jet variables--------------------------------------------
      
      _file = new TFile("/data/glamprinoudis/RIVET/Run2_WZ_364850_100k.root", "RECREATE");
      _tree = new TTree("events", "Rivet ntuple");
      // --Your mom is a fat jet 
      _tree->Branch("jet_fat_pt", &_jet_fat_pt);
      _tree->Branch("jet_fat_m", &_jet_fat_m);
      _tree->Branch("jet_fat_eta", &_jet_fat_eta);
      _tree->Branch("jet_fat_phi", &_jet_fat_phi);
      _tree->Branch("jet_fat_E", &_jet_fat_E);
      _tree->Branch("jet_fat_y", &_jet_fat_y);
      _tree->Branch("Jet_fat_m", &_Jet_fat_m);
      _tree->Branch("Jet_fat_pt", &_Jet_fat_pt);
      _tree->Branch("Jet_fat_eta", &_Jet_fat_eta);
      _tree->Branch("Jet_fat_phi", &_Jet_fat_phi);
      _tree->Branch("Jet_fat_E", &_Jet_fat_E);
      _tree->Branch("Jet_fat_y", &_Jet_fat_y);
      
      //---Small R Jets   
      _tree->Branch("jet_pt", &_jet_pt); 
      _tree->Branch("jet_m", &_jet_m);
      _tree->Branch("jet_eta", &_jet_eta);
      _tree->Branch("jet_phi", &_jet_phi);
      _tree->Branch("jet_E", &_jet_E);
      _tree->Branch("jet_y", &_jet_y);
      _tree->Branch("Jet_n", &_Jet_n);
      _tree->Branch("Jet1_m", &_Jet1_m);
      _tree->Branch("Jet1_pt", &_Jet1_pt);
      _tree->Branch("Jet1_eta", &_Jet1_eta);
      _tree->Branch("Jet1_phi", &_Jet1_phi);
      _tree->Branch("Jet1_E", &_Jet1_E);
      _tree->Branch("Jet1_y", &_Jet1_y);

      _tree->Branch("Jet2_m", &_Jet2_m);
      _tree->Branch("Jet2_pt", &_Jet2_pt);
      _tree->Branch("Jet2_eta", &_Jet2_eta);
      _tree->Branch("Jet2_phi", &_Jet2_phi);
      _tree->Branch("Jet2_E", &_Jet2_E);
      _tree->Branch("Jet2_y", &_Jet2_y);
      
      _tree->Branch("Jet3_m", &_Jet3_m);
      _tree->Branch("Jet3_pt", &_Jet3_pt);
      _tree->Branch("Jet3_eta", &_Jet3_eta);
      _tree->Branch("Jet3_phi", &_Jet3_phi);
      _tree->Branch("Jet3_E", &_Jet3_E);
      _tree->Branch("Jet3_y", &_Jet3_y);

      _tree->Branch("Jet4_m", &_Jet4_m);
      _tree->Branch("Jet4_pt", &_Jet4_pt);
      _tree->Branch("Jet4_eta", &_Jet4_eta);
      _tree->Branch("Jet4_phi", &_Jet4_phi);
      _tree->Branch("Jet4_E", &_Jet4_E);
      _tree->Branch("Jet4_y", &_Jet4_y);

      //---leptons
    
      _tree->Branch("lep_m", &_lep_m);
      _tree->Branch("lep_pt", &_lep_pt);
      _tree->Branch("lep_eta", &_lep_eta);
      _tree->Branch("lep_phi", &_lep_phi);
      _tree->Branch("lep_E", &_lep_E);
      _tree->Branch("lep_y", &_lep_y);
      _tree->Branch("lep_pid", &_lep_pid);      
 
      _tree->Branch("Lepton1_m", &_Lepton1_m);
      _tree->Branch("Lepton1_pt", &_Lepton1_pt);
      _tree->Branch("Lepton1_eta", &_Lepton1_eta);
      _tree->Branch("Lepton1_phi", &_Lepton1_phi);
      _tree->Branch("Lepton1_E", &_Lepton1_E);
      _tree->Branch("Lepton1_y", &_Lepton1_y);
      _tree->Branch("Lepton1_pid", &_Lepton1_pid);
      /* 
      _tree->Branch("Lepton2_m", &_Lepton2_m);
      _tree->Branch("Lepton2_pt", &_Lepton2_pt);
      _tree->Branch("Lepton2_eta", &_Lepton2_eta);
      _tree->Branch("Lepton2_phi", &_Lepton2_phi);
      _tree->Branch("Lepton2_E", &_Lepton2_E);
      _tree->Branch("Lepton2_y", &_Lepton2_y);
      _tree->Branch("Lepton2_pid", &_Lepton2_pid);
      */
      //Missing energy
      _tree->Branch("Miss_ET", &_Miss_ET);
      _tree->Branch("Invisible_ET", &_Invisible_ET);
      //Neutrino variables
      _tree->Branch("Neutrino_ET", &_Neutrino_ET);
      _tree->Branch("Neutrino_pt", &_Neutrino_pt);
      _tree->Branch("Neutrino_eta", &_Neutrino_eta);
      _tree->Branch("Neutrino_phi", &_Neutrino_phi);
      //Btagging
      _tree->Branch("Btags", &_Btags); 
      //==event variables
      _tree->Branch("evt", &_evt, "evt/L");

      // generator bookkeeping (unweighted sample -> genWeight = 1)
      _tree->Branch("genWeight", &_genWeight, "genWeight/D");
      _tree->Branch("genWeights", &_genWeights);  
      // cut bookkeeping
      _tree->Branch("cutMask", &_cutMask, "cutMask/l"); // 64-bit mask
      
      // Cutflow histogram (raw generator-weighted counts)
      _h_cutflow = new TH1D("cutflow_physics", "Cutflow (physics gen weights)", 14, 0, 14);
      _h_cutflow->Sumw2();
      
      _h_cutflow_raw = new TH1D("cutflow_raw", "Cutflow (raw mc events)", 14, 0, 14);
      _h_cutflow_raw->Sumw2();
      
      
      // Init counters
      _nEvents = 0;
      _sumGenWeights = 0.0;
      _sumGenWeights2 = 0.0;

      // User-supplied normalization metadata (set these correctly)
      //_crossSection_pb = 0.018687;      // pb   <-- CHANGE to your sample's cross-section
      _kFactor = 1.0;                // optional k-factor
      _targetLumi_fb = 300;         // fb^-1 <-- target lumi you want to normalize to
    }
    /// Perform the per-event analysis
    void analyze(const Event& event) {
      
      //const GenEvent* gev = event.hepmcEventPtr(); // Rivet -> HepMC pointer
      const GenEvent* gev = event.genEvent();
      if (!gev) {
          MSG_WARNING("No HepMC GenEvent pointer for this Rivet Event!");
        return;
        }
     // const GenRunInfo* runinfo = gev->run_Info();
      // Selection cuts
      //std::cout << "Run Attributes" << gev->run_info()->attributes() << std::endl;
      
      _evt = gev->event_number();
      
        _genWeights = gev->weights();            // nominal is usually index 0
        if (!_genWeights.empty()) _genWeight = _genWeights[0];
        else {
          // fallback: try legacy single-weight accessor
            _genWeight = gev->weight();
            _genWeights.clear();
            _genWeights.push_back(_genWeight);
            }
      
      // increment global bookkeeping
    _nEvents += 1;
    _sumGenWeights += _genWeight;
    _sumGenWeights2 += _genWeight * _genWeight;
    
    _h_cutflow->Fill(1.0,1.0);
    

      
       //const vector<DressedLepton>& dressedleptons = apply<LeptonFinder>(event, "DressedLeptons").dressedLeptons();
       const vector<DressedLepton>& dressedleptons = apply<DressedLeptons>(event, "DressedLeptons").dressedLeptons();
       const FinalState &fs = apply<FinalState>(event, "FS");
       cutflow_resolved->fill(0); // All events
      
       
       if (dressedleptons.size() < 1) vetoEvent;
       
       cutflow_resolved->fill(1);
       _h_cutflow->Fill(2.0, _genWeight);
       
       Jets jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 0.0*GeV && Cuts::abseta < 4.5);
       Jets jets_fat = apply<FastJets>(event, "Jets_fat").jetsByPt(Cuts::pT > 0.0*GeV && Cuts::abseta < 4.5);
      const Vector3& vet = apply<MissingMomentum>(event, "MET").vectorEt();
      const double etmiss = vet.perp(); 
      const Particles& neutrinos = apply<PromptFinalState>(event, "Neutrinos").particlesByPt();
      const FinalState& ifs = apply<VetoedFinalState>(event, "InvisibleFS");
      const FourMomentum metinvisible = sum(ifs.particles(), FourMomentum());
      double Invisible_et = metinvisible.pt();
      if(jets.size() < 4) vetoEvent;
      cutflow_resolved->fill(2);
      _h_cutflow->Fill(3.0, _genWeight);
    
      sort(jets.begin(), jets.end(), [](const Jet& a, const Jet& b) {
         return a.pT() > b.pT();
      });    
      sort(jets_fat.begin(), jets_fat.end(), [](const Jet& a, const Jet& b) {
         return a.pT() > b.pT();
      });
      
    Particles lep;
    for (size_t i = 0; i < dressedleptons.size(); ++i) {
          if(dressedleptons[i].pt() > 0*GeV){ lep +=dressedleptons[i];
          }
    }
           
    sort(lep.begin(), lep.end(), [](const Particle& a, const Particle& b) { 
        return a.pT() > b.pT(); 
    });
    
    if (lep.empty()) vetoEvent;     
    if (lep.size() < 1) vetoEvent;    
    
    cutflow_resolved->fill(3);
    _h_cutflow->Fill(4.0, _genWeight);

    //Vector filling
    
    _lep_pt.clear();
    _lep_m.clear();
    _lep_phi.clear();
    _lep_eta.clear();
    _lep_E.clear();
    _lep_y.clear();
    _lep_pid.clear();

    for (const Particle& lepton : lep) {
     _lep_pt.push_back(lepton.pt());
     _lep_m.push_back(lepton.mass());
     _lep_phi.push_back(lepton.phi());
     _lep_eta.push_back(lepton.eta());
     _lep_E.push_back(lepton.E());
     _lep_y.push_back(lepton.rap());
     _lep_pid.push_back(lepton.pid());
    }

    _jet_pt.clear();
    _jet_m.clear();
    _jet_phi.clear();
    _jet_eta.clear();
    _jet_E.clear();
    _jet_y.clear();
    
    int btags=0;
    for (const Jet& j : jets) {
     _jet_pt.push_back(j.pt());
     _jet_m.push_back(j.mass());
     _jet_phi.push_back(j.phi());
     _jet_eta.push_back(j.eta());
     _jet_E.push_back(j.E());
     _jet_y.push_back(j.rap());
     if (j.bTagged()) btags=btags+1;
    }
    
    _jet_fat_pt.clear();
    _jet_fat_m.clear();
    _jet_fat_phi.clear();
    _jet_fat_eta.clear();
    _jet_fat_E.clear();
    _jet_fat_y.clear();

    for (const Jet& j2 : jets_fat) {
     _jet_fat_pt.push_back(j2.pt());
     _jet_fat_m.push_back(j2.mass());
     _jet_fat_phi.push_back(j2.phi());
     _jet_fat_eta.push_back(j2.eta());
     _jet_fat_E.push_back(j2.E());
     _jet_fat_y.push_back(j2.rap());
    } 
   
    //Vector filling
    _Jet_fat_m = jets_fat[0].mass();
    _Jet_fat_pt = jets_fat[0].pt();
    _Jet_fat_eta = jets_fat[0].eta();
    _Jet_fat_phi = jets_fat[0].phi();
    _Jet_fat_E = jets_fat[0].E();
    _Jet_fat_y = jets_fat[0].rap();
    _Jet_n = jets.size();
    _Miss_ET = etmiss;
    _Invisible_ET = Invisible_et;
    _Neutrino_ET = neutrinos[0].Et();
    _Neutrino_pt = neutrinos[0].pt();
    _Neutrino_phi = neutrinos[0].phi();
    _Neutrino_eta = neutrinos[0].eta();
    _Btags = btags;
    _Jet1_m = jets[0].mass();
    _Jet1_pt = jets[0].pt();
    _Jet1_eta = jets[0].eta();
    _Jet1_phi = jets[0].phi();
    _Jet1_E = jets[0].E();
    _Jet1_y = jets[0].rap();
    _Jet2_m = jets[1].mass();
    _Jet2_pt = jets[1].pt();
    _Jet2_eta = jets[1].eta();
    _Jet2_phi = jets[1].phi();
    _Jet2_E = jets[1].E();
    _Jet2_y = jets[1].rap();
    _Jet3_m = jets[2].mass();
    _Jet3_pt = jets[2].pt();
    _Jet3_eta = jets[2].eta();
    _Jet3_phi = jets[2].phi();
    _Jet3_E = jets[2].E();
    _Jet3_y = jets[2].rap();
    _Jet4_m = jets[3].mass();
    _Jet4_pt = jets[3].pt();
    _Jet4_eta = jets[3].eta();
    _Jet4_phi = jets[3].phi();
    _Jet4_E = jets[3].E();
    _Jet4_y = jets[3].rap();
    _Lepton1_m = lep[0].mass();
    _Lepton1_pt = lep[0].pt();
    _Lepton1_eta = lep[0].eta();
    _Lepton1_phi = lep[0].phi();
    _Lepton1_E = lep[0].E();
    _Lepton1_y = lep[0].rap();
    _Lepton1_pid = lep[0].pid();
    /*
    _Lepton2_m = lep[1].mass();
    _Lepton2_pt = lep[1].pt();
    _Lepton2_eta = lep[1].eta();
    _Lepton2_phi = lep[1].phi();
    _Lepton2_E = lep[1].E();
    _Lepton2_y = lep[1].rap();
    _Lepton2_pid = lep[1].pid();
    */
    _tree->Fill();
    
    
    }
    
    /// Normalise histograms etc., after the run
    void finalize() {
       std::cout << "FINALIZE CALLED" << std::endl;
      //for (auto& item : _h) {
      //  scale(item, (crossSection()/picobarn)/sumOfWeights());
        // merge overflow bin into the last bin
      //}
      

      //normalize(_h, crossSection()/femtobarn);
      normalize(_h, 1.0);

// Scale histograms using femtobarn unit
      //scale(_h, crossSection()/femtobarn/sumW());
      cout << "Events before selection" << ": " << cutflow_resolved->bin(0).numEntries() << " entries" << endl;
      cout << "Events 1st selection : min lepton " << ": " << cutflow_resolved->bin(1).numEntries() << " entries" << endl;
      cout << "Events 2nd selection : min jets" << ": " << cutflow_resolved->bin(2).numEntries() << " entries" << endl;
      cout << "Events 3rd selection : lep opp charge same flavour" << ": " << cutflow_resolved->bin(3).numEntries() << " entries" << endl;
      
    double L_pb = _targetLumi_fb * 1000.0;
    double sumWeights = _sumGenWeights; // for unweighted = _nEvents
    double sigma_eff_pb = (crossSection()/picobarn) * _kFactor;
    _crossSection_pb = crossSection()/picobarn;
    double norm = 0.0;
    if (sumWeights > 0.0) norm = sigma_eff_pb * L_pb / sumWeights;
    
     _h_cutflow_raw->Fill(1.0, cutflow_resolved->bin(0).numEntries());
     _h_cutflow_raw->Fill(2.0, cutflow_resolved->bin(1).numEntries());
     _h_cutflow_raw->Fill(3.0, cutflow_resolved->bin(2).numEntries());
     _h_cutflow_raw->Fill(4.0, cutflow_resolved->bin(3).numEntries());
    
       
      _file->cd();
      _tree->Write();
      // write raw cutflow and a scaled copy (expected yields at target lumi)
    _h_cutflow->Write(); // physics counts (sum of gen weights)
    _h_cutflow_raw->Write();//raw counts
    TH1D* h_cutflow_scaled = (TH1D*)_h_cutflow->Clone("cutflow_scaled");
    h_cutflow_scaled->Scale(norm);
    h_cutflow_scaled->SetTitle("Cutflow (expected physics events at target lumi)");
    h_cutflow_scaled->Write();

    // Store global metadata in an Info tree or as TParameters
    TTree* info = new TTree("Info", "Sample info");
    info->Branch("N_gen", &_nEvents, "N_gen/L");
    info->Branch("sumGenWeights", &_sumGenWeights, "sumGenWeights/D");
    info->Branch("sumGenWeights2", &_sumGenWeights2, "sumGenWeights2/D");
    info->Branch("crossSection_pb", &_crossSection_pb, "crossSection_pb/D");
    info->Branch("kFactor", &_kFactor, "kFactor/D");
    info->Branch("targetLumi_fb", &_targetLumi_fb, "targetLumi_fb/D");
    double norm_for_tree = norm;
    info->Branch("normFactor", &norm_for_tree, "normFactor/D");
    info->Fill();
    info->Write();

    // also write norm as a TParameter (handy)
    TParameter<double>("normFactor", norm).Write();
      
      
      _file->Close();
      

    }

    /// @}
    
    


    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    map<string, Histo1DPtr> _d;
    //map<string, BinnedHistoPtr<int>> _d;
    Histo1DPtr cutflow_resolved;
    /// @}
    fastjet::Filter _trimmer;
    const double MZ_PDG = 91.1876;
    const double MW_PDG = 83.385;
    const double GammaZ_PDG = 2.4952;
    const double GammaW_PDG = 2.085;
    
    private:
    TFile* _file;
    TTree* _tree;
    TH1D* _h_cutflow;
    TH1D* _h_cutflow_raw;
    std::vector<double> _lep_m, _lep_pt , _lep_eta, _lep_phi,  _lep_E, _lep_y, _lep_pid;
    std::vector<double> _jet_fat_pt, _jet_fat_m, _jet_fat_eta, _jet_fat_phi, _jet_fat_E, _jet_fat_y;
    std::vector<double> _jet_pt, _jet_m, _jet_eta, _jet_phi, _jet_E, _jet_y;
    unsigned long long _Jet_n;
    double _Jet_fat_m, _Jet_fat_pt, _Jet_fat_eta, _Jet_fat_phi, _Jet_fat_E, _Jet_fat_y;
    double _Jet2_m, _Jet2_pt, _Jet2_eta, _Jet2_phi, _Jet2_E, _Jet2_y, _Jet1_m, _Jet1_pt, _Jet1_eta, _Jet1_phi, _Jet1_E, _Jet1_y,_Jet3_m, _Jet3_pt, _Jet3_eta, _Jet3_phi, _Jet3_E, _Jet3_y,_Jet4_m, _Jet4_pt, _Jet4_eta, _Jet4_phi, _Jet4_E, _Jet4_y;
    double _Lepton1_m, _Lepton1_pt , _Lepton1_eta, _Lepton1_phi, _Lepton1_px, _Lepton1_py, _Lepton1_pz, _Lepton1_E, _Lepton1_y, _Lepton2_m, _Lepton2_pt , _Lepton2_eta, _Lepton2_phi, _Lepton2_px, _Lepton2_py, _Lepton2_pz, _Lepton2_E, _Lepton2_y, _Lepton1_pid, _Lepton2_pid;
    double _Miss_ET,_Btags,_Invisible_ET;
    double _Neutrino_ET, _Neutrino_pt, _Neutrino_eta, _Neutrino_phi;
   // event branches
  unsigned long long _run;
  unsigned long long _evt;
  int _lumi;
  double _genWeight;
  std::vector<double> _genWeights;
  unsigned long long _cutMask;
  // bookkeeping
  unsigned long long _nEvents;
  double _sumGenWeights;
  double _sumGenWeights2;

  // normalization metadata (user must set correct values)
  double _crossSection_pb;
  double _kFactor;
  double _targetLumi_fb;
  
  };


  RIVET_DECLARE_PLUGIN(VBS_lv);

}
