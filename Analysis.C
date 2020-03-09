// Analysis.py
// Skeleton code in python provided for you
// In place of this comment you should write [your name] -- [the date] and update it as you go!
// Make sure to make backups and comment as you go along :)

// Header guard to ensure file is imported properly
#ifndef Analysis
#define Analysis

// Include the file that lets the program know about the data
#include "backend/CLoop.h"

#define COLOUR kRed

void CLoop::Book() {
    // This function is where you "book" your histograms
    // It is called once per data set
    // You will need to do this for each histogram you want to plot
    // The syntax is: histogram = new TH1F("histogram_name", "Title", number_of_bins, x_min, x_max); 
    // For example, booking a histogram to plot number of leptons per event:
    h_lep_n = new TH1F("lep_n","Number of leptons",10,-0.5,9.5);
    // A couple of other basic histograms
    h_lep_type = new TH1F("lep_type","Type of leptons",17,0.5,17.5);
    h_lep_pt = new TH1F("lep_pt","Transverse momentum of electrons",200,0,14e4);
    h_lep_ptm = new TH1F("lep_ptm","Transverse momentum of muons",200,0,14e4);
    h_lep_pt_positive = new TH1F("lep_pt_positive","Transverse momentum of positive leptons",200,0,14e4); 
    h_lep_pt_negative = new TH1F("lep_pt_negative","Transverse momentum of negative leptons",200,0,14e4);
    h_lep_eta = new TH1F("lep_eta","Pseudorapidity of the leptons",200,-3.14,3.14);
    h_lep_phi = new TH1F("lep_phi","Azimuthal angle of the leptons",200,-3.14,3.14);
    h_invariant_mass = new TH1F("invariant_mass","Invariant mass of the system (ee)",400,0,200e3);
    h_invariant_mass2 = new TH1F("invariant_mass2","Invariant mass of the system (mumu)",400,0,200e3);
    h_ptcone30_e = new TH1F("pt30e","PTCone30 (electrons)",50,0,2000);
    h_ptcone30_m = new TH1F("pt30m","PTCone30 (muons)",50,0,2000);
    h_etcone20_e= new TH1F("et20e","ETCone20 (electrons)",400,-7e3,7e3);
    h_etcone20_m= new TH1F("et20m","ETCone20 (muons)",400,-7e3,7e3);
    
}

void CLoop::Fill(double weight) {
    // This function is where you select events and fill your histograms
    // It is called once PER EVENT.
    
    // To fill a histogram you should write:
    // histogram->Fill(quantity,weight)

    // For example, filling a histogram with the number of leptons in each event
    h_lep_n -> Fill(lep_n, weight);
    
    int i_pos=-1, i_neg=-1;
    if (lep_type->size() > 2){
        float pt_pos=0.0, pt_neg=0.0;
        for(int i=0;i<lep_type->size();i++) {
            if(lep_charge->at(i)>0 && lep_pt->at(i)>pt_pos) {
            pt_pos=lep_pt->at(i);
            i_pos=i;
            } 
            else if (lep_charge->at(i)<0 && lep_pt->at(i)>pt_neg) {
            pt_neg=lep_pt->at(i);
            i_neg=i;
            }
        }
        if (i_pos<0 || i_neg<0 || lep_type->at(i_pos)!=lep_type->at(i_neg)) {
            return;
        }
    }else if (lep_type->size() == 2){ 
        if(lep_charge->at(0)*lep_charge->at(1)>0){
        return;
        }
	i_pos=0;
	i_neg=1;
    }
    
    //cout << "got" << i_pos << " " << i_neg << endl;
    
    float invariant_mass = sqrt(2*(lep_pt->at(i_pos)*lep_pt->at(i_neg))*(cosh((lep_eta->at(i_pos))-lep_eta->at(i_neg))-cos(lep_phi->at(i_pos)-lep_phi->at(i_neg))));
    if (lep_type->at(i_pos)==11) {
      h_invariant_mass -> Fill(invariant_mass,weight);
      h_ptcone30_e -> Fill(lep_ptcone30->at(i_pos),weight);
      h_ptcone30_e -> Fill(lep_ptcone30->at(i_neg),weight);
      h_etcone20_e -> Fill(lep_etcone20->at(i_pos),weight);
      h_etcone20_e -> Fill(lep_etcone20->at(i_neg),weight);
    } else if (lep_type->at(i_pos)==13) {
      h_invariant_mass2 -> Fill(invariant_mass,weight);
      h_ptcone30_m -> Fill(lep_ptcone30->at(i_pos),weight);
      h_ptcone30_m -> Fill(lep_ptcone30->at(i_neg),weight);
      h_etcone20_m -> Fill(lep_etcone20->at(i_pos),weight);
      h_etcone20_m -> Fill(lep_etcone20->at(i_neg),weight);
    }
    
    if (lep_ptcone30->at(i_pos)>100 || lep_ptcone30->at(i_neg)>100) {return;}
    
    // loop over leptons in the event
    for (size_t ilep=0; ilep<lep_type->size(); ilep++) {
	h_lep_eta -> Fill(lep_eta->at(ilep),weight);
	h_lep_phi -> Fill(lep_phi->at(ilep),weight);
        // This is where you will want to fill histograms with properties of individual leptons
        // You will need to use quantity->at(ilep) to refer to properties of the lepton

        // For example, filling a histogram with the type of each lepton
        // muons have type = 13, electrons have type = 11

        h_lep_type -> Fill(lep_type->at(ilep),weight);
        // Filling the transverse momentum of each lepton, but only if it is an electron
	
        if (lep_type->at(ilep) == 11) {
            h_lep_pt -> Fill(lep_pt->at(ilep),weight);
        }else if(lep_type->at(ilep) == 13) {
	    h_lep_ptm -> Fill(lep_pt->at(ilep),weight);
	}
	
	if (lep_charge->at(ilep) > 0){
	    h_lep_pt_negative -> Fill(lep_pt->at(ilep),weight);
	}else if (lep_charge->at(ilep) <0){
	    h_lep_pt_positive -> Fill(lep_pt->at(ilep),weight);
	}
     // End of loop over leptons
    }
}

void CLoop::Style() {
    // This function is where you can control the style elements of your histograms and write them to a file
    // It is called once per data set

    // For example, set some properties of the lep_n histogram
    h_lep_n->GetXaxis()->SetTitle("Number of leptons per event"); // label x axis
    h_lep_n->GetYaxis()->SetTitle("Number of entries"); // label y axis
    h_lep_n->SetLineColor(COLOUR); // set the line colour to red
    // For more information see https://root.cern.ch/root/htmldoc/guides/users-guide/Histograms.html
    h_lep_type->GetXaxis()->SetTitle("Type"); // label x axis
    h_lep_type->GetYaxis()->SetTitle("Number of entries"); // label y axis
    h_lep_type->SetLineColor(COLOUR); // set the line colour to green
    
    h_lep_pt->GetXaxis()->SetTitle("Transverse momentum of electrons"); // label x axis
    h_lep_pt->GetYaxis()->SetTitle("Number of entries"); // label y axis
    h_lep_pt->SetLineColor(COLOUR); // set the line colour to red
    // Write histograms to a file
    // This needs to be done for each histogram
    h_lep_ptm->GetXaxis()->SetTitle("Transverse momentum of muons"); // label x axis
    h_lep_ptm->GetYaxis()->SetTitle("Number of entries"); // label y axis
    h_lep_ptm->SetLineColor(COLOUR); // set the line colour to red
    
    h_lep_pt_negative->GetXaxis()->SetTitle("Transverse momentum of negative leptons"); // label x axis
    h_lep_pt_negative->GetYaxis()->SetTitle("Number of entries"); // label y axis
    h_lep_pt_negative->SetLineColor(COLOUR); // set the line colour to red
    
    h_lep_pt_positive->GetXaxis()->SetTitle("Transverse momentum of positive leptons"); // label x axis
    h_lep_pt_positive->GetYaxis()->SetTitle("Number of entries"); // label y axis
    h_lep_pt_positive->SetLineColor(COLOUR); // set the line colour to red
   
    
    
    h_lep_eta->GetXaxis()->SetTitle("Pseudorapidity of the leptons"); // label x axis
    h_lep_eta->GetYaxis()->SetTitle("Number of entries"); // label y axis
    h_lep_eta->SetLineColor(COLOUR); // set the line colour to red
    
    h_lep_phi->GetXaxis()->SetTitle("Azimuthal angle of the leptons"); // label x axis
    h_lep_phi->GetYaxis()->SetTitle("Number of entries"); // label y axis
    h_lep_phi->SetLineColor(COLOUR); // set the line colour to red 
    
    h_invariant_mass->GetXaxis()->SetTitle("Invariant Mass of System"); // label x axis
    h_invariant_mass->GetYaxis()->SetTitle("Number of entries"); // label y axis
    h_invariant_mass->SetLineColor(COLOUR); // set the line colour to red
    
    h_invariant_mass2->GetXaxis()->SetTitle("Invariant Mass of System"); // label x axis
    h_invariant_mass2->GetYaxis()->SetTitle("Number of entries"); // label y axis
    h_invariant_mass2->SetLineColor(COLOUR); // set the line colour to red
    
    h_ptcone30_e->GetXaxis()->SetTitle("30-cone Momentum Sum"); // label x axis
    h_ptcone30_e->GetYaxis()->SetTitle("Number of entries"); // label y axis
    h_ptcone30_e->SetLineColor(COLOUR); // set the line colour to red
    
    h_ptcone30_m->GetXaxis()->SetTitle("30-cone Momentum Sum"); // label x axis
    h_ptcone30_m->GetYaxis()->SetTitle("Number of entries"); // label y axis
    h_ptcone30_m->SetLineColor(COLOUR); // set the line colour to red
    
    h_etcone20_e->GetXaxis()->SetTitle("20-cone Energy Sum"); // label x axis
    h_etcone20_e->GetYaxis()->SetTitle("Number of entries"); // label y axis
    h_etcone20_e->SetLineColor(COLOUR); // set the line colour to red
    
    h_etcone20_m->GetXaxis()->SetTitle("20-cone Energy Sum"); // label x axis
    h_etcone20_m->GetYaxis()->SetTitle("Number of entries"); // label y axis
    h_etcone20_m->SetLineColor(COLOUR); // set the line colour to red
    
    h_lep_n->Write();

    h_lep_type->Write();
    h_lep_pt->Write();
    h_lep_ptm->Write();
    h_lep_pt_positive->Write();
    h_lep_pt_negative->Write();
    h_lep_eta -> Write();
    h_lep_phi -> Write();
    h_invariant_mass -> Write();
    h_invariant_mass2 -> Write();
    h_ptcone30_e -> Write();
    h_ptcone30_m -> Write();
    h_etcone20_e -> Write();
    h_etcone20_m -> Write();

    

    
}
#endif // End header guard
