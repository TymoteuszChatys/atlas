// Analysis.py
// Skeleton code in python provided for you
// In place of this comment you should write [your name] -- [the date] and update it as you go!
// Make sure to make backups and comment as you go along :)

// Header guard to ensure file is imported properly
#ifndef Analysis
#define Analysis

// Include the file that lets the program know about the data
#include "backend/CLoop.h"

#define TOTAL_WEIGHT_E 19630128.89
#define TOTAL_WEIGHT_MU 19631161.45
double selected_weight_sum_e{}, selected_weight_sum_mu{};
size_t number_of_selected_e{}, number_of_selected_mu;
double total_weight_sum_e{}, total_weight_sum_mu{};
size_t total_number_of_e{}, total_number_of_mu{};

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
    
    //PT cone and ET cone histograms
    h_ptcone30 = new TH1F("pt30e","PTCone30",50,0,2000);
    h_etcone20 = new TH1F("et20e","ETCone20",400,-7e3,7e3);
    
    //Invariant mass histograms
    h_invariant_mass = new TH1F("invariant_mass","Invariant mass of the system",400,0e3,400e3);
    //Now our estimate of background events we plot hisograms of the data we discarded from the above histograms
    invariant_mass_background = new TH1F("invariant_mass_background_ee","Invariant mass of the system",400,0e3,200e3);
    //Monte Carlo simulations only valid for 60MeV+, therefore define new invariant mass graphs between 60GeV-120MeV
    h_invariant_mass_100 = new TH1F("invariant_mass_60","Invariant mass of the system",50,100e3,150e3);

}

void CLoop::Fill(double weight) {
    // This function is where you select events and fill your histograms
    // It is called once PER EVENT.
    bool background{false};
    // For example, filling a histogram with the number of leptons in each event
    h_lep_n -> Fill(lep_n, weight);
    if (lep_type->size() != 4) return;
    
    int cs = 0, fs = 0;
    for (int i=0;i<4;i++) {
      cs += lep_charge->at(i);
      fs += lep_charge->at(i)*lep_type->at(i);
    }
    if ((cs != 0)  || (fs != 0)) {
        background = true;
    }
    
    double invariant_mass = 0.0;
    for(int i=0;i<4;i++) {
      for(int j=i+1;j<4;j++) {
	invariant_mass += lep_pt->at(i)*lep_pt->at(j)*(cosh(lep_eta->at(i)-lep_eta->at(j))-cos(lep_phi->at(i)-lep_phi->at(j)));}
    }
    invariant_mass = sqrt(2.0 * invariant_mass);
    h_invariant_mass -> Fill(invariant_mass,weight);

    if(invariant_mass < 150e3 && invariant_mass > 100e3){
        h_invariant_mass_100 ->Fill(invariant_mass,weight);
    }


    if (background == false) {
        /*if (lep_type->at(i_pos)==11) {
            h_invariant_mass_ee -> Fill(invariant_mass,weight);
            h_invariant_mass_ee_60 -> Fill(invariant_mass,weight);
            h_ptcone30_e -> Fill(lep_ptcone30->at(i_pos),weight);
            h_ptcone30_e -> Fill(lep_ptcone30->at(i_neg),weight);
            h_etcone20_e -> Fill(lep_etcone20->at(i_pos),weight);
            h_etcone20_e -> Fill(lep_etcone20->at(i_neg),weight);
                //number_of_selected_e++;
	            //selected_weight_sum_e += weight;
        } else if (lep_type->at(i_pos)==13) {
            h_invariant_mass_mumu -> Fill(invariant_mass,weight);
            h_invariant_mass_mumu_60 -> Fill(invariant_mass,weight);
            h_ptcone30_m -> Fill(lep_ptcone30->at(i_pos),weight);
            h_ptcone30_m -> Fill(lep_ptcone30->at(i_neg),weight);
            h_etcone20_m -> Fill(lep_etcone20->at(i_pos),weight);
            h_etcone20_m -> Fill(lep_etcone20->at(i_neg),weight);
	            //number_of_selected_mu++;
                //selected_weight_sum_mu += weight;
        }*/
        
        //if (lep_ptcone30->at(i_pos)>100 || lep_ptcone30->at(i_neg)>100) {return;}
        
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

    }/*else if(background == true){
        if (lep_type->at(i_pos)==11) {
            invariant_mass_background_ee -> Fill(invariant_mass,weight);
        } else if (lep_type->at(i_pos)==13) {
            invariant_mass_background_mumu -> Fill(invariant_mass,weight);
        }
    }*/
    
}

void CLoop::Style(int colourcode) {
    // This function is where you can control the style elements of your histograms and write them to a file
    // It is called once per data set

    // For example, set some properties of the lep_n histogram
    h_lep_n->GetXaxis()->SetTitle("Number of leptons per event"); // label x axis
    h_lep_n->GetYaxis()->SetTitle("Number of entries"); // label y axis
    h_lep_n->SetLineColor(colourcode); // set the line colour to red
    // For more information see https://root.cern.ch/root/htmldoc/guides/users-guide/Histograms.html
    h_lep_type->GetXaxis()->SetTitle("Type"); // label x axis
    h_lep_type->GetYaxis()->SetTitle("Number of entries"); // label y axis
    h_lep_type->SetLineColor(colourcode); // set the line colour to green
    
    h_lep_pt->GetXaxis()->SetTitle("Transverse momentum of electrons"); // label x axis
    h_lep_pt->GetYaxis()->SetTitle("Number of entries"); // label y axis
    h_lep_pt->SetLineColor(colourcode); // set the line colour to red
    // Write histograms to a file
    // This needs to be done for each histogram
    h_lep_ptm->GetXaxis()->SetTitle("Transverse momentum of muons"); // label x axis
    h_lep_ptm->GetYaxis()->SetTitle("Number of entries"); // label y axis
    h_lep_ptm->SetLineColor(colourcode); // set the line colour to red
    
    h_lep_pt_negative->GetXaxis()->SetTitle("Transverse momentum of negative leptons"); // label x axis
    h_lep_pt_negative->GetYaxis()->SetTitle("Number of entries"); // label y axis
    h_lep_pt_negative->SetLineColor(colourcode); // set the line colour to red
    
    h_lep_pt_positive->GetXaxis()->SetTitle("Transverse momentum of positive leptons"); // label x axis
    h_lep_pt_positive->GetYaxis()->SetTitle("Number of entries"); // label y axis
    h_lep_pt_positive->SetLineColor(colourcode); // set the line colour to red
   
    h_lep_eta->GetXaxis()->SetTitle("Pseudorapidity of the leptons"); // label x axis
    h_lep_eta->GetYaxis()->SetTitle("Number of entries"); // label y axis
    h_lep_eta->SetLineColor(colourcode); // set the line colour to red
    
    h_lep_phi->GetXaxis()->SetTitle("Azimuthal angle of the leptons"); // label x axis
    h_lep_phi->GetYaxis()->SetTitle("Number of entries"); // label y axis
    h_lep_phi->SetLineColor(colourcode); // set the line colour to red 
    
    //Invariant mass plots
    //Selected data
        //Monte Carlo valid for >60GeV therefore plot only valid data from histogram
    h_invariant_mass->GetXaxis()->SetTitle("Invariant Mass of System"); // label x axis
    h_invariant_mass->GetYaxis()->SetTitle("Number of entries"); // label y axis
    h_invariant_mass->SetLineColor(colourcode); // set the line colour to red
    
    h_invariant_mass_100->GetXaxis()->SetTitle("Invariant Mass of System"); // label x axis
    h_invariant_mass_100->GetYaxis()->SetTitle("Number of entries"); // label y axis
    h_invariant_mass_100->SetLineColor(colourcode); // set the line colour to red
    //Background (rejected) data
    invariant_mass_background->GetXaxis()->SetTitle("Invariant Mass of System"); // label x axis
    invariant_mass_background->GetYaxis()->SetTitle("Number of entries"); // label y axis
    invariant_mass_background->SetLineColor(colourcode); // set the line colour to red

    h_ptcone30->GetXaxis()->SetTitle("30-cone Momentum Sum"); // label x axis
    h_ptcone30->GetYaxis()->SetTitle("Number of entries"); // label y axis
    h_ptcone30->SetLineColor(colourcode); // set the line colour to red
    
    h_etcone20->GetXaxis()->SetTitle("20-cone Energy Sum"); // label x axis
    h_etcone20->GetYaxis()->SetTitle("Number of entries"); // label y axis
    h_etcone20->SetLineColor(colourcode); // set the line colour to red
    
    h_lep_n->Write();

    h_lep_type->Write();
    h_lep_pt->Write();
    h_lep_ptm->Write();
    h_lep_pt_positive->Write();
    h_lep_pt_negative->Write();
    h_lep_eta -> Write();
    h_lep_phi -> Write();
    h_ptcone30 -> Write();
    h_etcone20 -> Write();

    //invariant masses
    h_invariant_mass -> Write();
    h_invariant_mass_100 -> Write();
    invariant_mass_background -> Write();
    //Cross section calculations
    double INTEGRATED_LUMINOSITY = 10.064;
    double efficiency_e = selected_weight_sum_e/total_weight_sum_e;
    double efficiency_mu = selected_weight_sum_mu/total_weight_sum_mu;
    cout << "Efficiency - electrons: " << efficiency_e << endl;
    cout << "Efficiency - muons: " << efficiency_mu << endl << endl;
    double cross_section_e = ((2*number_of_selected_e)-total_number_of_e)/(INTEGRATED_LUMINOSITY*efficiency_e);
    double cross_section_mu = ((2*number_of_selected_mu)-total_number_of_mu)/(INTEGRATED_LUMINOSITY*efficiency_mu);
    cout << "Cross section - electrons: " << cross_section_e << endl;
    cout << "Cross section - muons: " << cross_section_mu << endl;

   
}
#endif // End header guard

