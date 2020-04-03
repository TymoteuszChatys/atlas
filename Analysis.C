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

int pt,et,mass,trans,phi,eta{};

double mass_Z{91.1876};

// CUTS APPLIED
double pt_cone_lower{2000};
double et_cone_higher{5000};
double mass_deviation{30};
double transverse_momentum_lower{30};
double phi_range{3.1415926};
double eta_range{2.8};


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
    h_ptcone30_e = new TH1F("pt30e","PTCone30 (electrons)",50,0,2000);
    h_ptcone30_m = new TH1F("pt30m","PTCone30 (muons)",50,0,2000);
    h_etcone20_e= new TH1F("et20e","ETCone20 (electrons)",400,-7e3,7e3);
    h_etcone20_m= new TH1F("et20m","ETCone20 (muons)",400,-7e3,7e3);
    
    //Invariant mass histograms
    h_invariant_mass_ee = new TH1F("invariant_mass","Invariant mass of the system (ee)",400,0,200e3);
    h_invariant_mass_mumu = new TH1F("invariant_mass_mumu","Invariant mass of the system (mumu)",400,0,200e3);
   
    //Now our estimate of background events we plot hisograms of the data we discarded from the above histograms
    invariant_mass_background_ee = new TH1F("invariant_mass_background_ee","Invariant mass of the system (ee)",400,0e3,200e3);
    invariant_mass_background_mumu = new TH1F("invariant_mass_background_mumu","Invariant mass of the system (mumu)",400,0e3,200e3);
    
    //Monte Carlo simulations only valid for 60MeV+, therefore define new invariant mass graphs between 60MeV-120MeV
    h_invariant_mass_ee_60 = new TH1F("invariant_mass_60","Invariant mass of the system (ee)",400,60e3,120e3);
    h_invariant_mass_mumu_60 = new TH1F("invariant_mass_mumu_60","Invariant mass of the system (mumu)",400,60e3,120e3);
    
}

void CLoop::Fill(double weight) {
    // This function is where you select events and fill your histograms
    // It is called once PER EVENT.
    
    // To fill a histogram you should write:
    // histogram->Fill(quantity,weight)
    bool background{false};
    bool rejected{false};
    // For example, filling a histogram with the number of leptons in each event
    h_lep_n -> Fill(lep_n, weight);
    
    int i_pos=-1, i_neg=-1;
    if (lep_type->size() > 2){
        float pt_pos=0.0, pt_neg=0.0;
        for(int i=0;i<lep_type->size();i++) {
            if(lep_charge->at(i)>0 && lep_pt->at(i)>pt_pos) {
                pt_pos=lep_pt->at(i);
                i_pos=i;
            }else if(lep_charge->at(i)<0 && lep_pt->at(i)>pt_neg) {
                pt_neg=lep_pt->at(i);
                i_neg=i;
            }
        }
        if (i_pos<0 || i_neg<0 || lep_type->at(i_pos)!=lep_type->at(i_neg)) {
            //If there is no positive charge leptons in the event add the previous or next lepton 
            //transverse momentums to the background and vice versa
            background = true;
            if(i_pos == -1){
                if(i_neg == 0){
                    //Fake i_pos, actually a negative charge
                    i_pos = 1;
                }else{
                    i_pos = i_neg-1;
                }                
            }else if(i_neg == -1){
                if(i_pos == 0){
                    i_neg = 1;
                }else{
                    i_pos = i_neg-1;
                }
            }
            return;
        }
    }else if (lep_type->size() == 2){ 
        if(lep_charge->at(0)*lep_charge->at(1)>0){
            background = true;
        }
	i_pos=0;
	i_neg=1;
    }

    double invariant_mass = sqrt(2*(lep_pt->at(i_pos)*lep_pt->at(i_neg))*(cosh((lep_eta->at(i_pos))-lep_eta->at(i_neg))-cos(lep_phi->at(i_pos)-lep_phi->at(i_neg))));
        //ptcone
        if (lep_ptcone30->at(i_pos)+lep_ptcone30->at(i_neg)>pt_cone_lower){
            rejected = true;
            pt++;
        }else{
        }

        //etcone
        if (lep_etcone20->at(i_pos)>et_cone_higher || lep_etcone20->at(i_neg)>et_cone_higher) {
            rejected = true;
            et++;
        }else{
        }

        //mass
        if (invariant_mass/1000 < mass_Z+mass_deviation && invariant_mass/1000 > mass_Z-mass_deviation) {
        }else{
            rejected = true;
            mass++; 
        }

        //transverse momentum
        if (lep_pt->at(i_pos)/1000<transverse_momentum_lower || lep_pt->at(i_neg)/1000<transverse_momentum_lower) {
            rejected = true; 
            trans++;
        }else{
        }

        //phi
        if (lep_phi->at(i_pos)<phi_range && lep_phi->at(i_neg)<phi_range && lep_phi->at(i_pos)>(phi_range*-1) && lep_phi->at(i_neg)>(phi_range*-1)) {
        }else{
            rejected = true;
            phi++;
        }

        //eta
        if (lep_eta->at(i_pos)<eta_range && lep_eta->at(i_neg)<eta_range && lep_eta->at(i_pos)>(eta_range*-1) && lep_eta->at(i_neg)>(eta_range*-1)) {
        }else{
            rejected = true;
            eta++;
        }

    if (rejected == false){
        if (lep_type->at(i_pos)==11 && lep_type->at(i_neg)==11) {
        total_number_of_e++;
        total_weight_sum_e+=weight;
        } else if (lep_type->at(i_pos)==13 && lep_type->at(i_neg)==13) {
        total_number_of_mu++;
        total_weight_sum_mu+=weight;
    }

    }
    
    //cout << "got" << i_pos << " " << i_neg << endl;
    if (rejected == false){
        if (background == false) {
            if (lep_type->at(i_pos)==11 && lep_type->at(i_neg)==11) {
                h_invariant_mass_ee -> Fill(invariant_mass,weight);
                h_invariant_mass_ee_60 -> Fill(invariant_mass,weight);
                h_ptcone30_e -> Fill(lep_ptcone30->at(i_pos),weight);
                h_ptcone30_e -> Fill(lep_ptcone30->at(i_neg),weight);
                h_etcone20_e -> Fill(lep_etcone20->at(i_pos),weight);
                h_etcone20_e -> Fill(lep_etcone20->at(i_neg),weight);
                    number_of_selected_e++;
                    selected_weight_sum_e += weight;
            } else if (lep_type->at(i_pos)==13 && lep_type->at(i_neg)==13) {
                h_invariant_mass_mumu -> Fill(invariant_mass,weight);
                h_invariant_mass_mumu_60 -> Fill(invariant_mass,weight);
                h_ptcone30_m -> Fill(lep_ptcone30->at(i_pos),weight);
                h_ptcone30_m -> Fill(lep_ptcone30->at(i_neg),weight);
                h_etcone20_m -> Fill(lep_etcone20->at(i_pos),weight);
                h_etcone20_m -> Fill(lep_etcone20->at(i_neg),weight);
                    number_of_selected_mu++;
                    selected_weight_sum_mu += weight;
            }
            
        
            
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

        }else if(background == true){
            if (lep_type->at(i_pos)==11) {
                invariant_mass_background_ee -> Fill(invariant_mass,weight);
            } else if (lep_type->at(i_pos)==13) {
                invariant_mass_background_mumu -> Fill(invariant_mass,weight);
            }
        }
    }
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
    h_invariant_mass_ee->GetXaxis()->SetTitle("Invariant Mass of System"); // label x axis
    h_invariant_mass_ee->GetYaxis()->SetTitle("Number of entries"); // label y axis
    h_invariant_mass_ee->SetLineColor(colourcode); // set the line colour to red
    
    h_invariant_mass_mumu->GetXaxis()->SetTitle("Invariant Mass of System"); // label x axis
    h_invariant_mass_mumu->GetYaxis()->SetTitle("Number of entries"); // label y axis
    h_invariant_mass_mumu->SetLineColor(colourcode); // set the line colour to red
    //Monte Carlo valid for >60GeV therefore plot only valid data from histogram
    h_invariant_mass_ee_60->GetXaxis()->SetTitle("Invariant Mass of System"); // label x axis
    h_invariant_mass_ee_60->GetYaxis()->SetTitle("Number of entries"); // label y axis
    h_invariant_mass_ee_60->SetLineColor(colourcode); // set the line colour to red
    
    h_invariant_mass_mumu_60->GetXaxis()->SetTitle("Invariant Mass of System"); // label x axis
    h_invariant_mass_mumu_60->GetYaxis()->SetTitle("Number of entries"); // label y axis
    h_invariant_mass_mumu_60->SetLineColor(colourcode); // set the line colour to red
    //Background (rejected) data
    invariant_mass_background_ee->GetXaxis()->SetTitle("Invariant Mass of System"); // label x axis
    invariant_mass_background_ee->GetYaxis()->SetTitle("Number of entries"); // label y axis
    invariant_mass_background_ee->SetLineColor(colourcode); // set the line colour to red

    invariant_mass_background_mumu->GetXaxis()->SetTitle("Invariant Mass of System"); // label x axis
    invariant_mass_background_mumu->GetYaxis()->SetTitle("Number of entries"); // label y axis
    invariant_mass_background_mumu->SetLineColor(colourcode); // set the line colour to red

    h_ptcone30_e->GetXaxis()->SetTitle("30-cone Momentum Sum"); // label x axis
    h_ptcone30_e->GetYaxis()->SetTitle("Number of entries"); // label y axis
    h_ptcone30_e->SetLineColor(colourcode); // set the line colour to red
    
    h_ptcone30_m->GetXaxis()->SetTitle("30-cone Momentum Sum"); // label x axis
    h_ptcone30_m->GetYaxis()->SetTitle("Number of entries"); // label y axis
    h_ptcone30_m->SetLineColor(colourcode); // set the line colour to red
    
    h_etcone20_e->GetXaxis()->SetTitle("20-cone Energy Sum"); // label x axis
    h_etcone20_e->GetYaxis()->SetTitle("Number of entries"); // label y axis
    h_etcone20_e->SetLineColor(colourcode); // set the line colour to red
    
    h_etcone20_m->GetXaxis()->SetTitle("20-cone Energy Sum"); // label x axis
    h_etcone20_m->GetYaxis()->SetTitle("Number of entries"); // label y axis
    h_etcone20_m->SetLineColor(colourcode); // set the line colour to red
    
    h_lep_n->Write();
    h_lep_type->Write();
    h_lep_pt->Write();
    h_lep_ptm->Write();
    h_lep_pt_positive->Write();
    h_lep_pt_negative->Write();
    h_lep_eta -> Write();
    h_lep_phi -> Write();
    h_invariant_mass_ee -> Write();
    h_invariant_mass_mumu -> Write();
    h_ptcone30_e -> Write();
    h_ptcone30_m -> Write();
    h_etcone20_e -> Write();
    h_etcone20_m -> Write();
    h_invariant_mass_ee_60 -> Write();
    h_invariant_mass_mumu_60 -> Write();
    invariant_mass_background_ee -> Write();
    invariant_mass_background_mumu -> Write();


    //Cross section calculations
    double INTEGRATED_LUMINOSITY = 10.064;
    double efficiency_e = selected_weight_sum_e/TOTAL_WEIGHT_E;
    double efficiency_mu = selected_weight_sum_mu/TOTAL_WEIGHT_MU;
    
    double cross_section_e = (((selected_weight_sum_e)-(total_weight_sum_e - selected_weight_sum_e))/(INTEGRATED_LUMINOSITY*efficiency_e))/1000;
    double cross_section_mu = (((selected_weight_sum_mu)-(total_weight_sum_mu - selected_weight_sum_mu))/(INTEGRATED_LUMINOSITY*efficiency_mu))/1000;
    

    
    std::ofstream outfile;
    auto end = std::chrono::system_clock::now();
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    std::cout << "Date: " << std::ctime(&end_time) << std::endl;
    std::cout << "parameters: " 
        << "ptcone: " << pt_cone_lower << " rejected: " << pt  << " "
        << "etcone: " << et_cone_higher << " rejected: " << et  << " "
        << "invariant mass: " << mass_deviation << " rejected: " << mass  << " "
        << "transverse mom: " << transverse_momentum_lower << " rejected: " << trans  << " "
        << "phi: " << phi_range << " rejected: " << phi  << " "
        << "eta: " << eta_range << " rejected: " << eta  << std::endl;

    std::cout << "Efficiency - electrons: " << efficiency_e << std::endl;
    std::cout << "Efficiency - muons: " << efficiency_mu << std::endl << std::endl;

    std::cout << "Cross section - electrons: " << cross_section_e << std::endl;
    std::cout << "Cross section - muons: " << cross_section_mu << std::endl << std::endl;

    std::cout << selected_weight_sum_e << ':' << total_weight_sum_e - selected_weight_sum_e << std::endl;
    std::cout << selected_weight_sum_mu << ':' << total_weight_sum_mu - selected_weight_sum_mu << std::endl << std::endl;

    outfile.open("log.txt", std::ios_base::app); // append instead of overwrite
    outfile << "Date: " << std::ctime(&end_time) << std::endl;
    outfile << "parameters: "  << std::endl
        << "ptcone: " << pt_cone_lower << " rejected: " << pt  << std::endl
        << "etcone: " << et_cone_higher << " rejected: " << et  << std::endl
        << "invariant mass: " << mass_deviation << " rejected: " << mass  << std::endl
        << "transverse mom: " << transverse_momentum_lower << " rejected: " << trans  << std::endl
        << "phi: " << phi_range << " rejected: " << phi  << std::endl
        << "eta: " << eta_range << " rejected: " << eta  << std::endl << std::endl; 

    outfile << "Efficiency - electrons: " << efficiency_e 
            << "  Efficiency - muons: " << efficiency_mu << std::endl << std::endl;
        
    outfile << "Cross section - electrons: " << cross_section_e
            << "  Cross section - muons: " << cross_section_mu << std::endl << std::endl;

    outfile << selected_weight_sum_e << ':' << total_weight_sum_e - selected_weight_sum_e << std::endl;
    outfile << selected_weight_sum_mu << ':' << total_weight_sum_mu - selected_weight_sum_mu << std::endl << std::endl;

    std::ofstream tempecsv;
    std::ofstream tempmucsv;
    tempecsv.open("tempe.csv", std::ios_base::app);
    tempmucsv.open("tempmu.csv", std::ios_base::app);
    tempecsv << selected_weight_sum_e << "," << total_weight_sum_e - selected_weight_sum_e << "\n";
    tempmucsv << selected_weight_sum_mu << "," << total_weight_sum_mu - selected_weight_sum_mu << "\n";
    tempecsv.close();
    tempmucsv.close();
}
#endif // End header guard
