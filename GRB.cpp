#include "TH1D.h"
#include "TFile.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TApplication.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TProfile.h"
#include<iostream>
using namespace std;

TFile* PULL(double pull_max=8,TString filename="~/Desktop/Uni/Master/Astrophysics_Lab_I/Automne/scrate.root"){
    
//Useful constant for the Program-----------------------------------------
    const double bin_size=1.0;//search bin size seconds
    const int jump=5; //jump in bins before and after for doing mean
    const int MeanWidth=25; // nbr of bins to do the mean estimation
    const int TShift=100000; // Time window in bins
    const double Time_GRB = 1488809240;
    
//File opening/creation---------------------------------------------------
    TFile * input= new TFile(filename,"READ");
    TFile *output = new TFile("out.root","RECREATE");
    
//Creation of alias for TTree Branches------------------------------------
    TTree * tree=(TTree*)input->Get("f1rate");
    double unix_time;
    tree->SetBranchAddress("unix_time",&unix_time);
    double rate[14];
    tree->SetBranchAddress("rate", rate);
    double rate_err[14];
    tree->SetBranchAddress("rate_err", rate_err);
    float fe_cosmic;
    tree->SetBranchAddress("fe_cosmic", &fe_cosmic);
    
    
//Creation of the histogram-----------------------------------------------
    
    TProfile *h = new TProfile("h","hist with right error", TShift/bin_size,-TShift/2,TShift/2);
    TH1D *h1 = new TH1D("h1","hist of the mean", TShift/bin_size,-TShift/2,TShift/2);
    TH1D *h2 = new TH1D("h2","hist data-mean", TShift/bin_size,-TShift/2,TShift/2);
    TH1D *h3 = new TH1D("h3","hist pull", TShift/bin_size,-TShift/2,TShift/2);
    TH1D *e = new TH1D("e","hist of the rate_err^2", TShift/bin_size,-TShift/2,TShift/2);
    TH1D *fe = new TH1D("fe","hist of fe_cosmic", TShift/bin_size,-TShift/2,TShift/2);
    TH1D *p = new TH1D("p","hist of pull values", 100,-10,40);
    
    
    const Long64_t first_entry=33000000;
    const Long64_t last_entry=34000000;//tree->GetEntries()

//Filling h profile histogram and e/fe histogram--------------------------
    
    for (Long64_t i=first_entry; i<last_entry; i++){
        tree->GetEntry(i);
        h->Fill(unix_time-Time_GRB, rate[13]);
        e->Fill(unix_time-Time_GRB, rate_err[13]*rate_err[13]);
        fe->Fill(unix_time-Time_GRB, fe_cosmic);
    }
    
//Setting e/fe value if there is bins in h--------------------------------
    for (int i=0 ;i<h->GetNbinsX();i++){
        int n_entry=h->GetBinEntries(i+1);
        if (n_entry!=0){
            e->SetBinContent(i+1,sqrt(e->GetArray()[i+1])/n_entry);
            fe->SetBinContent(i+1,fe->GetArray()[i+1]/n_entry);
        }
    }

//Calculation of the bins expected value----------------------------------
    double* mean=new double[3];
    double* t=new double[3];
    double timecheck=0;
    double tot_pull=0;
    double a,b,slope,val,pull;
    int na,nb;
    bool in_a_grb=false;
    
    for (int j=jump+MeanWidth;j<h->GetNbinsX()-jump-MeanWidth;j++){
        if (h->GetBinEntries(j)!=0){
            a =0;
            b =0;
            na=0;
            nb=0;
            
            for(int i=0;i<MeanWidth;i++){
                if (fe->GetArray()[j-jump-i+1]<2500) continue;
                if (h->GetBinEntries(j+jump+i+1)!=0){
                    a+=h->GetBinContent(j+jump+i+1);
                    na++;
                }
            }
            for(int i=0;i<MeanWidth;i++){
                if (fe->GetArray()[j-jump-i+1]<2500) continue;
                if (h->GetBinEntries(j-jump-i+1)!=0){
                    b+=h->GetBinContent(j-jump-i+1);
                    nb++;
                }
            }
            if (na<10) continue;
            if (nb<10) continue;
            if (na!=0){
                mean[0]=a/na;
            }
            if (nb!=0){
                mean[2]=b/nb;
            }

            t[0]=h->GetBinLowEdge(j+jump+MeanWidth/2);
            t[1]=h->GetBinLowEdge(j);
            t[2]=h->GetBinLowEdge(j-jump-MeanWidth/2);
            slope=(mean[2]-mean[0])/(t[2]-t[0]);
            mean[1]=slope*(t[1]-t[0])+mean[0];
            val=(h->GetBinContent(j+1)-mean[1]);
        
            pull=val;
            if (e->GetBinContent(j+1)!=0){
                pull/=e->GetBinContent(j+1);
            }
            h1->Fill(t[1],mean[1]);
            h2->Fill(t[1],val);
            h3->Fill(t[1],pull);
            p->Fill(pull);
            if (pull>pull_max){
                if (!in_a_grb) {
                    timecheck=t[1];
                    in_a_grb=true;
                }
                tot_pull+=pull;
            }
            else{
                if (in_a_grb){
                    if (t[1]-timecheck!=0) {
                    cout<<"anomaly at unix time: "<< timecheck << " duration "<<t[1]-timecheck<<" tot pull "<<tot_pull<< endl;
                    }
                    tot_pull=0;
                    in_a_grb=false;
                }
            }
        }
        else{
            if (in_a_grb){
                if (t[1]-timecheck!=0) {
                cout<<"anomaly at unix time: "<< timecheck << " duration "<<t[1]-timecheck<<" tot pull "<<tot_pull<< endl;
                }
                tot_pull=0;
                in_a_grb=false;
                
            }
        }
    }
    delete[] mean;
    delete[] t;
    
    output->cd();
    h->Write();
    e->Write();
    p->Write();
    h1->Write();
    h2->Write();
    h3->Write();
    fe->Write();
    return output;
}


int main(){
    TApplication *a = new TApplication("a", 0, 0);
    TFile* f=PULL();
    
    
    TCanvas *c1 = new TCanvas("c1","Test with fe<2500 and nbr mean 10 ",200,10,700,1200);
    TPad *pad1 = new TPad("pad1","The pad with the function",0.1,0.75,0.9,1);
    TPad *pad2 = new TPad("pad2","The pad with the histogram",0.1,0.5,0.9,0.75);
    TPad *pad3 = new TPad("pad3","The pad with the histogram",0.1,0.25,0.9,0.5);
    TPad *pad4 = new TPad("pad4","The pad with the histogram",0.1,0.01,0.9,0.25);
    pad1->Draw();
    pad2->Draw();
    pad3->Draw();
    pad4->Draw();
    TH1D* h=(TH1D*)f->Get("h");
    TH1D* h1=(TH1D*)f->Get("h1");
    TH1D* h2=(TH1D*)f->Get("h2");
    TH1D* h3=(TH1D*)f->Get("h3");
    TH1D* p=(TH1D*)f->Get("p");
    TH1D* e=(TH1D*)f->Get("e");
    /*
    h->GetXaxis()->SetRangeUser(-(TShift/2-jump-MeanWidth-1), TShift/2-jump-MeanWidth);
    h1->GetXaxis()->SetRangeUser(-(TShift/2-jump-MeanWidth-1), TShift/2-jump-MeanWidth);
    h2->GetXaxis()->SetRangeUser(-(TShift/2-jump-MeanWidth-1), TShift/2-jump-MeanWidth);
    h3->GetXaxis()->SetRangeUser(-(TShift/2-jump-MeanWidth-1), TShift/2-jump-MeanWidth);
    h->GetYaxis()->SetRangeUser(-(TShift/2-jump-MeanWidth-1), TShift/2-jump-MeanWidth);
    h1->GetYaxis()->SetRangeUser(-(TShift/2-jump-MeanWidth-1), TShift/2-jump-MeanWidth);
    h2->GetYaxis()->SetRangeUser(-(TShift/2-jump-MeanWidth-1), TShift/2-jump-MeanWidth);
    h3->GetYaxis()->SetRangeUser(-(TShift/2-jump-MeanWidth-1), TShift/2-jump-MeanWidth);
     */
    h2->GetYaxis()->SetRangeUser(-1500, 4000);
    h3->GetYaxis()->SetRangeUser(-20, 40);
    pad1->cd();
    h->Draw();
    pad2->cd();
    h1->Draw("hist");
    pad3->cd();
    h2->Draw("hist");
    pad4->cd();
    h3->Draw("hist");
    c1->Update();
    c1->Write();
    f->Close();
    delete a;
return 0;
}
//tree->Draw("unix_time-1488809240>>h","rate[10]/4","",1000000,33000000);
/*
TH1D *h = new TH1D("h","a trial histogram", 100, -5, 5);
for (Int_t i = 0; i < 10000; i++) h->Fill(gRandom->Gaus(0, 1));
cout<<"Hello"<<endl;
 
 
 */


//Look for fi lower than 2000
