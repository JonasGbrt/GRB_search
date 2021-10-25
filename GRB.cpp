#include "TH1D.h"
#include "TFile.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TApplication.h"
#include "TDirectory.h"
#include "TTree.h"
#include<iostream>
using namespace std;

void PULL(){
    TApplication *a = new TApplication("a", 0, 0);
    double bin=0.25;//quarter of second
    int jump=5; //jump before and after the bean evaluated (in second)
    int MeanWidth=25; // nbr of second used to make the mean for the estimation
    int TShift=1300; // Time shift for the histogram
    double timecheck=0; //used to check if it's one single event
    TFile * input= new TFile("~/Desktop/Uni/Master/Astrophysics_Lab_I/Automne/scrate.root","READ");
    TTree * tree=(TTree*)input->Get("f1rate");
    double unix_time;
    tree->SetBranchAddress("unix_time",&unix_time);
    double rate[14];
    tree->SetBranchAddress("rate", rate);
    double rate_err[14];
    tree->SetBranchAddress("rate_err", rate_err);
    TFile *output = new TFile("out.root","RECREATE");
    
    TCanvas *c1 = new TCanvas("c1","Histogram Drawing Options",200,10,700,1200);
    TPad *pad1 = new TPad("pad1","The pad with the function",0.1,0.75,0.9,1);
    TPad *pad2 = new TPad("pad2","The pad with the histogram",0.1,0.5,0.9,0.75);
    TPad *pad3 = new TPad("pad3","The pad with the histogram",0.1,0.25,0.9,0.5);
    TPad *pad4 = new TPad("pad4","The pad with the histogram",0.1,0.01,0.9,0.25);
    pad1->Draw();
    pad2->Draw();
    pad3->Draw();
    pad4->Draw();
    
    TH1D *h = new TH1D("h","hist with right error", TShift,-TShift/2,TShift/2);
    TH1D *h1 = new TH1D("h1","hist of the mean", TShift,-TShift/2,TShift/2);
    TH1D *h2 = new TH1D("h2","hist data-mean", TShift,-TShift/2,TShift/2);
    TH1D *h3 = new TH1D("h3","hist pull", TShift,-TShift/2,TShift/2);
    TH1D *e = new TH1D("e","hist of the rate_err^2", TShift,-TShift/2,TShift/2);
    TH1D *p = new TH1D("p","hist of pull values", 100,-10,40);
    for (int i=0; i< tree->GetEntries(); i++){
        if (i < 33000000) continue;
        if (i > 33000000 + 1000000) continue;
        tree->GetEntry(i);
        h->Fill(unix_time-1488809240, rate[13]*bin);
        e->Fill(unix_time-1488809240, rate_err[13]*rate_err[13]);
    }
    double* er=new double[e->GetNbinsX()+2];
    for (int i=0 ;i<TShift+2;i++) er[i]=sqrt(e->GetArray()[i])*bin;
    h->SetError(er);
    


    for (int j=jump+MeanWidth;j<TShift+2-MeanWidth-jump;j++){
        double a =0;
        double b =0;
        for(int i=0;i<MeanWidth;i++) a+=h->GetArray()[j+jump+i];
        for(int i=0;i<MeanWidth;i++) b+=h->GetArray()[j-jump-i];
        double* mean=new double[3];
        double* t=new double[3];
        mean[0]=a/MeanWidth;
        mean[2]=b/MeanWidth;
        t[0]=j+jump+MeanWidth/2;
        t[1]=j;
        t[2]=j-jump-MeanWidth/2;
        double slope=(mean[2]-mean[0])/(t[2]-t[0]);
        mean[1]=slope*(t[1]-t[0])+mean[0];
        double pull=(h->GetArray()[j]-mean[1])/er[j];
        h1->Fill(j-TShift/2+1,mean[1]);
        p->Fill(pull);
        delete[] mean;
        delete[] t;
        if (pull>8){
            if (timecheck==0) {
                timecheck=j;
                cout<<"anomaly at unix time: "<< j<< " with a pull > 8"<< endl;
            }
            if (timecheck-j>20){
                cout<<"anomaly at unix time: "<< j<< " with a pull > 8"<< endl;
                timecheck=j;
            }
                
            
            
        }
        
    }
    for (int i=0;i<h->GetNbinsX();i++){
        h2->Fill(i-(TShift/2+1),h->GetArray()[i]-h1->GetArray()[i]);
        }
    
    for (int i=0;i<h->GetNbinsX();i++){
        h3->Fill(i-(TShift/2+1),h2->GetArray()[i]/er[i]);
        }
    h->GetXaxis()->SetRangeUser(-(TShift/2-jump-MeanWidth-1), TShift/2-jump-MeanWidth);
    h1->GetXaxis()->SetRangeUser(-(TShift/2-jump-MeanWidth-1), TShift/2-jump-MeanWidth);
    h2->GetXaxis()->SetRangeUser(-(TShift/2-jump-MeanWidth-1), TShift/2-jump-MeanWidth);
    h3->GetXaxis()->SetRangeUser(-(TShift/2-jump-MeanWidth-1), TShift/2-jump-MeanWidth);
    pad1->cd();
    h->Draw();
    pad2->cd();
    h1->Draw("hist");
    pad3->cd();
    h2->Draw("hist");
    pad4->cd();
    h3->Draw("hist");
    
    output->cd();
    h->Write();
    e->Write();
    p->Write();
    h1->Write();
    c1->Write();
    output->Close();
    delete a ;
}


int main(){

    PULL();

return 0;
}
//tree->Draw("unix_time-1488809240>>h","rate[10]/4","",1000000,33000000);
/*
TH1D *h = new TH1D("h","a trial histogram", 100, -5, 5);
for (Int_t i = 0; i < 10000; i++) h->Fill(gRandom->Gaus(0, 1));
cout<<"Hello"<<endl;
 */
