#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <dirent.h>
#include <sys/stat.h>
#include <TROOT.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TF1.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TStyle.h>

#define NSAMPLE 150
#define NFILE 16

void displaycaen(TString filename){  //filenameはEASIROCで記録されたunixtime
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    std::string path = "/Users/naokiotani/NINJA/trackeranalysis/rayrawanalysis/cosmictestatKyoto/data/caen/"; //caenのデータの場所
    std::string ifname;
    ifname = path + std::string(filename) + "wave.root";
    TFile *f = TFile::Open(ifname.c_str());

    TTree *t = (TTree*) f->Get("t");

    std::string msg;

    // set x axis
    double x[NSAMPLE];
    for(int i=0; i<NSAMPLE; i++){
        x[i] = i;
    }

    // set branch address
    double waveform[NFILE][NSAMPLE];
    for(int i=0; i<NFILE; i++){
      t->SetBranchAddress(Form("waveform%d",i), waveform[i]);
    }

    int nEntry = t->GetEntries();

    for(int iEntry=0; iEntry<nEntry; iEntry++){
    //for(int iEntry=nEntry; iEntry>0; iEntry--){
        t->GetEntry(iEntry);
        //t->GetEntry(nEntry-iEntry-1);

	//cout << waveform[0][0] <<waveform[7][7]<<endl;
        double cons[NFILE]= {0.0};
        int counter[NFILE]= {30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30};
        double value[NFILE+1] = {0.0};
        //cout<<cons[0]<<cons[1]<<cons[2]<<endl;
	//cout << counter[3] <<endl;
        for(int i=0; i<NFILE; i++){
            for(int j=0; j<30; j++){
                cons[i] += waveform[i][j];
                if(fabs(waveform[i][j]-waveform[i][j+1])>5){
                    counter[i]=j+1;
                    break;
                }
            }

            cons[i]=cons[i]/(double)counter[i];
            //cout << "con[0]=" <<cons[0] << "counter[3]=" << counter[3] <<endl;
            //積分
            for(int j=35;j<75;j++){
                value[i]+=cons[i]-waveform[i][j];
            }
        }

        //if(value[21]<-20){
            std::cout << iEntry << std::endl;
            for(int i=0; i<NFILE; i++){
                std::cout << "ch"+std::to_string(i)+": " << value[i] << std::endl;
                //std::cout << cons[i] << std::endl;
                //std::cout << cons << std::endl;
            }
            //cout<<cons[21]<<endl;
            // for(int j=30;j<70;j++){
            //     //cout<<cons[21]-waveform[21][j]<<endl;
            // }
            TGraph *g[NFILE];
            for(int i=0; i<NFILE; i++){
                g[i] = new TGraph(NSAMPLE, x, waveform[i]);
                g[i]->SetTitle("CAEN;Sampling;ADC value");
                g[i]->SetMarkerColor(i);
                g[i]->SetMarkerStyle(8);
                //g[i]->Draw("AP");
            }
            g[10]->SetMarkerColor(21);

	    //cout << waveform[0][0] << " "  << waveform[7][7] << endl;

	    TMultiGraph *mg =new TMultiGraph();
        //mg->Add(g[0]);
	    for(int i=0;i<NFILE;i++){
            mg->Add(g[i]);
	    }
            mg->SetTitle(";Sampling;ADC value");
            //mg->Add(g[1]);
            mg->Draw("AP");

            c1->Update();


            TLegend *legend = new TLegend( 0.8, 0.5, 0.99, 0.99) ; //（）の中は位置の指定（左下の x , y 、右上の x , y ）
            //legend->SetTextSize(0.05);

	    legend->AddEntry( g[0], "ch0" , "p");
            legend->AddEntry( g[1], "ch1" , "p");
            legend->AddEntry( g[2], "ch2" , "p");
            legend->AddEntry( g[3], "ch3" , "p");
            legend->AddEntry( g[4], "ch4" , "p");
            legend->AddEntry( g[5], "ch5" , "p");
            legend->AddEntry( g[6], "ch6" , "p");
            legend->AddEntry( g[7], "ch7" , "p");
            legend->AddEntry( g[8], "ch8" , "p");
            legend->AddEntry( g[9], "ch9" , "p");
            legend->AddEntry( g[10], "ch10" , "p");
            legend->AddEntry( g[11], "ch11" , "p");
            legend->AddEntry( g[12], "ch12" , "p");
            legend->AddEntry( g[13], "ch13" , "p");
            legend->AddEntry( g[14], "ch14" , "p");
            legend->AddEntry( g[15], "ch15" , "p");


            legend->SetFillColor(0);
            legend->Draw() ;
            c1->Update();
            gSystem->ProcessEvents();

            std::cout << "Input q to quit." << std::endl;
            getline(std::cin, msg);
            if(msg == "q"){
                break;
            }
	    //}
    }
}
