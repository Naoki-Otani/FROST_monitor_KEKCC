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
#define INTMIN 0
#define INTRANGE 60

void calibcaen(TString filename){  //filenameは日付-no eg)20240708-1

    std::string path_caen = "/group/nu/ninja/work/otani/FROST_cosmictest_dataandoutput/data/caen/"; //caenのデータの場所
    TFile *f = new TFile(TString(path_caen + "caencalib" + std::string(filename) + "wave.root"));
    TTree *t = (TTree*) f->Get("t");

    TCanvas *c1 = new TCanvas("c1", "c1", 2000, 1400);
    c1->Divide(4,4);

    std::ofstream fout(TString(path_caen) + "../../calibration/caen/calibresult/caencalib"+ filename + ".csv");
    fout << "channel,m0,m1,integralrange\n";

    // set branch address
    double waveform[NFILE][NSAMPLE];
    for(int i=0; i<NFILE; i++){
      t->SetBranchAddress(Form("waveform%d",i), waveform[i]);
    }

    int nEntry = t->GetEntries();

    TH1D *h[NFILE];
    for(int i=0;i<NFILE;i++){
      h[i] = new TH1D(Form("h[%d]",i),Form("h[%d]",i),150,-100,200);
    }

    double cons[NFILE] = {1228.3,1241.1,1232.3,1295.0,
        1236.8,1239.9,1250.0,1248.9,
        1228.6,1237.9,1237.7,1224.7,
        1215.9,1168.2,1228.1,1225.8};

    const int windowSize = 20;
    const int startMax = 70; //0,5,...,,70
    const int slideStep = 5;

    int integralrange[NFILE] = {INTRANGE};

    for(int i=0; i<NFILE; i++){ //default integral range
        integralrange[i] = INTRANGE;
    }

    for(int iEntry=0; iEntry<nEntry; iEntry++){
        t->GetEntry(iEntry);

        for(int i=0; i<NFILE; i++){

            std::vector<double> conssliding;

            for(int start = 0; start<=startMax; start+=slideStep){
                double sum = 0.0;

                for(int j=start; j<start+windowSize; j++){
                    sum += waveform[i][j];
                }
                double ave = sum / (double)windowSize;
                conssliding.push_back(ave);
            }
            std::sort(conssliding.begin(), conssliding.end());

            cons[i] = conssliding[conssliding.size()-1];

            double value = 0.0;
            for(int j=INTMIN;j<INTMIN+integralrange[i];j++){
                value += cons[i] - waveform[i][j];
            }
            h[i]->Fill(value);
        }
    }


//////////////////////fitting範囲の決定////////////////////////////////////////////////
    double fitmin0[NFILE];
    double fitmax0[NFILE];
    double fitmin1[NFILE];
    double fitmax1[NFILE];
    for(int i = 0; i < NFILE; i++) {
        fitmin0[i]=h[i]->GetBinCenter(h[i]->GetMaximumBin())-5;
        fitmax0[i]=h[i]->GetBinCenter(h[i]->GetMaximumBin())+5;

        // int nbins = h[i]->GetNbinsX();
        int nbins = h[i]->GetNbinsX();
        double maxval = -1;
        int binmax = -1;

        // x=15以上のbinのみ対象
        for(int b = 1; b <= nbins; b++) {
            double binCenter = h[i]->GetXaxis()->GetBinCenter(b);
            if(binCenter < h[i]->GetBinCenter(h[i]->GetMaximumBin())+15) continue;

            double content = h[i]->GetBinContent(b);
            if(content > maxval) {
                maxval = content;
                binmax = b;
            }
        }
        // x=20以上のbinのみ対象 1 p.e.のピークが小さくて上の方法だとmaxbinがx~15になってしまう場合
        if((h[i]->GetXaxis()->GetBinCenter(binmax) < h[i]->GetBinCenter(h[i]->GetMaximumBin())+18) && (h[i]->GetBinContent(binmax) > h[i]->GetBinContent(binmax+1))){
            maxval = -1;
            binmax = -1;
            for(int b = 1; b <= nbins; b++) {
                double binCenter = h[i]->GetXaxis()->GetBinCenter(b);
                if(binCenter < h[i]->GetBinCenter(h[i]->GetMaximumBin())+20) continue;

                double content = h[i]->GetBinContent(b);
                if(content > maxval) {
                    maxval = content;
                    binmax = b;
                }
            }
        }
        // x=30以上のbinのみ対象 1 p.e.のピークが小さくて上の方法だとmaxbinがx~20になってしまう場合
        if((h[i]->GetXaxis()->GetBinCenter(binmax) < h[i]->GetBinCenter(h[i]->GetMaximumBin())+23) && (h[i]->GetBinContent(binmax) > h[i]->GetBinContent(binmax+1))){
            maxval = -1;
            binmax = -1;
            for(int b = 1; b <= nbins; b++) {
                double binCenter = h[i]->GetXaxis()->GetBinCenter(b);
                if(binCenter < h[i]->GetBinCenter(h[i]->GetMaximumBin())+30) continue;

                double content = h[i]->GetBinContent(b);
                if(content > maxval) {
                    maxval = content;
                    binmax = b;
                }
            }
        }

        // binmaxが見つかったか確認
        if(binmax == -1){
            // 見つからない場合はデフォルト設定（例：無理やり範囲を決めるなど）
            fitmin1[i] = 40;
            fitmax1[i] = 60;
        } else {
            // double peak = h[i]->GetXaxis()->GetBinCenter(binmax);
            double peak = h[i]->GetXaxis()->GetBinCenter(binmax);
            fitmin1[i] = peak - 10;
            fitmax1[i] = peak + 10;
        }
    }


/////////////////////////////////////////////////////////////////////////////////
    TF1 *Gaus0 = new TF1("Gaus0", "gaus", -10, 10);
    TF1 *Gaus1 = new TF1("Gaus1", "gaus", 20, 50);
    std::string title[NFILE];

    // gStyle->SetOptFit(); //Fitting結果の統計boxを表示

    for (int i = 0; i < NFILE; i++) {
        title[i] = "CAEN ch " + std::to_string(i);
        c1->cd(i+1);
        h[i]->SetFillColor(5);
        h[i]->SetTitle(title[i].c_str());
        h[i]->GetXaxis()->SetTitle("ADC Integral");
        h[i]->GetYaxis()->SetTitle("Number of Events");
        h[i]->Draw();


        // === Fit 1 ===
        TF1* f0 = new TF1(Form("f0_%d", i), "gaus", fitmin0[i], fitmax0[i]);
        f0->SetLineColor(kRed);
        // コピーしたヒストにフィット（これにより stats 分離）
        TH1F* hcopy0 = (TH1F*)h[i]->Clone(Form("hgaus0[%d]", i));
        hcopy0->Fit(f0, "", "", fitmin0[i], fitmax0[i]);

        // === Fit 2 ===
        TF1* f1 = new TF1(Form("f1_%d", i), "gaus", fitmin1[i], fitmax1[i]);
        f1->SetLineColor(kBlue);
        TH1F* hcopy1 = (TH1F*)h[i]->Clone(Form("hgaus1[%d]", i));
        hcopy1->Fit(f1, "", "", fitmin1[i], fitmax1[i]);
        // 元ヒスト再描画＋フィット関数・統計ボックス
        h[i]->Draw();
        f0->Draw("same");
        f1->Draw("same");


        // ===== 赤フィットの値と誤差 =====
        double c0  = f0->GetParameter(0);
        double c0e = f0->GetParError(0);
        double m0  = f0->GetParameter(1);
        double m0e = f0->GetParError(1);
        double s0  = f0->GetParameter(2);
        double s0e = f0->GetParError(2);

        // ===== 青フィットの値と誤差 =====
        double c1  = f1->GetParameter(0);
        double c1e = f1->GetParError(0);
        double m1  = f1->GetParameter(1);
        double m1e = f1->GetParError(1);
        double s1  = f1->GetParameter(2);
        double s1e = f1->GetParError(2);

        // ===== TLatex描画 =====
        TLatex *latex = new TLatex();
        latex->SetNDC();
        latex->SetTextSize(0.04);

        // 赤パラメータ
        latex->SetTextColor(kRed);
        latex->DrawLatex(0.45, 0.85, Form("Const = %.2f #pm %.2f", c0, c0e));
        latex->DrawLatex(0.45, 0.80, Form("Mean  = %.2f #pm %.2f", m0, m0e));
        latex->DrawLatex(0.45, 0.75, Form("Sigma = %.2f #pm %.2f", s0, s0e));

        // 青パラメータ（Y位置を下げる）
        latex->SetTextColor(kBlue);
        latex->DrawLatex(0.45, 0.70, Form("Const = %.2f #pm %.2f", c1, c1e));
        latex->DrawLatex(0.45, 0.65, Form("Mean  = %.2f #pm %.2f", m1, m1e));
        latex->DrawLatex(0.45, 0.60, Form("Sigma = %.2f #pm %.2f", s1, s1e));

        fout << i << "," << m0 << "," << m1 << "," << integralrange[i] << "\n";
    }

    c1->Update();
    c1->SaveAs(TString(path_caen) + "../../calibration/caen/plot/caencalib"+ filename + ".pdf");
    fout.close();

}
