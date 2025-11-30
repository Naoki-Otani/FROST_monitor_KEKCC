#define IntegralMin 1
#define IntegralMax 60
#define PED 510

#define IntegralMin 1
#define IntegralMax 60
#define PED 510

#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TSystem.h>
#include <TPad.h>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <map>

// -------------- グローバルch配列とタイトル情報を受けて描画する本体 --------------
void waveform_impl(TString filename,
                   const std::vector<int>& ch_list,
                   const std::vector<int>& rayraw_list,
                   const std::vector<int>& local_list,
                   int ev_start=0)
{
    // キャンバスを8分割
    TCanvas *c1 = new TCanvas("c1", "c1", 1200, 800);
    c1->Divide(4,2); // 最大8分割固定

    // ファイル・ツリーを開く
    TFile *f = TFile::Open(filename);
    if(!f || f->IsZombie()){ std::cerr << "Failed to open file: " << filename << std::endl; return; }
    TTree *t = (TTree*) f->Get("tree");
    if(!t){ std::cerr << "TTree 'tree' not found.\n"; return; }

    // waveformブランチ（vector<vector<double>>）
    std::vector<std::vector<double>> *waveform = nullptr;
    t->SetBranchAddress("waveform", &waveform);

    // サンプル数取得
    t->GetEntry(0);
    if(!waveform || waveform->empty()){ std::cerr << "Empty waveform at entry 0.\n"; return; }
    int NSAMPLE = (int)(waveform->at(0).size());

    // x軸
    std::vector<double> x(NSAMPLE);
    for (int i = 0; i < NSAMPLE; ++i) x[i] = i;

    int nEntry = t->GetEntries();
    std::string msg;

    // 実際に描画するch配列（最大8）
    std::vector<int> chs;
    std::vector<int> rrs;
    std::vector<int> lcs;
    for(size_t i=0; i<ch_list.size() && i<8; ++i){
        if(ch_list[i] >= 0){
            chs.push_back(ch_list[i]);
            rrs.push_back(rayraw_list[i]);
            lcs.push_back(local_list[i]);
        }
    }

    if(chs.empty()){
        std::cerr << "No valid channels specified.\n";
        return;
    }

    // イベントループ
    for (int iEntry = ev_start; iEntry < nEntry; ++iEntry) {
        std::cout << "==== Event " << iEntry << " ====" << std::endl;
        t->GetEntry(iEntry);

        // 各padをクリア
        for(int ipad=1; ipad<=8; ++ipad){
            c1->cd(ipad);
            gPad->Clear();
            gPad->SetGridx(); gPad->SetGridy();
        }

        // chごとに描画
        for(size_t idx=0; idx<chs.size() && idx<8; ++idx){
            c1->cd((int)idx+1);

            int ch_idx = chs[idx];
            int rr_id = rrs[idx];
            int lc_id = lcs[idx];

            if ((int)waveform->size() <= ch_idx || (int)(*waveform)[ch_idx].size() < NSAMPLE) {
                std::cerr << "Invalid waveform size at entry " << iEntry << " for ch " << ch_idx << std::endl;
                continue;
            }

            // ベースライン計算（立ち上がり検知）
            double cons = 0.0;
            int counter = 10;
            if(NSAMPLE <= 12) { counter = std::min(NSAMPLE, 10); }
            for (int j = 1; j < std::min(11, NSAMPLE-1); ++j) {
                cons += (*waveform)[ch_idx][j];
                if (std::fabs((*waveform)[ch_idx][j] - (*waveform)[ch_idx][j + 1]) > 5) {
                    counter = j;
                    break;
                }
            }
            if(counter > 0) cons /= counter;

            // Max ADC
            double max_adc = 0.0;
            for (int j = 0; j < NSAMPLE; ++j) {
                if ((*waveform)[ch_idx][j] > max_adc) max_adc = (*waveform)[ch_idx][j];
            }

            // 積分
            double value = 0.0;
            int imin = std::max(0, IntegralMin);
            int imax = std::min(NSAMPLE-1, IntegralMax);
            for (int j = imin; j <= imax; ++j) {
                value += (*waveform)[ch_idx][j] - cons;
            }

            // ログ出力
            std::cout << "[RAYRAW#" << rr_id << " ch" << lc_id << "] Max ADC: " << max_adc
                      << " | baseline: " << cons
                      << " | integral(" << imin << "-" << imax << "): " << value << std::endl;

            // グラフ描画
            TGraph *g = new TGraph(NSAMPLE, x.data(), (*waveform)[ch_idx].data());
            g->SetTitle(Form("RAYRAW#%d ch%d;Sampling;ADC value", rr_id, lc_id));
            g->SetMarkerColor(1);
            g->SetMarkerStyle(8);
            g->Draw("ALP");
            gPad->Update();
        }

        c1->Update();
        gSystem->ProcessEvents();

        std::cout << "Input q to quit (Enter to next event)." << std::endl;
        std::getline(std::cin, msg);
        if (msg == "q") break;
    }
}

// -------------- (RAYRAW, local_ch) → global_ch 対応表を作る --------------
static std::vector<std::vector<int>> mapRayrawLocalToGlobal(TTree* t) {
    std::vector<std::vector<int>> rayraw_to_global(12); // index 0..11: RAYRAW#=index+1

    std::vector<std::vector<double>> *waveform = nullptr;
    std::vector<int> *plane = nullptr;

    t->SetBranchAddress("waveform", &waveform);
    bool hasPlane = (t->GetBranch("plane") != nullptr);
    if (hasPlane) t->SetBranchAddress("plane", &plane);

    t->GetEntry(0);
    if (!waveform) return rayraw_to_global;

    int nChAll = (int) waveform->size();
    for (int gch = 0; gch < nChAll; ++gch) {
        int pid = (hasPlane && plane && (int)plane->size() > gch) ? plane->at(gch) : gch / 32;
        if (pid < 0 || pid >= 12) continue;
        rayraw_to_global[pid].push_back(gch);
    }
    return rayraw_to_global;
}

// -------------- 新ラッパ：RAYRAW番号とその中のch番号で指定する --------------
void waveformMultiRC(TString filename,
                     int ev_start=0,
                     int rayraw0=-1, int ch0=-1,
                     int rayraw1=-1, int ch1=-1,
                     int rayraw2=-1, int ch2=-1,
                     int rayraw3=-1, int ch3=-1,
                     int rayraw4=-1, int ch4=-1,
                     int rayraw5=-1, int ch5=-1,
                     int rayraw6=-1, int ch6=-1,
                     int rayraw7=-1, int ch7=-1
                    )
{
    TFile *f = TFile::Open(filename);
    if(!f || f->IsZombie()){ std::cerr << "Failed to open file: " << filename << std::endl; return; }
    TTree *t = (TTree*) f->Get("tree");
    if(!t){ std::cerr << "TTree 'tree' not found.\n"; f->Close(); return; }

    auto mapR2G = mapRayrawLocalToGlobal(t);
    f->Close();

    std::vector<std::pair<int,int>> pairs = {
        {rayraw0, ch0}, {rayraw1, ch1}, {rayraw2, ch2}, {rayraw3, ch3},
        {rayraw4, ch4}, {rayraw5, ch5}, {rayraw6, ch6}, {rayraw7, ch7}
    };

    std::vector<int> globals;
    std::vector<int> rayraws;
    std::vector<int> locals;

    for (size_t i=0;i<pairs.size();++i){
        int rr = pairs[i].first;
        int lc = pairs[i].second;
        if (rr < 1 || rr > 12 || lc < 0) continue;

        int pid = rr - 1;
        if ((int)mapR2G[pid].size() == 0) {
            std::cerr << "[WARN] RAYRAW#" << rr << " はこのファイルには存在しない可能性があります。\n";
            continue;
        }
        if (lc >= (int)mapR2G[pid].size()) {
            std::cerr << "[WARN] RAYRAW#" << rr << " の ch" << lc
                      << " は存在しません（この面の実チャンネル数は " << mapR2G[pid].size() << " です）。\n";
            continue;
        }
        int gch = mapR2G[pid][lc];
        globals.push_back(gch);
        rayraws.push_back(rr);
        locals.push_back(lc);

        std::cout << "Mapping: RAYRAW#" << rr << " ch" << lc << " -> global ch " << gch << std::endl;
    }

    if (globals.empty()) {
        std::cerr << "No valid (RAYRAW, ch) pairs mapped to global channels.\n";
        return;
    }

    waveform_impl(filename, globals, rayraws, locals, ev_start);
}

// -------------- 互換API（グローバルchで直接指定） --------------
void waveformMulti(TString filename,
              int ch0=-1,int ch1=-1,int ch2=-1,int ch3=-1,
              int ch4=-1,int ch5=-1,int ch6=-1,int ch7=-1,
              int ev_start=0)
{
    std::vector<int> chs;
    int arr[8] = {ch0,ch1,ch2,ch3,ch4,ch5,ch6,ch7};
    for(int i=0;i<8;++i) if(arr[i]>=0) chs.push_back(arr[i]);

    std::vector<int> dummy_rr(chs.size(), -1);
    std::vector<int> dummy_lc(chs.size(), -1);

    waveform_impl(filename, chs, dummy_rr, dummy_lc, ev_start);
}
