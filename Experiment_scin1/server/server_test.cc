#include <THttpServer.h>

void server_test(const char* jobname = "job1",Long64_t maxcnt = 0){
  THttpServer* serv = new THttpServer(Form("http:8080?top=%s",jobname));
  serv->SetReadOnly(kFALSE);
  gBenchmark->Start(jobname);
  
  Int_t title = 210924;
  TFile *file3 = new TFile(Form("../analysis_data/%d_correction.root",title),"read");
  TTree *data = (TTree*) file3->Get("data");

  Int_t TDC_a[6];
  Int_t ADC_a[6];

  data->SetBranchAddress("TDC",TDC_a);
  data->SetBranchAddress("ADC",ADC_a);

  int total1 = data ->GetEntries();

  TH1D *TDC1 = new TH1D("tdc1","tdc1",200,-100,100);
  for(int n=0;n<total1;n++){
    data->GetEntry(n);
    TDC1->Fill(TDC_a[5]);
  }
  //TCanvas *c1 = new TCanvas("c1","c1",1000,650);
  //c1->cd();
  //TDC1->Draw();

  gBenchmark->Show(jobname);

  

}
  
				      
