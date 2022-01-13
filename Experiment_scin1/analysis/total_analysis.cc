Double_t timewalkfit(Double_t *x, Double_t *par){
  return -par[0]+par[1]/std::sqrt(x[0]);
}


void total_analysis(){

  //-----------bring offset data-------------------
  Int_t title = 211021;   //have to change whenever you change the data file!!
  Int_t off_title = 210927;
  Int_t N = 6;  //the number of the channels

  TFile *file1 = new TFile(Form("../offset_data/offset_%d.root",off_title),"read");
  TTree *data1 = (TTree*)file1->Get("data");
  Int_t ADC_off[16];

  data1->SetBranchAddress("ADC",ADC_off);

  Double_t total_off = data1->GetEntries();

  //------------making root file for pictures------
  TFile *file_pic = new TFile(Form("../analysis_data/%d_pic.root",title),"recreate");
  
  //-----------pedestal fitting--------------------
  TH1D* hist_off[N];
  for(int i=0;i<N;i++){
    hist_off[i] = new TH1D(Form("hist_off%d",i),Form("ADC_off[%d]",i),300,0,300); 
  }

  for(int n=0;n<total_off;n++){
    data1->GetEntry(n);
    for(int j=0;j<N;j++){
      hist_off[j]->Fill(ADC_off[j]);
    }
  }
  TF1 *f_off[N];
  TCanvas *cv1 = new TCanvas("cv1","ADC pedestal",800,650);
  cv1->Divide(N,1);
  for(int i=0;i<N;i++){
    f_off[i] = new TF1(Form("f_off%d",i),"gaus(0)",0,300);
    cv1->cd(i+1);
    hist_off[i]->Fit(f_off[i],"Q","",0,300);
  }
    
  hist_off[0]->SetTitle("ADC_trig1-1_pedestal;ADC;n");
  hist_off[1]->SetTitle("ADC_trig1-2_pedestal;ADC;n");
  hist_off[2]->SetTitle("ADC_scin1_pedestal;ADC;n");
  hist_off[3]->SetTitle("ADC_scin2_pedestal;ADC;n");
  hist_off[4]->SetTitle("ADC_trig2-1_pedestal;ADC;n");
  hist_off[5]->SetTitle("ADC_trig2-2_pedestal;ADC;n");

  file_pic->cd();
  for(int i=0;i<N;i++){
    hist_off[i]-> Write();
  }
  

  //---------print pedestal mean and sigma-----------------
  Double_t ped[3*N];
  for(int i=0;i<N;i++){
    f_off[i]->GetParameters(&ped[i*3]);
  }

  //have to change whenever you change the setup
  std::cout<<"----------Trigger counter pedestal-------------"<<std::endl;
  std::cout<<"-------------------Trig1-1---------------------"<<std::endl;
  std::cout<<"mean : "<<ped[1]<<std::endl;
  std::cout<<"sigma : "<<ped[2]<<std::endl;
  std::cout<<"-------------------Trig1-2---------------------"<<std::endl;
  std::cout<<"mean : "<<ped[4]<<std::endl;
  std::cout<<"sigma : "<<ped[5]<<std::endl;
  std::cout<<"-------------------Scin1-----------------------"<<std::endl;
  std::cout<<"mean : "<<ped[7]<<std::endl;
  std::cout<<"sigma : "<<ped[8]<<std::endl;
  std::cout<<"-------------------Scin2-----------------------"<<std::endl;
  std::cout<<"mean : "<<ped[10]<<std::endl;
  std::cout<<"sigma : "<<ped[11]<<std::endl;
  std::cout<<"-------------------Trig2-1---------------------"<<std::endl;
  std::cout<<"mean : "<<ped[13]<<std::endl;
  std::cout<<"sigma : "<<ped[14]<<std::endl;
  std::cout<<"-------------------Trig2-2---------------------"<<std::endl;
  std::cout<<"mean : "<<ped[16]<<std::endl;
  std::cout<<"sigma : "<<ped[17]<<std::endl;

  
  
  //---------pedestal reduction----------------

  //---------bring experiment data-------------

  TFile *file2 = new TFile(Form("../data/%d.root",title),"read");
  TTree *data2 = (TTree*)file2->Get("data");

  Int_t ADC_data[16];
  Int_t TDC_data[16];
  Int_t unixtime_data;

  data2->SetBranchAddress("ADC",ADC_data);
  data2->SetBranchAddress("TDC",TDC_data);
  data2->SetBranchAddress("unixtime",&unixtime_data);

  TFile *file3 = new TFile(Form("../analysis_data/%d_reduction.root",title),"recreate");
  TTree *data3 = new TTree("data","data without pedestal");


  Int_t ADC[N];
  Int_t TDC[N];

  data3->Branch("ADC",&ADC,Form("ADC[%d]/I",N));  //Branch(branchname, &value);
  data3->Branch("TDC",&TDC,Form("TDC[%d]/I",N));

  file3->cd();
  Double_t total_data = data2->GetEntries();
  for(int n=0;n<total_data;n++){
    data2->GetEntry(n);
    if(TDC_data[0]!=-9999){
      for(int i=0;i<N;i++){
	ADC[i] = ADC_data[i]-ped[3*i+1];
	TDC[i] = TDC_data[i];
      }
      if(ADC[0]>0&&ADC[1]>0&&ADC[2]>0&&ADC[3]>0&&ADC[4]>0&&ADC[5]>0){
	//if(unixtime_data>1633800000&&unixtime_data<1634300000){
	  data3->Fill();
	  //}
      }
    }
  }
  data3->Write();
  file3->Close();

  
  //----------weird value suppression-------- (ADC values under 0)

  TFile *file8 = new TFile(Form("../analysis_data/%d_reduction.root",title),"read");
  TTree *data8 = (TTree*) file8->Get("data");

  Int_t ADC_red[N];
  Int_t TDC_red[N];

  data8->SetBranchAddress("ADC",ADC_red);
  data8->SetBranchAddress("TDC",TDC_red);

  Int_t total_red1 = data8->GetEntries();

  TH1D *adc_check[N];
  TF1 *adc_check_fit[N];
  for(int i=0;i<N;i++){
    adc_check[i] = new TH1D(Form("adc_check%d",i),"checking ADC",350,-500,3000);
  }
  for(int n=0;n<total_red1;n++){
    data8->GetEntry(n);
    if(TDC_red[1]!=-9999&&TDC_red[2]!=-9999&&TDC_red[3]!=-9999&&TDC_red[4]!=-9999){
      for(int i=0;i<N;i++){
	adc_check[i]->Fill(ADC_red[i]);
      }
    }
  }


  TCanvas *cv5 = new TCanvas("cv5","ADC after pedestal reduction",1000,650);
  cv5->Divide(N,1);
  for(int i=0;i<N;i++){
    adc_check_fit[i] = new TF1(Form("adc_check_fit[%d]",i),"landau(0)",-500,3000);
    cv5->cd(i+1);
    adc_check[i]->Fit(Form("adc_check_fit[%d]",i),"","",-500,3000);
  }
  Double_t adc_cons[3*N];
  for(int i=0;i<N;i++){
    adc_check_fit[i]->GetParameters(&adc_cons[3*i]);
  }
  
  //---------sorting and timewalk-------------

  TFile *file4 = new TFile(Form("../analysis_data/%d_reduction.root",title),"read");
  TTree *data4 = (TTree*) file4->Get("data");

  Int_t ADC_reduction[N];
  Int_t TDC_reduction[N];

  data4->SetBranchAddress("ADC",ADC_reduction);
  data4->SetBranchAddress("TDC",TDC_reduction);

  Int_t total_red = data4->GetEntries();
  TH2D *adc_tdc[N];
  TH2D *adc_tdc_c[N];

  
  //--------here you have to change all histogram name whenever you change your experimental setup------
  //--------making time-walk, time-walk correction histogram--------------------------------------------
  adc_tdc[0] = new TH2D("adc_tdc0","Trig1-1 adc_tdc time-walk",210,-100,2000,190,100,2000);
  adc_tdc_c[0] = new TH2D("adc_tdc_c0","Trig1-1 adc_tdc time-walk correction",210,-100,2000,100,-500,500);
  adc_tdc[1] = new TH2D("adc_tdc1","Trig1-2 adc_tdc time-walk",210,-100,2000,190,100,2000);
  adc_tdc_c[1] = new TH2D("adc_tdc_c1","Trig1-2 adc_tdc time-walk correction",210,-100,2000,100,-500,500);
  adc_tdc[2] = new TH2D("adc_tdc2","Scin1 adc_tdc time-walk",210,-100,2000,190,100,2000);
  adc_tdc_c[2] = new TH2D("adc_tdc_c2","Scin1 adc_tdc time-walk correction",210,-100,2000,100,-500,500);
  adc_tdc[3] = new TH2D("adc_tdc3","Scin2 adc_tdc time-walk",210,-100,2000,190,100,2000);
  adc_tdc_c[3] = new TH2D("adc_tdc_c3","Scin2 adc_tdc time-walk correction",210,-100,2000,100,-500,500);
  adc_tdc[4] = new TH2D("adc_tdc4","Trig2-1 adc_tdc time-walk",210,-100,2000,190,100,2000);
  adc_tdc_c[4] = new TH2D("adc_tdc_c4","Trig2-1 adc_tdc time-walk correction",210,-100,2000,100,-500,500);
  adc_tdc[5] = new TH2D("adc_tdc5","Trig2-2 adc_tdc time-walk",210,-100,2000,190,100,2000);
  adc_tdc_c[5] = new TH2D("adc_tdc_c5","Trig2-2 adc_tdc time-walk correction",210,-100,2000,100,-500,500);

  //new
  TH2D *adc_tdc1_c[N];
  adc_tdc1_c[0] = new TH2D("adc_tdc1_0","Trig1-1",210,-100,2000,100,-500,500);
  adc_tdc1_c[1] = new TH2D("adc_tdc1_1","Trig1-2",210,-100,2000,100,-500,500);
  adc_tdc1_c[2] = new TH2D("adc_tdc1_2","Scin1",210,-100,2000,100,-500,500);
  adc_tdc1_c[3] = new TH2D("adc_tdc1_3","Scin2",210,-100,2000,100,-500,500);
  adc_tdc1_c[4] = new TH2D("adc_tdc1_4","Trig2-1",210,-100,2000,100,-500,500);
  adc_tdc1_c[5] = new TH2D("adc_tdc1_5","Trig2-2",210,-100,2000,100,-500,500);
  
  for(int n=0;n<total_red;n++){
    data4->GetEntry(n);
    if(TDC_reduction[0]!=-9999&&TDC_reduction[1]!=-9999&&TDC_reduction[2]!=-9999&&TDC_reduction[3]!=-9999&&TDC_reduction[4]!=-9999){ //sorting condition
      if(ADC_reduction[2]>adc_cons[7]-3*adc_cons[8]&&ADC_reduction[2]<adc_cons[7]+20*adc_cons[8]){   //proper values of scintillator 1
	if(ADC_reduction[0]>adc_cons[1]-3*adc_cons[2]&&ADC_reduction[0]<adc_cons[1]+20*adc_cons[2]){   //proper values of trig1-1
	  for(int i=0;i<N;i++){
	    adc_tdc[i]->Fill(ADC_reduction[i],TDC_reduction[i]);
	  }
	}
      }
    }
  }
  TF1 *f_time[N];
  Double_t time[2*N];
  TCanvas *cv2 = new TCanvas("cv2","ADC-TDC 2D histogram for time-walk correction",1000,650);
  cv2->Divide(N,2);
  for(int i=0;i<N;i++){
    f_time[i] = new TF1(Form("f_time[%d]",i),timewalkfit,0,2000,2);
    cv2->cd(i+1);
    adc_tdc[i] -> Fit(Form("f_time[%d]",i),"","",adc_cons[3*i+1]-3*adc_cons[3*i+2],adc_cons[3*i+1]+20*adc_cons[3*i+2]);
    adc_tdc[i]->Draw("colz");
    adc_tdc[i]->SetTitle(Form("adc_tdc[%d];ADC;TDC",i));
    f_time[i]->GetParameters(&time[2*i]);
    //cout<<time[2*i]<<" and "<<time[2*i+1]<<endl;
					   
  }

  TFile *file5 = new TFile(Form("../analysis_data/%d_correction.root",title),"recreate");
  TTree *data5 = new TTree("data","data with time-walk correction");

  Int_t ADC_c[N];
  Int_t TDC_c[N];

  data5->Branch("ADC",&ADC_c,Form("ADC[%d]/I",N));
  data5->Branch("TDC",&TDC_c,Form("TDC[%d]/I",N));
  file4->cd();

  for(int n=0; n<total_red;n++){
    file4->cd();
    data4->GetEntry(n);
    if(TDC_reduction[0]!=-9999&&TDC_reduction[1]!=-9999&&TDC_reduction[2]!=-9999&&TDC_reduction[3]!=-9999&&TDC_reduction[4]!=-9999){ //sorting condition
      if(ADC_reduction[0]>adc_cons[1]-3*adc_cons[2]&&ADC_reduction[0]<adc_cons[1]+20*adc_cons[2]){
	if(ADC_reduction[2]>adc_cons[7]-3*adc_cons[8]&&ADC_reduction[2]<adc_cons[7]+20*adc_cons[8]){
	  for(int i=0;i<N;i++){
	    TDC_c[i] = TDC_reduction[i]+time[2*i]-(time[2*i+1]/std::sqrt(ADC_reduction[i]));
	    adc_tdc_c[i]->Fill(ADC_reduction[i],TDC_c[i]);
	    adc_tdc1_c[i]->Fill(ADC_reduction[i],TDC_c[1]); //new

	    ADC_c[i] = ADC_reduction[i];
	  }
	  file5->cd();
	  data5->Fill();
	  }
      }
    }
  }
  file5->cd();
  data5->Write();
  file5->Close();

  for(int i=0;i<N;i++){
    cv2->cd(i+1+N);
    adc_tdc_c[i]->Draw("colz");
    adc_tdc_c[i]->SetTitle(Form("adc_tdc_c[%d];ADC;TDC",i));
  }

    //new
    /*TCanvas *cv10 = new TCanvas("cv10","cv10",1000,650);
      cv10->Divide(N);
      for(int i=0;i<N;i++){
      cv10->cd(i+1);
      adc_tdc1_c[i]->Draw("colz");
      adc_tdc1_c[i]->SetTitle(Form("adc_tdc1_c[%d];ADC;TDC",i));
      }
    */
  
  
    file_pic->cd();
    for(int i=0;i<N;i++){
      adc_tdc[i]-> Write();
      adc_tdc_c[i]->Write();
    }
    
  
    //----------time resolution-----------------------
    TFile *file6 = new TFile(Form("../analysis_data/%d_correction.root",title),"read");
    TTree *data6 = (TTree*) file6->Get("data");

    Int_t TDC_correction[N];

    data6->SetBranchAddress("TDC",TDC_correction);

    Int_t total_cor = data6->GetEntries();

    //have to change whenever you change your setup
    TH1D *hist_tdc[N];
    hist_tdc[0] = new TH1D("hist_tdc1","tdc Trig1-1",100,-500,500);
    hist_tdc[1] = new TH1D("hist_tdc2","tdc Trig1-2",100,-500,500);
    hist_tdc[2] = new TH1D("hist_tdc3","tdc Scin1",100,-500,500);
    hist_tdc[3] = new TH1D("hist_tdc4","tdc Scin2",100,-500,500);
    hist_tdc[4] = new TH1D("hist_tdc5","tdc Trig2-1",100,-500,500);
    hist_tdc[5] = new TH1D("hist_tdc6","tdc Trig2-2",100,-500,500);

    TF1 *tdc_fit[N];

    for(int n=0;n<total_cor;n++){
      data6->GetEntry(n);
      for(int i=0;i<N;i++){
	hist_tdc[i]->Fill(TDC_correction[i]);  
      }
    }

    Double_t tdc_par[3*N];
    TCanvas *cv3 = new TCanvas("cv3","TDC after time-walk correction",1000,650);
    cv3->Divide(N,1);
    for(int i=0;i<N;i++){
      tdc_fit[i] = new TF1(Form("tdc_fit[%d]",i),"gaus(0)",-500,500);
      cv3->cd(i+1);
      hist_tdc[i]->Fit(Form("tdc_fit[%d]",i),"Q","",-500,500);
      tdc_fit[i]->GetParameters(&tdc_par[3*i]);
    }
  
    //print time resolution, tdc mean value
    //have to change whenever you change the setup
    std::cout<<"--------Trigger counter time resolution--------"<<std::endl;
    std::cout<<"------------------channel----------------------"<<std::endl;
    std::cout<<"-------------------Trig1-1---------------------"<<std::endl;
    std::cout<<"mean : "<<tdc_par[1]<<std::endl;
    std::cout<<"sigma : "<<tdc_par[2]<<std::endl;
    std::cout<<"-------------------Trig1-2---------------------"<<std::endl;
    std::cout<<"mean : "<<tdc_par[4]<<std::endl;
    std::cout<<"sigma : "<<tdc_par[5]<<std::endl;
    std::cout<<"-------------------Scin1-----------------------"<<std::endl;
    std::cout<<"mean : "<<tdc_par[7]<<std::endl;
    std::cout<<"sigma : "<<tdc_par[8]<<std::endl;
    std::cout<<"-------------------Scin2-----------------------"<<std::endl;
    std::cout<<"mean : "<<tdc_par[10]<<std::endl;
    std::cout<<"sigma : "<<tdc_par[11]<<std::endl;
    std::cout<<"-------------------Trig2-1---------------------"<<std::endl;
    std::cout<<"mean : "<<tdc_par[13]<<std::endl;
    std::cout<<"sigma : "<<tdc_par[14]<<std::endl;
    std::cout<<"-------------------Trig2-2---------------------"<<std::endl;
    std::cout<<"mean : "<<tdc_par[16]<<std::endl;
    std::cout<<"sigma : "<<tdc_par[17]<<std::endl;

  
    file_pic->cd();
    for(int i=0;i<N;i++){
      hist_tdc[i]-> Write();
    }

    
 

    //-----------ADC value fitting-----------------------------

    //have to change whenever you change your setup
  /*
  Int_t ADC_correction[N];

  data6->SetBranchAddress("ADC",ADC_correction);
  
  TH1D *hist_adc[5];
  hist_adc[0] = new TH1D("hist_adc1","adc Trig0",2500,0,2500);
  hist_adc[1] = new TH1D("hist_adc2","adc Trig1-1",2500,0,2500);
  hist_adc[2] = new TH1D("hist_adc3","adc Trig1-2",2500,0,2500);
  hist_adc[3] = new TH1D("hist_adc4","adc Trig2-1",2500,0,2500);
  hist_adc[4] = new TH1D("hist_adc5","adc Trig2-2",2500,0,2500);

  TF1 *adc_fit[5];

  for(int n=0;n<total_cor;n++){
    data6->GetEntry(n);
    for(int i=0;i<N;i++){
      hist_adc[i]->Fill(ADC_correction[i]);  
    }
  }

  Double_t adc_par[15];
  TCanvas *cv4 = new TCanvas("cv4","cv4",1000,650);
  cv4->Divide(N,1);
  for(int i=0;i<N;i++){
    adc_fit[i] = new TF1(Form("adc_fit[%d]",i),"landau(0)",0,2500);
    cv4->cd(i+1);
    hist_adc[i]->Fit(Form("adc_fit[%d]",i),"","",0,2500);
    adc_fit[i]->GetParameters(&adc_par[3*i]);
  }

  file_pic->cd();
  for(int i=0;i<N;i++){
    hist_adc[i]-> Write();
  }

  TFile *file7 = new TFile(Form("../trigger_data_1/%d_cor_adccut.root",title),"recreate");
  TTree *data7 = new TTree("data","trigger data with time-walk correction and adc cut");

  Int_t ADC_cut[N];
  Int_t TDC_cut[N];

  data7->Branch("ADC",&ADC_cut,Form("ADC[%d]/I",N));
  data7->Branch("TDC",&TDC_cut,Form("TDC[%d]/I",N));
  //file6->cd();

  for(int n=0; n<total_cor;n++){
    //file6->cd();
    data6->GetEntry(n);
    if(ADC_correction[0]>(adc_par[1]-adc_par[2]) && ADC_correction[0]<(adc_par[1]+3*adc_par[2])&&ADC_correction[1]>adc_par[4]-adc_par[5]&&ADC_correction[1]<adc_par[4]+3*adc_par[5]&&ADC_correction[2]>adc_par[7]-adc_par[8]&&ADC_correction[2]<adc_par[7]+3*adc_par[8]&&ADC_correction[3]>adc_par[10]-adc_par[11]&&ADC_correction[3]<adc_par[10]+3*adc_par[11]&&ADC_correction[4]>adc_par[13]-adc_par[14]&&ADC_correction[4]<adc_par[13]+3*adc_par[14]){ //ADC_cut condition
      for(int i=0;i<N;i++){
	TDC_cut[i] = TDC_correction[i];
	ADC_cut[i] = ADC_correction[i];
      }
    
      //file7->cd();
      data7->Fill();
      }
  }
  file7->cd();
  data7->Write();
  file7->Close();
  */
  
  //---time resolution after considering the hit position-----
  //here you have to change the code whenever you change the setup


  /*TH1D *hit[2];
  TF1 *f_pos[2];
  for(int i=0;i<2;i++){
    hit[i] = new TH1D(Form("hit%d",i),Form("hit position trig%d",i+1),200,-100,100);
    f_pos[i] = new TF1(Form("f_pos[%d]",i),"gaus(0)",-100,100);
  }

  TH2D *hit2 = new TH2D("hit2","hit position of trig 1 and 2",200,-100,100,200,-100,100);

  for(int n=0;n<total_cor;n++){
    data6->GetEntry(n);
    hit[0]->Fill(TDC_correction[2]-TDC_correction[1]);
    hit[1]->Fill(TDC_correction[4]-TDC_correction[3]);
    hit2->Fill(TDC_correction[2]-TDC_correction[1],TDC_correction[4]-TDC_correction[3]);
  }

  TCanvas *cv4 = new TCanvas("cv4","cv4",1000,650);
  cv4->Divide(3,1);
  cv4->cd(1);
  hit[0]->Fit("f_pos[0]","","",-100,100);
  cv4->cd(2);
  hit[1]->Fit("f_pos[1]","","",-100,100);
  cv4->cd(3);
  hit2->Draw();

  Double_t hit_par[6];
  f_pos[0]->GetParameters(&hit_par[0]);
  f_pos[1]->GetParameters(&hit_par[3]);

  TH1D *hist_tdc_pos[5];
  hist_tdc_pos[0] = new TH1D("hist_tdc_pos1","tdc considering hit position Trig0",400,-200,200);
  hist_tdc_pos[1] = new TH1D("hist_tdc_pos2","tdc considering hit position Trig1-1",400,-200,200);
  hist_tdc_pos[2] = new TH1D("hist_tdc_pos3","tdc considering hit position Trig1-2",400,-200,200);
  hist_tdc_pos[3] = new TH1D("hist_tdc_pos4","tdc considering hit position Trig2-1",400,-200,200);
  hist_tdc_pos[4] = new TH1D("hist_tdc_pos5","tdc considering hit position Trig2-2",400,-200,200);

  TF1 *tdc_fit_pos[5];

  Double_t time_dif_1;
  Double_t time_dif_2;
  
  for(int n=0;n<total_cor;n++){
    data6->GetEntry(n);
    time_dif_1 = TDC_correction[2]-TDC_correction[1];
    time_dif_2 = TDC_correction[4]-TDC_correction[3];
    if(time_dif_1>(hit_par[1]-hit_par[2]/2) && time_dif_1<(hit_par[1]+hit_par[2]/2) && time_dif_2>(hit_par[4]-hit_par[5]/2) && time_dif_2<(hit_par[4]+hit_par[5]/2)){
      for(int i=0;i<N;i++){
	hist_tdc_pos[i]->Fill(TDC_correction[i]);
      }
    }
  }

  TCanvas *cv5 = new TCanvas("cv5","cv5",1000,650);
  cv5->Divide(N,1);
  for(int i=0;i<N;i++){
    tdc_fit_pos[i] = new TF1(Form("tdc_fit_pos[%d]",i),"gaus(0)",-200,200);
    cv5->cd(i+1);
    hist_tdc_pos[i]->Fit(Form("tdc_fit_pos[%d]",i),"","",-200,200);
    //tdc_fit[i]->GetParameters(&tdc_par[3*i]);
  }
  */
  
  
  
  
    


 
  
  
  

  
  
  
  
  
  
  

  

  
}
  
