Double_t cut1 = 0.19;
Double_t cut2 = 0.24;

void hitpo(){
  Int_t N = 6; //number of the channels
  Int_t title = 211026;    //have to change whenever you change the data file!
  TFile *file = new TFile(Form("../analysis_data/%d_correction.root",title),"read");
  TTree *data = (TTree*) file->Get("data");
  Int_t TDC[N];
  Int_t ADC[N];

  data->SetBranchAddress("TDC",TDC);
  data->SetBranchAddress("ADC",ADC);
  

  TH1D *hit1[2];
  TH1D *hit2[2];
  TH1D *hit3[2];
  TF1 *hit_fit_ne[3];
  TF1 *hit_fit_po[3];
  for(int i=0;i<2;i++){
    hit1[i] = new TH1D(Form("hit1%d",i),"hit position trig1",280,-3500,3500);
    hit2[i] = new TH1D(Form("hit2%d",i),"hit position Scin",280,-3500,3500); 
    hit3[i] = new TH1D(Form("hit3%d",i),"hit position trig2",280,-3500,3500);
  }
  for(int i=0;i<3;i++){
    hit_fit_ne[i] = new TF1(Form("hit_fit_ne[%d]",i),"[0]*TMath::Erf((x-[1])/[2])+[3]",-3500,0);
    hit_fit_po[i] = new TF1(Form("hit_fit_po[%d]",i),"[0]*TMath::Erfc((x-[1])/[2])+[3]",0,3500);
  }

  int total = data ->GetEntries();

  for(int n=0;n<total;n++){
    data->GetEntry(n);
    for(int i=0;i<2;i++){
      hit1[i]->Fill((TDC[1]-TDC[0])*35./2.);
      hit2[i]->Fill((TDC[3]-TDC[2])*35./2.);
      hit3[i]->Fill((TDC[5]-TDC[4])*35./2.);
    }
  }

  //parameters
  ifstream hitparameters;
  hitparameters.open(Form("parameter/%d_pa.txt",title));
  
  string c;
  int in;
  double var1;
  double var2;
  double trig1_ne_parlim[8];
  double trig1_po_parlim[8];
  double trig1_range[4];

  double trig2_ne_parlim[8];
  double trig2_po_parlim[8];
  double trig2_range[4];

  double scin_ne_parlim[8];
  double scin_po_parlim[8];
  double scin_range[4];

  
  while(hitparameters.peek()!=EOF){
    hitparameters>>c>>in>>var1>>var2;
    if(c=="1-parlim"){
      trig1_ne_parlim[in*2]=var1;
      trig1_ne_parlim[in*2+1]=var2;
    }
    else if(c=="1-range"||c=="1+range"){
      trig1_range[2*in]=var1;
      trig1_range[2*in+1]=var2;
    }
    else if(c=="1+parlim"){
      trig1_po_parlim[in*2]=var1;
      trig1_po_parlim[in*2+1]=var2;
    }
    else if(c=="2-parlim"){
      scin_ne_parlim[in*2]=var1;
      scin_ne_parlim[in*2+1]=var2;
    }
    else if(c=="2-range"||c=="2+range"){
      scin_range[2*in]=var1;
      scin_range[2*in+1]=var2;
    }
    else if(c=="2+parlim"){
      scin_po_parlim[in*2]=var1;
      scin_po_parlim[in*2+1]=var2;
    }
    else if(c=="3-parlim"){
      trig2_ne_parlim[in*2]=var1;
      trig2_ne_parlim[in*2+1]=var2;
    }
    else if(c=="3-range"||c=="3+range"){
      trig2_range[2*in]=var1;
      trig2_range[2*in+1]=var2;
    }
    else if(c=="3+parlim"){
      trig2_po_parlim[in*2]=var1;
      trig2_po_parlim[in*2+1]=var2;
    }
  }
  hitparameters.close();

  //trigger counter 1
  for(int i=0;i<4;i++){
    hit_fit_ne[0]->SetParLimits(i,trig1_ne_parlim[i*2],trig1_ne_parlim[i*2+1]);
    hit_fit_po[0]->SetParLimits(i,trig1_po_parlim[i*2],trig1_po_parlim[i*2+1]);
  }

  hit1[0]->Fit("hit_fit_ne[0]","","",trig1_range[0],trig1_range[1]);
  hit1[1]->Fit("hit_fit_po[0]","","",trig1_range[2],trig1_range[3]);

  double hit_pa1[8];
  hit_fit_ne[0]->GetParameters(&hit_pa1[0]);
  hit_fit_po[0]->GetParameters(&hit_pa1[4]);

  double v1 = cut1/((hit_pa1[5]-hit_pa1[1])*pow(10,-12)); //unit : m/s
  cout<<hit_pa1[5]-hit_pa1[1]<<endl;
  cout<<"experimental velocity "<<v1<<endl;
  cout<<"theoretical velocity "<<299792458/1.58<<" m/s"<<endl;

  //scintillator
  for(int i=0;i<4;i++){
    hit_fit_ne[1]->SetParLimits(i,scin_ne_parlim[i*2],scin_ne_parlim[i*2+1]);
    hit_fit_po[1]->SetParLimits(i,scin_po_parlim[i*2],scin_po_parlim[i*2+1]);
  }

  hit2[0]->Fit("hit_fit_ne[1]","","",scin_range[0],scin_range[1]);
  hit2[1]->Fit("hit_fit_po[1]","","",scin_range[2],scin_range[3]);

  double hit_pa2[8];
  hit_fit_ne[1]->GetParameters(&hit_pa2[0]);
  hit_fit_po[1]->GetParameters(&hit_pa2[4]);

  double v2 = cut2/((hit_pa2[5]-hit_pa2[1])*pow(10,-12)); //unit : m/s
  cout<<hit_pa2[5]-hit_pa2[1]<<endl;
  cout<<"experimental velocity "<<v2<<endl;
  cout<<"theoretical velocity "<<299792458/1.58<<" m/s"<<endl;

  //trigger counter 2
  for(int i=0;i<4;i++){
    hit_fit_ne[2]->SetParLimits(i,trig2_ne_parlim[i*2],trig2_ne_parlim[i*2+1]);
    hit_fit_po[2]->SetParLimits(i,trig2_po_parlim[i*2],trig2_po_parlim[i*2+1]);
  }

  hit3[0]->Fit("hit_fit_ne[2]","","",trig2_range[0],trig2_range[1]);
  hit3[1]->Fit("hit_fit_po[2]","","",trig2_range[2],trig2_range[3]);

  double hit_pa3[8];
  hit_fit_ne[2]->GetParameters(&hit_pa3[0]);
  hit_fit_po[2]->GetParameters(&hit_pa3[4]);

  double v3 = cut1/((hit_pa3[5]-hit_pa3[1])*pow(10,-12)); //unit : m/s
  cout<<hit_pa3[5]-hit_pa3[1]<<endl;
  cout<<"experimental velocity "<<v3<<endl;
  cout<<"theoretical velocity "<<299792458/1.58<<" m/s"<<endl;

  
  

  TCanvas *c1 = new TCanvas("c1","c1",1000,650);
  c1->Divide(6,1);
  c1->cd(1);
  gPad->SetLogy();
  hit1[0]->Draw();

  c1->cd(2);
  gPad->SetLogy();
  hit1[1]->Draw();

  c1->cd(3);
  gPad->SetLogy();
  hit2[0]->Draw();

  c1->cd(4);
  gPad->SetLogy();
  hit2[1]->Draw();

  c1->cd(5);
  gPad->SetLogy();
  hit3[0]->Draw();

  c1->cd(6);
  gPad->SetLogy();
  hit3[1]->Draw();

  /*
  TFile *file2 = new TFile(Form("../analysis_data/%d_hitposition.root",title),"recreate");
  TTree *data2 = new TTree("data","hit position test");
  
  Int_t TDC_hit[N];
  Int_t ADC_hit[N];
  Double_t hitpos[2];

  Int_t t12;
  Int_t t34;
  data2->Branch("ADC",&ADC_hit,"ADC[6]/I");
  data2->Branch("TDC",&TDC_hit,"TDC[6]/I");
  data2->Branch("HitPosition",&hitpos,"HitPosition[2]/D");
  
  for(int n=0;n<total;n++){
    data->GetEntry(n);
    t12 = TDC[2]-TDC[1];
    t34 = TDC[4]-TDC[3];
    if(fabs(t12)*35*pow(10,-12)*v1<cut1){
      if(fabs(t34)*35*pow(10,-12)*v2<cut2){
        if(fabs(t12+t34*19/12)<10){
	  for(int i=0;i<5;i++){
	    TDC_hit[i] = TDC[i];
	    ADC_hit[i] = ADC[i];
	  }
	  hitpos[0] = t12*35*pow(10,-12)*v1*500;
	  hitpos[1] = t34*35*pow(10,-12)*v2*500;
	  data2->Fill();
	  }
      }
    }
  }
  data2->Write();
  file2->Close();
  */


  TH1D *tof[3];
  tof[0] = new TH1D("tof01","tof01",100,-100,100);
  tof[1] = new TH1D("tof02","tof02",100,-100,100);
  tof[2] = new TH1D("tof12","tof12",100,-100,100);

  TF1 *tof_fit[3];
  tof_fit[0] = new TF1("tof01_fit","gaus(0)",-100,100);
  tof_fit[1] = new TF1("tof02_fit","gaus(0)",-100,100);
  tof_fit[2] = new TF1("tof12_fit","gaus(0)",-100,100);

  for(int n=0;n<total;n++){
    data->GetEntry(n);
    if(v1*fabs(TDC[1]-TDC[0])*35.*pow(10,-12)<0.1&&v2*fabs(TDC[3]-TDC[2])*35.*pow(10,-12)<0.1&&v3*fabs(TDC[5]-TDC[4])*35.*pow(10,-12)<0.1){
    tof[0]->Fill((TDC[0]+TDC[1])/2-(TDC[2]+TDC[3])/2);
    tof[1]->Fill((TDC[0]+TDC[1])/2-(TDC[4]+TDC[5])/2);
    tof[2]->Fill((TDC[2]+TDC[3])/2-(TDC[4]+TDC[5])/2);
    }
  }

  TCanvas *c3 = new TCanvas("c3","c3",1000,650);
  c3->Divide(3,1);
  c3->cd(1);
  tof[0]->Fit(tof_fit[0],"","",-100,100);
  c3->cd(2);
  tof[1]->Fit(tof_fit[1],"","",-100,100);
  c3->cd(3);
  tof[2]->Fit(tof_fit[2],"","",-100,100);



  //-------------------------
  /*
  TFile *file3 = new TFile(Form("../../trigger_data_1/%d_hitposition.root",title),"read");
  TTree *data3 = (TTree*) file3->Get("data");

  Int_t TDC[5];
  Int_t ADC_a[5];
  Double_t hit_a[2];

  data3->SetBranchAddress("TDC",TDC);
  data3->SetBranchAddress("ADC",ADC_a);
  data3->SetBranchAddress("HitPosition",hit_a);

  TH1D *hitposition[2];
  for(int i=0;i<2;i++){
    hitposition[i] = new TH1D(Form("hitposition%d",i),"hit position",20,-100,100);
   }

  int total1 = data3 ->GetEntries();

  for(int n=0;n<total1;n++){
    data3->GetEntry(n);
    hitposition[0]->Fill(hit_a[0]);
    hitposition[1]->Fill(hit_a[1]);
  }

  TCanvas *c2 = new TCanvas("c2","c2",1000,650);
  c2->Divide(2,1);
  c2->cd(1);
  hitposition[0]->Draw();

  c2->cd(2);
  hitposition[1]->Draw();

  TH1D *tof[3];
  tof[0] = new TH1D("tof01","tof01",20,-100,100);
  tof[1] = new TH1D("tof02","tof02",20,-100,100);
  tof[2] = new TH1D("tof12","tof12",20,-40,40);

  TF1 *tof_fit[3];
  tof_fit[0] = new TF1("tof01_fit","gaus(0)",-100,100);
  tof_fit[1] = new TF1("tof02_fit","gaus(0)",-100,100);
  tof_fit[2] = new TF1("tof12_fit","gaus(0)",-40,40);

  for(int n=0;n<total1;n++){
    data3->GetEntry(n);
    if(fabs(hit_a[0])<50&&fabs(hit_a[1])<50){
      tof[0]->Fill(TDC[0]-(TDC[1]+TDC[2])/2);
      tof[1]->Fill(TDC[0]-(TDC[3]+TDC[4])/2);
      tof[2]->Fill((TDC[1]+TDC[2])/2-(TDC[3]+TDC[4])/2);
    }
  }

  TCanvas *c3 = new TCanvas("c3","c3",1000,650);
  c3->Divide(3,1);
  c3->cd(1);
  tof[0]->Fit(tof_fit[0],"","",-100,100);
  c3->cd(2);
  tof[1]->Fit(tof_fit[1],"","",-100,100);
  c3->cd(3);
  tof[2]->Fit(tof_fit[2],"","",-40,40);
  */
  

}
