//This program creates a macro for finding Average_pulse_integral and also to create pdf. This program also makes use of Form and shows how it's used.

void Average_pulse_integral( Int_t nrun=14464 ) {



  gStyle->SetCanvasColor(0);
  gStyle->SetFrameFillColor(0);

  

  TFile *f = new TFile(Form("shms_coin_replay_production_%d_50000.root", nrun)); 
  TTree *t = (TTree*) f->Get("T");
  
  Int_t NSampWaveFormP;
  t->SetBranchAddress("Ndata.P.cal.pr.adcPosSampWaveform",&NSampWaveFormP) ;
  Double_t SampWaveFormP[10000];
  t->SetBranchAddress("P.cal.pr.adcPosSampWaveform",&SampWaveFormP) ;
  
  Int_t NSampWaveFormN;
  t->SetBranchAddress("Ndata.P.cal.pr.adcNegSampWaveform",&NSampWaveFormN) ;
  Double_t SampWaveFormN[10000];
  t->SetBranchAddress("P.cal.pr.adcNegSampWaveform",&SampWaveFormN) ;


  const int nblocks = 14;
  auto h_hitpos= new TH1F("h_hit_pos","hitpos;PMT number; positive hits", 16,-0.5, 15.5);
  auto h_hitneg= new TH1F("h_hit_neg","hitneg;PMT number; negative hits", 16,-0.5, 15.5);


  Double_t pulse_integral_pos[nblocks];
  Double_t pulse_integral_neg[nblocks];
  for (Int_t n=0;n<nblocks;n++)
    {
      pulse_integral_pos[n]=0;
      pulse_integral_neg[n]=0;
    }
  
  Int_t SampPMT;
  Int_t NSamp;
  Double_t SampADC;

  Long64_t nentries = t->GetEntries();    
  cout<<" nentries="<<nentries<<endl;
  int inp;
  for (int i = 0; i < nentries; i++)
    {
      
      Bool_t verbose = kTRUE;
      //   if(i<10) verbose =kTRUE;
      
      if(verbose) cout<<"event="<<i<<endl;
      if((i%1000)==0) cout<<i<<":"<<nentries<<endl;
      
      t->GetEntry(i);
      
      Int_t nsp = 0;
      while (nsp < NSampWaveFormP)
	{
	  SampPMT=SampWaveFormP[nsp++];
	  h_hitpos->Fill(SampPMT);
	  NSamp=SampWaveFormP[nsp++];
	  if(verbose) cout<<" positive side event number ="<<i<<"  SampPMT="<<
		      SampPMT<<"   NSamp="<<NSamp<<endl;
	  for (Int_t n=0;n<NSamp;n++) {
	  SampADC=SampWaveFormP[nsp++];
      	  pulse_integral_pos[SampPMT-1]+=SampADC;
	}	
      }

      if(verbose)
	{
	  for(int i=0;i<14;i++)
	    cout<<"pmt ="<<i<<"pulse_intergal_pos ="<<pulse_integral_pos[i]<<endl;
	}

      Int_t nsn = 0;
      while (nsn < NSampWaveFormN)
	{
	  SampPMT=SampWaveFormN[nsn++];//which PMT
	  h_hitneg->Fill(SampPMT);
	  NSamp=SampWaveFormN[nsn++];//how much information we have for blocks
	  if(verbose) cout<<" negative side event number ="<<i<<"  SampPMT="<<SampPMT<<"   NSamp="<<NSamp<<endl;
	  for (Int_t n=0;n<NSamp;n++)
	    {
	      SampADC=SampWaveFormN[nsn++];//what information
	      pulse_integral_neg[SampPMT-1]+=SampADC;

	    }
	}	
      
    }

  auto h_averagepulsepos= new TH1F("avergae pulse pos","average pulse pos;PMT number;average pulse integral pos", 16,-0.5, 15.5);
  auto h_averagepulseneg= new TH1F("average pulse neg","avergae pulse neg; PMT number;average pulse integral neg", 16,-0.5, 15.5);
  
  for(int i=0;i<nblocks; i++)
    {
      pulse_integral_neg[i]=pulse_integral_neg[i]/h_hitneg->GetBinContent(i+2);
      pulse_integral_pos[i]=pulse_integral_pos[i]/h_hitpos->GetBinContent(i+2);
       cout<<"hitpositive"<<h_hitpos->GetBinContent(i+2)<<endl;
      //cout<<"pulse_integral_pos"<<pulse_integral_neg[i]<<endl;
      cout<<"hitnegative"<<h_hitneg->GetBinContent(i+2)<<endl;
       
      //      cout<<Form(" block i pos: neg  = %4.0f : %4.0f ",pulse_integral_pos[i],pulse_integral_neg[i])<<endl;
      h_averagepulsepos->SetBinContent(i+2, pulse_integral_pos[i]);
      h_averagepulseneg->SetBinContent(i+2, pulse_integral_neg[i]);

    }

 
  
  // for(int i=0;i<nblocks; i++)
  //cout<<i<<"   "<<pulse_integral_neg[i]<<"   "<<h_hitneg->GetBinContent(i+2)<<endl;



  
    for(int i=0;i<nblocks; i++) cout<<Form(" block i pos: neg  = %4.0f : %4.0f ",pulse_integral_pos[i],pulse_integral_neg[i])<<endl;
    
    h_averagepulsepos->SetStats(1);
    h_averagepulseneg->SetStats(1);                
    gStyle->SetOptStat(111111);  
  

    TCanvas *c1= new TCanvas();
    c1->Divide(1,2);
    c1->cd(1);
    h_averagepulsepos->SetMinimum(0);
    h_averagepulsepos->SetFillColor(kBlue-7);
    h_averagepulsepos->SetStats(1);
    gStyle->SetOptStat(111111);

  //TFile f("hsimple.root","READ");
   //c1 = new TCanvas("c1","sub data",200,10,700,500);
   // hpx->Draw();
    h_averagepulsepos->Draw();
    c1->Print("h9.pdf(","pdf");
    //  h_averagepulsepos->Draw();
    //TCanvas *c2 = new TCanvas();
    c1->cd(2);
    h_averagepulseneg->SetMinimum(0);
    h_averagepulseneg->SetFillColor(kRed);
    h_averagepulseneg->SetStats(1);
    gStyle->SetOptStat(111111);
    h_averagepulseneg->Draw(); 
    
    c1->Print("h9.pdf)","pdf");
    //h_averagepulseneg->Draw(); 
    

    
  // --------------------------------------------------------------------------------
  
}
