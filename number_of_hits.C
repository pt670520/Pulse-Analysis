void number_of_hits( Int_t nrun=14464 ) {

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetOptStat(0);

  TString pdfname=Form("number_of_hit_%d.pdf",nrun);
  
  TFile *f = new TFile(Form("shms_coin_replay_production_%d_50000.root", nrun)); 
  TTree *t = (TTree*) f->Get("T");

  const int nblock_row=14;
  const int nblock_column= 2;
  const int data_per_block=177; 
  
  Int_t NSampWaveFormP;
  t->SetBranchAddress("Ndata.P.cal.pr.adcPosSampWaveform",&NSampWaveFormP) ;
  Double_t SampWaveFormP[data_per_block*nblock_row+10]; // this value 4956 is a large number, on block hit is theoritically 177 so data for one row of 14 is 177*14. the additional +10 is precautionary
  t->SetBranchAddress("P.cal.pr.adcPosSampWaveform",&SampWaveFormP) ;
  
  Int_t NSampWaveFormN;
  t->SetBranchAddress("Ndata.P.cal.pr.adcNegSampWaveform",&NSampWaveFormN) ;
  Double_t SampWaveFormN[data_per_block*nblock_row+10];
  t->SetBranchAddress("P.cal.pr.adcNegSampWaveform",&SampWaveFormN) ;

  TH1F *hit_per_column[nblock_column];
  for(int icol=0;icol<nblock_column;icol++)
    hit_per_column[icol]=new TH1F(Form("hits_in_column_%d",icol),
				  Form("hits in column %d for run %d; row number; hits",icol, nrun),
				  nblock_row+2, -0.5, nblock_row+1.5);
  auto h= new TH2F("hit on the calo", Form("Recorded hits for run %d;column number;row number",nrun),
		   nblock_column+2, -0.5, nblock_column+1.5,
		   nblock_row+2, -0.5, nblock_row+1.5);
  auto h_numberofPMTS = new TH1F("h_numberofPMTS",
				 Form("Number of Block hit per event for run %d;number of PMTS hit;",nrun),
				 10,-0.5,9.5);

  Int_t PMT_number;
  Int_t SampPMT_N;
  Int_t NSamp;
  Double_t SampADC;
  
  Long64_t nentries = t->GetEntries();    
  cout<<" nentries="<<nentries<<endl;
  int inp;
  for (int i = 0; i < nentries; i++)
    {
      Bool_t verbose = kFALSE;
      // if(i<10)
      // 	{verbose =kTRUE;}
      if(verbose)
	{
	  cout<<" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
	  cout<<"event="<<i<<endl;
	}

      if((i%1000)==0) cout<<i<<":"<<nentries<<endl;
      
      t->GetEntry(i);
      Int_t Block_hit_per_event_per_column[nblock_column];
      for(int icol= 0; icol<nblock_column;icol++)
	Block_hit_per_event_per_column[icol]=0; 
      
      // this is for the first column (or column P in the HMS parlance)
     Int_t nsp = 0;
     while (nsp < NSampWaveFormP)
       {
	 PMT_number=SampWaveFormP[nsp++]; // also the pmt number (or row) starts at 1 not 0
	 hit_per_column[0]->Fill(PMT_number); // the positive column is the first column
	 NSamp=SampWaveFormP[nsp++];
	 nsp+=NSamp; //skip all the info about the traces
	 Block_hit_per_event_per_column[0]+=1;   
	 if(verbose)
	   cout<<" positive side event number ="<<i<<"  PMT number="<<PMT_number <<"   NSamp="<<NSamp<<"   nsp ="<<nsp
	       <<"  Block_hit_per_event_per_column[0]="<<Block_hit_per_event_per_column[0]<< endl;
       }

     // this is for the second column (or column N in the HMS parlance)
      Int_t nsn = 0;
      while (nsn < NSampWaveFormN)
	{
	PMT_number=SampWaveFormN[nsn++];//which PMT
	hit_per_column[1]->Fill(PMT_number);
	NSamp=SampWaveFormN[nsn++];//how much information we have for blocks
	nsn+=NSamp; //skip all the info about the traces
	Block_hit_per_event_per_column[1]+=1;
       	if(verbose)
	  cout<<" negative side event number ="<<i<<"  PMT number="<<PMT_number <<"   NSamp="<<NSamp<<"   nsn ="<<nsn
	       <<"  Block_hit_per_event_per_column[1]="<<Block_hit_per_event_per_column[1]<< endl;
      }

      Int_t total_number_of_hit_for_this_event=0;
      for(int ib= 0; ib<nblock_column;ib++)
	total_number_of_hit_for_this_event+=Block_hit_per_event_per_column[ib];
      h_numberofPMTS->Fill(total_number_of_hit_for_this_event);
    }
  ////////// end of the event loop //////////
  // construct the 2D map of the hits here
  for(int icol=0;icol<nblock_column;icol++)
    for(int irow=0;irow<nblock_row;irow++)
      {
	h->Fill(icol+1,irow+1,hit_per_column[icol]->GetBinContent(irow+2));
      }
  
  //// now onto display
  TCanvas* c1 = new TCanvas();
  h->Draw("colz");
  c1->Print(pdfname+"(");
  
  TCanvas *c2 = new TCanvas();
  h_numberofPMTS->SetMinimum(0);
  h_numberofPMTS->SetFillColor(kGreen);
  h_numberofPMTS->SetStats(1);
  h_numberofPMTS->Draw();
    c2->Print(pdfname);

  
  TCanvas *c3= new TCanvas();
  c3->Divide(2,1);
  for(int icol=0;icol<nblock_column;icol++)
    {
      c3->cd(icol+1);
      hit_per_column[icol]->SetMinimum(0);             
      hit_per_column[icol]->SetFillColor(kGreen);
      hit_per_column[icol]->Draw();
    }
      c3->Print(pdfname+")");

 
    
  // --------------------------------------------------------------------------------
  
}
