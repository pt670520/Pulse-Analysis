//This scripts creates bunch of events(In this case it's 25 events) for PMT 8 for positive side PMTs.
void specific_PMT_NPS_EEL108( Int_t nrun = 55, Int_t PMT=22, Int_t istart=0)
{
  gStyle->SetTitleOffset(1.1, "X");                                                
  gStyle->SetTitleOffset(1.08, "Y");                                               
  gStyle->SetTitleFont(42,"X");                                                    
  gStyle->SetTitleFont(42,"Y");                                                    
  gStyle->SetTitleSize(0.055,"X");                                                 
  gStyle->SetTitleSize(0.055,"Y");

  gStyle->SetTitleColor(kRed);
  
  gStyle->SetLabelOffset(0.01, "X");                                               
  gStyle->SetLabelOffset(0.01, "Y");                                               
  gStyle->SetLabelFont(42,"X");                                                    
  gStyle->SetLabelFont(42,"Y");                                                    
  gStyle->SetLabelSize(0.075,"X");                                                 
  gStyle->SetLabelSize(0.075,"Y");
  
  Int_t NROWS = 16;
  Int_t NCOLS = 15;
  Int_t sample_size = 100;
  const Int_t nblocks = NROWS*NCOLS;
  auto h_integralvsmaxamp = new TH2F("h_integralvsmaxamp",
				     Form("Number of Block hit per event for run %d;charge integral;maximum amplitude;",nrun),
				     100,0,900,
				     100,0,150);

  TFile *f = new TFile(Form("nps_eel108_%d.root", nrun)); 
  TTree *t = (TTree*) f->Get("T");
  Int_t NSampWaveForm;
  t->SetBranchAddress("Ndata.NPS.cal.fly.adcSampWaveform",&NSampWaveForm) ;
  Double_t SampWaveForm[nblocks*NROWS*NCOLS];//FIXME
  t->SetBranchAddress("NPS.cal.fly.adcSampWaveform",&SampWaveForm) ;
  
  Long64_t nentries = t->GetEntries();
 
  // --------------------------------------------------------------------------------
 const int event_display = 25;
 const int edis=nentries;
 TH1F *h_SampWaveForm[edis];
  
 TString eventnumberlabel[edis];

 for (Int_t n=0;n<edis;n++)
    {
      h_SampWaveForm[n] = new TH1F(Form("h_SampWaveForm_%d",n),"",sample_size,0,4*sample_size);
      h_SampWaveForm[n]->SetLineWidth(2);
      h_SampWaveForm[n]->SetLineColor(2);
      h_SampWaveForm[n]->GetXaxis()->SetTitle("Time");
      h_SampWaveForm[n] ->GetYaxis()->SetTitle("FADC voltage");
    }


  Int_t SampPMT;
  Int_t NSamp;
  Double_t SampADC;
  
  
  Int_t nevent=0;  
  cout<<" nentries="<<nentries<<endl;
  int count = 0;
  
  for (int i = istart ; i < nentries; i++)
    {
      if(nevent>=edis) i=nentries+1;
      Bool_t verbose = kTRUE;
      t->GetEntry(i);
      Int_t ns=0;
       
      while (ns < NSampWaveForm)
	{
	  SampPMT=SampWaveForm[ns++];
	  //  cout<<"SampPMT"<<SampPMT<<endl;
	  NSamp=SampWaveForm[ns++];
	  if(SampPMT==PMT && nevent<edis)
	    {
	      // if(verbose) cout<<"positive side event number =" << i << "NSamp="<<NSamp<<endl;
	      Double_t threshold_lo=5;  // pick these to focus on events in one of the 2 peaks in the ampl. histo
	      Double_t threshold_hi=99999;
	      Bool_t  pass_threshold=kFALSE;
	      Double_t amp=0;
	      for (Int_t n=0;n<NSamp;n++)
		{
		  SampADC=SampWaveForm[ns++];
		  h_SampWaveForm[nevent]->SetBinContent(n, SampADC);
		  amp+=SampADC;
		  
		  
		}
	      eventnumberlabel[nevent]=Form("Event %d (Th:%5.1f -- %5.1f)",i, threshold_lo, threshold_hi);
	      count = count+1;
	      if(amp >= 50)
		{
		    
		  if ((h_SampWaveForm[nevent]->GetMaximum()) > 60 && ( h_SampWaveForm[nevent]->GetMaximum())  < 140 )
		    {
		      cout <<"event_number = "<<  i  <<" Amplitude = " << h_SampWaveForm[nevent]->GetMaximum()<<endl;
		     }
	       h_integralvsmaxamp-> Fill(amp, h_SampWaveForm[nevent]->GetMaximum());
	       if( ( h_SampWaveForm[nevent]->GetMaximum() > threshold_lo ) &&
		  ( h_SampWaveForm[nevent]->GetMaximum() < threshold_hi )) 
		 
		 {
		  pass_threshold = kTRUE;
		 }
	      
	       if(pass_threshold)
		nevent+=1;
		}
	    }
	  else
	    ns+=NSamp; // skip to the next value of SampWaveForm  to find the number of the next PMT hit (if any)
	}
      
      
    }
  cout<< "count" << count<<endl;
   
   
  gStyle->SetPadBorderMode(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TCanvas *c1= new TCanvas();
  c1->SetTitle(Form("PMT_%d Positive side histogram",PMT));//put the values for cut for the reminder.
  c1->cd();
  TPad *padtitle = new TPad("padtitle", "padtitle",0.01,0.95,0.99,0.99);
  padtitle->Draw();
  padtitle->cd();
  padtitle->SetFillStyle(0);
  
  auto tex = new TLatex(0.5,0.5,Form("Run %d PMT %d",nrun,PMT));
  tex->SetTextAlign(22);
  tex->SetTextSize(0.5);
  tex->Draw();

  c1->cd();
  TPad *padmain = new TPad("padmain", "padmain",0.04,0.04,0.99,0.94);
  padmain->Draw();
  padmain->cd();
  padmain->SetFillStyle(0);
  
  
  int nca=sqrt( event_display);
  if((nca-sqrt( event_display))!=0)nca=sqrt( event_display)+1;
  padmain->Divide(nca,nca,0,0);

  Double_t maxdisplay=0;
  
  TLatex *tt[ event_display];
  TLatex *ttt=new TLatex();
  ttt->SetNDC();
  ttt->SetTextSize(0.1);
  for(int i=0;i<25;i++)
    {
      tt[i]=new TLatex();
      tt[i]->SetNDC();
      tt[i]->SetTextSize(0.1);
      
      if(h_SampWaveForm[i]->GetMaximum()> maxdisplay)
	maxdisplay= h_SampWaveForm[i]->GetMaximum();
    }
  for(int i= 0;i< event_display;i++)
    {
      padmain->cd(i+1);      
      gStyle->SetOptStat(0);
      h_SampWaveForm[i]->SetMaximum(maxdisplay*1.1);
 
      h_SampWaveForm[i]->Draw();
      tt[i]->DrawLatex(0.15, 0.8, eventnumberlabel[i]);
    }


  TString pdfname= Form("%d_run_%d_pmt__%d_istart.pdf",nrun,PMT,istart); 
  c1->Print(pdfname,"pdf");
  TCanvas *c2= new TCanvas(pdfname+" average amplitudes", pdfname+" average amplitudes");
  h_integralvsmaxamp->Draw("colz");

  // --------------------------------------------------------------------------------

 }
 

