//This scripts creates bunch of events(In this case it's 25 events) for PMT  for positive side PMTs.
void specific_Pmt( Int_t nrun=14464, Int_t PMT=10, char side ='P', Int_t istart=500)
// PMT can be between 1 and 14
// side can be "P" or "N"
// i start is the event number at which we start loking for a signal  
//, String_t side="P")
{
  // gStyle->SetOptTitle(0);                                                         
  //gStyle->SetOptStat(0);                                                            
  //gStyle->SetPalette(1);                                                            

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameFillColor(0);

  gStyle->SetPadTopMargin(.0001);                                                    
  //gStyle->SetPadLeftMargin(.01);
  // gStyle->SetTitle("PMT number");
  gStyle->SetPadRightMargin(.01);                                                  
  //gStyle->SetPadBottomMargin(.01);                                                 
  
  gStyle->SetTitleOffset(1.1, "X");                                                
  gStyle->SetTitleOffset(1.08, "Y");                                               
  gStyle->SetTitleFont(42,"X");                                                    
  gStyle->SetTitleFont(42,"Y");                                                    
  gStyle->SetTitleSize(0.055,"X");                                                 
  gStyle->SetTitleSize(0.045,"Y");                                                 

  gStyle->SetLabelOffset(0.01, "X");                                               
  gStyle->SetLabelOffset(0.01, "Y");                                               
  gStyle->SetLabelFont(42,"X");                                                    
  gStyle->SetLabelFont(22,"Y");                                                    
  gStyle->SetLabelSize(0.05,"X");                                                 
  gStyle->SetLabelSize(0.05,"Y");                                                 
  
  //gStyle->SetNdivisions(101,"X");                                                  
  //gStyle->SetNdivisions(101,"Y");                                                  
  
  gStyle->SetStripDecimals(kFALSE);                                                
 



// this example is going to be for both positive and negative pmts only
{

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
  const int edis=25;
  TH1F *h_SampWaveForm[edis];

  for (Int_t n=0;n<edis;n++)
    {
      h_SampWaveForm[n] = new TH1F(Form("h_SampWaveForm_%d",n),"",175,0,175*4);
      h_SampWaveForm[n]->SetLineWidth(2);
      h_SampWaveForm[n]->SetLineColor(2);
      h_SampWaveForm[n]->GetXaxis()->SetTitle("Time");
      h_SampWaveForm[n] ->GetYaxis()->SetTitle("Voltage");
    }


  Int_t SampPMT;
  Int_t NSamp;
  Double_t SampADC;

  Long64_t nentries = t->GetEntries();
 
  // --------------------------------------------------------------------------------
  Int_t nevent=0;  
  cout<<" nentries="<<nentries<<endl;
  int inp;

  for (int i = istart ; i < nentries; i++)
    {
      if(nevent>=edis) i=nentries+1;
      Bool_t verbose = kTRUE;
      if((i%1000)==0) cout<<i<<":"<<nentries<<endl;
      
      t->GetEntry(i);
      
      Int_t nsp=0;
      Int_t nsn=0;
      switch (side)
	{
	case 'P':
	  
	  //	   Int_t nsp = 0;
	  while (nsp < NSampWaveFormP)
	    {
	      SampPMT=SampWaveFormP[nsp++];
	      NSamp=SampWaveFormP[nsp++];
	      if(SampPMT==PMT && nevent<edis)
		{
		  // if(verbose) cout<<"positive side event number =" << i << "NSamp="<<NSamp<<endl;
		  for (Int_t n=0;n<NSamp;n++)
		    {
		      SampADC=SampWaveFormP[nsp++];
		      h_SampWaveForm[nevent]->SetBinContent(n, SampADC);
		    }
		  h_SampWaveForm[nevent]->SetTitle(Form("Event %d",i));
		  h_SampWaveForm[nevent]->SetTitleSize(0.90);
		  nevent+=1;
		}
	      else
		nsp+=NSamp; // skip to the next value of SampWaveForm  to find the number of the next PMT hit (if any)
	    }
	  break;
	case 'N':
	  //Int_t nsn = 0;
          while (nsn < NSampWaveFormN)
            {
              SampPMT=SampWaveFormN[nsn++];
              NSamp=SampWaveFormN[nsn++];
              if(SampPMT==PMT && nevent<edis)
                {
		  //  if(verbose) cout<<"negative side event number ="<<i<< "NSamp="<<NSamp<<endl;                                                               
                  for (Int_t n=0;n<NSamp;n++)
                    {
                      SampADC=SampWaveFormN[nsn++];
                      h_SampWaveForm[nevent]->SetBinContent(n, SampADC);
                    }
		  //  h_SampWaveForm[nevent]->SetTitle(Form("Event %d",i));
		  //h_SampWaveForm[nevent
		  nevent+=1;
                }
              else
                nsn+=NSamp; // skip to the next value of SampWaveForm  to find the number of the next PMT hit (if any)                                                    
            }
	  break;
	}
    }
      
  TCanvas *c1= new TCanvas();
  c1->SetTitle(Form("PMT_%d Positive side historgam",PMT));
	       // c1->SetTitle("(Form("PMT_%d",PMT))", "Positive side histogram");
  
  int nca=sqrt(edis);
  if((nca-sqrt(edis))!=0)nca=sqrt(edis)+1;
  c1->Divide(nca,nca);
  //h_SampWaveForm->SetMinimum(0);
  // h_SampWaveForm->SetFillColor(kBlue-7);
  // h_SampWaveForm->SetStats(1);
  // gStyle->SetOptStat(111111);

  Double_t maxdisplay=0;
  for(int i= 0;i<edis;i++)
    {
      //   cout<< "i="<<i<<" maxdisaply="<<maxdisplay<<" max histo="<<h_SampWaveForm[i]->GetMaximum()<<endl;
      if(h_SampWaveForm[i]->GetMaximum()> maxdisplay) maxdisplay= h_SampWaveForm[i]->GetMaximum();
    }
      

  TString pdfname= Form("%d_run %d_pmt %c_side %d_istart.pdf",nrun,PMT,side,istart); 
 
  for(int i= 0;i<edis;i++)
    {
      c1->cd(i+1);
      //h_SampWaveForm->SetMinimum(0);
      // h_SampWaveform->SetFillColor(kBlue-7);
      //h_SampWaveForm->SetStats(1);
      //gStyle->SetOptStat(111111);
      h_SampWaveForm[i]->SetMaximum(maxdisplay*1.1);
      h_SampWaveForm[i]->Draw();
    }
  c1->Print(pdfname,"pdf");
      //     h_SampWaveForm[i]->SetMaximum(40); //* H_SampWaveForm[i]->GetMaximum());     
    // cout<<i<<"->"<<  h_SampWaveForm[i]->GetMaximum()<<endl;
    



    
  // --------------------------------------------------------------------------------
    
 }
}
