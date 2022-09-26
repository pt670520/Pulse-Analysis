//This scripts creates bunch of events(In this case it's 25 events) for PMT 8 for positive side PMTs.
void specific_PMT_new( Int_t nrun=14464, Int_t PMT=8, char side ='P', Int_t istart=500)
// PMT can be between 1 and 14
// side can be "P" or "N"
// i start is the event number at which we start loking for a signal  
//, String_t side="P")
{
  // // gStyle->SetOptTitle(0);                                                         
  // //gStyle->SetOptStat(1111);                                                            
  // //gStyle->SetPalette(1);                                                            

  // gStyle->SetCanvasColor(0);
  // gStyle->SetFrameFillColor(0);
  // gStyle->SetPadTopMargin(.0001);                                               
  // gStyle->SetPadLeftMargin(.01);
  // gStyle->SetTitle("PMT number");
  // gStyle->SetPadRightMargin(.01);
  // gStyle->SetPadBottomMargin(.01);                                             
  
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
                                                 
  
  // gStyle->SetStripDecimals(kFALSE);                                                
  

// this example is going to be for both positive and negative pmts only


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
  TString eventnumberlabel[edis];

  for (Int_t n=0;n<edis;n++)
    {
      h_SampWaveForm[n] = new TH1F(Form("h_SampWaveForm_%d",n),"",175,0,175*4);
      h_SampWaveForm[n]->SetLineWidth(2);
      h_SampWaveForm[n]->SetLineColor(2);
      h_SampWaveForm[n]->GetXaxis()->SetTitle("Time");
      h_SampWaveForm[n] ->GetYaxis()->SetTitle("FADC voltage");
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
		  eventnumberlabel[nevent]=Form("Event %d",i);
		  //		  	  h_SampWaveForm[nevent]->SetTitle(Form("Event %d",i));
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
                  h_SampWaveForm[nevent]->SetTitle(Form("Event %d",i));
                  nevent+=1;
                }
              else
                nsn+=NSamp; // skip to the next value of SampWaveForm  to find the number of the next PMT hit (if any)                                                    
            }
	  break;
	}
    }

  
  gStyle->SetPadBorderMode(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TCanvas *c1= new TCanvas();
  c1->SetTitle(Form("PMT_%d Positive side histogram",PMT));
  c1->cd();
  TPad *padtitle = new TPad("padtitle", "padtitle",0.01,0.95,0.99,0.99);
  padtitle->Draw();
  padtitle->cd();
  padtitle->SetFillStyle(0);

  auto tex = new TLatex(0.5,0.5,Form("Run %d PMT %d-%s",nrun,PMT,&side));
  tex->SetTextAlign(22);
  tex->SetTextSize(0.5);
  tex->Draw();

  c1->cd();
  TPad *padmain = new TPad("padmain", "padmain",0.04,0.04,0.99,0.94);
  padmain->Draw();
  padmain->cd();
  padmain->SetFillStyle(0);
  
  
  int nca=sqrt(edis);
  if((nca-sqrt(edis))!=0)nca=sqrt(edis)+1;
  padmain->Divide(nca,nca,0,0);

  Double_t maxdisplay=0;
 
  TLatex *tt[edis];
  TLatex *ttt=new TLatex();
  ttt->SetNDC();
  ttt->SetTextSize(0.1);
  for(int i=0;i<edis;i++)
    {
      tt[i]=new TLatex();
      tt[i]->SetNDC();
      tt[i]->SetTextSize(0.1);

      if(h_SampWaveForm[i]->GetMaximum()> maxdisplay)
	maxdisplay= h_SampWaveForm[i]->GetMaximum();
    }
     
  for(int i= 0;i<edis;i++)
    {
      padmain->cd(i+1);      
      gStyle->SetOptStat(0);
      h_SampWaveForm[i]->SetMaximum(maxdisplay*1.1);
      h_SampWaveForm[i]->Draw();
      tt[i]->DrawLatex(0.3,  0.8,  eventnumberlabel[i]);
    }


  TString pdfname= Form("%d_run_%d_pmt_%c_side_%d_istart.pdf",nrun,PMT,side,istart); 
  c1->Print(pdfname,"pdf");

     
  // --------------------------------------------------------------------------------
    
 }
 

