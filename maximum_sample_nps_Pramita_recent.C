#define SHOWER_NROWS 36
#define SHOWER_NCOLS 7
#define SHOWER_NSAMPLES 100  // PMTindex + n_Samp + [Samples]
#define COLUMN_NROWS 6       // rows, cols for per-column canvas
#define COLUMN_NCOLS 4
const int nPedSamples = 4;    // number of samples used to evaluate pedestal
  

Int_t row(int PMT, int nrownumber )
{
  return  (PMT-1)%nrownumber+1;
}

Int_t column(int PMT, int nrownumber)
{
  return  int ((PMT-1) / nrownumber)+1;
}

tuple<int,int> testcolrow(int PMT)
{

  // local variable definition
  int resultrow=row(PMT,SHOWER_NROWS);
  //  cout<<"row ="<< resultrow<<endl;
  int resultcol=column(PMT,SHOWER_NROWS);
  //cout<<" column="<<resultcol<<endl;
  return {resultrow,resultcol} ;
}

// nblock_column == 1 -- SHOWER_NCOLS

void maximum_sample_nps_Pramita_recent( int nblock_column=-1, Int_t nrun=24 )
{

  // Check the argument(s)
  if((nblock_column != -1) && 
     ( (nblock_column < 1) || (nblock_column > SHOWER_NCOLS) ) ) {
    cout
      << Form("  First argument (nblock_column=%d) must be between 1 and %d; or -1 for all columns", nblock_column, SHOWER_NCOLS) 
      << endl;
    return;  // exit macro
  }
  
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetOptStat(11111);

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
  
  TString pdfname=Form("number_of_hits_%d.pdf",nrun);

  TFile *f = new TFile(Form("nps_eel108_%d_162515.root", nrun));
  TTree *t = (TTree*) f->Get("T");
  
  Int_t NSampWaveForm;
  t->SetBranchAddress("Ndata.NPS.cal.fly.adcSampWaveform",&NSampWaveForm) ;
  cout <<"NSampWaveform"<< NSampWaveForm<<endl;;
  Double_t SampWaveForm[(1+1+SHOWER_NSAMPLES)*SHOWER_NCOLS*SHOWER_NROWS]; // maximum possible size of waveform data array from ROOT file
                                                                          // PMT_num, nSamples, [ samples ... ]  => 1+1+SHOWER_NSAMPLES
  t->SetBranchAddress("NPS.cal.fly.adcSampWaveform",&SampWaveForm) ;  auto h= new TH2F("hit on the calo",
										       Form("Recorded hits for run %d;column number;row number",nrun),
										       SHOWER_NCOLS+2, -0.5, SHOWER_NCOLS+1.5,
										       SHOWER_NROWS+2, -0.5, SHOWER_NROWS+1.5);
  
  auto h_numberofPMTS = new TH1F("h_numberofPMTS",
                                 Form("Number of Block hit per event for run %d;number of PMTS hit;",nrun),
				 21,-0.5,20.5);
  auto h_mean_for_column_vs_PMT_number = new TH1F("h_mean_for_column_vs_column_number",
						  Form("h_mean_for_column vs column_number %d;column number;",nrun),
						  16,-0.5,15.5);
  TH1F *h_amp[SHOWER_NCOLS*SHOWER_NROWS];              // indexed by PMT number
  for(int i=0;i<SHOWER_NCOLS*SHOWER_NROWS;i++)
    {
      h_amp[i]=new TH1F(Form("h_amp_%d",i),Form("h_amp_%d run %d",i,nrun), 50, 50, 500);
    }
  
  auto h_average_pulse= new TH2F( "average_pulse_amp", "Average pulse amplitude in PMT",
				  SHOWER_NCOLS+2, -0.5, SHOWER_NCOLS+1.5,
				  SHOWER_NROWS+2, -0.5, SHOWER_NROWS+1.5);
 
  Int_t PMT_number;
  Int_t NSamp;
  Double_t SampADC;
  Double_t pedVal = 0; // pedestal value
  
  Long64_t nentries = t->GetEntries();
  int inp;
  nentries=83000;
  for (int i = 0; i < nentries; i++)
    {
      Bool_t verbose = kTRUE;
      if((i%10000)==0) cout<<i<<":"<<nentries<<endl;
      t->GetEntry(i);
      Int_t ns = 0;
      for (Int_t n=0; n<nPedSamples; n++) { // loop over 1st 'nPedSamples' samples to get pedestal data for 'SampPMT'
        SampADC=SampWaveForm[ns++];
	pedVal += SampADC;
      }
      pedVal = pedVal/nPedSamples;
      ns = ns - nPedSamples; // reset sample index 'ns' back to start of this hit's waveform array

      while (ns < NSampWaveForm)
        {
	  
	  PMT_number=SampWaveForm[ns++]; // also the pmt number (or row) starts at 1 not 0
	  NSamp=SampWaveForm[ns++];  // number of samples that follow

	  Int_t column_n = column(PMT_number,SHOWER_NROWS)-1;
	  Int_t row_n    =    row(PMT_number,SHOWER_NROWS)-1;
	 
	  if(column_n > SHOWER_NCOLS-1) {
	    ns+=NSamp;
	    continue;   // not a column we care about; skip to next block
	  }
	 
	  Double_t amp=0;
	  Double_t max=0;
	  for (Int_t n=0;n<NSamp;n++)
	    {
	      SampADC=SampWaveForm[ns++];//  extract a single sample value
	      amp+= SampADC-pedVal;
	      //if(SampADC>max)max=SampADC; //better to have integral amplitude more sensitive to noise
	    }
	  if(amp>=5)
	    
	    h_amp[PMT_number]->Fill(amp);
	} // end of loop over waveform sample array (all hits in this event)
    } // end of loop over events
  

  // creating a 2D histo to display all the means at once
  // also creating 14 1D histo to display the mean of one column at once
  TH1F *h_per_column[SHOWER_NCOLS];
  TString name_histo;
  for (int i=0;i<SHOWER_NCOLS;i++)
    {
      name_histo=Form("per colum %d",i+1);
      h_per_column[i]=new TH1F(name_histo, name_histo, 36, -0.5, 35.5);
    }
  
 for(int irow=0;irow<SHOWER_NROWS;irow++)
    {
      Double_t sum = 0;
      
      for(int icol=0;icol<SHOWER_NCOLS;icol++)
      {
	int ipmt=irow+icol*SHOWER_NROWS;
      
	h_per_column[icol]->Fill(irow+1, h_amp[ipmt]->GetEntries());

	h_average_pulse->Fill(icol+1,irow+1,h_amp[ipmt]->GetMean());
	
	sum = sum+h_amp[ipmt]->GetMean();	
      }
      cout<<" irow = " <<irow<< "average ="<<sum/SHOWER_NCOLS<<endl;
      h_mean_for_column_vs_PMT_number->Fill(irow,sum/SHOWER_NCOLS);
    }
  



 //h_amp is the list(14 columns*16 rows= 224) of 1 dimensional histogram and GetMean gives the list of means of each list.
 TCanvas *c1= new TCanvas(pdfname+" average amplitudes", pdfname+" average amplitudes");
 // h_average_pulse->SetStats(0);
 gStyle->SetPaintTextFormat("4.1f");
 h_average_pulse->Draw("colz text");
 h_average_pulse->SetMinimum(0);
 h_average_pulse->SetFillColor(kBlue-7);
 h_average_pulse->SetStats(0);
 c1->Print(pdfname+"(");
 // Open PDF file and leave open

 // Plot amplitudes in requested column; or plot all amplitudes if nblock_column==-1
 //  - generalize the plotting code below so it can work in a for() loop
 //    - need to be able to make new canvas in the loop (if needed)
 //  - then wrap the whole thing in a for() loop that loops over all columns if nblock_column=-1,
 //    OR just does a single loop from nblock_column to nblock_column if nblock_column>-1

TCanvas *c2= new TCanvas(pdfname+" average amplitudes 2", pdfname+" average amplitudes 2");
 c2->Divide(2,4);
 for(int icol=0;icol<SHOWER_NCOLS;icol++)
   {
     c2->cd(icol+1);
     h_per_column[icol]->Draw("HIST");
     h_per_column[icol]->SetMinimum(0);
     h_per_column[icol]->SetFillColor(kBlue-7);
     h_per_column[icol]->SetStats(0);

   }
 c2->Print(pdfname);


 
  int nStart=0;
  int nEnd=SHOWER_NCOLS;

  if(nblock_column > -1) {
    nStart=nblock_column-1;
    nEnd=nblock_column;
  }

  TCanvas *colCanvas[nEnd - nStart];

  for(int col=nStart; col<nEnd; col++) {
    colCanvas[col]=new TCanvas(
			       Form("Column %d",col),
			       Form("Amplitudes in Column %d (run %d)",col+1, nrun)
      );
    colCanvas[col]->Divide(6,6);
    
    int subplot=1;
    int PMT_start=(col)*SHOWER_NROWS;   // first PMT number in the column we are plotting
 
    for(int nPMT=PMT_start; nPMT<PMT_start+SHOWER_NROWS; nPMT++) {
      if(col==0)cout<<"col="<<col<<" subplot="<<subplot<<" nPMT=" << nPMT<<endl;
      colCanvas[col]->cd(subplot++);
      
      h_amp[nPMT]->Draw();
      h_amp[nPMT]->SetMinimum(0);
     
      h_amp[nPMT]->SetFillColor(kBlue-7);
      h_amp[nPMT]->SetStats(0);
	}
    if(col<nEnd-2)
      colCanvas[col]->Print(pdfname); // print to (previously opened) PDF file
    else
	colCanvas[col]->Print(pdfname+")"); // print to (previously opened) PDF file
      
  }


 
}    
