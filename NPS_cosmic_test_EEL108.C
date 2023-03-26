//This programs outputs the good time corresponding to maximum amplitiude and average energy deposited in the crystals. Additionally, the relative energy resolution of the crystals is also determined. It also has a text file using average charge integral to perform gain matching.

void NPS_cosmic_test_EEL108( Int_t nrun=55, Int_t nruns = 56, Double_t threshold=5.0, Bool_t plotAllHits=kFALSE) {
   
  const int maxSamp = 1+1+100;    // Size of sample array for single hit: Ch# + NSamp + <samples> (run 22,  PTW=404 (100(+1) samples)
  const int nrows = 16;
  const int ncols = 15;         // 15 FADCs/Crate in NPS setup
  const int nPedSamples = 4;    // To get average pedestal for 4 smaples.
  const int NPS_rows = 36;
  const int NPS_cols = 30;
  // number of samples used to evaluate pedestal
  
  int nblocks = nrows*ncols;
  Bool_t loopOverAllEvents = kFALSE;
  
  
  //########### Tchain adds two rootfiles#########################3
  
   TChain *t= new TChain("T");   // No need to get root file while using tchain. acts as a pointer.
   t->Add(Form("nps_eel108_%d.root ", nrun));
   t->Add(Form("nps_eel108_%d.root ", nruns));
   //####################################       //ch.Merge("all.root");
   
   Int_t NSampWaveForm;
   t->SetBranchAddress("Ndata.NPS.cal.fly.adcSampWaveform",&NSampWaveForm) ;  // Total size of the FADC sample array (all hits)
   Double_t SampWaveForm[91502];//830000];91502 = maxSamp * nblocks];// maxSamp * nblocks//// FIXME: Allocating max possible size here (based on constants in this script). Really == NSampWaveForm; dynamically allocate instead?
   t->SetBranchAddress("NPS.cal.fly.adcSampWaveform",&SampWaveForm) ;//&SampWaveForm: pointer to the variable and copy of data from root leaf to SampWaveForm. (contains the waveforms for all the blocks)
  TH1F *h_SampWaveForm[nblocks];
  TH1F *h_SampWaveFormHitAve[nblocks];
  TH1F *h_SampWaveFormHitMax[nblocks];
  auto h_Cal_Time_PMT = new TH2F(Form("h_Cal_Time%d", nrun), "Cal Time",
				 110,0,110,  // PMT count: Range should be equal to bin count if filling with integers
				 maxSamp,0,maxSamp*4);  // put y-axis in nanosecs
  
    auto *h_time_block_number = new TH2F("h_time_block_number ","average time ",20,0,110
					 ,100,95,130);
    TH1F *h_amp[nblocks];  //Defining array of histogram over whole run
    TH1F *h_time[nblocks];
    Int_t cellHitCount[nblocks];
    Double_t fit_means[nblocks]; // will hold the mean values from the fits
    Double_t fit_mean[nblocks];
    Int_t SampPMT_array[nblocks];
    for (Int_t n=0;n<nblocks;n++)
      {
	cellHitCount[n] = 0;
	fit_means[n] = 0;
      h_SampWaveForm[n] = new TH1F(Form("h_SampWaveForm_%d",n),"",maxSamp,0,maxSamp*4); // *4 to put x-axis in nanosecs
      h_SampWaveForm[n]->SetLineWidth(1);
      h_SampWaveForm[n]->SetLineColor(2);
      //h_SampWaveForm[n]->GetYaxis()->SetRangeUser(-1,1);
      h_SampWaveFormHitAve[n] = new TH1F(Form("h_SampWaveFormHitAve_%d",n),"",maxSamp,0,maxSamp*4); // *4 to put x-axis in nanosecs
      h_SampWaveFormHitAve[n]->SetLineWidth(1);
      h_SampWaveFormHitAve[n]->SetLineColor(4);//line color 4 is blue, 2 is red. 
      h_SampWaveFormHitMax[n] = new TH1F(Form("h_SampWaveFormHitMax_%d",n),"",maxSamp,0,maxSamp*4); // *4 to put x-axis in nanosecs
      h_amp[n]=new TH1F(Form("h_amp_%d",n) ,Form("h_amp_%d run %d",n,nrun), 50, 0, 500);
      h_time[n] = new TH1F("h_time ","average time ",50,0,250);//FIXME? What time should i consider?
    }
    
    Int_t SampPMT;
    Int_t NSamp;
    Double_t SampADC;
  
    Long64_t nentries = t->GetEntries();
    // cout << "nentries" << nentries<<endl;
    TCanvas *CanSamp= new TCanvas("canSamp", 
				  Form("NPS Calorimeter FADC Sample (EEL108) -- Run %d", nrun), 
				  3200,1500);  // 4k monitor size
    TCanvas *CanSampHitAve= new TCanvas("canSampHitAve", 
					Form("NPS Calorimeter FADC Averaged Hits over %4.1f mV threshold (EEL108) -- Run %d", threshold, nrun),
					3200,1500);  // 4k monitor size
    // --------------------------------------------------------------------------------
    Int_t cosmicEventCount = 0;
  Long64_t evNum;
  //nentries = 83000;
  for (evNum = 0; evNum < nentries; evNum++)//this a for loop that loops over all events in the root tree
    {  
      t->GetEntry(evNum);
      if (evNum%1000==0) cerr << endl << "---> Event = " << evNum << endl;
      if (evNum%10==0)   cerr << "." ;
      Int_t cellHitMap[nblocks];
      for (Int_t n=0;n<nblocks;n++)
	{
	  cellHitMap[n]=0;
	  for (Int_t nn=0;nn<maxSamp;nn++)  // recall that the ROOT histogram array is also #bins + underflow + overflow, which happens to equal maxSamp == 1 + 1 + nsamples.
	    { 
	      h_SampWaveForm[n]->SetBinContent(nn,0); //zeroing out the histogram incase it has the older memory. Its mportant to re initialize.
	    }
	}
      
      Int_t ns = 0;
      Bool_t cosmicHit=kFALSE; 
      Double_t amp[3*NPS_rows]={0}; /// 3* coz only summed from three columns with  NPS_rows.
      Double_t time[3*NPS_rows]={0}; /// 3* coz only summed from three columns with  NPS_rows.
      //cout << "NSampWaveform" << NSampWaveForm <<endl;
      while (ns < NSampWaveForm) //one gigantic long array with all hits for a single events.
	{ 
	  SampPMT=SampWaveForm[ns++];
	  //SampPMT_array[ns] = SampPMT;
	  NSamp = SampWaveForm[ns++];    // NSamp = sample size.
	  Double_t sampMin= -9999;
	  Double_t sampMax=9999;
	  Bool_t cellHit=kFALSE;
	  
	  // ##################### Block removes Background Noise #############################################
	  Double_t pedVal = 0; // pedestal value
	  for (Int_t n=0; n<nPedSamples; n++)  // loop over 1st 'nPedSamples' samples to get pedestal data for 'SampPMT'
	    {
	      SampADC=SampWaveForm[ns++];
	      pedVal += SampADC;
	    }
	  pedVal = pedVal/nPedSamples;         // This gets the average pedestal
	  ns = ns - nPedSamples; // reset sample index 'ns' back to start of this hit's waveform array
	  
	  //###########################ENDS HERE ##############################################################
	  
	  // #######################Block calculates the Pulse/Charge Integral #####################
	  //if ()
	  
	  for (Int_t n=0;n<NSamp;n++)  // loop over waveform data for 'SampPMT'
	    {
	      SampADC=SampWaveForm[ns++];
	      h_SampWaveForm[SampPMT]->SetBinContent(n, SampADC);  // should probably be n+1 here
	      
	      amp[SampPMT]+=SampADC-pedVal;// subtract pedVal from SampADC ?
	      
	      
	      
	      if(SampADC > sampMax) sampMax = SampADC;
	      if(SampADC < sampMin) sampMin = SampADC;
	      
	      
	      // Declare a cell to be 'hit' if sample is over threshold
	      if( (SampADC > threshold) && (SampPMT < 108) ) {
		cellHit = kTRUE; // cells above 108 are not hooked up yet
	      }
	      
	    }
	  
	  // If we found a hit in cell 
	  if( cellHit ) {
	    cellHitCount[SampPMT]++; // counts how many hits were seen in eac cell index over the whole run
	    cellHitMap[SampPMT]=1;   // contains which cells were hit in this event
	    h_SampWaveFormHitAve[SampPMT]->Add( h_SampWaveForm[SampPMT] );

	    Double_t hit_time = 4* h_SampWaveForm[SampPMT]->GetMaximumBin();  // 4 *returns bin number for peak/biggest pulse
	    // cout<< "hit_time" << hit_time<<endl;
	    
	    time[SampPMT]=hit_time;
	    
	    cout<< "hit time or bin number corresponding to pulse  maximum  = " << time[SampPMT]<< "PMT number =" << SampPMT<< endl;
	    
	    h_time_block_number-> Fill(SampPMT, time[SampPMT]);//hit_time );
	    
	    if( h_SampWaveForm[SampPMT]->GetMaximum() > h_SampWaveFormHitMax[SampPMT]->GetMaximum() ) {
	      h_SampWaveFormHitMax[SampPMT]->Delete();
	      h_SampWaveFormHitMax[SampPMT] = (TH1F*) h_SampWaveForm[SampPMT]->Clone();
	      
	      
	    }
	    
	  }
	  
	} // end of loop over waveform blob/array
      
      for (Int_t j= 0; j< nblocks/NPS_rows; j++) {
	if ((cellHitMap[NPS_rows*j+1] == 1) && (cellHitMap[NPS_rows*(j+1)-1] == 1)) { //+1 in first cellindex as cell index 36 is broken in run 24.
	  for (Int_t jPMT = NPS_rows*j+1; jPMT<NPS_rows*(j+1);jPMT++)  {
	    h_amp[jPMT]->Fill(amp[jPMT]);    // integrated amplitude histo
	    if (time[jPMT] < 10) continue;
	    h_time[jPMT]->Fill(time[jPMT]);  // time histo
	   
	  }
	  
	  if (
	      ((cellHitMap [0+1] == 1) && (cellHitMap[ 35] == 1))  ||  // column 0
	      ((cellHitMap[36+1] == 1) && (cellHitMap[ 71] == 1))  ||  // column 1
	      ((cellHitMap[72+1] == 1) && (cellHitMap[107] == 1))      // column 2
	      ) 
	    
	    // h_amp[j]->Fill(amp[j]);
	    cosmicHit = kTRUE;
	  cout<< endl << " cosmic in column " <<j<< endl;
	  
	}
	

	
	
      }
    
      // See if we have column hits''
      //column 0-PMT 0-35
      //column 1-PMT 36-71
      //column 2- PMT 72-107
      //You have a cosmic that you can use if the top and bottom block on that column give signal above threshold i.e., if PMT is in colum 0, block 0 and block 35 are above threshold, Then histogram the integral of that PMT)
      
      if(cosmicHit) cosmicEventCount++;
      
      
    // Skip plotting waveforms unless
    //   - we have asked to plot all hits
    //     - OR -
    //   - we found a hit
      if( (! cosmicHit ) && (! plotAllHits) ) continue;  // process next event
      if (loopOverAllEvents)                  continue;  // process next event
    
      // Update per-event histo Canvas
      CanSamp->Clear();
      CanSamp->Update();
      
      CanSamp->Divide(ncols, nrows+1, 0.005, 0.0001);//,
      
      CanSamp->cd(1);
      TLatex* tex = new TLatex(0.1,0.6,"EEL108"); //defines the x and y value for the text in subcanvas
      tex->SetNDC(1); // NDC sets like the coordinate system
      tex->SetTextFont(42);
      tex->SetTextColor(1);
      tex->SetTextSize(0.3);
      tex->Draw();
      
      tex = new TLatex(0.1,0.1,"NPS Crate 1");
      tex->SetNDC(1);
      tex->SetTextFont(42);
      tex->SetTextColor(1);
      tex->SetTextSize(0.3);
      tex->Draw();
      
      CanSamp->cd(2);
    
      tex = new TLatex(0.1,0.7,Form("Event %lld",evNum));
      tex->SetNDC(1);
      tex->SetTextFont(42);
      tex->SetTextColor(2);
      tex->SetTextSize(0.3);
      tex->Draw();

      tex = new TLatex(0.1,0.4,"Amplitude (mV)");
      tex->SetNDC(1);
      tex->SetTextFont(42);
      tex->SetTextColor(1);
      tex->SetTextSize(0.3);
      tex->Draw();
      
      tex = new TLatex(0.1,0.1,"Time (ns)");
      tex->SetNDC(1);
      tex->SetTextFont(42);
      tex->SetTextColor(1);
      tex->SetTextSize(0.3);
      tex->Draw();
      cerr << "got here 3" << endl;
      for (Int_t n=0;n<108;n++) {
      Int_t canv_ind = ncols*(ncols+1-n%nrows) + (n/nrows) + 1;
      CanSamp->cd(canv_ind);
      if(canv_ind < (nrows-1+1)*ncols + 1) { gPad->SetBottomMargin(0); } else { gPad->SetBottomMargin(0.25); }
      Int_t n_color = 1;
      if(n > 35 && n < 72) n_color = 2;
      if(n > 71 && n < 108) n_color = 3; //hv is off for colum 3
      
       h_SampWaveForm[n]->SetLineColor(n_color);
       h_SampWaveForm[n]->Draw();
      }

      CanSamp->Update();
      
      Int_t inp=1;
      cout << endl
	   << Form("---> Processed event %lld.  Enter", evNum) << endl
	   << "           1 for next event," << endl
	   << "           N to jump to event 'N'" << endl
	   << "          -1 to loop over all remaining events " << endl
	   << "          -N to loop over next 'N' events" << endl
	   << "           0 to quit" << endl; 
      cout << "   ? : ";
      cin >> inp;
      cout << endl;

      if (inp ==  0) break; // exit event loop
      if (inp ==  1) continue; // next-event
      if (inp == -1) {  // loop over all events in file
	loopOverAllEvents=kTRUE;
	continue;
      }
      if (inp < -1) {  // loop over next '-inp' events (starting at this event 'evNum')
	loopOverAllEvents=kTRUE;
	inp = -inp + evNum;
	Long64_t tmp = nentries;
	nentries = nentries < inp ? nentries : inp;  // set to smaller of {nentries, inp}
	cout << Form("---> Looping until event: %lld  (of %lld)", nentries, tmp) << endl;
	continue;
      }
      if (inp > 1) {  // make next event to be loaded to 'inp'
	cout << Form("---> Loading event: %d", inp) << endl;
	evNum=inp-1;
	continue;
      }
      
    }
  cout << endl << Form("---> Processed %lld events", evNum) << endl;


  { // Plot "Average Hit" Canvas
    CanSampHitAve->Divide(ncols, nrows+1, 0.005, 0.0001);
    
    CanSampHitAve->cd(1);
    TLatex* tex = new TLatex(0.1,0.6,"EEL108");
    tex->SetNDC(1); tex->SetTextFont(42); tex->SetTextColor(1); tex->SetTextSize(0.3);
    tex->Draw();

    tex = new TLatex(0.1,0.2,"NPS Crate 1");
    tex->SetNDC(1); tex->SetTextFont(42); tex->SetTextColor(1); tex->SetTextSize(0.3);
    tex->Draw();
    
    CanSampHitAve->cd(2);
    tex = new TLatex(0.1,0.6,"Amplitude (mV)");
    tex->SetNDC(1); tex->SetTextFont(42); tex->SetTextColor(1); tex->SetTextSize(0.3);
    tex->Draw();
    
    tex = new TLatex(0.1,0.2,"Time (ns)");
    tex->SetNDC(1); tex->SetTextFont(42); tex->SetTextColor(1); tex->SetTextSize(0.3);
    tex->Draw();
    
    CanSampHitAve->cd(3);
    tex = new TLatex(0.1,0.6,Form("Events Processed: "));
    tex->SetNDC(1); tex->SetTextFont(42); tex->SetTextColor(2); tex->SetTextSize(0.3);
    tex->Draw();
    
    tex = new TLatex(0.1,0.2,Form("  %lld", evNum));
    tex->SetNDC(1); tex->SetTextFont(42); tex->SetTextColor(1); tex->SetTextSize(0.3);
    tex->Draw();

    CanSampHitAve->cd(4);
    tex = new TLatex(0.1,0.6,Form("Cosmic Hits:"));
    tex->SetNDC(1); tex->SetTextFont(42); tex->SetTextColor(2); tex->SetTextSize(0.3);
    tex->Draw();

    tex = new TLatex(0.1,0.2,Form("  %d", cosmicEventCount));
    tex->SetNDC(1); tex->SetTextFont(42); tex->SetTextColor(1); tex->SetTextSize(0.3);
    tex->Draw();
    
    for (Int_t n=0;n<108;n++) {
      Int_t canv_ind = ncols*(ncols+1-n%nrows) + (n/nrows) + 1;
      CanSampHitAve->cd(canv_ind);
      if(canv_ind < (nrows-1+1)*ncols + 1) { gPad->SetBottomMargin(0); } else { gPad->SetBottomMargin(0.25); }//tweaks the vertical packing of histogram and squish it. 
      
      if(cellHitCount[n] > 0) {
        h_SampWaveFormHitAve[n]->Scale(1./cellHitCount[n], "nosw2");//"" doesn't change the weights/error values in the histo
      }
      //h_SampWaveFormHitAve[n]->GetYaxis()->SetRangeUser(-10,40);
      h_SampWaveFormHitAve[n]->Draw();
      tex = new TLatex(0.3,0.6,Form("%d",(Int_t) cellHitCount[n]));
      tex->SetNDC(1); tex->SetTextFont(42); tex->SetTextColor(1); tex->SetTextSize(0.2);
      tex->Draw();

      //h_SampWaveFormHitMax[n]->GetYaxis()->SetRangeUser(-1,1);
      h_SampWaveFormHitMax[n]->SetLineWidth(1);
      h_SampWaveFormHitMax[n]->SetLineColor(2);
      h_SampWaveFormHitMax[n]->SetLineStyle(3);
      h_SampWaveFormHitMax[n]->SetLineStyle(3);
      h_SampWaveFormHitMax[n]->Draw("same");  // FIXME; roughed in here -- need some logic for correct y-axis size

    }

    CanSampHitAve->Update();
    
  }
  Double_t rms_amplitude[3*NPS_rows];
  Double_t mean_good_pulse[3*NPS_rows];
  Double_t good_amp[3*NPS_rows];
  Double_t  rms_time[3*NPS_rows];
  Double_t  relative_resolution[3*NPS_rows];
  Double_t ey[3*NPS_rows];
  Double_t ez[3*NPS_rows];
   
  TCanvas *C10= new TCanvas("c10","Pulse integral vs count",3200,1500);
  C10->Divide(ncols,nrows, 0, 0);
  for(Int_t n=0;n<3*NPS_rows;n++)

    {
     Int_t canv_ind = ncols*(ncols-n%nrows) + (n/nrows) + 1;
     C10->cd(canv_ind);//changes focus to a particular sub-canvas
     Int_t n_color = 4;
     if(n > (NPS_rows-1) && n < (2*NPS_rows)) n_color = 1;
     if(n >(2*NPS_rows)-1 && n < (3*NPS_rows)) n_color = 3;//hv is off for colum 3
       
     h_amp[n]->SetLineColor(n_color);
      
     h_amp[n]->SetMinimum(0);
   
     h_amp[n]->Draw();
     cout<< " PMT number ="<< n  <<"number of entries"<<  h_amp[n]->GetEntries() <<endl;
     TF1  *f1 = new TF1(Form("f1_%d",n),"gaus+ expo(3)"); // 5 params ==> 3 for gaus(p0, p1, p2), 2 for exp(p3 + p4*x)
     f1->SetParameters(h_amp[n]->GetMaximum(), h_amp[n]->GetMean(), h_amp[n]->GetRMS(), -10, -10 );  // exp seeds reflect no background (small); will need to fix for data with an exp background
     f1->FixParameter(3, -10);  // fix exp params to keep them from screwing up the fit for now
     f1->FixParameter(4, -10);
     
     cout<<"   Maximum = "<<h_amp[n]->GetMaximum() <<"  mean = "<< h_amp[n]->GetMean() <<"   RMS =  "<<h_amp[n]->GetRMS() <<endl;
     good_amp[n] =f1->GetParameter(1);
     
     h_amp[n]->Fit(Form("f1_%d",n),"LR","",  h_amp[n]->GetMean()-( h_amp[n]->GetRMS()),   h_amp[n]->GetMean()+   (h_amp[n]->GetRMS()));//L” Use log likelihood method (default is chi-square method). To be used when the histogram represents counts
     //“R” Use the range specified in the function range
     
     rms_amplitude[n] =f1->GetParameter(2);
     
     ey[n] = f1->GetParError(2);
  
     fit_mean[n]=f1->GetParameter(1);  // get the mean from fit function
    
     //fit_means[n]=f1->GetParError(1);
     TLatex *tex = new TLatex(0.6,0.7,Form("%4.0f", fit_mean[n]));
     tex->SetNDC(1); tex->SetTextFont(42); tex->SetTextColor(2); tex->SetTextSize(0.2);
     tex->Draw();
     ofstream myfile;
     myfile.open (Form("Run%d_cell_numbers_means_parrots.txt", nrun));
     if( !myfile ) 
       {
	 cerr << "Can't open output file!" << endl;
	}
     else
       {
	 myfile << Form("# Cell amplitude means for Run %d", nrun) << endl;
	 myfile << "# Cell_idx   mean   HV  v_new   gain_change   " << endl;  // just for testing (true for runs up to 55)
	 for (int idx=0; idx <nblocks;idx++) 
	   {
	     Int_t HV = rand()%99+700;
	     double mean = fit_mean[idx];
	     double v_old = HV;
	     double power = 0.167; //1./(a*n);
	     double g_new = 100;
	     double g_old = mean;
	     double v_new;
	     double gain_change;
	     if (g_old <= 0 ) { // fixme for negative too|| g_old < 0
	       gain_change = 0;
	       v_new = 0;
	     } else {
	       v_new = v_old *pow(gain_change,power);
	       gain_change = g_new/g_old;	}
	     
	     //cout <<" fit means = "<< b  << "    " <<" cell number = " <<n<< endl; //print to the screen
	     myfile << Form("%6d  %8.1f   %8d   %8.1f %8.1f" , idx,  mean , HV, v_new,  gain_change ) << endl;
		
	   }
	 myfile.close();
       }
    }
  
  TCanvas *C100= new TCanvas("c100","Time vs count",3200,1500);
  C100->Divide(ncols,nrows, 0, 0);
 
  for(Int_t n=0;n<3*NPS_rows;n++)
    {
      Int_t canv_ind = ncols*(ncols-n%nrows) + (n/nrows) + 1;
      C100->cd(canv_ind);//changes focus to a particular sub-canvas
      h_time[n]->Draw();
      cout<< " PMT number ="<< n  <<"number of entries"<<  h_time[n]->GetEntries() <<endl;
      
      TF1  *f5 = new TF1(Form("f5_%d",n),"gaus");
      f5->SetParameters(h_time[n]->GetMaximum(), h_time[n]->GetMean(), h_time[n]->GetRMS());//, -10, -10 );  // exp seeds reflect no background (small); will need to fix for data with an exp background
      
      cout<<"   Maximum = "<<h_time[n]->GetMaximum() <<"  mean = "<<  h_time[n]->GetMean()  <<"   RMS =  "<<h_time[n]->GetRMS() <<endl;
      mean_good_pulse[n] =  h_time[n]->GetMean();
     
      h_time[n]->Fit(Form("f5_%d",n),"Q","",  h_time[n]->GetMean()-  ( h_time[n]->GetRMS()),h_time[n]->GetMean()+ ( h_time[n]->GetRMS()));//FITMIN, FITMAX);//L” Use log likelihood method (default is chi-square method). To be used when the histogram represents counts
      //“R” Use the range specified in the function range
      fit_means[n]=f5->GetParameter(1);  // get the mean from fit function
      rms_time[n]=f5->GetParameter(2);
      ez[n]=f5->GetParError(2);
      cout<< " mean =" << f5->GetParameter(1) <<" Maximum = "<<f5->GetParameter(0)<< " RMS = "<<f5->GetParameter(2)<<endl;
      TLatex *tex = new TLatex(0.6,0.7,Form("%4.0f", fit_means[n]));
      tex->SetNDC(1); tex->SetTextFont(42); tex->SetTextColor(2); tex->SetTextSize(0.2);
      tex->Draw();
      Double_t x[3*NPS_rows];
      for (Int_t n_i=0; n_i<3*NPS_rows; n_i++) {
	
       x[n_i]=n_i;
       h_time_block_number->SetBinContent( n_i, h_time[n]->GetBinContent(n+2),h_SampWaveForm[SampPMT]->GetBinContent(n_i));
				
      }
      TCanvas *C7= new TCanvas("c7"," h_time_vs_block_number",2000,1500);
   
      h_time_block_number->Draw("colz");
     
     {
       Int_t m =3*NPS_rows;
       
       Double_t x[m], y[m];
       for (Int_t i=0;i< m;i++) {
      x[i] = i ;
      if( fit_means[i]<0)  {
	fit_means[i]=0;
      }
      else
	{
	  y[i] =   fit_means[i];// h_time[n]->GetMean();
	}
   }
       
   TGraph *gr  = new TGraph(m,x,y);
   TCanvas *c8 = new TCanvas("c8","Good pulse average time for bloc i (after gaussian fit of the pulse time spectrum)",

			    				      200,10,600,400);
   gr->SetTitle("Good pulse average time for bloc i (after gaussian fit of the pulse time spectrum)");
   
   
   gr->Draw("AP*");
   
     }
     {
       Int_t m =3*NPS_rows;
       
       Double_t x[m], y[m];
       for (Int_t i=0;i< m;i++) {
	 x[i] = i ;
	 if( fit_mean[i] <= 0)  {
	   fit_mean[i]=0;
	 }
	 else
	   {
	    y[i] =   fit_mean[i];
	    relative_resolution[i] = (rms_amplitude[i]/4) /(fit_mean[i]);
	  }
       }
       
   TGraph *gra  = new TGraph(m,x,y);
   TCanvas *c11 = new TCanvas("c11","Good pulse average amplitude for bloc i (after gaussian fit of the pulse amplitude spectrum)",
			      200,10,600,400);
   gra->SetTitle("Good pulse average amplitude for bloc i (after gaussian fit of the pulse amplitude spectrum");
   gra->Draw("AP*");//AC*
   TGraph *graphs  = new TGraph(m,x,rms_time);
   TCanvas *c12 = new TCanvas("c12","RMS pulse time amplitude for bloc i (after gaussian fit of the pulse amplitude spectrum)");
   graphs->SetTitle("rms pulse time for bloc i (after gaussian fit of the pulse amplitude spectrum)");
   
   graphs->Draw("AP*");
   TGraph *graphss  = new TGraph(m,x,rms_amplitude);
   TCanvas *c120 = new TCanvas("c120","rms pulse average amplitude for bloc i (after gaussian fit of the pulse amplitude spectrum)",
			      200,10,600,400);
   graphss->SetTitle("rms pulse average amplitude for bloc i (after gaussian fit of the pulse amplitude spectrum");
  
   graphss->Draw("AP*");
   
   TGraph *graphsss  = new TGraph(m,x, relative_resolution);
   TCanvas *c123 = new TCanvas("c123","Relative energy resolution for bloc i (after gaussian fit of the pulse amplitude spectrum)",
			      200,10,600,400);
   graphsss->SetTitle("Relative energy resolution for bloc i (after gaussian fit of the pulse amplitude spectrum)");
  
   graphsss->Draw("AP*");
   TCanvas*  cp = new TCanvas("cp","A Simple Graph with error bars",200,10,700,500);
   cp->SetGrid();
   
   TGraphErrors*  gr = new TGraphErrors(n,x,rms_amplitude, 0,ey);
   gr->SetTitle("A Graph with  error bars in energy resolution");
   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(21);
   gr->Draw("AP");
   cp->Update();
   TCanvas*  cp1 = new TCanvas("cp1","A Graph with  error bars in time resolution",200,10,700,500);
   cp1->SetGrid();
   
   TGraphErrors*  grs = new TGraphErrors(n,x,rms_time,0,ez);
   grs->SetTitle("A Graph with  error bars in time resolution");
   grs->SetMarkerColor(4);
   grs->SetMarkerStyle(21);
   grs->Draw("AP");
   cp1->Update();
    
     }
    
    }
}   

