Double_t FittingF_GE (Double_t *x, Double_t *par) {
  //par[0] - gauss ampl
  //par[1] - gauss mean
  //par[2] - gause sigmna
  //par[3] - exp start t
  //par[4] - exp ampl
  //par[5] - exp decay time
  return par[0]*std::exp(-0.5*std::pow((x[0]-par[1])/par[2], 2))
      + (x[0] > par[3] ? par[4]*std::exp(-1.0*x[0]/par[5]) : 0);
}

Double_t AnodeCurrentF (Double_t *x, Double_t *par) {
  //par[0] - gauss ampl
  //par[1] - exp start t
  //par[2] - exp ampl
  //par[3] - exp decay time
  return (x[0] > par[1] ? par[2]*std::exp(-1.0*(x[0]-par[1])/par[3]) : par[0]);
}

Double_t ElectronF (Double_t *x, Double_t *par) {
  //par[0] - exp finish t
  //par[1] - exp ampl
  //par[2] - exp decay time
  return (x[0] < par[0] ? par[1]*std::exp((x[0]-par[0])/par[2]) : 0.0);
}

void convolution(void) {
  gStyle->SetStatY(0.9);
  gStyle->SetStatX(0.9);

  double Dtime_left = 0, Dtime_right = 60;
  double dt = 2e-2; //Integration time in us
  unsigned int Np = (Dtime_right - Dtime_left)/dt;
  Double_t *Is, *Qs, *Sig, *ts;
  ts = new Double_t[Np];
  Is = new Double_t[Np];
  Qs = new Double_t[Np];
  Sig = new Double_t[Np];

  TF1 *fc = new TF1("AnodeCurrentF", AnodeCurrentF, Dtime_left, Dtime_right, 4);
  fc->SetParNames("Amplitude", "Slow t0[s]", "Slow amplitude", "Slow tau[s]");
  fc->SetParameter(0, 5.2);
  fc->SetParameter(1, 3.05);
  fc->SetParameter(2, 0.0);
  fc->SetParameter(3, 6.841e-1);

  double ComptonL = 1.8; //milimeters
  double ComptonTau = ComptonL/1.85e6;//in s

  TF1 *fe = new TF1("ElectronF", ElectronF, Dtime_left, Dtime_right, 3);
  fe->SetParNames("t0[s]", "Amplitude", "tau[s]");
  fe->SetParameter(0, 26);//in us, 48mm LAr, 7.0Td id gas in our detector
  fe->SetParameter(1, 10.0);
  fe->SetParameter(2, ComptonTau*1e6);

  double Int = 0;
  for (unsigned int i = 0; i<Np; ++i) {
	  ts[i] = Dtime_left + i*(Dtime_right - Dtime_left)/(Np-1);
    Is[i] = (*fc)(ts[i]);
    Qs[i] = (*fe)(ts[i]);
    Sig[i] = 0;
    for (unsigned int j = 0; j<Np; ++j) {
	    double t2 = Dtime_left + j*(Dtime_right - Dtime_left)/(Np-1);
      Sig[i] += (*fc)(t2)*(*fe)(ts[i]-t2)*dt;
    }
  }
  TCanvas *c1 = new TCanvas("Charge Signal","Charge Signal", 900, 700);
	TGraph *grSig = new TGraph(Np, ts, Sig);
  TGraph *grIA = new TGraph(Np, ts, Is);
  TGraph *grQe = new TGraph(Np, ts, Qs);
	grSig->Draw("AL");
  TCanvas *c2 = new TCanvas("Anode Current","Anode Current", 900, 700);
  grIA->Draw("AL");  
  TCanvas *c3 = new TCanvas("Electrons","Electrons", 900, 700);
  grQe->Draw("AL");
  delete [] ts;
  delete [] Sig;
  delete [] Is;
  delete [] Qs;
}

 


