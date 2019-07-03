// ROOT script for generating dilepton plots
// run via "root -l -q -b analysis_root.cc"

void draw_spectrum(TChain &ch, TString var, TString file, TString label, TString unit = "") {
  // draw a spectrum with several components, given a TChain and a variable to be plotted
  TCanvas *c = new TCanvas();
  c->SetLogy();

  ch.SetLineColor(1); ch.Draw(var, "weight");
  ch.SetLineColor(2); ch.Draw(var, "weight*(pdgcode[0]==111  && nout==3)", "same");  // π⁰ Dalitz
  ch.SetLineColor(3); ch.Draw(var, "weight*(pdgcode[0]==221  && nout==3)", "same");  // η Dalitz
  ch.SetLineColor(4); ch.Draw(var, "weight*(pdgcode[0]==113  && nout==2)", "same");  // ρ⁰
  ch.SetLineColor(5); ch.Draw(var, "weight*(pdgcode[0]==223  && nout==2)", "same");  // ω
  ch.SetLineColor(6); ch.Draw(var, "weight*(pdgcode[0]==223  && nout==3)", "same");  // ω Dalitz
  ch.SetLineColor(7); ch.Draw(var, "weight*(pdgcode[0]==2214 && nout==3)", "same");  // Δ⁺ Dalitz
  ch.SetLineColor(8); ch.Draw(var, "weight*(pdgcode[0]==2114 && nout==3)", "same");  // Δ⁰ Dalitz

  TH1F *h = (TH1F*)gPad->GetPrimitive("htemp");
  if (unit != "")
    h->GetXaxis()->SetTitle(label + " [" + unit + "]");
  else
    h->GetXaxis()->SetTitle(label);
  h->GetYaxis()->SetTitle("dN/d" + label);
  h->SetTitle("");
  c->Update();
  c->SaveAs(file + ".pdf");
}


int analysis() {
  gROOT->Reset();

  // open several files and add all of them to a TChain
  TChain chain("collisions");
  for (int i=1; i<=10; i++) {
    TString nr;
    nr.Form("%i", i);
    chain.Add("data/" + nr + "/DileptonOutput.root");
  }

  // mass spectrum
  TString inv_mass = "sqrt((p0[1]+p0[2])^2-(px[1]+px[2])^2-(py[1]+py[2])^2-(pz[1]+pz[2])^2)";
  draw_spectrum(chain, inv_mass, "mass", "m_{ee}", "GeV");

  // pT spectrum
  TString pT = "sqrt((px[1]+px[2])^2+(py[1]+py[2])^2)";
  draw_spectrum(chain, pT, "pt", "p_{T}", "GeV");

  // rapidity spectrum
  TString rap = "0.5 * log((p0[1]+p0[2] + (pz[1]+pz[2])) / (p0[1]+p0[2] - (pz[1]+pz[2])))";
  draw_spectrum(chain, rap, "y",  "y");

  return 0;
}
