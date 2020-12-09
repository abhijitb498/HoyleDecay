#include "Riostream.h"
{
    TH1F *seq = new TH1F("seq", "Sequential", 500, 0.0, 0.25);
    //TH1F *ddl = new TH1F("ddl", "DDL", 100, 0.05, 0.25);
    //TH1F *dde = new TH1F("dde", "DDE", 100, 0.05, 0.25);
    //TH1F *ddphi = new TH1F("ddphi", "DDphi", 100, 0.05, 0.25);
    //TH1F *combined = new TH1F("combined", "Combined", 500, 0.05, 0.25);
    ifstream file;

    Float_t a,b,c,d,e;
    int i;
    i = 0;
    //b = 0.;
    
    file.open("dalitz_seq.txt");

    while(i<4000000)
    {
        file >> a;// >> b >> c;
        seq->Fill(a);
        cout << a << endl;
	i = i+1;
 	if (!file.good()) break;
    }
    file.close();

    TCanvas *c1 = new TCanvas();
    //c1->SetLogy();
    seq->GetXaxis()->SetTitle("Relative energy in MeV");
    seq->GetYaxis()->SetTitle("Entries");
    //seq->GetYaxis()->SetRangeUser(0,9000);
    seq->SetTitle("Relative Energy of 8Be like pair");
    seq->Draw();
}
