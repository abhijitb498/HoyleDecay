{
    TH1F *seq = new TH1F("seq", "Sequential", 500, 0.0, 0.25);
    TH1F *ddl = new TH1F("ddl", "DDL", 100, 0.0, 0.25);
    TH1F *dde = new TH1F("dde", "DDE", 100, 0.0, 0.25);
    TH1F *ddphi = new TH1F("ddphi", "DDphi", 100, 0.0, 0.25);
    TH1F *combined = new TH1F("combined", "Combined", 500, 0.0, 0.25);
    fstream file;

    double a,b,c,d,e;
    
    file.open("dalitz_seq.txt");

    while(1)
    {
        file >> a ;// >> b >> c;
        seq->Fill(a);
        //cout << (a);
	if (!file.good()) break;

    }
    file.close();

    file.open("dalitz_ddl.txt");

    while(1)
    {
        file >> b ;// >> b >> c;
        ddl->Fill(b);
        //cout << (a);
	if (!file.good()) break;

    }
    file.close();

    file.open("dalitz_dde.txt");

    while(1)
    {
        file >> c ;// >> b >> c;
        dde->Fill(c);
        //cout << (a);
	if (!file.good()) break;

    }
    file.close();


    file.open("dalitz_ddphi.txt");

    while(1)
    {
        file >> d ;// >> b >> c;
        ddphi->Fill(d);
        //cout << (a);
	if (!file.good()) break;

    }
    file.close();


    file.open("dalitz_all.txt");

    while(1)
    {
        file >> e ;// >> b >> c;
        combined->Fill(e);
        //cout << (a);
	if (!file.good()) break;

    }
    file.close();


    
    TCanvas *c1 = new TCanvas();
    c1->SetLogy();
    seq->GetXaxis()->SetTitle("E_rms MeV");
    seq->GetYaxis()->SetTitle("Entries");
    //seq->GetYaxis()->SetRangeUser(0,9000);
    seq->SetTitle("rms Energy Plot");
    seq->SetStats(0);
    seq->SetLineColor(2);
    seq->Draw();

    ddl->SetLineColor(3);
    ddl->Draw("SAME");

    dde->SetLineColor(6);
    dde->Draw("SAME");

    ddphi->SetLineColor(7);
    ddphi->Draw("SAME");
    //hist.Fit("gaus");

    combined->SetLineColor(1);
    combined->Draw("SAME");

    TLegend *leg = new TLegend(0.6,0.7,0.9,0.9);
    leg->AddEntry(seq,"Sequential","l");
    leg->AddEntry(ddl,"DDL","l");
    leg->AddEntry(dde,"DDE","l");
    leg->AddEntry(ddphi,"DDPhi","l");
    leg->AddEntry(combined,"Combined","l");
    leg->Draw();
}
