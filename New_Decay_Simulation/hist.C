{
    TH1F *seq = new TH1F("seq", "Sequential", 500, 0.0, .5);
    fstream file;

    double a,b,c,d,e;
    
    file.open("dalitz_seq.txt");

while(1)
    {
        file >> a  >> b ;//>> endl;
        seq->Fill(a);
	//seq->Fill(b);
        //cout << (a);
	if (!file.good()) break;

    }
    file.close();




    
    TCanvas *c1 = new TCanvas();
    //c1->SetLogy();
    seq->Draw();

}
