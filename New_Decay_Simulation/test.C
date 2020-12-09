{
    TH1F *seq = new TH1F("seq", "Sequential", 500, 0.0, 0.4);
    TH1F *ddl = new TH1F("ddl", "DDL", 500, 0.0, 0.4);
    TH1F *dde = new TH1F("dde", "DDE", 500, 0.0, 0.4);
    TH1F *ddphi = new TH1F("ddphi", "DDphi", 500, 0.0, 0.4);
    TH1F *data = new TH1F("combined", "Combined", 500, 0.0, 0.4);

    // retrieve histograms

    fstream file;

    double a,b,c,d,e;
    
    file.open("dalitz_seq.txt");

    while(1)
    {
        file >> a  >> b ;//>> endl;
        //seq->Fill(a);
	seq->Fill(b);
        //cout << (a);
	if (!file.good()) break;

    }
    file.close();

    file.open("dalitz_ddl.txt");

    while(1)
    {
        file >> a  >> b ;//>> endl;
        //ddl->Fill(a);
	ddl->Fill(b);
        //cout << (a);
	if (!file.good()) break;

    }
    file.close();

    file.open("dalitz_dde.txt");

    while(1)
    {
        file >> a  >> b ;//>> endl;
        //dde->Fill(a);
	dde->Fill(b);
        //cout << (a);
	if (!file.good()) break;

    }
    file.close();


    file.open("dalitz_ddphi.txt");

    while(1)
    {
        file >> a  >> b ;//>> endl;
        //ddphi->Fill(a);
	ddphi->Fill(b);
        //cout << (a);
	if (!file.good()) break;

    }
    file.close();


    file.open("dalitz_all.txt");

    while(1)
    {
        file >> a  >> b ;//>> endl;
        //data->Fill(a);
	data->Fill(b);
        //cout << (a);
	if (!file.good()) break;

    }
    file.close();


   TObjArray *mc = new TObjArray(4);        // MC histograms are put in this array
   mc->Add(seq);
   mc->Add(ddl);
   mc->Add(dde);
   mc->Add(ddphi);


   TFractionFitter* fit = new TFractionFitter(data, mc); // initialise
   fit->Constrain(0,0.99,1.0);               // constrain fraction 1 to be between 0 and 1
   fit->Constrain(1,0.00001,0.001);
   fit->Constrain(2,0.0,0.0001);
   fit->Constrain(3,0.0,0.0001);                // constrain fraction 2 to be between 0 and 1
   //fit->SetRangeX(0,400);                    // use only the first 15 bins in the fit
   Int_t status = fit->Fit();               // perform the fit
   std::cout << "fit status: " << status << std::endl;
   if (status == 0) {                       // check on fit status
     TH1F* result = (TH1F*) fit->GetPlot();
     data->Draw("Ep");
     result->Draw("same");
   }
}
