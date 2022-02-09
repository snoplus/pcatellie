void fit(){
  int ipw [40] = { 0,250,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000,3250,3500,3750,4000,4250,4500,4750,5000,5250,5500,5750,6000,6250,6500,6750,7000,7250,7500,7750,8000,8250,8500,8750,9000,9250,9500,9750 };
  int pin [40] = { 65535,65535,65535,65535,65535,65535,65535,65535,65535,65535,65535,65535,65535,65535,64556,60096,55728,51299,46795,42133,37290,32337,27362,22479,17830,13475,9465,6025,3248,1379,747,720,721,721,715,726,716,723,726,726};

  TCanvas* c1 = new TCanvas("c1","c1",1024, 768);
  TGraph *plot  = new TGraph();

  for (int i = 0; i < 41; i++){
    plot->SetPoint(i,ipw[i],pin[i]);
  }

  plot->SetName("IPW vs PIN");
  plot->SetTitle("IPW vs PIN");
  plot->GetXaxis()->SetTitle("ipw");
  plot->GetYaxis()->SetTitle("PIN");
  plot->SetMarkerStyle(5);
  plot->SetMarkerSize(2);
  //plot->GetXaxis()->SetRangeUser(0.0,fSelectAngle+0.1);
  //plot->GetYaxis()->SetRangeUser(bmin - 0.01,bmax + 0.01);
  plot->Draw("AP");

}
