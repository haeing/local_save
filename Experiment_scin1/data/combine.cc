void combine(){
  TChain ch("data");
  ch.Add("210930.root");
  ch.Add("211001.root");


  ch.Merge("211005.root");
}
