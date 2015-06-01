void rootlogon()
{
  printf("Start Loading RooUnfold\n");
  gSystem->Load("/uscms/home/zlesko/ZFinder/Dev/ZScripts/ZFinder-Analysis-Scripts/ZFinder/RooUnfold-1.1.1/libRooUnfold");
  gSystem->AddIncludePath("-I/uscms/home/zlesko/ZFinder/Dev/ZScripts/ZFinder-Analysis-Scripts/ZFinder/RooUnfold-1.1.1/src");
  gStyle->SetPalette(1);
  printf("Finished Loading RooUnfold\n");
}
