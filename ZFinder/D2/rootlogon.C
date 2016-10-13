void rootlogon()
{
  gSystem->Load("/home/user1/lesko/work/NicoleStuff/ZFinder/RooUnfold-1.1.1/libRooUnfold");
  gSystem->AddIncludePath("-I/home/user1/lesko/work/NicoleStuff/ZFinder/RooUnfold-1.1.1/src");
  gStyle->SetPalette(1);
  printf("Finished Loading RooUnfold\n");
}
