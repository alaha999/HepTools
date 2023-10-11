#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>


void anaSkim(int sample=1)
{
  const char *hstfilename, *sumfilename, *skimfilename;

  TChain *chain = new TChain("Events");
  NanoSkim m_selec;

  if(sample==1){
    chain->Add("/home/sdube/Sourabh/Work/FullRun2/testSample/tree_50.root");
    hstfilename = "hst_DYJetsToLLM50_RunIISummer16.root";
    sumfilename = "sum_DYJetsToLLM50_RunIISummer16.txt";
    skimfilename = "skimmed_tree_50.root";
    m_selec.SetData(0); //0 - running over MC, 1 - running over Data
    m_selec.SetYear(2018);
    //m_selec.SetMCwt(1);
  }
  
  std::cout<<"Output files are "<<hstfilename<<" and "<<sumfilename<<std::endl;
  m_selec.SetHstFileName(hstfilename);
  m_selec.SetSumFileName(sumfilename);
  m_selec.SetSkimFileName(skimfilename);
  m_selec.SetVerbose(10);//set verbosity level for output.

  chain->Process(&m_selec);


}

int main(int argc, char *argv[])
{

  if(argc<2){
    std::cout<<" Please give one integer argument "<<std::endl;
    return 0;
  }
  int sample_id = atoi(argv[1]);

  anaSkim(sample_id);
  return 0;
}
