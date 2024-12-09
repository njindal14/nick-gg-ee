// this macro runs analyze.cxx, which takes as input an Ascii starlight output
// file, slight.out, and creates a standard set of histograms, which are stored
// in histograms.root  
// just modifying a file for svn -jb 7.11.17
#include "SLAnalyze.cxx"


void SLAna()
{
  //uranium
  Analyze a("/Users/Nick/STAR/starlight/utils/slight.out", 10000000);

  //gold 
  //Analyze a("/Users/Nick/STAR/docker_mount/JobFilesAu1M/slight.out", 1000000);
  a.doAnalysis();
}
