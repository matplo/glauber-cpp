#include "tglaubermc.cxx"
#include <fstream>
#include <iomanip>
void reference_fast(unsigned trials=15000, unsigned seed=72911, const char* path="fast.csv") {
  gRandom->SetSeed(seed);
  TGlauberMC g("Pbpnrw","Pbpnrw",67.6);
  g.SetBmin(0);
  g.SetBmax(20);
  g.SetMinDistance(0.4);
  g.SetNodeDistance(0.);
  std::ofstream out(path);
  out << "trial,b_fm,npart,ncoll\n" << std::setprecision(12);
  for(unsigned i=0;i<trials;++i)
    if(g.NextEvent()) out << i+1 << ',' << g.GetB() << ',' << g.GetNpart() << ',' << g.GetNcoll() << '\n';
  out.flush();
}
