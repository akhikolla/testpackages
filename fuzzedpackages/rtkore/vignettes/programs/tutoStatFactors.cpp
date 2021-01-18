#include "STKpp.h"
using namespace STK;
int main(int argc, char** argv)
{
  CArray<char, 13, 3> fac;
  fac <<'b','a','a', 'c','b','a', 'b','a','a', 'c','a','a', 'd','a','a'
       ,'b','b','a', 'd','c','b', 'c','c','b', 'b','c','b', 'd','b','b'
       ,'b','a','b', 'd','a','b', 'c','c','b';

  Stat::Factor<CArrayVector<char, 13> > f1d(fac.col(1));
  f1d.run();
  stk_cout << "nbLevels= "  << f1d.nbLevels() << "\n";
  stk_cout << "asInteger=\n"<< f1d.asInteger().transpose() << "\n";
  stk_cout << "Levels: "     << f1d.levels().transpose();
  stk_cout << "Levels counts: " << f1d.counts().transpose();

  Stat::MultiFactor<CArray<char, 13, 3> > f2d(fac);
  f2d.run();
  stk_cout << "nbLevels= "  << f2d.nbLevels() << "\n";
  stk_cout << "asInteger=\n"<< f2d.asInteger().transpose() << "\n";
  for (int i=f2d.levels().begin(); i<f2d.levels().end(); ++i)
    stk_cout <<"Levels "<< i <<": "<< f2d.levels()[i].transpose();
  stk_cout << "\n";
  for (int i=f2d.levels().begin(); i<f2d.levels().end(); ++i)
    stk_cout <<"Counts "<< i << ": "<< f2d.counts()[i].transpose();
  return 0;
}
