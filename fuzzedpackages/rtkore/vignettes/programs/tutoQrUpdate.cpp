#include <STKpp.h>
using namespace STK;
/** @ingroup tutorial */
int main(int argc, char** argv)
{
  ArrayXX a(5,4), QR;
  a << 0, 1, 2, 3,
       2, 3, 4, 5,
       2, 1, 6, 7,
       0, 3,-1, 2,
       3,-1, 1, 1;
  stk_cout << "lapack QR decomposition:\n";
  lapack::Qr lq(a); lq.run();
  // remove column
  lq.eraseCol(2);
  stk_cout << "R matrix:\n";
  stk_cout << lq.R();
  QR = lq.Q() * lq.R();
  stk_cout << "QR matrix:\n";
  stk_cout << QR;
  // insert constant column with value 1 in column 2
  lq.insertCol(Const::Vector<Real>(5), 2);
  stk_cout << "\nR matrix:\n";
  stk_cout << lq.R();
  QR = lq.Q() * lq.R();
  stk_cout << "QR matrix:\n";
  stk_cout << QR;
  return 0;
}
