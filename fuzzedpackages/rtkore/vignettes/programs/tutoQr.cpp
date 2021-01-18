#include <STKpp.h>
using namespace STK;
/** @ingroup tutorial */
int main(int argc, char** argv)
{
  Array2D<Real> a(5,4);
  a << 0, 1, 2, 3,
       2, 3, 4, 5,
       2, 1, 6, 7,
       0, 3,-1, 2,
       3,-1, 1, 1;
  stk_cout << _T("STK++ QR decomposition:\n");
  stk_cout << _T("-----------------------\n");
  Qr q(a); q.run();
  stk_cout << _T("R matrix:\n");
  stk_cout << q.R();
  ArrayXX QR = q.R();
  applyLeftHouseholderArray(QR, q.Q());
  stk_cout << _T("QR matrix:\n");
  stk_cout << QR;
  stk_cout << _T("\nlapack QR decomposition:\n");
  stk_cout << _T("--------------------------\n");
  lapack::Qr lq(a); lq.run();
  stk_cout << _T("R matrix:\n");
  stk_cout << lq.R();
  QR = lq.R();
  applyLeftHouseholderArray(QR, lq.Q());
  stk_cout << _T("QR matrix:\n");
  stk_cout << QR;
  return 0;
}
