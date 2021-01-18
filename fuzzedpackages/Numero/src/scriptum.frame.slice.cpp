/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "scriptum.local.h"

/*
 *
 */
bool
Frame::slice(const mdreal x, const mdreal y, const mdreal rA,
	     const mdreal rB, const mdreal phi0, const mdreal phiB) {
  mdreal rlnan = medusa::rnan();
  FrameBuffer* p = (FrameBuffer*)buffer;

  /* Check if unusable coordinates. */
  if(x == rlnan) return false;
  if(y == rlnan) return false;
  if(rA == rlnan) return false;
  if(rB == rlnan) return false;
  if(phi0 == rlnan) return false;
  if(phiB == rlnan) return false;
  if(rA < 0.0) return false; /* zero allowed */
  if(rB <= rA) return false;

  /* Check if circle. */
  double rotat = 0.0;
  double delta = fabs(phiB - phi0);
  delta = modf(delta/360.0, &rotat);
  if(delta < 1e-6) {
    sprintf(p->f(), "\n<circle ");
    sprintf(p->f(), "cx=\"%.2f\" cy=\"%.2f\" ", x, y);
    sprintf(p->f(), "r=\"%.3f\"\n", rB);
    sprintf(p->f(), "%s", (p->linestycode).c_str());
    sprintf(p->f(), "/>\n");
    return true;
  }

  /* Skip full rotations. */
  double phiA = (phi0 + rotat*360.0);
  delta -= rotat*360.0;

  /* Determine end—points. */
  vector<mdreal> xp(6, rlnan);
  vector<mdreal> yp(6, rlnan);
  xp[0] = rA*cos(3.141593*phiB/180.0);
  yp[0] = rA*sin(3.141593*phiB/180.0);
  xp[1] = rA*cos(3.141593*phiA/180.0);
  yp[1] = rA*sin(3.141593*phiA/180.0);
  xp[2] = rB*cos(3.141593*phiA/180.0);
  yp[2] = rB*sin(3.141593*phiA/180.0);
  xp[3] = rB*cos(3.141593*phiB/180.0);
  yp[3] = rB*sin(3.141593*phiB/180.0);

  /* Set position. */
  if((x != 0.0) || (y != 0.0))
    sprintf(p->f(), "\n<path transform=\"translate(%.2f,%.2f)\"\n", x, y);
  else
    sprintf(p->f(), "\n<path\n");

  /* Draw slice. */
  sprintf(p->f(), "d=\"M %.4f %.4f\n", xp[0], yp[0]);
  sprintf(p->f(), "A %.4f %.4f 0 0 0\n", rA, rA);
  sprintf(p->f(), "%.4f %.4f\n", xp[1], yp[1]);
  sprintf(p->f(), "L %.4f %.4f\n", xp[2], yp[2]);
  sprintf(p->f(), "A %.4f %.4f 0 0 1\n", rB, rB);
  sprintf(p->f(), "%.4f %.4f Z\"\n", xp[3], yp[3]);
  sprintf(p->f(), "%s", (p->linestycode).c_str());
  sprintf(p->f(), "/>\n");

  /* Update limits. */
  for(mdsize i = 0; i < xp.size(); i++) {
    if(xp[i] != rlnan) xp[i] += x;
    if(yp[i] != rlnan) yp[i] += y;
  }
  (p->limits).first.update(xp, p->style);
  (p->limits).second.update(yp, p->style);
  return true;
}
