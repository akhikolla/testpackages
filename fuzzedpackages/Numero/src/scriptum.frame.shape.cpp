/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "scriptum.local.h"

static string make_path(const mdreal, const mdreal, const mdreal,
			const string&);


/*
 *
 */
bool
Frame::shape(const mdreal x, const mdreal y, const mdreal r,
	     const string& name) {
  mdreal rlnan = medusa::rnan();
  FrameBuffer* p = (FrameBuffer*)buffer;

  /* Check if unusable coordinates. */
  if(x == rlnan) return false;
  if(y == rlnan) return false;
  if(r == rlnan) return false;
  if(r <= 0.0) return false;

  /* Outline. */
  string path = make_path(x, y, r, name);
  if(path.size() < 1) return false;
  p->append(path);
  
  /* Set style. */
  p->append(p->linestycode);
  p->append("/>\n");

  /* Update limits. */
  vector<mdreal> xp(2, x);
  vector<mdreal> yp(2, y);
  xp[0] -= r; xp[1] += r;
  yp[0] -= r; yp[1] += r;
  (p->limits).first.update(xp, p->style);
  (p->limits).second.update(yp, p->style);
  return true;
}

/*
 *
 */
string
make_path(const mdreal x, const mdreal y, const mdreal r,
	  const string& name) {
  char buf[512];
  char* ptr = buf;
  vector<mdreal> px, py;
  if(name.size() < 1) return "";

  /* Create a circle. */
  if(name.substr(0, 6) == "circle") {
    ptr += sprintf(ptr, "\n<circle ");
    ptr += sprintf(ptr, "cx=\"%.2f\" cy=\"%.2f\" ", x, y);
    ptr += sprintf(ptr, "r=\"%.3f\"\n", r);
    return string(buf);
  }

  /* Determine symbol rotation. */
  string s = (name + "|0.0 ");
  const char* array = s.c_str();
  char* pos = (char*)strchr(array, '|'); pos++;
  double theta = atof(pos)*3.1416/180;

  /* Create polygonal shapes. */
  switch(s[0]) {
  case 'c':
    if(s.substr(0, 6) == "clover") {
      px.resize(12); py.resize(12);
      px[0]  = -0.5664; py[0]  = -0.3270;
      px[1]  = -0.7616; py[1]  = -0.7798;
      px[2]  =  0.0000; py[2]  = -1.1445;
      px[3]  =  0.7616; py[3]  = -0.7798;
      px[4]  =  0.5664; py[4]  = -0.3270;
      px[5]  =  1.0561; py[5]  = -0.2697;
      px[6]  =  0.9912; py[6]  =  0.5722;
      px[7]  =  0.2945; py[7]  =  1.0494;
      px[8]  =  0.0000; py[8]  =  0.6540;
      px[9]  = -0.2945; py[9]  =  1.0495;
      px[10] = -0.9911; py[10] =  0.5723;
      px[11] = -1.0561; py[11] = -0.2696;
    }
    if(s.substr(0, 5) == "cross") {
      px.resize(12); py.resize(12);
      px[0]  = -0.5810; py[0]  = -0.5810;
      px[1]  = -1.0294; py[1]  = -0.4289;
      px[2]  = -1.0294; py[2]  =  0.4289;
      px[3]  = -0.5810; py[3]  =  0.5810;
      px[4]  = -0.4289; py[4]  =  1.0294;
      px[5]  =  0.4289; py[5]  =  1.0294;
      px[6]  =  0.5810; py[6]  =  0.5810;
      px[7]  =  1.0294; py[7]  =  0.4289;
      px[8]  =  1.0294; py[8]  = -0.4289;
      px[9]  =  0.5810; py[9]  = -0.5810;
      px[10] =  0.4289; py[10] = -1.0294;
      px[11] = -0.4289; py[11] = -1.0294;
    }
    break;
  case 's':
    if(s.substr(0, 6) == "square") {
      px.resize(4); py.resize(4);
      px[0]  = -0.9400; py[0]  = -0.9400;
      px[1]  =  0.9400; py[1]  = -0.9400;
      px[2]  =  0.9400; py[2]  =  0.9400;
      px[3]  = -0.9400; py[3]  =  0.9400;
     }
    if(s.substr(0, 4) == "star") {
      px.resize(12); py.resize(12);
      px[0]  =  0.8125; py[0]  =  0.0000;
      px[1]  =  1.0825; py[1]  =  0.6250;
      px[2]  =  0.4063; py[2]  =  0.7036;
      px[3]  =  0.0000; py[3]  =  1.2500;
      px[4]  = -0.4062; py[4]  =  0.7036;
      px[5]  = -1.0825; py[5]  =  0.6250;
      px[6]  = -0.8125; py[6]  =  0.0000;
      px[7]  = -1.0825; py[7]  = -0.6250;
      px[8]  = -0.4063; py[8]  = -0.7036;
      px[9]  = -0.0000; py[9]  = -1.2500;
      px[10] =  0.4063; py[10] = -0.7036;
      px[11] =  1.0825; py[11] = -0.6250;
    }
    break;
  case 't':
    if(s.substr(0, 8) == "triangle") {
      px.resize(6); py.resize(6);
      px[0]  =  0.1991; py[0]  =  1.1715;
      px[1]  = -0.1991; py[1]  =  1.1715;
      px[2]  = -1.1683; py[2]  = -0.5071;
      px[3]  = -0.9691; py[3]  = -0.8520;
      px[4]  =  0.9691; py[4]  = -0.8520;
      px[5]  =  1.1683; py[5]  = -0.5071;
    }
    break;
  case 'p':
    if(s.substr(0, 8) == "pentagon") {
      px.resize(5); py.resize(5);
      px[0]  =  1.0937; py[0]  =  0.3539;
      px[1]  =  0.0000; py[1]  =  1.1485;
      px[2]  = -1.0937; py[2]  =  0.3539;
      px[3]  = -0.6760; py[3]  = -0.9319;
      px[4]  =  0.6760; py[4]  = -0.9319;
    }
    break;
  default:
    return "";
  }
  
  /* Rotate. */
  for(mdsize i = 0; i < px.size(); i++) {
    double phi = 0.0;
    double rad = sqrt(px[i]*px[i] + py[i]*py[i]);
    if(px[i] > 0) phi = atan(py[i]/(px[i] + FLT_EPSILON));
    else phi = (3.1416 + atan(py[i]/(px[i] - FLT_EPSILON)));
    px[i] = rad*cos(phi + theta);
    py[i] = rad*sin(phi + theta);
  }

  /* Draw polygon. */
  ptr += sprintf(ptr, "\n<polygon points=\"");
  for(mdsize i = 0; i < px.size(); i++) {
    ptr += sprintf(ptr, "\n\t%.2f,", (x + r*px[i]));
    ptr += sprintf(ptr, "%.2f", (y + r*py[i]));
  }
  ptr += sprintf(ptr, "\"\n");
  return string(buf);
}
