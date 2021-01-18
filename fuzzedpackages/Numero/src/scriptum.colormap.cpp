/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "scriptum.local.h"

/*
 *
 */
class RGB {
public:
  mdreal red;
  mdreal green;
  mdreal blue;
public:
  RGB() {};
  RGB(mdsize r, mdsize g, mdsize b) {
    this->red = r/255.0;
    this->green = g/255.0;
    this->blue = b/255.0;
  };
  ~RGB() {};
};

/*
 *
 */
class Palette {
private:
  vector<RGB> hues;
  vector<mdreal> qpoints;
public:
  Palette() {};
  Palette(const string& name);
  ~Palette() {};
  Color color(const mdreal) const;
};

static unordered_map<string, Palette> PaletteCache;

/*
 * Return red, green and blue components.
 */
Color
scriptum::colormap(const mdreal q, const string& name) {
  if(PaletteCache.size() < 1) {
    PaletteCache["grey"] = Palette("grey");
    PaletteCache["fire"] = Palette("fire");
    PaletteCache["jungle"] = Palette("jungle");
    PaletteCache["miami"] = Palette("miami");
    PaletteCache["rhodo"] = Palette("rhodo");
    PaletteCache["tan"] = Palette("tan");
    PaletteCache["default"] = Palette("");
  }
  if(PaletteCache.count(name) < 1)
    return PaletteCache["default"].color(q);
  return PaletteCache[name].color(q);
}

/*
 *
 */
Palette::Palette(const string& name) {

  /* Grey-scale. */
  if(name == "grey") {
    hues.push_back(RGB(  0,   0,   0));
    hues.push_back(RGB(127, 127, 127)); 
    hues.push_back(RGB(255, 255, 255));
  }

  /* Red palette. */
  if(name == "fire") {
    hues.push_back(RGB(  0,   0,   0));
    hues.push_back(RGB(150,   5,   0));
    hues.push_back(RGB(245,  90,   0)); 
    hues.push_back(RGB(255, 170,  20)); 
    hues.push_back(RGB(255, 220, 100)); 
    hues.push_back(RGB(255, 240, 150));
    hues.push_back(RGB(255, 255, 255));
  }

  /* Green palette. */
  if(name == "jungle") {
    hues.push_back(RGB(240, 255, 150));
    hues.push_back(RGB(220, 245, 100)); 
    hues.push_back(RGB(170, 225,  20)); 
    hues.push_back(RGB( 90, 205,   0)); 
    hues.push_back(RGB(  5, 120,   0)); 
  }

  /* Saturated but balanced rainbow hues. */
  if(name == "miami") {
    hues.push_back(RGB(255,  71, 189));
    hues.push_back(RGB(255,  65,  50));
    hues.push_back(RGB(255, 126,  25));
    hues.push_back(RGB(220, 213,   0));
    hues.push_back(RGB( 76, 245,  50));
    hues.push_back(RGB( 10, 190, 213));
    hues.push_back(RGB( 35, 130, 255));
  }

  /* Red-blue hues. */
  if(name == "rhodo") {
    hues.push_back(RGB( 40,  10, 220));
    hues.push_back(RGB(  0, 100, 255));
    hues.push_back(RGB( 80, 150, 255));
    hues.push_back(RGB(120, 210, 255));
    hues.push_back(RGB(255, 255, 255));
    hues.push_back(RGB(255, 200, 120));
    hues.push_back(RGB(255, 140,  80));
    hues.push_back(RGB(255,  90,   0));
    hues.push_back(RGB(220,  15,  40));
  }

  /* Dark-bright hues. */
  if(name == "tan") {
    hues.push_back(RGB(  0,   0,   0));
    hues.push_back(RGB(125,  85,  65)); 
    hues.push_back(RGB(205, 140, 100)); 
    hues.push_back(RGB(255, 210, 150)); 
    hues.push_back(RGB(255, 255, 255));
  }

  /* Default  colors. */
  if(hues.size() < 1) {
    hues.push_back(RGB(255,   0,   0));
    hues.push_back(RGB(  0, 255,   0)); 
    hues.push_back(RGB(  0,   0, 255));
  }
  
  /* Determine quantile positions. */
  mdsize n = hues.size();
  (this->qpoints).resize(n);
  for(mdsize i = 0; i < n; i++)
    qpoints[i] = 1.0*i/(n - 1.0);
}

/*
 *
 */
Color
Palette::color(const mdreal q0) const {
  mdsize sznan = medusa::snan();
  mdreal rlnan = medusa::rnan();
  if(q0 == rlnan) return Color();

  /* Clip extreme values. */
  mdreal q = q0;
  if(q < 0.0) q = 0.0;
  if(q > 1.0) q = 1.0;

  /* Find quantile segment. */
  Site slot = binsearch(qpoints, q);
  mdsize a = slot.bounds.first;
  mdsize b = slot.bounds.second;
  if(a == sznan) panic("Inconsistent state.", __FILE__, __LINE__);
  if(b == sznan) panic("Inconsistent state.", __FILE__, __LINE__);

  /* Interpolate color. */
  Color c;
  mdreal wA = slot.weights.first;
  mdreal wB = slot.weights.second;
  c.red = (wA*(hues[a].red) + wB*(hues[b].red));
  c.green = (wA*(hues[a].green) + wB*(hues[b].green));
  c.blue = (wA*(hues[a].blue) + wB*(hues[b].blue));
  c.opacity = 1.0;
  return c;
}
