/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "scriptum.local.h"

static void update_color(Color&, const Color&);

/*
 *
 */
void
Frame::stylize(const Style& st) {
  FrameBuffer* p = (FrameBuffer*)buffer;
  mdreal rlnan = medusa::rnan();

  /* Pointability. */
  Style& base = p->style;
  base.pointable = st.pointable;
  
  /* Text-anchor. */
  if(st.anchor == "start") base.anchor = st.anchor;
  if(st.anchor == "middle") base.anchor = st.anchor;
  if(st.anchor == "end") base.anchor = st.anchor;

  /* Rotation angle in degrees [-180, 180]. */
  if(st.angle == 0.0) base.angle = 0.0;
  if(st.angle != rlnan) {
    long nrot = (long)(fabs(st.angle)/360.0);
    if(st.angle < 0.0) base.angle = (st.angle + nrot*360.0);
    if(st.angle > 0.0) base.angle = (st.angle - nrot*360.0);
    if(base.angle < -180.0) base.angle += 360.0;
    if(base.angle > 180.0) base.angle -= 360.0;
  }

  /* Fill color. */
  update_color(base.fillcolor, st.fillcolor);

  /* Font family. */
  if(st.fontfamily.size() > 0) base.fontfamily = st.fontfamily;

  /* Font size. */
  if(st.fontsize != rlnan) base.fontsize = st.fontsize;
  if(base.fontsize < 0.0) base.fontsize = 0.0;
  if(base.fontsize > 1e3) base.fontsize = 1e3;
    
  /* Font weight. */
  if((st.fontweight >= 100) && (st.fontweight <= 900))
    base.fontweight = st.fontweight;

  /* Origin of the element coordinate system. */
  if((st.origin).size() >= 2) {
    if(st.origin[0] != rlnan) base.origin[0] = st.origin[0];
    if(st.origin[1] != rlnan) base.origin[1] = st.origin[1];
    if(base.origin[0] < -1e4) base.origin[0] = -1e4;
    if(base.origin[0] > 1e4) base.origin[0] = 1e4;
    if(base.origin[1] < -1e4) base.origin[1] = -1e4;
    if(base.origin[1] > 1e4) base.origin[1] = 1e4;
  }

  /* Minimum distance to canvas margin. */
  if(st.padding != rlnan) base.padding = st.padding;
  if(base.padding < 0.0) base.padding = 0.0;
  if(base.padding > 1e3) base.padding = 1e3;

  /* Stroke color. */
  update_color(base.strokecolor, st.strokecolor);

  /* Stroke width. */
  if(st.strokewidth != rlnan) base.strokewidth = st.strokewidth;
  if(base.strokewidth < 0.0) base.strokewidth = 0.0;
  if(base.strokewidth > 1e3) base.strokewidth = 1e3;

  /* Identity and value. */
  base.identity = medusa::string2safe(st.identity, 255);
  base.values.resize(st.values.size());
  for(mdsize i = 0; i < st.values.size(); i++)
    base.values[i] = medusa::string2safe(st.values[i], 255);
  
  /* Apply changes. */
  scriptum_local::style2code(p->linestycode, p->textstycode, base);
}

/*
 *
 */
void
update_color(Color& c, const Color& c0) {
  mdreal rlnan = medusa::rnan();
  if(c0.red != rlnan) c.red = c0.red;
  if(c0.green != rlnan) c.green = c0.green;
  if(c0.blue != rlnan) c.blue = c0.blue;
  if(c0.opacity != rlnan) c.opacity = c0.opacity;

  if(c.red < 0.0) c.red = 0.0;
  if(c.green < 0.0) c.green = 0.0;
  if(c.blue < 0.0) c.blue = 0.0;
  if(c.opacity < 0.0) c.opacity = 0.0;

  if(c.red > 1.0) c.red = 1.0;
  if(c.green > 1.0) c.green = 1.0;
  if(c.blue > 1.0) c.blue = 1.0;
  if(c.opacity > 1.0) c.opacity = 1.0;
}
