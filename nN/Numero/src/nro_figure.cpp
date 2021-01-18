/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "nro.h"

/*
 *
 */
class SVGFrame : public Frame {
private:
  string code;
  vector<mdreal> bbox;
public:
  SVGFrame() : Frame() {};
  SVGFrame(const string& c, const vector<mdreal>& b) : Frame() {
    this->code = c;
    this->bbox = b;
  };
  string flush() {return code;};
  pair<mdreal, mdreal> horizontal() const {
    pair<mdreal, mdreal> res(0.0, 0.0);
    if(bbox.size() != 4) return res;
    res.first = bbox[0];
    res.second = bbox[2];
    return res;
  };
  pair<mdreal, mdreal> vertical() const {
    pair<mdreal, mdreal> res(0.0, 0.0);
    if(bbox.size() != 4) return res;
    res.first = bbox[1];
    res.second = bbox[3];
    return res;
  };
};

/*
 *
 */
RcppExport SEXP
nro_figure(SEXP fname_R, SEXP data_R, SEXP bbox_R, SEXP script_R) {
  string fname = as<string>(fname_R);
  vector<string> data = as<vector<string> >(data_R);
  string script = as<string>(script_R);

  /* Make sure bounding box is the right size. */
  vector<mdreal> bbox = nro::vector2reals(bbox_R);
  bbox.resize(4, 0.0);
  
  /* Open output file. */
  Artist art(fname);

  /* Use the derived class to pass data to renderer. */
  for(mdsize i = 0; i < data.size(); i++) {
    SVGFrame frame(data[i], bbox);
    art.paint(frame);
  }
  
  /* Return file size. */
  List output;
  unsigned long nbytes = art.close(script);
  output.push_back(long2string(nbytes), "nbytes");
  output.push_back(long2text(nbytes), "text");
  return output;
}
