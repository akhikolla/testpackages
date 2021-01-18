/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#ifndef scriptum_INCLUDED
#define scriptum_INCLUDED

#include <string>
#include <vector>
#include "medusa.h"

namespace scriptum {

  /* Color definition. */
  class Color {
  public:
    medusa::mdreal red;
    medusa::mdreal green;
    medusa::mdreal blue;
    medusa::mdreal opacity;
  public:
    Color();
    ~Color();

    /* Create a color from an XML string (with or without the '#'). */
    Color(const std::string&);

    /* Estimate perceived contrast between two colors.
       If the input is darker, returns a negative value. 
       If the input is lighter, returns a positive value. */
    medusa::mdreal contrast(const Color&) const;

    /* Return RGB color as a six-character hexadecimal code. */
    std::string hex() const;
  };

  /* Color functions. */
  extern Color colormap(const medusa::mdreal, const std::string&);
  extern std::vector<Color> colorize(const std::vector<medusa::mdreal>&,
				     const medusa::mdreal,
				     const std::string&);

  /*
   * Parameters for graphics primitives.
   */
  class Style {
  public:

    /* Determine if object can be pointed (e.g. by mouse). */
    bool pointable;
    
    /* Text-anchor: 'start', 'middle', 'end'. */
    std::string anchor;

    /* Rotation angle about the element origin in degrees[0, 360]. */
    medusa::mdreal angle;

    /* Fill color in RGB[0, 1]. */
    Color fillcolor;

    /* Font family. */
    std::string fontfamily;

    /* Font size (pixels). */
    medusa::mdreal fontsize;

    /* Font weight: 100 | 200 | 300 | 400 | 500 | 600 | 700 | 800 | 900. */
    medusa::mdsize fontweight;

    /* Identifier (max 255 characters). */
    std::string identity;

    /* Origin of the element coordinate system (pixels). */
    std::vector<medusa::mdreal> origin;

    /* Minimum distance to canvas margin. */
    medusa::mdreal padding;

    /* Stroke color in RGB[0, 1]. */
    Color strokecolor;

    /* Stroke width (pixels). */
    medusa::mdreal strokewidth;

    /* Optional values to attach to an element (max 255 characters
       each). These will appear as v0="content", v1="content",
       v2="content" et cetera in the SVG element. */
    std::vector<std::string> values;
    
  public:
    Style();
    ~Style();
  };

  /*
   * Create Scalable Vector Graphics code.
   */
  class Frame {
  private:
    void* buffer;
  public:
    Frame();

    /* Copy the contents from the argument.  */
    Frame(const Frame&);
    void operator=(const Frame&);

    /* Free resources. */
    ~Frame();

    /* Expand the bounding box to include a point. Style attributes
       are not considered. Returns false if point was unusable. */
    bool box(const medusa::mdreal, const medusa::mdreal);
    
    /* Draw curve or shape. The 1st argument sets the horizontal
       and the 2nd the vertical coordinates. If the first and last
       points are equal, the curve is drawn closed. */
    bool curve(const std::vector<medusa::mdreal>&,
	       const std::vector<medusa::mdreal>&);

    /* Draw a quadratic Bezier curve. The first two inputs set the
       starting point. The third and fourth inputs set the horizontal
       and vertical positions for the control point. The last two
       inputs set the end point. */
    bool curve(const medusa::mdreal&, const medusa::mdreal&,
	       const medusa::mdreal&, const medusa::mdreal&,
	       const medusa::mdreal&, const medusa::mdreal&);

    /* Return current graphics code and remove it from the object. */
    virtual std::string flush();

    /* Close the current group. Returns the number of remaining groups. */
    medusa::mdsize group();
    
    /* Open a new element group with the specified identity.
       Returns the current number of groups. */
    medusa::mdsize group(const std::string&);

    /* Coordinate range that contains elements of the frame. */
    virtual std::pair<medusa::mdreal, medusa::mdreal> horizontal() const;

    /* Draw a predefined shape. The first two arguments set the position
       and the third sets the radius. The fourth argument defines the shape
       itself and the optional rotation angle (e.g. 'circle', 'cross|50',
       'square|120', 'star|130' or 'triangle|200'). */
    bool shape(const medusa::mdreal, const medusa::mdreal,
	       const medusa::mdreal, const std::string&);

    /* Draw an elliptical slice. The first two arguments set the ellipse
       center (X and Y). The 3rd and 4th arguments set the smaller and larger
       radii for the slice. The 5th and 6th arguments set the smaller and
       larger angle (in degrees) for the slice. */
    bool slice(const medusa::mdreal, const medusa::mdreal,
	       const medusa::mdreal, const medusa::mdreal,
	       const medusa::mdreal, const medusa::mdreal);

    /* Return current style. */
    Style style() const;

    /* Set current style. Empty style data is reverted to previous values
       or system defaults if history is not available. */
    void stylize(const Style&);

    /* Write a string of text. */
    bool text(const medusa::mdreal, const medusa::mdreal,
	      const std::string&);

    /* Coordinate range that contains elements of the frame. */
    virtual std::pair<medusa::mdreal, medusa::mdreal> vertical() const;

    /* List of shape names available. */
    static std::vector<std::string> shapes();
  };

  /*
   *
   */
  class Artist {
  private:
    void* buffer;
  public:
    Artist();

    /* Create a new object with the specified output file. */
    Artist(const std::string&);

    /* Copy the contents from the argument. Objects in active
       rendering cannot be copied.  */
    Artist(const Artist&);
    void operator=(const Artist&);

    /* Free resources. */
    ~Artist();

    /* Set background color. */
    void background(const Color&);
    
    /* Finish rendering and return the number of bytes. The background
       will be filled according to the current style. The input is any
       custom inline content that will be added within the 'svg' or
       the HTML 'body' depending on the file format. */
    unsigned long close(const std::string&);

    /* Close the current group. Returns the number of remaining groups. */
    medusa::mdsize group();
    
    /* Open a new element group with the specified identity.
       Returns the current number of groups. */
    medusa::mdsize group(const std::string&);

    /* Coordinate range that contains rendered elements. */
    std::pair<medusa::mdreal, medusa::mdreal> horizontal() const;
    
    /* Render graphics from a frame. */
    bool paint(Frame&);
    
    /* Coordinate range that contains rendered elements. */
    std::pair<medusa::mdreal, medusa::mdreal> vertical() const;
  };
  
  /* Version information. */
  extern std::string version();
}

#endif /* scriptum_INCLUDED */
