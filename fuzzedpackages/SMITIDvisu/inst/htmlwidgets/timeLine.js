/* This file is part of SMITIDvisu package.
* Copyright (C) 2018-2019 Jean-François Rey <jean-francois.rey@inra.fr>

* SMITIDvisu is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* SMITIDvisu is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with SMITIDvisu. If not, see <https://www.gnu.org/licenses/>.
*/

HTMLWidgets.widget({

  name: 'timeLine',

  type: 'output',

  factory: function(el, width, height) {

    // a TimeLine object
    var tl = null;
    // options
    var tlOptions = [];

    return {

      renderValue: function(x) {
        
        tlOptions["width"] = width;
        tlOptions["height"] = height;
        if( x.options.color !== undefined) tlOptions["color"] =x.options.color;
        
        var mouseclick = null;
        if (HTMLWidgets.shinyMode) {
          
          mouseclick = function(d) {
            console.log( d);
            Shiny.onInputChange(
                "timeline_click",
                {label: d.label, id: d.ID, time: d.timestart}
                //{id: d, nonce: Math.random()}
              );
          }
          
          tlOptions["mouseclick"] = mouseclick;
        }
        
        if(tl === null) {
          tl = new TimeLine(el, tlOptions);
          tl.loadJSON(x.data);
          tl.draw();
        }
        else {
          tl.hide();
          tl.loadJSON(x.data);
          tl.redraw(width, height);
        }
        
        tl.show();

      },

      resize: function(width, height) {

        if(tl) { tl.redraw(width,height);  }

      },
      
      // expose timeLine to the outside
      getTimeLine: function() {
        return tl;
      }
      
    };
  }
});

/** getTimeLine
 * get timeLine htmlwidgetinstance
 * @apram id a timeline instance id
 */
function getTimeLine(id){
  
  var htmlWidgetsObj = HTMLWidgets.find("#" + id);
  var tl = htmlWidgetsObj.getTimeLine();
  return(tl);
}

/** updatetl
 * add message handler to update tiemline
 */
if(Shiny !== undefined) {
  Shiny.addCustomMessageHandler("updatetl", function(message) {
    var tl = getTimeLine(message.id);
    tl.hide();
    tl.loadJSON(message.data);
    tl.redraw();
  });
}
