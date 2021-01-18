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

  name: 'mstVariant',

  type: 'output',

  factory: function(el, width, height) {
    
    // a ForceGraph object
    var mstV = null;
    // otpions
    var mstVOptions = {};

    return {

      renderValue: function(x) {
      
        mstVOptions["width"] = width;
        mstVOptions["height"] = height;
        mstVOptions["data"] = x.data;

        if(mstV === null){
          mstV = new ForceGraph(el,mstVOptions);
          mstV.draw();
        }
        else {
          mstV.hide();
          mstV.setJSON(x.data);
          mstV.redraw(width, height);
        }
        
        mstV.show();
      },

      resize: function(width, height) {
        if(mstV){ mstV.redraw( width, height);}
      },
      
      getMstVariant: function() {
        return(mstV);
      }
    };
  }
});

/** getMstVariant
* return object from id
* @param id : widget name (output name in server)
*/
function getMstVariant(id){
  
  var htmlWidgetsObj = HTMLWidgets.find("#" + id);
  var mstV = htmlWidgetsObj.getMstVariant();
  return(mstV);
}

/** updatemstV
 * add a custom message handler to update mstVariant
 * client side
 */
if(Shiny !== undefined) {
  Shiny.addCustomMessageHandler("updatemstV", function(message) {
    var mstV = getMstVariant(message.id);
    mstV.hide();
    //console.log(HTMLWidgets.dataframeToD3(message.data));
    mstV.setJSON(message.data);
    mstV.redraw();
  });
}
