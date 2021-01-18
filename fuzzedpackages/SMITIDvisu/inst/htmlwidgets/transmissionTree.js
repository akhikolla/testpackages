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

  name: 'transmissionTree',

  type: 'output',

  factory: function(el, width, height) {

    var tt = null;

    return {

      renderValue: function(x) {
        
        // color function for nodes
        if( x.options.nodes_color !== undefined) {
          var tt_nodes_color = function(value) {
            
            if( x.options.nodes_color[value] !== undefined ) {
              return x.options.nodes_color[value];
            }
            
            return x.options.nodes_color["default"] !== undefined ? x.options.nodes_color["default"] : "black";
          }
        }
        
        var mouseclick = null;
        if (HTMLWidgets.shinyMode) {

          mouseclick = function(d) {
            //console.log("click " + d);
            Shiny.onInputChange(
                "transmissionTree_node_click", d
                //{id: d, nonce: Math.random()}
              );
          }
        }
        
        if(tt === null) {
          tt = new TimeTree(el, {width: width, height: height, nodes_color: tt_nodes_color, mouseclick: mouseclick});
          tt.loadJSON(x.data);
          tt.draw();
        }else{
          tt.loadJSON(x.data);
          tt.loadOptions({nodes_color: tt_nodes_color, mouseclick: mouseclick});
          tt.redraw();
        }
        tt.show();
      },

      resize: function(width, height) {

        if(tt) { tt.redraw(width,height);  }

      },

      getTransmissionTree: function() {
        return(tt);
      }
    };
  }
});

/** getTransmissionTree
* return object from id
* @param id : widget name (output name in server)
*/
function getTransmissionTree(id){

  var htmlWidgetsObj = HTMLWidgets.find("#" + id);
  var tt = htmlWidgetsObj.getTransmissionTree();
  return(tt);
}

/** updateTT
 * add a custom message handler to update transmissionTree
 * client side
 */
if(Shiny !== undefined) {
  Shiny.addCustomMessageHandler("updateTT", function(message) {
    var tt = getTransmissionTree(message.id);
    tt.hide();
    //console.log(HTMLWidgets.dataframeToD3(message.data));
    tt.loadJSON(message.data);
    tt.loadOptions(message.options);
    tt.redraw();
  });
}

