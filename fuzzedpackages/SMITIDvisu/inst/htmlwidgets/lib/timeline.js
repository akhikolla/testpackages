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

// TimeLine object constructor
// will draw a time line 
function TimeLine(el, options) {
  
  if( ! this instanceof TimeLine) {
    return new TimeLine(el, options);
  }
  
  this.name = "timeline";

  // main element
  this.el = d3.select(el);
  if(options.width) this.el.style("width", options.width+"px");
  else this.el.style("width", "100%");
  if(options.height) this.el.style("height", options.height+"px");
  else this.el.style("height", "100%");
  
  this.el.classed("timeline-base", true);
  // data (json)
  this.root = null;
  // level array and collision by id (same timestart -> x)
  this.levelCollision = null;
  
  // Legend
  this.legend = null;
  this.legend_height = 50;
  
  // width and height
  this.width = options.width;
  this.height = options.height - this.legend_height;
  
  if( options.mouseclick === null ) { this.mouseclick = null;}
  else { this.mouseclick = options.mouseclick;}
  
  // margin
  this.margin = 30;
  
  // number of level to display
  this.nbLevels = 1;
  // default level height
  this.levelHeight = 30;
  // position of the top of level 0
  this.level0_top = null;
  
  // a Scale function for Date
  this.xScale = null;
  // x Axis constructor
  this.xAxis = null;
  // estimated axis height size
  this.xAxisSize = 40;
  
  // color function 
  this.color = (options.color === undefined) ? {"default":"black"} : options.color;
  
  // pointer on a div tooltip
  this.tooltip = null;
  
  // Day format
  this.timeFormat = d3.utcFormat("%Y-%m-%dT%H:%M:%SZ");
  
  // is time timestamp
  this.time_timestamp = true;
  
  

}


// getColor
// return a color for an associated name
TimeLine.prototype.getColor = function(name) {
  
  return this.color[name] !== undefined ? this.color[name] : this.color["default"];
}


// loadJSON 
// load a json data format
//
// { "type" : "timeline",
//    mintine : "Date",
//    maxtime: "Date",
//    nblevels: "nb",
//    { "level": "nb",
//      "label": "alabel",
//      "timeline" : [
//          { "label": "alabel",
//              "ID": "ID",
//              "timestart": "Date"
//              "timeend": "Date"
//            },
//          { "label": "alabel",
//              "ID": "ID",
//              "timestart": "Date"
//              "timeend": "Date"
//            }
//          ]
//    },
//    etc...
// }
TimeLine.prototype.loadJSON = function(data) {

  _this = this;
  this.root = data;
  
  if( data.length === 0 ){return(true);}
  
  //console.log(data);
  
  if( parseInt(new Date(parseInt(data[0].mintime)).getFullYear()) < 1971) {
    
    console.log("Julian Day time step On");
    this.time_timestamp = false;
  }

  this.minTime = parseInt(this.root[0].mintime);
  this.maxTime = parseInt(this.root[0].maxtime);
  this.nbLevels = parseInt(this.root[0].nblevels);
  
  // levels Height depend of the number of levels
  if(this.nbLevels < 2) this.levelHeight = (this.height - this.xAxisSize - this.margin) / 2;
  else this.levelHeight = (this.height - this.xAxisSize - this.margin) / this.nbLevels;
  
  // levels position depend of the number of levels
  this.level0_top = 0;
  if( this.nbLevels === 3 ) this.level0_top = this.levelHeight;
  if( this.nbLevels === 2 && (parseInt(this.root[1].level[0]) === -1 || parseInt(this.root[2].level[0]) === -1) ) {
    this.level0_top = this.levelHeight;
  }
  //this.level0_top = (this.height - this.xAxisSize)/2 - this.levelHeight/2;
  
  // create an array 2D by level to mark timeline with the same timestart
  this.levelCollision = {};
  this.root.forEach( function(line) {
    if(line.level === undefined) return(false);
    _this.levelCollision[line.level] = {};
    line.timeline.forEach( function(tl) {
      _this.levelCollision[line.level][tl._row] = {}
      line.timeline.forEach( function(tl2) {
        if(tl.timestart === tl2.timestart) _this.levelCollision[line.level][tl._row][tl2._row] = 1;
      })
    })
  })
}

// Add color legend to an element
TimeLine.prototype.drawLegend = function(el) {
  
  var _this = this;
  
  var divL = "";
  Object.keys(_this.color).forEach( function(k,i) {
    divL += '<div><span class="tl-legend-square" style="background: '+_this.color[k]+'"></span><span class="tl-legend-text">'+k+'</span></div>';
  })
  
  el.html(divL);
}

// Draw the object timeline
// from initial parameters and json data
TimeLine.prototype.draw = function() {
  
  var _this = this;
  
  if(_this.root.length === 0) { return(false);}
  
  _this.legend = _this.el.append("div").attr("id","tl-legend");
  _this.drawLegend(_this.legend);

  _this.svg = _this.el.append("svg")
                      .attr("id","tl-svg")
                      .style("height", _this.height + "px")
                      .style("width", _this.width + "px");
                      //.style("top", _this.legend_height + "px");

  // init scale with mintime and maxtime
  _this.xScale = d3.scaleTime()
                   .domain([_this.minTime, _this.maxTime])
                   .range([_this.margin,_this.width - _this.margin]);

  // init axis
  if(_this.time_timestamp) {
    _this.xAxis = d3.axisBottom(_this.xScale)
                    .tickFormat(d3.timeFormat("%Y-%m-%d"))
                    .ticks(d3.timeWeek.every(_this.width/(_this.maxTime - _this.minTime)));
  } else {
    _this.xAxis = d3.axisBottom(_this.xScale)
                    .tickFormat(d3.format(".2f"))
                    .ticks(_this.width/(_this.maxTime - _this.minTime));
  }
  
  // add a div for tooltip	
  this.tooltip = this.el.append("div")
            .classed("tl-tooltip",true)
            .style("opacity",0);

  // add defs for arrows
  this.svg.append("defs")
    .append("marker")
		.attr("id","tl-arrow")
		.attr("viewBox","0 0 10 10")
		.attr("refX",0)
		.attr("refY",5)
		.attr("markerUnits","strokeWidth")
		.attr("markerWidth",4)
		.attr("markerHeight",4)
		.attr("orient","auto")
		.append("path")
		.attr("d", "M 0 0 L 10 5 L 0 10 z");

  // draw timeline
  var middle = _this.svg.selectAll("g")
                    .data(_this.root)
                    .enter()
                    .filter( function(d) { return d.level !== undefined ;})
                    .append("g")
                    .each( function(level) {
                      elmt = this;
                      // draw ech level
                      level.timeline.forEach( function(tl){ _this.drawThings(elmt, level.level, tl); });
                      // Add title from level0 label
                      if(parseInt(level.level) === 0) {
                        _this.el.select("div#tl-legend").insert("div","div")
                                                     .html("Host " + level.label)
                                                     .classed("tl-legend-title",true);
                      }
                    });


  // Add axis and rotate label
  _this.svg.append("g")
    .attr("transform","translate(0,"+ (_this.height - _this.xAxisSize - _this.margin) +")")
    .call(_this.xAxis)
    .selectAll("text")	
    .style("text-anchor", "end")
    .attr("dx", "-.8em")
    .attr("dy", ".15em")
    .attr("transform", function(d) {
        return "rotate(-65)" 
    });

}

// redraw everything
// using new width and height
TimeLine.prototype.redraw = function(width, height) {

  var _this = this;
  this.hide();
  if( width !== undefined ) this.width = width;
  if( height !== undefined ) this.height = height - _this.legend_height;
  
  if(this.nbLevels < 2) this.levelHeight = (this.height - this.xAxisSize - this.margin) / 2;
  else this.levelHeight = (this.height - this.xAxisSize - this.margin) / this.nbLevels;
  
  this.level0_top = 0;
  if( this.nbLevels === 3 ) this.level0_top = this.levelHeight;
  if( this.nbLevels === 2 && (parseInt(this.root[1].level[0]) === -1 || parseInt(this.root[2].level[0]) === -1) ) {
    this.level0_top = this.levelHeight;
  }
  
  if(_this.legend) _this.legend.remove();
  if(_this.svg) _this.svg.remove();
  if(_this.tooltip) _this.tooltip.remove();

  _this.draw();
  _this.show();
}

// draw an items of the timeline (circle or rectangle)
TimeLine.prototype.drawThings = function(el, level, line) {
  
  var _this = this;
  
  // give Y of top level position
  var getLevelY = function(level) {
    if(parseInt(level) === 0 ) { return _this.level0_top; }
    return _this.level0_top + (_this.levelHeight * parseInt(level));
  }
  
  // draw rectangle
  if( line.timeend !== undefined && line.timeend !== "") {
    var x = _this.xScale(line.timestart);
    var y = getLevelY(level);
    
    var large;
    if( line.timeend !== "Inf" ) {large = _this.xScale(line.timeend) - x;}
    else {large = _this.xScale(_this.maxTime) - x; }

    d3.select(el)
      .append("rect")
      .attr("x", x)
      .attr("y", y )
      .attr("width", large )
      .attr("height",_this.levelHeight)
      .style("fill",_this.getColor(line.label));
    d3.select(el)
      .append("text")
      .attr("x",x+(large/2))
      .attr("y",y+_this.levelHeight-20);
      
    if( line.timeend !== "Inf") {
      // add line at the end of the rect
      d3.select(el)
        .append("line")
        .attr("x1",x+large)
        .attr("x2",x+large)
        .attr("y1",y)
        .attr("y2",getLevelY(level) + _this.levelHeight)
        .attr("stroke-dasharray","15,10,5,10")
        .classed("tl-line-rect",true)
        .on("mouseover", function() {
          if( _this.time_timestamp ) {
            _this.tooltip.html("time: " + _this.timeFormat(new Date(line.timeend)) + "<br/>" + line.label + "<br/>" + line.ID)
          }
          else {
            _this.tooltip.html("time: " + line.timeend + "<br/>" + line.label + "<br/>" + line.ID)
          }
          _this.tooltip.style("left", (d3.select(this).attr("x1")) + "px")
                       .style("top", (d3.select(this).attr("y1")) + "px")
                       .transition()
                       .duration(200)
                       .style("opacity",0.9);
        })
        .on("mouseout", function() {
          _this.tooltip.transition()
             .duration(500)
             .style("opacity",0);
        });
    }
  }
  //draw circle
  else {
    var x = _this.xScale(line.timestart);
    var y = getLevelY(level) + _this.levelHeight/2;
    var r = _this.levelHeight/6;
  
    
    // check for collision -> vertical align circles
    if( _this.levelCollision[level][line._row] !== undefined && Object.keys(_this.levelCollision[level][line._row]).length > 1) {
      r = _this.levelHeight/(3 * Object.keys(_this.levelCollision[level][line._row]).length);
      y = getLevelY(level) + ((Object.keys(_this.levelCollision[level][line._row]).indexOf(line._row)+1)*2)*r;
    }

    // an arc for textpath
    var c_arc = d3.arc()
                .innerRadius(r/2)
                .outerRadius(r/2)
                .startAngle(-Math.PI)
                .endAngle(Math.PI);
    
    // items info
    d3.select(el)
      .append("circle")
      .attr("cx",x)
      .attr("cy",y)
      .attr("r",r)
      .style("fill",_this.getColor(line.label))
      .on("mouseover", function() {
        if( _this.time_timestamp ) {
          _this.tooltip.html("time: " + _this.timeFormat(new Date(line.timestart)) + "<br/>" + line.label + "<br/>" + line.ID)
        } else {
          _this.tooltip.html("time: " + line.timestart + "<br/>" + line.label + "<br/>" + line.ID)
          
        }
        _this.tooltip.style("left",d3.select(this).attr("cx") + "px")
                     .style("top", d3.select(this).attr("cy") + "px");
        _this.tooltip.transition()
           .duration(200)
           .style("opacity",0.9);
      })
      .on("mouseout", function() {
        _this.tooltip.transition()
           .duration(500)
           .style("opacity",0);
      })
      .on("click", function() { _this.mouseclick(line)});
      
    d3.select(el)
      .append("path")
      .attr("id","tl-path-text"+line.ID)
      .attr("d", c_arc())
      .attr("transform","translate("+x+","+y+")")
      .classed("tl-path-text",true);
      
    d3.select(el)
      .append("text")
      .classed("tl-text",true)
      .style("font-size",r/2+"px")
      .append("textPath")
	    .attr("text-anchor","middle")
	    .attr("startOffset","25%")
      .attr("xlink:href", "#tl-path-text"+line.ID)
      .text(line.ID);
      
    
    // add lines and arrows depends of the circle position in the level
    if(parseInt(level) !== 0) {
      if( (parseInt(level) === -1 && (Object.keys(_this.levelCollision[level][line._row]).indexOf(line._row) === 0 ))
      || (parseInt(level) === 1 && (Object.keys(_this.levelCollision[level][line._row]).indexOf(line._row) === (Object.keys(_this.levelCollision[level][line._row]).length - 1) ) ) ) {
      d3.select(el)
        .append("line")
        .attr("x1",x)
        .attr("x2",x)
        .attr("y1",y+r)
        .attr("y2",getLevelY(level) + _this.levelHeight - 10)
        .attr("transform",parseInt(level) < 0 ? "" : "translate(0,-" +( y+r - getLevelY(level)) +")")
        .attr("marker-end","url(#tl-arrow)")
        .classed("tl-line-event",true);
      }       
    }
    
  }
}

// show
TimeLine.prototype.show = function(){
  this.el.style("display",'block');
}

// hide
TimeLine.prototype.hide = function(){
  this.el.style("display",'none');
}
