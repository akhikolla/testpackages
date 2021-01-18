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

// TimeTree object constructor
// Will draw an svg container as nodes in the perimeter of a circle
// and edge between node
// using time dimension
function TimeTree(el,options) {

  //do options scoping with default attributes
  
  this.name = "timetree";
  
  // main element
  this.el = d3.select(el);
  if(options.width) this.el.style("width", options.width+"px")
  else this.el.style("width", "100%")
  if(options.height) this.el.style("height", options.height+"px")
  else this.el.style("height", "100%")

  this.el.classed("timetree-base",true);
  
  // date slider object 
  this.slider_date = null;
  
  // player control
  this.player_control = null;
  // animation timer
  this.timer = null;
  
  // svg elelment
  this.svg = null;
  
  // data (json)
  this.root = null;
  
  // list nodes ID
  this.nodes = [];
  
  // font size 
  this.nodes_fsize = "20";
  this.nodes_fsize_hover = "30px";
  
  
  if( options.nodes_color === undefined) { this.nodes_color = default_nodes_color;}
  else { this.nodes_color = options.nodes_color; }
  
  if( options.mouseclick === null ) { this.mouseclick = null;}
  else { this.mouseclick = options.mouseclick;}
  
  // offsprings colors
  this.edges_colors_offsprings_hex = d3.quantize(d3.interpolateCool,22)
  this.edges_colors_offsprings_hex.forEach( function(d,i,a){a[i]=d3.color(d).hex()})
  this.edges_colors_offsprings =  d3.scaleQuantize([0,1], this.edges_colors_offsprings_hex);
  // infected by colors
  this.edges_colors_infectedby_hex =  d3.quantize(d3.interpolateRgb('orange','red'),22);
  this.edges_colors_infectedby_hex.forEach(function(d,i,a){a[i]=d3.color(d).hex()})
  this.edges_colors_infectedby =  d3.scaleQuantize([0,1], this.edges_colors_infectedby_hex);
  
  // svg nodes group
  this.nodes_group = null;
  
  // svg g time graphs
  this.time_graphs = null;
  
  // time beginning
  this.starttime = 0;
  
  // time ending
  this.endtime = -1;
  
  // list all times available
  this.times = [];
  // is times is timestamp format
  this.time_timestamp = true;
  
  // svg width and height
  this.default_width = options.width;
  this.default_height = options.height;
  this.width = options.width;
  this.height = options.height;
  //console.log(options.width + " " + options.height)
  
  // circle center coordinates
  this.centerx = this.width/2;
  this.centery = this.height/2;
  
  // angle between nodes
  this.stepangle = 0;
  // raduis of the circle
  if(this.width < this.height) {
    this.radius = this.width/2.6;
  }
  else {
    this.radius = this.height/2.6;
  }
  
}

// default function to color nodes
var default_nodes_color = function(value) {
  
  if(value !== undefined) { return "red";}
  
  return "black";
}

// TimeTree loadOptions
TimeTree.prototype.loadOptions = function(options) {
  if( options.nodes_color === undefined) { this.nodes_color = default_nodes_color;}
  else { this.nodes_color = options.nodes_color; }
  
  if( options.mouseclick === null ) { this.mouseclick = null;}
  else { this.mouseclick = options.mouseclick;}
}

// TimeTree laod a json data object 
// data should be a json format object
// format :
// { "graphs" : {
//      "graph": {
//          "nodes" : [
//           { 
//              "id":0,
//              "status": "Inf"
//            },
//            { 
//              "id":1,
//              "status": "dob"
//            }
//          ],
//          "edges" : [
//            {
//              "id":1,
//              "source": 0,
//              "target": 1
//            }
//          ],
//          "time" : "adate"
//     }
//     ...
//   }
// }
//
TimeTree.prototype.loadJSON = function(data) {

  if(data !== undefined){
    delete this.root;
    this.root = data;
  }
  else { console.log("TimeTree data empty"); return; }

  var _this = this;
  
  //console.log(data);
  
  if( parseInt(new Date(parseInt(data[0].time)).getFullYear()) <= 1971) {
    console.log("Julian Day time step On");
    _this.time_timestamp = false;
  }
  else {_this.time_timestamp = true;}
  
  // set begining and ending time in milliseconde (json graph Have to be sort by time)
  this.starttime = parseInt(data[0].time);
  this.endtime = parseInt(data[data.length - 1].time);
  
  // load times array and nodes array
  delete _this.times;
  delete _this.nodes;
  _this.times = [];
  _this.nodes = [];
  data.forEach( function(graph) {
    _this.times.push(parseInt(graph.time));
    graph.nodes.forEach( function(node){
      if( _this.nodes.indexOf(node.ID) === -1) {
        _this.nodes.push(node.ID);
      }
    })
  })
  
  // stepangle depend of the number of nodes
  this.stepangle = (Math.PI*2) / this.nodes.length;
  
}

// createPlayerControl
// create a player control and animation over times 
// need slider_date and times to be set before.
// @param el : an element id for this player control
TimeTree.prototype.createPlayerControl = function(el) {
  var _this = this;
  
  // Add control icon as svg path icons
  if(d3.select("tt-icons").empty()) {
    _this.el.append("div")
      .html(
      '<svg id="timetree-svg-icons" display="none" version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="176" height="32" viewBox="0 0 176 32">'+
      '  <defs>'+
      '    <g id="icon-play"><path class="path1" d="M26.717 15.179l-13.698-8.486c-0.998-0.654-1.814-0.171-1.814 1.072v16.474c0 1.243 0.818 1.725 1.814 1.070l13.699-8.486c0 0 0.486-0.342 0.486-0.822-0.002-0.478-0.488-0.821-0.488-0.821z"></path></g>'+
      '    <g id="icon-pause"><path class="path1" d="M21.6 4.8c-1.59 0-2.88 0.49-2.88 2.080v18.24c0 1.59 1.29 2.080 2.88 2.080s2.88-0.49 2.88-2.080v-18.24c0-1.59-1.29-2.080-2.88-2.080zM10.4 4.8c-1.59 0-2.88 0.49-2.88 2.080v18.24c0 1.59 1.29 2.080 2.88 2.080s2.88-0.49 2.88-2.080v-18.24c0-1.59-1.29-2.080-2.88-2.080z"></path></g>'+
      '    <g id="icon-first"><path class="path1" d="M11.976 16c0 0.413 0.419 0.707 0.419 0.707l11.64 7.31c0.862 0.565 1.565 0.149 1.565-0.92v-14.195c0-1.070-0.702-1.486-1.565-0.922l-11.64 7.312c0 0.002-0.419 0.294-0.419 0.707zM6.4 8.571v14.858c0 1.421 0.979 1.856 2.4 1.856s2.4-0.435 2.4-1.854v-14.859c0-1.422-0.979-1.858-2.4-1.858s-2.4 0.437-2.4 1.858z"></path></g>'+
      '  </defs>'+
      '</svg>');
    
  }
  
  // Add the player control interface
  _this.player_control = _this.el.append("div").attr('class', 'player-control-container')
              .attr("id",el);
  _this.player_control.html(
      "<div class='step-back-control'><svg class='icon' width='32' height='32' viewBox='0 0 32 32'><use xlink:href='#icon-first'></use></svg></div>"+
      "<div class='play-back-control'><svg class='icon' width='32' height='32' viewBox='0 0 32 32'><g transform='rotate(180, 16, 16)'><use xlink:href='#icon-play'></use></g></svg></div>"+
      "<div class='pause-control'><svg class='icon' width='32' height='32' viewBox='0 0 32 32'><use xlink:href='#icon-pause'></use></svg></div>"+
      "<div class='play-forward-control'><svg class='icon' width='32' height='32' viewBox='0 0 32 32'><use xlink:href='#icon-play'></use></svg></div>"+
      "<div class='step-forward-control'><svg class='icon' width='32' height='32' viewBox='0 0 32 32'><g transform='rotate(180, 16, 16)'><use xlink:href='#icon-first'></use></g></svg></div>"
    );

  // Animation function
  // @param step : position step move start time forward, negative step move end time backward
  var playAnimation = function(step) {
    if(_this.timer !== null && _this.timer !== undefined) {
      _this.timer.stop();
      delete _this.timer;
    }
    
    // set the timer
    _this.timer = d3.interval(function(el) {
      var start = parseInt(_this.slider_date.get()[0]);
      var end = parseInt(_this.slider_date.get()[1]);
      
      // when buttons are at the same position
      //if( start === end) {
        if(step > 0 && end+step <= (_this.times.length - 1)) _this.slider_date.set([start+step,end+step]); 
        else if(step < 0 && start+step >= 0) _this.slider_date.set([start+step,end+step]);
        else _this.timer.stop();
        return(0);
      //}
      
      // change behavior : move button to overlap each other   
     /*if( step > 0 && (start + step) > end) {
        _this.timer.stop();
        return(0);
      } else if(step > 0) {
        _this.slider_date.set([start+step,null]);  
      }
      
      if( step < 0 && (end + step ) < start ) {
        _this.timer.stop();
        return(0);
      } else if(step < 0){
        _this.slider_date.set([null, end+step]);  
      }*/  
    },1000)
  }
  
  // stop animation
  var endAnimation = function() {
    if(_this.timer !== null && _this.timer !== undefined) {
      _this.timer.stop();
      delete _this.timer;
    }
  }
  
  // move start or end time to the min or max range
  // @param value : positive move end time, negative move start time
  var stepAnimation = function(value) {
    if(_this.timer !== null && _this.timer !== undefined) {
      _this.timer.stop();
      delete _this.timer;
    }
    if(value > 0) _this.slider_date.set([(_this.times.length - 1),(_this.times.length - 1)]);
    if(value < 0) _this.slider_date.set([0,0]);
  }
  
  // set event on player control
  _this.player_control.select('.step-back-control').on('click', function() { stepAnimation(-1); });
  _this.player_control.select('.play-back-control').on('click', function() { playAnimation(-1); });
  _this.player_control.select('.pause-control').on('click', function() { endAnimation(); });
  _this.player_control.select('.play-forward-control').on('click', function() { playAnimation(1); });
  _this.player_control.select('.step-forward-control').on('click', function() { stepAnimation(1); });
}

// Create a slider with two times handler
// @param el : element id name
TimeTree.prototype.createSliderDate = function(el) {
  
  var _this = this;
  
  // add div element + id
  _this.el.insert("div","svg").attr("id",el);
  
  _this.slider_date = noUiSlider.create(document.getElementById(el), {
    animate: true,
    animationDuration: 500,
    range: {
        min: [ 0 ],
        max: [ (_this.times.length - 1) ]
    },
    start: [ 0, (_this.times.length - 1) ],
    step: 1,
    pips: {
		  mode: 'range',
		  density: 1,
	    format: {
	      to: function(value) {
	        return _this.time_timestamp ? (new Date(_this.times[parseInt(value)]).toDateString()) : _this.times[parseInt(value)] ;
	      }
	    }
    },
    tooltips: [
      {
        to: function(value) {
          return _this.time_timestamp ? (new Date(_this.times[parseInt(value)]).toDateString()) : _this.times[parseInt(value)] ;
        }
      },{
        to: function(value) {
          return _this.time_timestamp ? (new Date(_this.times[parseInt(value)]).toDateString()) : _this.times[parseInt(value)] ;
        }
      }
      ]
  });
  
  // event update on slider -> update starttime and endtime
  // + throw event update on svg
  _this.slider_date.on("update", function(){
  
    _this.starttime = _this.times[parseInt(_this.slider_date.get()[0])];
    _this.endtime = _this.times[parseInt(_this.slider_date.get()[1])];
    _this.svg.dispatch("update");
    
  });
  
};


// Nodep mouseover 
// show node source and target
// show edges between them
TimeTree.prototype.mouseovernode = function(data) {
  
  var _this = this;
  
  color_sources = d3.interpolateReds();
  
  this.nodes_group.selectAll("text")
    .style("opacity",0.2)
    .attr("dx", function(d) { return (parseInt(d3.select(this).attr("dx")) <= 0) ? "-10px" : "10px"; })
    .style("font-size",_this.nodes_fsize+"px");

  this.nodes_group.select("#node-"+data)
    .style("opacity",1)
    .transition()
    .duration(1000)
    .attr("dx", function(d) { return (parseInt(d3.select(this).attr("dx")) <= 0) ? "-15px" : "15px"; })
    .style("font-size",_this.nodes_fsize_hover);
 
  this.root.forEach( function(graph) {
    if( parseInt(graph.time) >= _this.starttime && parseInt(graph.time) <= _this.endtime) {
      graph.edges.forEach( function(edge){
        if( edge.source == data ) {
          _this.time_graphs.select("#edge-"+edge.ID)
            .style("stroke", _this.edges_colors_offsprings(parseFloat(edge.weight)))
            .attr("marker-end", function(d){ return("url(#arrow_" + _this.edges_colors_offsprings(parseFloat(edge.weight)).substring(1) + ")");})
            .style("opacity",1)
            .raise();
            
          _this.nodes_group.select("#node-"+edge.target)
            .style("opacity",1)
            .transition()
            .duration(1000)
            .attr("dx", function(d) { return (parseInt(d3.select(this).attr("dx")) <= 0) ? "-15px" : "15px"; })
            .style("font-size",_this.nodes_fsize_hover);
          
          _this.time_graphs.selectAll("#tt-graph-"+ graph.time).raise();
        }
        if(edge.target == data) {
          _this.time_graphs.select("#edge-"+edge.ID)
            .style("stroke", _this.edges_colors_infectedby(parseFloat(edge.weight)))
            .attr("marker-end", function(d){ return("url(#arrow_" + _this.edges_colors_infectedby(parseFloat(edge.weight)).substring(1) + ")");})
            .style("opacity",1)
            .raise();
            
          _this.nodes_group.select("#node-"+edge.source)
            .style("opacity",1)
            .transition()
            .duration(1000)
            .attr("dx", function(d) { return (parseInt(d3.select(this).attr("dx")) <= 0) ? "-15px" : "15px"; })
            .style("font-size",_this.nodes_fsize_hover);
          
          _this.time_graphs.selectAll("#tt-graph-"+ graph.time).raise();
        }
        if(edge.source != data && edge.target != data) {
          _this.time_graphs.select("#edge-"+edge.ID)
            .attr("marker-end", function(d){ return("url(#arrow_4682b4)"); })
            .style("stroke","steelblue")
            .style("opacity",0.05);
          _this.time_graphs.selectAll("#tt-graph-"+ graph.time).lower();
        }
      })
    }
    else {
          _this.time_graphs.selectAll("#tt-graph-"+ graph.time).lower();

    }
  }) 
};

// mode mouse out 
// show default nodes and edges
TimeTree.prototype.mouseoutnode = function(data) {
  
  var _this = this;
  
  this.nodes_group.selectAll("text")
    .style("opacity",1)
    .transition()
    .attr("dx", function(d) { return (parseInt(d3.select(this).attr("dx")) <= 0) ? "-10px" : "10px"; })
    .style("font-size",_this.nodes_fsize+"px");
 
  _this.time_graphs.selectAll("path")
    .attr('marker-end', function(d){ return('url(#arrow_4682b4)'); })
    .style("stroke","steelblue")
    .style("opacity",1);
};

// show
TimeTree.prototype.show = function(){
  this.el.style("display",'block');
}

// hide
TimeTree.prototype.hide = function(){
  this.el.style("display",'none');
}


// Draw everythings
// From initial parameters and json object
// create date slider + svg geometry
TimeTree.prototype.draw = function() {
  
  var _this = this;

  _this.svg = _this.el.append("svg")
              .classed("timetree-svg",true);
  
  if( ! _this.slider_date ) {
    _this.createSliderDate("timetree-slider-date");
  }
  
  if( ! _this.player_control ) {
    _this.createPlayerControl("timetree-player-control");
  }
  
  // remove top height of SLideDate from svg height
  _this.height = _this.height - document.getElementById("timetree-slider-date").offsetHeight - 50;
  _this.svg.style("height",_this.height+50+"px")
            .style("width",_this.width+"px");
  
  // set center coordinates
  _this.centerx = _this.width/2;
  _this.centery = _this.height/2;
  
  // circle radius
  if(this.width < this.height) {
    this.radius = this.width/2.6;
  }
  else {
    this.radius = this.height/2.6;
  }
  
  _this.nodes_fsize = Math.min(Math.ceil( _this.radius * Math.sin(_this.stepangle)), 18);
  _this.nodes_fsize_hover = _this.nodes_fsize + _this.nodes_fsize/2 + "px";
  
  // Append a g with all nodes
  _this.nodes_group = this.svg.append("g")
                    .attr("id","nodes");
  
  _this.nodes_group.selectAll("g")
          .data(_this.nodes)
          .enter()
          .each( function(d,i) {
            
            // Node angle
            angle = _this.stepangle*i;
            
            // Node coordinates
            x = _this.centerx + (_this.radius)*Math.cos(angle);
            y = _this.centery + (_this.radius)*Math.sin(angle);
            // Add Node Id as text and set position
            d3.select(this)  
              .append("text")
              .attr("id", function(d) { return "node-" + d; })
              .on("mouseover", function(d) { _this.mouseovernode(d);})
              .on("mouseout", function(d) { _this.mouseoutnode(d);})
              .on("click", function(d) { _this.mouseclick(d);})
              .attr("x",function(d){ return x;})
              .attr("y",function(d){ return y;})
              .attr("dx", function(d) { return (angle > Math.PI/2 && angle < (3*Math.PI)/2) ? "-10px" : "10px"; })
              .attr("transform", function(d) {
                return "rotate(" + ((angle*180)/Math.PI) + "," + x + "," + y + ")"
                + ((angle > Math.PI/2 && angle < (3*Math.PI)/2)  ? "rotate(180," + x + "," + y + ")" : "")
              })
              .attr("text-anchor", function(d) { return (angle > Math.PI/2 && angle < (3*Math.PI)/2) ? "end" : "start"; })
              .text(function(d) { return d})
              .style("font-size",_this.nodes_fsize+"px")
              .classed("timetree-node",true);
          });
  
  // Create a line generator
  var line = d3.line()
    .curve(d3.curveBundle.beta(0.8))
    .x(function(d){return d[0];})
    .y(function(d){return d[1];});
  
  //create a path from source to target nodes trough the center
  var createPath = function(edges) {
    path = []
    
    sourcex = _this.nodes_group.select("#node-"+edges.source).attr("x");
    sourcey = _this.nodes_group.select("#node-"+edges.source).attr("y");
    
    targetx = _this.nodes_group.select("#node-"+edges.target).attr("x");
    targety = _this.nodes_group.select("#node-"+edges.target).attr("y");

    path.push([sourcex,sourcey])
    path.push([_this.centerx,_this.centery])
    path.push([targetx,targety])
    
    return path;
  }

  // add defs and marker for arrows colors
  defs = _this.svg.append("defs");
    defs.selectAll("marker")
    .data( ["#4682b4"].concat(_this.edges_colors_offsprings_hex, _this.edges_colors_infectedby_hex) )
    .enter()
    .append("marker")
		.attr("id", function(d){ return("arrow_"+d.substring(1)) })
		.attr("viewBox","0 0 10 10")
		.attr("refX",0)
		.attr("refY",5)
		.attr("markerUnits","strokeWidth")
		.attr("markerWidth",4)
		.attr("markerHeight",4)
		.attr("orient","auto")
		.append("path")
		.attr("d", "M 0 0 L 10 5 L 0 10 z")
		.attr("fill", function(d){ return(d); });

  // create a g for each time to display
  _this.time_graphs = _this.svg.selectAll("g")
          .enter()
          .data(_this.root)
          .enter()
          .append("g")
          .attr("id", function(d) {return "tt-graph-"+ d.time;});
          
  // create edge fior this times        
  _this.time_graphs.selectAll('path')
              .data(function(d){return d.edges;})
              .enter()
              .append("path")
              .each( function(d,i){
                d3.select(this)
                .attr("id",function(d){return "edge-" + d.ID;})
                .attr("d", line(createPath(d)))
                .attr("marker-end","url(#arrow_4682b4)")
                .classed("timetree-edge",true)
              });

  
  // update color and visible edges
  _this.update(); 
   
  _this.svg.on("update", function(){
    _this.update();
  });

};

// update svg container
// active nodes and edges between start and end time
TimeTree.prototype.update = function() {
  
  var _this = this;
  
  // default text color
  _this.nodes_group.selectAll("text")
    .style("fill",_this.nodes_color("default"))
    .style("pointer-events","none");
  
  // color nodes status from begining until endtime
  this.root.forEach( function(graph){
    if( parseInt(graph.time) <= _this.endtime ) {
      graph.nodes.forEach( function(node){
        _this.nodes_group.select("#node-"+node.ID).style("fill",_this.nodes_color(node.status))
                                                  .style("pointer-events","auto");
      })
    }
    
  });
  
  // display edges between starttime and endtime  
  _this.times.forEach( function(t) {
    if( _this.starttime <= t && t <= _this.endtime ) {
      _this.svg.select("#tt-graph-"+t).style("display","block");
      
    }
    else{
      _this.svg.select("#tt-graph-"+t).style("display","none");
    }
  });
}

// Redraw an existing timetree with new width and height
TimeTree.prototype.redraw = function(width, height) {
  var _this = this;
  _this.hide();
  
  if( width !== undefined) _this.width = width;
  else _this.width =  _this.default_width;
  if( height !== undefined) _this.height = height;
  else _this.height = _this.default_height;
  
  // set default size to the new size
  _this.default_width = _this.width;
  _this.default_height = _this.height;
  

  _this.el.style("width", _this.default_width+"px");
  _this.el.style("height", _this.default_height+"px");
  
  delete this.slider_date;
  _this.slider_date = null;
  d3.select("#timetree-slider-date").remove();

  // player control
  delete this.player_control;
  _this.player_control = null;
  if(_this.timer) _this.timer.stop();
  delete _this.timer;
  d3.select("#timetree-player-control").remove();
  d3.select("#timetree-svg-icons").remove();

  
  this.nodes_group.remove();
  this.time_graphs.remove();
  _this.svg.remove();

  _this.draw();
  _this.show();
}
