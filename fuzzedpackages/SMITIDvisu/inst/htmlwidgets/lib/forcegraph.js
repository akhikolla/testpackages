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

// ForceGraph
// Draw a graph using force simulation
// add tiny node on links using link weight
// node radius depend of weight
function ForceGraph(el,options) {

  // main frame
  this.el = d3.select(el);
  if(options.width) this.el.style("width", options.width+"px");
  else this.el.style("width", "100%");
  if(options.height) this.el.style("height", options.height+"px");
  else this.el.style("height", "100%");


  // data to draw as json format                             
  this.root = options.data;
  //console.log(this.root);

  // svg container
  this.svg = null;
  
  // main container
  this.container = null;
  
  // control container
   this.control = null
  
  // element width and height
  this.width = options.width;
  this.height = options.height;
  
  // radius max size
  if(this.width < this.height) {
    this.radius = this.width/4;
  }
  else {
    this.radius = this.height/4;
  }
  // raduis min value
  //this.minRadius = 10.0;
  
  // array containing nodes and link information 
  // use by force simulator
  this.graph = {};
  this.graph.nodes = [];
  this.graph.links = [];
  
  // link (line) selector
  this.link = null;
  
  // node (circle) selector
  this.node = null;

  // force simulation object
  this.simulation = null;
  
  // bool to active simulation (default true)
  this.simulationIsActive = true;
  
  // tiny node on link selector
  this.linkNode = null;
  
  // the zoom
  this.zoom = null;
  this.zoomTranslate = [0,0];
  this.zoomScale = 1;
  
  ////////
  //Pie
  this.havePie = false;
  // Pie container
  this.nodePie = null;
  
  // Time informations
  this.minTime = Infinity;
  this.maxTime = 0;
  // Day format
  this.timeFormat = d3.timeFormat("%Y-%m-%d");
  // is time timestamp
  this.time_timestamp = false;
  
  // fill graph with root data
  if(this.root !== undefined) {
    makeGraph(this.graph, this.root);
    
    if(this.root.node_prop !== undefined && this.root.node_prop !== null) this.addPie();
  }
  
}

// setJSON set root attribut with json data 
// load data
// {
//  edges: [
//      {"source": sourceId,
//        "target": targetId,
//        "weight":weight
//      },
//      ...
//    ],
//    nodes: [
//      {"id": id,
//        "weight": weight
//      }
//      ...
//    ],
//    node_prop: [
//      { "time" : []
//        "value" : []
//      }
//      ...
//    ]
//}
ForceGraph.prototype.setJSON = function(data) {
  
  if(data !== undefined){
    delete this.root;
    this.root = data;
  } 
  if(this.root === undefined) return;
  
  this.havePie = false;
  delete this.graph.nodes;
  delete this.graph.links;
  delete this.graph;
  this.graph = {};
  this.graph.nodes = [];
  this.graph.links = [];
  
  if(this.root !== undefined && this.root !== null) makeGraph(this.graph, this.root);
  if(this.root.node_prop !== undefined && this.root.node_prop !== null) this.addPie();
  
};

// createControl
// create a control container
// @param el : an element id for this control
ForceGraph.prototype.createControl = function(el) {
  var _this = this;
  
  // Add control icon as svg path icons
  if(d3.select("fg-icons").empty()) {
    _this.el.append("div")
      .html(
      '<svg id="forcegraph-svg-icons" display="none" version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="176" height="32" viewBox="0 0 176 32">'+
      '  <defs>'+
      '    <g id="icon-play"><path d="M26.717 15.179l-13.698-8.486c-0.998-0.654-1.814-0.171-1.814 1.072v16.474c0 1.243 0.818 1.725 1.814 1.070l13.699-8.486c0 0 0.486-0.342 0.486-0.822-0.002-0.478-0.488-0.821-0.488-0.821z"></path></g>'+
      '    <g id="icon-pause"><path d="M21.6 4.8c-1.59 0-2.88 0.49-2.88 2.080v18.24c0 1.59 1.29 2.080 2.88 2.080s2.88-0.49 2.88-2.080v-18.24c0-1.59-1.29-2.080-2.88-2.080zM10.4 4.8c-1.59 0-2.88 0.49-2.88 2.080v18.24c0 1.59 1.29 2.080 2.88 2.080s2.88-0.49 2.88-2.080v-18.24c0-1.59-1.29-2.080-2.88-2.080z"></path></g>'+
      '    <g id="icon-camera"><path d="M10,6.536c-2.263,0-4.099,1.836-4.099,4.098S7.737,14.732,10,14.732s4.099-1.836,4.099-4.098S12.263,6.536,10,6.536M10,13.871c-1.784,0-3.235-1.453-3.235-3.237S8.216,7.399,10,7.399c1.784,0,3.235,1.452,3.235,3.235S11.784,13.871,10,13.871M17.118,5.672l-3.237,0.014L12.52,3.697c-0.082-0.105-0.209-0.168-0.343-0.168H7.824c-0.134,0-0.261,0.062-0.343,0.168L6.12,5.686H2.882c-0.951,0-1.726,0.748-1.726,1.699v7.362c0,0.951,0.774,1.725,1.726,1.725h14.236c0.951,0,1.726-0.773,1.726-1.725V7.195C18.844,6.244,18.069,5.672,17.118,5.672 M17.98,14.746c0,0.477-0.386,0.861-0.862,0.861H2.882c-0.477,0-0.863-0.385-0.863-0.861V7.384c0-0.477,0.386-0.85,0.863-0.85l3.451,0.014c0.134,0,0.261-0.062,0.343-0.168l1.361-1.989h3.926l1.361,1.989c0.082,0.105,0.209,0.168,0.343,0.168l3.451-0.014c0.477,0,0.862,0.184,0.862,0.661V14.746z"></path></g>'+
      '  </defs>'+
      '</svg>');
    
  }
  
  // Add the control interface
  _this.control = _this.el.append("div").attr('class', 'fg-control-container')
              .attr("id",el);
  _this.control.html(
    "<div class='fg-save-control'><svg class='icon' width='32' height='32' viewBox='0 0 20 20'><use xlink:href='#icon-camera'></use></svg></div><br/>"+
      "<div class='fg-pause-control'><svg class='icon' width='32' height='32' viewBox='0 0 32 32'><use xlink:href='#icon-pause'></use></svg></div>"+
      "<div class='fg-play-control' style='display:none;'><svg class='icon' width='32' height='32' viewBox='0 0 32 32'><use xlink:href='#icon-play'></use></svg></div>");

  
  // Save SVG as Image
  saveSVG = function() {
    try {
        var isFileSaverSupported = !!new Blob();
    } catch (e) {
      _this.control.select('.fg-save-control').style('display', 'none');
        //alert("blob not supported");
    }
    
    function save( dataBlob ){
		  saveAs( dataBlob, 'mstVariant.png' ); // FileSaver.js function
	  }

    var svgString = getSVGString(_this.svg.node());
	  svgString2Image( svgString, _this.width*1.5, _this.height*1.5, 'png', save ); // passes Blob and filesize String to the callback
  }

  /*** Fuction To convert and save css value eto png ***/
  // From http://bl.ocks.org/Rokotyan/0556f8facbaf344507cdc45dc3622177
  // Below are the functions that handle actual exporting:
  // getSVGString ( svgNode ) and svgString2Image( svgString, width, height, format, callback )
  function getSVGString( svgNode ) {
  	svgNode.setAttribute('xlink', 'http://www.w3.org/1999/xlink');
  	var cssStyleText = getCSSStyles( svgNode );
  	appendCSS( cssStyleText, svgNode );

  	var serializer = new XMLSerializer();
  	var svgString = serializer.serializeToString(svgNode);
  	svgString = svgString.replace(/(\w+)?:?xlink=/g, 'xmlns:xlink='); // Fix root xlink without namespace
  	svgString = svgString.replace(/NS\d+:href/g, 'xlink:href'); // Safari NS namespace fix
  
  	return svgString;

  	function getCSSStyles( parentElement ) {
  		var selectorTextArr = [];
  
  		// Add Parent element Id and Classes to the list
  		selectorTextArr.push( '#'+parentElement.id );
  		for (var c = 0; c < parentElement.classList.length; c++)
  				if ( !contains('.'+parentElement.classList[c], selectorTextArr) )
  					selectorTextArr.push( '.'+parentElement.classList[c] );
  
  		// Add Children element Ids and Classes to the list
  		var nodes = parentElement.getElementsByTagName("*");
  		for (var i = 0; i < nodes.length; i++) {
  			var id = nodes[i].id;
  			if ( !contains('#'+id, selectorTextArr) )
  				selectorTextArr.push( '#'+id );
  
  			var classes = nodes[i].classList;
  			for (var c = 0; c < classes.length; c++)
  				if ( !contains('.'+classes[c], selectorTextArr) )
  					selectorTextArr.push( '.'+classes[c] );
  		}
  
  		// Extract CSS Rules
  		var extractedCSSText = "";
  		for (var i = 0; i < document.styleSheets.length; i++) {
  			var s = document.styleSheets[i];
  			
  			try {
  			    if(!s.cssRules) continue;
  			} catch( e ) {
  		    		if(e.name !== 'SecurityError') throw e; // for Firefox
  		    		continue;
  		    	}
  
  			var cssRules = s.cssRules;
  			for (var r = 0; r < cssRules.length; r++) {
  				if ( contains( cssRules[r].selectorText, selectorTextArr ) )
  					extractedCSSText += cssRules[r].cssText;
  			}
  		}
  		
  
  		return extractedCSSText;
  
  		function contains(str,arr) {
  			return arr.indexOf( str ) === -1 ? false : true;
  		}
  
  	}
  
  	function appendCSS( cssText, element ) {
  		var styleElement = document.createElement("style");
  		styleElement.setAttribute("type","text/css"); 
  		styleElement.innerHTML = cssText;
  		var refNode = (element.hasChildNodes() && element.children[0] !== undefined) ? element.children[0] : null;
  		element.insertBefore( styleElement, refNode );
  	}
  	
  }

  function svgString2Image( svgString, width, height, format, callback ) {
  	var format = format ? format : 'png';
  
  	var imgsrc = 'data:image/svg+xml;base64,'+ btoa( unescape( encodeURIComponent( svgString ) ) ); // Convert SVG string to data URL
  
  	var canvas = document.createElement("canvas");
  	var context = canvas.getContext("2d");
  
  	canvas.width = width;
  	canvas.height = height;
  
  	var image = new Image();
  	image.onload = function() {
  		context.clearRect ( 0, 0, width, height );
  		context.drawImage(image, 0, 0, width, height);
  
  		canvas.toBlob( function(blob) {
  			if ( callback ) callback( blob);
  		});
  
  		
  	};
  
  	image.src = imgsrc;
  }
  
  /*** end unction to convert svg to png ***/

  // enable force simmulator
  play = function() {
    _this.simulation.alpha(0.2).restart();
    _this.simulationIsActive = true;
    _this.control.select(".fg-pause-control").style("display","block");
    _this.control.select(".fg-play-control").style("display", "none");
  }
  
  // disable force simulator
  stop = function() {
    _this.simulation.stop();
    _this.simulationIsActive = false;
    _this.control.select(".fg-play-control").style("display","block");
    _this.control.select(".fg-pause-control").style("display","none");
  }

  // set event on control
  _this.control.select('.fg-save-control').on('click', function() { saveSVG(); });
  _this.control.select('.fg-pause-control').on('click', function() { stop(); });
  _this.control.select('.fg-play-control').on('click', function() { play();});
  
}

// draw everythings
ForceGraph.prototype.draw = function() {

  var _this = this;
  
  // add a div for tooltip	
  _this.tooltip = this.el.append("div")
            .classed("fg-tooltip",true)
            .style("opacity",0);
            
  if( ! _this.control ) _this.createControl("forcegraph-control");
  
  // zoom transform main container
  _this.zoom = d3.zoom()
                .on("zoom", function() {
                      _this.container.attr("transform",d3.event.transform);
                      _this.zoomTranslate[0] = parseFloat(d3.event.transform.x);
                      _this.zoomTranslate[1] = parseFloat(d3.event.transform.y);
                      _this.zoomScale = parseFloat(d3.event.transform.k);
                });

  // svg container with zoom function
  _this.svg = _this.el.append("svg")
                    .attr("id", "fgid")
                    .style("width",_this.width + "px")
                    .style("height",_this.height + "px")
                    .style("pointer-events", "all")
                    .call(_this.zoom);
  // main container                
  _this.container = this.svg.append("g")
                          .style("width","100%")
                          .style("height","100%")
                          .attr("transform","translate(0,0)")
                          .style("pointer-events", "all")
                          .style("cursor","move");
  
  if( _this.root === null || _this.root === undefined) return(0);
  // draw links (line)                  
  _this.link = this.container.append("g")
              .attr("class", "fg-links")
              .selectAll("line")
              .data(_this.graph.links)
              .enter()
              .append("line")
              .attr("stroke", "green")
              .attr("stroke-width", 3)
              .attr("fill", "green");
  
  // draw nodes (circles)                  
  _this.node = this.container.append("g")
                .attr("class", "fg-nodes")
                .selectAll("circle")
                .data(_this.graph.nodes)
                .enter()
                .append("circle")
                .attr("r", function(n) { return Math.max(n.weight*250, 5)+10})//Math.max(_this.radius*parseFloat(n.weight),_this.minRadius)})
                .attr("fill", "green")
                .attr("stroke","darkgreen")
                .attr("stroke-width", 3)
                // add mouse over listener for tool-tip
                .on("mouseover", function(n) {
                   _this.tooltip.html(n.id + "<br/>prop.= " + n.weight + "<br/>count= " +n.count)
                      // OK
                      //.style("left", ( (d3.mouse(this)[0] + ( _this.zoomTranslate[0] * _this.zoomScale)) + "px"))
                      //.style("top", ( (d3.mouse(this)[1] + ( _this.zoomTranslate[1] * _this.zoomScale)) + "px"))
                      // OK but with jquery
                      //.style("left", d3.event.pageX - ($('#mstvariant svg').offset().left) + "px")
                      //.style("top", (d3.event.pageY - ($('#mstvariant svg').offset().top)) + "px")
                      .style("left",10 + "px")
                      .style("top", 30 + "px")
                      .transition()
                      .duration(200)
                      .style("opacity",0.9);
                })
                // add mouse out listener for tool-tip
                .on("mouseout", function() {
                  _this.tooltip.transition()
                  .duration(500)
                  .style("opacity",0);
                })
                // double click and zoom listener
                // to center the selected node
                .on("dblclick.zoom", function(d) {
                  d3.event.stopPropagation();
                  //_this.container.attr("transform", "translate("+ (_this.width/2 - d.x) + "," + (_this.height/2 - d.y) + ")");
	                //_this.zoom.translateTo(_this.container, d.x, d.y);
	                _this.zoom.translateTo(_this.svg, d.x, d.y);
	              })
	              .call(d3.drag()
	                        .on("start",function(e){e.x=d3.event.x; e.y=d3.event.y; ticked();})
	                        .on("drag",function(e){e.x=d3.event.x; e.y=d3.event.y; ticked();})
	                        .on("end",function(e){e.x=d3.event.x; e.y=d3.event.y; ticked();
	                                              if( _this.simulationIsActive ) _this.simulation.alpha(0.1).restart();})
	              );
	
	if(_this.havePie) {
	   
	 // arc generator using pie info a and nodes info d           
	  var arc = d3.arc()
	          .startAngle(function(d,a){ return a.startAngle})
	          .endAngle(function(d,a){ return a.endAngle;})
	          .innerRadius(0)
	          .outerRadius(function(d,a){return Math.max(d.weight*250,5)+10}) //Math.max(_this.radius*d.weight,_this.minRadius)});
	
	  var colors_pie = d3.scaleOrdinal(d3.schemeCategory10)
	                  .domain([_this.minTime, _this.maxTime]);
	                  //d3.scaleLinear()
	                  //.domain([_this.minTime, (_this.maxTime-_this.minTime)/2 , _this.maxTime])
                    //.range(["red","green","purple"]);
                
                  
  
  // for all node
  _this.nodePie = this.container.append("g")
                    .attr("class", "fg-nodes-pies")
                    .selectAll("g")
                    .data(_this.graph.nodes)
                    .enter()
                    .filter(function(d) { return (typeof d.value !== 'undefined' && d.value.length > 0);})
                    .append("g")
                    .call(d3.drag()
	                        .on("start",function(e){e.x=d3.event.x; e.y=d3.event.y; ticked();})
	                        .on("drag",function(e){e.x=d3.event.x; e.y=d3.event.y; ticked();})
	                        .on("end",function(e){e.x=d3.event.x; e.y=d3.event.y; ticked();
	                                              if( _this.simulationIsActive ) _this.simulation.alpha(0.1).restart();})
	                  );
  
  // for each node
  // and each pie values
  // draw path of a pie
  _this.nodePie
              .each( function(d,i){
                d3.select(this)
                  .selectAll("path")
                  .data(d.value)
                  .enter()
                  .append('path')
                  .attr('d', function(a,j){ return arc(d,a); })
                  .attr("fill", function(a,j){return colors_pie(d.time[a.index]);})
                  .attr("stroke", function(a,j){return colors_pie(d.time[a.index]);})
                  .on("mouseover", function(n) {
                    if( _this.timestamp ) {
                    // pie index from value and time array from graph.nodes
                    _this.tooltip.html(d.id + "<br/>time : " + _this.timeFormat(new Date(d.time[n.index])) + "<br/>prop.= " + n.value + "<br/>count= " +Math.round(n.value*d.count));
                    }
                    else _this.tooltip.html(d.id + "<br/>time : " + d.time[n.index] + "<br/>prop.= " + n.value + "<br/>count= " +Math.round(n.value*d.count));
                   _this.tooltip
                      // OK
                      //.style("left", ( (d3.mouse(this)[0] + ( _this.zoomTranslate[0] * _this.zoomScale)) + "px"))
                      //.style("top", ( (d3.mouse(this)[1] + ( _this.zoomTranslate[1] * _this.zoomScale)) + "px"))
                      // OK but with jquery
                      //.style("left", d3.event.pageX - ($('#mstvariant svg').offset().left) + "px")
                      //.style("top", (d3.event.pageY - ($('#mstvariant svg').offset().top)) + "px")
                      .style("left",10 + "px")
                      .style("top", 30 + "px")
                      .transition()
                      .duration(200)
                      .style("opacity",0.9);
                })
                // add mouse out listener for tool-tip
                .on("mouseout", function() {
                  _this.tooltip.transition()
                  .duration(500)
                  .style("opacity",0);
                })
              });
	}
	              
  // the tiny node over link                
  _this.linkNode = this.container.append("g")
                    .attr("class", "fg-nodelink")
                    .selectAll("circle")
                    .data(_this.graph.links)
                    .enter()
                    // for each links
                    .each( function(e,i) {
                      // for link weight from 1 to weight-1
                      // add a new circle
                      for(var iw=1; iw < e.weight; iw++) {
                        d3.selectAll(".fg-nodelink").append("circle")
                        .attr("id", "node-link-"+i+"-"+iw)
                        .attr("cx", 0)
                        .attr("cy",0)
                        .attr("r",3)
                        .attr("fill", "black")
                        .attr("stroke","black")
                        .attr("stroke-width", 1);
                      }
                    });
  
  // force simulator
  _this.simulation = d3.forceSimulation()
    .nodes(_this.graph.nodes)
    .on("tick", ticked)
    .force("link", d3.forceLink(_this.graph.links).id(function(d){ return d.id;}).distance( function(l){return l.weight*250}))
    .force("charge", d3.forceManyBody().strength(function(n) {return -(1/n.count)*_this.radius ;}))
    /*.force("link", d3.forceLink(_this.graph.links).id(function(d){ return d.id;}).distance( function(l){return _this.radius/2*l.weight}))
    .force("charge", d3.forceManyBody().strength(function(n) { return  - _this.radius/n.weight ;}).distanceMax(_this.radius*3))*/
    .force("center", d3.forceCenter(_this.width / 2, _this.height / 2));
    /*.on("end", function(){
      //simulation.force("link", d3.forceLink().distance( function() { return _this.radius + _this.move; }))
      
      //simulator.force("charge", d3.forceManyBody().strength(function(n) { return  - _this.radius/n.weight ;}))
      _this.simulation.alpha(0.1).restart();
      //simluator.restart();
    })
    .alphaMin(0.05);*/
    
    if( ! _this.simulationIsActive ) _this.simulation.stop();
  
  // how simulator change attributes
  function ticked() {
    _this.link
        .attr("x1", function(d) { return d.source.x; })
        .attr("y1", function(d) { return d.source.y; })
        .attr("x2", function(d) { return d.target.x; })
        .attr("y2", function(d) { return d.target.y; });
    
    // for the tiny node (between each real nodes)
    _this.link.each( function(d) {
      
      if(d.weight > 1) {

        // get nodes radius
        var radiussource = parseFloat(_this.graph.nodes.filter( function(e){ return e.id === d.source.id;})[0].weight) * 250 + 10; //_this.radius;
        var radiustarget = parseFloat(_this.graph.nodes.filter( function(e){ return e.id === d.target.id;})[0].weight) * 250 + 10; //_this.radius;
        
        // hypotenus 
        var linelength = Math.sqrt(Math.pow(d.source.x - d.target.x,2) + Math.pow(d.source.y - d.target.y,2) );
        // get cos and sin angle
        var cosangle = (d.target.x - d.source.x) / linelength;
        var sinangle = (d.target.y - d.source.y) / linelength;
        // step length for each new node
        var steplength = (linelength - radiussource - radiustarget) / d.weight;
        
        // for 1 to weight-1 
        for(var i=1; i<d.weight; i++) {
          _this.linkNode.select("#node-link-"+d.index+"-"+i)
            .attr("transform", "translate("+ (d.source.x + (cosangle * (radiussource+(steplength*i)))) + "," + (d.source.y + (sinangle * (radiussource+(steplength*i))))  + ")");
        }
      }
      
    });

    _this.node
        .attr("cx", function(d) { return d.x; })
        .attr("cy", function(d) { return d.y; });
    
     // update pie position    
     if(_this.havePie) _this.nodePie.attr("transform", function(d){ return "translate(" +d.x + ","+d.y+")" });
  }
  
  // node stroke animation
  _this.signe = 1;
  d3.interval(function(elapsed) {
    _this.signe = -1*_this.signe
    _this.node.transition()
        .duration(1200)
        .attr("stroke-width", 3+_this.signe*2);
  },1000);
}

// forcegraph redraw everythings
ForceGraph.prototype.redraw = function(width, height) {
  var _this = this;
  this.hide();
  
  if( width !== undefined ) this.width = width;
  if( height !== undefined ) this.height = height;
  
  if(this.width < this.height) {
    this.radius = this.width/4;
  }
  else {
    this.radius = this.height/4;
  }
  
  if(_this.control) _this.control.remove();
  _this.control = null;
  if(_this.nodePie) _this.nodePie.remove();
  _this.node.remove();
  _this.linkNode.remove();
  _this.link.remove();
  _this.linkNode.remove();
  _this.container.remove();
  _this.svg.remove();
  _this.tooltip.remove();

  _this.draw();
  _this.show();
}

// show
ForceGraph.prototype.show = function(){
  this.el.style("display",'block');
}

// hide
ForceGraph.prototype.hide = function(){
  this.el.style("display",'none');
}

// addPie
// Active and set PieCharts into graph nodes
// need data prop_node
// An array indexed by node id (rows)
// that contain time and values in columns
ForceGraph.prototype.addPie = function(){
  
  _this = this;
  
  if( _this.root.node_prop === undefined || _this.root.node_prop === null || _this.root.node_prop.length === undefined) return false;
  
  //console.log(_this.root.node_prop);
  
  this.havePie = true;
  
  // generate pie info
  var pie = d3.pie()
	            .value(function(d) {return d})
	            .sort(null);
  
  
  this.minTime = Infinity;
  this.maxTime = 0;
  this.root.node_prop.forEach(function(el,i) {
    
    if(Object.assign !== undefined) Object.assign(_this.graph.nodes[i], {"value" : pie(el.value), "time" : el.time });
    else {
      _this.graph.nodes[i]["value"] = pie(el.value);
      _this.graph.nodes[i]["time"] = el.time;
    }
    
    el.time.forEach(function(el,i){
      if(_this.minTime > el) _this.minTime = el;
      if(_this.maxTime < el) _this.maxTime = el;
    });
  })
  //console.log(_this.minTime + " " + _this.maxTime);
  //console.log(_this.graph.nodes[0]);
}

// makeGraph
// fill graph with data
// graph.nodes = [{"id":id, "weight": weight},...]
// graph.edges = [{"source":sourceID, "target": targetId, "weight": weight},...]
var makeGraph = function(graph, data) {

  function makeNode(id, weight, count) {
    return {"id": id.replace(/\./g,'-'), "weight": parseFloat(weight), "count": parseInt(count)};
  }
  
  function makeLink(sourceId, targetId, weight) {
    return {"source": sourceId.replace(/\./g,'-'), "target": targetId.replace(/\./g,'-'), "weight": parseInt(weight)};
  }

  if(data.nodes instanceof Array)
  data.nodes.forEach( function(n){
    graph.nodes.push(makeNode(String(n.ID), n.weight, n.count));
  });
  
  if(data.edges instanceof Array)
  data.edges.forEach( function(e){
    graph.links.push(makeLink(String(e.source), String(e.target), e.weight));
  });
}


