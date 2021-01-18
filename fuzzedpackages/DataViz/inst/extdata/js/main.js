var slider = document.getElementById("myRange");
USER_SPEED_SLIDE = 5000/slider.value; // Display the default slider value

// Update the current slider value (each time you drag the slider handle)
slider.oninput = function() {
  USER_SPEED_SLIDE = 5000/this.value;
}

Object.size = function(obj) {
    var size = 0, key;
    for (key in obj) {
        if (obj.hasOwnProperty(key)) size++;
    }
    return size;
};

// Get the size of an object
var sizeSchedule = Object.size(data[0]);
var sizeElements = Object.size(data);
var TRUE_MINUTE = 0;


var arrayPreDistinct= [];

for (var i = 0; i < sizeElements; i++) {
  for (var j = 1; j < sizeSchedule+1; j++) {
    arrayPreDistinct.push(data[i][j])
  }
}

var unique_ids = [...new Set(arrayPreDistinct)];
var sizeUniqueItems = Object.size(unique_ids);

var POPULATIONMAP = new Map();
var shapes = [];

for (var i = 0; i < sizeUniqueItems ; i++) {
  POPULATIONMAP.set(unique_ids[i],{"popu" : 0, "coordX" :0, "coordY": 0, "index" :i, "shape" : d3.symbolSquare});
  // unique_ids[i].push(d3.symbolSquare);
  // shapes.push(d3.symbolSquare)
  shapes.push({"label" : unique_ids[i], "shapey" : d3.symbolSquare})
}

//*
//* Actuellement l'application refait le calcul pour le cercle pour chaque boule,
//* vaut mieux fier l'élément "ancre" en dur et donc juste vérifier dans le quel
//* il doit aller, ça sera plus facile pour permettre à l'utilisateur de bouger les points
//*

var pause = true;
var tmpPop ={};


var width = $(window).width(),
        height = $(window).height(),
        widthmiddle = width/2,
        heightmiddle = height/2;

var curr_minute = 0;

// Activity to put in center of circle arrangement
var center_act = "NA",
        center_pt = { "x": widthmiddle, "y": heightmiddle };

// Coordinates for activities

POPULATIONMAP.forEach(function(d, i) {
    var theta = 2 * Math.PI / (sizeUniqueItems);
    tmpPop = d;
    tmpPop.coordX = 250 * Math.cos(d.index * theta)+widthmiddle;
    tmpPop.coordY = 250 * Math.sin(d.index * theta)+heightmiddle;
    POPULATIONMAP.set(i, tmpPop);
});

// Start the SVG
var svg = d3.select("#chart").append("svg")
        .attr("width", width)
        .attr("height", height);

var svg2 = d3.select("#bar-chart").append("svg")
        .attr("width", width)
        .attr("height", height);

// A node for each person's schedule
var nodes = data.map(function (o, i) {
    var act2 = o[TRUE_MINUTE+1];
    tmpPop = POPULATIONMAP.get(act2);
    tmpPop.popu += 1;
    POPULATIONMAP.set(act2,tmpPop);
    return {
        activity : act2,
        radius: 3,
        x: POPULATIONMAP.get(act2).coordX,
        y: POPULATIONMAP.get(act2).coordY,
        color: "#00cdc0",
        schedule: o
    }
});


var FORCE = d3.forceSimulation(nodes)
      .force('x', d3.forceX().x(function(d) {
          return d.x;
        }))
        .force('y', d3.forceY().y(function(d) {
            return d.y;
          }))
          .force('collision', d3.forceCollide(4))

        .on("tick", tick);

var circle = svg.selectAll("circle")
        .data(nodes)
        .enter().append("circle")
        .attr("r", function (d) {
            return d.radius;
        })
        .attr("class","dot");


//--------------------------------------------------
        // Activity labels
        var label = svg.selectAll("text")
                .data(unique_ids)
                .enter().append("text")
                .attr("class", "actlabel")
                .attr("draggable","true")
                .attr("x", function (d, i){
                  return POPULATIONMAP.get(d).coordX;

                })
                .attr("y", function (d, i) {
                  return POPULATIONMAP.get(d).coordY;
                });

        label.append("tspan")
                // .attr("x", function () {
                //     return d3.select(this.parentNode).attr("x");
                // })
                .attr("text-anchor", "middle")
                .text(function (d) {
                    return d;
                });
        label.append("tspan")
                .attr("dy", "14")
                .attr("x", function () {
                    return d3.select(this.parentNode).attr("x");
                })
                .attr("text-anchor", "middle")
                .attr("class", "actpct")
                .text(function (d,i) {
                    return POPULATIONMAP.get(d).popu;
                });








                var svg_dx = 800,
                    svg_dy = 400,
                    margin_x = 100;

                // var shapes = [ d3.symbolSquare];

                svg.append("svg")
                            .attr("width", svg_dx)
                            .attr("height", svg_dy);


                var symbol = d3.symbol().size([3000])
                var drag_behavior = d3.drag()
                                      .on("start", dragstarted)
                                      .on("drag", dragged);
                  var squares = svg.append("g")
                   .selectAll("path")
                   .data(shapes);

                   var squaresEnter = squares.enter()
                   .append("path")
                   // .append("text")
                   .attr("d", symbol.type(shape => shape.shapey))
                   .attr("id", (shape, i) =>  shapes[i].label)
                   .attr("transform", (shape, i) => "translate(" + POPULATIONMAP.get(shapes[i].label).coordX + ", " + POPULATIONMAP.get(shapes[i].label).coordY + ")")
                   .style("fill", (shape, i) => '#red')
                   .style("opacity", 0.5)
                   .call(drag_behavior)
                   .transition()
                   .duration((shape, i) => i * 800)



                   squares.append("text")
                   .attr("dy", "-35")
                   // .attr("opacity", "1")
                   .style("fill","red")
                   .text("text");
                   // .attr("transform", (shape, i) => "translate(" + x(i) + "," + (svg_dy / 2) + ")");

                function dragstarted() {
                  d3.select(this).raise();
                }

                function dragged(shape) {

                    var dx = d3.event.sourceEvent.offsetX,
                        dy = d3.event.sourceEvent.offsetY;
                    tmpPop = POPULATIONMAP.get(this.id);
                    tmpPop.coordX =dx;
                    tmpPop.coordY = dy;
                    POPULATIONMAP.set(i, tmpPop);

                    d3.select(this)
                      .attr("transform", shape => "translate(" + dx + "," + dy + ")");
                }










//--------------------------------------------------

// Update nodes based on activity and duration
function timer() {
    d3.range(nodes.length).map(function (i) {
        var curr_node = nodes[i];
                    tmpPop = POPULATIONMAP.get(curr_node.activity);
                    tmpPop.popu -= 1;
                    POPULATIONMAP.set(curr_node.activity,tmpPop);
                    curr_node.activity =  curr_node.schedule[TRUE_MINUTE+1];
                    tmpPop = POPULATIONMAP.get(curr_node.activity);
                    tmpPop.popu += 1;
                    curr_node.x = POPULATIONMAP.get(curr_node.activity).coordX;
                    curr_node.y = POPULATIONMAP.get(curr_node.activity).coordY;
                    POPULATIONMAP.set(curr_node.activity,tmpPop);
                    FORCE
                    .force('x', d3.forceX().x(function(d) {
                        return d.x;
                      }))
                      .force('y', d3.forceY().y(function(d) {
                          return d.y;
                        }))
                        FORCE.alpha(1).restart();
    });

    curr_minute += 1;

    // Update percentages
    label.selectAll("tspan.actpct")
            .text(function (d,i) {
              return POPULATIONMAP.get(d).popu;
            });

    // Update times
    TRUE_MINUTE = curr_minute % sizeSchedule;
    var showMinute = TRUE_MINUTE == 0 ? sizeSchedule : TRUE_MINUTE;
    d3.select("#current_time").text(units[showMinute]);

    if (!pause) {
        setTimeout(timer, USER_SPEED_SLIDE);
    }
}

// Time Control
d3.select("#play").style("display", "initial").on("click", function () {
    setTimeout(timer,USER_SPEED_SLIDE);
    pause = false;
    d3.select("#play").style("display", "none");
    d3.select("#pause").style("display", "initial");
});
d3.select("#pause").on("click", function () {
    pause = true;
    d3.select("#play").style("display", "initial");
    d3.select("#pause").style("display", "none");
});

function tick(e) {

var u = d3.select('svg')
   .selectAll('circle')
   .data(nodes)

 u.enter()
   .append('circle')
   .attr('r', function(d) {
     return d.radius
   })
   .merge(u)
   .transition().duration(USER_SPEED_SLIDE*10/sizeSchedule).ease(d3.easeLinear)
   .attr('cx', function(d) {
     return d.x
   })
   .attr('cy', function(d) {
     return d.y
   })

 u.exit().remove()

}
