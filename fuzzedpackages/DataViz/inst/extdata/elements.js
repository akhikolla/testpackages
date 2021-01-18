/**
 * Getting the current size of the window in order to scale everything
 */
let margin = { top: 0, right: 0, bottom: 0, left: 0 },
    width = $(window).width(),
    height = $(window).height(),
    halfHeight = height / 2;

let tooltip = d3.select('#tooltip');

/**
 * Setting the windows through d3
 */
let svg = d3.select('.container')
  .append('svg')
  .attr('width', width)
  .attr('height', height );

//$(".container").css({ left: width / 4 , right: 3*width/4});
