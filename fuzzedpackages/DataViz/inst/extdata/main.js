  /**
   * Using D3's curveCardinal for curving the pyramid, all three points are part of the curve
   */
  let lineGenerator = d3.line()
    .curve(d3.curveCardinal);


  let g = svg.selectAll('.line')
    .data(data)
    .enter()

  /**
   * activePath is null, it is the token of the activeCurve to color it or no
   */
  let activePath = null;

  /**
   * Drawing the lines
   */
  g.append('path')
    .attr('d', (d) => lineGenerator(d.points))
    .attr('class', 'line')
    .attr('stroke', (d) => d.col4)
    .attr('stroke-width',(d) => d.col5)

    .on('mouseover', function (d)
    {
      activePath = d;
      /**
       * If the mouse is on the curve, highlight it in red
       */
      let color = 'red';

      $(this).attr('stroke', color);
      $(this).attr('stroke-width', 5.5);

      /**
       * TEST
       */
      tooltip
        .select('#col3')
        .text(d.col3);
      /**
       * TEST
       */
       tooltip
         .select('#col4')
         .text(names1);
       tooltip
         .select('#col5')
         .text(names2);

       /**
        * Tooltips for the values of the active curves.
        */
      tooltip
        .select('#col1')
        .text(d.col1);
      tooltip
        .select('#col2')
        .text(d.col2);

      tooltip
      .select('#height')
      .text(d.col2-d.col1);

    })

    /**
     * If not selecting the curve anymore, put back normal colour
     */
    .on('mouseout', function (d)
    {
        if (activePath === d)
        {
        d3.selectAll('.line').attr('stroke', (d) =>d.col4)
        .attr('stroke-width',(d) => d.col5)
      }
    });

  var rectLeft = svg.append("rect")
      .attr("x", 0)
      .attr("y", 0)
      .attr("height", height)
      .attr("width", width / 4)
      .style("stroke", 'none')
      .style("fill", "white");
    var rectRight = svg.append("rect")
      .attr("x", 3*width/4)
      .attr("y", 0)
      .attr("height", height)
      .attr("width", width / 4)
      .style("stroke", 'none')
      .style("fill", "white");
    if (offSet < 0 || (yMax <= 0 && yMin < 0))
    {
        var rectTop = svg.append("rect")
          .attr("x",0)
          .attr("y", 0)
          .attr("height", halfHeight)
          .attr("width", width)
          .style("stroke", 'none') 
          .style("fill", "white");
    }
    if (offSet > 0 || (yMin >=0 && yMax > 0))
    {
        var rectBot = svg.append("rect")
          .attr("x", 0)
          .attr("y", halfHeight)
          .attr("height", halfHeight)
          .attr("width", width)
          .style("stroke", 'none')
          .style("fill", "white");
    }

    svg.append('g')
      .attr('transform', 'translate(' + 0 + ',' + halfHeight + ')')
      .call(d3.axisBottom(x)
        .tickSizeOuter(0))
$(window).resize(function(){location.reload();});
