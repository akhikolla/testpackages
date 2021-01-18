/* This file is part of SMITIDvisu package.
* Copyright (C) 2018-2019 Jean-François Rey <jean-francois.rey@inra.fr>
*                         Julien Boge <julien.boge.u@gmail.com>
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

/*****************************************************************************
 * ▼ Map Transmission Tree
 * ---------------------------------------------------------------------------
 * Perform animations using D3.js and Leaflet.TimeDimension libraries on a
 * Leaflet map.
 * ---------------------------------------------------------------------------
 * Notes:
 *  - GeoJSON coordinates and Leaflet latLng are inversed.
 *  - `this._currentTimeIndex` corresponds to the running animation.
 *    It may differ from `this._timeDimension.getCurrentTimeIndex()`.
 *  - Be carreful with the `this` variable inside d3 functions.
 *    Inside arrow functions, it refers to this class context (MapTT).
 *    Otherwise, it refers to the svg element handled by d3.
 *    If we need both contexts, we usually declare a variable `that` which
 *    refers to this class context (to summarize: `that` = this class,
 *    `this` = the d3 element).
 *****************************************************************************/
var log = {
	fired: true,
	update: true,
	feature: true,
	anim: true,
}

L.MapTT = L.TimeDimension.Layer.GeoJson.extend({

	_defaultOptions: {
		circleRadius: 6,
		nodeColorByState: {},
		defaultNodeColor: 'steelblue',
		moveEdgeColor: 'royalblue',
		color1: 'green',
		color2: 'red',
		nbColors: 10,
		minWeight: 0,
		maxWeight: 1,
		weight1: 0,
		weight2: 0,
		autoFocus: true,
		keepOldFeatures: true,
		onSelectHost: null,
		player: null,
	},

	//---------------------------------------------------------------------------
	// * Basic functions
	//---------------------------------------------------------------------------
	// Initialize a new MapTT instance
	initialize: function(layer, options) {
		console.log("initialize");
		L.TimeDimension.Layer.GeoJson.prototype.initialize.call(this, layer, options);

		this._setOptions(options);

		this._oldFeatures = [];
		this._currentFeatures = [];
		this._historic = null;
		this._skipAnimation = false;

		this._highlightedHost = null;
	},
	// Set the options. It fixs the default `L.setOptions` to ensure that all the
	// options have the necessary values to not break anything.
	_setOptions: function(options) {
		// Set the player
		this._player = options.player || null;
		if (this._player == null) {
			throw('Missing TimeDimension Player. You must provide it through the options.');
		}

		// Reset null options with defaut non-null values
		const nonNullDefaultOptions =
			Object.keys(this._defaultOptions)
			.filter(opt => this._defaultOptions[opt] != null);

		for (const opt of nonNullDefaultOptions) {
			if (this.options[opt] == null) {
				this.options[opt] = this._defaultOptions[opt];
			}
		}

		// Set the fraction duration of each animation step
		this._setFractionsDurations(options.durations ? options.durations : {});
	},
	// Set each step duration percentage
	_setFractionsDurations: function(options) {
		let fractions = {};

		fractions.fractionNodes = options.nodes || (2/4);
		fractions.fractionInfections = options.infections || (1/4);
		fractions.fractionPause = options.pause || (1/4);

		this._durations = fractions;

		this._setStepsDurations();
	},
	// Set each step duration
	_setStepsDurations: function() {
		const totalDuration = this._player.getTransitionTime();
		console.log("totalDuration:", totalDuration);

		this._durations.nodes      = parseInt(totalDuration * this._durations.fractionNodes);
		this._durations.infections = parseInt(totalDuration * this._durations.fractionInfections);
		this._durations.pause      = parseInt(totalDuration * this._durations.fractionPause);

		console.log("this._durations:", this._durations);
	},
	// Initialize SVG on the map
	_initializeSVG: function() {
		let svg = L.svg(); // Create a Leaflet SVG layer to display our features
		this._map.addLayer(svg);
		this._currentLayer = svg;

		// Allows pointer interactions (mouseover, click...) with the layer's elements
		this._currentLayer._container.setAttribute("pointer-events", "visible");

		// Bind a Leaflet Tooltip to this layer
		this.bindTooltip(L.tooltip());

		// Select the SVG group we will use to display d3js animations
		this._svg = d3.select(svg._rootGroup); // We can not work directly on the Leaflet SVG layer, so we must select the _rootGroup attribute
		this._svg.classed("leaflet-zoom-hide", true); // Hide SVG elements when zooming

		// Initialize the animation variables
		this._resetAnimation();
		this._oldTimeIndex = this._timeDimension.getCurrentTimeIndex();

		// Add a <defs> section to the svg
		this._defs = this._svg.append('svg:defs');

		this.setEdgeColorScheme();
	},
	// Define a scheme of colors for the edges
	setEdgeColorScheme: function() {
		const opt = this.options;

		const interpolate = d3.interpolateRgb(opt.color1, opt.color2),
		      cs = d3.quantize(interpolate, opt.nbColors)
		             .map(c => d3.color(c).hex());

		const domain1 = [opt.minWeight, opt.weight1],
		      domain2 = [opt.weight1, opt.weight2],
		      domain3 = [opt.weight2, opt.maxWeight];


		let r1 = cs.slice(0, opt.minColor),
		    r2 = cs.slice(opt.minColor, opt.maxColor),
		    r3 = cs.slice(opt.maxColor, opt.nbColors);

		if (r1.length === 0) r1 = [cs[0]];
		if (r2.length === 0) r2 = [cs[opt.minColor-1]];
		if (r3.length === 0) r3 = [cs[opt.nbColors-1]];

		const ecs1 = d3.scaleQuantize().domain(domain1).range(r1),
		      ecs2 = d3.scaleQuantize().domain(domain2).range(r2),
		      ecs3 = d3.scaleQuantize().domain(domain3).range(r3);

		// Add arrows markers for each color in the scheme
		this._addArrowMarkers(cs);
		this._addArrowMarkers(['#1e90ff']);

		// Return the color corresponding to the infection probability
		this._edgesColorScheme = function(n) {
			if (n < domain1[1]) return ecs1(n);
			if (n >= domain2[0] && n < domain2[1]) return ecs2(n);
			return ecs3(n);
		}

		// Re-draw infection edges if they are not being drawn
		if (!this._isAnimationRunning()) {
			this._updateInfectionsEdgesColor();
		}
	},
	// Called when the instance is attached to the Leaflet map
	onAdd: function(map) {
		L.TimeDimension.Layer.GeoJson.prototype.onAdd.call(this,map);
		this._map.on("zoomend", this._onZoom, this);
		this._timeDimension.on("loop", this._onLoop, this);
		this._timeDimension.on("timedivert", this._onTimeDivert, this);
		this._player.on("speedchange", this._onSpeedChange, this);
		this._focusOnFeatures();
		if (!this._player.isPlaying()) {
			this._player.start();
		}
	},
	// Called when the instance is detached from the Leaflet map
	onRemove: function(map) {
		this._map.removeLayer(this._currentLayer);
		this._map.off("zoomend", this._onZoom, this);
		this._map.off("click", this._reinitializeHighlightedHost, this);
		this._timeDimension.off("loop", this._onLoop, this);
		this._timeDimension.off("timedivert", this._onTimeDivert, this);
		this._player.off("speedchange", this._onSpeedChange, this);
		L.TimeDimension.Layer.GeoJson.prototype.onRemove.call(this,map);
	},
	// Add an SVG arrow marker for any new color in the given array
	_addArrowMarkers: function(colors) {
		const markers = this._defs.selectAll('marker');
		const currentColors = markers.data();
		colors = currentColors.concat(colors);
		colors = [...new Set(colors)]; // keep distinct colors

		markers
			.data(colors)
			.enter()
				.append('svg:marker')
				.attr('id', d => 'marker_' + d.substring(1))
				.attr('markerHeight', 5)
				.attr('markerWidth', 5)
				.attr('markerUnits', 'strokeWidth')
				.attr('orient', 'auto')
				.attr('refX', 0)
				.attr('refY', 0)
				.attr('viewBox', "-5 -5 10 10")
				.append('svg:path')
					.attr('d', "M 0,0 m -5,-5 L 5,0 L -5,5 Z")
					.attr('fill', d => d);
	},
	// Call the user's function when a host is selected
	_selectHost: function(hostId) {
		if (typeof this.options.onSelectHost === "function") {
			this.options.onSelectHost(hostId);
		}
	},

	//---------------------------------------------------------------------------
	// * On Events
	//---------------------------------------------------------------------------
	// Reposition nodes if needed, then resume animation
	_onZoom: function(ev) {
		if (log.fired) console.log("-------------- Zoom -------------");
		try {
			this._repositionFeatures();
			// Resume animation after it has been broken.
			// Step 0 (updateMapView) did not break animation, so we don't run anew
			if (this._step > 0) {
				this._runAnimation();
			}
		} catch (err) {
			console.log("[zoom] Error caught!");
			console.log(err);
			console.trace();
		}
	},
	// Called when the user clicked on the timebar
	// Display immediatly the current time
	_onTimeDivert: function(data) {
		if (log.fired) console.log("----------- TimeDivert -----------");

		if (log.fired) console.log("[timedivert] newIndex:", data.newIndex);
		if (log.fired) console.log("[timedivert] this._currentTimeIndex:", this._currentTimeIndex);
		if (log.fired) console.log("[timedivert] this._oldTimeIndex:", this._oldTimeIndex);

		this._currentFeatures = []; // Use only the new features to set the map view
		this._skipAnimation = true;
	},
	// Reset variables so that the next animation is performed in the appropriate
	// order
	_onLoop: function(data) {
		if (log.fired) console.log("-------------- Loop -------------");

		this._currentFeatures = []; // Use only the new features to set the map view
		this._resetAnimation();

		// Set oldTimeIndex to perform the animation steps in the appropriate order
		if (data.lowerLimit) {
			this._oldTimeIndex = data.newIndex - 1;
		} else {
			this._oldTimeIndex = data.newIndex + 1;
		}

		if (log.fired) console.log("[loop] oldTimeIndex:", this._oldTimeIndex);
	},
	// Update animation steps duration
	_onSpeedChange: function() {
		this._setStepsDurations();
	},

	//---------------------------------------------------------------------------
	// * Animation main methods
	//---------------------------------------------------------------------------
	// Rewrite _update function
	_update: function() {
		if (log.update) console.log("============= UPDATE =============");

		if (!this._map || !this._loaded) return;

		if (!this._svg) {
			this._initializeSVG();
		}

		if (this._timeDimension.getCurrentTimeIndex() < 0) return;

		if (this._historic != null) {
			// move edges are the only elements we need to clear,
			// animation will clear the rest automatically
			this._clearMovesEdges(0);
			this._historic = null;
			this._moves = null;
		}

		console.log("[update]          oldTimeIndex:", this._oldTimeIndex);
		console.log("[update] this.currentTimeIndex:", this._currentTimeIndex);
		console.log("[update]   td.currentTimeIndex:", this._timeDimension.getCurrentTimeIndex());
		console.log("[update] animationRunning:", this._isAnimationRunning());
		if (this._isAnimationRunning()) {
			this._finishRunningAnimation();
		} else if (this._timeDimension.getCurrentTimeIndex() != this._oldTimeIndex) {
			try {
				if (!this._player.isPlaying()) {
					this._skipAnimation = true;
				}
				this._initializeAnimation();
				this._runAnimation();
			} catch(err) {
				console.log("[update] >> Animation interrupted <<");
				console.log(err);
				console.trace();
			}
		}
	},
	// Initialize the animation for the TimeDimension's current time
	_initializeAnimation: function() {
		// Determine the current time features
		const currentFeatures = this._getCurrentFeatures();

		// Set the current features
		this._oldFeatures = this._currentFeatures;
		this._currentFeatures = currentFeatures;

		// Select nodes elements
		this._nodes = this._selectNodes(currentFeatures);
		this._newNodes = null;

		// Determine infections lines coordinates from nodes' data
		this._infectionsData = this._determineInfections(currentFeatures);

		// Initialize variables for the new animation
		this._step = 0;
		this._currentTimeIndex = this._timeDimension.getCurrentTimeIndex();

		// Set the related hosts at the current time
		if (this._highlightedHost != null) {
			const host = this._currentFeatures.find(f => this._highlightedHost.id === f.properties.id);
			this._highlightedHost.relatedHosts = this._getRelatedHosts(host);
		}
	},
	// Run or resume an animation at the appropriate step.
	// If runNextAnimation is true, call the next animation after this one ended
	_runAnimation: async function(runNextAnimation = true) {
		if (log.anim) console.log("** ANIMATION **");

		if (runNextAnimation && this._isAnimationFinished()) {
			if (this._player.isPlaying()) {
				console.log("[tick] Tic & Tac!");
				this._player._tick();
			}
			return;
		}

		// Initialize a nonce that will break this animation instance
		// if a newer animation is running
		const animationNonce = this._animationNonce = new Object();

		if (log.anim) console.log("[anim] this._skipAnimation:", this._skipAnimation);
		if (log.anim) console.log("[anim] this._step:", this._step);
		if (log.anim) console.log("[anim] >>>");

		const steps = this._getSteps();

		if (runNextAnimation) {
			// Add a pause step at the end of the animation
			steps.push({"methods": [this._pause],
			            "duration": this._durations.pause,
			            "label": "pause" });
		}

		for ( ; this._step < steps.length ; this._step++) {
			await this._runAnimationStep(steps[this._step]);
			if (animationNonce != this._animationNonce) return;
		}

		console.log("[anim] steps ended successfully");

		// Reset animation variables
		this._step = null;
		this._oldTimeIndex = this._currentTimeIndex;
		this._currentTimeIndex = null;
		this._skipAnimation = false;

		if (runNextAnimation && this._player.isPlaying()) {
			console.log("[tick] Tic & Tac!");
			this._player._tick();
		}
	},
	// Finish the running animation
	_finishRunningAnimation: async function() {
		// Interrupt any transitions and perform all the steps of the running
		// animation immediatly
		try {
			if (log.anim) console.log("[finish] start");
			// const goingForward = this._oldTimeIndex <= this._currentTimeIndex;
			const goingForward = this._currentTimeIndex <= this._timeDimension.getCurrentTimeIndex();
			const playerGoingForward = this._player.getSteps() > 0;
			const goingSameDirection = goingForward === playerGoingForward;

			// Finish the running animation immediatly
			this._skipAnimation = true;
			await this._runAnimation(false);

			// Perform the next animation (which is the one TimeDimension is currently at)
			if (!this._player.isPlaying() || !goingSameDirection) {
				// If we diverge from the player « direction » (ie. player is going forward
				// while this animation is backward), then perform the animation immediatly
				this._skipAnimation = true;
			}
			this._initializeAnimation();
			this._runAnimation();
			if (log.anim) console.log("[finish] stop");
		} catch(err) {
			console.log("[finish] >> CATCH <<");
			console.log("err:", err);
			console.trace();
		}
	},
	// Reset animation : displayed features, SVG groups, variables, etc.
	_resetAnimation: function() {
		if (log.anim) console.log("[anim] reset");

		// Ensure that the current animation is broken
		this._animationNonce = new Object();

		// Reset D3JS selections
		this._nodes = null;
		this._newNodes = null;

		this.closeTooltip(); // Close the bound Leaflet Tooltip
		this._featureOnMouseover = null; // id of the node on which the pointer is

		// Remove all SVG groups
		this._svg.selectAll("g").remove();

		// Add SVG groups
		this._svg.append("g").attr("id", "points");
		this._svg.append("g").attr("id", "move");
		this._svg.append("g").attr("id", "infection");

		// Reset animation variables
		this._currentTimeIndex = null;
		this._step = null;
		this._transitionDuration = null;
		this._transitionStartTime = null;

		this._oldTimeIndex = null;
	},
	// Return a promise waiting for a certain amount of time if player.isPlaying
	_pause: function(duration) {
		return this._player.isPlaying() ?
			new Promise(resolve => setTimeout(resolve, this._durations.pause)) :
			null;
	},
	// Return the animation steps (order & duration) depending if animation
	// is played forward of backward
	_getSteps: function() {
		let steps = [];
		const goingForward = this._oldTimeIndex <= this._currentTimeIndex;

		if (this.options.autoFocus) {
			steps.push({
				"methods": [this._focusOnFeatures],
				"label": "waitmapview"
			});
		}

		if (goingForward) {
			steps.push({"methods": [this._clearInfectionsEdges],
			            "duration": this._durations.infections,
			            "label": "clearinfect" });
			steps.push({"methods": [this._addNewNodes,
			                        this._removeOldNodes,
			                        this._updateNodesPosition],
			            "duration": this._durations.nodes,
			            "label": "update"
			          });
			steps.push({"methods": [this._updateNodesColor,
			                        this._drawInfectionsEdges],
			            "duration": this._durations.infections,
			            "label": "infect" });
		} else {
			steps.push({"methods": [this._clearInfectionsEdgesReversed],
			            "duration": this._durations.infections,
			            "label": "clearinfect" });
			steps.push({"methods": [this._updateNodesColor],
			            "duration": this._durations.infections,
			            "label": "infect" });
			steps.push({"methods": [this._addNewNodes,
			                        this._removeOldNodes,
			                        this._updateNodesPosition],
			            "duration": this._durations.nodes,
			            "label": "update"
			          });
			steps.push({"methods": [this._drawInfectionsEdgesReversed],
			            "duration": this._durations.infections,
			            "label": "infect" });
		}

		return steps;
	},

	//---------------------------------------------------------------------------
	// * Animation utils methods
	//---------------------------------------------------------------------------
	// Perform an animation step, waiting for d3js transitions to finish.
	// It can also resume transitions (stop & restart at the right time) for the
	// given step
	_runAnimationStep: function(step) {
		console.log(`    [${step.label}] >>>`);
		let duration = 0,
		    elapsed = 0;

		if (!this._skipAnimation) {
			if (this._transitionDuration > 0) {
				duration = this._transitionDuration;
			} else if (step.duration > 0) {
				duration = step.duration;
			}

			// Resume animation by re-defining the elapsed time and the remaining duration
			if (this._transitionStartTime != null) {
				const elapsedTime = Date.now() - this._transitionStartTime;

				elapsed = Math.min(elapsedTime / duration, 1);
				duration = Math.max(duration - elapsedTime, 0);
			}
		}

		this._transitionDuration = duration;
		const methods = step.methods.map(m => this[m.name](duration, elapsed));

		return Promise.all(methods)
			.then(() => {
				console.log(`    [${step.label}] <<<`);
				// Reset animation step variables
				this._transitionStartTime = null;
				this._transitionDuration = null;
			});
	},
	// Check if an animation is running
	_isAnimationRunning: function() {
		return this._currentTimeIndex != null;
	},
	// Check if the animation is finished
	_isAnimationFinished: function() {
		return !this._isAnimationRunning() &&
		       this._oldTimeIndex === this._timeDimension.getCurrentTimeIndex();
	},
	// Return the ease function for a given time + the proportion of elapsed time
	// It is useful if we have to « resume » the transition
	// Read documentation :
	// https://github.com/d3/d3-transition/blob/master/README.md#transition_ease
	_easeFunction: function(t, elapsed) {
		t = t + elapsed - t * elapsed;
		return d3.easeCubic(t);
	},
	// Set the transition start time if it is not defined
	_setTransitionStartTime: function() {
		if (this._transitionStartTime == null) {
			this._transitionStartTime = Date.now();
		}
	},
	// Reposition the appropriate features if we are at the appropriate step
	_repositionFeatures: function() {
		this._repositionNodes(this._nodes.exit());

		const goingForward = this._oldTimeIndex <= this._currentTimeIndex;

		// Reposition nodes unless they are already moving
		if ((goingForward && this._step != 2) ||
		    (!goingForward && this._step != 3)) {
			this._repositionNodes(this._nodes);
		}

		// Reposition edges unless they are being drawn
		if ((goingForward && this._step != 2) ||
		    (!goingForward && this._step != 3)) {
			this._repositionEdges(this._infections);
		}

		if (this._historic != null) {
			this._repositionEdges(this._moves);
		}
	},

	//---------------------------------------------------------------------------
	// * Nodes methods
	//---------------------------------------------------------------------------
	// Add undisplayed nodes elements according to features data, and merge them
	// to `this._nodes`. Also add the mouseover event to display node's info
	_addNewNodes: function(duration, elapsed) {
		console.log("    [update][   add] start");

		if (this._nodes.enter().empty()) {
			console.log("    [update][   add] enter().empty():", true);
			return null;
		}

		if (!this._newNodes) {
			this._newNodes = this._nodes.enter()
				.append("circle")
				.attr("r", this.options.circleRadius)
				.style("fill", d => this._getNodeFillColor(d))
				.style("stroke", d => this._getNodeStrokeColor(d))
				.style("stroke-width", 2)
				.style("opacity", 0)
				.each(function(d) { this._oldCoordinates = d.geometry.coordinates; }) // here : this = node
				.on("mouseover", d => {
					// After a short time period, display node's informations
					// within the Layer tooltip
					this._mouseoverTimeout = setTimeout(() => {
						// Always place the tooltip at the node destination to prevent
						// strange positionning when the node is moving
						const latLng = d.geometry.coordinates.slice(0).reverse();
						this._featureOnMouseover = d.properties.id;
						this.setTooltipContent(this._getNodeInfo(d))
						    .openTooltip(latLng);
					}, 500);
				})
				.on("mouseout", () => {
					clearTimeout(this._mouseoverTimeout);
					this.closeTooltip();
					this._featureOnMouseover = null;
				})
				.on("click", d => {
					L.DomEvent.stopPropagation(d3.event);
					this._highlightHost(d.properties.id);
					this.closeTooltip();
					this._featureOnMouseover = null;
				})
				.on("dblclick", d => {
					L.DomEvent.stopPropagation(d3.event);
					this._displayHostHistoric(d.properties.id);
				});
		}

		this._repositionNodes(this._newNodes);

		const transition = this._newNodes
			.transition()
			.duration(duration)
			.ease(t => this._easeFunction(t, elapsed))
			.on("start", this._setTransitionStartTime())
			.style("opacity", d => this._getNodeOpacity(d));

		console.log("    [update][   add] transition.empty():", transition.empty());

		return transition.end()
			.then(() => {
				this._nodes = this._nodes.merge(this._newNodes);
				this._newNodes = null;
			});
	},
	// Remove nodes that doesn't match features data
	_removeOldNodes: function(duration, elapsed) {
		console.log("    [update][remove] start");

		let transition = this._nodes.exit()
			.transition()
			.duration(duration)
			.on("start", this._setTransitionStartTime())
			.style("opacity", 0)
			.each(d => {
				if (d.properties.id === this._featureOnMouseover) {
					this.closeTooltip();
					this._featureOnMouseover = null;
				}
			})
			.remove();

		console.log("    [update][remove] transition.empty():", transition.empty());

		return transition.empty() ? null : transition.end();
	},
	// Update nodes' position, by either moving or repositionning immediatly
	_updateNodesPosition: function(duration, elapsed) {
		console.log("    [update][smooth] start");
		const that = this;

		// Reposition nodes whose map coordinates didn't change
		const nodesNotMoving = this._nodes
			.select(function(d) { return that._nodeCoordinatesChanged(this, d) ? null : this; });
		this._repositionNodes(nodesNotMoving);

		// let edgesData = [];

		// Move nodes whose map coordinates changed
		const transition = this._nodes
			// Select the nodes that (still) have to move
			.select(function(d) { return that._nodeCoordinatesChanged(this, d) ? this : null; })
			// .each(function(d) {
			// 	edgesData.push([this._oldCoordinates, d.geometry.coordinates]);
			// })
			.transition()
			.duration(duration)
			.ease(t => this._easeFunction(t, elapsed))
			.on("start", this._setTransitionStartTime())
			.on("end", function(d) { this._oldCoordinates = d.geometry.coordinates; })
			.attrTween("cx", function(d) {
				return d3.interpolateNumber(that._projectPoint(this._oldCoordinates).x,
				                            that._projectPoint(d.geometry.coordinates).x);
			})
			.attrTween("cy", function(d) {
				return d3.interpolateNumber(that._projectPoint(this._oldCoordinates).y,
				                            that._projectPoint(d.geometry.coordinates).y);
			});

		console.log("    [update][smooth] transition.empty():", transition.empty());

		// const edges = this._selectMoveEdges();
		// this._drawEdges(edges, duration, elapsed);

		return transition.empty() ? null : transition.end();
	},
	// Update nodes' colors
	_updateNodesColor: function(duration) {
		console.log("    [infect][color] start");

		const transition = this._nodes
			.transition()
			.duration(duration)
			.on("start", this._setTransitionStartTime())
			.style("fill", d => this._getNodeFillColor(d))
			.style("stroke", d => this._getNodeStrokeColor(d))
			.style("opacity", d => this._getNodeOpacity(d));

		console.log("    [infect][color] transition.empty():", transition.empty());

		return transition.empty() ? null : transition.end();
	},
	// Reposition nodes to their old coordinates
	_repositionNodes: function(nodes) {
		const that = this;
		nodes
			.attr("cx", function() { return that._projectPoint(this._oldCoordinates).x })
			.attr("cy", function() { return that._projectPoint(this._oldCoordinates).y });
	},
	// Return the color of a node depending on its state
	_getNodeColor: function(node) {
		const state = node.properties.status;

		return this.options.nodeColorByState.hasOwnProperty(state) ?
			this.options.nodeColorByState[state] :
			this.options.defaultNodeColor;
	},
	// Return the fill color of a node whether it is an old feature
	_getNodeFillColor: function(node) {
		return node.properties._oldFeature ?
			'transparent' :
			this._getNodeColor(node);
	},
	// Return the stroke's color of a node whether it is an old feature
	_getNodeStrokeColor: function(node) {
		return node.properties._oldFeature ?
			this._getNodeColor(node) :
			this._historic === node.properties.id ?
				this.options.moveEdgeColor :
				'transparent';
	},
	// Return the node opacity depending it is related to the highlighted host
	_getNodeOpacity: function(node) {
		if (this._highlightedHost == null) return 1;
		if (this._highlightedHost.relatedHosts == null) return 0.25;

		return this._isHostRelatedToHighligtedHost(node) ? 1 : 0.25;
	},
	// Return formated node's info from data
	_getNodeInfo: function(node) {
		const date = new Date(node.properties._time);
		const time = `${date.toLocaleDateString()} ${date.toLocaleTimeString()} (${date.getTime()})`;

		let nodeInfo = '';
		nodeInfo += `Host: ${node.properties.id}<br>`;
		nodeInfo += `Time: ${time}<br>`;
		nodeInfo += `Coordinnates: ${node.geometry.coordinates}<br>`;
		nodeInfo += `Status: ${node.properties.status}<br>`;
		nodeInfo += `Infected: ${node.properties.infected}<br>`;
		nodeInfo += `Infected by: ${node.properties.infectedby}<br>`;
		nodeInfo += '';

		return nodeInfo;
	},
	// Return true if the given node's coordinates have changed ; false otherwise.
	// `node` is the DOM element, `host` is the data hold by d3js
	_nodeCoordinatesChanged: function(node, host) {
		return !this._sameCoordinates(node._oldCoordinates, host.geometry.coordinates);
	},
	// Check if two arrays of coordinates are the same
	_sameCoordinates: function(c1, c2) {
		return c1[0] === c2[0] && c1[1] === c2[1];
	},

	//---------------------------------------------------------------------------
	// * Infections edges methods
	//---------------------------------------------------------------------------
	// Add infections edges, set stroke and marker colors and set up the tooltip
	// to display the edge's info on mouseover event.
	_addInfectionsEdges: function() {
		this._infections = this._addEdges(this._infectionsData, "g #infection");

		this._infections
			.style("stroke", d => this._getInfectionEdgeColor(d))
			.attr('marker-end', d => {
				// Set the edge's end-marker with the appropriate color
				const color = this._getInfectionEdgeColor(d).substring(1);
				return `url(#marker_${color})`;
			})
			.on("mouseover", d => {
				// After a short time period, display edge's informations
				// within the Layer tooltip
				// Tooltip coordinates are set at the mouse position
				const latLng = this._map.mouseEventToLatLng(d3.event);
				this._mouseoverTimeout = setTimeout(() => {
					const date = new Date(d.time);
					let content = `<span style="color: red">${d.src.properties.id}</span>`;
					content += ` infected <span style="color: green">${d.dst.properties.id}</span>`
					if (typeof d.probability === 'number') {
						content += ` with ${d.probability}% of probability`;
					}
					content += ` on ${date.toLocaleDateString()}`;

					this
						.setTooltipContent(content)
						.openTooltip(latLng);
					this._featureOnMouseover = -1;
				}, 500);
			})
			.on("mouseout", () => {
				clearTimeout(this._mouseoverTimeout);
				this.closeTooltip();
				this._featureOnMouseover = null;
			});
	},
	// Draw infections edges normally (point => line)
	_drawInfectionsEdges: function(...args) {
		this._addInfectionsEdges();
		return this._drawEdges(this._infections, ...args);
	},
	// Draw infections edges in the reversed order (lines appear gradually)
	_drawInfectionsEdgesReversed: function(...args) {
		this._addInfectionsEdges();
		return this._drawEdgesReversed(this._infections, ...args);
	},
	// Clear infections edges normally (lines disappear gradually)
	_clearInfectionsEdges: function(...args) {
		return this._clearEdges(this._infections, ...args);
	},
	// Clear infections edges in the reversed order (line => point)
	_clearInfectionsEdgesReversed: function(...args) {
		return this._clearEdgesReversed(this._infections, ...args);
	},
	// Update infections edges color immediatly
	_updateInfectionsEdgesColor: function() {
		if (!this._infections || this._infections.empty()) return null;

		this._infections
			.style("stroke", d => this._getInfectionEdgeColor(d))
			.attr('marker-end', d => {
				// Set the edge's end-marker with the appropriate color
				const color = this._getInfectionEdgeColor(d).substring(1);
				return `url(#marker_${color})`;
			});
	},
	// Return the edge color depending on it's info
	_getInfectionEdgeColor: function(d) {
		return this._edgesColorScheme(d.probability || 1);
	},

	//---------------------------------------------------------------------------
	// * Moves edges methods
	//---------------------------------------------------------------------------
	// Add infections edges, set stroke and marker colors and set up the tooltip
	// to display the edge's info on mouseover event.
	_addMovesEdges: function() {
		console.log("addMovesEdges");
		console.log("this._movesData:", this._movesData);
		if (!this._movesData || this._movesData.length === 0) return null;
		this._moves = this._addEdges(this._movesData, "g #move");
		console.log("this._moves:", this._moves);
		this._moves
			.style("stroke", 'dodgerBlue')
			.attr('marker-end', 'url(#marker_1e90ff)')
			.on("mouseover", d => {
				// After a short time period, display edge's informations
				// within the Layer tooltip
				// Tooltip coordinates are set at the mouse position
				const latLng = this._map.mouseEventToLatLng(d3.event);
				this._mouseoverTimeout = setTimeout(() => {
					const date = new Date(d.time);
					const content = `<span style="color:green">#${d.src.properties.id}</span> moved at <span style="color:dodgerBlue">[${d.dst.geometry.coordinates}]</span> on ${date.toLocaleDateString()}`;

					this
						.setTooltipContent(content)
						.openTooltip(latLng);
					this._featureOnMouseover = -1;
				}, 500);
			})
			.on("mouseout", () => {
				clearTimeout(this._mouseoverTimeout);
				this.closeTooltip();
				this._featureOnMouseover = null;
			});;
	},
	// Draw infections edges normally (point => line)
	_drawMovesEdges: function(...args) {
		this._addMovesEdges();
		return this._drawEdges(this._moves, ...args);
	},
	// Draw moves edges in the reversed order (lines appear gradually)
	_drawMovesEdgesReversed: function(...args) {
		this._addMovesEdges();
		return this._drawEdgesReversed(this._moves, ...args);
	},
	// Clear moves edges normally (lines disappear gradually)
	_clearMovesEdges: function(...args) {
		return this._clearEdges(this._moves, ...args);
	},
	// Clear moves edges in the reversed order (line => point)
	_clearMovesEdgesReversed: function(...args) {
		return this._clearEdgesReversed(this._moves, ...args);
	},

	//---------------------------------------------------------------------------
	// * Edges methods
	//---------------------------------------------------------------------------
	// Add edges features to the given SVG group and return a selection on both
	// current and new edges
	_addEdges: function(data, group) {
		console.log("    [append-edges][start]");

		let edges = this._svg
			.select(group)
			.selectAll("path")
			.data(data);

		if (!edges.enter().empty()) {
			console.log("    [append-edges] Appending edges!");
			const newEdges = edges.enter().append("path");

			edges = edges
				.merge(newEdges)
				.data(data); // We must reset data because it was erased during the merge
		}

		return edges;
	},
	// Draw the given edges normally (point => line)
	_drawEdges: function(edges, duration, elapsed) {
		console.log("    [draw-edges] start");
		if (!edges || edges.empty()) return null;

		const transition = edges
			.style("opacity", d => this._getEdgeOpacity(d))
			.transition()
			.duration(duration)
			.ease(t => this._easeFunction(t, elapsed))
			.on("start", this._setTransitionStartTime())
			.attrTween('d', d => {
				// Gradually draw a path from the infector to the infected node.
				const paths = this._getEdgePaths(d);
				return d3.interpolateString(paths.point, paths.line);
			});

		return transition.empty() ? null : transition.end();
	},
	// Draw the given edges in reversed order (lines appear gradually)
	_drawEdgesReversed: function(edges, duration, elapsed) {
		console.log("    [draw-edges-r] start");
		if (!edges || edges.empty()) return null;

		const transition = edges
			.attr('d', d => this._getEdgePaths(d).line)
			.style('opacity', 0)
			.transition()
			.duration(duration)
			.ease(t => this._easeFunction(t, elapsed))
			.on("start", this._setTransitionStartTime())
			.style('opacity', d => this._getEdgeOpacity(d));

		return transition.empty() ? null : transition.end();
	},
	// Clear the given edges normally (lines disappear gradually)
	_clearEdges: function(edges, duration, elapsed) {
		console.log("    [clear-edges] start");
		if (!edges || edges.empty()) return null;

		if (duration === 0) {
			edges.remove();
			return null;
		}

		const transition = edges
			.transition()
			.duration(duration)
			.ease(t => this._easeFunction(t, elapsed))
			.on("start", this._setTransitionStartTime())
			.style("opacity", 0)
			.remove();

		return transition.empty() ? null : transition.end();
	},
	// Clear the given edges in the reversed order (line => point)
	_clearEdgesReversed: function(edges, duration, elapsed) {
		console.log("    [clear-edges-r] start");
		if (!edges || edges.empty()) return null;

		if (duration === 0) {
			edges.remove();
			return null;
		}

		const transition = edges
			.transition()
			.duration(duration)
			.ease(t => this._easeFunction(t, elapsed))
			.on("start", this._setTransitionStartTime())
			.attrTween('d', d => {
				// Draw a path from the infector to the infected nodes and reduce it
				// gradually to a point at the infector position
				const paths = this._getEdgePaths(d);
				return d3.interpolateString(paths.line, paths.point);
			})
			.remove();

		return transition.empty() ? null : transition.end();
	},
	// Reposition edges immediatly
	_repositionEdges: function(edges) {
		if (!edges || edges.empty()) return;
		edges.attr('d', d => this._getEdgePaths(d).line);
	},
	// Return the starting and ending paths of an edge between a source host
	// and the destination one.
	// The starting path is equivalent to a point at the origin position ;
	// the second one goes from the origin host to the destination one.
	// The path is shorten to not overlap on nodes.
	_getEdgePaths: function(edge) {
		// Project geographical coordinates into pixel coordinates
		let src = this._projectPoint(edge.src.geometry.coordinates);
		    dst = this._projectPoint(edge.dst.geometry.coordinates);

		// Calculate the length of the final path
		const length = Math.sqrt(Math.pow(src.x - dst.x, 2) +
		                         Math.pow(src.y - dst.y, 2));

		// Shift the starting point of the path to not overlap the source node
		const fraction1 = this.options.circleRadius / length,
		      x1 = src.x + (dst.x - src.x) * fraction1,
		      y1 = src.y + (dst.y - src.y) * fraction1;

		// Shift the destination point of the path to not overlap the destination node
		const fraction2 = (length - 2 * this.options.circleRadius) / length,
		      x2 = src.x + (dst.x - src.x) * fraction2,
		      y2 = src.y + (dst.y - src.y) * fraction2;

		src = {x: x1, y: y1};
		dst = {x: x2, y: y2};

		// Starting path = point at the source node
		let pointPath = d3.path();
		pointPath.moveTo(src.x, src.y);
		pointPath.lineTo(src.x, src.y);
		pointPath = pointPath.toString();

		// Ending path = line between source and destination nodes
		let linePath = d3.path();
		linePath.moveTo(src.x, src.y);
		linePath.lineTo(dst.x, dst.y);
		linePath = linePath.toString();

		return {point: pointPath, line: linePath};
	},
	// Check if the edge is related to the highlighted host
	_isEdgeRelatedToHighlightedHost: function(edge) {
		return edge.src.properties.id === this._highlightedHost.id ||
		       edge.dst.properties.id === this._highlightedHost.id;
	},
	// Return the edge opacity depending on the highlighted host
	_getEdgeOpacity: function(edge) {
		if (this._highlightedHost == null) return 1;
		return this._isEdgeRelatedToHighlightedHost(edge) ? 1 : 0.25;
	},

	//---------------------------------------------------------------------------
	// * Features data methods
	//---------------------------------------------------------------------------
	// Return the features to display at the TimeDimension's current time
	_getCurrentFeatures: function() {
		const currentTime = this._timeDimension.getCurrentTime();

		let currentFeatures = [];
		this._baseLayer.eachLayer(layer => {
			let feature = this._lineStringToPoint(layer.feature, currentTime);
			if (feature) currentFeatures.push(feature);
		});

		return currentFeatures;
	},
	// Transform a LineString feature to a Point and return it.
	// The feature's fields are filtered according to the given time.
	// The original feature properties are not modified.
	// Return null if the feature is invalid or has no data for the
	// current time *yet* (i.e. the first time at which we have info
	// on the feature is *after* the current time).
	// If `keepOldFeatures` option is set to true and we don't have info
	// *anymore* on the feature (i.e. the last known info is *before*
	// the current time), returns the feature at the last known time.
	// Otherwise, returns null.
	_lineStringToPoint: function(feature, time) {
		if (feature.geometry.type === 'Point' &&
		    feature.hasOwnProperty('property') &&
		    feature.properties._time <= time) {
			return feature;
		}
		if (feature.geometry.type != 'LineString') return null;

		const featureTimes = this._getFeatureTimes(feature);
		if (featureTimes.length === 0) return null;

		let index = null,
		    len = featureTimes.length,
		    oldFeature = false;

		if (featureTimes[0] > time) return null;
		if (!this.options.keepOldFeatures && featureTimes[len-1] < time) return null;

		if (this.options.keepOldFeatures && featureTimes[len-1] < time) {
			// We don't have info for this feature *anymore*, so we use the last known info
			index = len - 1;
			oldFeature = true;
		} else {
			// Pick the most recent info compared to the given time
			index = -1;
			for (const t of featureTimes) {
				if (t > time) break;
				index++;
			}
		}

		if (index < 0) return null;

		const c = feature.geometry.coordinates[index];
		const p = JSON.parse(JSON.stringify(feature.properties));

		p.status = p.hasOwnProperty('status') ? p.status[index] : 'default';
		p.infected = [];

		if (oldFeature || time > featureTimes[index]) {
			p.infectedby = [];
		} else if (p.hasOwnProperty('infectedby') &&
		           p.infectedby.hasOwnProperty(index)) {
			p.infectedby = p.infectedby[index];
			if (p.hasOwnProperty('probabilities') &&
			    p.probabilities.hasOwnProperty(index)) {
				p.probabilities = p.probabilities[index];
			}
		}

		p._time = time;
		p._oldFeature = oldFeature;

		return {
			type: 'Feature',
			properties: p,
			geometry: {
				type: 'Point',
				coordinates: c
			}
		};
	},
	// Determine the infections edges' infos from GeoJSON features
	// Also adds the infected host's id to the infector's `infected` list.
	// Return a list formated like this:
	// [
	//  	{
	//  		'probability': float,
	//  		'src': infector host data
	//  		'dst': infected host data
	//  	},
	//  	...
	// ]
	_determineInfections: function(features) {
		console.log("    [infect][getInfections] start");

		// Returns the most recent host instance compared to the given time
		function getMostRecentHost(id, time) {
			const hosts = features.filter(f => f.properties.id === id);

			if (hosts.length === 0) return null;

			let mostRecentHost = null;

			hosts.forEach(function(f) {
				if (f.properties._time <= time &&
				    (mostRecentHost == null ||
				     f.properties._time > mostRecentHost.properties._time)) {
					mostRecentHost = f;
				}
			});

			return mostRecentHost;
		}

		// Add an infection relation between two hosts
		function addInfection(fp, infectorId, i = null) {
			let infection = {};

			infection.src = infectorId;
			infection.dst = fp.id;
			infection.time = fp._time;

			if (typeof fp.probabilities === 'number') {
				infection.probability = fp.probabilities;
			} else if (Number.isInteger(i) && fp.probabilities) {
				infection.probability = fp.probabilities[i];
			}

			infection.src = getMostRecentHost(infectorId, infection.time);
			infection.dst = getMostRecentHost(fp.id, infection.time);

			if (infection.src != null && infection.dst != null) {
				infections.push(infection);
				infection.src.properties.infected.push(fp.id);
			} else {
				if (infection.src == null) {
					console.warn(`Could not find infector host '${infectorId}'`);
				}
				if (infection.dst == null) {
					console.warn(`Could not find infected host '${fp.id}'`);
				}
			}
		}

		let infections = [], // the list we want to build
		    hostsIds = [];   // list of the ids of both infected and infectors hosts
		                     // It is not necessary but used for convenience.

		features.forEach(function(feature) {
			const fp = feature.properties;

			if (fp._oldFeature ||
			    !fp.hasOwnProperty('infectedby') ||
			    fp.infectedby.length === 0) {
				return;
			}

			if (Number.isInteger(fp.infectedby)) {
				addInfection(fp, fp.infectedby);
			} else if (fp.infectedby.length) {
				fp.infectedby.forEach(function(infectorId, i) {
					addInfection(fp, infectorId, i);
				});
			}
		});

		console.log("    [infect][getInfections] infections:", infections);

		return infections;
	},
	// Select the already existing nodes elements and set data
	_selectNodes: function(geojsonFeatures) {
		const data = geojsonFeatures.filter(d => d.geometry.type === "Point");

		return this._svg
			.select("g #points")
			.selectAll("circle")
			.data(data, d => d.properties.id);
	},

	//---------------------------------------------------------------------------
	// * Host highlight
	//---------------------------------------------------------------------------
	// Highlight host relations.
	// Return a message of validity (empty if no error)
	_highlightHost: function(hostId) {
		if (this._historic != null) return `Can't highlight a host during the historic view.`;

		let features = this._baseLayer.getLayers().map(l => l.feature);

		// Check if the host exists in the entire time
		const host = features.find(f => f.properties.id === hostId);
		if (!host) {
			const warningMessage = `Could not find a host with id ${hostId}`;
			console.warn(warningMessage);
			return warningMessage;
		}

		this._selectHost(hostId);

		const currentHost = this._currentFeatures.find(f => f.properties.id === hostId);

		this._highlightedHost = {
			id: hostId,
			relatedHosts: this._getRelatedHosts(currentHost)
		};

		this._highlightRelatedHosts();
		this._highlightRelatedInfections();

		this.fire('highlight', {
			hostId: hostId
		});

		this._map.once("click", this._reinitializeHighlightedHost, this);

		return '';
	},
	// Partially hide all the nodes then highlight the ones related to
	// the selected one
	_highlightRelatedHosts: function() {
		if (this._nodes == null || this._nodes.empty()) return;

		this._nodes.style('opacity', 0.25);

		this._nodes
			.filter(d => this._isHostRelatedToHighligtedHost(d))
			.style('opacity', 1);
	},
	// Partially hide all the edges then highlight the ones related to
	// the selected host
	_highlightRelatedInfections: function() {
		if (this._infections == null || this._infections.empty()) return;

		this._infections.style('opacity', 0.25);

		this._infections
			.filter(d => this._isEdgeRelatedToHighlightedHost(d))
			.style('opacity', 1);
	},
	// Reinitialize highlighted host
	_reinitializeHighlightedHost: function() {
		if (this._highlightedHost == null) return;

		this._highlightedHost = null;

		if (this._nodes != null && !this._nodes.empty()) {
			this._nodes.style('opacity', 1);
		}

		if (this._infections != null && !this._infections.empty()) {
			this._infections.style('opacity', 1);
		}

		this.fire('highlight', {
			hostId: null
		});
	},
	// Check if the host is related to the highlighted one
	_isHostRelatedToHighligtedHost: function(host) {
		return this._highlightedHost.relatedHosts.includes(host.properties.id);
	},

	//---------------------------------------------------------------------------
	// * Host historic
	//---------------------------------------------------------------------------
	// Display host historic
	_displayHostHistoric: function(hostId) {
		console.log(`----- Display node historic of ${hostId}`);

		if (hostId === '' || this._historic === hostId) return;

		if (this._player.isPlaying()) {
			this._player.stop();
		}

		this._resetAnimation();

		this._historic = hostId;
		this._selectHost(hostId);

		this._reinitializeHighlightedHost();

		const historicFeatures = this._getHostHistoricFeatures(hostId);
		if (!historicFeatures.length) {
			const warningMessage = `Could not find host with id ${hostId}`;
			console.warn(warningMessage);
			return warningMessage;
		}

		this._historic = hostId;

		this._nodes = this._selectHistoricNodes(historicFeatures);
		const hostNodes = historicFeatures.filter(d => d.properties.id === hostId);
		this._movesData = this._getHostMovesEdges(hostNodes);
		this._infectionsData = this._determineInfections(historicFeatures);

		this._updateMapView(historicFeatures);

		this._addNewNodes(0);

		this._clearMovesEdges(0);
		this._drawMovesEdges(0);

		this._clearInfectionsEdges(0);
		this._drawInfectionsEdges(0);

		this.fire('historic', {
			hostId: hostId
		});

		return '';
	},
	// Clear host historic
	_clearHostHistoric: function() {
		if (this._historic == null) return;

		this._historic = null;
		this.fire('historic', {
			hostId: null
		});

		this._resetAnimation();
		this._initializeAnimation();
		this._skipAnimation = true;
		this._clearMovesEdges(0);
		this._runAnimation();
	},
	// Return all the host features and their related hosts
	_getHostHistoricFeatures: function(hostId) {
		let features = this._baseLayer.getLayers().map(l => l.feature);

		const hostFeature = features.find(f => hostId === f.properties.id);
		if (!hostFeature) {
			console.warn(`Warning: could not find host ${hostId}`);
			return [];
		}

		// Add the selected host features for each distinct time we have info on him
		const hostFeatures = hostFeature.featureTimes
			.map(t => this._lineStringToPoint(hostFeature, t));

		const hp = hostFeature.properties;

		// Add the hosts at each time they are related to the selected host
		let hosts = features.map(f => {
				const fp = f.properties;
				// Check when the host is related to the selected host
				const relatedF = fp.infectedby.map(arr => arr.includes(hp.id)),
				      relatedH = hp.infectedby.map(arr => arr.includes(fp.id));

				if (!relatedF.includes(true) && !relatedH.includes(true)) {
					return null;
				}

				// Filter times where the host is related to the selected host
				const timesF = f.featureTimes.filter((t,i) => relatedF[i]),
				      timesH = hp.linestringTimestamps.filter((t,i) => relatedH[i]),
				      times = [...new Set(timesF.concat(timesH))];

				// Retrieve the host info filtered at each distinct time
				return times.map(t => this._lineStringToPoint(f, t));
			})
			.flat();

		// Add time inside each feature's id
		hosts.map(f => {
			if (!f) return;
			const fp = f.properties;
			// This new field will be used as the key to help d3 identify nodes
			fp.timeId = fp.id + '-' + fp._time;
			// Filter the infectedby field to have only the selected host
			fp.infectedby = fp.infectedby.filter(id => id === hp.id);
		});

		hosts = hosts.concat(hostFeatures);
		hosts = hosts.filter(h => h != null);

		console.log("hosts:", hosts);

		return hosts;
	},
	// Select the already existing nodes elements and set data
	_selectHistoricNodes: function(geojsonFeatures) {
		const data = geojsonFeatures.filter(function(d) {
			return d.geometry.type === "Point"
		});
		return this._svg
			.select("g #points")
			.selectAll("circle")
			.data(data, d => d.properties.timeId);
	},
	// Return the move edges for each host coordinates changing
	_getHostMovesEdges: function(hostNodes) {
		let edges = [],
		    previousHost = null;

		hostNodes.map(d => {
			if (previousHost == null) {
				previousHost = d;
			} else if (edges.length === 0 || !this._sameCoordinates(previousHost.geometry.coordinates, d.geometry.coordinates)) {
				const edge = {
					src: previousHost,
					dst: d,
					time: d.properties._time
				};
				edges.push(edge);
				previousHost = d;
			}
		});

		return edges;
	},
	// Return the related hosts' ids of a host
	_getRelatedHosts: function(host) {
		if (host == null) return null;

		return [host.properties.id]
			.concat(host.properties.infected)
			.concat(host.properties.infectedby);
	},

	//---------------------------------------------------------------------------
	// * Leaflet utils
	//---------------------------------------------------------------------------
	// Project GeoJSON coordinates onto the Leaflet map and return the
	// corresponding coordinates.
	_projectPoint: function(c) {
		// Be careful: GeoJSON and Leaflet coordinates are inversed
		return this._map.latLngToLayerPoint(new L.LatLng(c[1], c[0]));
	},
	// Update the map view to cover both old and new features
	_focusOnFeatures: function(duration, elapsed) {
		const features = this._oldFeatures.concat(this._currentFeatures);
		this._updateMapView(features);
	},
	// Update the map position & zoom to cover the given features.
	// Return a Promise which resolves when the map view has been
	// updated, or `null` if it doesn't need to be updated.
	_updateMapView: function(features) {
		const newBounds = this._getFeaturesBounds(features),
		      zoomDifference = this._map.getZoom() - this._map.getBoundsZoom(newBounds);

		// Don't reset the map view if the new bounds are more or less equivalent
		// (i.e. if the view clearly shows all the features)
		if (this._map.getBounds().contains(newBounds) &&
		    Math.abs(zoomDifference) < 1) {
			return;
		}

		const promise = new Promise(resolve => this._map.once("moveend", resolve));

		if (newBounds.getNorthEast().equals(newBounds.getSouthWest())) {
			this._map.panTo(newBounds.getNorthEast());
		} else {
			// Fit the new bounds with some padding in pixels to prevent the features
			// to be overlayed by controls
			this._map.fitBounds(newBounds, {
    		paddingTopLeft: [50, 30],
    		paddingBottomRight: [100, 100]
  		});
		}

		return promise;
	},
	// Return a `latLngBounds` containing all the GeoJSON features given.
	// If bounds are invalid (if there are no features or invalid ones, for
	// example), return the current map bounds.
	_getFeaturesBounds: function(features) {
		const bounds = L.geoJson({
			type: "FeatureCollection",
			features: features
		}).getBounds();

		return bounds.isValid() ? bounds : this._map.getBounds();
	},
});

L.maptt = function(layer, options) {
	return new L.MapTT(layer, options);
};
