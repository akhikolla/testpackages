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

	name: "maptt",

	type: "output",

	factory: function(el, width, height) {

		console.log("Inside factory");
		console.log("el:", el);
		console.log("width:", width);
		console.log("height:", height);

		//=====================================================
		// * Initialize Leaflet Map
		//=====================================================

		var map = L.map(el.id, {
			zoomSnap: 0.5
		}).setView([43.91436, 4.88260], 17);
		el.setAttribute("height", "100%");
		el.setAttribute("width", "100%");

		// Ajoute la carte en fond à partir d'OpenStreetMap
		L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
			attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors'
		}).addTo(map);

		// // Ajoute un polygone
		// var polygon = L.polygon([
		// 	[43.911853, 4.879099],
		// 	[43.911853, 4.885933],
		// 	[43.917085, 4.885933],
		// 	[43.917085, 4.879099],
		// ]).addTo(map);

		//=====================================================
		// * Initialize Leaflet.TimeDimension
		//=====================================================

		// Initialize TimeDimension
		var timeDimension = new L.MyTimeDimension({
			period: 'P1D',
		});
		// Share TimeDimension object between all the layers
		map.timeDimension = timeDimension;

		// Initialize TimeDimension Player
		var player = new L.TimeDimension.MyPlayer({
			transitionTime: 700,
			loop: true,
			startOver: false
		}, timeDimension);

		// Set TimeDimension controller options
		var timeDimensionControlOptions = {
			player: player,
			timeDimension: timeDimension,
			position: 'bottomleft',
			autoPlay: true,
			minSpeed: 0.1,
			maxSpeed: 3.0,
			timeSliderDragUpdate: true,
			playReverseButton: true,
		};

		// Initialize TimeDimension Controller
		var timeDimensionControl = new L.Control.MyTimeDimension(timeDimensionControlOptions);
		// Add TimeDimension Controller to the map
		map.addControl(timeDimensionControl);

		//=====================================================
		// * Initialize MapTT Controls
		//=====================================================

		const mapttOptions = L.control.maptt.options(),
		      mapttGradient = L.control.maptt.gradient(),
		      mapttLegend = L.control.maptt.legend();

		//=====================================================
		// * Initialize MapTT
		//=====================================================
		var maptt = null;
		function initializeMaptt(args) {

			// Ensure that there is only one MapTT
			if (maptt != null) {
				console.log("Removing MapTT from map!");
				maptt.remove();
			}

			console.log("args:", args);

			geoJsonFeatures = JSON.parse(args.geoJson);
			geoJsonFeatures['features'] = geoJsonFeatures['features'];

			console.log("Formated geoJsonFeatures:", geoJsonFeatures);

			let onSelectHost = null;
			if (HTMLWidgets.shinyMode) {
				onSelectHost = function(hostId) {
					Shiny.onInputChange('mapttSelectedHost', hostId);
				}
			}

			var geoJsonLayer = L.geoJson(geoJsonFeatures);
			console.log("geoJsonLayer:", geoJsonLayer);

			// Replace TimeDimension's times with this layer's times
			args.updateTimeDimension = true;
			args.updateTimeDimensionMode = 'replace';

			args.player = player;
			args.onSelectHost = onSelectHost;

			maptt = L.maptt(geoJsonLayer, args);
			maptt.addTo(map);

			if (args.optionsControl) {
				mapttOptions.setMaptt(maptt);
				mapttOptions.addTo(map);
			} else {
				mapttOptions.remove();
			}

			if (args.gradientControl) {
				mapttGradient.setMaptt(maptt);
				mapttGradient.addTo(map);
			} else {
				mapttGradient.remove();
			}

			if (args.legend) {
				mapttLegend.setMaptt(maptt);
				mapttLegend.addTo(map);
			} else {
				mapttLegend.remove();
			}
		}

		//=====================================================
		// * Return to Shiny
		//=====================================================

		return {
			renderValue: (geoJsonFeatures) => initializeMaptt(geoJsonFeatures),
			resize: (width, height) => {},
			getMap: () => map,
			getTd: () => timeDimension,
			getPlayer: () => player,
			getTdControl: () => timeDimensionControl,
			getMaptt: () => maptt
		};
	}
});

/** getMaptt
* return object from id
* @param id : widget name (output name in server)
*/
function getMaptt(id) {
  const widget = HTMLWidgets.find("#" + id);
  const maptt = widget.getMaptt();
  return(maptt);
}

/** selectHost
 * add a custom message handler to select a MapTT host
 * client side
 */
if(typeof Shiny !== 'undefined' && Shiny != null) {
  Shiny.addCustomMessageHandler("selectHost", function(message) {
    const maptt = getMaptt(message.id);
    const hostId = message.hostId[0];
    maptt._displayHostHistoric(hostId);
  });
}
