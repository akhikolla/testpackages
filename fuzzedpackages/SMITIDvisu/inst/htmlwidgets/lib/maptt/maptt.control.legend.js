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

L.Control.MapTT.Legend = L.Control.MapTT.extend({

	options: {
		position: 'bottomright'
	},

	//---------------------------------------------------------------------------
	// * Basic functions
	//---------------------------------------------------------------------------
	initialize: function(options) {
		L.Control.MapTT.prototype.initialize.call(this, options);
	},
	onAdd: function(map) {
		L.Control.MapTT.prototype.onAdd.call(this, map);

		return this._maptt ? this._addLegendControl() : null;
	},

	//---------------------------------------------------------------------------
	// * Legend control
	//---------------------------------------------------------------------------
	_addLegendControl: function() {
		let container = L.DomUtil.create('div', 'maptt-control-content');
		L.DomEvent.disableClickPropagation(container);

		let title = L.DomUtil.create('h4', '', container);
		title.innerHTML = 'Legend';

		this._addStatusSection(container);
		this._addEdgeColorsSection(container);

		return this._wrapInControlContainer(container);
	},
	_addStatusSection: function(container) {
		const stateColors = this._maptt.options.nodeColorByState;

		if (!stateColors || Object.keys(stateColors).length === 0) return;

		let section = L.DomUtil.create('div', 'options-section', container);

		let title = L.DomUtil.create('h4', '', section);
		title.innerHTML = "Hosts' status";

		let div = L.DomUtil.create('div', 'legend', section);
		for (var state in stateColors) {
			div.innerHTML += `<i style="background:${stateColors[state]}"></i>${state}<br>`;
		}
	},
	_addEdgeColorsSection: function(container) {
		let section = L.DomUtil.create('div', 'options-section', container);

		let title = L.DomUtil.create('h4', '', section);
		title.innerHTML = "Edges' colors";

		const opt = this._maptt.options;
		let div = L.DomUtil.create('div', 'legend', section);
		div.innerHTML += `<i style="background:${opt.moveEdgeColor}"></i>move<br>`;
		div.innerHTML += `<i style="background:${opt.color1}"></i>infection probability: 0%<br>`;
		div.innerHTML += `<i style="background:${opt.color2}"></i>infection probability: 100%<br>`;
	},
});

L.control.maptt.legend = function(options) {
	return new L.Control.MapTT.Legend(options);
};
