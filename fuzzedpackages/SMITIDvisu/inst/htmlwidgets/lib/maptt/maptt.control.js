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

L.Control.MapTT = L.Control.extend({
	initialize: function(options) {
		console.log("maptt control: initialize");
		L.Control.prototype.initialize.call(this, options);
		if (!options) {
			options = {};
		}
		if (options.maptt != null) this.setMaptt(options.maptt);
	},
	onAdd: function(map) {
		this._map = map;
	},
	_wrapInControlContainer: function(content) {
		let wrapper = L.DomUtil.create('div', 'maptt-control leaflet-control');
		L.DomEvent.disableClickPropagation(wrapper);
		L.DomEvent.disableScrollPropagation(wrapper);
		L.DomEvent.on(wrapper, 'click', () => {
			L.DomUtil.addClass(wrapper, 'maptt-control-expanded');
		});

		let toggle = L.DomUtil.create('div', 'maptt-control-toggle', wrapper);

		wrapper.appendChild(content);

		let closeButton = L.DomUtil.create('button', 'maptt-control-close-button', content);
		closeButton.type = 'button';
		closeButton.innerHTML = '&times;';
		L.DomEvent.on(closeButton, 'click', (event) => {
			event.stopPropagation();
			L.DomUtil.removeClass(wrapper, 'maptt-control-expanded');
		});

		return wrapper;
	},
	setMaptt(maptt) {
		this._maptt = maptt;
	},
});

L.control.maptt = function(options) {
	return new L.Control.MapTT(options);
};
