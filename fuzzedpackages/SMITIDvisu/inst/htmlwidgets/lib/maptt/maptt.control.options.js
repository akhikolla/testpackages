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

L.Control.MapTT.Options = L.Control.MapTT.extend({

	options: {
		position: 'topright'
	},

	//---------------------------------------------------------------------------
	// * Basic functions
	//---------------------------------------------------------------------------
	initialize: function(options) {
		L.Control.MapTT.prototype.initialize.call(this, options);
	},
	onAdd: function(map) {
		L.Control.MapTT.prototype.onAdd.call(this, map);

		return this._maptt ? this._addOptionsControls() : null;
	},
	setMaptt: function(maptt) {
		L.Control.MapTT.prototype.setMaptt.call(this, maptt);
		this._maptt.on("highlight", this._onHighlight, this);
		this._maptt.on("historic", this._onHistoric, this);
	},
	_addCheckbox: function(name, checked, container, onClick) {
		let holder = L.DomUtil.create('label', 'checkbox', container),
		    input = L.DomUtil.create('input', '', holder),
		    label = L.DomUtil.create('span', '', holder);

		label.innerHTML = name;

		input.type = 'checkbox';
		input.checked = checked;
		L.DomEvent.on(input, 'click', onClick, input);
	},

	//---------------------------------------------------------------------------
	// * Events
	//---------------------------------------------------------------------------
	_onHighlight: function(data) {
		if (data.hostId == null) {
			this._highlightInput.value = '';
		} else {
			this._highlightInput.value = data.hostId;
		}
	},
	_onHistoric: function(data) {
		if (data.hostId == null) {
			this._highlightInput.disabled = false;
		} else {
			this._historicInput.value = data.hostId;
			this._highlightInput.disabled = true;
		}
	},

	//---------------------------------------------------------------------------
	// * Options control
	//---------------------------------------------------------------------------
	// Add the options controls to the map
	_addOptionsControls: function() {
		const that = this;

		let container = L.DomUtil.create('div', 'maptt-control-content');
		L.DomEvent.disableClickPropagation(container);

		this._addOptionsSection(container);
		this._addHighlightSection(container);
		this._addHistoricSection(container);

		return this._wrapInControlContainer(container);
	},
	_addOptionsSection: function(container) {
		const that = this;

		let section = L.DomUtil.create('div', 'options-section', container);

		let title = L.DomUtil.create('h4', '', section);
		title.innerHTML = 'Options';

		this._addCheckbox('Auto-focus', this._maptt.options.keepOldFeatures, section, function() {
			that._maptt.options.autoFocus = this.checked;
		});
		this._addCheckbox('Old features', this._maptt.options.keepOldFeatures, section, function() {
			that._maptt.options.keepOldFeatures = this.checked;

			if (!that._maptt._isAnimationRunning() && that._maptt._historic == null) {
				that._maptt._initializeAnimation();
				that._maptt._runAnimation();
			}
		});
	},
	// Add historic section controls
	_addHistoricSection: function(container) {
		let historicSection = L.DomUtil.create('div', 'options-section', container);

		let historicLabel = L.DomUtil.create('h4', '', historicSection);
		historicLabel.innerHTML = 'Host historic';

		let line = L.DomUtil.create('div', '', historicSection);
		line.style.marginTop = '4px';
		line.style.marginBottom = '8px';

		let historicInput = L.DomUtil.create('input', '', line);

		const askHostHistoric = () => {
			const hostId = historicInput.value;

			if (hostId == '') {
				this._maptt._clearHostHistoric();
			} else {
				const validityMsg = this._maptt._displayHostHistoric(hostId);
				historicInput.setCustomValidity(validityMsg);
			}
		}

		historicInput.placeholder = 'Host ID';
		historicInput.style.width = '80px';
		L.DomEvent.on(historicInput, 'change', askHostHistoric);
		this._historicInput = historicInput;

		let ok = L.DomUtil.create('button', '', line);
		ok.innerHTML = 'OK';
		ok.style.marginLeft = '4px';
		L.DomEvent.on(ok, 'click', askHostHistoric);

		let reset = L.DomUtil.create('button', '', historicSection);
		reset.innerHTML = 'Reset';
		L.DomEvent.on(reset, 'click', (event) => {
			historicInput.value = '';
			this._maptt._clearHostHistoric();
		});
	},
	// Add highlight section controls
	_addHighlightSection: function(container) {
		let section = L.DomUtil.create('div', 'options-section', container);

		let label = L.DomUtil.create('h4', '', section);
		label.innerHTML = 'Highlighted host';

		let line = L.DomUtil.create('div', '', section);
		line.style.marginTop = '4px';
		line.style.marginBottom = '8px';

		let highlightInput = L.DomUtil.create('input', '', line);

		const askHostHighlight = () => {
			const hostId = highlightInput.value;
			if (hostId == '') {
				highlightInput.value = '';
				this._maptt._reinitializeHighlightedHost();
			} else {
				highlightInput.setCustomValidity(this._maptt._highlightHost(hostId));
			}
		}

		highlightInput.placeholder = 'Host ID';
		highlightInput.style.width = '80px';
		L.DomEvent.on(highlightInput, 'change', askHostHighlight);
		this._highlightInput = highlightInput;

		let ok = L.DomUtil.create('button', '', line);
		ok.innerHTML = 'OK';
		ok.style.marginLeft = '4px';
		L.DomEvent.on(ok, 'click', askHostHighlight);

		let reset = L.DomUtil.create('button', '', section);
		reset.innerHTML = 'Reset';
		L.DomEvent.on(reset, 'click', (event) => {
			this._maptt._reinitializeHighlightedHost();
		});
	},
});

L.control.maptt.options = function(options) {
	return new L.Control.MapTT.Options(options);
};
