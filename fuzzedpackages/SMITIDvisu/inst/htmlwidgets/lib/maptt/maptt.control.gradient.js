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


L.Control.MapTT.Gradient = L.Control.MapTT.extend({

	options: {
		position: 'topright'
	},

	//---------------------------------------------------------------------------
	// * Basic functions
	//---------------------------------------------------------------------------
	initialize: function(options) {
		L.Control.MapTT.prototype.initialize.call(this, options);

		if (this._maptt != null) {
			this.setGradient();
		}
	},
	onAdd: function(map) {
		L.Control.MapTT.prototype.onAdd.call(this, map);
		return this._maptt ? this._addGradientControls() : null;
	},
	setMaptt: function(maptt) {
		L.Control.MapTT.prototype.setMaptt.call(this, maptt);
		this.gradient = this._maptt.options;
		this.setGradient();
	},
	setGradient: function() {
		this.gradient.minColor = parseInt(this.gradient.minColor) || 1;
		this.gradient.maxColor = parseInt(this.gradient.maxColor) || 10;
		this.gradient.nbColors = parseInt(this.gradient.nbColors) || 10;
		this.gradient.minWeight = parseFloat(this.gradient.minWeight) || 0;
		this.gradient.maxWeight = parseFloat(this.gradient.maxWeight) || 1;
		this.gradient.weightStep = parseFloat(this.gradient.weightStep) || 0.1;
	},

	//---------------------------------------------------------------------------
	// * Gradient controls
	//---------------------------------------------------------------------------
	_addGradientControls: function() {
		let container = L.DomUtil.create('div', 'maptt-control-content');
		this._addColorSchemeSection(container);
		this._addWeightsSection(container);

		return this._wrapInControlContainer(container);
	},
	_addColorSchemeSection: function(container) {
		let section = L.DomUtil.create('div', 'options-section', container);

		let label = L.DomUtil.create('h4', '', section);
		label.innerHTML = 'Color Scheme';

		let colorSchemeContainer = L.DomUtil.create('div', 'doubleRangeContainerr', section);

		let minColor = L.DomUtil.create('input', 'value', colorSchemeContainer);
		minColor.type = 'number';
		minColor.placeholder = 'min';
		minColor.min = 1;
		minColor.max = this.gradient.nbColors;
		minColor.value = 1;
		L.DomEvent.on(minColor, 'change', this._onMinColorChange, this);

		let doubleRange = L.DomUtil.create('div', 'doubleRange', colorSchemeContainer);
		const color1 = this.gradient.color1,
		      color2 = this.gradient.color2,
		      bg = `linear-gradient(0.25turn, ${color1}, ${color2})`;
		doubleRange.style.background = bg;

		let colorRange1 = L.DomUtil.create('input', 'slider', doubleRange);
		colorRange1.type = 'range';
		colorRange1.min = 1;
		colorRange1.max = this.gradient.nbColors;
		colorRange1.value = 1;
		L.DomEvent.on(colorRange1, 'input', this._onColorRangeUpdate, this);

		let colorRange2 = L.DomUtil.create('input', 'slider', doubleRange);
		colorRange2.type = 'range';
		colorRange2.min = 1;
		colorRange2.max = this.gradient.nbColors;
		colorRange2.value = this.gradient.nbColors;
		L.DomEvent.on(colorRange2, 'input', this._onColorRangeUpdate, this);

		let maxColor = L.DomUtil.create('input', 'value', colorSchemeContainer);
		maxColor.type = 'number';
		maxColor.placeholder = 'max';
		maxColor.min = 1;
		maxColor.max = this.gradient.nbColors;
		maxColor.value = this.gradient.nbColors;
		L.DomEvent.on(maxColor, 'change', this._onMaxColorChange, this);

		let nbColorsLabel = L.DomUtil.create('label', '', section);
		nbColorsLabel.innerHTML = 'Number of colors';

		let nbColors = L.DomUtil.create('input', 'nbColors', nbColorsLabel);
		nbColors.type = 'number';
		nbColors.min = 2;
		nbColors.value = this.gradient.nbColors;
		L.DomEvent.on(nbColors, 'change', this._onNbColorsChange, this);

		this._minColor = minColor;
		this._colorRange1 = colorRange1;
		this._colorRange2 = colorRange2;
		this._maxColor = maxColor;
		this._nbColors = nbColors;
	},
	_addWeightsSection: function(container) {
		let section = L.DomUtil.create('div', 'options-section', container);

		let label = L.DomUtil.create('h4', '', section);
		label.innerHTML = 'Weights';

		let weightContainer = L.DomUtil.create('div', 'doubleRangeContainerr', section);

		let minWeight = L.DomUtil.create('input', 'value', weightContainer);
		minWeight.type = 'number';
		minWeight.placeholder = 'min';
		minWeight.min = this.gradient.minWeight;
		minWeight.max = this.gradient.maxWeight;
		minWeight.step = this.gradient.weightStep;
		minWeight.value = this.gradient.weight1;
		L.DomEvent.on(minWeight, 'change', this._onMinWeightChange, this);

		let doubleRange = L.DomUtil.create('div', 'doubleRange', weightContainer);

		let weightRange1 = L.DomUtil.create('input', 'slider', doubleRange);
		weightRange1.type = 'range';
		weightRange1.min = this.gradient.minWeight;
		weightRange1.max = this.gradient.maxWeight;
		weightRange1.step = this.gradient.weightStep;
		weightRange1.value = this.gradient.weight1;
		L.DomEvent.on(weightRange1, 'input', this._onWeightRangesUpdate, this);

		let weightRange2 = L.DomUtil.create('input', 'slider', doubleRange);
		weightRange2.type = 'range';
		weightRange2.min = this.gradient.minWeight;
		weightRange2.max = this.gradient.maxWeight;
		weightRange2.step = this.gradient.weightStep;
		weightRange2.value = this.gradient.weight2;
		L.DomEvent.on(weightRange2, 'input', this._onWeightRangesUpdate, this);

		let maxWeight = L.DomUtil.create('input', 'value', weightContainer);
		maxWeight.type = 'number';
		maxWeight.placeholder = 'max';
		maxWeight.min = this.gradient.minWeight;
		maxWeight.max = this.gradient.maxWeight;
		maxWeight.step = this.gradient.weightStep;
		maxWeight.value = this.gradient.weight2;
		L.DomEvent.on(maxWeight, 'change', this._onMaxWeightChange, this);

		this._minWeight = minWeight;
		this._weightRange1 = weightRange1;
		this._weightRange2 = weightRange2;
		this._maxWeight = maxWeight;
	},
	_onGradientChange: function() {
		this.gradient.minColor = parseInt(this._minColor.value);
		this.gradient.maxColor = parseInt(this._maxColor.value);
		this.gradient.nbColors = parseInt(this._nbColors.value);
		this.gradient.weight1 = parseFloat(this._minWeight.value);
		this.gradient.weight2 = parseFloat(this._maxWeight.value);

		this._maptt.setEdgeColorScheme();
	},

	//---------------------------------------------------------------------------
	// * Color scheme
	//---------------------------------------------------------------------------
	_onColorRangeUpdate: function() {
		const v1 = parseInt(this._colorRange1.value),
		      v2 = parseInt(this._colorRange2.value),
		      min = Math.min(v1, v2),
		      max = Math.max(v1, v2);

		this._minColor.value = min;
		this._maxColor.value = max;

		this._onGradientChange();
	},
	_updateColorRanges: function() {
		const min = parseInt(this._minColor.value),
		      max = parseInt(this._maxColor.value);

		this._colorRange1.value = min;
		this._colorRange2.value = max;

		this._onGradientChange();
	},
	_onMinColorChange: function() {
		const minV = parseInt(this._minColor.value),
		      maxV = parseInt(this._maxColor.value),
		      min = parseInt(this._minColor.min);

		// Prevent forbidden values
		if (minV > maxV) {
			this._minColor.value = maxV;
		} else if (isNaN(minV) || minV < min) {
			this._minColor.value = min;
		}

		this._updateColorRanges();
	},
	_onMaxColorChange: function() {
		const minV = parseInt(this._minColor.value),
		      maxV = parseInt(this._maxColor.value),
		      max = parseInt(this._maxColor.max);

		// Prevent forbidden values
		if (maxV < minV) {
			this._maxColor.value = minV;
		} else if (isNaN(maxV) || maxV > max) {
			this._maxColor.value = max;
		}

		this._updateColorRanges();
	},
	_onNbColorsChange: function() {
		const max = parseInt(this._nbColors.value),
		      minV = parseInt(this._minColor.value),
		      maxV = parseInt(this._maxColor.value);

		if (minV > max) {
			this._minColor.value = max;
		}
		if (maxV > max) {
			this._maxColor.value = max;
		}

		this._colorRange1.max = this._colorRange2.max = this._minColor.max = this._maxColor.max = max;

		this._gradient.max = max;
		this._onGradientChange();
	},

	//---------------------------------------------------------------------------
	// * Weights
	//---------------------------------------------------------------------------
	_onWeightRangesUpdate: function() {
		const v1 = parseFloat(this._weightRange1.value),
		      v2 = parseFloat(this._weightRange2.value),
		      min = Math.min(v1,v2),
		      max = Math.max(v1,v2);

		this._minWeight.value = min;
		this._maxWeight.value = max;

		this._onGradientChange();
	},
	_updateWeightRanges: function() {
		const min = parseFloat(this._minWeight.value),
		      max = parseFloat(this._maxWeight.value);

		this._weightRange1.value = min;
		this._weightRange2.value = max;

		this._minColor.title = `Indice of the first color to allocate to the weights between ${min} and ${max}`;
		this._maxColor.title = `Indice of the last color to allocate to the weights between ${min} and ${max}`;

		this._onGradientChange();
	},
	_onMinWeightChange: function() {
		const minV = parseFloat(this._minWeight.value),
		      maxV = parseFloat(this._maxWeight.value),
		      min = parseFloat(this._minWeight.min);

		// Prevent forbidden values
		if (minV > maxV) {
			this._minWeight.value = maxV;
		} else if (isNaN(minV) || minV < min) {
			this._minWeight.value = min;
		}

		this._updateWeightRanges();
	},
	_onMaxWeightChange: function() {
		const minV = parseFloat(this._minWeight.value),
		      maxV = parseFloat(this._maxWeight.value),
		      max = parseFloat(this._maxWeight.max);

		// Prevent forbidden values
		if (maxV < minV) {
			this._maxWeight.value = minV;
		} else if (isNaN(maxV) || maxV > max) {
			this._maxWeight.value = max;
		}

		this._updateWeightRanges();
	},
});

L.control.maptt.gradient = function(options) {
	return new L.Control.MapTT.Gradient(options);
};
