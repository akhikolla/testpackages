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

L.Control.MyTimeDimension = L.Control.TimeDimension.extend({
	// --------------------------------------------------------------------------
	// Edit function
	// Doesn't call setCurrentTimeIndex function directly.
	// It means that user is using the time bar, so we now that we need to
	// update display accordingly.
	// --------------------------------------------------------------------------
	_sliderTimeValueChanged: function(newIndex) {
		this._timeDimension.updateCurrentTimeIndex(newIndex);
	},
	_buttonBackwardClicked: function(...args) {
		console.log("------------ Backward ------------");
		L.Control.TimeDimension.prototype._buttonBackwardClicked.call(this, args);
	},
	_buttonForwardClicked: function(...args) {
		console.log("------------ Forward -------------");
		L.Control.TimeDimension.prototype._buttonForwardClicked.call(this, args);
	},
});
