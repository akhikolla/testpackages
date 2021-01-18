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

L.MyTimeDimension = L.TimeDimension.extend({
	// --------------------------------------------------------------------------
	// New function
	// Fires 'timedivert', meaning that currentTimeIndex variable is not
	// changing naturally (forward/backward), but can jump to different values,
	// like when the user is using the time bar.
	// --------------------------------------------------------------------------
	updateCurrentTimeIndex: function(newIndex) {
		this.fire('timedivert', {
			newIndex: newIndex
		});
		this.setCurrentTimeIndex(newIndex);
	},

	// --------------------------------------------------------------------------
	// Edit function
	// Fires 'loop' when looping
	// --------------------------------------------------------------------------
	nextTime: function (numSteps, loop) {
		if (!numSteps) {
			numSteps = 1;
		}
		var newIndex = this._currentTimeIndex;
		var upperLimit = this._upperLimit || this._availableTimes.length - 1;
		var lowerLimit = this._lowerLimit || 0;
		if (this._loadingTimeIndex > -1) {
			newIndex = this._loadingTimeIndex;
		}
		newIndex = newIndex + numSteps;
		if (newIndex > upperLimit) {
			if (!!loop) {
				newIndex = lowerLimit;
				this.fire('loop', {
					lowerLimit: true,
					newIndex: newIndex
				});
			} else {
				newIndex = upperLimit;
			}
		}
		// loop backwards
		if (newIndex < lowerLimit) {
			if (!!loop) {
				newIndex = upperLimit;
				this.fire('loop', {
					lowerLimit: false,
					newIndex: newIndex
				});
			} else {
				newIndex = lowerLimit;
			}
		}
		this.setCurrentTimeIndex(newIndex);
	},
});
