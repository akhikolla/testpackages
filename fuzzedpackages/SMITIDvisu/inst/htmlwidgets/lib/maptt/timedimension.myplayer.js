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

/*jshint indent: 4, browser:true*/
/*global L*/


/*
 * L.TimeDimension.Player
 */
//'use strict';
L.TimeDimension.MyPlayer = (L.TimeDimension.Player).extend({
    start: function(numSteps) {
        if (this._intervalID) return;
        this._steps = numSteps || 1;
        this._waitingForBuffer = false;
        if (this.options.startOver){
            if (this._timeDimension.getCurrentTimeIndex() === this._getMaxIndex()){
                 this._timeDimension.setCurrentTimeIndex(this._timeDimension.getLowerLimitIndex() || 0);
            }
        }
        this.release();
        this._intervalID = true;
        this._tick();
        this.fire('play');
        this.fire('running');
    },

    stop: function() {
        if (!this._intervalID) return;
        this._intervalID = null;
        this._waitingForBuffer = false;
        this.fire('stop');
    },
});
