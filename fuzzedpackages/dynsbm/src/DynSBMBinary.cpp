/*
    This file is part of dynsbm.

    dysbm is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    dynsbm is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with dynsbm.  If not, see <http://www.gnu.org/licenses/>
 */
#include<DynSBMBinary.h>
namespace dynsbm{
  void DynSBMBinary::updateTheta(int*** const Y){// M-step
    DynSBMBinaryAddEventFunctor addEventFunctor(*this);
    updateThetaCore<DynSBMBinaryAddEventFunctor>(Y, addEventFunctor);
  }
  void DynSBMBinary::updateFrozenTheta(int*** const Y){// M-step
    DynSBMBinaryAddEventFunctor addEventFunctor(*this);
    updateFrozenThetaCore<DynSBMBinaryAddEventFunctor>(Y, addEventFunctor);
  }
}
