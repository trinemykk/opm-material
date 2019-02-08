// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Opm::PengRobinsonParams
 */
#ifndef OPM_PENG_ROBINSON_PARAMS_HPP
#define OPM_PENG_ROBINSON_PARAMS_HPP

#include <opm/material/common/Valgrind.hpp>

namespace Opm
{
/*!
 * \brief Stores and provides access to the Peng-Robinson parameters
 *
 * See:
 *
 * R. Reid, et al.: The Properties of Gases and Liquids, 4th edition,
 * McGraw-Hill, 1987, pp. 43-44
 */
template <class Scalar>
class PengRobinsonParams
{
public:
    /*!
     * \brief Returns the attractive parameter 'a' of the
     *        Peng-Robinson fluid.
     */
    Scalar a(int phaseIdx=0) const
    { return a_; }

    /*!
     * \brief Returns the repulsive parameter 'b' of the Peng-Robinson
     *        fluid.
     */
    Scalar b(int phaseIdx=0) const
    { return b_; }

    /*!
     * \brief If run under valgrind, this method produces an warning
     *        if the parameters where not determined correctly.
     */
    void checkDefined() const
    {
#ifndef NDEBUG
        Valgrind::CheckDefined(a_);
        Valgrind::CheckDefined(b_);
#endif
    }

    /*!
     * \brief Set the attractive parameter 'a' of the Peng-Robinson
     *        fluid.
     */
    void setA(Scalar value)
    { a_ = value; }

    /*!
     * \brief Set the repulsive parameter 'b' of the Peng-Robinson
     *        fluid.
     */
    void setB(Scalar value)
    { b_ = value; }

    Scalar temperature(int phaseIdx = 0) const
    { return temperature_; };

    Scalar pressure(int phaseIdx = 0) const
    { return pressure_; };

    Scalar molarVolume(int phaseIdx = 0) const
    { return molarVolume_; };

    void setTemperature(Scalar v)
    { temperature_ = v; };

    void setPressure(Scalar v)
    { pressure_ = v; };

    void setMolarVolume(Scalar v)
    { molarVolume_ = v; };

protected:
    Scalar a_;
    Scalar b_;
    Scalar temperature_;
    Scalar pressure_;
    Scalar molarVolume_;
};

} // namespace Opm

#endif
