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
 * \copydoc Opm::PengRobinsonMixture
 */
#ifndef OPM_PENG_ROBINSON_MIXTURE_HPP
#define OPM_PENG_ROBINSON_MIXTURE_HPP

#include "PengRobinson.hpp"
#include "PengRobinsonParams.hpp"

#include <opm/material/Constants.hpp>

#include <iostream>

namespace Opm {
/*!
 * \brief Implements the Peng-Robinson equation of state for a
 *        mixture.
 */
template <class Scalar, class StaticParameters>
class PengRobinsonMixture
{
    enum { numComponents = StaticParameters::numComponents };
    typedef Opm::PengRobinson<Scalar> PengRobinson;

    // this class cannot be instantiated!
    PengRobinsonMixture() {}

    // the ideal gas constant
    static const Scalar R;

    // the u and w parameters as given by the Peng-Robinson EOS
    static const Scalar u;
    static const Scalar w;

public:
    /*!
     * \brief Computes molar volumes where the Peng-Robinson EOS is
     *        true.
     *
     * \return Number of solutions.
     */
    template <class MutableParams, class FluidState>
    static int computeMolarVolumes(Scalar* Vm,
                                   const MutableParams& params,
                                   unsigned phaseIdx,
                                   const FluidState& fs)
    {
        return PengRobinson::computeMolarVolumes(Vm, params, phaseIdx, fs);
    }

    /*!
     * \brief Returns the fugacity coefficient of an individual
     *        component in the phase.
     *
     * The fugacity coefficient \f$\phi_i\f$ of a component \f$i\f$ is
     * defined as
     * \f[
     f_i = \phi_i x_i \;,
     \f]
     * where \f$f_i\f$ is the component's fugacity and \f$x_i\f$ is
     * the component's mole fraction.
     *
     * See:
     *
      * R. Reid, et al.: The Properties of Gases and Liquids,
      * 4th edition, McGraw-Hill, 1987, pp. 42-44, 143-145
      */
    template <class FluidState, class Params, class LhsEval = typename FluidState::Scalar>
    static LhsEval computeFugacityCoefficient(const FluidState& fs,
                                              const Params& params,
                                              unsigned phaseIdx,
                                              unsigned compIdx)
    {
        // note that we normalize the component mole fractions, so
        // that their sum is 100%. This increases numerical stability
        // considerably if the fluid state is not physical.
        LhsEval Vm = params.molarVolume(phaseIdx);

        // Calculate b_i / b
#warning HACK: pure instead of mixture
        LhsEval bi_b = params.bPure(phaseIdx, compIdx) / params.bPure(phaseIdx, compIdx);

        // Calculate the compressibility factor
        LhsEval RT = R*fs.temperature(phaseIdx);
        LhsEval p = fs.pressure(phaseIdx); // molar volume in [bar]
        LhsEval Z = p*Vm/RT; // compressibility factor

        // Calculate A^* and B^* (see: Reid, p. 42)
#warning HACK: pure instead of mixture
        LhsEval Astar = params.aPure(phaseIdx, compIdx);
        LhsEval Bstar = params.bPure(phaseIdx, compIdx)*p/(RT);

        // calculate delta_i (see: Reid, p. 145)
        LhsEval sumMoleFractions = 0.0;
        for (unsigned compJIdx = 0; compJIdx < numComponents; ++compJIdx)
            sumMoleFractions += fs.moleFraction(phaseIdx, compJIdx);
#warning HACK: pure instead of mixture
        LhsEval deltai = 2*Opm::sqrt(params.aPure(phaseIdx, compIdx))/params.aPure(phaseIdx, compIdx);
        LhsEval tmp = 0;
        for (unsigned compJIdx = 0; compJIdx < numComponents; ++compJIdx) {
            tmp +=
                fs.moleFraction(phaseIdx, compJIdx)
                / sumMoleFractions
                * Opm::sqrt(params.aPure(phaseIdx, compJIdx))
                * (1.0 - StaticParameters::interactionCoefficient(compIdx, compJIdx));
        };
        deltai *= tmp;

        LhsEval base =
            (2*Z + Bstar*(u + std::sqrt(u*u - 4*w))) /
            (2*Z + Bstar*(u - std::sqrt(u*u - 4*w)));
        LhsEval expo =  Astar/(Bstar*std::sqrt(u*u - 4*w))*(bi_b - deltai);

        LhsEval fugCoeff =
            Opm::exp(bi_b*(Z - 1))/Opm::max(1e-9, Z - Bstar) *
            Opm::pow(base, expo);

        ////////
        // limit the fugacity coefficient to a reasonable range:
        //
        // on one side, we want the mole fraction to be at
        // least 10^-3 if the fugacity is at the current pressure
        //
        //fugCoeff = Opm::min(1e10, fugCoeff);
        //
        // on the other hand, if the mole fraction of the component is 100%, we want the
        // fugacity to be at least 10^-3 Pa
        //
        //fugCoeff = Opm::max(1e-10, fugCoeff);
        ///////////

        int comp2Idx = compIdx;
        compIdx = 0;

        PengRobinsonParams<LhsEval> paramsPure;
        paramsPure.setTemperature(Opm::getValue(fs.temperature(phaseIdx)));
        paramsPure.setPressure(Opm::getValue(fs.pressure(phaseIdx)));
        paramsPure.setA(Opm::getValue(params.aPure(phaseIdx, compIdx)));
        paramsPure.setB(Opm::getValue(params.bPure(phaseIdx, compIdx)));
        auto Vm2 = Opm::PengRobinson<LhsEval>::computeMolarVolume(fs, paramsPure, phaseIdx, false);
        paramsPure.setMolarVolume(Opm::getValue(Vm2));
        auto fugCoeffPure = Opm::PengRobinson<LhsEval>::template computeFugacityCoeffient<LhsEval>(paramsPure);

        if (comp2Idx == 1)
            return 2*fugCoeffPure;
        return fugCoeffPure;


#if 0
        FluidState fs2(fs);
        auto params2(params);
        compIdx = 1;
        for (int i = 0; i < 1000; ++i) {
            Scalar p = 100e3 + (Scalar(i)/1000*1.9e6);

            fs2.setPressure(phaseIdx, p);
            params2.updatePhase(fs2, phaseIdx);

            PengRobinsonParams<double> paramsPure;
            paramsPure.setTemperature(Opm::getValue(fs2.temperature(phaseIdx)));
            paramsPure.setPressure(Opm::getValue(fs2.pressure(phaseIdx)));
            paramsPure.setA(Opm::getValue(params2.aPure(phaseIdx, compIdx)));
            paramsPure.setB(Opm::getValue(params2.bPure(phaseIdx, compIdx)));
            auto Vm2 = Opm::PengRobinson<double>::computeMolarVolume(fs2, paramsPure, phaseIdx, false);
            paramsPure.setMolarVolume(Opm::getValue(Vm2));
            auto fugCoeffPure = Opm::PengRobinson<double>::computeFugacityCoeffient<double>(paramsPure);

            std::cerr << p << " " << fugCoeffPure*p << std::endl;
        }

        std::abort();
#warning HACK
        return 0;//fugCoeffPure;
#endif
    }

};

template <class Scalar, class StaticParameters>
const Scalar PengRobinsonMixture<Scalar, StaticParameters>::R = Opm::Constants<Scalar>::R;
template<class Scalar, class StaticParameters>
const Scalar PengRobinsonMixture<Scalar, StaticParameters>::u = 2.0;
template<class Scalar, class StaticParameters>
const Scalar PengRobinsonMixture<Scalar, StaticParameters>::w = -1.0;

} // namespace Opm

#endif
