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
    static LhsEval computeFugacityCoefficient(const FluidState& fs2,
                                              const Params& params2,
                                              unsigned phaseIdx,
                                              unsigned compIdx)
    {
#if 0
#warning HACK
        auto fs = fs2;
        auto params = params2;
        fs.setMoleFraction(0, 0, 1.0);
        fs.setMoleFraction(0, 1, 0.0);
        fs.setMoleFraction(0, 2, 0.0);
        params.updatePhase(fs, 0);
#else
#warning also a HACK
        auto fs = fs2;
        double sumx = 0.0;
        for (int i = 0; i < FluidState::numComponents; ++i)
            sumx += Opm::scalarValue(fs2.moleFraction(phaseIdx, i));
        if (sumx < 0.95)  {
            double alpha = 0.95/sumx;
            std::cerr << "normalize: " << sumx
                      << " alpha: " << alpha << "\n";
            for (int i = 0; i < FluidState::numComponents; ++i)
                fs.setMoleFraction(phaseIdx, i, alpha*fs2.moleFraction(phaseIdx, i));
        }

        /*
#warning more of a HACK
        fs.setMoleFraction(phaseIdx, 0, 0.0); // octane
        fs.setMoleFraction(phaseIdx, 1, 0.0); // co2
        fs.setMoleFraction(phaseIdx, 2, 0.95); // water
        */

        auto params = params2;
        params.updatePhase(fs, phaseIdx);
#endif
     
        // note that we normalize the component mole fractions, so
        // that their sum is 100%. This increases numerical stability
        // considerably if the fluid state is not physical.
        LhsEval Vm = params.molarVolume(phaseIdx);

        // Calculate b_i / b
        LhsEval bi_b = params.bPure(phaseIdx, compIdx) / params.b(phaseIdx);

        // Calculate the compressibility factor
        LhsEval RT = R*fs.temperature(phaseIdx);
        LhsEval p = fs.pressure(phaseIdx); // molar volume in [bar]
        LhsEval Z = p*Vm/RT; // compressibility factor

        // Calculate A^* and B^* (see: Reid, p. 42)
        LhsEval Astar = params.a(phaseIdx)*p/(RT*RT);
        LhsEval Bstar = params.b(phaseIdx)*p/(RT);

        // calculate delta_i (see: Reid, p. 145)
        LhsEval sumMoleFractions = 0.0;
        for (unsigned compJIdx = 0; compJIdx < numComponents; ++compJIdx)
            sumMoleFractions += fs.moleFraction(phaseIdx, compJIdx);
        LhsEval deltai = 2*Opm::sqrt(params.aPure(phaseIdx, compIdx))/params.a(phaseIdx);
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

        std::cout << "Z - Bstar:" << Z - Bstar << "\n";
#if 1
        LhsEval fugCoeff =
            Opm::exp(bi_b*(Z - 1))/Opm::max(1e-9, Z - Bstar) *
            Opm::pow(base, expo);
#else
#warning HACK
        LhsEval fugCoeff =
            Opm::exp(bi_b*(Z - 1))/Opm::max(1e-9, Z - Bstar);
#endif

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

        auto fugCoeffMix = fugCoeff;
        
#if 1 // -> use fugacity coefficient of mixture
        return fugCoeffMix;
#else // -> use pure fugacity coefficient
        PengRobinsonParams<LhsEval> paramsPure;
        paramsPure.setTemperature(Opm::getValue(fs.temperature(phaseIdx)));
        paramsPure.setPressure(Opm::getValue(fs.pressure(phaseIdx)));
        paramsPure.setA(Opm::getValue(params.aPure(phaseIdx, compIdx)));
        paramsPure.setB(Opm::getValue(params.bPure(phaseIdx, compIdx)));
        auto Vm2 = Opm::PengRobinson<LhsEval>::computeMolarVolume(fs, paramsPure, phaseIdx, false);
        paramsPure.setMolarVolume(Opm::getValue(Vm2));
        auto fugCoeffPure = Opm::PengRobinson<LhsEval>::template computeFugacityCoeffient<LhsEval>(paramsPure);

        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            PengRobinsonParams<LhsEval> paramsPure;
            paramsPure.setTemperature(Opm::getValue(fs.temperature(phaseIdx)));
            paramsPure.setPressure(Opm::getValue(fs.pressure(phaseIdx)));
            paramsPure.setA(Opm::getValue(params.aPure(phaseIdx, compIdx)));
            paramsPure.setB(Opm::getValue(params.bPure(phaseIdx, compIdx)));
            auto Vm2 = Opm::PengRobinson<LhsEval>::computeMolarVolume(fs, paramsPure, phaseIdx, false);
            paramsPure.setMolarVolume(Opm::getValue(Vm2));
            auto fugCoeffPure = Opm::PengRobinson<LhsEval>::template computeFugacityCoeffient<LhsEval>(paramsPure);

            std::cout << "fugCoeffPure_" << phaseIdx
                      << "^" << compIdx << " = "
                      << fugCoeffPure << "\n";
            
        }

        std::cout << fugCoeffMix << " vs: " << fugCoeffPure << "\n";
        return fugCoeffPure;
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
