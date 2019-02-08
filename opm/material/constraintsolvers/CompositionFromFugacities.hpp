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
 * \copydoc Opm::CompositionFromFugacities
 */
#ifndef OPM_COMPOSITION_FROM_FUGACITIES_HPP
#define OPM_COMPOSITION_FROM_FUGACITIES_HPP

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/common/Exceptions.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <limits>

namespace Opm {

/*!
 * \brief Calculates the chemical equilibrium from the component
 *        fugacities in a phase.
 */
template <class Scalar, class FluidSystem, class Evaluation = Scalar>
class CompositionFromFugacities
{
    enum { numComponents = FluidSystem::numComponents };

public:
    typedef Dune::FieldVector<Evaluation, numComponents> ComponentVector;

    /*!
     * \brief Guess an initial value for the composition of the phase.
     */
    template <class FluidState>
    static void guessInitial(FluidState& fluidState,
                             unsigned phaseIdx,
                             const ComponentVector& /*fugVec*/)
    {
        if (FluidSystem::isIdealMixture(phaseIdx))
            return;

        // Pure component fugacities
        if (phaseIdx == FluidSystem::oilPhaseIdx) {
            //std::cout << f << " -> " << mutParams.fugacity(phaseIdx, i)/f << "\n";
            fluidState.setMoleFraction(phaseIdx,
                                       FluidSystem::OctaneIdx,
                                       0.02);
            fluidState.setMoleFraction(phaseIdx,
                                       FluidSystem::BrineIdx,
                                       0.25);
            fluidState.setMoleFraction(phaseIdx,
                                       FluidSystem::CO2Idx,
                                       0.02);
        }
        else
            assert(false);
    }

    /*!
     * \brief Calculates the chemical equilibrium from the component
     *        fugacities in a phase.
     *
     * The phase's fugacities must already be set.
     */
    template <class FluidState>
    static void solve(FluidState& fluidState,
                      typename FluidSystem::template ParameterCache<typename FluidState::Scalar>& paramCache,
                      unsigned phaseIdx,
                      const ComponentVector& targetFug)
    {
        // use a much more efficient method in case the phase is an
        // ideal mixture
        if (FluidSystem::isIdealMixture(phaseIdx)) {
            solveIdealMix_(fluidState, paramCache, phaseIdx, targetFug);
            return;
        }



#if 1
        // save initial composition in case something goes wrong
        Dune::FieldVector<Evaluation, numComponents> xInit;
        for (unsigned i = 0; i < numComponents; ++i) {
            xInit[i] = fluidState.moleFraction(phaseIdx, i);
        }

        /////////////////////////
        // Newton method
        /////////////////////////

        // Jacobian matrix
        Dune::FieldMatrix<Evaluation, numComponents, numComponents> J;
        // solution, i.e. phase composition
        Dune::FieldVector<Evaluation, numComponents> x;
        // right hand side
        Dune::FieldVector<Evaluation, numComponents> b;

        paramCache.updatePhase(fluidState, phaseIdx);

        // maximum number of iterations
        const int nMax = 151;
        for (int nIdx = 0; nIdx < nMax; ++nIdx) {
            paramCache.updateAll(fluidState);

            const Evaluation& rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
            fluidState.setDensity(phaseIdx, rho);

            // calculate Jacobian matrix and right hand side
            linearize_(J, b, fluidState, paramCache, phaseIdx, targetFug);
            Valgrind::CheckDefined(J);
            Valgrind::CheckDefined(b);

#if 0
            std::cout << FluidSystem::phaseName(phaseIdx) << "Phase composition: ";
            for (unsigned i = 0; i < FluidSystem::numComponents; ++i)
                std::cout << fluidState.moleFraction(phaseIdx, i) << " ";
            std::cout << "\n";
            std::cout << FluidSystem::phaseName(phaseIdx) << "Phase phi: ";
            for (unsigned i = 0; i < FluidSystem::numComponents; ++i)
                std::cout << fluidState.fugacityCoefficient(phaseIdx, i) << " ";
            std::cout << "\n";
#endif

            // Solve J*x = b
            x = 0.0;
            try { J.solve(x, b); }
            catch (const Dune::FMatrixError& e)
            { throw Opm::NumericalIssue(e.what()); }

            //std::cout << "original delta: " << x << "\n";

            Valgrind::CheckDefined(x);

#if 0
            std::cout << FluidSystem::phaseName(phaseIdx) << "Phase composition: ";
            for (unsigned i = 0; i < FluidSystem::numComponents; ++i)
                std::cout << fluidState.massFraction(phaseIdx, i) << " ";
            std::cout << "\n";
            std::cout << "J: " << J << "\n";
            std::cout << "rho: " << fluidState.density(phaseIdx) << "\n";
            std::cout << "delta: " << x << "\n";
            std::cout << "defect: " << b << "\n";

            std::cout << "J: " << J << "\n";

            std::cout << "---------------------------\n";
#endif

            // update the fluid composition. b is also used to store
            // the defect for the next iteration.
            Scalar relError = update_(fluidState, paramCache, x, b, phaseIdx, targetFug);

            if (relError < 1e-10) {
                paramCache.updatePhase(fluidState, phaseIdx);
                const Evaluation& rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
                fluidState.setDensity(phaseIdx, rho);

                for (int i = 0; i < numComponents; ++i) {
                    const Evaluation& phi = FluidSystem::fugacityCoefficient(fluidState, paramCache, phaseIdx, i);
                    fluidState.setFugacityCoefficient(phaseIdx, i, phi);
                }

                //std::abort();

                //std::cout << "num iterations: " << nIdx << "\n";
                return;
            }
        }

        std::ostringstream oss;
        oss << "Calculating the " << FluidSystem::phaseName(phaseIdx)
            << "Phase composition failed. Initial {x} = {"
            << xInit
            << "}, {fug_t} = {" << targetFug << "}, p = " << fluidState.pressure(phaseIdx)
            << ", T = " << fluidState.temperature(phaseIdx);
        throw Opm::NumericalIssue(oss.str());
#endif
    }


protected:
    // update the phase composition in case the phase is an ideal
    // mixture, i.e. the component's fugacity coefficients are
    // independent of the phase's composition.
    template <class FluidState>
    static void solveIdealMix_(FluidState& fluidState,
                               typename FluidSystem::template ParameterCache<typename FluidState::Scalar>& paramCache,
                               unsigned phaseIdx,
                               const ComponentVector& fugacities)
    {
        paramCache.updatePhase(fluidState, phaseIdx);
        for (unsigned i = 0; i < numComponents; ++ i) {
            const Evaluation& phi = FluidSystem::fugacityCoefficient(fluidState,
                                                                     paramCache,
                                                                     phaseIdx,
                                                                     i);
            const Evaluation& gamma = phi * fluidState.pressure(phaseIdx);
            Valgrind::CheckDefined(phi);
            Valgrind::CheckDefined(gamma);
            Valgrind::CheckDefined(fugacities[i]);
            fluidState.setFugacityCoefficient(phaseIdx, i, phi);
            fluidState.setMoleFraction(phaseIdx, i, fugacities[i]/gamma);
        };

        paramCache.updatePhase(fluidState, phaseIdx);

        const Evaluation& rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
        fluidState.setDensity(phaseIdx, rho);
        return;
    }

#warning HACK
#if 1
    template <class FluidState>
    static void linearize_(Dune::FieldMatrix<Evaluation, numComponents, numComponents>& J,
                           Dune::FieldVector<Evaluation, numComponents>& defect,
                           const FluidState& fluidStateIn,
                           typename FluidSystem::template ParameterCache<typename FluidState::Scalar>& paramCache,
                           unsigned phaseIdx,
                           const ComponentVector& targetFug)
    {
        // reset jacobian
        J = 0;

        FluidState fluidState(fluidStateIn);

        // normalize composition of input fluid state
#if 0
        Evaluation sumx = 0.0;
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            sumx += fluidState.moleFraction(phaseIdx, compIdx);
        sumx = Opm::max(sumx, 1e-7);

        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            Evaluation  xi = fluidState.moleFraction(phaseIdx, compIdx)/sumx;
            fluidState.setMoleFraction(phaseIdx, compIdx, xi);
        }
#endif

        // assemble jacobian matrix of the constraints for the composition
        typedef Opm::DenseAd::Evaluation<Evaluation, numComponents> InnerEval;
        typename FluidSystem::template ParameterCache<InnerEval> myParamCache;
        typedef Opm::CompositionalFluidState<InnerEval, FluidSystem> CFS;
        CFS myFluidState;
        myFluidState.setPressure(phaseIdx, fluidState.pressure(phaseIdx));
        myFluidState.setTemperature(fluidState.temperature(phaseIdx));
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            myFluidState.setMoleFraction(phaseIdx, compIdx, fluidState.moleFraction(phaseIdx, compIdx));
        }
        myParamCache.updatePhase(myFluidState, phaseIdx);
        for (unsigned i = 0; i < numComponents; ++ i) {
            // deviate the mole fraction of the i-th component
            myFluidState.setMoleFraction(phaseIdx, i, Opm::variable<InnerEval>(fluidState.moleFraction(phaseIdx, i), i));
            myParamCache.updateAll(myFluidState);
            //myFluidState.setDensity(phaseIdx, FluidSystem::density(myFluidState, myParamCache, phaseIdx));

            // compute new defect and derivative for all component
            // fugacities
            for (unsigned j = 0; j < numComponents; ++j) {
                // compute the j-th component's fugacity coefficient ...
                const auto& phi = FluidSystem::fugacityCoefficient(myFluidState, myParamCache, phaseIdx, j);
                // ... and its fugacity ...
                const auto& f =
                    phi *
                    myFluidState.pressure(phaseIdx) *
                    myFluidState.moleFraction(phaseIdx, j);

#if 0
#warning HACK
                if (j == i)
                    assert(std::abs(f.derivative(i) - phi.value()*fluidState.pressure(phaseIdx)) < 1e-2);
#endif

                // use forward differences to calculate the defect's
                // derivative
                J[j][i] = - f.derivative(j);
            }

            myFluidState.setMoleFraction(phaseIdx, i, fluidState.moleFraction(phaseIdx, i));

            // compute the j-th component's fugacity coefficient ...
            const auto& phi = FluidSystem::fugacityCoefficient(myFluidState, myParamCache, phaseIdx, i);
            // ... and its fugacity ...
            const auto& f =
                phi *
                myFluidState.pressure(phaseIdx) *
                myFluidState.moleFraction(phaseIdx, i);

            defect[i] = targetFug[i] - f.value();


            // end forward differences
            ////////

        }
    }

    template <class FluidState>
    static Scalar update_(FluidState& fluidState,
                          typename FluidSystem::template ParameterCache<typename FluidState::Scalar>& paramCache,
                          Dune::FieldVector<Evaluation, numComponents>& x,
                          Dune::FieldVector<Evaluation, numComponents>& /*b*/,
                          unsigned phaseIdx,
                          const Dune::FieldVector<Evaluation, numComponents>& targetFug)
    {
        // change composition
        Scalar relError = 0.0;
        for (unsigned i = 0; i < numComponents; ++i) {
            Evaluation newx = fluidState.moleFraction(phaseIdx, i) - x[i];
            fluidState.setMoleFraction(phaseIdx, i, newx);
     //       if (targetFug[i] == 0)
     //           std::cout << targetFug[i] << std::endl;

            relError = std::max(relError, std::abs(Opm::getValue(x[i])));
        }

        paramCache.updateComposition(fluidState, phaseIdx);

        return relError;
    }

    template <class FluidState>
    static Scalar calculateDefect_(const FluidState& params,
                                   unsigned phaseIdx,
                                   const ComponentVector& targetFug)
    {
        Scalar result = 0.0;
        for (unsigned i = 0; i < numComponents; ++i) {
            // sum of the fugacity defect weighted by the inverse
            // fugacity coefficient
            result += std::abs(
                (targetFug[i] - params.fugacity(phaseIdx, i))
                /
                params.fugacityCoefficient(phaseIdx, i) );
        };
        return result;
    }
#endif
}; // namespace Opm

} // end namespace Opm

#endif
