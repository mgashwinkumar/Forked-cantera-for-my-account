/**
 *  @file SpeciesThermoInterpType.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/SpeciesThermoInterpType.h"
#include "cantera/thermo/PDSS.h"

namespace Cantera
{

SpeciesThermoInterpType::SpeciesThermoInterpType() :
    m_lowT(0.0),
    m_highT(0.0),
    m_Pref(0.0)
{
}

SpeciesThermoInterpType::SpeciesThermoInterpType(double tlow,
                                                 double thigh,
                                                 double pref) :
    m_lowT(tlow),
    m_highT(thigh),
    m_Pref(pref)
{
}

void SpeciesThermoInterpType::updateProperties(const doublereal* tempPoly,
        doublereal* cp_R, doublereal* h_RT, doublereal* s_R) const
{
    double T = tempPoly[0];
    updatePropertiesTemp(T, cp_R, h_RT, s_R);
}

void SpeciesThermoInterpType::updatePropertiesTemp(const double temp,
        double* cp_R, double* h_RT, double* s_R) const
{
    throw CanteraError("SpeciesThermoInterpType::updatePropertiesTemp",
                       "Not implemented");
}

void SpeciesThermoInterpType::reportParameters(size_t& index, int& type,
        double& minTemp, double& maxTemp, double& refPressure,
        double* const coeffs) const
{
    throw CanteraError("SpeciesThermoInterpType::reportParameters",
                       "Not implemented");
}

doublereal SpeciesThermoInterpType::reportHf298(doublereal* const h298) const
{
    throw CanteraError("SpeciesThermoInterpType::reportHf298",
                       "Not implemented");
}

void SpeciesThermoInterpType::modifyOneHf298(const size_t k,
                                             const doublereal Hf298New)
{
    throw CanteraError("SpeciesThermoInterpType::modifyOneHf298",
                       "Not implemented");
}

}
