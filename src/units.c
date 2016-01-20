/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "units.h"

/* Includes. */
#include "const.h"
#include "error.h"
#include "units.h"

/**
 * @brief Initialises the UnitSystem structure with the constants given in
 * const.h
 * @param us The UnitSystem to initialize
 */

void initUnitSystem(struct UnitSystem* us) {
  us->UnitMass_in_cgs = const_unit_mass_in_cgs;
  us->UnitLength_in_cgs = const_unit_length_in_cgs;
  us->UnitTime_in_cgs = 1. / ((double)const_unit_velocity_in_cgs /
                              ((double)const_unit_length_in_cgs));
  us->UnitCurrent_in_cgs = 1.;
  us->UnitTemperature_in_cgs = 1.;
}

/**
 * @brief Returns the base unit conversion factor for a given unit system
 * @param us The UnitSystem used
 * @param baseUnit The base unit
 */
double getBaseUnit(struct UnitSystem* us, enum BaseUnits baseUnit) {
  switch (baseUnit) {
    case UNIT_MASS:
      return us->UnitMass_in_cgs;
    case UNIT_LENGTH:
      return us->UnitLength_in_cgs;
    case UNIT_TIME:
      return us->UnitTime_in_cgs;
    case UNIT_CURRENT:
      return us->UnitCurrent_in_cgs;
    case UNIT_TEMPERATURE:
      return us->UnitTemperature_in_cgs;
    default:
      error("Invalid base Unit");
  }
  return 0.0;
}

/**
 * @brief Returns the base unit symbol
 * @param baseUnit The base unit
 */
const char* getBaseUnitSymbol(enum BaseUnits baseUnit) {
  switch (baseUnit) {
    case UNIT_MASS:
      return "U_M";
    case UNIT_LENGTH:
      return "U_L";
    case UNIT_TIME:
      return "U_t";
    case UNIT_CURRENT:
      return "U_I";
    case UNIT_TEMPERATURE:
      return "U_T";
    default:
      error("Invalid base Unit");
  }
  return "";
}

/**
 * @brief Returns the base unit symbol in the cgs system
 * @param baseUnit The base unit
 */
const char* getBaseUnitCGSSymbol(enum BaseUnits baseUnit) {
  switch (baseUnit) {
    case UNIT_MASS:
      return "g";
    case UNIT_LENGTH:
      return "cm";
    case UNIT_TIME:
      return "s";
    case UNIT_CURRENT:
      return "A";
    case UNIT_TEMPERATURE:
      return "K";
    default:
      error("Invalid base Unit");
  }
  return "";
}

void getBaseUnitExponantsArray(float baseUnitsExp[5],
                               enum UnitConversionFactor unit) {
  switch (unit) {
    case UNIT_CONV_NO_UNITS:
      break;

    case UNIT_CONV_MASS:
      baseUnitsExp[UNIT_MASS] = 1.f;
      break;

    case UNIT_CONV_LENGTH:
      baseUnitsExp[UNIT_LENGTH] = 1.f;
      break;

    case UNIT_CONV_TIME:
      baseUnitsExp[UNIT_TIME] = 1.f;
      break;

    case UNIT_CONV_FREQUENCY:
      baseUnitsExp[UNIT_TIME] = -1.f;
      break;

    case UNIT_CONV_DENSITY:
      baseUnitsExp[UNIT_MASS] = 1.f;
      baseUnitsExp[UNIT_LENGTH] = -3.f;
      break;

    case UNIT_CONV_SPEED:
      baseUnitsExp[UNIT_LENGTH] = 1.f;
      baseUnitsExp[UNIT_TIME] = -1.f;
      break;

    case UNIT_CONV_ACCELERATION:
      baseUnitsExp[UNIT_LENGTH] = 1.f;
      baseUnitsExp[UNIT_TIME] = -2.f;
      break;

    case UNIT_CONV_FORCE:
      baseUnitsExp[UNIT_MASS] = 1.f;
      baseUnitsExp[UNIT_LENGTH] = 1.f;
      baseUnitsExp[UNIT_TIME] = -2.f;
      break;

    case UNIT_CONV_ENERGY:
      baseUnitsExp[UNIT_MASS] = 1.f;
      baseUnitsExp[UNIT_LENGTH] = 2.f;
      baseUnitsExp[UNIT_TIME] = -2.f;
      break;

    case UNIT_CONV_ENERGY_PER_UNIT_MASS:
      baseUnitsExp[UNIT_LENGTH] = 2.f;
      baseUnitsExp[UNIT_TIME] = -2.f;
      break;

    case UNIT_CONV_ENTROPY:
      baseUnitsExp[UNIT_MASS] = 1.f - const_hydro_gamma;
      baseUnitsExp[UNIT_LENGTH] = 3.f * const_hydro_gamma - 1.f;
      baseUnitsExp[UNIT_TIME] = -2.f;
      break;

    case UNIT_CONV_ENTROPY_PER_UNIT_MASS:
      baseUnitsExp[UNIT_MASS] = -const_hydro_gamma;
      baseUnitsExp[UNIT_LENGTH] = 3.f * const_hydro_gamma - 1.f;
      baseUnitsExp[UNIT_TIME] = -2.f;
      break;

    case UNIT_CONV_POWER:
      baseUnitsExp[UNIT_MASS] = 1.f;
      baseUnitsExp[UNIT_LENGTH] = 2.f;
      baseUnitsExp[UNIT_TIME] = -3.f;
      break;

    case UNIT_CONV_PRESSURE:
      baseUnitsExp[UNIT_MASS] = 1.f;
      baseUnitsExp[UNIT_LENGTH] = -1.f;
      baseUnitsExp[UNIT_TIME] = -2.f;
      break;

    case UNIT_CONV_ELECTRIC_CHARGE:
      baseUnitsExp[UNIT_TIME] = 1.f;
      baseUnitsExp[UNIT_CURRENT] = 1.f;
      break;

    case UNIT_CONV_ELECTRIC_VOLTAGE:
      baseUnitsExp[UNIT_MASS] = 1.f;
      baseUnitsExp[UNIT_LENGTH] = 2.f;
      baseUnitsExp[UNIT_TIME] = -3.f;
      baseUnitsExp[UNIT_CURRENT] = -1.f;
      break;

    case UNIT_CONV_ELECTRIC_CAPACITANCE:
      baseUnitsExp[UNIT_MASS] = -1.f;
      baseUnitsExp[UNIT_LENGTH] = -2.f;
      baseUnitsExp[UNIT_TIME] = 4;
      baseUnitsExp[UNIT_CURRENT] = 2.f;
      break;

    case UNIT_CONV_ELECTRIC_RESISTANCE:
      baseUnitsExp[UNIT_MASS] = 1.f;
      baseUnitsExp[UNIT_LENGTH] = 2.f;
      baseUnitsExp[UNIT_TIME] = -3.f;
      baseUnitsExp[UNIT_CURRENT] = -2.f;
      break;

    case UNIT_CONV_ELECTRIC_CONDUCTANCE:
      baseUnitsExp[UNIT_MASS] = -1.f;
      baseUnitsExp[UNIT_LENGTH] = -2.f;
      baseUnitsExp[UNIT_TIME] = 3.f;
      baseUnitsExp[UNIT_CURRENT] = 2.f;
      break;

    case UNIT_CONV_MAGNETIC_FLUX:
      baseUnitsExp[UNIT_MASS] = 1.f;
      baseUnitsExp[UNIT_LENGTH] = 2.f;
      baseUnitsExp[UNIT_TIME] = -2.f;
      baseUnitsExp[UNIT_CURRENT] = -1.f;
      break;

    case UNIT_CONV_MAGNETIC_FIELD:
      baseUnitsExp[UNIT_MASS] = 1.f;
      baseUnitsExp[UNIT_TIME] = -2.f;
      baseUnitsExp[UNIT_CURRENT] = -1.f;
      break;

    case UNIT_CONV_MAGNETIC_INDUCTANCE:
      baseUnitsExp[UNIT_MASS] = 1.f;
      baseUnitsExp[UNIT_LENGTH] = 2.f;
      baseUnitsExp[UNIT_TIME] = -2.f;
      baseUnitsExp[UNIT_CURRENT] = -2.f;
      break;

    case UNIT_CONV_TEMPERATURE:
      baseUnitsExp[UNIT_TEMPERATURE] = 1.f;
  }
}

/**
 * @brief Returns the conversion factor for a given unit in the chosen unit
 * system
 * @param us The system of units in use
 * @param unit The unit to convert
 */
double conversionFactor(struct UnitSystem* us, enum UnitConversionFactor unit) {
  float baseUnitsExp[5] = {0.f};

  getBaseUnitExponantsArray(baseUnitsExp, unit);

  return generalConversionFactor(us, baseUnitsExp);
}

/**
 * @brief Returns the h factor exponentiation for a given unit
 * @param us The system of units in use
 * @param unit The unit to convert
 */
float hFactor(struct UnitSystem* us, enum UnitConversionFactor unit) {
  float baseUnitsExp[5] = {0.f};

  getBaseUnitExponantsArray(baseUnitsExp, unit);

  return generalhFactor(us, baseUnitsExp);
}

/**
 * @brief Returns the scaling factor exponentiation for a given unit
 * @param us The system of units in use
 * @param unit The unit to convert
 */
float aFactor(struct UnitSystem* us, enum UnitConversionFactor unit) {
  float baseUnitsExp[5] = {0.f};

  getBaseUnitExponantsArray(baseUnitsExp, unit);

  return generalaFactor(us, baseUnitsExp);
}

/**
 * @brief Returns a string containing the exponents of the base units making up
 * the conversion factors
 */
void conversionString(char* buffer, struct UnitSystem* us,
                      enum UnitConversionFactor unit) {
  float baseUnitsExp[5] = {0.f};

  getBaseUnitExponantsArray(baseUnitsExp, unit);

  generalConversionString(buffer, us, baseUnitsExp);
}

/**
 * @brief Returns the conversion factor for a given unit (expressed in terms of
 * the 5 fundamental units) in the chosen unit system
 * @param us The unit system used
 * @param baseUnitsExponants The exponent of each base units required to form
 * the desired quantity. See conversionFactor() for a working example
 */
double generalConversionFactor(struct UnitSystem* us,
                               float baseUnitsExponants[5]) {
  double factor = 1.;
  int i;

  for (i = 0; i < 5; ++i)
    if (baseUnitsExponants[i] != 0)
      factor *= pow(getBaseUnit(us, i), baseUnitsExponants[i]);
  return factor;
}

/**
 * @brief Returns the h factor exponentiation for a given unit (expressed in
 * terms of the 5 fundamental units)
 * @param us The unit system used
 * @param baseUnitsExponants The exponent of each base units required to form
 * the desired quantity. See conversionFactor() for a working example
 */
float generalhFactor(struct UnitSystem* us, float baseUnitsExponants[5]) {
  float factor_exp = 0.f;

  factor_exp += -baseUnitsExponants[UNIT_MASS];
  factor_exp += -baseUnitsExponants[UNIT_LENGTH];
  factor_exp += -baseUnitsExponants[UNIT_TIME];

  return factor_exp;
}

/**
 * @brief Returns the scaling factor exponentiation for a given unit (expressed
 * in terms of the 5 fundamental units)
 * @param us The unit system used
 * @param baseUnitsExponants The exponent of each base units required to form
 * the desired quantity. See conversionFactor() for a working example
 */
float generalaFactor(struct UnitSystem* us, float baseUnitsExponants[5]) {
  float factor_exp = 0.f;

  factor_exp += baseUnitsExponants[UNIT_LENGTH];

  return factor_exp;
}

/**
 * @brief Returns a string containing the exponents of the base units making up
 * the conversion factors (expressed in terms of the 5 fundamental units)
 * @param buffer The buffer in which to write (The buffer must be long enough,
 * 140 chars at most)
 * @param us The UnitsSystem in use.
 * @param baseUnitsExponants The exponent of each base units required to form
 * the desired quantity. See conversionFactor() for a working example
 */
void generalConversionString(char* buffer, struct UnitSystem* us,
                             float baseUnitsExponants[5]) {
  char temp[14];
  double a_exp = generalaFactor(us, baseUnitsExponants);
  double h_exp = generalhFactor(us, baseUnitsExponants);
  int i;

  /* Check whether we are unitless or not */
  char isAllNonZero = 1;
  for (i = 0; i < 5; ++i)
    if (baseUnitsExponants[i] != 0.) isAllNonZero = 0;

  if (isAllNonZero) {
    sprintf(buffer, "[ - ] ");
    return;
  }

  /* Add a-factor */
  if (a_exp == 0)
    sprintf(buffer, " ");
  else if (a_exp == 1)
    sprintf(buffer, "a ");
  else if (remainder(a_exp, 1.) == 0)
    sprintf(buffer, "a^%d ", (int)a_exp);
  else
    sprintf(buffer, "a^%7.4f ", a_exp);

  /* Add h-factor */
  if (h_exp == 0)
    sprintf(temp, " ");
  else if (h_exp == 1)
    sprintf(temp, "h ");
  else if (remainder(h_exp, 1.) == 0)
    sprintf(temp, "h^%d ", (int)h_exp);
  else
    sprintf(temp, "h^%7.4f ", h_exp);
  strncat(buffer, temp, 12);

  /* Add conversion units */
  for (i = 0; i < 5; ++i)
    if (baseUnitsExponants[i] != 0) {
      if (baseUnitsExponants[i] == 0.)
        sprintf(temp, " ");
      else if (baseUnitsExponants[i] == 1.)
        sprintf(temp, "%s ", getBaseUnitSymbol(i));
      else if (remainder(baseUnitsExponants[i], 1.) == 0)
        sprintf(temp, "%s^%d ", getBaseUnitSymbol(i),
                (int)baseUnitsExponants[i]);
      else
        sprintf(temp, "%s^%7.4f ", getBaseUnitSymbol(i), baseUnitsExponants[i]);
      strncat(buffer, temp, 12);
    }

  /* Add CGS units */
  strncat(buffer, " [ ", 3);

  for (i = 0; i < 5; ++i) {
    if (baseUnitsExponants[i] != 0) {
      if (baseUnitsExponants[i] == 0.)
        continue;
      else if (baseUnitsExponants[i] == 1.)
        sprintf(temp, "%s ", getBaseUnitCGSSymbol(i));
      else if (remainder(baseUnitsExponants[i], 1.) == 0)
        sprintf(temp, "%s^%d ", getBaseUnitCGSSymbol(i),
                (int)baseUnitsExponants[i]);
      else
        sprintf(temp, "%s^%7.4f ", getBaseUnitCGSSymbol(i),
                baseUnitsExponants[i]);
      strncat(buffer, temp, 12);
    }
  }

  strncat(buffer, "]", 2);
}
