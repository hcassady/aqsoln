import numpy as np
import pandas as pd
import scipy as sp
from scipy.interpolate import LinearNDInterpolator
from scipy.optimize import fsolve
import scipy.optimize
import sqlite3
import molmass



class Calculator:

  def __init__(self, salt, check_bounds=True):
    self.salt = salt
    self.solute_molar_mass = molmass.Formula(salt).mass
    self._load_data()
    self._initialize_interpolator()


  def _load_data(self):
    query = f"""
      --begin-sql
      SELECT
        weight_percent,  -- (g solute / g total) x 100%
        temperature,     -- (deg C)
        density,         -- (g total / mL total)
        source
      FROM
        aqueous_solution_densities
      WHERE (salt = ?
        OR salt = 'None')
      AND temperature <= 100;
    """
    conn = sqlite3.connect('physical_data.db')
    self.raw_data = pd.read_sql(
      query, conn, params=(self.salt,))
    conn.close()

    self.sources = self.raw_data.source.unique()

  def _initialize_interpolator(self):
    """ Generate interpolator function.

    This creates an interpolator function of the LinearNDInterpolator type,
    which takes the data points (weight_percent, temperature) tuples as an input
    and outputs (density).
    """
    x = self.raw_data.weight_percent.to_numpy()
    y = self.raw_data.temperature.to_numpy()
    z = self.raw_data.density.to_numpy()
    self._interpolate = LinearNDInterpolator(
      list(zip(x, y)), z
      )

  def molar_to_weight_percent(self, molarity, temperature=20):
    """Molarity (mol/L) to weight percent.

    Input a molarity (in mol/L) and a temperature, and output the
    weight percent.
    """
    molarity_per_mL = molarity / 1000    # mol/mL
    def zero_function(weight_percent):
      zero = (
        ((molarity_per_mL * self.solute_molar_mass) /
          self._interpolate(weight_percent, temperature)) - weight_percent / 100
      )
      return zero
    # Initial guess based on density of 1 g/mL
    initial_guess = molarity_per_mL * self.solute_molar_mass * 100
    weight_percent = fsolve(zero_function, [initial_guess])
    return weight_percent

  def molar_to_mass_fraction(self, molarity, temperature=20):
    return self.molar_to_weight_percent(molarity, temperature) / 100

  def molar_to_density(self, molarity, temperature=20):
    """Molar concentration to density.

    Returns the solution density (g /mL) of the solution for a given molarity
    (mol solute / L solution) and temperature (deg C).
    """
    density = self._interpolate(
      self.molar_to_weight_percent(molarity, temperature), temperature
    )
    return density
  
  def molar_to_molal(self, molarity, temperature=20):
    """Molar concentration to molality.

    Returns the molality (mol solute / kg solvent) of the solution for a given
    molarity (mol solute / L solution) and temperature (deg C).
    """
    w = molar_to_mass_fraction(molarity, temperature)
    Mw = self.solute_molar_mass / 1000  # kg/mol
    molality = w / ((1 - w) *  Mw)      # mol solute / kg solvent
    return molality
