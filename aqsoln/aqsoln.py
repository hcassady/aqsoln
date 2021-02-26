

import numpy as np
import pandas as pd
import scipy as sp
import scipy.interpolate
import sqlite3
import molmass



class Solution:

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
    self._interpolate = sp.interpolate.LinearNDInterpolator(
      list(zip(x, y)), z
    )

  # def molar_TO_mass_concentration(self, molarity):
  #   mass_concentration = molarity * self.solute_molar_mass
  #   return mass_concentration

  # # def molar_TO_mass_fraction(self, molarity):

  # def solve_density_function(self, function, temperature):
  #   x_ref = self.raw_data.weight_percent
  #   y_ref = self._interpolate(x_ref, temperature)

