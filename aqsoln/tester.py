#%%
import aqsoln as aq
from importlib import reload
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# %%
# Generate NaCl tester
NaCl_test = aq.Solution("NaCl")
# %%
fig = plt.figure()
x_plot = np.linspace(0, 26, 6969)
def add_temp_to_plot(T):
  Ts = np.full(len(x_plot), T)
  y_plot = NaCl_test._interpolate(x_plot, Ts)
  plt.plot(x_plot, y_plot, '-')
add_temp_to_plot(22)
add_temp_to_plot(30)
add_temp_to_plot(40)
add_temp_to_plot(50)
fig.show()
# plt.plot(x, y, 'o', color='bo', markersize=1);

# %%
