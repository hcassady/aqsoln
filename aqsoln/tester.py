#%%
import aqsoln as aq
from importlib import reload
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go

# %%
# Generate NaCl tester

# %%
reload(aq)
NaCl_test = aq.Solution("NaCl")

# Test the interpolator by plotting:
x_reg = np.linspace(0, 26, 420)
y_reg = np.linspace(0, 100, 420)
X_reg, Y_reg = np.meshgrid(x_reg, y_reg)
Z_reg = NaCl_test._interpolate(X_reg, Y_reg)

fig = go.Figure(data=[
  go.Surface(x=X_reg, y=Y_reg, z=Z_reg),
  go.Scatter3d(
    x = NaCl_test.raw_data.weight_percent,
    y = NaCl_test.raw_data.temperature,
    z = NaCl_test.raw_data.density,
    mode='markers',
    marker=dict(size=2, color='red')
    )
  ])
fig.show()
# %%
