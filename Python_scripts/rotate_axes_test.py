from mpl_toolkits.axisartist.floating_axes import GridHelperCurveLinear, FloatingSubplot

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from mpl_toolkits.axisartist.grid_finder import FixedLocator, MaxNLocator, DictFormatter

data = np.random.normal(0, 1, 1000)

# Setting up the extremes for the grid (x_min, x_max, y_min, y_max)
extremes = (-4, 4, 0, 60)  # Limits for the x and y axes

# Create a GridHelper with no rotation
transform = mtransforms.Affine2D()
grid_helper_original = GridHelperCurveLinear(transform, extremes=extremes)

# Create a GridHelper with the rotation transform and resize the axes to half the size
transform = mtransforms.Affine2D().scale(1, 0.25).rotate_deg(-45)
grid_helper_rotated = GridHelperCurveLinear(transform, extremes=extremes)

# Create the figure and floating axes
fig = plt.figure(figsize=(6, 6))
original_ax = FloatingSubplot(fig, 111, grid_helper=grid_helper_original)
rotated_ax = FloatingSubplot(fig, 111, grid_helper=grid_helper_rotated)

fig.add_subplot(original_ax)
fig.add_subplot(rotated_ax)

rotated_ax.set_position([0.5, 0.5, 0.3, 0.3])

# Create an auxiliary axes object
aux_ax = rotated_ax.get_aux_axes(transform)

# Plot the original histogram
aux_ax.hist(data, bins=30, alpha=0.7, color='blue', edgecolor='black')

# Plot the smaller, rotated histogram over the first
small_hist_data = np.histogram(data, bins=30)
# aux_ax.plot(small_hist_data[1][1:], small_hist_data[0], color='green')

# plot the original histogram (not rotated)
original_ax.hist(data, bins=30, alpha=0.7, color='blue', edgecolor='black')

plt.show()
