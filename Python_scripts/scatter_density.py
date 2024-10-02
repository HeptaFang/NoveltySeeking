import matplotlib.pyplot as plt
import numpy as np

# Example data
x = np.random.rand(10000)
y = np.random.rand(10000)

plt.figure(figsize=(8, 6))

# Create 2D histogram
hist, xedges, yedges, image = plt.hist2d(x, y, bins=200, cmap='gist_yarg', clim=[0, 10])

# Add color bar to indicate density levels
cbar = plt.colorbar(image, label='Density')

plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('2D Histogram Density Plot')

# add axes
plt.axhline(0, color='black', linewidth=1.5)
plt.axvline(0, color='black', linewidth=1.5)
plt.tight_layout()
plt.gca().set_aspect('equal', adjustable='box')

plt.show()
