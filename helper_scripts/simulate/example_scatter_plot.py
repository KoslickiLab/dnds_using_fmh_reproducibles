#import pandas as pd
import matplotlib.pyplot as plt
#import scipy.stats
import numpy as np
from matplotlib import rcParams

# Set font to Times New Roman
rcParams['font.family'] = 'serif'

#Create random values
x = np.random.uniform(0, 8, 100)
y = x + np.random.normal(0, 0.5, 100)  # Add small noise to make them close to 1:1


# Create a figure and axis
fig, ax = plt.subplots()

#Plot values
plt.scatter(x, y, alpha=0.7, color='black')
plt.plot([0, 8], [0, 8], color='gray', linestyle='--', label='1:1 Line')  # Add a 1:1 line

# Shading regions as specified
plt.fill_betweenx(y=np.linspace(1, 8, 100), x1=1, x2=8, color='#ff7f0e', alpha=0.3)  # Light red region
plt.fill_betweenx(y=np.linspace(0, 1, 100), x1=0, x2=1, color='#add8e6', alpha=0.3)  # Light blue region

# Set x and y axis limits
ax.set_xlim(-2, 8)
ax.set_ylim(-2, 8)

# Set specific ticks for x and y axes
ax.set_xlim(0, 10)
ax.set_xticks([0, 1, 2, 5, 8, 10])

ax.set_ylim(0, 9)
ax.set_yticks([0, 1, 2, 5, 8, 10])

# Set labels for axes and title
fs=18
dnds = r"$\mathrm{FMH}\ d_{\mathrm{N}}/d_{\mathrm{S}}$"
ax.set_xlabel('FracMinHash '+dnds,fontsize=fs)
ax.set_ylabel('Ground Truth '+dnds,fontsize=fs)
ax.set_title("Sequence Simulation Sample",fontsize=fs)

# Adjust layout
fig.figure.savefig(f"../..//manuscript_figures/updated_pdf/example_scatter_plot.pdf",bbox_inches='tight') 
