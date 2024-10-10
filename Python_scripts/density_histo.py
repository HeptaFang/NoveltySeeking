import numpy as np
import scipy.io
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as mcolors

def main():
    all_diffs = scipy.io.loadmat('GLM_data/all_diffs.mat')
    all_Js = scipy.io.loadmat('GLM_data/all_Js.mat')
    # all_Js = scipy.io.loadmat('GLM_data/all_Js_abs.mat')
    fig, axs = plt.subplots(3, 3, figsize=(15, 15))
    
    # cmap = 'viridis'
    # cmap = 'gray'
    # cmap = 'gist_yarg'

    # define a colormap from white to red
    cmap = mcolors.LinearSegmentedColormap.from_list('mycmap', ['white', 'red'], N=256)
    
    # Prepare for the colorbar by tracking the max log(density)
    max_log_density = -np.inf

    area_names = ['ACC', 'Thalamus', 'VLPFC']
    
    for i in range(3):
        for j in range(3):
            ax = axs[i, j]
            fast_Js = all_Js['all_Js'][i, j, 0].flatten()
            slow_Js = all_Js['all_Js'][i, j, 1].flatten()
            
            # Create 2D histogram
            hist, xedges, yedges = np.histogram2d(fast_Js, slow_Js, bins=100, range=[[-2, 2], [-2, 2]])
            hist_log = np.log1p(hist)  # Use log1p to avoid log(0) issues
            max_log_density = max(max_log_density, hist_log.max())  # Track the maximum log(density)
            
            # Plot the 2D histogram
            image = ax.imshow(hist_log.T, origin='lower', aspect='auto', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap=cmap)
            ax.plot([-2, 2], [-2, 2], 'k--', linewidth=1)
            ax.set_aspect('equal', adjustable='box')
            ax.set_xlim(0, 2)
            ax.set_ylim(0, 2)
            ax.set_xticks([-2, 0, 2])
            ax.set_yticks([-2, 0, 2])

            # add gray axes at 0
            ax.axhline(0, color='gray', linewidth=1)
            ax.axvline(0, color='gray', linewidth=1)

            # Add histograms along x and y axes
            divider = make_axes_locatable(ax)
            ax_hist_x = divider.append_axes("bottom", 0.5, pad=0.3, sharex=ax)
            ax_hist_y = divider.append_axes("left", 0.5, pad=0.3, sharey=ax)
            
            # Plot the x histogram
            ax_hist_x.hist(fast_Js, bins=50, range=(-2, 2), color='black')
            # ax_hist_x.set_xticks([])
            ax_hist_x.invert_yaxis()

            # Plot the y histogram
            ax_hist_y.hist(slow_Js, bins=50, range=(-2, 2), color='black', orientation='horizontal')
            # ax_hist_y.set_yticks([])
            ax_hist_y.invert_xaxis()

            # Label x and y axes
            ax_hist_x.set_xlabel('J_fast')
            ax_hist_y.set_ylabel('J_slow')

            # Set a title
            ax.set_title(f'{area_names[j]} to {area_names[i]}')

    # Add a global colorbar for the entire figure
    cbar_ax = fig.add_axes([0.92, 0.1, 0.02, 0.8])  # Position the colorbar beside the plots
    norm = mcolors.Normalize(vmin=0, vmax=max_log_density)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # Required for colorbar
    cbar = plt.colorbar(sm, cax=cbar_ax)
    cbar.set_label('Density')
    cbar.set_ticks(np.log1p([0, 1, 10, 100, 1000]))
    cbar.set_ticklabels([0, 1, 10, 100, 1000])

    # Add a global title
    fig.suptitle('Density plots of J_fast vs J_slow', fontsize=16)

    plt.tight_layout(rect=[0, 0, 0.95, 0.95], pad=3.5)  # Adjust layout to make space for the colorbar and title
    plt.savefig('figures/kernel_comp_density.png', dpi=200)
    # plt.show()
    
    plt.close(fig)

if __name__ == "__main__":
    main()
