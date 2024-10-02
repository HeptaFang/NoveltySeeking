import numpy as np
import scipy.io
from matplotlib import pyplot as plt

def main():
    all_diffs = scipy.io.loadmat('GLM_data/all_diffs.mat')
    all_Js = scipy.io.loadmat('GLM_data/all_Js.mat')
    fig, axs = plt.subplots(3, 3, figsize=(15, 15))
    
    # cmap = 'gist_yarg'
    cmap = 'viridis'
    
    for i in range(3):
        for j in range(3):
            ax = axs[i, j]
            fast_Js = all_Js['all_Js'][i, j, 0].flatten()
            slow_Js = all_Js['all_Js'][i, j, 1].flatten()
            
            # axs[i, j].plot(fast_Js, slow_Js, 'o')
            hist, xedges, yedges = np.histogram2d(fast_Js, slow_Js, bins=100, range=[[0, 2], [0, 2]])
            hist = np.log1p(hist)  # Use log1p to avoid log(0) issues
            image = ax.imshow(hist.T, origin='lower', aspect='auto', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap=cmap)
            ax.plot([0, 2], [0, 2], 'k-')
            
            
            ax.set_aspect('equal', adjustable='box')
            ax.set_xlim(0, 2)
            ax.set_ylim(0, 2)
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()