import numpy as np
import scipy.io
from matplotlib import pyplot as plt

def main():
    all_diffs = scipy.io.loadmat('GLM_data/all_diffs.mat')
    all_Js = scipy.io.loadmat('GLM_data/all_Js.mat')
    fig, axs = plt.subplots(3, 3, figsize=(15, 15))
    
    for i in range(3):
        for j in range(3):
            fast_Js = all_Js['all_Js'][i, j, 0]
            slow_Js = all_Js['all_Js'][i, j, 1]
            
            axs[i, j].plot(fast_Js, slow_Js, 'o')
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()