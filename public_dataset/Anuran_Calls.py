import matplotlib.pyplot as plt
import numpy as np

def draw_scatter(s):
    """
    param s : Scatter size
    return : None
    """
    data1 = np.loadtxt('/home/ljz/CLionProjects/npca/pc.pca', encoding='utf-8', delimiter='\t')
    data2 = np.loadtxt('/home/ljz/CLionProjects/npca/pc.npca', encoding='utf-8', delimiter='\t')
	# Get the abscissa x by slicing
    x1 = data1[0:672, 0]            # Leptodactylidae
    x2 = data1[673:1214, 0]         # Dendrobatidae
    x3 = data1[1215:4692, 0]        # Leptodactylidae
    x4 = data1[4693:6595, 0]        # Hylidae
    x5 = data1[6596:6865, 0]        # Leptodactylidae
    x6 = data1[6866:6979, 0]        # Hylidae
    x7 = data1[6980:7047, 0]        # Bufonidae
    x8 = data1[7048:7195, 0]        # Hylidae
    nx1 = data2[0:672, 0]           # Leptodactylidae
    nx2 = data2[673:1214, 0]        # Dendrobatidae
    nx3 = data2[1215:4692, 0]       # Leptodactylidae
    nx4 = data2[4693:6595, 0]       # Hylidae
    nx5 = data2[6596:6865, 0]       # Leptodactylidae
    nx6 = data2[6866:6979, 0]       # Hylidae
    nx7 = data2[6980:7047, 0]       # Bufonidae
    nx8 = data2[7048:7195, 0]       # Hylidae
	# Get ordinate y by slice
    y1 = data1[0:672, 1]            # Leptodactylidae
    y2 = data1[673:1214, 1]         # Dendrobatidae
    y3 = data1[1215:4692, 1]        # Leptodactylidae
    y4 = data1[4693:6595, 1]        # Hylidae
    y5 = data1[6596:6865, 1]        # Leptodactylidae
    y6 = data1[6866:6979, 1]        # Hylidae
    y7 = data1[6980:7047, 1]        # Bufonidae
    y8 = data1[7048:7195, 1]        # Hylidae
    ny1 = data2[0:672, 1]           # Leptodactylidae
    ny2 = data2[673:1214, 1]        # Dendrobatidae
    ny3 = data2[1215:4692, 1]       # Leptodactylidae
    ny4 = data2[4693:6595, 1]       # Hylidae
    ny5 = data2[6596:6865, 1]       # Leptodactylidae
    ny6 = data2[6866:6979, 1]       # Hylidae
    ny7 = data2[6980:7047, 1]       # Bufonidae
    ny8 = data2[7048:7195, 1]       # Hylidae
	# Create a drawing window
    fig = plt.figure(figsize=(10, 5))
	# Divide the drawing window into 1 row and 2 columns
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)
	# set title
    ax1.set_title('pca')
    ax2.set_title('npca')
	# Set the abscissa name
    ax1.set_xlabel('pc1')
    ax2.set_xlabel('pc1')
	# Set the ordinate name
    ax1.set_ylabel('pc2')
    ax2.set_ylabel('pc2')
	# Draw a scatterplot
    ax1.scatter(x1, y1, s = s, c='k', label='Leptodactylidae')
    ax1.scatter(x2, y2, s = s, c='b', label='Dendrobatidae')
    ax1.scatter(x3, y3, s = s, c='k')
    ax1.scatter(x4, y4, s = s, c='r', label='Hylidae')
    ax1.scatter(x5, y5, s = s, c='k')
    ax1.scatter(x6, y6, s = s, c='r')
    ax1.scatter(x7, y7, s = s, c='g', label='Bufonidae')
    ax1.scatter(x8, y8, s = s, c='r')
    ax2.scatter(nx1, ny1, s = s, c='k')
    ax2.scatter(nx2, ny2, s = s, c='b')
    ax2.scatter(nx3, ny3, s = s, c='k')
    ax2.scatter(nx4, ny4, s = s, c='r')
    ax2.scatter(nx5, ny5, s = s, c='k')
    ax2.scatter(nx6, ny6, s = s, c='r')
    ax2.scatter(nx7, ny7, s = s, c='g')
    ax2.scatter(nx8, ny8, s = s, c='r')
    # add legend
    ax1.legend(markerscale=6, loc='upper right', ncol=4, bbox_to_anchor=(1.9, 1.15))
	# show
    plt.show()
 
 
# main
if __name__ == "__main__":
	# run
	draw_scatter(s=0.5)
