import imageio as iio
import os
import numpy as np
import matplotlib.pyplot as plt

def draw_scatter(s, num, data_path, image_path):
    """
    param s : Scatter size
    image_path : image path
    return : None
    """
    data = np.loadtxt(data_path, encoding='utf-8', delimiter='\t')
	# Get the abscissa x by slicing
    nx1 = data[0:672, 0]           # Leptodactylidae
    nx2 = data[673:1214, 0]        # Dendrobatidae
    nx3 = data[1215:4692, 0]       # Leptodactylidae
    nx4 = data[4693:6595, 0]       # Hylidae
    nx5 = data[6596:6865, 0]       # Leptodactylidae
    nx6 = data[6866:6979, 0]       # Hylidae
    nx7 = data[6980:7047, 0]       # Bufonidae
    nx8 = data[7048:7195, 0]       # Hylidae
	# Get ordinate y by slice
    ny1 = data[0:672, 1]           # Leptodactylidae
    ny2 = data[673:1214, 1]        # Dendrobatidae
    ny3 = data[1215:4692, 1]       # Leptodactylidae
    ny4 = data[4693:6595, 1]       # Hylidae
    ny5 = data[6596:6865, 1]       # Leptodactylidae
    ny6 = data[6866:6979, 1]       # Hylidae
    ny7 = data[6980:7047, 1]       # Bufonidae
    ny8 = data[7048:7195, 1]       # Hylidae
	# Create a drawing window
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(1,1,1)
	# set title
    ax.set_title('npca')
	# Set the abscissa name
    ax.set_xlabel('pc1')
	# Set the ordinate name
    ax.set_ylabel('pc2')
	# Draw a scatterplot
    ax.scatter(nx1, ny1, s = s, c='k')
    ax.scatter(nx2, ny2, s = s, c='b')
    ax.scatter(nx3, ny3, s = s, c='k')
    ax.scatter(nx4, ny4, s = s, c='r')
    ax.scatter(nx5, ny5, s = s, c='k')
    ax.scatter(nx6, ny6, s = s, c='r')
    ax.scatter(nx7, ny7, s = s, c='g')
    ax.scatter(nx8, ny8, s = s, c='r')
	# Adjust the upper and lower bounds of the horizontal and vertical coordinates
    plt.xlim(xmax=15, xmin=-15)
    plt.ylim(ymax=10, ymin=-20)
	# show
    plt.savefig('%s/%s.png' % (image_path, num))
    plt.close()

def draw_multipng(data_path, image_path):
    '''
    Draw multiple results in the npca dimensionality reduction process
    data_path : result path
    '''
    pathDir = os.listdir(data_path)
    for eachdata in pathDir:
        eachDir = os.path.join('%s/%s' % (data_path, eachdata))
        draw_scatter(s=0.5, num=eachdata, data_path=eachDir, image_path=image_path)

def png_to_gif(image_path, gif_name, duration):
    '''
    image_path : Image save path
    gif_name : gif save path name
    duration : Intervals
    '''
    # frames used to save the picture list
    frames = []
    # Get all pictures in the picture path
    image_list = os.listdir(image_path)
    image_list.sort(key=lambda x:int(x[:-4]))
    # Imageio reads in pictures and saves them in frames list 
    for image in image_list:
        image = image_path + '/' + image
        frames.append(iio.v2.imread(image))
    # create gif
    iio.mimsave(gif_name, frames, 'GIF', duration=duration)



def main():
    data_path = '/home/ljz/VSCode_PY/npca/code/png_to_gif'
    image_path = '/home/ljz/VSCode_PY/npca/code/png'
    gif_name = '/home/ljz/VSCode_PY/npca/code/test2.gif'
    draw_multipng(data_path, image_path)
    png_to_gif(image_path, gif_name, 0.5)

if __name__ == '__main__':
    main()
