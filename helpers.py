# import pltting library
from matplotlib import pyplot as plt
import numpy as np

def draw_ball_trajectory(all_x, all_y): 
    
    print(f"Throwing distance: {all_x[-1]:0.2f} m")
    
    plt.figure() 
    plt.plot(all_x, all_y, '-o')
    # plot floor
    all_x = np.array(all_x)
    all_y = np.array(all_y)
    xmin = min(-0.5, all_x.min() - 0.5)
    xmax = max(10., all_x.max() + 0.5)
    ymin = min(-0.5, all_y.min() - 0.5)
    ymax = max(6., all_y.max() + 0.5)
    # approx equal axes 
    if ymax < 0.6 * xmax: 
        ymax = 0.6*xmax
    else: 
        xmax = ymax / 0.6
    

    plt.plot([xmin, xmax],[0., 0.], color='black')
    plt.axis([xmin, xmax, ymin, ymax])
    plt.xlabel('x')
    plt.ylabel('y')


def plot_wave(ax, cell_centers, height, cell_width, time):
    ax.bar(cell_centers, height, width=cell_width)
    ax.set_xlabel('x')
    ax.set_title(f't = {time:0.2f}')