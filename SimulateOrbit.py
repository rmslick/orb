import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = 4,3
from matplotlib.animation import FuncAnimation

#r = 1 # radius of circle
def circle(phi,r):
    return np.array([r*np.cos(phi), r*np.sin(phi)])
def show_orbit(r,radian_frames):
    # create a figure with an axes
    fig, ax = plt.subplots()
    # set the axes limits
    ax.axis([-1.5,1.5,-1.5,1.5])
    # set equal aspect such that the circle is not shown as ellipse
    ax.set_aspect("equal")
    # create a point in the axes
    point, = ax.plot(0,1, marker="o")

    # Updating function, to be repeatedly called by the animation
    def update(phi):
        # obtain point coordinates 
        x,y = circle(phi,r)
        # set point's coordinates
        point.set_data([x],[y])
        return point,

    # create animation with 10ms interval, which is repeated,
    # provide the full circle (0,2pi) as parameters
    ani = FuncAnimation(fig, update, interval=10, blit=True, repeat=True,
                        frames=radian_frames )

    plt.show()