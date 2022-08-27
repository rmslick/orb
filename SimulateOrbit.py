import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = 40,30
from matplotlib.animation import FuncAnimation

#r = 1 # radius of circle
def circle(phi,r):
    return np.array([r*np.cos(phi), r*np.sin(phi)])
def show_orbit(orbit_radius, radian_frames, full_radius,orbital_object_radius=6371000):
    # create a figure with an axes
    fig, ax = plt.subplots()
    # set the axes limits
    ax.axis([-full_radius,full_radius,-full_radius,full_radius])
    # set equal aspect such that the circle is not shown as ellipse
    ax.set_aspect("equal")
    Drawing_colored_circle = plt.Circle(( 4 , 4 ), 100 )
    ax.add_artist( Drawing_colored_circle )
    orbital_center_toorbital_object = plt.Line2D((0,0),(1,1),color="red")
    ax.add_artist(orbital_center_toorbital_object)
    ln, = plt.plot([], [], 'ro-', animated=True)
    ground_station_line, = plt.plot([], [], 'ro-', animated=True)
    #plt.plot([[0,0]],[[x,y]],'ro-')
    # create a point in the axes
    point, = ax.plot(0,0.5, marker="o")
    earth_center, = ax.plot(0,0.5, marker="x")
    xfixdata = 0,0
    yfixdata = 0,0
    plt.plot([xfixdata], [yfixdata], 'bo', ms=10)

    # Updating function, to be repeatedly called by the animation
    def update(phi):
        # obtain point coordinates 
        x,y = circle(phi,orbit_radius)
        # set point's coordinates
        point.set_data([x],[y])

        # update line
        Drawing_colored_circle.set_radius(orbital_object_radius)

        ln.set_data([0,x], [0,y])
        ground_station_line.set_data([0,x],[orbital_object_radius,y])
        #plt.show()
        earth_center.set_data([int(0)],[int(0)])
        return point,earth_center,Drawing_colored_circle,orbital_center_toorbital_object,ln,ground_station_line,

    # create animation with 10ms interval, which is repeated,
    # provide the full circle (0,2pi) as parameters
    ani = FuncAnimation(fig, update, interval=10, blit=True, repeat=True,
                        frames=radian_frames )

    plt.show()