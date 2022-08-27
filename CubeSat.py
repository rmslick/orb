from math import asin, atan, cos, pi, sin, sqrt
import numpy as np
from SimulateOrbit import show_orbit
class OrbitalPosition:
    def __init__(self, phi, aoe, aoe_rc):
        # angle from ground station
        self.phi = phi
        # angle of elevation
        self.aoe = aoe
        # angle of elevation rate of change
        self.aoe_rc = aoe_rc

class Orbit:
    def __init__(self,P=100,d_p=1,d_a=1,r=1):
        # Period of rotation (mins)
        self.P = P
        # perihelion
        self.d_p = d_p
        # aphelion
        self.d_a=d_a
        # radius r (https://openstax.org/books/university-physics-volume-1/pages/13-5-keplers-laws-of-planetary-motion#:~:text=areal%20velocity%20%3D%20%CE%94%20A%20%CE%94,is%20exactly%20Kepler's%20second%20law.)
        self.r=r
        # Angular velocity - Rate of change (https://math.stackexchange.com/questions/1057115/angular-velocity-around-ellipse) dangle/dt
        # (radians/(mins*meters*meters))*meters*meters
        # radians/min
        self.o_dot = (pi/(self.P*self.r**2)) *(self.d_a+self.d_p)*sqrt(self.d_a*self.d_p)
class CubeSat:
    def __init__(self):
        pass
# single axis sube sat
# Dev note: Factor out pitch specific to a maneuver class
class SingleAxisCubeSat:
    def __init__( self, J , Orbit , A = np.array( [[0,1],[0,0]] ) ) :
        # Inertial matrix
        self.J_m = None
        # Moment of intertia
        self.J = J
        # Pitch transition matrix
        self.A = A
        self.torque_transition = np.array([0,1/self.J])

        '''
            Orbital Dynamics: Position properties - Angle of elevation and pitch angle
        '''
        self.orbit = Orbit
        # Angle of elevation
        self.e = 0
        self.e_dot = 0
        # Angle between the ground station and the satellite from the center of the Earth, 𝜑.
        self.phi = 0
        # Supporting properties
        # radius of the earth (m)
        self.Re = 6371000 
        # altitude of the satellite (m)
        self.h = 400*1000
        # rate of change (assumed constant)
        #self.theta = None
        self.pitch = 0
    
        '''
            Drag properties
        '''
        # Atmospheric density -  (kg/m^3)
        # accurate for orbit at 400km orbit
        self.p = 2.8*(10**(-1*12))
        # Drag coefficient
        self.C_d = 2.5
        # Velocity
        self.V = None
        # Dimensional constant
        self.summa_R_A = None
        # gravitational constant of eaarth (9.807 m/s²)
        self.u = 9.807
    # Equations of angle
    def getSingleAxisMotion(self,ang_rate, ang_vel, torque_command, disturbance_torque):
        return self.A * np.array([ang_rate,ang_vel])+self.torque_transition*torque_command+ self.torque_transition * disturbance_torque
    # Equations of position
    def getAngleOfElevationROC(self):
        return ( (self.Re+self.h)*sin(self.e+ self.phi) * self.orbit.o_dot ) / (self.Re*sin(self.e) - (self.Re+self.h)*sin(self.e+self.phi) ) 
    def getAngleOfElevation(self):
        return atan( ( -1*self.Re+ (self.Re+self.h) * cos(self.phi) ) ) / ((self.Re +self.h)*sin( self.phi )) 
    '''
    def getPitchAngle(self):
        return asin(self.Re/(self.Re+self.h)*cos(self.e))
    def getPosition(self):
        return self.getAngleOfElevation(), self.getPitchAngle()
    # Equations of Drag
    def getAtmosphericDrag(self):
        return -1*0.5*self.C_d*self.p*self.V**2*self.summa_R_A
    def getMagneticDisturbance(self):
        pass
    def getGravityGradient(self, R):
        return (3*self.u/R**3)(R*self.J*R)
    '''
    def simulate(self,start_time,end_time,interval):
        # change angle between satellite and ground station
        t = start_time
        all_positions = []
        while t <= end_time:
            self.phi = t*self.orbit.o_dot
            
            if self.phi > 2*pi:
                self.phi -= 2*pi
            
            # set angle of elevation
            self.e = self.getAngleOfElevation()
            self.e_dot = self.getAngleOfElevationROC()
            #self.pitch = asin((self.Re+(self.Re+self.h)*cos(self.e)))
            
            #op = OrbitalPosition(self.phi, self.e, self.e_dot)
            all_positions.append(self.phi)
            t+= interval
        return all_positions


# Define an orbit
o = Orbit(100, 6371000+400*1000 , 6371000+400*1000, 6371000+400*1000)
print(o.o_dot)
s = SingleAxisCubeSat(5,o)
pos = s.simulate(1,100,1)

#show_orbit(1,np.linspace(0,2*np.pi,360, endpoint=False))
show_orbit(1,pos)