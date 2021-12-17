# To get an elliptical orbit, we can change the velocity
# Implement RK4
# Check the formula for velocity without centripetal motion

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import time

#************************************** SIMULATION PARAMETERS ***************************************#
dt = 1000              # Timestep in seconds
N = 100000             # Number of timesteps to be calulated. dt=86400(1 Day) and N=365(1 year) => 1 year
framerate = 2000        # Number of frames
#****************************************************************************************************#

#*********************************** DEFINING UNIVERSAL CONSTANTS ***********************************#
G = 6.674e-11                   # Gravitation Constant in SI unit (m3 kg-1 s-2)
M = 1.989e30                    # Mass of the Sun in kg
AU = 1.496e11                   # 1 Astronomical Unit in m
#****************************************************************************************************#

#************************************* DEFINING VALUES FOR EARTH ************************************#
Mass_E = 5.972e24                          # Mass of Earth
Dist_E = 1*AU                              # Distance of Earth from Sun
# Pos_E = np.array([AU,0])                   # Initial position of Earth
Pos_E = np.array([0.98329*AU,0])                   # Initial position of Earth
# Vel_E = np.array([0,-np.sqrt(G*M/Pos_E[0])])  # Initial velocity of Earth
Vel_E = np.array([0,-30290])  # Initial velocity of Earth
x_earth = []                               # This array contains all x-components of Earth's position
y_earth = []                               # This array contains all y-component of Earth's position
#****************************************************************************************************#

#*********************************** DEFINING VALUES FOR JUPITER ************************************#
Mass_J = 1.898e27                          # Mass of Jupiter
Rad_J = 69911e3                            # Radius of Jupiter in m

'''Angle phi is calculated such that Jupiter meets the rocket at the apogee'''
phi = -1.72                                # Angle between initial position of Earth and Jupiter
Dist_J = 5.2044*AU                          # Distance of Jupiter from Sun

'''Initial position of Jupiter is calculated by splitting the Dist_J in components using angle phi'''
Pos_J = np.array([Dist_J*np.cos(phi),Dist_J*np.sin(phi)])  # Initial position of Jupiter

'''Inital velocity of jupiter calculated using the orbital speed and then dividing it into components'''
Vel_J_Mod = np.sqrt(G*M/np.sqrt(Pos_J.dot(Pos_J))) # Orbital speed of Jupiter
Vel_J = np.array([np.cos(phi-1.5708)*Vel_J_Mod, np.sin(phi-1.5708)*Vel_J_Mod]) # Initial velocity of Jupiter

x_jupiter = []                             # x-component of Jupiter's position
y_jupiter = []                             # y-component of Jupiter's position
#****************************************************************************************************#

#*********************************** DEFINING VALUES FOR SATURN ************************************#
Mass_S = 5.683e26                          # Mass of Saturn

'''Angle phi is calculated such that Jupiter meets the rocket at the apogee'''
phi1 = -2.75                               # Angle between initial position of Earth and Saturn
Dist_S = 9.5857*AU                         # Distance of Saturn from Sun

'''Initial position of Jupiter is calculated by splitting the Dist_J in components using angle phi'''
Pos_S = np.array([Dist_S*np.cos(phi1),Dist_S*np.sin(phi1)])  # Initial position of Saturn

'''Inital velocity of jupiter calculated using the orbital speed and then dividing it into components'''
Vel_S_Mod = np.sqrt(G*M/np.sqrt(Pos_S.dot(Pos_S))) # Orbital speed of Jupiter
Vel_S = np.array([np.cos(phi1-1.5708)*Vel_S_Mod, np.sin(phi1-1.5708)*Vel_S_Mod]) # Initial velocity of Jupiter

x_saturn = []                             # x-component of Saturn's position
y_saturn = []                             # y-component of Saturn's position
#****************************************************************************************************#

#************************************ DEFINING VALUES FOR ROCKET ************************************#
rR2 = 5*Rad_J                              # Distance of rocket from the center of Jupiter in m
'''The following formula needs to be revised'''
# vR = np.array([0, -np.sqrt(G*M*2*(Dist_J+rR2)/(Pos_E[0]*(Pos_E[0]+(Dist_J+rR2))))]) # Initial velocity of rocket
vR = np.array([0,-np.sqrt(G*M/Pos_E[0])+3000])
x_rocket = []                              # x-component of rocket's position
y_rocket = []                              # y-component of rocket's position
#****************************************************************************************************#


class planet:

    def __init__(self, name, m, r, v):
        self.name = name
        self.r = r
        self.r_mag_array = []
        self.r_array = []
        self.m = m
        self.F = np.array([0, 0])
        self.F_array = []
        self.a = np.array([0, 0])
        self.a_array = []
        self.v = v
        self.v_array = []

        t = 0
        while t < N:
            # Calculating force
            r_mag = np.sqrt(self.r.dot(self.r))
            if self.name == 'Earth':
                self.r_mag_array.append(r_mag)
            r_hat = self.r/r_mag
            self.F = -G*M*self.m*r_hat/self.r.dot(self.r)
            self.F_array.append(self.F)

            # Calculating acceleration
            self.a = self.F/self.m          # ?
            self.a_array.append(self.a)

            # Calculating velocity
            self.v = self.v + self.a*dt     # ?
            self.v_array.append(self.v)

            # Calculating position
            self.r = self.r + self.v*dt     # ?
            self.r_array.append(self.r)
            if self.name == 'Earth':
                x_earth.append(self.r[0])
                y_earth.append(self.r[1])
            if self.name == 'Jupiter':
                x_jupiter.append(self.r[0])
                y_jupiter.append(self.r[1])
            if self.name == 'Saturn':
                x_saturn.append(self.r[0])
                y_saturn.append(self.r[1])
            t = t+1

class rocket:

    def __init__(self, mR, vR):
        i = 0
        closest_approach_J = 10e15
        closest_approach_S = 10e15
        self.rR2 = 0
        self.rR2_array = []
        self.rR3 = 0
        self.rR3_array = []
        self.r = Pos_E # + rR1 # Distance from sun
        self.r_array = []
        self.m = mR
        self.F1 = np.array([0, 0])
        self.F2 = np.array([0, 0])
        self.F3 = np.array([0, 0])
        self.F4 = np.array([0, 0])
        self.F = np.array([0, 0])
        self.F_array = []
        self.a = np.array([0, 0])
        self.a_array = []
        self.v = vR
        self.v_array = []
        self.speed_array = []

        t = 0
        while t < N:
            # Calculating force
            r_mag = np.sqrt(self.r.dot(self.r))
            r_hat = self.r/r_mag
            self.F1 = -G*M*self.m*r_hat/self.r.dot(self.r) #Force due to sun
            
            #Force due to jupiter
            self.rR2 = self.r-jupiter.r_array[t]
            rR2_mag = np.sqrt(self.rR2.dot(self.rR2))
            self.rR2_array.append(rR2_mag)
            if rR2_mag < closest_approach_J : closest_approach_J = rR2_mag
            rR2_hat = self.rR2/rR2_mag
            self.F3 = -G*Mass_J*self.m*rR2_hat/self.rR2.dot(self.rR2)

            #Force due to saturn
            self.rR3 = self.r-saturn.r_array[t]
            rR3_mag = np.sqrt(self.rR3.dot(self.rR3))
            self.rR3_array.append(rR3_mag)
            if rR3_mag < closest_approach_S : closest_approach_S = rR3_mag
            rR3_hat = self.rR3/rR3_mag
            self.F4 = -G*Mass_S*self.m*rR3_hat/self.rR3.dot(self.rR3)

            self.F = self.F1 + self.F2 + self.F3 + self.F4
            self.F_array.append(self.F)

            # Calculating acceleration
            self.a = self.F/self.m
            self.a_array.append(self.a)

            # Calculating velocity
            self.v = self.v + self.a*dt
            # Extra velocity added at Jupiter
            # if r_mag > Dist_J and i==0:
            #     print(self.v)
            #     self.v = self.v + [-1000, 6000]
            #     print(self.v)
            #     i = 1

            # Extra velocity added at Saturn
            # if r_mag > Dist_S and i==0:
            #     print(self.v)
            #     self.v = self.v + [1000, -500]
            #     print(self.v)
            #     i = 1
            self.speed_array.append(np.sqrt(self.v[0]**2 + self.v[1]**2))
            self.v_array.append(self.v)

            # Calculating position
            self.r = self.r + self.v*dt
            self.r_array.append(self.r)
            x_rocket.append(self.r[0])
            y_rocket.append(self.r[1])

            t = t+1
        
        print("Distance of closest approach to Jupiter:", closest_approach_J)
        print("Distance of closest approach to Saturn:", closest_approach_S)

earth = planet('Earth', Mass_E, Pos_E, Vel_E)
jupiter = planet('Jupiter', Mass_J, Pos_J, Vel_J)
saturn = planet('Saturn', Mass_S, Pos_S, Vel_S)
rocket = rocket(722, vR)

# Making graph
fig, ax = plt.subplots()
range = 10*AU
plt.xlim([-range, range])
plt.ylim([-range, range])
ax.plot([-range, range], [0,0], lw=0.5)
rocket_speed = ax.text(0.4, 0.9, "Speed: " + str(int(rocket.speed_array[0])), transform=ax.transAxes)
rocket_dist = ax.text(0.4, 0.85, "Dist. to Jupiter: " + str(int(rocket.rR2_array[0])), transform=ax.transAxes)

Earth, = ax.plot([], [], lw=2, marker='o', label='Earth', color='blue')
EarthOrbit, = ax.plot([], [], lw=0.5, color='blue')
Jupiter, = ax.plot([], [], lw=2, marker='o', label='Jupiter', color='orange')
JupiterOrbit, = ax.plot([], [], lw=0.5, color='orange')
Saturn, = ax.plot([], [], lw=2, marker='o', label='Saturn', color='brown')
SaturnOrbit, = ax.plot([], [], lw=0.5, color='brown')
Rocket, = ax.plot([], [], lw=2, marker='o', label='Rocket', color='red')
RocketOrbit, = ax.plot([], [], lw=0.5, color='red')
x1, y1 = [], []
x2, y2 = [], []
x3, y3 = [], []
x4, y4 = [], []
def animate(i):
    j = int(N * i / framerate)
    if rocket.rR2_array[j] < 2e10: time.sleep(0.4)
    if rocket.rR3_array[j] < 2e10: time.sleep(0.4)
    x, y = earth.r_array[j]
    x1.append(x)
    y1.append(y)
    Earth.set_data(x, y)
    EarthOrbit.set_data(x1, y1)
    x, y = jupiter.r_array[j]
    x2.append(x)
    y2.append(y)
    Jupiter.set_data(x, y)
    JupiterOrbit.set_data(x2, y2)
    x, y = saturn.r_array[j]
    x3.append(x)
    y3.append(y)
    Saturn.set_data(x, y)
    SaturnOrbit.set_data(x3, y3)
    x, y = rocket.r_array[j]
    x4.append(x)
    y4.append(y)
    Rocket.set_data(x, y)
    RocketOrbit.set_data(x4, y4)
    rocket_speed.set_text("Speed: " + str(int(rocket.speed_array[j])))
    rocket_dist.set_text("Dist. to Jupiter: " + str(int(rocket.rR2_array[j])))
    return Earth, EarthOrbit, Jupiter, JupiterOrbit, Saturn, SaturnOrbit, Rocket, RocketOrbit, rocket_speed, rocket_dist,
anim = FuncAnimation(fig, animate, frames = framerate, interval = 1, blit = True)

plt.plot(x_earth, y_earth, lw=0.5, color='lightgray')
plt.plot(x_jupiter, y_jupiter, lw=0.5, color='lightgray')
plt.plot(x_saturn, y_saturn, lw=0.5, color='lightgray')
plt.plot(x_rocket, y_rocket, lw=0.5, color='lightgray')
plt.legend()
plt.show()