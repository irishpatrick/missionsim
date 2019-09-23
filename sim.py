#!/usr/bin/python3

import numpy as np
from scipy.integrate import odeint, simps, trapz, cumtrapz, quad
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

class Rocket:
    def __init__(self):
        self.motor = Motor()
        self.mass = 10.0 # kilograms
        self.cd = 0.7 # guess
        self.area = 0.018 # meters^2

    def load(self, fn):
        with open(fn, "r") as fp:
            line = fp.readline()
            while line:
                if line.startswith("mass="):
                    self.mass = float(line.split("=")[1][:-1])
                elif line.startswith("area="):
                    self.area = float(line.split("=")[1][:-1])
                elif line.startswith("cd="):
                    self.cd = float(line.split("=")[1][:-1])
                elif line.startswith("motor="):
                    self.motor.load("./motors/" + line.split("=")[1][:-1])
                
                line = fp.readline()

    def get_mass(self, t):
        return self.mass + self.motor.get_mass(t)

class Motor:
    def __init__(self):
        self.thrust = []
        self.mass = []

    def load(self, fn):
        with open(fn, "r") as fp:
            line = fp.readline()

            while line:

                line = fp.readline()

    def get_mass(self, t):
        return 2 

    def get_thrust(self, t):
        return 200

class World:
    def __init__(self):
        self.g = 9.81 # m/s/s
        self.atm = 0.133322 # kPa

    def get_pressure(self, height, tmp):
        boltzmann = 1.38064852e-23
        return self.atm * np.exp((-1.0 * 4.81069412e-26 * self.g * height) / (boltzmann * tmp))

def acceleration(v, t, rocket, world):
    thrust = 0
    if (t < 2):
        thrust = rocket.motor.get_thrust(t)

    at = thrust / rocket.get_mass(t)
    ag = world.g
    ad = (rocket.cd * rocket.area * world.atm * np.power(v, 2)) / rocket.get_mass(t)

    dvdt = at - (ag + ad)

    return dvdt

if __name__ == "__main__":
    rocket = Rocket()
    world = World()
    rocket.load("./rocket.txt")

    v0 = 0
    t = np.linspace(0, 5, num=1000)
    v = odeint(acceleration, v0, t, args=(rocket, world,))

    x = cumtrapz(np.swapaxes(v,0,1)[0], t, initial=0)
    
    plt.plot(t, v, "black", t, x, "red")
    plt.xlabel("time")
    plt.ylabel("v(t)")

    plt.show()
