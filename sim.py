#!/usr/bin/python3
# see https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html

import numpy as np
from scipy.integrate import odeint, simps, trapz, cumtrapz, quad
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

class Rocket:
    def __init__(self):
        self.motor = Motor()
        self.mass = 10.0 # kilograms
        self.cd = 1.2 # guess
        self.area = 0.018 # meters^2

    def load(self, fn):
        with open(fn, "r") as fp:
            line = fp.readline()
            while line:
                line = line[:-1]
                if line.startswith("mass="):
                    self.mass = float(line.split("=")[1])
                elif line.startswith("area="):
                    self.area = float(line.split("=")[1])
                elif line.startswith("cd="):
                    self.cd = float(line.split("=")[1])
                elif line.startswith("motor="):
                    self.motor.load("./motors/" + line.split("=")[1])
                
                line = fp.readline()

    def get_mass(self, t):
        return self.mass + self.motor.get_mass(t)

class Motor:
    def __init__(self):
        self.thrust = []
        self.time = []
        self.mass_casing = 0
        self.mass_prop = 0

    def load(self, fn):
        print("loading motor ", fn)
        with open(fn, "r") as fp:
            line = fp.readline()
            read_thrust_data = False

            while line:
                line = line[:-1]
                parts = line.split("=")
                name = parts[0]
                value = ""
                if len(parts) > 1:
                    value = parts[1]

                if line.startswith("end thrust"):
                    read_thrust_data = False

                if read_thrust_data and line != "":
                    force = float(line.split(" ")[0])
                    time = float(line.split(" ")[1])

                    self.time.append(time)
                    self.thrust.append(force)

                if line.startswith("begin thrust"):
                    read_thrust_data = True

                if name == "mass_casing":
                    self.mass_casing = float(value)
                elif name == "mass_prop":
                    self.mas_prop = float(value)

                line = fp.readline()

    def get_mass(self, t):
        return 2 

    def get_thrust(self, t):
        low = self.time[0]
        high = self.time[-1]
        if t < low:
            return 0
        if t > high:
            return 0

        return interp1d(self.time, self.thrust, kind="cubic")(t)

class World:
    def __init__(self):
        self.g = 9.81 # m/s/s
        self.atm = 0.133322 # kPa

    def get_pressure(self, height, tmp):
        boltzmann = 1.38064852e-23
        return self.atm * np.exp((-1.0 * 4.81069412e-26 * self.g * height) / (boltzmann * tmp))

def acceleration(v, t, rocket, world):
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
    t = np.linspace(0, 34)
    v = odeint(acceleration, v0, t, args=(rocket, world,))

    x = cumtrapz(np.swapaxes(v,0,1)[0], t, initial=0)
    
    plt.plot(t, v, "blue", t, x, "black")
    plt.xlabel("time")
    plt.ylabel("v(t)")

    plt.show()
