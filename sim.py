#!/usr/bin/python3
# see https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html

import numpy as np
from scipy.integrate import odeint, cumtrapz
from scipy.optimize import *
from scipy.interpolate import interp1d, barycentric_interpolate
import matplotlib.pyplot as plt

class Rocket:
    def __init__(self):
        self.motor = Motor()
        self.mass = 0 # kilograms
        self.cd = 0 # guess
        self.area = 0 # meters^2

    def load(self, fn):
        print("loading rocket", fn)
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
        print("loading motor", fn)
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
                    force = float(line.split(" ")[1])
                    time = float(line.split(" ")[0])

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
        s = [self.mass_casing + self.mass_prop, self.mass_prop]
        low = self.time[0]
        high = self.time[-1]
        if t < low:
            return s[0]
        elif t > high:
            return s[-1]

        # linear interpolation
        # this assumes the rate of propellant burn is constant
        # which it is not
        return (1 - t) * s[0] + t * s[1]

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
        self.atm = 101.3 # kPa
        self.p = 1.29

    def get_pressure(self, height, tmp):
        boltzmann = 1.38064852e-23
        return self.atm * np.exp((-1.0 * 4.8156e-26 * self.g * height) / (boltzmann * tmp))

def func(x, y):
    return interp1d(x, y, kind="cubic", fill_value="extrapolate")

# the second order diffeq
def model(v, t, rocket, world):
    x, dx = v

    thrust = rocket.motor.get_thrust(t)

    at = thrust / rocket.get_mass(t)
    ag = world.g
    #ad = (rocket.cd * rocket.area * world.get_pressure(x, 300) * np.power(dx, 2)) / rocket.get_mass(t)
    ad = np.sign(dx) * (0.5 * rocket.cd * rocket.area * world.p * np.power(dx, 2)) / rocket.get_mass(t)
    dvdt = [dx, at - ag - ad]

    return dvdt

if __name__ == "__main__":
    rocket = Rocket()
    world = World()
    rocket.load("./rocket.txt")

    v0 = [0, 0]
    t = np.linspace(0, 50, 1000)
    solution = odeint(model, v0, t, args=(rocket, world,))
    x = solution[:, 0]
    v = solution[:, 1]

    end = leastsq(func(t, x), 30)[0][0]
    v_terminal = func(t, v)(end)

    print("Terminal Velocity:", v_terminal)

    plt.plot(t, x, color="black", label="x(t)")
    plt.plot(t, v, color="red", label="v(t)")
    plt.xlabel("t")
    plt.xlim([0, end])
    plt.ylim(-500)
    plt.grid()

    plt.show()
