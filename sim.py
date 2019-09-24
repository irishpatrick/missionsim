#!/usr/bin/python3
# see https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html

import numpy as np
from scipy.integrate import *
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
        self.impulse = None
        self.impulse_total = 0.0
        self.mass_casing = 0
        self.mass_prop = 0
        self.mass_prop_cur = 0

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
                    self.mass_prop = float(value)
                    self.mass_prop_cur = self.mass_prop

                line = fp.readline()

        self.calculate()

    def calculate(self):
        self.impulse = cumtrapz(self.thrust, self.time, initial=0)
        #self.impulse_total = trapz(self.thrust, self.time)
        self.impulse_total = 2760.8899
        #print(self.impulse_total)

    def get_mass(self, t):
        if t < self.time[0]:
            return self.mass_casing + self.mass_prop
        elif t > self.time[-1]:
            return self.mass_casing

        J = func(self.time, self.impulse)
        prop = self.mass_prop - ((J(t) / self.impulse_total) * self.mass_prop)
        return self.mass_casing + prop

    def get_thrust(self, t):
        low = self.time[0]
        high = self.time[-1]
        if t < low:
            return self.thrust[0]
        if t > high:
            return 0

        return func(self.time, self.thrust)(t)

class World:
    def __init__(self):
        self.g = 9.8 # m/s/s
        self.atm = 101300 # Pa
        self.p = 1.225 # density
        self.mm = 0.0289644
        #self.base = 472.1352
        self.base = 0
        self.R = 8.3144598

    def get_kelvin(self, temp):
        #return (temp - ((height / 304.8) * 2)) + 273.15
        return temp + 273.15

    def get_density(self, height, temp):
        return self.p * np.exp((-self.g * self.mm * (height - self.base)) / (self.R * self.get_kelvin(temp)))

def func(x, y):
    return interp1d(x, y, kind="quadratic", fill_value="extrapolate")

def plot(f, x):
    out = []
    for i in range(len(x)):
        out.append(f(x[i]))

    return out
    
# the second order diffeq
def model(v, t, rocket, world):
    x, dx = v

    thrust = rocket.motor.get_thrust(t)

    at = thrust / rocket.get_mass(t)
    ag = world.g
    ad = np.sign(dx) * (0.5 * rocket.cd * rocket.area * world.get_density(x, 20) * np.power(dx, 2)) / rocket.get_mass(t)
    dvdt = [dx, at - ag - ad]
    return dvdt

if __name__ == "__main__":
    rocket = Rocket()
    world = World()
    rocket.load("./rocket.txt")

    v0 = [0, 0]
    t = np.linspace(0, 50, 1000)
    solution = odeint(model, v0, t, args=(rocket, world))
    x = solution[:, 0]
    v = solution[:, 1]

    end = leastsq(func(t, x), 30)[0][0]

    plt.plot(t, x, color="black", label="x(t)")
    plt.plot(t, v, color="blue", label="v(t)")
    #plt.plot(t, plot(, t), color="red", label="a(t)")
    plt.xlabel("t")
    plt.xlim([0, end])
    plt.ylim(-500)
    plt.grid()

    """
    fig, ax1 = plt.subplots()
    
    ax1.set_xlabel("time")
    ax1.set_ylabel("x(t)", color="black")
    ax1.plot(t, x, color="black")
    ax1.tick_params(axis="y", labelcolor="black")
    ax1.set_xlim([0, end])
    ax1.set_ylim(0)

    ax2 = ax1.twinx()

    ax2.set_ylabel("v(t)", color="red")
    ax2.plot(t, v, color="red")
    ax2.tick_params(axis="y", labelcolor="red")

    fig.tight_layout()
    """
    plt.show()
