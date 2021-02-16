"""
The nmo_correction function was created by Uieda, L. and published
on the article: Step-by-step NMO correction, The Leading Edge,
doi:10.1190/tle36020179.1 and can be found on
https://github.com/pinga-lab/nmo-tutorial
"""

import numpy as np 
from scipy.interpolate import CubicSpline

def nmo_correction(cmp, dt, offsets, velocities): 
    nmo = np.zeros_like(cmp) 
    nsamples = cmp.shape[0] 
    times = np.arange(0, nsamples*dt, dt) 
    for i, t0 in enumerate(times): 
        for j, x in enumerate(offsets): 
            t = reflection_time(t0, x, velocities[i])
            amplitude = sample_trace(cmp[:, j], t, dt)
            if amplitude is not None: 
                nmo[i, j] = amplitude 
    return nmo

def reflection_time(t0, x, vnmo): 
    t = np.sqrt(t0**2 + x**2/vnmo**2) 
    return t

def sample_trace(trace, time, dt): 
    before = int(np.floor(time/dt)) 
    N = trace.size 
    a = before - 1
    b = before + 3
    samples = np.arange(a, b) 
    if any(samples < 0) or any(samples >= N): 
        amplitude = None 
    else: 
        times = dt*samples 
        amps = trace[samples] 
        interpolator = CubicSpline(times, amps) 
        amplitude = interpolator(time) 
    return amplitude

def calc_velocity_spectrum(cmp, sample_rate, offsets, velocities): 
    cvs = np.zeros_like(cmp)
    nsamples = cvs.shape[0] 
    velocity_spectrum = np.zeros([nsamples, len(velocities)])
    times = np.arange(0, nsamples*sample_rate, sample_rate)
    i = 0
    while i <= 10:
        for a, t0 in enumerate(times): 
            for j, x in enumerate(offsets): 
                t = reflection_time(t0, x, velocities[i])
                amplitude = sample_trace(cmp[:, j], t, sample_rate)
                if amplitude is not None: 
                    cvs[a, j] = amplitude
            row_value = sum(cvs[a, :])
            if row_value < 0:
                velocity_spectrum[a, i] = 0
            else:
                velocity_spectrum[a, i] = row_value
        i += 1
    return velocity_spectrum

def create_velocity_function(min_vel, max_vel, ns, dt, vel_picks):
    velocity_coordinates = np.arange(min_vel, max_vel, (max_vel-min_vel)/336)
    times_coordinates = np.arange(0, ns*dt, ns*dt/302)
    vel_function = np.zeros(ns)
    velocity_picks = np.zeros((len(vel_picks), 2))
    for id1, (value_x, value_y) in enumerate(vel_picks):
        for id2, vel_value in enumerate(velocity_coordinates, 548):
            for id3, time_value in enumerate(times_coordinates, 48):
                if value_x == id2 and value_y == id3:
                    velocity_picks[id1, :] = [vel_value, time_value]
                    vel_function[int(round(time_value/dt))] = vel_value
    id5 = 0
    for id4, value in enumerate(vel_function):
        if id4 <= round(velocity_picks[0,1]/dt):
            vel_function[id4] = velocity_picks[0,0]
        elif id4 >= round(velocity_picks[len(velocity_picks)-1,1]/dt):
            vel_function[id4] = velocity_picks[len(velocity_picks)-1,0]
        else:
            if id4 > round(velocity_picks[id5 + 1, 1]/dt):
                id5 += 1
            a1 = velocity_picks[id5, 0]
            b1 = round(velocity_picks[id5, 1]/dt)
            a2 = velocity_picks[id5 + 1, 0]
            b2 = round(velocity_picks[id5 + 1, 1]/dt)
            v = [a1, a2]
            t = [b1, b2]
            interpolator = CubicSpline(t, v) 
            vel_function[id4] = interpolator(id4)
    return vel_function, velocity_picks
def stretch_mute(velocity, ns, dt, stretch_mute):
    time = np.arange(0, ns*dt, dt)
    Offset_max = np.multiply(velocity, time) * pow(np.square(1 + stretch_mute/100) - 1, 0.5)
    return Offset_max

def apply_stretch_mute(mute_function, nmo, fold):
    for i, x in enumerate (mute_function):
        j = int(np.floor(x/fold))
        if j > fold:
            j = fold
        nmo[i, j:] = 0
    return nmo
