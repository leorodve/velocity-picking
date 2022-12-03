import numpy as np 
from scipy.interpolate import CubicSpline

def nmo_correction(cmp, dt, offsets, velocities): 
    nmo = np.zeros_like(cmp) 
    nsamples = cmp.shape[0] 
    times = np.arange(0, nsamples*dt, dt)
    t = np.sqrt(np.reshape(times**2, (times.size, 1)) + np.multiply.outer(1/(velocities**2), offsets**2))
    for i in range(0, len(offsets)):
        interpolator = CubicSpline(times, cmp[:,i])
        amplitude = interpolator(t[:,i])
        if amplitude is not None:
            nmo[:, i] = amplitude 
    return nmo

def calc_velocity_spectrum(cmp, sample_rate, offsets, velocities): 
    cvs = np.zeros_like(cmp)
    nsamples = cvs.shape[0] 
    velocity_spectrum = np.zeros([nsamples, len(velocities)])
    times = np.arange(0, nsamples*sample_rate, sample_rate)
    i = 0
    while i <= 10:
        t = np.sqrt(np.reshape(times**2, (times.size, 1)) + np.multiply.outer(1/(velocities[i]**2), offsets**2))
        for j in range(0, len(offsets)):
            interpolator = CubicSpline(times, cmp[:,j])
            amplitude = interpolator(t[:,j])
            if amplitude is not None:
                cvs[:, j] = amplitude 
        velocity_spectrum[:,i] = np.sum(cvs, axis=1)
        velocity_spectrum.clip(0, np.max(velocity_spectrum), out=velocity_spectrum)
        i += 1
    return velocity_spectrum

def create_velocity_function(min_vel, max_vel, ns, dt, vel_picks):
    vel_function = np.zeros(ns)
    velocity_picks = np.array(vel_picks, dtype=float)
    velocity_picks[:,0] =  min_vel + (velocity_picks[:,0] - 548)*(max_vel-min_vel)/336
    velocity_picks[:,1] = (velocity_picks[:,1] - 48)*(ns*dt)/302
    times = np.arange(0, ns*dt, dt)
    interpolator = CubicSpline(velocity_picks[:,1], velocity_picks[:,0])
    vel_function[:] = interpolator(times)
    return vel_function, velocity_picks

def regression_coordinates(velocity_prediction, minimum_velocity, 
                           maximum_velocity, ns):
    coord = np.argwhere(velocity_prediction >= minimum_velocity)
    ml_picks_coord = np.zeros((len(coord),2))
    ml_vel_coord = np.round(548 + (velocity_prediction[coord] - minimum_velocity) * 336 / (maximum_velocity - minimum_velocity))
    ml_time_coord = np.round(48 + (coord * 302 / ns))
    ml_picks_coord[:,:1] = ml_vel_coord
    ml_picks_coord[:,1:] = ml_time_coord

    return ml_picks_coord

def stretch_mute(velocity, ns, dt, stretch_mute):
    time = np.arange(0, ns*dt, dt)
    Offset_max = np.multiply(velocity, time) * pow(np.square(1 
                            + stretch_mute/100) - 1, 0.5)
    return Offset_max

def apply_stretch_mute(mute_function, nmo, fold):
    for i, x in enumerate (mute_function):
        j = int(np.floor(x/fold))
        if j > fold:
            j = fold
        nmo[i, j:] = 0
    return nmo
