import yt
import numpy as np
import matplotlib.pyplot as plt
import sys

### define problem name
problem_name = 'CRSodShockTubeAMR_ZEUS_'


### define simulation output directory and filename base
output_dir_base = 'DD'
datafile_base = 'data'


### define output to be plotted
dumpid = '0001'

problem_name += dumpid

### construct input filename
filename = './' + output_dir_base + dumpid + '/' + datafile_base + dumpid
print("Plotting output file %s\n" % filename)


### some more filenames
exact_solution_filename = 'analytic.txt'
png_filename = './' + problem_name + '.png'



### load data
ds = yt.load(filename)


### define CR pressure
def _CRPressure(field,data):
    return (1.0/3.0)*data[('enzo','CREnergyDensity')] * yt.units.erg/yt.units.cm**3

ds.add_field('cosmic_ray_pressure', function=_CRPressure, units='erg/cm**3', sampling_type='cell')

### extract an ortho_ray (1D solution vector)
ray = ds.ortho_ray(0, [0.5, 0.5])


### read exact solution
exact = np.genfromtxt( exact_solution_filename, names=('x', 'density', 'velocity_x', 'pressure') )

### make plot

plt.figure(1, figsize=(8,7))


# density Plot
plt.subplot(221)
plt.plot((0,250,250,500),(1.0,1.0,0.2,0.2),color='k',linestyle='dotted')
plt.plot(ray['x'], ray['density'], color='k', marker='.', ls='none')
plt.plot(exact['x'], exact['density'],color='r', label='Analytic')
plt.legend()

#plt.axis([0,500,0,1.1])
plt.xlabel('Position')
plt.ylabel('Density')


## Velocity Plot
plt.subplot(222)
plt.plot(ray['x'],ray['velocity_x'], color='k', marker='.', ls='none')
plt.plot(exact['x'],exact['velocity_x'],color='r', label='Analytic')
plt.legend()

plt.axis([0,500,0,500])
plt.xlabel('Position')
plt.ylabel('Velocity')


# pressure Plot
plt.subplot(223)
plt.plot((0,250,250,500),(2.0,2.0,0,0),color='k',linestyle='dotted')
plt.plot(ray['x'],ray['pressure']/100000, 
         color='g', marker='.', linestyle='none', label='Gas')
plt.plot(ray['x'],ray['cosmic_ray_pressure']/100000,
         color='b', marker='.', linestyle='none', label='CRs')
plt.plot(ray['x'],ray['cosmic_ray_pressure']/100000 + ray['pressure']/100000,
         color='k', marker='.', ls='none', label='Total')
plt.plot(exact['x'],exact['pressure'],color='r', label='Analytic')
plt.legend()

plt.axis([0,500,0,2.2])
plt.xlabel('Position')
plt.ylabel('Pressure/100,000')


plt.savefig(png_filename)
