from cProfile import label
import os
import numpy as np
import matplotlib.pyplot as plt

np.set_printoptions(suppress = True)
filepath = os.path.abspath('.')

# imu_circle   imu_spline
error_integration = []
error_integration = np.loadtxt(filepath + '/data_gen/integration_error.txt', usecols = (0,1))

error_preintegration = []
error_preintegration = np.loadtxt(filepath + '/data_gen/preintegration_error.txt', usecols = (0,1))

x_axis = range(len(error_integration))

r_error = plt.subplot(1,2,1)
plt.plot(x_axis,error_integration[:,0])
plt.plot(x_axis,error_preintegration[:,0])
plt.ylabel('r_error')

t_error = plt.subplot(1,2,2)
plt.plot(x_axis,error_integration[:,1])
plt.plot(x_axis,error_preintegration[:,1])
plt.ylabel('t_error')

plt.show()
