import os
import numpy as np
import matplotlib.pyplot as plt

np.set_printoptions(suppress = True)
filepath = os.path.abspath('.')
tx_index = 5

position = []
position = np.loadtxt(filepath + '/data_gen/imu_pose_noise.txt', usecols = (5,6,7))

position1 = []
position1 = np.loadtxt(filepath + '/data_gen/vision_pose.txt',usecols=(4,5,6))

position2 = []
position2 = np.loadtxt(filepath + '/data_gen/integration_pose.txt',usecols=(4,5,6))

position3 = []
position3 = np.loadtxt(filepath + '/data_gen/preintegration_pose.txt',usecols=(4,5,6))

fig = plt.figure()
ax = fig.gca(projection='3d')

ax.plot(position[:,0], position[:,1], position[:,2], label='gt')
ax.plot(position1[:,0], position1[:,1], position1[:,2], label='vision')
ax.plot(position2[:,0], position2[:,1], position2[:,2], label='integration')
ax.plot(position3[:,0], position3[:,1], position3[:,2], label='preintegration')
ax.plot([position[0,0]], [position[0,1]], [position[0,2]], 'r.', label='start')

ax.legend()
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()
