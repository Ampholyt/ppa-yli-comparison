import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

soa =np.array( [[0,0,0,1,2,0], 
                [0,0,0,1,1,0],
                [0,0,0,2,1,0]])

# Glucose, Oxygen, Growth rate
# (2.2599, 3.2739, 0.2261) # EFM 1
# (0.0218, 0.0315, 0.0022) # EFM 2 
# (0.0515, 0.0746, 0.0052) # EFM 3
# (0.0969, 0.1404, 0.0097) # EFM 4
efm1 = [2.2599, 3.2739, 0.2261]
efm2 = [0.0218, 0.0315, 0.0022]
efm3 = [0.0515, 0.0746, 0.0052]
efm4 = [0.0969, 0.1404, 0.0097]
start = [0, 0, 0]
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.quiver(start[0], start[1], start[2], efm1[0], efm1[1], efm1[2], color='b')
ax.quiver(start[0], start[1], start[2], efm2[0], efm2[1], efm2[2], color='r')
ax.quiver(start[0], start[1], start[2], efm3[0], efm3[1], efm3[2], color='g')
ax.quiver(start[0], start[1], start[2], efm4[0], efm4[1], efm4[2], color='y')
ax.set_xlim([0, 4])
ax.set_ylim([0, 4])
ax.set_zlim([0, 4])
# Set labels for x, y, and z axes
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_ylim([0, 0.7])
ax.set_ylim([0, 0.3])
ax.set_zlim([0, 0.1])
ax.view_init(45, 215)
plt.show()

# Vector origin location
# EFM 1, 2, 3, 4, 5, 6, 7, 8
# X = [0] # 
# Y = [0]
# Z = [0]
  
# # Directional vectors
# U = [3]  
# V = [1]  
  
# # Creating plot
# plt.quiver(X, Y, U, V, color='b', units='xy', scale=1)
# plt.title('Single Vector')
  
# # x-lim and y-lim
# plt.xlim(-1, 3)
# plt.ylim(-1, 3)
  
# # Show plot with grid
# plt.grid()
# plt.show()