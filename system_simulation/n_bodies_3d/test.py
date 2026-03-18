import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


fig = plt.figure()
ax = fig.add_subplot(1, 2, 1, projection='3d')
ax.scatter([0, 1, 2], [2, 3, 4], [5, 50, 2])
bx = fig.add_subplot(1, 2, 2, projection=None)



plt.show()