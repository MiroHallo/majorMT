#! /usr/bin/env python

############################################################
# 
# Plot trinity of beach-balls from Major and Minor DC decomposition
#
############################################################
#
# Hallo,M., Asano,K., Gallovic,F. (2017): Bayesian inference and interpretation
#      of centroid moment tensors of the 2016 Kumamoto earthquake sequence, 
#      Kyushu, Japan, Earth, Planets and Space, 69:134.
# Hallo,M., Oprsal,I., Asano,K., Gallovic,F. (2019): Seismotectonics of the 2018
#      Northern Osaka M6.1 earthquake and its aftershocks: joint movements
#      on strike-slip and reverse faults in inland Japan, Earth,
#      Planets and Space, 71:34.
#
# Code author: Miroslav Hallo
# Charles University in Prague, Faculty of Mathematics and Physics
# Web: http://geo.mff.cuni.cz/~hallo/
# E-mail: hallo@karel.troja.mff.cuni.cz
# Revision 2/2017: The first version of the function.
# Revision 12/2018: Enhanced version.
#
# Copyright (C) 2017,2018  Miroslav Hallo
#
# This program is published under the GNU General Public License (GNU GPL).
#
# This program is free software: you can modify it and/or redistribute it
# or any derivative version under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.
#
# This code is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY. We would like to kindly ask you to acknowledge the authors
# and don't remove their names from the code.
#
# You should have received copy of the GNU General Public License along
# with this program. If not, see <http://www.gnu.org/licenses/>.
#
############################################################

import matplotlib.pyplot as plt
import numpy as np
from obspy.imaging.mopad_wrapper import beach as beach2

############################################################

# INPUT

# Moment tensors for Osaka M6.1 mainshock
# [M11, M22, M33, M12, M13, M23]
mt  = np.array([[1.1067,       1.5367,     -2.6433,       0.26,       0.08,       -0.720]])  # Input non-DC MT (Harvard)
mt1 = np.array([[0.21921,      1.4434,     -1.6626,    0.58336,   -0.035697,    -0.52161]]) # Major DC MT (Harvard)
mt2 = np.array([[0.88745,    0.093313,    -0.98077,   -0.32336,      0.1157,    -0.19839]]) # Minor DC MT (Harvard)
Ra1 = 0.632 # Scalar seismic moment ratio of Major DC MT
Ra2 = 0.368 # Scalar seismic moment ratio of Minor DC MT

# END INPUT

############################################################
############################################################

outfile = "input_MT.png"
facecolor='red'

fig = plt.figure(figsize=(5,5))
ax = plt.axes()
plt.axis('off')
ax.axes.get_xaxis().set_visible(False)
ax.axes.get_yaxis().set_visible(False)
lw=2
plt.xlim(-100-lw/2, 100+lw/2)
plt.ylim(-100-lw/2, 100+lw/2)

Na=200
full = beach2(mt[0,0:6], xy=(0, 0), linewidth=lw, facecolor=facecolor, edgecolor='black', zorder=1, width=Na)
ax.add_collection(full)

plt.savefig(outfile, bbox_inches='tight', pad_inches=0)
plt.clf()

#------------------
RaN = Ra1 + Ra2
Na1 = Ra1/RaN * 200
Na2 = Ra2/RaN * 200

outfile = "trinity_MT.png"
  
fig = plt.figure(figsize=(15,5))
ax = plt.axes()
plt.axis('off')
ax.axes.get_xaxis().set_visible(False)
ax.axes.get_yaxis().set_visible(False)
lw=2
plt.xlim(-100-lw/2, 500+lw/2)
plt.ylim(-100-lw/2, 100+lw/2)

full = beach2(mt1[0,0:6], xy=(0, 0), linewidth=lw, facecolor='black', edgecolor='black', zorder=1, width=Na1)
ax.add_collection(full)

full = beach2(mt2[0,0:6], xy=(200, 0), linewidth=lw, facecolor='black', edgecolor='black', zorder=1, width=Na2)
ax.add_collection(full)

full = beach2(mt1[0,0:6]+mt2[0,0:6], xy=(400, 0), linewidth=lw, facecolor=facecolor, edgecolor='black', zorder=1)
ax.add_collection(full)

plt.savefig(outfile, bbox_inches='tight', pad_inches=0)
plt.clf()
plt.close()


