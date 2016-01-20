###############################################################################
 # This file is part of SWIFT.
 # Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 #                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 # 
 # This program is free software: you can redistribute it and/or modify
 # it under the terms of the GNU Lesser General Public License as published
 # by the Free Software Foundation, either version 3 of the License, or
 # (at your option) any later version.
 # 
 # This program is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU General Public License for more details.
 # 
 # You should have received a copy of the GNU Lesser General Public License
 # along with this program.  If not, see <http://www.gnu.org/licenses/>.
 # 
 ##############################################################################

import h5py
import random
from numpy import *

# Generates a swift IC file for the Gresho Vortex in a periodic box

# Parameters
periodic= 1      # 1 For periodic box
factor = 3
boxSize = [ 1.0 , 1.0, 1.0/factor ]
L = 120           # Number of particles along one axis
gamma = 5./3.     # Gas adiabatic index
rho = 1           # Gas density
P0 = 0.           # Constant additional pressure (should have no impact on the dynamics)
fileName = "greshoVortex.hdf5" 
vol = boxSize[0] * boxSize[1] * boxSize[2]


#---------------------------------------------------

numPart = L**3 / factor
mass = boxSize[0]*boxSize[1]*boxSize[2] * rho / numPart

#Generate particles
coords = zeros((numPart, 3))
v      = zeros((numPart, 3))
m      = zeros((numPart, 1))
h      = zeros((numPart, 1))
u      = zeros((numPart, 1))
ids    = zeros((numPart, 1), dtype='L')

partId=0
for i in range(L):
    for j in range(L):
        for k in range(L/factor):
            index = i*L*L/factor + j*L/factor + k
            x = i * boxSize[0] / L + boxSize[0] / (2*L)
            y = j * boxSize[0] / L + boxSize[0] / (2*L)
            z = k * boxSize[0] / L + boxSize[0] / (2*L)
            r2 = (x - boxSize[0] / 2)**2 + (y - boxSize[1] / 2)**2
            r = sqrt(r2)
            coords[index,0] = x
            coords[index,1] = y
            coords[index,2] = z
            v_phi = 0.
            if r < 0.2:
                v_phi = 5.*r
            elif r < 0.4:
                v_phi = 2. - 5.*r
            else:
                v_phi = 0.
            v[index,0] = -v_phi * (y - boxSize[0] / 2) / r
            v[index,1] =  v_phi * (x - boxSize[0] / 2) / r
            v[index,2] = 0.
            m[index] = mass
            h[index] = 2.251 * boxSize[0] / L
            P = P0
            if r < 0.2:
                P = P + 5. + 12.5*r2
            elif r < 0.4:
                P = P + 9. + 12.5*r2 - 20.*r + 4.*log(r/0.2)
            else:
                P = P + 3. + 4.*log(2.)
            u[index] = P / ((gamma - 1.)*rho)
            ids[index] = partId
            partId = partId + 1



#File
file = h5py.File(fileName, 'w')

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = boxSize
grp.attrs["NumPart_Total"] =  [numPart, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]

#Runtime parameters
grp = file.create_group("/RuntimePars")
grp.attrs["PeriodicBoundariesOn"] = periodic

#Particle group
grp = file.create_group("/PartType0")
ds = grp.create_dataset('Coordinates', (numPart, 3), 'd')
ds[()] = coords
ds = grp.create_dataset('Velocities', (numPart, 3), 'f')
ds[()] = v
ds = grp.create_dataset('Masses', (numPart,1), 'f')
ds[()] = m
ds = grp.create_dataset('SmoothingLength', (numPart,1), 'f')
ds[()] = h
ds = grp.create_dataset('InternalEnergy', (numPart,1), 'f')
ds[()] = u
ds = grp.create_dataset('ParticleIDs', (numPart,1), 'L')
ds[()] = ids[:]

file.close()


