import numpy as np
import numpy.random as rand
import random
from matplotlib import collections  as mc
import math
import matplotlib.pyplot as plt

"""An implementation of a model for a mycelial network.
Based off the one formulated in "The Development of Fungal Networks in Complex
Environments" by Graeme P. Boswell, et al. (2007) https://link.springer.com/content/pdf/10.1007/s11538-005-9056-6.pdf

This model has 5 key member arrays: 
   -- active [0 .. 63]
      -- Tracks the active hyphae in a cell
      -- as a binary number 
   -- inactive [0 .. 63]
      -- Tracks the inactive hyphae in a cell
   -- internal (0, 1)
      -- Holds internal nutrient values for a cell
   -- external (0, 1)
      -- holds external nutrient values for a cell
   -- tips [0 1]
      -- says whether or not a cell holds a tip

And has 5 distinct functions:
   -- grow
      -- causes tips to move according to a probability density
   -- branch
      -- causes tips to branch according to a probability
   -- maintain
      -- get rid of unused hyphae, merge colliding tips
   -- uptake
      -- pick up nutrients from the environment
   -- translocate
      -- transfer nutrients throughout the network
"""

sqrt3half = math.sqrt(3)/2

# vectors at angles 0, π/3, ... , 5π/3 with magnitude 0.5
direcs = (np.array([0,1])/2, np.array([-sqrt3half,0.5])/2,np.array([-sqrt3half,-0.5])/2,
np.array([0,-1])/2,np.array([sqrt3half,-0.5])/2,
np.array([sqrt3half,0.5])/2)

# cartesian displacement vectors for the adjacent hexagons corresponding to the above directions
#           ((even rows), (odd rows))
hexDirecs = ((np.array([0,1]),np.array([-1,0]), np.array([-1,-1]),np.array([0,-1]),np.array([1,-1]),np.array([1,0])),
(np.array([0,1]), np.array([-1,1]), np.array([-1,0]), np.array([0,-1]), np.array([1,0]), np.array([1,1])))

class Model:

   def __init__(self,n=30):
      self.size = n
      self.lines = []
      self.makeLines()
      self.dx = 1
      # self.Dp = 1e4
      # self.v = 1e5
      # self.b = 1e6
      # self.c1 = 10
      # self.c2 = 1e-7
      # self.c3 = 1e2
      # self.c4 = 1e-8
      # self.c5 = 1e-9
      self.Dp = 0.3
      self.Di = 0.1
      self.Da = 0.1
      self.v = 1
      self.b = 2
      self.c1 = 0.03
      self.c2 = 0.01
      self.c3 = 0.05
      self.c4 = 0.001
      self.c5 = 0.001
      self.zero = np.zeros((n+4,n+4))
      self.one = np.ones((n+4,n+4))
      self.external_init = 1.5*math.sqrt(3)*1e-5#1e-11
      self.omega = 1e-13
      self.precision = np.float16
      self.active = np.zeros((n+4,n+4),dtype=np.int8)
      self.active[1+(n//2),1+(n//2)] = 1
      self.inactive = np.zeros((n+4,n+4),dtype=np.int8)
      self.tips = np.zeros((n+4,n+4),dtype=np.int8)
      self.tips[(n//2)+1,(n//2)+1] = 1
      self.orientation = np.where(self.tips,np.mod(np.log2(self.active)+3,6),-1)
      print(self.orientation)
      self.inCell = np.zeros((n,n),dtype=np.int8)
      self.internal = np.ones((n+4,n+4),dtype=self.precision)
      self.external = np.ones((n+4,n+4),dtype=self.precision)*self.external_init
      self.dt = 0.5

   # from https://stackoverflow.com/questions/8871204/count-number-of-1s-in-binary-representation
   def count_ones(self, n):
      count=np.zeros(n.shape)
      while np.sum(self.one[(remaining:=np.where(n!=0))])>0:
         n[remaining] = n[remaining]&(n[remaining]-1)
         count[remaining]+=1
      return count

   def update(self):
      model.grow()
      model.maintain()
      model.translocate()

   def grow(self):

      # calculate movement probabilities from internal substrate concentrations
      movementDensity = np.array([diffusion:=self.Dp*self.internal*self.dt*self.dx**(-2), 2*diffusion+(convection:=self.v*self.internal*self.dt*self.dx**(-1)), 3*diffusion+convection])
      
      # generate random numbers
      probs = rand.rand(self.size+4,self.size+4)
      
      # figure out which cells are moving in each direction
      # then move them

      # hyphae which randomly turn to the right
      eta=np.where(np.logical_and(probs<movementDensity[0,:,:],self.tips==1))
      # lambd = np.where(lambda_cond:=np.logical_and(np.logical_and(movementDensity[0,:,:]<=probs,probs<movementDensity[1,:,:]),self.tips))
      
      # hyphae which randomly turn to the left
      mu = np.where(np.logical_and(np.logical_and(movementDensity[1,:,:]<=probs,probs<movementDensity[2,:,:]),self.tips==1))
      
      # hyphae which don't move
      self.orientation = np.where(np.logical_and(probs>=movementDensity[2,:,:],self.tips==1), -1, self.orientation)

      # change cell orientation where we are turning
      self.orientation[eta] -= 1
      self.orientation[mu] += 1

      # for each hypha orientation
      for theta in range(6):

         # mask out all other cells
         tips_masked = np.ma.masked_where(self.orientation!=theta, self.tips)

         # print("tips @ 72")
         # print(self.tips)
         
         # for each cell if there's a tip in the current orientation grow in that direction
         for j in range(2, self.size+2):
            for i in range(2, self.size+2):

               # tip in cell i,j at orientation theta
               if tips_masked[i,j]==1:

                  # add a tip to the adjacent cell in the theta direction
                  self.tips[i+hexDirecs[i&1][theta][0],j+hexDirecs[i&1][theta][1]] += 1

                  self.internal[i,j] -= self.c2 * self.dx # use up nutrient

         # self.tips[2 + hexDirecs[i][0] : self.size + 2 + hexDirecs[theta][0],
         # 2 + hexDirecs[theta][1]:  self.size + 2 + hexDirecs[theta][1]] += tips_masked[2:-2,2:-2]

         # tips that just grew are no longer where they used to be
         self.tips = np.where(tips_masked==1, self.tips-1, self.tips)

         # add a hypha in the theta direction
         self.active[np.where(tips_masked==1)] |= (2**theta)
         
         # branch
         self.branch(tips_masked)

         # for each cell that grew add a hypha in the cell it grew into in the theta direction
         for j in range(2, self.size+2):
            for i in range(2, self.size+2):
               if tips_masked[i,j] == 1:
                  self.active[i+hexDirecs[i&1][theta][0], j+hexDirecs[i&1][theta][1]] |= (othernum:=(2**((theta+3)%6)))
      

   def maintain(self):
      
      # more than one tip should be resolved to one tip
      self.tips = np.where(cond:=self.tips>1, 1, self.tips)

      # where there are tips get the orientation of those tips, no tip -> -1
      self.orientation = np.where(self.tips,np.mod(np.log2(self.active)+3,6),-1)

      # it costs nutrients to be a hypha! deduct for each in cell.
      # self.internal[cond] -= self.c5*self.count_ones(self.active[(cond:=np.where(self.active!=0))])


      self.internal[(new_inactive:=np.where(self.internal<self.omega))] = self.omega

      self.inactive[new_inactive] = self.active[new_inactive]

      self.active[new_inactive] = 0
      # # anastosmosis
      # # for i in range(6):
      # #    self.tips = np.where(np.any(self.active[np.where(self.orientation==i)]==np.setdiff1d(2**np.arange(7),2**i), axis=0),0,self.tips)


   def branch(self, tips_masked):
      options = []
      probs = rand.rand(self.size+4,self.size+4)

   
      # for each cell
      for j in range(2, self.size+2):
         for i in range(2, self.size+2):

            # if our random variable is less than b Si dt and there is a tip
            if probs[i,j]<self.b*self.internal[i,j]*self.dt and (tips_masked[i,j]==1):
               print("i: %d, j = %d" % (i,j))

               # enumerate possible branch directions
               for option in (-1,0,1):
                  num = int(2**((self.orientation[i,j]+option)%6))
                  if (self.active[i,j] | num) != num:
                     options.append(option)
               if len(options)==0: # nowhere to branch
                  continue
               print(options)

               # get growth angle
               theta = int((self.orientation[i,j]+random.choice(options)) % 6)

               # make new tip in adjacent hex in that direction
               self.tips[i+hexDirecs[i&1][theta][0],j+hexDirecs[i&1][theta][1]] += 1

               # make hypha in this cell and adjacent cell in that direction
               self.active[i,j] |= (2**theta)
               self.active[i+hexDirecs[i&1][theta][0],j+hexDirecs[i&1][theta][1]] |= (2**((theta+3)%6))
   
   def uptake(self):
      # absorb amount of nutrients proportional to length of hyphae in cell and amount of external nutrients
      self.internal *= (1+self.c1*(num:=(self.count_ones(self.active)+self.count_ones(self.inactive))*self.external*self.dt))
      
      # subtract the nutrients we absorb from the environment
      self.external -= self.c3 * self.internal * num


   def translocate(self):

      # for each direction of adjacent cells
      for theta in range(6):
         for j in range(2, self.size+2):
            for i in range(2, self.size+2):

               # if there's hyphae connecting this cell to the adjacent cell in direction theta
               if (self.active[i,j] | 2**theta==2**theta) and (self.active[i+hexDirecs[i&1][theta][0],j+hexDirecs[i&1][theta][1]] | 2**((theta+3)%6)==2**((theta+3)%6)):

                  # passive translocation
                  self.internal[i,j]+=self.Di*(self.internal[i+hexDirecs[i&1][theta][0],j+hexDirecs[i&1][theta][1]]-self.internal[i,j])*self.dx**(-2)

                  # active translocation
                  if self.tips[i,j]-self.tips[i+hexDirecs[i&1][theta][0],j+hexDirecs[i&1][theta][1]]>0:
                     m = self.internal[i+hexDirecs[i&1][theta][0],j+hexDirecs[i&1][theta][1]] * self.Da * self.dt * self.dx**(-2)
                     self.internal[i+hexDirecs[i&1][theta][0],j+hexDirecs[i&1][theta][1]] -= m * self.c4
                  elif self.tips[i,j]-self.tips[i+hexDirecs[i&1][theta][0],j+hexDirecs[i&1][theta][1]]<0:
                     m = self.internal[i,j] * self.Da * self.dt * self.dx**(-2)
                     self.internal[i,j] -= m * self.c4
                  else:
                     m = 0
                  self.internal[i,j]+=m

   # create hash table
   def makeLines(self):
      for m in range(64):
         line = []
         for i in range(6):
            if m&2**i==2**i:
               line.append(direcs[i])
         self.lines.append(line)

   # display
   def display(self, ax):
      for i in range(2, self.size+2):
         for j in range(2,self.size+2):
            pt = np.array([i*sqrt3half, (j+0.5) if i%2==1 else j])
            lst = [(pt,pt+direc) for direc in self.lines[self.active[i,j]]]
            # max = np.max(self.internal)
            # grey = (self.internal[i,j]/max,self.internal[i,j]/max,self.internal[i,j]/max,1)
            # colors=((grey for direc in self.lines[self.active[i,j]])),
            ax.add_collection(mc.LineCollection(lst, linewidths=2))

      # for making hex dot grid
      # x,y = np.mgrid[0:self.size+4,0:self.size+4]
      # x = np.where(y%2==1,x+0.5,x)
      # y = y*math.sqrt(3)/2
      # ax.scatter(y,x)

      ax.autoscale()
      ax.margins(0.1)
      # plt.draw()

model = Model()
fig, ax = plt.subplots()
plt.ion()
model.display(ax)
for i in range(10):
   plt.pause(0.2)
   model.update()
   model.display(ax)



