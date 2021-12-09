import model as m
import numpy as np
import time
import matplotlib.pyplot as plt

def make_display(params=m.Model.defaultParams,filename=None, n=15,timesteps=20):
   model = m.Model(params,n=n, timeSteps=timesteps)
   fig, ax = plt.subplots()
   plt.ion()
   model.display(ax)
   for i in range(timesteps):
      plt.pause(0.1)
      model.update()
      model.display(ax)

def timeRuns(params,n):
   increment = n
   
   timeTable = np.zeros(n)

   for q in range(n):
      print(n)
      start = time.time()

      model = m.Model(params, n)
      fig, ax = plt.subplots()
      plt.ion()
      model.display(ax)
      for i in range(20):
         plt.pause(0.2)
         model.update()
         model.display(ax)

      end = time.time()

      timeTable[q] = (end - start)
      plt.close()
      n += increment

   plt.plot(timeTable)
   plt.xlabel("Number of Runs")
   plt.ylabel("Time(seconds)")
   plt.title("Runtime")
   plt.show()
   plt.pause(1.5)