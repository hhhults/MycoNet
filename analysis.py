import Model as m
import numpy as np
import matplotlib.pyplot as plt
import time



def make_display(params,filename=None, n=15,timesteps=20):
   model = m.Model(params,n=n, timeSteps=timesteps)
   fig, ax = plt.subplots()
   plt.ion()
   model.display(ax)
   for i in range(timesteps):
      plt.pause(0.1)
      model.update()
      model.display(ax)

def plot_avg_length_vs_time(params,iterations=10,n=15,timesteps=40):
   data = np.zeros((iterations,timesteps))
   for it in range(iterations):
      model = m.Model(params,n=n, timeSteps=timesteps)
      data[it,0] = model.get_total_length()
      for t in range(1,timesteps):
         model.update()
         data[it,t] = model.get_total_length()
   plt.plot(np.arange(timesteps), np.mean(data,axis=0))
   plt.show()

   
#        Dp v Di Da b c1 c2 c3 c4 c5
params = [1e4, 1e5, 2, 1, 1e6, 10, 1e-7, 0.001, 1e-8, 1e-9]

params = [0.3, 1, 0.1,0.1,0.6,0.015,0.2,0.05,0.05,0.01]

plot_avg_length_vs_time(params)
# make_display(params,n=30, timesteps=30)

def timeRuns(n):
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

   plt.plot([5,10,15,20,25],timeTable,)
   plt.xlabel("Number of Runs")
   plt.ylabel("Time(seconds)")
   plt.title("Runtime")
   plt.show()
   plt.pause(1.5)

def testDpDi():
   Dp = .3  # coefficient of diffusive movement
   v = 1  # coefficient of active movement
   Di = .1  # coefficient of diffusive transport
   Da = .1  # coefficient of active transport
   b = .6  # branching coefficient
   c1 = .015  # coefficient for internal gain of nutrient through uptake
   c2 = .2  # growth cost coefficient
   c3 = .05  # coefficient for environmental loss of nutrient through uptake
   c4 = .05  # active translocation cost coefficient
   c5 = .01

   parameters = [Dp, v,Di, Da, b,c1,c2,c3,c4,c5]
   density = np.zeros(5)

   for i in range(5):

      testModel = m.Model(parameters, 15)

      fig, ax = plt.subplots()
      plt.ion()
      testModel.display(ax)
      for i in range(20):
         plt.pause(0.2)
         testModel.update()
         testModel.display(ax)


      plt.close()


test = timeRuns(5)



