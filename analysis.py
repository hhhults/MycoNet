import model as m
import numpy as np
import matplotlib.pyplot as plt



def make_display(params,filename=None, n=15,timesteps=20):
   model = m.Model(params,n=n, timeSteps=timesteps)
   fig, ax = plt.subplots()
   plt.ion()
   model.display(ax)
   for i in range(timesteps):
      plt.pause(0.1)
      model.update()
      model.display(ax)

def plot_avg_length_vs_time(params,iterations=10,n=15,timesteps=20):
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
      

