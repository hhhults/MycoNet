import model as m
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import time

labels = ["Dp","v","Di","Da","b","c1","c2","c3","c4","c5"]

def plot_avg_length_vs_time(params,iterations=10,n=20,timesteps=50,ax=None):
   data = np.zeros((iterations,timesteps))
   for it in range(iterations):
      model = m.Model(params,n=n, timeSteps=timesteps)
      data[it,0] = get_total_length(model)
      for t in range(1,timesteps):
         model.update()
         data[it,t] = get_total_length(model)
   if ax is None:
      plt.plot(np.arange(timesteps), np.mean(data,axis=0))
      plt.show()
   else:
      ax.plot(np.arange(timesteps), np.mean(data,axis=0))

def plot_avg_density_vs_time(params,iterations=10,n=15,timesteps=20,ax=None):
   data = np.zeros((iterations,timesteps))
   for it in range(iterations):
      model = m.Model(params,n=n, timeSteps=timesteps)
      data[it,0] = get_mean_length(model)
      for t in range(1,timesteps):
         model.update()
         data[it,t] = get_mean_length(model)
   if ax is None:
      plt.plot(np.arange(timesteps), np.mean(data,axis=0))
      plt.title("Average Density vs Time")
      plt.show()
   else:
      ax.plot(np.arange(timesteps), np.mean(data,axis=0))

def plot_avg_degree_vs_time(params,iterations=10,n=15,timesteps=20,ax=None):
   data = np.zeros((iterations,timesteps))
   for it in range(iterations):
      model = m.Model(params,n=n, timeSteps=timesteps)
      data[it,0] = get_mean_degree(model)
      for t in range(1,timesteps):
         model.update()
         data[it,t] = get_mean_degree(model)
   if ax is None:
      plt.plot(np.arange(timesteps), np.mean(data,axis=0))
      plt.show()
   else:
      ax.plot(np.arange(timesteps), np.mean(data,axis=0))

def plot_avg_nutrients_vs_time(params,iterations=10,n=15,timesteps=40,ax=None):
   data = np.zeros((iterations,timesteps))
   for it in range(iterations):
      model = m.Model(params,n=n, timeSteps=timesteps)
      data[it,0] = get_mean_nutrients(model)
      for t in range(1,timesteps):
         model.update()
         data[it,t] = get_mean_nutrients(model)
   if ax is None:
      plt.plot(np.arange(timesteps), np.mean(data,axis=0))
      plt.title("Average Nutrients vs Time")
      plt.show()
   else:
      ax.plot(np.arange(timesteps), np.mean(data,axis=0))

def do_sensitivity_analysis(filename, timesteps,iter, xvals, yvals,xindex,yindex):
   fig1, ax1 = plt.subplots(yvals.shape[0],xvals.shape[0],sharex=True,sharey=True)
   fig1.suptitle("Total Length vs Time When Varying "+labels[xindex]+" and "+labels[yindex])
   fig2, ax2 = plt.subplots(yvals.shape[0],xvals.shape[0],sharex=True,sharey=True)
   fig2.suptitle("Density vs Time When Varying "+labels[xindex]+" and "+labels[yindex])
   fig3, ax3 = plt.subplots(yvals.shape[0],xvals.shape[0],sharex=True,sharey=True)
   fig3.suptitle("Mean Degree vs Time When Varying "+labels[xindex]+" and "+labels[yindex])

   for i in range(yvals.shape[0]):
      params = defaultParams
      params[yindex] = yvals[i]
      for j in range(xvals.shape[0]):
         params[xindex] = xvals[j]
         data1 = np.zeros((iter,timesteps))
         data2 = np.zeros((iter,timesteps))
         data3 = np.zeros((iter,timesteps))

         for it in range(iter):
            model = m.Model(params, n=20,timeSteps=timesteps)
            data1[it,0] = get_total_length(model)
            data2[it,0] = get_mean_length(model)
            data3[it,0] = get_mean_degree(model)
            for t in range(1,timesteps):
               model.update()
               data1[it,t] = get_total_length(model)
               data2[it,t] = get_mean_length(model)
               data3[it,t] = get_mean_degree(model)
               ax1[i,j].plot(np.arange(timesteps), np.mean(data1,axis=0))
               ax2[i,j].plot(np.arange(timesteps), np.mean(data2,axis=0))
               ax3[i,j].plot(np.arange(timesteps), np.mean(data3,axis=0))
   fig1.savefig(filename+"-TotalLength")
   fig2.savefig(filename+"-Density")
   fig3.savefig(filename+"-MeanDegree")

def compute_all_analyses(timesteps,iter,xvals,yvals):
   indices = (0,1,2,3,4)
   for i in indices:
      for j in indices:
         if(i>j):
            do_sensitivity_analysis(labels[j]+"-vs-"+labels[i],timesteps,iter,xvals,yvals,j,i)


def get_total_length(model):
   return np.sum(model.count_ones(model.active))

def get_mean_length(model):
   count = model.count_ones(model.active)
   return np.sum(count[np.where(count>0)])/np.sum(np.where(count>0,1,0))
   # return ma.mean(ma.masked_less(model.count_ones(model.active),1))

def get_mean_nutrients(model):
   return np.mean(model.internal[np.where(model.internal>0)])

def get_mean_degree(model):
   count = model.count_ones(model.active)
   if (sum:=np.sum(np.where(count>2,1,0)))>0:
      return np.sum(count[np.where(count>2)])/sum
   else:
      return 0
   # return ma.mean(ma.masked_less(model.count_ones(model.active),3))
   
#        Dp v Di Da b c1 c2 c3 c4 c5
params = [1e4, 1e5, 2, 1, 1e6, 10, 1e-7, 0.001, 1e-8, 1e-9] # from paper`

defaultParams = [0.3, 1, 0.1,0.1,0.6,0.015,0.2,0.05,0.05,0.01]

xvals = 10.**np.array([-3,-1,1,3,5])
yvals = 10.**np.array([-3,-1,1,3,5])
# compute_all_analyses(50,10,xvals,yvals)#plot_avg_length_vs_time(params)
# make_display(defaultParams,n=50, timesteps=50)
plot_avg_density_vs_time(defaultParams,iterations=100)



         







