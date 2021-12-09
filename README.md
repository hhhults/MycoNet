# MycoNet

README / USER's MANUAL

MycoNet is an implementation of a model for a mycelial network described in "The Development of Fungal Networks in Complex
Environments" by Graeme P. Boswell, et al. (2007) (https://link.springer.com/content/pdf/10.1007/s11538-005-9056-6.pdf)

Included in this repository is the model, a set of analysis code, and some tests.

If you just want to see the network grow, use the make_display() function in testing.py to make a live visualization of the network. It gets pretty slow around n = 50.

If you want to do analysis, look to analysis.py where the function do_sensitivity_analysis allows you to compare different metrics of the model as you vary multiple parameters.

Besides the best testing which is watching it and seeing what it does, there is also the function timeRuns() in testing.py which allows the user to time a specific number of simulations and benchmark their speed.
This would be useful when optimizing this model, which is the obvious next step after creation so that it may be used more widely.