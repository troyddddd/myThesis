

def makegraph(x,y,name):
	import matplotlib.pyplot as plt
	fig = plt.figure(name)
	plt.plot(x,y)
	fig.savefig(name)
