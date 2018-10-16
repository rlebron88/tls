import numpy as np

def color(meth):
	'''
	usage: color(meth)
	
	input:
	meth (float): percentage of reads that show methylation at this position in the genome [0, 100]
	
	output:
	color (str): red [0, 255], green [0, 255], blue [0, 255] bedMethyl color values
	'''
	#bins = np.array([5, 15, 25, 35, 45, 55, 65, 75, 85, 95, 100])
	bins = np.array([-0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9, 1])
	meth_bin = np.digitize(np.array([meth]), bins, right = True)[0]
	G_bin = np.clip(meth_bin, 0, 5)
	R_bin = np.clip(10 - meth_bin, 0, 5)
	G = 50 * G_bin + 5 if G_bin != 0 else 0
	R = 50 * R_bin + 5 if R_bin != 0 else 0
	return ",".join(map(str, [R, G, 0]))
