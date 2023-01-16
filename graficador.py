#Graficador parte 1

import numpy as np
import matplotlib.pyplot as plt
PC= 3.085677588e+18

archivo = "10MSUN-"      #input("Escribe el nombre del archivo out: ")


#poner nombre de input

def grafica(i,x,rho,u,P):
	n = 1.0 
	K = 1.38e-16
	# T = P/(n*K)	
	fig, axs = plt.subplots(3, 1, constrained_layout=True, figsize=(5, 10))
	fig.suptitle(str(i*1000)+"yr", fontsize=16)

	axs[0].plot(x, rho)
	axs[0].set_title('Densidad')
	# axs[0].set_ylim([0,1.8e-23])
	axs[0].set_xlabel("[pc]")
	axs[0].set_ylabel("g/cm^3")

	
	axs[1].plot(x, u, 'tab:orange')
	axs[1].set_title('Velocidad')
	axs[1].set_ylim([0,8.5e8])
	axs[1].set_xlabel("[pc]")
	axs[1].set_ylabel("m/s")
	
	axs[2].plot(x, P, 'tab:green')
	axs[2].set_title('Presion')
	axs[2].set_ylim([-5e-7,6.5e-6])
	axs[2].set_xlabel("[pc]")
	axs[2].set_ylabel("dyn/cm^2")
	
	# axs[1, 1].plot(x, T, 'tab:red')
	# axs[1, 1].set_title('Temperatura')
	# axs[1, 1].set_ylim([-14,4.5e10])
	if i<10:
		fig.savefig(archivo+"0"+str(i)+".png")
	else:
		fig.savefig(archivo+ str(i)+".png")
	# fig.show()
	
	

#for ax in axs.flat:
 #   ax.set(xlabel='x-label', ylabel='y-label')

# Hide x labels and tick labels for top plots and y ticks for right plots.
#for ax in axs.flat:
 #   ax.label_outer()

for i in range(50):
	if i < 10: 
		data = np.loadtxt(archivo+"0"+str(i)+".dat")
		grafica(i, x=data[:,0], rho=data[:,1], u=data[:,2], P=data[:,3])
	else:
		data = np.loadtxt(archivo+ str(i) +".dat")
		grafica(i, x=data[:,0], rho=data[:,1], u=data[:,2], P=data[:,3])


	
# data = np.loadtxt("HLL/10MSUN-00.dat")
# grafica(i, x=data[:,0], rho=data[:,1], u=data[:,2], P=data[:,3])

