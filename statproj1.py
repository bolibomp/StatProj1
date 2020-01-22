import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance
import time
   
Data_set_1 = np.loadtxt(r"C:\Users\Robin\Downloads\skol\P1data01.txt",usecols=(0, 1))
Data_set_2 = np.loadtxt(r"C:\Users\Robin\Downloads\skol\P1data02.txt",usecols=(0, 1))
Data_set_3 = np.loadtxt(r"C:\Users\Robin\Downloads\skol\P1data03.txt",usecols=(0, 1))
Data_set_4 = np.loadtxt(r"C:\Users\Robin\Downloads\skol\P1data04.txt",usecols=(0, 1))
Data_set_5 = np.loadtxt(r"C:\Users\Robin\Downloads\skol\P1data05.txt",usecols=(0, 1))
Data_set_6 = np.loadtxt(r"C:\Users\Robin\Downloads\skol\P1data06.txt",usecols=(0, 1))
Data_list=[Data_set_1 ,Data_set_2 ,Data_set_3 ,Data_set_4, Data_set_5, Data_set_6] #List of the data sets so to make it able to loop through

def Nicetime(time): #Small function to get the time it took to run the code in minutes and seconds
    temp_1=int(time)
    temp_2=time-temp_1
    print('The code took', temp_1, 'minutes and', round(temp_2/0.0166666667,2),'seconds to complete')

def corr_func(data, index, minimize_stat_fluct=10):
    """
    Returns plots of the Natrual, Peebles, Landy & Szalay and Hamilton estimators for the given data.

    Parameters
    ----------
    data : 2D array
        Data for which the estimatores are to be plotted from.

    index: int
        For keeping track of the number of times the function have been called
        
    minimize_stat_fluct: int
        Number of times the distances in the random data set should be averaged.
        Default value set to 10
        
    -------
    """
    print('Data set ', index+1) #Print which data set is currently being processed
    Bins = np.linspace(0,1000,101)  #The bins for the for the histogram. Bins span from 0 to 1000 with witdh 10
    N_r=len(data) #Parameter for the estimators
    N=len(data) #Parameter for the estimators
    DR = 0 
    for i in range(minimize_stat_fluct): #Finding the distance 'min_stat_fluct'-times to minimize statistical fluctuations in DR
        print('DR: ',i+1, ' out of ',minimize_stat_fluct)
        a=(np.random.rand(data.shape[0], data.shape[1])*1000).round(2) #Making random data in the same size as the data
        temp_DR = np.asarray(distance.cdist(data,a, 'euclidean')) #Finding the distance betyween the data and the random data
        DR= DR + np.histogram(temp_DR, bins=Bins)[0] #making 'min_stat_fluct'-times histograms of the distance data
    DR=DR/minimize_stat_fluct  

    print('DR histogram done')
    RR = 0 
    for i in range(minimize_stat_fluct):
        print('RR: ',i+1, ' out of ',minimize_stat_fluct)
        a=(np.random.rand(data.shape[0], data.shape[1])*1000).round(2)
        temp_RR = np.asarray(distance.cdist(a,a, 'euclidean'))
        RR = RR + np.histogram(temp_RR, bins=Bins)[0]
    RR=RR/minimize_stat_fluct

    RR[0] = RR[0]-data.shape[0] #removing the zeros from the distance between two identical random data sets
    RR=RR/2 #removing double counts
    #RR_dist=0
    print('RR histogram done')
    DD_dist = distance.cdist(data, data, 'euclidean')
    DD = np.histogram(DD_dist, bins=Bins)[0]
    DD[0] = DD[0]-data.shape[0]
    DD=DD/2
    DD_dist=0
    print('DD histogram done')
    corr_1 = DD / RR - 1 #Natrual estimator
    corr_2 = 2*N_r/(N-1)*(DD)/(DR)-1 #Peebles estimator
    corr_3 = (N_r*(N_r-1))/(N*(N-1))*(DD)/(RR)-(N_r-1)/(N)*(DR)/(RR)+1 #Landy & Szalay estimator
    corr_4 = (4*N*N_r)/((N-1)*(N_r-1))*(DD)/(DR)*(RR)/(DR)-1 #Hamilton estimator
    Est_1.append(corr_1) #appedning to use for plots
    Est_2.append(corr_2)
    Est_3.append(corr_3)
    Est_4.append(corr_4)
    DD=0#reseting values
    RR=0
    DR=0
    corr_1=0
    corr_2=0
    corr_3=0
    corr_4=0
start_time = time.clock()
Est_1=[]
Est_2=[]
Est_3=[]
Est_4=[]
for idx, Data in enumerate(Data_list): #Running the function for each of the data sets  
    corr_func(Data,idx)
Nicetime((time.clock() - start_time)*0.0166666667)

for i in range(len(Est_1)):
    if i==1:
        plt.plot(np.linspace(5,995,100),Est_1[i], marker='o', markersize=1, label='Data set '+str(i+1))
    else:
        plt.plot(np.linspace(5,995,100),Est_1[i], label='Data set '+str(i+1))
plt.title('All the data sets with the natrual estimator' )
plt.ylabel(r'$w(\theta)$')
plt.xlabel(r'$\theta$')
plt.ylim(-0.06,0.03)
plt.legend(loc=8)
plt.show()

for i in range(len(Est_2)):
    if i==1:
        plt.plot(np.linspace(5,995,100),Est_2[i], marker='o', markersize=1, label='Data set '+str(i+1))
    else:
        plt.plot(np.linspace(5,995,100),Est_2[i], label='Data set '+str(i+1))
plt.title('All the data sets with the Davis and Peebles estimator' )
plt.ylabel(r'$w(\theta)$')
plt.xlabel(r'$\theta$')
plt.ylim(-0.06,0.03)
plt.legend(loc=8)
plt.show()


for i in range(len(Est_3)):
    if i==1:
        plt.plot(np.linspace(5,995,100),Est_3[i], marker='o', markersize=1, label='Data set '+str(i+1))
    else:
        plt.plot(np.linspace(5,995,100),Est_3[i], label='Data set '+str(i+1))
plt.title('All the data sets with the Hamilton estimator' )
plt.ylabel(r'$w(\theta)$')
plt.xlabel(r'$\theta$')
plt.ylim(-0.06,0.03)
plt.legend(loc=8)
plt.show()

for i in range(len(Est_4)):
    if i==1:
        plt.plot(np.linspace(5,995,100),Est_4[i], marker='o', markersize=1, label='Data set '+str(i+1))
    else:
        plt.plot(np.linspace(5,995,100),Est_4[i], label='Data set '+str(i+1))
plt.title('All the data sets with the Landy & Szalay estimator' )
plt.ylabel(r'$w(\theta)$')
plt.xlabel(r'$\theta$')
plt.ylim(-0.06,0.03)
plt.legend(loc=8)
plt.show()