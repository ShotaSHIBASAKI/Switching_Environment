import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
def data_generate2(nu):
    #for two consumer model
    time=200
    trial=1000

    resource=np.zeros([time+1, trial])
    consumer1=np.zeros([time+1, trial])
    consumer2=np.zeros([time+1, trial])
    toxin=np.zeros([time+1, trial])
    fold_name=str('./nu10^%d' %(np.log10(nu)))
    os.chdir(fold_name)
    
    for i in range(trial):
        fname=str('Convergence_RTC_model1_switching10^%d_trial%d.csv'%(np.log10(nu), i))
        data=np.loadtxt(fname, delimiter=',', skiprows=1)
        # data at t=0
        resource[0,i]=data[0,1]
        consumer1[0,i]=data[0,2]
        consumer1[0,i]=data[0,3]
        toxin[0,i]=data[0,4]
        # data at t=1,2,...199
        counter=1
        j=1
        while counter<200 and j < np.size(data,0)-1:
            if data[j,0]<counter and data[j+1,0]>=counter:
                del_neg=counter-data[j,0]
                del_pos=data[j+1,1]-counter
                resource[counter, i]=(del_pos*data[j, 1]+del_neg*data[j+1, 1])/(del_pos+del_neg)
                consumer1[counter, i]=(del_pos*data[j, 2]+del_neg*data[j+1, 2])/(del_pos+del_neg)
                consumer2[counter, i]=(del_pos*data[j, 3]+del_neg*data[j+1, 3])/(del_pos+del_neg)
                toxin[counter, i]=(del_pos*data[j, 4]+del_neg*data[j+1, 4])/(del_pos+del_neg)
                counter+=1
            j+=1
        # data at t=200
        resource[-1,i]=data[-1,1]
        consumer1[-1,i]=data[-1,2]
        consumer1[-1,i]=data[-1,3]
        toxin[-1,i]=data[-1,4]
    np.savetxt('TwoConsumer_resource_hist.csv', resource, delimiter=',', fmt='%.3f')
    np.savetxt('TwoConsumer_consumer1_hist.csv', consumer1, delimiter=',', fmt='%.3f')
    np.savetxt('TwoConsumer_consumer2_hist.csv', consumer2, delimiter=',', fmt='%.3f')
    np.savetxt('TwoConsumer_toxin_hist.csv', toxin, delimiter=',', fmt='%.3f')
    return [resource, consumer1, consumer2, toxin]



def GifAnimation(data, number_of_frames, col_list, xlabel, gname):
    """
    data=[data1, data2]: data set to generate animation of histgram
    num_of_frames: time step (from 0 to time)
    col_list: color of each data inhistgram
    xlabel; 'resource','species', toxin''
    gname: str of graph name
    """
    
    def update_hist(num, data):
        plt.cla()
        data1=data[0]
        data2=data[1]
        plt.hist(data1[num], bins=20, color=col_list[0], range=(0, np.max(data)),
        density=True, alpha=0.5, label='species 1')
        plt.hist(data2[num], bins=20, color=col_list[1], range=(0, np.max(data)),
        density=True, alpha=0.5, label='species 2')
        plt.xlim(left=0, right=180)
        plt.xlabel(xlabel, fontsize=16)
        plt.ylabel('count', fontsize=16)
        plt.ylim(bottom=0, top=0.15)
        plt.title(str('time %d' %(num)))
        plt.legend(loc='upper right')
    fig = plt.figure()
    ani = animation.FuncAnimation(fig, update_hist, number_of_frames, fargs=(data,) )
    ani.save(gname+'_hist.gif', writer="pillow", dpi=300)




#data_generate2(0.1)
sp1=np.loadtxt('TwoConsumer_consumer1_hist.csv', delimiter=',', skiprows=0)
sp2=np.loadtxt('TwoConsumer_consumer2_hist.csv', delimiter=',', skiprows=0)
data=[sp1,sp2]
GifAnimation(data, 201, ['#8dd3c7','#b3de69'], 'abundance', 'TwoSpecies')