# fonction search_for_signal(string)
# return : tempsgps et SNR
# filename.find('truc') = -1 s'il ne le trouve pas 
from detections import search_for_signal, compare_times

data_L = ['L1_1359826218_800.gwf','L1_1359827018_800.gwf','L1_1359827818_800.gwf','L1_1359828618_800.gwf','L1_1359829418_800.gwf','L1_1359830218_800.gwf']
data_H = ['H1_1359826218_800.gwf','H1_1359827018_800.gwf','H1_1359827818_800.gwf','H1_1359828618_800.gwf','H1_1359829418_800.gwf','H1_1359830218_800.gwf']
data_V = ['V1_1359826218_800.gwf','V1_1359827018_800.gwf','V1_1359827818_800.gwf','V1_1359828618_800.gwf','V1_1359829418_800.gwf','V1_1359830218_800.gwf']    

for i in range(len(data_L)):
    compare_times(data_L[i],data_H[i],data_V[i])

