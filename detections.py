import numpy as np
import pycbc.noise
import pycbc.psd
import pycbc.filter
import pycbc.waveform
from pycbc import frame
import matplotlib.pyplot as plt
import scipy
import style

plt.style.use(style.style1)

def search_for_signal(filename):

    flow = 20.0
    delta_f = 1.0 / 800
    flen = int(2048 / delta_f) + 1
    
    best_snr = 0
    best_time = 0
    best_mass = 0
    
    if filename.find('H1') != -1:
    # Get Ligo PSD 
        psd = pycbc.psd.read.from_txt('psd_data/aLIGO_O4_high_psd.txt', 
                                  length=flen, delta_f=delta_f,
                                  low_freq_cutoff = flow, 
                                  is_asd_file=False)
        strain = frame.read_frame(filename, 'H1:LOSC-STRAIN')
        
    elif filename.find('L1') != -1:
        # Get Ligo PSD 
        psd = pycbc.psd.read.from_txt('psd_data/aLIGO_O4_high_psd.txt', 
                                      length=flen, delta_f=delta_f,
                                      low_freq_cutoff = flow, 
                                      is_asd_file=False)
        strain = frame.read_frame(filename, 'L1:LOSC-STRAIN')
                                              
    elif filename.find('V1') != -1:                              
        # Get Virgo PSD
        psd = pycbc.psd.read.from_txt('psd_data/AdV_psd.txt',
                                            length=flen,
                                            delta_f=delta_f,
                                            low_freq_cutoff = flow,
                                            is_asd_file=False)
        strain = frame.read_frame(filename, 'V1:LOSC-STRAIN')

    # Filter low frequencies bellow 20 Hz
    strain = pycbc.filter.highpass(strain, 20.0) 
    stilde = strain.to_frequencyseries()
  
    for mass in [1.4, 10, 30, 50] :
    
        # Use a waveform as a matched filter
        hp, hc = pycbc.waveform.get_fd_waveform(approximant='IMRPhenomD',
                                                mass1=mass,
                                                mass2=mass,
                                                f_lower=flow,
                                                delta_f=stilde.delta_f)
        hp.resize(len(stilde))
    
        # match-filter signal_noise
        snr = pycbc.filter.matched_filter(hp,
                                          stilde,
                                          psd=psd,
                                          low_frequency_cutoff=flow)
        peak = abs(snr).numpy().argmax()
        snrp = snr[peak]
        time = snr.sample_times[peak]
        
        if abs(snrp) > best_snr :
            best_snr = abs(snrp)
            best_time = time
            best_mass = mass
            
            plt.close()
            plt.figure()
            plt.title("SNR of data from " + filename[7:19] +"\n(template "+ str(best_mass) + "M$_{\odot}$)\n")
            plt.plot(snr.sample_times, abs(snr))
            plt.ylabel('signal-to-noise ratio')
            plt.xlabel('time (s)')
            plt.savefig('./plots/'+filename[filename.find('1_')-1:-4]+'.png')
            plt.close()
            
            
            plt.figure()
            plt.title('Histogram of SNR distribution for '+ filename[7:19] +"\n(template "+ str(best_mass) + "M$_{\odot}$)\n")
            plt.hist(abs(snr), 100)
            plt.ylabel('number of counts')
            plt.xlabel('Signal-to-noise ratio')
            plt.savefig('./plots/distribution_'+filename[filename.find('1_')-1:-4]+'.png')
            plt.close()
   
   
    return best_snr, best_time, best_mass

def compare_times(file_L,file_H,file_V):
    
    pathfile = 'group4/'
    
    if file_L[2:] != file_H[2:] or file_L[2:] != file_V[2:]:
        print("\nWARNING : You are comparing different events.\n", file_L, "\n", file_H,"\n", file_V)
        exit()
        
    snr_ligo_L, time_L, mass_L = search_for_signal(pathfile + file_L)
    snr_ligo_H, time_H, mass_H = search_for_signal(pathfile + file_H)
    snr_virgo, time_V, mass_V = search_for_signal(pathfile + file_V)

    c = 3*10**5 #km/s
    d_HL = 3001 #km
    d_HV = 8181 #km

    tc_HL = d_HL/c
    tc_HV = d_HV/c

    t_HL = abs(time_H-time_L)
    t_HV = abs(time_H-time_V)

    print("\nSIGNAL : ",file_H[3:13])
    if t_HL < tc_HL:
        if t_HV < tc_HV:
            print("Coincidence confirmed between Hanford, Livingston and Virgo")
            print("The SNR of Hanford is : ", snr_ligo_H, "\nThe SNR of Livingston is : ", snr_ligo_L, "\nLe SNR de Virgo is :", snr_virgo)
            print("The mass found by Hanford is : ", mass_H, "\nThe mass found by Livingston is : ", mass_L, "\nThe mass found by Virgo is : ", mass_V)
            print("The arriving time mesured by Hanford is : ", time_H, "\nThe arriving time mesured by Livingston is : ", time_L, "\nThe arriving time mesured by Virgo is : ", time_V)
        
        else :
            print("Coincidence only confirmed between Hanford and Livinstone")
            print("The SNR of Hanford is : ", snr_ligo_H, "\nThe SNR of Livingston is : ", snr_ligo_L)
            print("The mass found by Hanford is : ", mass_H, "\nThe mass found by Livingston is : ", mass_L)
            print("The arriving time mesured by Hanford is : ", time_H, "\nThe arriving time mesured by Livingston is : ", time_L)
    else:
        if t_HV < tc_HV:
            print("Coincidence only confirmed between Hanford and Virgo")
            print("The SNR of Hanford is : ", snr_ligo_H, "\nThe SNR of Livingston is : ", snr_virgo)
            print("The mass found by Hanford is : ", mass_H, "\nThe mass found by Virgo is : ", mass_V)
            print("The arriving time mesured by Hanford is : ", time_H, "\nThe arriving time mesured by Virgo is : ", time_V)
        
        else:
            print("NO coincidence confirmed there is no signal")
            print("The SNR of Hanford is : ", snr_ligo_H, "\nThe SNR of Livingston is : ", snr_ligo_L, "\nLe SNR de Virgo is :", snr_virgo)
            print("The mass found by Hanford is : ", mass_H, "\nThe mass found by Livingston is : ", mass_L, "\nThe mass found by Virgo is : ", mass_V)
            print("The arriving time mesured by Hanford is : ", time_H, "\nThe arriving time mesured by Livingston is : ", time_L, "\nThe arriving time mesured by Virgo is : ", time_V)
        
        
        
        
        
