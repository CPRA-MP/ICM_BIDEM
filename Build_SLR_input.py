#Build the SLR input file for BIMODE
#input = TideData.csv, output = SLR_recordformodulation.txt
#should run once at the being of the 50 year simulation

import numpy as np

tide_input_directory = r'C:\Users\madel\Coastal Hydro Dropbox\Madeline Foster-Martinez\MRF_BarrierIslands'
tide_file = r'%s\TideData.csv' % tide_input_directory
tide_data = np.genfromtxt(tide_file,skip_header=1,delimiter=',',dtype='str')

ESLR_out_file_path = r'C:\Users\madel\Coastal Hydro Dropbox\Madeline Foster-Martinez\MRF_BarrierIslands'
ESLR_out_file_name = r'SLR_record4modulation.txt'
ESLR_out_file = r'%s\%s' % (ESLR_out_file_path, ESLR_out_file_name)

first_year = int(tide_data[0][0][0:4])
last_year = first_year+50 # alternative if it's not always 50 years: last_year = int(tide_data[-1][0][0:4])+1
years = range(first_year, last_year)

annual_mean_mm = [] 
for year in range(first_year,last_year):
    n = 0
    annual_total = 0
    for r in tide_data:
        if int(r[0][0:4]) == year:
            annual_total += float(r[3]) #4th column in TideData is Amerada Pass, LA
            n+=1 #8760 hours or 8784 hours if leap year      
    annual_mean_mm.append((annual_total/n)*1000) #take the mean of the hourly data and convert from m to mm
p = np.polyfit(years, np.asarray(annual_mean_mm),2) #fit a quadratic to the annual mean data
ESLR_rate_mmyr = (p[0]*2*years)+(p[1]) #take the first derivative and plug in years to get the rate
            
with open(ESLR_out_file, mode='w') as outf: 
    for i in range(0,len(years)):
        outf.write("%s %s\n" % (years[i],ESLR_rate_mmyr[i]))
                            
outf.close()     