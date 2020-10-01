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

with open(ESLR_out_file, mode='w') as outf:
    year_previous_mean_mm = 0
    for year in range(first_year,last_year):
        n = 0
        annual_total = 0
        for r in tide_data:
            if int(r[0][0:4]) == year:
                annual_total += float(r[3]) #4th column in TideData is Amerada Pass, LA
                n+=1 #8760 hours or 8784 hours if leap year 
        year_mean_mm = ((annual_total/n)*1000) #take the mean of the hourly data and convert from m to mm
        year_ESLR_mm = year_mean_mm-year_previous_mean_mm #subtract previous year to get the rate
        year_previous_mean_mm = year_mean_mm
        if year == first_year: #extrapolate the SLR rate from year 2 to 3 to year 1 to 2; otherwise the SLR in year 1 is off
            pass
        elif year == first_year +1:
            year2_ESLR_mm = year_ESLR_mm
        elif year == first_year +2:
            delta_SLR = year_ESLR_mm - year2_ESLR_mm
            year1_ESLR_mm = year2_ESLR_mm - delta_SLR
            outf.write('%s %s\n' % (first_year, year1_ESLR_mm))
            outf.write('%s %s\n' % (first_year+1, year2_ESLR_mm))                
            outf.write('%s %s\n' % (year, year_ESLR_mm))
        else:
            outf.write('%s %s\n' % (year, year_ESLR_mm))