# MoralWordEEG
Experiment, Materials, and Analysis code for Moral Word EEG project
"The time-course of moral perception: An ERP investigation of the moral pop-out effect"
Contributors: Ana P. Gantman Jay J Van Bavel Peter Mende-Siedlecki Kyle Elliott Mathewson

https://osf.io/5jmze/

Raw data is in EEG and BEH folder on OSF, 
Preproccessed segments are in the EEG folder in a subfolder 4conds_rej3
This code uses eeglab version 13_3_2b

download and change paths in the following script:
1) Analysis_Wrapper_MWP.m - This preprocesses the raw data and saves the preprocessed data segments
  -Preprocessing_MWP.m - subfunction

2) PlotERP_MWP.m - This loads the preprocessed segments plots the ERP's and topographies to make figures 3 and 4
  -Load_data_MWP.m - subfunction to load data

3) SingleTrial_MWP.m - This loads the behavioural data and creates a list of each trial to use for analysis, outputs a testout.csv file of the data that can be pasted into the columns of GEE_Data_MWP.xlsx.

4) SingleTrial_ERP_MWP.m - This loads the behavioural data and EEG data together and creates a list of each trial to use for analysis, outputs a testout.csv file of the data that can be pasted into the columns of GEE_Data_ERP_MWP.xlsx.

subfunctions - folder of functions used in the code
EOG-electrode-locs.ced - locations of electrodes for eeglab topographies
4conds_rej3_Settings_MWP.mat - saved preprocessing settings made by Analysis_Wrapper_MWP.m
