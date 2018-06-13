# MoralWordEEG
Experiment, Materials, and Analysis code for Moral Word EEG project
"The time-course of moral perception: An ERP investigation of the moral pop-out effect"
Contributors: Ana P. Gantman Jay J Van Bavel Peter Mende-Siedlecki Kyle Elliott Mathewson

https://osf.io/5jmze/

Raw data is in EEG and BEH folder on OSF, 
Preproccessed segments are in the EEG folder in a subfolder 4conds_rej3
This code uses eeglab version 13_3_2b

download and change paths in the following script:
Analysis_Wrapper_MWP.m - This preprocesses the raw data and saves the preprocessed data segments
  -Preprocessing_MWP.m - subfunction

PlotERP_MWP.m - This loads the preprocessed segments plots the ERP's and topographies to make figures 3 and 4
  -Load_data_MWP.m - subfunction to load data

SingleTrial_MWP.m - This loads the behavioural data and EEG data together and creates a list of each trial and the ERP magnitude on that trial to use for the GMM analysis, outputs a testout.csv file of the erp magnitudes that can be pasted into GEE_Data_MWP.xlsx after being resorted as indicated in the comments

subfunctions - folder of functions used in the code
EOG-electrode-locs.ced - locations of electrodes for eeglab topographies
4conds_rej3_Settings_MWP.mat - saved preprocessing settings made by Analysis_Wrapper_MWP.m
