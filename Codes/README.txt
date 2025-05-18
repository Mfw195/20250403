
This software package includes two folders, "Model_Forward" and "Real_Seismic_data_forward"：
The "Model_Forward" folder includes the rock physics modeling section and the forward synthesis of seismic data. 
The "Real_Seismic_data_forward" folder is mainly used for reading and writing real seismic data and comparing it with synthetic data.

These two demos folders contain script files for rock physical models and a series of script files for synthesizing seismic data using the reflectance method
, reproducing all figures presented in the paper. The description of each testing script can be found below:

In the folder named "Model_Forward" :
"Fig2_Aspect_ratio_distribution_app":  It is used to obtain the microscopic pore structure of the core.
"Fig3_to_5_Two_scale_model_app":  It is used to obtain the dispersion attenuation of the two-scale petrophysical model.
"Fig_12_Forward_Sw_AVA_app" : AVA seismic data synthesized at different saturation levels.
"Fig_13_Forward_perm_AVA_app":  The synthesized AVA seismic data at different permeability levels.
"Fig_14_Forward_Sw_app":  Analyze the influence of saturation on seismic data when a certain Angle is fixed.
"Fig_15_Forward_perm_app"  : Analyze the influence of permeability on seismic data when a certain Angle is fixed
This folder contains two files titled 'model_data' and 'Experimental_data'. Other .m files are the required sub-function programs . Keep these two files and script above in a single folder.

In the folder named "Real_Seismic_data_forward" :
"Fig_16_and_17_seismic_data_app":  Used for reading pre-stack and post-stack seismic data.
"Fig_18_a_forward_log_app":  Wellbore bypass seismic records are synthesized from well logging data.
"Fig_19_and_20_a_forward_4_app":  The dispersion velocity obtained from the petrophysical model is used to replace the logging velocity to synthesize the seismic data of the wellbore bypass.
This folder contains two files titled 'log_data_new' and 'seismic_data'. Other .m files are the required sub-function programs . Keep these two files and script above in a single folder.
