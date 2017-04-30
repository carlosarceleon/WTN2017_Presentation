from data_cleaning_routines       import \
        read_davis_tecplot_folder_and_rotate_to_serration_surface
from checklist                    import checklist
from raw_data_processing_routines import run_raw_data_collection
from boundary_layer_routines      import run_bl_analysis
from data_extraction_routines     import run_data_extraction
from data_analysis_routines       import do_the_time_resolved_analysis
#from publish                      import publish

def average_tecplot_dat_to_pandas_p(folder, 
                                    plot                = False,
                                    ignore_checklist    = False,
                                    rotate_to_airfoil   = True ,
                                    rotate_to_mean_flow = False,
                                    break_after_one     = True
                                   ):
    import os
    for tecplot_folder in [f for f in os.listdir( folder )]:
        if not checklist[tecplot_folder] or ignore_checklist:
            if 'STE' in tecplot_folder or 'z10' in tecplot_folder:
                rotate_to_airfoil_single   = True
                rotate_to_mean_flow_single = False
            else:
                rotate_to_airfoil_single   = False
                rotate_to_mean_flow_single = rotate_to_mean_flow

            print "  Processing {0}".format(tecplot_folder)
            read_davis_tecplot_folder_and_rotate_to_serration_surface(
                tecplot_folder      = os.path.join(folder,tecplot_folder),
                plot                = plot,
                rotate_to_airfoil   = rotate_to_airfoil_single,
                rotate_to_mean_flow = False,
            )
            if rotate_to_mean_flow_single:
                read_davis_tecplot_folder_and_rotate_to_serration_surface(
                    tecplot_folder      = os.path.join(folder,tecplot_folder),
                    plot                = plot,
                    rotate_to_airfoil   = False,
                    rotate_to_mean_flow = rotate_to_mean_flow_single,
                )

            if break_after_one:
                return

rotate_to_airfoil   = False
rotate_to_mean_flow = False

output_root = '/media/carlos/6E34D2CD34D29783/' +\
        '2015-02_SerrationPIV/TR_Data_Location_Calibrated_Article3'

root = '/media/carlos/6E34D2CD34D29783/' +\
        '2015-02_SerrationPIV/TR_Data_Processing_2016-05'

# 1) ###########################################################################
# Extract the averaged time-resolved results from the TECPLOT files from DaVis
# to a pickled pandas dataframe pickle #########################################
#
#average_tecplot_dat_to_pandas_p(
#    '/home/carlos/Documents/PhD/Articles/Article_3/Scripts/time_resolved'+\
#    '/tecplot_data_avg',
#    plot                = True,
#    ignore_checklist    = False,
#    rotate_to_airfoil   = rotate_to_airfoil,
#    rotate_to_mean_flow = rotate_to_mean_flow,
#    break_after_one     = False
#)
# ##############################################################################

# 2) ###########################################################################
# Do the boundary layer analysis on the averaged time-resolved pickles and put
# the boundary layer parameters in a pickle of its own #########################
#
#run_bl_analysis()
#
# ##############################################################################

# 3) ###########################################################################
# The pre-processing of the data, from raw TECPLOT to an aligned data frame ####
#
#run_raw_data_collection( 
#    only_for = [
#        'Sr20R21_a0_p0_U20_z05_tr', 
#        #'Sr20R21_a12_p0_U20_z05_tr', 
#        #'Sr20R21_a-12_p0_U20_z05_tr', 
#        'Sr20R21_a12_p6_U20_z05_tr', 
#        'Sr20R21_a-12_p6_U20_z05_tr', 
#        'STE_a-12_p0_U20_z00_tr', 
#        #'STE_a12_p0_U20_z00_tr', 
#    ], 
#    overwrite           = False ,
#    rotate_to_airfoil   = rotate_to_airfoil,
#    rotate_to_mean_flow = rotate_to_mean_flow,
#    output_root         = output_root,
#    plot_first          = True,
#    root                = root
#)
#
# ##############################################################################

# 4) ###########################################################################
# The data extraction of the interesting coordinate timeseries from the aligned
# data frame created in the previous step ######################################
#
#run_data_extraction( overwrite = False )
#
# ##############################################################################

# 5) ###########################################################################
# Do the data analysis #########################################################
#
do_the_time_resolved_analysis()
#
# ##############################################################################

#publish()
