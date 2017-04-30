import misc_functions as func
#from shutil import copy
import analyze_crossovers as xover
#import acoustic_functions as afunc
import bl_functions as bl
import os
#from os.path import split
#from publish import publish

raw_data_root = '/media/carlos/6E34D2CD34D29783/2015-07_BL/STE_BL_Data/'
try:
    raw_data_folders = [f for f in os.listdir(raw_data_root) \
                    if os.path.isdir(os.path.join(raw_data_root,f))]
except OSError:
    raw_data_folders = []


#bl.plot_article_deltas( overwrite = False )

#func.plot_interesting_cases( phi = 6, U = 30 )
##func.plot_interesting_cases( phi = 0, U = 30 )
#func.plot_interesting_cases( phi = 6, U = 35 )
#func.plot_interesting_cases( phi = 0, U = 35 )
#func.plot_interesting_cases( phi = 6, U = 40 )
##func.plot_interesting_cases( phi = 0, U = 40 )
#func.plot_article_relative_cases(alpha=0  , phi=6, article=True)
##func.plot_article_relative_cases(alpha=0  , phi=0, article=True)
##func.plot_article_relative_cases(alpha=6  , phi=0, article=True)
#func.plot_article_relative_cases(alpha=6  , phi=6, article=True)
##func.plot_article_relative_cases(alpha=12 , phi=0, article=True)
#func.plot_article_relative_cases(alpha=12 , phi=6, article=True)
##
#xover.plot_bl_crossover_relation( )

xover.predict_f_at_twenty()

#publish()
