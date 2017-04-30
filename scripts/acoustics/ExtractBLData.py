import os
import pandas as pd
import TELocations as teloc
import numpy as np
import re

file_dict = {
    'Avg_Vx'                            : "B00001.dat",
    'Avg_Vy'                            : "B00002.dat",
    'Length_of_Avg_V'                   : "B00003.dat",
    'Standard_deviation_of_Vx'          : "B00004.dat",
    'Standard_deviation_of_Vy'          : "B00005.dat",
    'Length_of_Standard_deviation_of_V' : "B00006.dat",
    'Turbulent_kinetec_energy'          : "B00007.dat",
    'Reynold_stress_XY'                 : "B00008.dat",
    'Reynold_stress_XX'                 : "B00009.dat",
    'Reynold_stress_YY'                 : "B00010.dat",
}

davis_dict = {
    'Avg_Vx'                            : 'Avg_Vx'                   ,
    'Avg_Vy'                            : 'Avg_Vy'                   ,
    'Length_of_Avg_V'                   : 'Length_of_Avg_V'          ,
    'Standard_deviation_of_Vx'          : 'RMS_Vx'                   ,
    'Standard_deviation_of_Vy'          : 'RMS_Vy'                   ,
    'Length_of_Standard_deviation_of_V' : 'Length_of_RMS_V'          ,
    'Turbulent_kinetec_energy'          : 'Turbulent_kinetec_energy' ,
    'Reynold_stress_XY'                 : 'Reynold_stress_XY'        ,
    'Reynold_stress_XX'                 : 'Reynold_stress_XX'        ,
    'Reynold_stress_YY'                 : 'Reynold_stress_YY'        ,
}

need_to_stitch_cases = [
    'STE_A6_U20_closed_SS',
    'STE_SS_a12_U20',
    'STE_SS_a12_U30',
    'STE_SS_a12_U35',
    'STE_SS_a12_U40',
]

BL_pct = 95

root = './Data'
raw_data_root = '/media/carlos/6E34D2CD34D29783/2015-07_BL/STE_BL_Data/'

data_folders = [f for f in os.listdir(root) \
                if os.path.isfile(os.path.join(root,f))
               and f.endswith(".p")]
try:
    raw_data_folders = [f for f in os.listdir(raw_data_root) \
                    if os.path.isdir(os.path.join(raw_data_root,f))]
except OSError:
    raw_data_folders = []

def read_data(root,case,variable):
    """ Reads the data. Prioritize reading an existing pickle
    of the data"""
    import os
    import pandas as pd

    if os.path.isfile(
        os.path.join(root,case)
    ):
        return pd.read_pickle(
            os.path.join(root,case)
        )
    else:
        tecplot_file = os.path.join(root,case,file_dict[variable])
        return read_tecplot(tecplot_file)

def stitch_cases(frame1_df,frame2_df,case,plot=False):
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.__version__

    te_location = teloc.TELocations[
        recognize_case(case)[0]
    ][1]

    y1 = find_nearest(float(te_location[1]),frame1_df.y.unique())
    y2 = find_nearest(float(te_location[1]),frame2_df.y.unique())

    frame1_df_TE = frame1_df[(frame1_df.y == y1)]
    frame2_df_TE = frame2_df[(frame2_df.y == y2)]

    near_frame_max = frame1_df_TE.Length_of_Avg_V.max()

    frame2_shift = frame2_df_TE[
        frame2_df_TE.Length_of_Avg_V == find_nearest(
            near_frame_max,
            frame2_df_TE.Length_of_Avg_V.values)
    ].x.values[-1]

    frame2_df_TE = frame2_df_TE[frame2_df_TE.x<=frame2_shift]
    frame2_df_TE.x = frame2_df_TE.x - \
            frame2_shift + frame1_df_TE.x.min()
    
    if plot:
        fig = plt.figure()

        plt.plot(
            frame1_df_TE.x,
            frame1_df_TE.Length_of_Avg_V,
        )
        plt.plot(
            frame2_df_TE.x,
            frame2_df_TE.Length_of_Avg_V,
        )
        plt.title(case)
        plt.savefig('test.png')

    frame2_df_TE = frame2_df_TE[frame2_df_TE.x<frame1_df_TE.x.min()]

    frame2_df.x = frame2_df.x - frame2_shift+frame1_df_TE.x.min()
    stitched_frames = frame1_df_TE.append(
        frame2_df_TE
    )
    return stitched_frames.sort('x')


def plot_surface(case,variable='Avg_Vy'):
    from matplotlib import pyplot as plt
    import seaborn as sns

    sns.set(context="notebook", style="whitegrid",
        rc={"axes.axisbelow": False,'image.cmap': 'YlOrRd'})

    df = read_data(root,case,variable)
    X,Y = np.meshgrid(df.x.unique(),df.y.unique())
    Z = df[variable].reshape(X.shape)
    te_location = teloc.TELocations[
        recognize_case(case)[0]
    ][1]

    bl_data,points  = get_bl(case=case,variable=variable)
    delta_BL,vel_BL = find_bl(case=case,variable=variable)
    points   = -points+te_location[0]
    delta_BL = -delta_BL+te_location[0]

    levels = np.linspace(float(Z.min()),float(Z.max())+1,30)
    fig = plt.figure()
    ax = plt.subplot(111,aspect=1)
    ax.contourf(X,Y,Z,levels=levels)
    C = ax.contour(X, Y, Z, levels=levels,
                       colors = ('k',),
                       linewidths = (1,),
              )
    ax.clabel(C, inline=1, fontsize=10,color='w')
    ax.scatter(points,[te_location[1]]*len(points),s=10,color='k')
    ax.scatter(delta_BL,te_location[1],s=40,color='k')
    ax.scatter(delta_BL,te_location[1],marker='x',s=80,color='k')
    plt.savefig('images/Surface_{0}.png'.format(case))
    fig.clear()

def read_tecplot(root,case_folder,variables_to_read):
    """Reads in the tecplot file, stitches the 2 frames
    and returns a pandas data frame with the requested variables

    Input:
        tecplot formated file
    Output:
        pandas data frame

    """
    import os

    for var in variables_to_read:
        tecplot_file = os.path.join(root,case_folder,file_dict[var])
        # Get available variables
        f = open(tecplot_file,'ro')

        variables = []
        # Do two things:
            # 1) Grab the important info from the header
            # 2) See where the second frame info starts so that 
            #       it passes it later to the pandas reader
        var_string = 0
        end_line   = 0
        final_line = 0
        stop_frame_count = False
        for line in f:
            if not stop_frame_count:
                end_line+=1
            if 'Frame 2' in line:
                stop_frame_count = True
            if not var_string:
                var_string = re.findall("^VARIABLES[ _A-Za-z0-9,\"=]+",line)
            if var_string:
                variables = [
                    v.replace(' ','_').replace("\"","") \
                    for v in var_string[0].replace("VARIABLES = ",'').\
                    split(", ")
                ]
                variables = [v for v in variables if len(v)]
            final_line += 1
        f.close()

        lines_to_skip = range(0,3)+range(end_line-1,final_line)

        if var == variables_to_read[0]:
            # Put the first frame data into a data frame
            data_frame1 = pd.read_table(
                    tecplot_file,
                    skiprows  = lines_to_skip,
                    names     = variables,
                    sep       = '[ \t]+',
                    index_col = False,
                    dtype     = np.float
            )

            # Put the second frame data into a data frame
            data_frame2 = pd.read_table(
                    tecplot_file,
                    skiprows  = range(0,end_line),
                    names     = variables,
                    sep       = '[ \t]+',
                    index_col = False,
                    dtype     = np.float
            )
        else:
            # Put the first frame data into a data frame
            df1_tmp = pd.read_table(
                    tecplot_file,
                    skiprows  = lines_to_skip,
                    names     = variables,
                    sep       = '[ \t]+',
                    index_col = False,
                    dtype     = np.float
            )

            # Put the second frame data into a data frame
            df2_tmp = pd.read_table(
                    tecplot_file,
                    skiprows  = range(0,end_line),
                    names     = variables,
                    sep       = '[ \t]+',
                    index_col = False,
                    dtype     = np.float
            )
            if not var in df2_tmp.columns or not var in df1_tmp.columns:
                data_frame1[var] = df1_tmp[davis_dict[var]]
                data_frame2[var] = df2_tmp[davis_dict[var]]
            else:
                data_frame1[var] = df1_tmp[var]
                data_frame2[var] = df2_tmp[var]
        
    # Crop the data
    data_frame1 = data_frame1[
        (data_frame1.x < data_frame1.x.max()*0.90) &\
        (data_frame1.x > data_frame1.x.min()*1.10) &\
        (data_frame1.y < data_frame1.y.max()*0.90) &\
        (data_frame1.y > data_frame1.y.min()*1.10) 
    ]
    data_frame2 = data_frame2[
        (data_frame2.x < data_frame2.x.max()*0.90) &\
        (data_frame2.x > data_frame2.x.min()*1.10) &\
        (data_frame2.y < data_frame2.y.max()*0.90) &\
        (data_frame2.y > data_frame2.y.min()*1.10) 
    ]

    data = stitch_cases(data_frame1, data_frame2,case_folder)

    return data

def pickle_all_data(root,case_name):
    """ Meant to be used only once... pickles the (relevant) 
        TECPLOT data into a single file

    Input:
        TECPLOT file folder
    """

    variables_to_read = [
        'Avg_Vx',
        'Avg_Vy',
        'Length_of_Avg_V',
        'Length_of_Standard_deviation_of_V'
    ]

    df = read_tecplot(raw_data_root,case_name,variables_to_read)

    # Only extract the trailing edge wall normal line
    te_location = teloc.TELocations[
        recognize_case(case_name)[0]
    ][1]
    x = find_nearest(float(te_location[0]),df.x.unique())
    y = find_nearest(float(te_location[1]),df.y.unique())

    df = df[
        (df.y == y) &\
        (df.x < x)
    ]


    df.to_pickle(os.path.join(root,case_name+'.p'))

def find_nearest(to_point,from_array):
   """ Finds the nearest available value in a array to a given value

   Inputs:
      to_point: value to find the nearest to in the array
      from_array: array of available values 
   Returns:
      The nearest value found in the array
   """
   deltas = np.ones(len(from_array))*1000
   for v,i in zip(from_array,range(len(from_array))):
       deltas[i] = abs(float(to_point) - float(v))

   return from_array[np.argmin(deltas)]

def recognize_case(case_name):
    """ Separates the case name folder into its distinctive parameters
    used in this campaign

    Input: folder name
    Output: 
        the key of the TELocations dictionary it belongs to
        case parameters [alpha,side,test_section]
    """

    alpha = int(re.findall('[Aa][0-9][0-9]?',case_name)[0]\
            .replace('A','').replace('a',''))
    try:
        side = re.findall('PS',case_name)[0]
    except:
        side = 'SS'
    try:
        test_section = re.findall('closed',case_name)[0]
    except:
        test_section = 'open'

    # A complicated search for equal terms in the dictionary keys and
    # the case parameters
    # (that's what happens when you don't use standard nomenclature)
    case_key = ''
    for keys,values in zip(
        teloc.TELocations.keys(),
        teloc.TELocations.values()
    ):
        if test_section == values[0]:
            if alpha == int(re.findall('[Aa][0-9][0-9]?',keys)[0]\
               .replace('A','').replace('a','')):
                if alpha: # Hey, alpha = 0 has no side
                    if side == re.findall('[PS]S',keys)[0]:
                        case_key = keys
                        break
                elif alpha==0:
                    case_key = keys
                    break

    case_parameters = [alpha,side,test_section]
    return case_key, case_parameters

def get_bl(case,variable='Length_of_Avg_V'):
    """ Get the TE boundary layer information and return it as an array

    Input:
        tecplot file location
    Output:
        array of flow velocity values
        locations of those velocity vectors
    """

    from numpy import sin,cos,tan
    from math import atan

    root = './Data'
    te_location = teloc.TELocations[
        recognize_case(case)[0]
    ][1]

    df = read_data(root,case,variable)

    # Get the angular value of the flow in the "freestream"
    
    freestream_location_min = df.x.min()*0.95
    freestream_location_max = df.x.min()*0.80

    vy_in_the_freestream = df[
        (df.x < freestream_location_max) & \
        (df.x > freestream_location_min) 
    ].Avg_Vy.mean()

    vx_in_the_freestream = df[
        (df.x < freestream_location_max) & \
        (df.x > freestream_location_min) 
    ].Avg_Vx.mean()

    deviation_angle = atan(vx_in_the_freestream/vy_in_the_freestream)

    df.Avg_Vy = df.Avg_Vy / cos(-deviation_angle)
    df.Avg_Vx = df.Avg_Vy * tan(-deviation_angle)

    bl_data =   np.array(map(float,df.Avg_Vy.values))

    points  = -(np.array(map(float,df['x'].values))-te_location[0])
    df.x = - ( df.x - te_location[0])
    df = df.sort('x')

    df = remove_outliers(df)

    bl_data = moving_average(bl_data,n=50)
    df = get_averaged_data(df,n=50)

    return bl_data,points,df

def remove_outliers(df):
    from scipy import stats
    return df[(np.abs(stats.zscore(df)) < 5).all(axis=1)]

def get_averaged_data(df,n=50):
    import pandas as pd

    variables = df.columns
    new_df = pd.DataFrame(columns=variables)
    for v in variables:
        new_df[v] = moving_average(df[v].values,n=n)

    return new_df

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def plot_bl(case,variable='Avg_Vy'):
    from matplotlib import pyplot as plt
    import seaborn as sns
    current_palette = sns.color_palette()

    x,y,df = get_bl(case,variable)

    fig = plt.figure()
    ax = plt.subplot(111)
    ax.scatter(
        df.Avg_Vy.values/df.Avg_Vy.max(),
        df.x,
        marker='x',
        c=current_palette[0]
        )
    ax.scatter(
        df.Length_of_Standard_deviation_of_V.values/\
        df.Length_of_Standard_deviation_of_V.max(),
        df.x,
        marker='x',
        c=current_palette[1]
        )
    loc_BL,vel_BL = find_bl(case)
    loc_BL = float(loc_BL)
    ax.axhline(y=loc_BL,ls='--',color='r',lw=2)
    ax.text(0.8*ax.get_xlim()[1],
            loc_BL,
            '$\\delta_{{{1}}} ={0:.2f} $ mm'.format(loc_BL,BL_pct),
            ha='right',va='bottom'
           )
    ax.text(0.87*ax.get_xlim()[1],
            df.x.max(),
            '$U ={0:.2f} $ m/s'.format(
                df[df.x==df.x.max()].Avg_Vy.values[0]
            ),
            ha='right',va='center',rotation=-90
           )
    ax.set_xlabel(
        '$U/U_\\mathrm{max}$ (blue), $U_\\mathrm{rms}/U_\\mathrm{rms,max}$ (green)'
    )
    ax.set_ylabel('$y$ [mm]')
    ax.set_ylim(0,45)
    ax.set_xlim(0,1)
    plt.title(case.replace('.p',''))
    plt.savefig('images/BL_{0}.png'.format(case.replace('.p','')))
    fig.clear()
    

def find_bl(case,variable='Length_of_Avg_V'):

    vel,loc,df = get_bl(case,variable)

    vel_BL = df.Avg_Vy.max()*BL_pct/100.
    delta_BL = df[df.Avg_Vy==find_nearest(vel_BL,df.Avg_Vy.values)].x
    #vy_BL = df[
    #    df.Avg_Vy==find_nearest(vel_BL,df.Avg_Vy.values)
    #].Avg_Vy.values[0]
    #vx_BL = df[
    #    df.Avg_Vy==find_nearest(vel_BL,df.Avg_Vy.values)
    #].Avg_Vx.values[0]
    #for v,l in zip(vel[::-1],loc[::-1]):
    #    if v>vel_BL:
    #        delta_BL = l
    #        break

    return delta_BL,vel_BL

def make_csv(out_file="BL_Data_Info.csv"):
    import pandas as pd
    from re import findall
    from os.path import join

    bl_info_DF = pd.DataFrame(
        columns = [
            'U_inf',
            'alpha',
            'Delta_BL',
            'U_BL',
            'Side',
            'Test_section'
        ])

    variable = 'Length_of_Avg_V'
    for case in data_folders:
        alpha = findall("_[aA][0-9][0-9]?",case)[0].\
                replace("_a","").\
                replace("_A","")
        U_inf = findall("_U[0-9][0-9]?",case)[0].replace("_U","")
        if len(findall("[PS]S",case)):
            side  = findall("[PS]S",case)[0]
        else: side = "NA"
        if len(findall("closed",case)):
            test_section  = findall("closed",case)[0]
        else: test_section = "open"
        Delta_BL,U_BL = find_bl(case,variable=variable)

        bl_info_DF = bl_info_DF.append({
            'U_inf'        : float(U_inf),
            'alpha'        : float(alpha),
            'Delta_BL'     : float(Delta_BL),
            'U_BL'         : float(U_BL),
            'Side'         : side,
            'Test_section' : test_section
        },ignore_index=True)

    if out_file:
        bl_info_DF.to_csv(join("outputs",out_file))
    return bl_info_DF

def plot_all_deltas(out_file="All_Deltas.png"):
    import matplotlib.pyplot as plt 
    import seaborn as sns
    from numpy import argmin,abs
    import os
    sns.__version__

    bl_info_DF = make_csv(out_file='')
    bl_info_DF = bl_info_DF.sort("U_inf")

    cmap_SS = sns.color_palette("Reds_r",len(bl_info_DF.U_inf.unique()))
    cmap_PS = sns.color_palette("Blues_r",len(bl_info_DF.U_inf.unique()))

    # Marker defines if open or closed
    marker = {
        'open'   : 'o',
        'closed' : 'x',
    }
    fig,ax = plt.subplots(1,1)
    for row_index, row in bl_info_DF.iterrows():
        if row.Side == 'SS' or row.Side == 'NA':
            cmap = cmap_SS
        else:
            cmap = cmap_PS

        if row.Side=='NA':
            label_side = ''
        else: label_side = row.Side
        label = "$U_\\infty = {0}$ m/s {1}".\
        format(row.U_inf,label_side)

        ax.scatter(
            row.alpha,
            row.Delta_BL,
            color = cmap[
                argmin(abs(bl_info_DF.U_inf.unique()-row.U_inf))
            ],
            marker = marker[row.Test_section],
            s=75,
            label = label
        )
        ax.annotate(label,
                    xy=(row.alpha,row.Delta_BL),
                    xytext=(row.alpha+1,row.Delta_BL+0.5),
                    arrowprops=dict(
                        facecolor='black', 
                        arrowstyle='->',
                        lw=2
                    ),
                   )
    ax.set_xticks([0,6,12])
    ax.set_xlim(-2,16)
    ax.set_xlabel("$\\alpha_g$ [deg]")
    ax.set_ylabel("$\\delta_{{{0}}}$ [mm]".format(BL_pct))
    #ax.legend(loc='best')
    plt.savefig(os.path.join('images',out_file))


#variable = 'Length_of_Standard_deviation_of_V'
##variable = 'Length_of_Avg_V'
#import os
#variables_to_read = [
#    'Avg_Vx',
#    'Avg_Vy',
#    'Length_of_Avg_V',
#    'Length_of_Standard_deviation_of_V'
#]
#
##for case in [raw_data_folders[2]]:
##    read_tecplot(raw_data_root,case,variables_to_read)
##    pickle_all_data(raw_data_root,case)
##for case in data_folders:
##    plot_bl(case,variable)
#    #plot_surface(case,variable)
plot_all_deltas()
#
#def plot_all_BLs():
#    for case in data_folders:
#        plot_bl(case,variable)
#
def pickle():
    for case in raw_data_folders:
        pickle_all_data(raw_data_root,case)
#
#pickle()
#plot_all_BLs()


