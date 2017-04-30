import os
import pandas as pd
import TELocations as teloc
import numpy as np
import re

file_dict = {
    'Avg_Vx':                            "B00001.dat",
    'Avg_Vy':                            "B00002.dat",
    'Avg_Vz':                            "B00003.dat",
    'Length_of_Avg_V':                   "B00004.dat",
    'RMS_Vx':                            "B00005.dat",
    'RMS_Vy':                            "B00006.dat",
    'Length_of_Standard_deviation_of_V': "B00007.dat",
    'Turbulent_kinetec_energy':          "B00008.dat",
    'Reynold_stress_XY':                 "B00009.dat",
    'Reynold_stress_XX':                 "B00010.dat",
    'Reynold_stress_YY':                 "B00011.dat",
}

davis_dict = {
    'Avg_Vx':                            'Avg_Vx',
    'Avg_Vy':                            'Avg_Vy',
    'Length_of_Avg_V':                   'Length_of_Avg_V',
    'Standard_deviation_of_Vx':          'RMS_Vx',
    'Standard_deviation_of_Vy':          'RMS_Vy',
    'RMS_Vx':                            'RMS_Vx',
    'RMS_Vy':                            'RMS_Vy',
    'Length_of_Standard_deviation_of_V': 'Length_of_RMS_V',
    'Turbulent_kinetec_energy':          'Turbulent_kinetec_energy',
    'Reynold_stress_XY':                 'Reynold_stress_XY',
    'Reynold_stress_XX':                 'Reynold_stress_XX',
    'Reynold_stress_YY':                 'Reynold_stress_YY',
}

need_to_stitch_cases = [
    'STE_A6_U20_closed_SS',
    'STE_SS_a12_U20',
    'STE_SS_a12_U30',
    'STE_SS_a12_U35',
    'STE_SS_a12_U40',
]

BL_pct = 99

root = '/home/carlos/Documents/PhD/Articles/Article_3/Scripts/acoustics_and_bl/PlaneData'
#raw_data_root = '/media/carlos/6E34D2CD34D29783/2015-07_BL/STE_BL_Data/'
raw_data_root = '/home/carlos/Documents/PhD/Articles/Article_3/Scripts/acoustics_and_bl/BLData_new'

data_folders = sorted( [os.path.join( root, f ) for f in os.listdir(root) \
                if os.path.isfile(os.path.join(root,f))
               and f.endswith(".p")] )
try:
    raw_data_folders = [f for f in os.listdir(raw_data_root) \
                    if os.path.isdir(os.path.join(raw_data_root,f))]
except OSError:
    raw_data_folders = []

def rename_DaVis_case_folders( raw_data_root = raw_data_root ):
    from re      import findall
    from shutil  import move
    from os      import listdir
    from os.path import join

    case_folders = [ f for f in listdir( raw_data_root ) ]

    for cf in case_folders:
        if 'closed' in cf: closed = True
        else: closed = False
        if 'PS' in cf: side = 'PS'
        elif 'SS' in cf: side = 'SS'
        else: side = ''

        try:
            base_name = findall( 
                'STE_[S]?[P]?[S]?[_]?[aA][0-9][0-9]?_[p]?[0]?' + \
                '[_]?[S]?[S]?[P]?[_]?U[0-9][0-9]', 
                cf 
            )[0]
            base_name = base_name.\
                    replace( 'PS_', '' ).replace( 'SS_', '' )
            if side:
                base_name = base_name + '_' + side
            if closed:
                base_name = base_name + '_closed'
            print base_name
            move( 
                join( raw_data_root, cf ), 
                join( raw_data_root, base_name )
            )
        except IndexError:
            print cf


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

def stitch_cases_2( frame1, frame2, case ):
    from numpy             import min, max, mgrid
    from scipy.interpolate import griddata
    from pandas            import DataFrame, concat

    stitch_dict = {
        #'STE_a12_U20_SS' : ( 16.5, 3.3 )
        'STE_a12_U20_SS' : ( 15.7, 3.0),
        'STE_a12_U30_SS' : ( 16.5, 3.0),
        'STE_a12_U35_SS' : ( 17.0, 1.5),
        'STE_a12_U40_SS' : ( 17.0, 2.0),
    }

    frame1 = frame1.copy()[ frame1.x > frame1.x.min() * 1.05 ]
    frame2 = frame2.copy()[ frame2.x < frame2.x.max() * 0.95 ]

    print '   Stitching by translating frame 2 by:'
    print '     x = {0}'.format( stitch_dict[ case ][0] )
    print '     y = {0}'.format( stitch_dict[ case ][1] )
    frame2.x = frame2.x - stitch_dict[ case ][0]
    frame2.y = frame2.y - stitch_dict[ case ][1]

    lower_x_limit = min( [ frame2.x.min(), frame1.x.min() ] )
    upper_x_limit = max( [ frame2.x.max(), frame1.x.max() ] )
    lower_y_limit = max( [ frame2.y.min(), frame1.y.min() ] )
    upper_y_limit = min( [ frame2.y.max(), frame1.y.max() ] )


    grid_x, grid_y = mgrid[ 
        lower_x_limit  : upper_x_limit  : 200j,
        lower_y_limit  : upper_y_limit  : 100j
    ]

    frame2 = frame2[ frame2.x <= frame1.x.min() ]

    stitched_frames = frame1.append(
        frame2
    )
    stitched_frames = stitched_frames[ stitched_frames.Avg_Vx != 0 ]

    interpolated_df = DataFrame()

    for var in [v for v in stitched_frames.columns if v != 'x' and v != 'y' ]:

        interpolated_df = concat( 
            [ interpolated_df, 
             DataFrame( data = {
                 'x' : grid_x.ravel(),
                 'y' : grid_y.ravel(),
                 var : griddata(
                     ( stitched_frames.x.values, stitched_frames.y.values ),
                     stitched_frames[var].values, 
                     (grid_x, grid_y), 
                     method='nearest'
                 ).ravel()
             } )
            ], axis = 1
        )

    interpolated_df = interpolated_df.T.drop_duplicates().T

    interpolated_df = interpolated_df\
            .sort_values( by = ['y','x'] )\
            .dropna()\
            .reset_index( drop = True )

    return interpolated_df

def stitch_cases(frame1_df,frame2_df,case,plot=False):
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.__version__

    stitching_parameter = 'RMS_Vy'

    frame1_df_TE = frame1_df 
    frame2_df_TE = frame2_df 

    near_frame_max = frame1_df_TE[ stitching_parameter ].max()

    frame2_shift = frame2_df_TE[
        frame2_df_TE[ stitching_parameter ] == find_nearest(
            near_frame_max,
            frame2_df_TE[ stitching_parameter ].values)
    ].x.values[-1]

    frame2_df_TE = frame2_df_TE[frame2_df_TE.x<=frame2_shift]
    frame2_df_TE.x = frame2_df_TE.x - \
            frame2_shift + frame1_df_TE.x.min()
    
    if plot:
        fig = plt.figure()

        plt.plot(
            frame1_df_TE.x,
            frame1_df_TE[ stitching_parameter ],
        )
        plt.plot(
            frame2_df_TE.x,
            frame2_df_TE[ stitching_parameter ],
        )
        plt.title(case)
        plt.savefig('PlaneData/{0}_Stitched.png'.format( case ) )
        plt.close( fig )

    frame2_df_TE = frame2_df_TE[frame2_df_TE.x<frame1_df_TE.x.min()]

    frame2_df.x = frame2_df.x - frame2_shift+frame1_df_TE.x.min()
    stitched_frames = frame1_df_TE.append(
        frame2_df_TE
    )
    return stitched_frames.sort_values( by = ['y', 'x'] )

def extract_case_details_from_name(case_name):
    """ Returns a dictionary of the case details that are extracted
    from the passed name

    Input: case data file name
    Output: a dictionary
    """
    from re import findall

    alpha = int(findall("[Aa][0-9][0-9]?",case_name)[0]\
            .replace('A','')\
            .replace('a',''))
    speed = int(findall("U[0-9][0-9]",case_name)[0].replace("U",''))
    if 'PS' in case_name or "SS" in case_name:
        side  = findall("[PS]S",case_name)[0]
    else:
        side = "zero alpha"
    if "closed" in case_name:
        test_section = 'closed'
    elif "open" in case_name:
        test_section = 'open'
    else:
        test_section = 'unknown'

    case_details = {
        'alpha'        : alpha,
        'speed'        : speed,
        'side'         : side,
        'test_section' : test_section
    }
    return case_details

def build_plot_case_label_from_dict(case_details):
    """ Turns the case details dictionary and returns a label that
    can be used for plotting

    Input: case details dictionary
    Output: label
    """

    label = "$\\alpha = {{{0}}}^\\circ$, $U_\\infty = {{{1}}}\\,\\mathrc{{m/s}}$".\
            format(case_details['alpha'],case_details['speed'])
    return label


def plot_article_bls(fig_name = 'BoundaryLayers.png'):
    from matplotlib import pyplot as plt
    import pandas as pd
    import seaborn as sns
    from matplotlib import rc
    from numpy import argmin,array
    rc('text',usetex=True)

    sns.set_context('paper')
    sns.set_style("whitegrid")
    sns.set(font='serif',font_scale=2.5,style='whitegrid')
    rc('font',family='serif', serif='cm10')

    #line_styles = ['--','-.',':']
    markers = [
            u'o', u'v', u'^', u'<', u'>', u'8', u's', u'p', u'*', 
        u'h', u'H', u'D', u'd'
    ]

    colors = [
            '#013F70',
            '#70A288',
            '#D5896F',
            '#BB9F06',
            '#DAB785',
            ]

    bls_df = pd.DataFrame()
    for case in data_folders:
        x,y,df = get_bl(case,variable='Avg_Vy')
        df['case'] = case
        bls_df = bls_df.append(df)

    cases = bls_df.case.unique()
    alphas = array([0,6,12])

    fig,ax = plt.subplots(1,1)
    for case in [c for c in cases \
                 if not 'closed' in c and not 'PS' in c]:

        case_dict = extract_case_details_from_name( case )       

        data = bls_df[bls_df.case==case]
        vel  = data.Avg_Vy / data.Avg_Vy.max()
        y    = data.x

        ax.scatter(vel,y,
                   label = build_plot_case_label_from_dict(case_dict),
                   marker = markers[
                       argmin(abs(
                           case_dict['alpha'] - alphas
                       ))
                   ],
                   color = colors[
                       argmin(abs(
                           case_dict['alpha'] - alphas
                       ))
                   ],
                   alpha = 0.2
                  )

    plt.savefig(fig_name)

    return bls_df

def plot_surface(pickled_case, variable='u', plot_name = ''):
    from matplotlib import pyplot as plt
    import seaborn as sns
    from pandas import read_pickle
    from numpy import linspace

    sns.set(context="notebook", style="whitegrid",
        rc={"axes.axisbelow": False,'image.cmap': 'YlOrRd'})

    try:
        df = read_pickle( pickled_case )
    except TypeError:
        df = pickled_case

    df = df.drop_duplicates()

    df_pivot = df.pivot( 'y', 'x', variable )
    X   = df_pivot.columns.values
    Y   = df_pivot.index.values
    x,y = np.meshgrid(X, Y)
    Z   = df_pivot.values
    df_pivot = df.pivot( 'y', 'x', 'u' )
    U   = df_pivot.values
    df_pivot = df.pivot( 'y', 'x', 'v' )
    V   = df_pivot.values

    #levels = np.linspace( float(Z.min()), float(Z.max())+1,30 )

    fig, ax  = plt.subplots( 1, 2 )

    if variable == 'u':
        levels = linspace( 0, df.u.max(), 10 )
    elif variable == 'v':
        levels = linspace( 0, df.v.max(), 13 )
    else:
        levels = linspace( -250000,10000 , 30 )
    ax[0].contourf(x, y, Z,  levels=levels)
    #C = ax.contour(
    #    X, Y, Z, 
    #    levels     = levels,
    #    colors     = ('k',),
    #    linewidths = (1,),
    #)

    stride = 10
    ax[0].quiver( 
        x[::stride, ::stride], 
        y[::stride, ::stride], 
        U[::stride, ::stride], 
        V[::stride, ::stride],
    )

    ax[1].plot( 
        df[ df.x == find_nearest( 0, df.x ) ].u,
        df[ df.x == find_nearest( 0, df.x ) ].y,
    )

    ax[0].set_aspect( 'equal' )

    #ax.clabel(C, inline=1, fontsize=10, color='w')

    if not plot_name:
        plt.show()
    else:
        plt.savefig( plot_name ) 
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
        tecplot_file = os.path.join( 
            root, 
            case_folder, 
            file_dict[ davis_dict[ var ] ] 
        )
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
                end_line += 1
            if 'Frame 2' in line:
                stop_frame_count = True
            if not var_string:
                var_string = re.findall(
                    "^VARIABLES[ _A-Za-z0-9,\"=]+",line
                )
            if var_string:
                variables = [
                    v.replace( ' ', '_' ).replace( "\"", "" ) \
                    for v in var_string[0].replace( "VARIABLES = ", '' ).\
                    split( ", " )
                ]
                variables = [v for v in variables if len(v)]
            final_line += 1
        f.close()

        lines_to_skip = range( 0, 3 ) + range( end_line-1, final_line )

        if not variables[-1] in variables_to_read:
            variables[-1] = davis_dict[ variables[-1] ]
            if not variables[-1] in variables_to_read:
                continue
        if var == variables_to_read[0]:
            # Put the first frame data into a data frame
            data_frame1 = pd.read_table(
                tecplot_file,
                skiprows  = lines_to_skip,
                names     = variables,
                sep       = '[ \t]+',
                index_col = False,
                engine    = 'python'
            )

            # Put the second frame data into a data frame
            data_frame2 = pd.read_table(
                tecplot_file,
                skiprows  = range(0,end_line),
                names     = variables,
                sep       = '[ \t]+',
                index_col = False,
                engine    = 'python'
            )
        else:
            # Put the first frame data into a data frame
            df1_tmp = pd.read_table(
                tecplot_file,
                skiprows  = lines_to_skip,
                names     = variables,
                sep       = '[ \t]+',
                index_col = False,
                engine    = 'python'
            )

            # Put the second frame data into a data frame
            df2_tmp = pd.read_table(
                tecplot_file,
                skiprows  = range(0,end_line),
                names     = variables,
                sep       = '[ \t]+',
                index_col = False,
                engine    = 'python'
            )
            if not var in df2_tmp.columns or not var in df1_tmp.columns:
                data_frame1[ var ] = df1_tmp[ davis_dict[ var ] ]
                data_frame2[ var ] = df2_tmp[ davis_dict[ var ] ]
            else:
                data_frame1[ var ] = df1_tmp[ var ]
                data_frame2[ var ] = df2_tmp[ var ]
        
    ## Crop the data
    #data_frame1 = data_frame1[
    #    ( data_frame1.x < data_frame1.x.max() * 0.95 ) &\
    #    ( data_frame1.x > data_frame1.x.min() * 1.05 ) &\
    #    ( data_frame1.y < data_frame1.y.max() * 0.95 ) &\
    #    ( data_frame1.y > data_frame1.y.min() * 1.05 ) 
    #].reset_index( drop = True )
    #data_frame2 = data_frame2[
    #    (data_frame2.x < data_frame2.x.max()*0.90) &\
    #    (data_frame2.x > data_frame2.x.min()*1.10) &\
    #    (data_frame2.y < data_frame2.y.max()*0.90) &\
    #    (data_frame2.y > data_frame2.y.min()*1.10) 
    #]


    to_stitch_cases = [
        'STE_a12_U20_SS',
        'STE_a12_U30_SS',
        'STE_a12_U35_SS',
        'STE_a12_U40_SS',
    ]

    if case_folder in to_stitch_cases:
        data = stitch_cases_2( data_frame1, data_frame2, case_folder )
    else:
        data = data_frame1
    return data

def rename_df_columns_from_DaVis_to_standard(df):
    DaVis_naming_dict= {
          "x":               "x",
          "y":               "y",
          "Avg_Vx":          "u",
          "Avg_Vy":          "v",
          "Avg_Vz":          "w",
          "Length_of_Avg_V": "U",
          "RMS_Vx":          "urms",
          "RMS_Vy":          "vrms",
          }

    for c in df.columns:
        if not c in DaVis_naming_dict.keys():
            df = df.drop( c, 1 )

    df.columns = [
        DaVis_naming_dict[col] for col in df.columns
    ]

    return df

def pickle_all_data(root,case_name,output_folder):
    """ Meant to be used only once... pickles the (relevant) 
        TECPLOT data into a single file

    Input:
        TECPLOT file folder
    """
    from os.path import join

    variables_to_read = [
        'Avg_Vx',
        'Avg_Vy',
        'RMS_Vx',
        'RMS_Vy',
        'Length_of_Avg_V',
        #'Length_of_Standard_deviation_of_V'
    ]

    print case_name
    df = read_tecplot( raw_data_root, case_name, variables_to_read )

    df = rename_df_columns_from_DaVis_to_standard( df ) 

    x = df.x.copy()
    y = df.y.copy()
    u = df.u.copy()
    v = df.v.copy()

    te_location = teloc.TELocations[
        recognize_case( case_name )[0]
    ][1]

    df.x =   y - te_location[1]
    df.y = - x + te_location[0]
    df.u =   v
    df.v = - u

    pickle_name = join( output_folder , case_name+'.p' )

    df.to_pickle( pickle_name )

    plot_surface( pickle_name, 
                 plot_name = join( output_folder , case_name+'.png' ) ,
                 variable = 'v'
                )


def find_nearest(to_point,from_array):
   """ Finds the nearest available value in a array to a given value

   Inputs:
      to_point: value to find the nearest to in the array
      from_array: array of available values 
   Returns:
      The nearest value found in the array
   """
   deltas = np.ones( len( from_array ) ) * 1000
   for v,i in zip( from_array, range( len( from_array ) ) ):
       deltas[i] = abs( float( to_point ) - float( v ) )

   return from_array[ np.argmin( deltas ) ]

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

def remove_outliers(df):
    #def reject_outliers(data, m=2):
    #    return data[abs(data - np.mean(data)) < m * np.std(data)]
    from scipy import stats
    return df[(np.abs(stats.zscore(df)) < 3).all(axis=1)]
    #df.Avg_Vy = df.Avg_Vy

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

def plot_bl(case, variable='Avg_Vy'):
    from matplotlib import pyplot as plt
    import seaborn as sns
    from pandas import read_pickle
    from numpy import abs, isnan
    from os.path import split

    current_palette = sns.color_palette()

    df = read_pickle( case )

    df = get_vorticity_2( df )

    df = df[ df.x == find_nearest( 0 , df.x.values ) ]

    loc_BL, delta_displacement, delta_momentum, vel_BL, U_max = find_bl( 
        read_pickle( case ), 
        split( case )[1].replace('.p','') 
    )

    if isnan( loc_BL ):
        return 0

    fig = plt.figure()
    ax  = plt.subplot(111)

    ax.scatter(
        df.u.values / U_max,
        df.y,
        marker='x',
        c=current_palette[0]
        )

    ax.scatter(
        abs( df.vorticity_xy.values/\
        df.vorticity_xy.mean() ),
        df.y,
        marker='+',
        c=current_palette[1]
        )

    ax.axhline(y=loc_BL,ls='--',color='r',lw=2)
    ax.text(0.8*ax.get_xlim()[1],
            loc_BL,
            '$\\delta_{{{1}}} ={0:.2f} $ mm'.format(loc_BL,BL_pct),
            ha='right',va='bottom'
           )
    ax.text(0.87*ax.get_xlim()[1],
            df.x.max(),
            '$U ={0:.2f} $ m/s'.format(
                df[df.x==df.x.max()].v.values[0]
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
    plt.savefig('images/{0}.png'.format( 
        os.path.split(case)[1].replace('.p','')) 
    )
    fig.clear()
    
def get_edge_velocity(
    df, 
    #condition = 'vorticity_integration_rate_of_change',
    condition = 'u_rate_of_change',
    threshold = 0.10,
    Ue_file   = 'Ue_data.csv',
    plot      = False
):
    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy.integrate import simps
    from numpy           import diff, nan, abs
    from copy            import copy
    from pandas          import rolling_mean

    def approximate_Uedge( df, threshold, condition ):

        getting_close = 0
        U_edge        = 0
        y_bl          = 0

        for y_loc in df.iterrows():
            #if abs( y_loc[1][condition] / df[ condition ].max() ) < threshold :

            if abs( y_loc[1][condition] ) < threshold :
                getting_close += 1

                if getting_close == 1:
                    U_edge = y_loc[1]['u']
                    y_bl   = y_loc[1]['y']

                if getting_close > 3 and y_loc[1].y > 5 :
                    return U_edge, y_bl

            else:
                getting_close = 0

        return U_edge, y_bl

    sns.__version__
    
    U_edge = 0 

    df = df.sort_values( by = 'y', ascending = True )
    df = df[ df.y > 2 ]
    df = df.reset_index( drop = True )

    df = wallnormal_average_data( df )
    df = df.dropna().reset_index( drop = True )

    # Get the integrated vorticity #############################################
    integration_result = []
    for i in range(len(df.y)):
        integration_result.append(
            - simps(
                df.vorticity_xy.ix[:i], 
                x    = df.y.ix[:i] / 1000.,
                even = 'avg'
                 )
        )
    if not len( integration_result ):
        return nan

    df['vorticity_integration'] = integration_result

    rate_of_change = diff( df.vorticity_integration ) / diff( df.y  )
    df = df.sort_values( by = 'y' ).reset_index( drop = True )
    if not len( rate_of_change ):
        return nan
    rate_of_change = list(rate_of_change) +  [rate_of_change[-1]]
    df['vorticity_integration_rate_of_change'] = abs( rate_of_change )

    rate_of_change = diff( df.u ) / diff( df.y  )
    df = df.sort_values( by = 'y' ).reset_index( drop = True )
    if not len( rate_of_change ):
        return nan
    rate_of_change = list(rate_of_change) +  [rate_of_change[-1]]
    df['u_rate_of_change'] = abs( rate_of_change )

    df = rolling_mean( df , window = 5 )

    df = df.dropna().reset_index( drop = True )

    threshold0    = copy( threshold )
    
    U_edge, y_bl = approximate_Uedge( df, threshold, condition )

    while not U_edge:
        threshold += threshold0
        if threshold > 2:
            print "      Could not find a U_e for x = {0:.2f}".format(
                df.x.unique()[0]
            )
            return nan
        U_edge, y_bl = approximate_Uedge( df, threshold, condition )
    
    if plot:
        fig, ax = plt.subplots(1, 1)
        ax.plot( df['y'], df[condition], )
        ax.axhline( threshold ) 
        ax2 = ax.twinx()
        ax2.plot( df['y'], df['u'] / U_edge, )
        ax2.plot( df['y'], df['urms'] / df.urms.max(), )
        ax2.axvline( y_bl, ls = '--' , c = 'k' ) 
        ax2.axvline( 
            df[ df.u == find_nearest( U_edge * 0.99, df.u.values )].y.values[0],
            ls = '-' , c = 'k' 
        ) 
        ax2.set_ylim( 0, 1.05 ) 
        plt.savefig( plot )


    return U_edge, y_bl

def find_bl( df , case_name = '', U_e = 0):

    def get_delta_99( df, U_e, y_e ):
        from numpy import isnan,nan

        if isnan(U_e):
            return nan

        y_bl = df[ df.y < y_e ]
        for y in y_bl.y.unique():
            if y_bl[ y_bl.y == y ].u.values[0] >= 0.99 * U_e:
                break

        return y

    def get_delta_momentum( df, U_e ):
        from scipy.integrate import simps
        from numpy import isnan,nan
        if isnan(U_e):
            return nan
        return simps(
            ( df.u / U_e ) * ( 1 - df.u / U_e ), 
            x = df.y , 
            even='avg'
        )

    def get_delta_displacement( df, U_e ):
        from scipy.integrate import simps
        from numpy import isnan,nan
        if isnan(U_e):
            return nan
        return simps(
            1 - df.u / U_e , 
            x = df.y , 
            even = 'avg'
        )

    #write_tecplot( 'test', df, 'Concat_tecplots.dat', interpolate=False )

    df = df[ df.y > 0 ]
    df = df[ df.y < df.y.max() * 0.95 ]
    df = df.sort_values( by = 'y' )

    if not U_e:

        df = df.dropna().reset_index( drop = True )
        df = get_vorticity_2( df )
        df = df.dropna().reset_index( drop = True )
        df = streamwise_average_data( df )
        #plot_surface( df , 'vorticity_xy' )

        df = df.dropna().reset_index( drop = True )
        df_x0 = df[ df.x == find_nearest( 0, df.x.values ) ]
        U_e, y_e = get_edge_velocity( 
            df_x0 ,
            plot = 'images/{0}_Profile.png'.format( case_name )
        )


    delta_99           = get_delta_99( df_x0, U_e, y_e )
    df_x0 = df_x0[ ( df_x0.y > 1.5 ) & ( df_x0.y < 1.1 * delta_99 ) ]
    delta_displacement = get_delta_displacement(df_x0,U_e)
    delta_momentum     = get_delta_momentum(df_x0,U_e)

    print '\t\t',
    print "Ue = {0:.2f}\tdelta_99 = {1:.2f}\tdelta* = {2:.2f}\ttheta = {3:.2f}"\
            .format( U_e, delta_99, delta_displacement, delta_momentum )

    return delta_99, delta_displacement, delta_momentum, U_e*0.99, U_e

def streamwise_average_data( df ):
    from pandas import rolling_mean

    for y in df.y.unique():
        df.loc[ df.y == y ] = rolling_mean( df[ df.y == y ] , window = 10 )

    return df

def wallnormal_average_data( df ):
    from pandas import rolling_mean

    for x in df.x.unique():
        df.loc[ df.x == x, 'vorticity_xy' ] = \
                rolling_mean( df[ df.x == x ].vorticity_xy , window = 10 )

    return df

def get_vorticity(df):
    from numpy import shape,zeros
    from sys import exit

    if "vorticity_xy" in df.columns:
        # Do nothing and return the same DF
        return df

    df = df.sort_values( by = ['y','x'] )
    # Get shape of 2D plane
    nx = len(df['x'].unique())
    ny = len(df['y'].unique())
    Ux = df['u'].values.reshape((ny,nx))
    Uy = df['v'].values.reshape((ny,nx))
    ax = df['x'].values.reshape((ny,nx))/1000. # [mm] -> [m]
    ay = df['y'].values.reshape((ny,nx))/1000. # [mm] -> [m]

    i,j = shape(Ux)

    # Sanity check:
    if i != shape(Uy)[0] or i != shape(ax)[0] or i != shape(ay)[0]:
        exit("   The shape of the matrices while getting the "+\
             "vorticity is not the same!")

    duy_dax = zeros((i,j))
    dux_day = zeros((i,j))
    for ii in range(1,i-1):
        for jj in range(1,j-1):
            duy_dax[ii,jj] = (Uy[ii,jj+1]-Uy[ii,jj-1])\
                    /(ax[ii,jj+1]-ax[ii,jj-1])
    for ii in range(1,i-1):
        for jj in range(1,j-1):
            dux_day[ii,jj] = (Ux[ii+1,jj]-Ux[ii-1,jj])\
                    /(ay[ii+1,jj]-ay[ii-1,jj])

    vorticity = duy_dax - dux_day

    df['vorticity_xy'] = vorticity.ravel()

    return df

def get_vorticity_2( df ):
    from numpy import gradient, diff

    df = df.sort_values( by = ['x','y'] )
    df_pivot = df.pivot( 'x', 'y', 'u' )
    U   = df_pivot.values
    df_pivot = df.pivot( 'x', 'y', 'v' )
    V   = df_pivot.values

    dx = diff( df.x ).mean() / 1000.
    dy = diff( df.y ).mean() / 1000.

    dv_dy, dv_dx = gradient( V, dy, dx )
    du_dx, du_dy = gradient( U, dx, dy )

    vorticity = dv_dx - du_dy

    df[ 'vorticity_xy' ] = vorticity.ravel()

    return df


def make_csv( out_file="BL_Data_Info.csv" ):
    import pandas as pd
    from re      import findall
    from os.path import join, split

    bl_info_DF = pd.DataFrame(
        columns = [
            'U',
            'alpha',
            'Delta_BL',
            'U_BL',
            'Side',
            'Test_section'
        ])

    #variable = 'Length_of_Avg_V'
    for case in data_folders:
        print "    Processing {0}".format( case )
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

        Delta_BL, delta_disp, delta_momentum, U_BL, U_max = find_bl( 
            pd.read_pickle( case ), split( case )[1].replace('.p','')
        )

        bl_info_DF = bl_info_DF.append({
            'U':                  float(U_inf),
            'U_max':              U_max,
            'alpha':              float(alpha),
            'Delta_BL':           float(Delta_BL),
            'Delta_displacement': float(delta_disp),
            'Delta_momentum':     float(delta_momentum),
            'U_BL':               float(U_BL),
            'Side':               side,
            'Test_section':       test_section
        },ignore_index=True)

    print bl_info_DF
    if out_file:
        bl_info_DF.to_csv(join("outputs",out_file))
        bl_info_DF.to_pickle(join("outputs",out_file.replace('.csv','.p')))
    return bl_info_DF

def plot_article_deltas( out_file = "Article_Deltas.png" , overwrite = False ):
    import matplotlib.pyplot as plt 
    import matplotlib.lines as mlines
    import seaborn as sns
    import os
    from os.path import isfile
    from pandas import read_pickle
    from matplotlib import rc
    sns.__version__

    rc('text',usetex=True)
    sns.set_context('paper')
    sns.set_style("whitegrid")
    sns.set(font='serif',font_scale=2.0,style='whitegrid')
    rc('font',family='serif', serif='cm10')

    bl_pickle = 'outputs/BL_Data_Info.p'
    if isfile( bl_pickle ) and not overwrite:
        bl_info_DF = read_pickle( bl_pickle )
    else:
        bl_info_DF = make_csv(out_file='')

    markeredgewidth = 3
    markerfacecolor = 'none'
    markersize      = 12
    mew             = 4 # Marker edge width

    bl_info_DF = bl_info_DF[ bl_info_DF.Test_section == 'open' ]

    palette = sns.color_palette("deep")

    # Marker defines if open or closed
    marker_dict = {
        '20': '+',
        '30': 'x',
        '35': '1',
        '40': '2',
    }

    marker_legend = [ mlines.Line2D(
        range(1),range(1),
        color='w', 
        marker          = marker_dict['20'],
        markersize      = markersize,
        markeredgecolor = 'k',
        markerfacecolor = 'none',
        mew             = 3,
        label           = '20 m/s'
    ) ]
    marker_legend.append( mlines.Line2D(
        range(1),range(1),
        color           = 'w',
        marker          = marker_dict['30'],
        label           = '30 m/s',
        markersize      = markersize,
        markeredgecolor = 'k',
        markerfacecolor = 'none',
        mew             = 3,
    ) )
    marker_legend.append( mlines.Line2D(
        range(1),range(1),
        color           = 'w',
        marker          = marker_dict['35'],
        label           = '35 m/s',
        markersize      = markersize,
        markeredgecolor = 'k',
        markerfacecolor = 'none',
        mew             = 3,
    ) )
    marker_legend.append( mlines.Line2D(
        range(1),range(1),
        color           = 'w',
        marker          = marker_dict['40'],
        label           = '40 m/s',
        markersize      = markersize,
        markeredgecolor = 'k',
        markerfacecolor = 'none',
        mew             = 3,
    ) )

    fig,ax = plt.subplots( 1, 3, figsize = ( 10, 4 ) )

    for row_index, row in bl_info_DF.iterrows():
        if row.Side == 'NA':
            color = palette[0]
        elif row.Side == 'SS':
            color = palette[1]
        elif row.Side == 'PS':
            color = palette[2]

        marker = marker_dict[ "{0}".format( int( row.U ) ) ]

        plot_config = {
            'ls'              : '',
            'marker'          : marker,
            'markeredgewidth' : markeredgewidth,
            'markerfacecolor' : markerfacecolor,
            'markeredgecolor' : color,
            'markersize'      : markersize,
            'mew'             : mew,
            'color'           : color,
        }

        ax[0].plot(
            row.alpha,
            row.Delta_BL,
            **plot_config
        )

        ax[1].plot(
            row.alpha,
            row.Delta_displacement,
            **plot_config
        )

        ax[2].plot(
            row.alpha,
            row.Delta_momentum,
            **plot_config
        )
    for a in ax:
        a.set_xlabel("$\\alpha$")
        a.set_xticks([0,6,12])
        a.set_xticklabels([r'$0^\circ$',r'$6^\circ$',r'$12^\circ$'])
        a.set_xlim(-2,14)
    ax[0].set_ylabel("$\\delta_{{{0}}}$ [mm]".format(BL_pct))
    ax[1].set_ylabel("$\\delta^*$ [mm]".format(BL_pct))
    ax[2].set_ylabel("$\\theta$ [mm]".format(BL_pct))
    #ax.legend(loc='best')
    ax[0].legend(handles = marker_legend,
               bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=4, borderaxespad=0., numpoints = 1
              )
    plt.tight_layout()
    plt.savefig(os.path.join('article_images',out_file), bbox_inches='tight')

def write_tecplot(case_name,df,outfile,interpolate=False):
   from numpy import nan
   import os

   # Construct the TECPLOT header for the combined data file
   header1 = "TITLE = \"test\""
   header2 = "VARIABLES = "
   for var in df.columns:
           header2 = header2+"\"{0}\"".format(var)
           if var != df.columns[-1]:
                   header2 = header2+", "
   header3 = "ZONE T=\"Frame 1\" I={0}, J={1}, F=POINT\n".\
           format(len(df['x'].unique()),len(df['y'].unique()))
   df = df.replace(nan,"NaN")
   df = df.sort(['y','x'])
   if os.path.isfile(outfile):
           os.remove(outfile)
   f = open(outfile,'w')
   f.write(header1+"\n")
   f.write(header2+"\n")
   f.write(header3+"\n")
   df.to_csv(f,sep=" ",header=False,index=False)
   f.close()

def plot_all_deltas(out_file="All_Deltas.png"):
    import matplotlib.pyplot as plt 
    import seaborn as sns
    from numpy import argmin,abs
    import os
    sns.__version__

    bl_info_DF = make_csv(out_file='')
    bl_info_DF = bl_info_DF.sort("U_BL")

    cmap_SS = sns.color_palette("Reds_r",len(bl_info_DF.U_BL.unique()))
    cmap_PS = sns.color_palette("Blues_r",len(bl_info_DF.U_BL.unique()))

    # Marker defines if open or closed
    marker = {
        'open'   : 'o',
        'closed' : 'x',
    }
    fig,ax = plt.subplots( 1, 3, figsize = ( 30, 10 ) )
    for row_index, row in bl_info_DF.iterrows():
        if row.Side == 'SS' or row.Side == 'NA':
            cmap = cmap_SS
        else:
            cmap = cmap_PS

        if row.Side=='NA':
            label_side = ''
        else: label_side = row.Side
        label = "$U_\\infty = {0:.1f}$ m/s {1}".\
        format(row.U,label_side)

        ax[0].scatter(
            row.alpha,
            row.Delta_BL,
            color = cmap[
                argmin(abs(bl_info_DF.U_BL.unique()-row.U_BL))
            ],
            marker = marker[row.Test_section],
            s=75,
            label = label
        )
        ax[0].annotate(label,
                    xy=(row.alpha,row.Delta_BL),
                    xytext=(row.alpha+1,row.Delta_BL+0.5),
                    arrowprops=dict(
                        facecolor='black', 
                        arrowstyle='->',
                        lw=2
                    ),
                   )

        ax[1].scatter(
            row.alpha,
            row.Delta_displacement,
            color = cmap[
                argmin(abs(bl_info_DF.U_BL.unique()-row.U_BL))
            ],
            marker = marker[row.Test_section],
            s=75,
            label = label
        )
        ax[1].annotate(label,
                    xy=(row.alpha,row.Delta_displacement),
                    xytext=(row.alpha+1,row.Delta_displacement+0.5),
                    arrowprops=dict(
                        facecolor='black', 
                        arrowstyle='->',
                        lw=2
                    ),
                   )
        ax[2].scatter(
            row.alpha,
            row.Delta_momentum,
            color = cmap[
                argmin(abs(bl_info_DF.U_BL.unique()-row.U_BL))
            ],
            marker = marker[row.Test_section],
            s=75,
            label = label
        )
        ax[2].annotate(label,
                    xy=(row.alpha,row.Delta_momentum),
                    xytext=(row.alpha+1,row.Delta_momentum+0.5),
                    arrowprops=dict(
                        facecolor='black', 
                        arrowstyle='->',
                        lw=2
                    ),
                   )
    for a in ax:
        a.set_xlabel("$\\alpha_g$ [deg]")
        a.set_xticks([0,6,12])
    ax[0].set_xlim(-2,16)
    ax[0].set_ylabel("$\\delta_{{{0}}}$ [mm]".format(BL_pct))
    ax[1].set_ylabel("$\\delta^*$ [mm]".format(BL_pct))
    ax[2].set_ylabel("$\\theta$ [mm]".format(BL_pct))
    #ax.legend(loc='best')
    plt.savefig(os.path.join('images',out_file))

def write_tecplot(case_name,df,outfile,interpolate=False):
   from numpy import nan
   import os

   # Construct the TECPLOT header for the combined data file
   header1 = "TITLE = \"test\""
   header2 = "VARIABLES = "
   for var in df.columns:
           header2 = header2+"\"{0}\"".format(var)
           if var != df.columns[-1]:
                   header2 = header2+", "
   header3 = "ZONE T=\"Frame 1\" I={0}, J={1}, F=POINT\n".format(len(df['x'].unique()),len(df['y'].unique()))
   df = df.replace(nan,"NaN")
   df = df.sort(['y','x'])
   if os.path.isfile(outfile):
           os.remove(outfile)
   f = open(outfile,'w')
   f.write(header1+"\n")
   f.write(header2+"\n")
   f.write(header3+"\n")
   df.to_csv(f,sep=" ",header=False,index=False)
   f.close()


def plot_all_BLs():
    for case in data_folders:
        plot_bl(case,'u')

def pickle():
    for case in sorted( raw_data_folders ):
        if case == 'STE_a12_U30_SS':
            pickle_all_data(raw_data_root,case,root)

#pickle()
#make_csv()
#plot_all_deltas()
