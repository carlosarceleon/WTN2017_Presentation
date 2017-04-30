import acoustic_functions as acfunc
import bl_functions as blfunc
import os

acoustic_root = \
        './AcousticDataRoberto/'
acoustic_root = '/home/carlos/Documents/PhD/Articles/Article_3/'\
            + 'Scripts/acoustics_and_bl/AcousticDataRoberto/'
acoustic_campaign = 'MarchData'

acoustic_data_path = os.path.join(acoustic_root)

devices = ['Sr20R21']
speeds = [30,35,40]
alphas = [0,6,12]
#phis   = [0,6]
phis   = [6]

def get_case_name(device,speed,alpha,phi,campaign):
    """ Given a case information, build the case name as it was
    given for the dataset (the folder that contains the case)

    Input:
        device,speed,alpha,phi,campaign
    Output:
        string
    """

    if "March" in campaign:
        return "psd_{0}_a{1:02d}_p{2}_U{3}.mat".\
                format(
                    device,
                    alpha,
                    phi,
                    speed
                )
    if "July" in campaign:
        return "psd_{0}_a{1}_U{2}.mat".\
                format(
                    device,
                    alpha,
                    speed
                )


def get_crossovers():
    """ Gets the crossover of the given case permutations and
    returns a pandas Data Frame with all the information

    Output:
        pandas Data Frame
    """
    import pandas as pd

    crossover_df = pd.DataFrame(
        columns = [
            'campaign',
            'device',
            'alpha',
            'phi',
            'U',
            'crossover'
        ]
    )

    campaign = acoustic_campaign
    for dev in devices:
        for speed in speeds:
            for alpha in alphas:
                for phi in phis:
                    case_name = get_case_name(
                        dev,speed,alpha,phi,campaign
                    )
                    ste_case_name = case_name\
                            .replace(dev,"STE")\
                            .replace("_p{0}".format(phi),'')

                    crossover = acfunc.get_crossover(
                        root_folder = os.path.join(
                            acoustic_root),
                        case        = case_name,
                        relative_to = ste_case_name,
                        #campaign='MarchData'
                    )

                    crossover_dict = {
                        'campaign'  : campaign,
                        'device'    : dev,
                        'alpha'     : alpha,
                        'phi'       : phi,
                        'U'         : speed,
                        'crossover' : crossover
                    }
                    crossover_df = crossover_df\
                            .append(crossover_dict,ignore_index=1)
    return crossover_df

def get_boundary_layers():
    return blfunc.make_csv()

def Strouhal(f,delta,U,side='NA'):
    Strouhal = f*delta/float(U)/1000.
    if side=='SS':
        print "{0},{1},{2},{3}".format(Strouhal,f,delta,float(U))
    return Strouhal

def plot_bl_crossover_relation( article=False, real_Uinf=True ):
    import matplotlib.pyplot as plt 
    import matplotlib.lines as mlines
    import seaborn.apionly as sns
    from matplotlib  import rc
    from os.path     import isfile
    from pandas      import read_pickle, DataFrame
    from scipy.stats import linregress
    from numpy       import array

    import matplotlib as mpl

    sns.set(font_scale=2.0,style='ticks',
            rc={"axes.axisbelow": False,'image.cmap': 'OrRd'})

    rc('text',usetex=True)
    rc('font',family='sans-serif', serif='sans-serif')

    mpl.rcParams['text.latex.preamble'] = [
        r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
        r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
        r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
        r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
    ]


    crossover_df = get_crossovers()
    BL_Data_File = 'outputs/BL_Data_Info.p'

    if isfile( BL_Data_File ):
        bl_df = read_pickle( BL_Data_File )
    else:
        bl_df = get_boundary_layers()
    
    colormap = sns.color_palette(
        'deep',
    )

    marker_dict = {
        '20': '+',
        '30': 'x',
        '35': '1',
        '40': '2',
    }

    color_dict = {
        '20': colormap[0],
        '30': colormap[1],
        '35': colormap[2],
        '40': colormap[3],
    }

    BL_measures = [
        'Delta_BL',
        'Delta_displacement',
        'Delta_momentum'
    ]

    delta_labels = [
        r'\delta_{99}',
        r'\delta^*',
        r'\theta',
    ]

    markersize      = 12
    mew             = 3 # Marker edge width

    marker_legend = [ mlines.Line2D(
        range(1),range(1),
        color           = 'w',
        marker          = marker_dict['30'],
        label           = '30 m/s',
        markersize      = markersize,
        markeredgecolor = color_dict['30'],
        markerfacecolor = 'none',
        mew             = mew,
    ) ]
    marker_legend.append( mlines.Line2D(
        range(1),range(1),
        color           = 'w',
        marker          = marker_dict['35'],
        label           = '35 m/s',
        markersize      = markersize,
        markeredgecolor = color_dict['35'],
        markerfacecolor = 'none',
        mew             = mew,
    ) )
    marker_legend.append( mlines.Line2D(
        range(1),range(1),
        color           = 'w',
        marker          = marker_dict['40'],
        label           = '40 m/s',
        markersize      = markersize,
        markeredgecolor = color_dict['40'],
        markerfacecolor = 'none',
        mew             = mew,
    ) )


    fig_ss,axes_ss = plt.subplots(1,len(BL_measures),sharex=True, 
                                  figsize = ( 10, 4 ) )
    fig_ps,axes_ps = plt.subplots(1,len(BL_measures),sharex=True, 
                                  figsize = ( 10, 4 ) )

    crossover_df = crossover_df[ crossover_df.crossover!=0 ]
    crossover_df = crossover_df[ crossover_df.phi==6 ]

    st_df = DataFrame( columns = [ 'alpha', 'u', 'st', 'measure_param' ])

    for ix,acoustic_case in crossover_df.iterrows():

        delta = bl_df[
            (bl_df.U            == acoustic_case.U)     & \
            (bl_df.alpha        == acoustic_case.alpha) & \
            (bl_df.Test_section == 'open')
        ]

        for BL_measure, ax_ss, ax_ps in zip( BL_measures, axes_ss, axes_ps ):
            
            delta_ps   = delta[ delta.Side == 'PS' ][ BL_measure ]
            delta_ss   = delta[ delta.Side == 'SS' ][ BL_measure ]
            delta_zero = delta[ delta.Side == 'NA' ][ BL_measure ]

            color  = color_dict[  '{0}'.format( int( delta.U.values[0] ) ) ]
            marker = marker_dict[ '{0}'.format( int( delta.U.values[0] ) ) ]

            marker_size = 175

            plot_args = {
                'edgecolors' : color,
                's'          : marker_size,
                'facecolors' : color,
                'linewidth'  : 3,
                'marker'     : marker,
            }

            if not delta_zero.empty:
                U_infty = delta[delta.Side=='NA'].U_BL.values[0]
                U_th    = delta[delta.Side=='NA'].U.values[0]
                for ax in [ ax_ss, ax_ps ]:
                    st = Strouhal(
                            acoustic_case.crossover,
                            delta_zero.values[0],
                            U_infty
                        )
                    ax.scatter(
                        acoustic_case.alpha,
                        st,
                        **plot_args
                    )

                    st_df = st_df.append(
                        DataFrame(
                            data = {
                                'alpha':         0,
                                'u':             U_th,
                                'st':            st,
                                'measure_param': BL_measure
                            }, index = [0]
                        ), ignore_index = True
                    )

            else:
                ax       = ax_ss
                U_infty  = delta[ delta.Side == 'SS' ].U_BL.values[0]
                U_th     = delta[ delta.Side == 'SS' ].U.values[0]
                delta_th = delta_ss.values[0]

                st = Strouhal(
                    acoustic_case.crossover,
                    delta_th,
                    U_infty,
                )
                ax.scatter(
                    acoustic_case.alpha,
                    st,
                    **plot_args
                )

                st_df = st_df.append(
                    DataFrame(
                        data = {
                            'alpha':         acoustic_case.alpha,
                            'u':             U_th,
                            'st':            st,
                            'measure_param': BL_measure
                        }, index = [0]
                    ), ignore_index = True
                )

                ax       = ax_ps
                U_infty  = delta[ delta.Side == 'PS' ].U_BL.values[0]
                U_th     = delta[ delta.Side == 'PS' ].U.values[0]
                delta_th = delta_ps.values[0]

                st = Strouhal(
                    acoustic_case.crossover,
                    delta_th,
                    U_infty,
                )
                ax.scatter(
                    acoustic_case.alpha,
                    st,
                    **plot_args
                )

                st_df = st_df.append(
                    DataFrame(
                        data = {
                            'alpha':         -acoustic_case.alpha,
                            'u':             U_th,
                            'st':            st,
                            'measure_param': BL_measure
                        }, index = [0]
                    ), ignore_index = True
                )

    for axes, side_ix in zip( [ axes_ss, axes_ps ], [ 1, -1 ] ):
        for ax, delta_label, bl_meas in zip( axes, delta_labels, BL_measures ):
            ax.set_xticks([0,6,12])
            ax.set_xlim(-2, 14)
            ax.set_xticklabels([r'$0^\circ$',r'$6^\circ$',r'$12^\circ$'])
            ax.set_xlabel("$\\alpha$")
            ax.set_ylabel(
                r"$\mathrm{{St}}_\mathrm{{c}}=f_\mathrm{{c}}{0}/u_e$"\
                .format( delta_label ),
            )
            st_vals = st_df[ 
                ( side_ix * st_df.alpha >= 0 ) & \
                ( st_df.measure_param == bl_meas ) 
            ].st.values
            alpha_vals = side_ix * st_df[ 
                ( side_ix * st_df.alpha >= 0 ) & \
                ( st_df.measure_param == bl_meas ) 
            ].alpha.values

            slope, intercept, rvalue, pvalue, stderr = linregress( 
                alpha_vals, 
                st_vals 
            )
            
            ax.plot( 
                [ 0, 6, 12 ], 
                slope * array([ 0, 6, 12 ]) + intercept,
                '--',
                color = 'darkgray'
            )

            ax.text(
                0.95, 0.95,
                r'$ r^2 = {0:.2f}$'.format( rvalue**2 ),
                ha = 'right', va = 'top',
                fontsize = 20,
                transform=ax.transAxes,
                bbox=dict( facecolor='white', edgecolor='none', alpha = 0.6 )
            )

    #ax[1].legend(fontsize=18,loc='upper right')
    axes_ss[0].legend(handles = marker_legend,
               bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=4, borderaxespad=0., numpoints = 1
              )
    axes_ps[0].legend(handles = marker_legend,
               bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=4, borderaxespad=0., numpoints = 1
              )

    fig_ss.tight_layout()
    fig_ps.tight_layout()
    fig_ss.savefig(
        './article_images/Crossover_Strouhal_and_AoA_SS.png',
        bbox_inches='tight')
    fig_ps.savefig(
        './article_images/Crossover_Strouhal_and_AoA_PS.png',
        bbox_inches='tight')

    st_df = st_df.drop_duplicates()
    st_df.to_pickle( 'st_results.p' )

    plt.close()

def predict_f_at_twenty():
    from pandas import read_pickle
    import matplotlib.pyplot as plt 
    from os.path import isfile
    from matplotlib  import rc
    import seaborn as sns
    from scipy.stats import linregress
    from numpy import array

    st_df = read_pickle( 'st_results.p' )
    sns.set_context('paper')
    sns.set_style("whitegrid")
    sns.set(font='serif',font_scale=2.0,style='ticks')
    rc('text',usetex=True)
    rc('font',family='serif', serif='cm10')
    color  = 'darkgray'
    marker = 'x'

    colormap = sns.color_palette(
        'deep',
    )

    marker_dict = {
        '20': '+',
        '30': 'x',
        '35': '1',
        '40': '2',
    }

    color_dict = {
        '20': colormap[0],
        '30': colormap[1],
        '35': colormap[2],
        '40': colormap[3],
    }

    marker_size = 175

    BL_Data_File = 'outputs/BL_Data_Info.p'

    if isfile( BL_Data_File ):
        bl_df = read_pickle( BL_Data_File )
    else:
        return 0

    bl_df = bl_df[ 
        ( bl_df.U == 20 ) & \
        ( bl_df.alpha == 12 ) & \
        ( bl_df.Side == 'PS' ) 
    ]

    # Get boundary layer thickness parameters
    deltas = st_df.measure_param.unique()

    fig_ss,axes_ss = plt.subplots(1,1,sharex=True, 
                                  figsize = ( 4, 3 ) )


    deltas_dict = {
        'Delta_BL':           [0, r'$\textrm{St}_{\textrm{c},\delta_{99}}$'],
        'Delta_displacement': [1, r'$\textrm{St}_{\textrm{c},\delta^*}$'],
        'Delta_momentum':     [2, r'$\textrm{St}_{\textrm{c},\theta}$'],
    }

    f_twenty = []

    for d,cnt in zip( deltas, range(len(deltas)) ):
        print " # {0} ".format( d )

        st_df_d = st_df[ 
            ( st_df.measure_param == d ) & \
            ( st_df.alpha == -12 ) 
        ]

        f_c_vec = []
        u_vec = []
        for st, u in zip( st_df_d.st.values, st_df_d.u.values ):
            
            color = color_dict[ "{0:d}".format(int(u))]
            marker = marker_dict["{0:d}".format(int(u))]
            
            plot_args = {
                'edgecolors' : color,
                's'          : marker_size,
                'facecolors' : color,
                'linewidth'  : 3,
                'marker'     : marker,
            }

            f_c = st * bl_df.U_BL.values 
            f_c = f_c / ( bl_df[ d ].values / 1000. )

            axes_ss.scatter( 
                deltas_dict[d][0], 
                [f_c],
                **plot_args
            )

            f_c_vec.append( f_c[0] )
            u_vec.append( u )

            if d == deltas[-1]:
                axes_ss.annotate(
                    r'$\mathrm{{St}}_\textrm{{c}}'\
                    r'\left( u_e \approx {0}\, \textrm{{m/s}}\right)$'.\
                    format( int(u) ),
                    xy         = (2.1, f_c[-1]),
                    xycoords   = 'data',
                    xytext     = (2.7, f_c[-1]),
                    textcoords = 'data',
                    arrowprops = dict(
                        arrowstyle="->",
                        connectionstyle=\
                        "arc,angleA=-20,armA=00,angleB=30,armB=30,rad=15",
                        linewidth = 1.5
                    ),
                    fontsize = 18
                )

        slope, intercept, rvalue, pvalue, stderr = linregress( 
            u_vec, 
            f_c_vec 
        )

        color = color_dict[ '20' ]
        marker = marker_dict[ '20' ]
        
        plot_args = {
            'edgecolors' : color,
            's'          : marker_size,
            'facecolors' : color,
            'linewidth'  : 3,
            'marker'     : marker,
        }

        f_twenty.append( [slope * 20 + intercept] )
        axes_ss.scatter(
            deltas_dict[d][0], 
            f_twenty[-1],
            **plot_args
        )
        print f_twenty[-1][0]
        print bl_df.U_BL.values[0]
        print bl_df.Delta_displacement.values[0]/1000.
        print f_twenty[-1][0] * (bl_df.Delta_displacement.values[0]/1000.) / bl_df.U_BL.values[0]

    axes_ss.annotate(
        r'$\mathrm{{St}}_\textrm{{c}}\left( u_e \approx {0}\, '.format( 20 )\
        +r'\textrm{{m/s}} \right)$',
        xy         = (2.1, f_twenty[-1][0]),
        xycoords   = 'data',
        xytext     = (2.7, f_twenty[-1][0]),
        textcoords = 'data',
        arrowprops = dict(
            arrowstyle="->",
            connectionstyle="arc,angleA=-20,armA=00,angleB=30,armB=30,rad=15",
            linewidth = 1.5
        ),
        fontsize = 18
    )

    axes_ss.set_xticks( [0, 1, 2] )
    axes_ss.set_xticklabels( [ deltas_dict[f][1] for f in deltas ] )
    axes_ss.set_title( 
        r'$f_{\textrm{c}} = \textrm{St}_\textrm{c}u_e /\delta$',
        fontsize = 20
    )
    axes_ss.set_ylabel( 
        r'$f_{\textrm{c}}\left( u_e \approx 20\,\textrm{m/s}\right) $ [Hz]' 
    )

    axes_ss.set_yticks( array( [ 500, 750, 1000, array( f_twenty ).mean() ] ) )
    axes_ss.axhline( array( f_twenty ).mean(), xmax =  0.9, 
                    ls = '--', color = color
                   )
    axes_ss.text(
        -0.4, array( f_twenty ).mean(),
        r'mean',
        ha = 'left', va = 'bottom',
        fontsize = 20,
        color= 'k'
    )

    plt.savefig( 
        'article_images/fc_prediction_twenty.png', 
        bbox_inches = 'tight' 
    )

