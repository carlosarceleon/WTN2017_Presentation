output_root = 'Results'
def return_bl_parameters( case, x_series ):
    from pandas import read_pickle, DataFrame
    from scipy.optimize import curve_fit

    def quad_func( x, a , b, c):
        return a*x**2 + b*x + c

    bl_pickle_file = '/home/carlos/Documents/PhD/Articles/Article_3/' + \
            'Scripts/time_resolved/BLData.p'
    bl_df = read_pickle( bl_pickle_file )

    if not case in bl_df.case.unique():
        print bl_df.case.unique()
        print case

    if 'z00' in case and not "STE" in case:
        limits = (0, 43)
    elif 'z05' in case:
        limits = (0, 25)
    else:
        limits = (-8, 10)

    variables = [ 'Ue', 'delta_99', 'delta_displacement', 
                'delta_momentum' ]

    bl_param_df = DataFrame()

    for x in x_series:
        data = {}
        for var in variables:

            case_df = bl_df[ bl_df.case == case ][ ['x', var] ].dropna()

            try:
                popt, pvar = curve_fit( 
                    quad_func, 
                    case_df[ ( case_df.x <= limits[1] ) & \
                            ( case_df.x >= limits[0] )].x, 
                    case_df[ ( case_df.x <= limits[1] ) & \
                            ( case_df.x >= limits[0] )][var] 
                )
            except TypeError:
                print x, case
                raise

            data[var] = quad_func( x, popt[0], popt[1], popt[2] ) ,

        bl_param_df = bl_param_df.append(
            DataFrame( data = data, index = [0] ),
            ignore_index = False
        )

    return bl_param_df

def remove_angle_jumps(df):
    from math import pi
    from copy import copy

    if df.empty:
        return df

    def remove_large_diff_jumps( df, var ):

        df[ 'diff_phi' ] = 0
        
        diff_phi = df.ix[ ixs[0] : ixs[-2] ][ var ].values - \
                df.ix[ ixs[1] : ixs[-1] ][ var ].values

        df['diff_phi'].ix[ : ixs[-2]] = copy( diff_phi )

        ix = df[ df.diff_phi > pi ].index

        if len( ix ):
            df.ix[ ix[0] + 1 : ].phi = df.ix[ ix[0] + 1 : ].phi + 2 * pi
        
        df = df.drop( 'diff_phi' , axis = 1)

        return df

    def remove_negative_values( df, var ):
        df.loc[ 
            ( df[var] < 0) & ( df.f < 3000 ), var 
        ] = df[ ( df[var] < 0) & ( df.f < 3000 ) ][ var ] + 2 * pi

        return df



    variables_with_phi = [v for v in df.columns if 'phi' in v]

    ixs = copy(df.index)

    for var in variables_with_phi:

        df = remove_negative_values( df, var )
        df = remove_large_diff_jumps( df, var )
        df = remove_large_diff_jumps( df, var )
        df = remove_large_diff_jumps( df, var )


    return df

def plot_frequency_spectra( hdf_cases , root = '', var = 'v'):
    import matplotlib.pyplot as plt
    from scipy.signal                 import welch
    from pandas                       import read_hdf
    from raw_data_processing_routines import find_nearest
    from os.path                      import split, join
    from numpy                        import log10, arange, array, nan
    from numpy                        import append
    from matplotlib                   import rc
    from math                         import degrees, atan
    from matplotlib.cbook             import get_sample_data

    #schematic = '/home/carlos/Documents/PhD/Articles/Article_2/'+\
            #'Figures/measurement_locations_TE_m2_wo05.png'

    fontsize = 20
    import matplotlib as mpl

    rc('text',usetex=True)
    rc('font',family='sans-serif', serif='sans-serif')

    mpl.rcParams['text.latex.preamble'] = [
        r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
        r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
        r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
        r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
    ]



    if not root:
        root = '/home/carlos/Documents/PhD/Articles/Article_3/' + \
                'Scripts/time_resolved/Results_v2'

    y_locations_f = array( [ 1.5, 5, 8, 12 ] )

    fig_f, axes_f = plt.subplots( 2 , 2 , 
                             figsize = ( 10, 8 ), sharex = True )

    for hdf_name in [ join( root, f ) for f in hdf_cases ]:

        if 'z00' in hdf_name and not 'STE' in hdf_name:
            x_bl_loc = 38
        elif 'z05' in hdf_name:
            x_bl_loc = 18
        elif 'z10' in hdf_name or 'STE' in hdf_name:
            x_bl_loc = -1

        case = split( hdf_name )[1].replace('.hdf5','').replace( '_mfr', '' )

        case_coords = read_hdf( hdf_name , where = ( 't=0' ),
                              columns = [ 'x', 'y' , 'near_x_2h', 
                                         'near_y_delta', 'u', 'v'],
                              )

        print "   Performing the wavenumber spectra analysis of {0}"\
                .format(case)

        delta = (case_coords.y / case_coords.near_y_delta).mean()

        case_coords = case_coords.sort_values( 
            by = [ 'y', 'x' ] 
        )

        case_coords = sort_block_coordinates_by_height_and_streamwise( 
            case_coords , case = case, 
            plot = "CoordExtraction_FreqPlots_"+case+".png", truncate = False
        )

        available_x_ix = case_coords[ 
            case_coords.x == find_nearest( x_bl_loc, case_coords.x.values)
        ].stream_ix.values[0]

        color, marker, cmap = get_color_and_marker( 
            case
        )

        plot_config = {
            'marker'          : marker,
            'markeredgewidth' : markeredgewidth,
            'markerfacecolor' : markerfacecolor,
            'markeredgecolor' : color,
            'markersize'      : markersize,
            'mew'             : mew,
            'color'           : color,
            'lw'              : 4
        }

        for y_req, ix, ax_f in zip( 
            y_locations_f[::-1], range( len( y_locations_f ) ) , axes_f.ravel()
        ):

            nearest_y_bl = find_nearest( y_req / delta, case_coords[ 
                case_coords.stream_ix == available_x_ix
            ].near_y_delta.values )

            ref_coord = case_coords[ 
                ( case_coords.near_y_delta         == nearest_y_bl      ) & \
                ( case_coords.stream_ix == available_x_ix )
                ]

            print y_req / delta

            case_ts = read_hdf(
                hdf_name, 
                where = ( 
                    'near_x_2h = {0} & near_y_delta = {1}'.format( 
                        ref_coord.near_x_2h.values[0],
                        ref_coord.near_y_delta.values[0] 
                    )
                ), 
                columns = [
                    'near_x_2h',
                    'near_y_delta',
                    't',
                    'u',
                    'v',
                    'x',
                    'y'
                ]
            ).sort_values( by = 't' ).reset_index( drop = True )

            if case_ts.empty:
                print "  Warning: found the following case empty: {0}".\
                        format( case )
                print "     at the location y = {0}".format(
                    ref_coord.near_y_delta.values[0]
                )

            case_ts[ 'case' ] = case

            case_ts.loc[ 
                is_outlier( case_ts, 'u', thresh = outlier_thresh) , 'u'
            ] = nan 
            case_ts.loc[ 
                is_outlier( case_ts, 'v', thresh = outlier_thresh) , 'v'
            ] = nan 

            case_ts[ 'u' ] = case_ts['u'].interpolate( 
                method='spline', order=3, s=0.
            )
            case_ts[ 'v' ] = case_ts['v'].interpolate( 
                method='spline', order=3, s=0.
            )

            print "        resulting mean flow angle at ({0:.1f},{1:.1f}) is:".\
                    format( 
                        case_ts.x.unique()[0], 
                        case_ts.y.unique()[0]
                    ),;
            print " {0:.2f} degrees"\
                    .format(
                        degrees( atan( case_ts.v.mean() / case_ts.u.mean() ) )
                    )

            freq, Pxx_u = welch(
                x       = case_ts[ 'u' ].dropna(),
                nperseg = nperseg,
                fs      = fs,
                scaling = 'spectrum',
            )
            freq, Pxx_v = welch(
                x       = case_ts[ 'v' ].dropna(),
                nperseg = nperseg,
                fs      = fs,
                scaling = 'spectrum',
            )
            

            #ax_f.plot(
            #    freq[:-1] / 1000.,
            #    10 * log10( Pxx_u )[:-1],
            #    **plot_config
            #)
            ax_f.plot(
                freq[:-1] / 1000.,
                10 * log10( Pxx_v )[:-1],
                **plot_config
            )

            kolmogorov_law_curve = get_kolmogorov_law_curve( 
                x_lim = (1, 3) 
            )
            kolmogorov_law_curve[1] = \
                    kolmogorov_law_curve[1] + 3

            if ix == 0 or ix == 1:
                y_shift = -20 
            elif ix == 2 or ix == 3:
                y_shift = -15 
            ax_f.plot( 
                kolmogorov_law_curve[0] , 
                kolmogorov_law_curve[1] + y_shift, 
                '-',
                color = 'k' ,
                lw    = 3,
            )

            ax_f.set_xlim( 0.2, 5   )
            ax_f.set_xscale('log')
            ax_f.tick_params(axis='both', which='major', 
                                      labelsize=fontsize)
            #ax_f.set_yticks( arange( -20, 5, 5 ) )
            #ax_f.set_yticks( arange( -20, 0, 5 ) )
            #axes[ix].set_yticks( arange( -20, 0, 5 ) )
            ax_f.set_xticks( append( [0.2], arange( 1, 6, 1) ) )
            ax_f.set_xticks( append( [0.2], arange( 1, 6, 1) ) )
            ax_f.set_xticklabels( [0.2,1,2,3,4,5] )
            ax_f.set_xticklabels( [0.2,1,2,3,4,5] )
            t = ax_f.text( 
                0.05, 0.05, r'$y = {0:.1f}$ mm'.format(y_req),
                ha = 'left', fontsize = fontsize,
                transform = ax_f.transAxes
            )
            t.set_bbox(
                dict(color='white', alpha=0.5, edgecolor='white', zorder = 300)
            )

    axes_f[0][0].annotate(
        r'$ 10\log_{10}\,f^{-5/3}$',
        xy                  = ( 1.300, -20 ),
        xytext              = ( 0.25,  -30 ),
        horizontalalignment = 'left',
        verticalalignment   = 'center',
        fontsize            = fontsize,
        arrowprops          = dict(
            arrowstyle      = "->",
            connectionstyle = "arc3, rad=-0.3"
        )
    )

    axes_f[-1][0].set_xlabel( r'$f$ [kHz]',
                       fontsize = fontsize )
    axes_f[-1][1].set_xlabel( r'$f$ [kHz]',
                       fontsize = fontsize )

    axes_f[0][0].set_ylabel( 
        r"$10\log_{10}\Phi_{vv}$ [dB]",
                      fontsize = fontsize )
    axes_f[0][0].set_ylim( -45, -5 )
    axes_f[0][1].set_ylim( -45, -5 )
    axes_f[1][0].set_ylabel( 
        r"$10\log_{10}\Phi_{vv}$ [dB]",
                      fontsize = fontsize )
    axes_f[1][0].set_ylim( -25, -5 )
    axes_f[1][1].set_ylim( -25, -5 )

    schematic = '/home/carlos/Documents/PhD/Articles/Article_3/Scripts/'\
            +'time_resolved/LegendNo00.png'

    im = plt.imread( get_sample_data( schematic  ) )

    newax = fig_f.add_axes([0.585, 0.18, 0.16, 0.16], 
                             anchor = 'SE', zorder=100)

    newax.imshow(im)
    newax.axis('off')

    fig_f.subplots_adjust(wspace = 0.3)

    fig_f.savefig( join( output_root, 'f_plot.png' ), bbox_inches = 'tight' )

def plot_wavenumber_spectra( hdf_cases , root = '', var = 'v'):
    import matplotlib.pyplot as plt
    from scipy.signal                    import welch
    from pandas                          import read_hdf, read_pickle
    from raw_data_processing_routines    import find_nearest
    from os.path                         import split, join
    from numpy                           import log10, arange, array, nan
    from numpy                           import append
    from matplotlib                      import rc
    from math                            import degrees, atan
    from matplotlib.cbook                import get_sample_data

    #schematic = '/home/carlos/Documents/PhD/Articles/Article_2/'+\
            #'Figures/measurement_locations_TE_m2_wo05.png'

    #rc('font',family='serif', serif='Linux Libertine', size=40)
    #rc('text',usetex=True)
    fontsize = 30
    import matplotlib as mpl

    rc('text',usetex=True)
    rc('font',family='sans-serif', serif='sans-serif')

    mpl.rcParams['text.latex.preamble'] = [
        r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
        r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
        r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
        r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
    ]



    if not root:
        root = '/home/carlos/Documents/PhD/Articles/Article_3/' + \
                'Scripts/time_resolved/Results_v2'

    y_locations_f = array( [ 1.5, 5, 8, 12 ] )
    y_locations_k = array( [ 1.5, 5, 8, 12 ] )

    fig_f, axes_f = plt.subplots( len( y_locations_f ) , 2 , 
                             figsize = ( 14, 11 ), sharex = True, sharey = True)
    fig_k, axes_k = plt.subplots( 2 , 2 , 
                             figsize = ( 10, 8 ), sharex = True, sharey = True)
    fig_s, axes_s = plt.subplots( 1 , 2 , 
                             figsize = ( 10, 5 ), sharex = False, sharey = True)

    for hdf_name in [ join( root, f ) for f in hdf_cases ]:

        if 'z00' in hdf_name and not 'STE' in hdf_name:
            x_bl_loc = 38
        elif 'z05' in hdf_name:
            x_bl_loc = 18
        elif 'z10' in hdf_name or 'STE' in hdf_name:
            x_bl_loc = -1

        case = split( hdf_name )[1].replace('.hdf5','').replace( '_mfr', '' )

        case_coords = read_hdf( hdf_name , where = ( 't=0' ),
                              columns = [ 'x', 'y' , 'near_x_2h', 
                                         'near_y_delta', 'u', 'v'],
                              )

        Uc_df = read_pickle( join( root, 'Uc_Spline_{0}.p'.format( case ) ) )

        print "   Performing the wavenumber spectra analysis of {0}"\
                .format(case)

        delta = (case_coords.y / case_coords.near_y_delta).mean()

        case_coords = case_coords.sort_values( 
            by = [ 'y', 'x' ] 
        )

        case_coords = sort_block_coordinates_by_height_and_streamwise( 
            case_coords , case = case, 
            plot = "CoordExtraction_FreqPlots_"+case+".png", truncate = False
        )

        available_x_ix = case_coords[ 
            case_coords.x == find_nearest( x_bl_loc, case_coords.x.values)
        ].stream_ix.values[0]

        color, marker, cmap = get_color_and_marker( 
            case
        )

        plot_config = {
            'marker'          : marker,
            'markeredgewidth' : markeredgewidth,
            'markerfacecolor' : markerfacecolor,
            'markeredgecolor' : color,
            'markersize'      : markersize,
            'mew'             : mew,
            'color'           : color,
            'lw'              : 4
        }

        f_cnt = 0
        for y_req, ix, ax_k in zip( 
            y_locations_k[::-1], range( len( y_locations_k ) ) , axes_k.ravel()
        ):

            nearest_y_bl = find_nearest( y_req / delta, case_coords[ 
                case_coords.stream_ix == available_x_ix
            ].near_y_delta.values )

            Uc = Uc_df[ 
                Uc_df.y_bl == find_nearest( y_req, Uc_df.y_bl.values ) 
            ].Uc_spline.values[0]

            ref_coord = case_coords[ 
                ( case_coords.near_y_delta == nearest_y_bl ) & \
                ( case_coords.stream_ix == available_x_ix )
                ]

            case_ts = read_hdf(
                hdf_name, 
                where = ( 
                    'near_x_2h = {0} & near_y_delta = {1}'.format( 
                        ref_coord.near_x_2h.values[0],
                        ref_coord.near_y_delta.values[0],
                    )
                ), 
                columns = [
                    'near_x_2h',
                    'near_y_delta',
                    't',
                    'u',
                    'v',
                    'x',
                    'y'
                ]
            ).sort_values( by = 't' ).reset_index( drop = True )

            if case_ts.empty:
                print "  Warning: found the following case empty: {0}".\
                        format( case )
                print "     at the location y = {0}".format(
                    ref_coord.near_y_delta.values[0]
                )

            case_ts[ 'case' ] = case

            case_ts.loc[ 
                is_outlier( case_ts, 'u', thresh = outlier_thresh) , 'u'
            ] = nan 
            case_ts.loc[ 
                is_outlier( case_ts, 'v', thresh = outlier_thresh) , 'v'
            ] = nan 

            case_ts[ 'u' ] = case_ts['u'].interpolate( 
                method='spline', order=3, s=0.
            )
            case_ts[ 'v' ] = case_ts['v'].interpolate( 
                method='spline', order=3, s=0.
            )

            print "        resulting mean flow angle at ({0:.1f},{1:.1f}) is:".\
                    format( 
                        case_ts.x.unique()[0], 
                        case_ts.y.unique()[0]
                    ),;
            print " {0:.2f} degrees"\
                    .format(
                        degrees( atan( case_ts.v.mean() / case_ts.u.mean() ) )
                    )

            axes_s[0].plot(
                [ case_ts.u.mean() ],
                [ case_ts.y.unique()[0] ],
                **plot_config
            )
            axes_s[1].plot(
                [ case_ts.v.mean() ],
                [ case_ts.y.unique()[0] ],
                **plot_config
            )

            freq, Pxx_u = welch(
                x       = case_ts[ 'u' ].dropna(),
                nperseg = nperseg,
                fs      = fs,
                scaling = 'spectrum',
            )
            freq, Pxx_v = welch(
                x       = case_ts[ 'v' ].dropna(),
                nperseg = nperseg,
                fs      = fs,
                scaling = 'spectrum',
            )
            

            k = freq / Uc

            ax_k.plot(
                k[:-1] / 100.,
                k[:-1] * Pxx_u[:-1],#10 * log10( Pxx_v )[:-1],
                **plot_config
            )

            ax_k.set_xlim( 0.1,   3.50 )
            ax_k.set_xlim( 0.1,   3.50 )
            #ax_k.set_xscale('log')
            ax_k.tick_params(axis='both', which='major', 
                                      labelsize=fontsize)
            ax_k.tick_params(axis='both', which='major', 
                                      labelsize=fontsize)


            ax_k.text( 0.11, 10, r'$y = {0:.1f}$ mm'.format(y_req),
                             ha = 'left', fontsize = fontsize
                            )

            ax_k.set_xticks( [ 0.1, 1, 2, 3 ] )
            ax_k.set_xticklabels( [ 0.1, 1, 2, 3 ] )

            if y_req in y_locations_f:
                
                axes_f[f_cnt][0].plot(
                    freq[:-1] / 1000.,
                    10 * log10( Pxx_u )[:-1],
                    **plot_config
                )
                axes_f[f_cnt][1].plot(
                    freq[:-1] / 1000.,
                    10 * log10( Pxx_v )[:-1],
                    **plot_config
                )

                kolmogorov_law_curve = get_kolmogorov_law_curve( 
                    x_lim = (1, 3) 
                )
                kolmogorov_law_curve[1] = \
                        kolmogorov_law_curve[1] + 3

                axes_f[f_cnt][0].plot( 
                    kolmogorov_law_curve[0] , 
                    kolmogorov_law_curve[1] - 7, 
                    '-',
                    color = 'k' ,
                    lw    = 4,
                )
                axes_f[f_cnt][1].plot( 
                    kolmogorov_law_curve[0] , 
                    kolmogorov_law_curve[1] - 10, 
                    '-',
                    color = 'k' ,
                    lw    = 4,
                )

                axes_f[f_cnt][0].set_xlim( 0.2, 5   )
                axes_f[f_cnt][1].set_xlim( 0.2, 5   )
                axes_f[f_cnt][0].set_xscale('log')
                axes_f[f_cnt][1].set_xscale('log')
                axes_f[f_cnt][0].tick_params(axis='both', which='major', 
                                          labelsize=fontsize)
                axes_f[f_cnt][1].tick_params(axis='both', which='major', 
                                          labelsize=fontsize)
                axes_f[f_cnt][0].set_yticks( arange( -20, 5, 5 ) )
                axes_f[f_cnt][1].set_yticks( arange( -20, 0, 5 ) )
                #axes[f_cnt][2].set_yticks( arange( -20, 0, 5 ) )
                axes_f[f_cnt][0].set_xticks( append( [0.2], arange( 1, 6, 1) ) )
                axes_f[f_cnt][1].set_xticks( append( [0.2], arange( 1, 6, 1) ) )
                axes_f[f_cnt][0].set_xticklabels( [0.2,1,2,3,4,5] )
                axes_f[f_cnt][1].set_xticklabels( [0.2,1,2,3,4,5] )
                axes_f[f_cnt][1].text( 
                    1.0, -5, r'$y = {0:.2f} mm$'.format(y_req),
                    ha = 'left', fontsize = fontsize
                )
                axes_f[f_cnt][0].text( 1.800, -6, r'$f^{-5/3}$',
                                 ha = 'left', fontsize = fontsize
                                )
                axes_f[f_cnt][1].text( 1.800, -10, r'$f^{-5/3}$',
                                 ha = 'left', fontsize = fontsize
                                )
                f_cnt += 1


    #axes_k[0][0].set_ylim( 0, 11 )
    axes_k[0][0].set_ylabel( r'$k_u \Phi_{vv}(k_u)$' , fontsize = fontsize )
    axes_k[1][0].set_ylabel( r'$k_u \Phi_{vv}(k_u)$' , fontsize = fontsize )

    axes_f[-1][0].set_xlabel( r'$f$ [kHz]',
                       fontsize = fontsize )
    axes_f[-1][1].set_xlabel( r'$f$ [kHz]',
                       fontsize = fontsize )
    axes_k[-1][0].set_xlabel( r'$k_u$ [1/m] $\times 10^2$',
                       fontsize = fontsize )
    axes_k[-1][1].set_xlabel( r'$k_u$ [1/m] $\times 10^2$',
                       fontsize = fontsize )
    axes_f[0][0].set_title( 
        r"$10\log_{10}\left[\Phi_{uu} \left( f\right)\right]$ [dB]",
                      fontsize = fontsize )
    axes_f[0][1].set_title( 
        r"$10\log_{10}\left[\Phi_{vv} \left( f\right)\right]$ [dB]",
                      fontsize = fontsize )

    schematic = '/home/carlos/Documents/PhD/Articles/Article_3/Scripts/'\
            +'time_resolved/LegendNo00.png'

    im = plt.imread( get_sample_data( schematic  ) )

    newax = fig_k.add_axes([0.15, 0.75, 0.16, 0.16], 
                             anchor = 'SE', zorder=100)

    newax.imshow(im)
    newax.axis('off')

    fig_k.subplots_adjust( wspace = -1.0 )

    #fig_f.tight_layout()
    #fig_f.savefig( join( output_root, 'f_plot.png' ), bbox_inches = 'tight' )
    fig_k.tight_layout()
    fig_k.savefig( join( output_root, 'k_plot.png' ), bbox_inches = 'tight' )
    fig_s.tight_layout()
    fig_s.savefig( join( output_root, 's_plot.png' ), bbox_inches = 'tight' )



def plot_pickled_Uc( pickled_Uc_list , root = '', print_integration = False):
    import matplotlib.pyplot as plt
    from pandas                          import read_pickle, DataFrame, concat
    from pandas                          import Series
    #from matplotlib.cbook import get_sample_data
    from scipy.optimize                  import curve_fit
    from matplotlib                      import rc
    from scipy.interpolate               import splprep, splev
    from scipy.integrate                 import simps
    from numpy                           import linspace, arange
    from numpy                           import diff as npdiff
    from os.path                         import join, isfile
    from raw_data_processing_routines    import find_nearest
    from re                              import findall

    #schematic = '/home/carlos/Documents/PhD/Articles/Article_2/'+\
            #'Figures/measurement_locations_TE_m2.png'

    if isfile( 'Um.p' ):
        U_mean_df = read_pickle( 'Um.p' )
    else:
        U_mean_df = DataFrame()

    import matplotlib as mpl

    rc('text',usetex=True)
    rc('font',family='sans-serif', serif='sans-serif')

    mpl.rcParams['text.latex.preamble'] = [
        r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
        r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
        r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
        r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
    ]


    fontsize = 20

    if not root:
        root = '/home/carlos/Documents/PhD/Articles/Article_3/' + \
                'Scripts/time_resolved/Results_v2'

    def log_law( y, y0, alpha ):
        from numpy import log

        return alpha * log( y / y0 )


    fig, axes = plt.subplots(1, len( pickled_Uc_list ), 
                           sharex = True, sharey = True , figsize = (20, 8) )

    fig_all, ax_all = plt.subplots(1,1,  figsize = (8, 6) )
    fig_shear, ax_shear = plt.subplots(1,1,  figsize = (8, 6) )
    fig_v_c, ax_v_c = plt.subplots(1,1,  figsize = (8, 6) )

    for pickled_Uc, ax in zip( pickled_Uc_list, axes ):

        f = pickled_Uc.\
                replace( 'Uc_data_Values_', '').replace('.p','').\
                replace('_mfr','')

        if 'z00' in f and not 'STE' in f:
            x_bl_loc = 38
        elif 'z05' in f:
            x_bl_loc = 18
        elif 'z10' in f or 'STE' in f:
            x_bl_loc = -2

        color, marker, cmap = get_color_and_marker( 
            f
        )

        plot_config = {
            'marker'          : marker,
            'markeredgewidth' : markeredgewidth,
            'markerfacecolor' : markerfacecolor,
            'markeredgecolor' : color,
            'markersize'      : markersize*1.5,
            'mew'             : mew,
            'color'           : color,
        }

        bl_parameters = return_bl_parameters( f, [x_bl_loc] )

        Uc_df = read_pickle( join( root, pickled_Uc ) )
        Um_df = read_pickle( join( root, pickled_Uc.\
                                  replace( 
                                      'Uc_data_Values_', 
                                      'WallNormalCorrelation_Values_'
                                  ) 
                                 )
                           )[ ['Um', 'y1_bl'] ].drop_duplicates()

        Um_df = Um_df.sort_values( by = 'y1_bl' )

        Uc_df = Uc_df[ 
            ( Uc_df.Uc != 0 ) & ( Uc_df.y_bl < 1. ) & \
            ( Uc_df.pct_of_series > 25 )
        ]

        popt, pcov = curve_fit( 
            log_law, 
            Uc_df.y_real, 
            Uc_df.Uc,
        )

        x_law = log_law( 
                Uc_df.y_real, popt[0], popt[1] 
            ) / bl_parameters.Ue.values[0]
        y_law = Uc_df.y_real / bl_parameters.delta_99.values[0]

        mean_Uc = Uc_df.groupby( Uc_df.y_bl ).Uc.mean().values / \
                bl_parameters.Ue.values[0] 
        std_Uc = Uc_df.groupby( Uc_df.y_bl ).Uc.std().values / \
                bl_parameters.Ue.values[0] 

        y_bl = Uc_df.y_bl.unique() 

        int_Uc = simps(
            mean_Uc,
            Uc_df.y_real.unique()
        )

        if print_integration:
            print Uc_df.case.unique(), int_Uc, Uc_df.y_real.min()


        # spline parameters
        s    = 1.0 # smoothness parameter
        k    = 5   # spline order
        nest = -1  # estimate of number of knots needed (-1 = maximal)

        # find the knot points
        tckp,u = splprep([mean_Uc,y_bl],s=s,k=k,nest=nest)

        # evaluate spline, including interpolated points
        xnew,ynew = splev(
            linspace(0,1,400),
            tckp
        )

        df_spline = DataFrame(
            data = {
                'Uc_spline_normed': xnew,
                'Uc_spline':        xnew * bl_parameters.Ue.values[0],
                'y_bl':             ynew,
            } 
        )

        df_spline.to_pickle( join( root, 'Uc_Spline_{0}.p'.format( f ) ) )

        if pickled_Uc == pickled_Uc_list[0]:
            for um, y in zip( x_law, y_law ):

                xd = df_spline[ 
                    df_spline.y_bl == find_nearest( y, df_spline.y_bl.values ) 
                ].Uc_spline_normed.values[0]
                yd = df_spline[ 
                    df_spline.y_bl == find_nearest( y, df_spline.y_bl.values ) 
                ].y_bl.values[0]

                diff = um - xd

                ax.annotate(
                    '',
                    xy = ( xd - diff * 8 , yd ),
                    xytext = ( xd, yd ),
                    horizontalalignment = 'center',
                    verticalalignment   = 'center',
                    arrowprops          = dict(
                        facecolor = 'k',
                        edgecolor = 'k',
                        shrink    = 0.2,
                        width     = 1.0,
                        headwidth = 8.0,
                    )
                )

        ax.plot(
            mean_Uc,
            y_bl,
            ls    = '',
            zorder = 100,
            **plot_config
        )
        ax_all.plot(
            mean_Uc,
            y_bl,
            ls    = '',
            zorder = 100,
            **plot_config
        )

        ax.fill_betweenx(y_bl, 
                        mean_Uc - std_Uc,
                        mean_Uc + std_Uc,
                        facecolor='gray', 
                        alpha=0.3,
                       )

        #ax.plot(
        #    xnew,
        #    ynew,
        #    c = 'k',
        #    zorder = 800,
        #    lw = 4
        #)

        ax.plot(
            x_law, 
            y_law,
            ls = '--',
            lw     = 7.0,
            color  = 'w',
            zorder = 200,
        )

        ax.plot(
            log_law( 
                Uc_df.y_real, popt[0], popt[1] 
            ) / bl_parameters.Ue.values[0], 
            Uc_df.y_real / bl_parameters.delta_99.values[0],
            zorder = 400,
            ls = '--',
            lw     = 4.0,
            color  = 'k'
        )

        z_loc = findall( 'z[0-9][0-9]',pickled_Uc )[0].replace( 'z', '' )
        dev = findall( 'S[TER201r]+',pickled_Uc )[0]
        if not U_mean_df.empty:
            cols = U_mean_df.columns
            for col in cols.levels[0]:
                if z_loc in col and dev in col:
                    ax.plot(
                        U_mean_df[ col ][ 'u' ],
                        U_mean_df[ col ][ 'ybl' ],
                        zorder = 100,
                        c = 'k',
                        ls = '-',
                        lw = 3
                    )

                    ax.plot(
                        U_mean_df[ col ][ 'u' ],
                        U_mean_df[ col ][ 'ybl' ],
                        zorder = 50,
                        c = 'w',
                        ls = '-',
                        lw = 6
                    )

                    dudy = ( npdiff( U_mean_df[ 
                            U_mean_df[col]['ybl'] > 0.2 
                        ][ col ][ 'u_real' ] ) / \
                        npdiff( U_mean_df[ 
                            U_mean_df[col]['ybl'] > 0.2 
                        ][ col ][ 'y_real' ] ) )**2
                    
                    ax_shear.plot(
                        dudy ,
                        U_mean_df[ 
                            U_mean_df[col]['ybl'] > 0.2 
                        ][ col ][ 'ybl' ][:-1],
                        ls = '',
                        **plot_config
                    )

                    int_shear = simps(
                        dudy,
                        U_mean_df[ 
                            U_mean_df[col]['ybl'] > 0.2 
                        ][ col ][ 'y_real' ][:-1],
                    )

                    print "Shear", col, int_shear

                    
                    v_uc_df = concat( [ 
                        U_mean_df[col][['ybl','v_real']].dropna().\
                        set_index('ybl'),
                        DataFrame( data = { 
                            'uc':mean_Uc * bl_parameters.Ue.values[0] 
                        }, index = y_bl )
                    ], axis = 1 )

                    v_uc_df = v_uc_df.apply( Series.interpolate ).dropna()

                    ax_v_c.plot(
                        (v_uc_df.v_real / v_uc_df.uc )[::3]**2 * 100,
                        v_uc_df.index.values[::3],
                        ls = '',
                        **plot_config
                    )

                    int_v_uc = simps(
                        (v_uc_df.v_real / v_uc_df.uc )[::3]**2 ,
                        v_uc_df.index.values[::3],
                    )
                    print "v_uc", col, int_v_uc

        ax.set_xlabel( r'$ u_c / u_e $' , fontsize = fontsize)

        if pickled_Uc == pickled_Uc_list[0]:
            ax.set_xticks( arange( 0.0, 1.2, 0.2 ) )
        else:
            ax.set_xticks( arange( 0, 1.2, 0.2 ) )

        ax.tick_params(axis='both', which='major', labelsize=fontsize)
        ax.set_ylim( 0 , 1 )
        ax.set_xlim( 0.4 , 1. )

    ax_all.tick_params(axis='both', which='major', labelsize=fontsize)
    ax_all.set_xlabel( r'$ u_c / u_e $' , fontsize = fontsize)
    ax_all.set_ylabel( r'$ y / \delta $' , fontsize = fontsize)
    ax_all.set_ylim( 0, 1 )
    ax_all.set_yticks( [0, 0.2, 0.4, 0.6, 0.8, 1] )

    ax_shear.set_xlabel( 
        r'$ \left[ \partial \overline{u} / \partial y \right]^2$'\
        +r'   [$1/\textrm{s}^2$]' , 
        fontsize = fontsize
    )
    ax_shear.set_ylabel( r'$ y / \delta $' , fontsize = fontsize)
    ax_shear.set_ylim( 0, 1 )
    ax_shear.set_yticks( [0, 0.2, 0.4, 0.6, 0.8, 1] )
    ax_shear.tick_params(axis='both', which='major', labelsize=fontsize)

    ax_v_c.set_xlabel( r'$ \left[ \overline{v} / u_c\right]^2\,\times10^{-2} $', 
                      fontsize = fontsize)
    ax_v_c.set_ylabel( r'$ y / \delta $' , fontsize = fontsize)
    ax_v_c.set_ylim( 0, 1 )
    ax_v_c.set_xlim( 0, 1.3 )
    ax_v_c.set_yticks( [0, 0.2, 0.4, 0.6, 0.8, 1] )
    ax_v_c.tick_params(axis='both', which='major', labelsize=fontsize)

    #im = plt.imread( get_sample_data( schematic  ) )
    #newax = fig_all.add_axes([0.2, 0.5, 0.3, 0.3], 
    #                         anchor = 'SW', zorder=100)
    #newax.imshow(im)
    #newax.axis('off')

    #newax = fig_v_c.add_axes([0.6, 0.15, 0.3, 0.3], 
    #                         anchor = 'SW', zorder=100)
    #newax.imshow(im)
    #newax.axis('off')

    axes[0].set_ylabel( r'$ y / \delta $' , fontsize = fontsize)
    axes[0].text( 0.7, 0.55, r"$\overline{u} / u_e$", ha = 'right', 
                 fontsize = fontsize )
    fig.tight_layout()
    fig.savefig( join( output_root, 'Uc.png' )) #, bbox_inches = 'tight' )
    fig_all.savefig( join( output_root, 'Uc_all.png' ), bbox_inches = 'tight' )
    fig_shear.savefig( join( output_root, 'U_shear.png' ), bbox_inches = 'tight' )
    fig_v_c.savefig( join( output_root, 'v_uc.png' ), bbox_inches = 'tight' )

def is_outlier(df, var, thresh=3.5):
    """
    Returns a boolean array with True if points are outliers and False 
    otherwise.

    Parameters:
    -----------
        points : An numobservations by numdimensions array of observations
        thresh : The modified z-score to use as a threshold. Observations with
            a modified z-score (based on the median absolute deviation) greater
            than this value will be classified as outliers.

    Returns:
    --------
        mask : A numobservations-length boolean array.

    References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor. 
    """
    from numpy import sum, sqrt
    from numpy import median as npmedian
    #from pandas import read_pickle, DataFrame
    #from os.path import isfile

    points = df[ var ].values

    if len(points.shape) == 1:
        points = points[:,None]

    median            = npmedian(points, axis=0)
    diff              = sum((points - median)**2, axis=-1)
    diff              = sqrt(diff)
    med_abs_deviation = npmedian(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    #if isfile( 'outliers.p' ):
    #    df_outliers = read_pickle( 'outliers.p')
    #else:
    #    df_outliers = DataFrame()

    data = {}
    for col in df.columns:
        if len( df[ col ].unique() ) == 1:
            data[ col ] = df[ col ].unique()[0] 
    data[ 'n_outliers' ] = sum( modified_z_score > thresh )
    data[ 'n_time_steps' ] = len( df.t.unique() )

    #df_outliers = df_outliers.append(
    #    DataFrame( data, index = [0] ), ignore_index = True
    #)

    #df_outliers.to_pickle( 'outliers.p' )

    return modified_z_score > thresh

def calculate_Uc(df):
    from math        import pi
    from scipy.stats import linregress

    r_value      = 0
    truncated_df = df[ ['f', 'phi'] ].copy()
    truncated_df = truncated_df.sort_values( by = [ 'f' ] )
    truncated_df = truncated_df.reset_index( drop = True )
    consider     = truncated_df.index[-2]

    while r_value**2 < 0.97:

        truncated_df = truncated_df.sort_values( by = [ 'f' ] ).\
                ix[:consider].\
                reset_index(drop=True)

        slope, intercept, r_value, p_value, std_err = linregress(
            truncated_df['phi'],
            truncated_df['f'],
        )

        consider -= 1
        if len(truncated_df) < 3:
            return 0, 0

    used_pct_of_series = 100 * len( truncated_df ) / len( df )
    
    Uc = 2 * pi * slope * df.xi.unique()[0] / 1000.

    return Uc, used_pct_of_series

def sort_block_coordinates_by_height_and_streamwise( coord_df , plot = '',
                                                   truncate = False ,
                                                   case = '' ):
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import cm
    from numpy             import linspace, meshgrid, zeros, diff, arange
    from numpy             import round
    from math              import degrees, radians
    from os.path           import split

    def find_angle( xy1, xy2 ):
        from math import atan, pi

        angle = atan(
            ( xy2[1] - xy1[1] ) / ( xy2[0] - xy1[0] )
        )
        if angle < 0:
            angle += pi / 2.

        return angle

    def find_grid_rotation( df ):
        from numpy                        import sqrt, array
        from raw_data_processing_routines import find_nearest
        from math                         import pi

        std     = 1
        epsilon = 0
        while std > 0.01:
            # Get central point ################################################
            x_centr = find_nearest( 
                ( df.x.max() + df.x.min() ) / ( 2. + epsilon ), df.x.values 
            )

            y_centr = df[ df.x == x_centr ].y.values[0]

            # Get distance to other points #####################################
            df[ 'dist' ] = sqrt( ( df.x - x_centr )**2 + ( df.y - y_centr )**2 )

            df = df.sort_values( by = 'dist' )
            df_trunc = df.iloc[ 0:5 ]

            angles = []
            xy1 = ( df_trunc.iloc[0].x, df_trunc.iloc[0].y )
            for i in range(1,5):
                xy2 = ( df_trunc.iloc[i].x, df_trunc.iloc[i].y )
                angles.append(
                    find_angle( xy1, xy2 )
                )

            angles = array( angles )
            angle = angles.mean()

            std = angles.std()
            
            epsilon -= 0.01

        if abs( angle ) > pi / 4.:
            angle = - angle + pi / 2.

        return angle

    def rotate_df( df, angle ):
        from math import sin, cos

        x = df.x
        y = df.y

        df.x = round( x * cos( angle ) - y * sin( angle ), 4 )
        df.y = round( x * sin( angle ) + y * cos( angle ), 4 )

        return df

    def find_indexes( coord_df , rotation_angle, correction = 0 ):

        coord_df = rotate_df( coord_df, rotation_angle + correction )
        coord_df[ 'height_ix' ]      = 0
        coord_df[ 'stream_ix' ]      = 0
        coord_df[ 'outside_fov' ]    = False
        coord_df[ 'rotation_angle' ] = rotation_angle

        height_cnt = 1
        stream_cnt = 0

        coord_df = coord_df.sort_values( by = ['x','y'] )
        unique_x = coord_df.x.unique()
        
        for x,dx in zip( unique_x[:-1], diff( unique_x ) ):
            coord_df.loc[ coord_df.x == x , 'stream_ix' ] = stream_cnt
            if dx > 0.09:
                stream_cnt -= 1

        coord_df = coord_df.sort_values( by = ['y','x'] )
        unique_y = coord_df.y.unique()

        import numpy
        numpy.set_printoptions(threshold=numpy.nan)
        for y,dy in zip( unique_y[:-1], diff( unique_y ) ):
            coord_df.loc[ coord_df.y == y , 'height_ix' ] = height_cnt
            if abs( dy ) > 0.020:
                height_cnt -= 1

        coord_df.height_ix = coord_df.height_ix - coord_df.height_ix.min()
        coord_df.stream_ix = coord_df.stream_ix - coord_df.stream_ix.min()

        return coord_df

    def truncate_grid( coord_df ):

        coord_df[ 'stream_cnt' ] = 0
        for h_ix in sorted( coord_df.height_ix.unique() ):
            coord_df.loc[ coord_df.height_ix == h_ix, 'stream_cnt' ] = \
                   coord_df[ coord_df.height_ix == h_ix].stream_ix.shape[0]

        max_stream_cnt = coord_df.stream_cnt.value_counts().index[0]

        coord_df = coord_df[
            coord_df.stream_cnt >= max_stream_cnt - 3
        ]
        
        coord_df.height_ix = coord_df.height_ix - coord_df.height_ix.min()
        coord_df.stream_ix = coord_df.stream_ix - coord_df.stream_ix.min()

        return coord_df


    angles_dict = {
        'Sr20R21_a-12_p6_U20_z05_tr.hdf5' : -6
    }

    case_name = split( case )[1]
    if case_name in angles_dict.keys():
        rotation_angle = radians( angles_dict[ case_name ] )
    else:
        rotation_angle = find_grid_rotation( coord_df )

    print "     rotating grid by {0:.2f} degrees to align points"\
            .format( degrees( - rotation_angle ) )

    coord_df  = find_indexes( 
        coord_df , 
        - rotation_angle 
    )

    if truncate == True:
        coord_df = truncate_grid( coord_df )

    stream_ixs = coord_df.stream_ix.unique()
    height_ixs = coord_df.height_ix.unique()

    if plot:
        fig, axes = plt.subplots( 1, 3 )

        color=cm.flag(linspace(0,1,max(stream_ixs)+1))
        for stream_ix, c_ix in zip( stream_ixs, range( len( stream_ixs ) ) ):
            axes[0].scatter(
                coord_df[ ( coord_df.stream_ix == stream_ix ) & \
                         ~coord_df.outside_fov ].x,
                coord_df[ ( coord_df.stream_ix == stream_ix ) & \
                         ~coord_df.outside_fov ].y,
                c = color[ c_ix ]
            )

        color=cm.flag(linspace(0,1,max(height_ixs)+1))
        for height_ix, c_ix in zip( height_ixs, range( len( height_ixs ) ) ):
            axes[1].scatter(
                coord_df[ ( coord_df.height_ix == height_ix ) & \
                         ~coord_df.outside_fov ].x,
                coord_df[ ( coord_df.height_ix == height_ix ) & \
                         ~coord_df.outside_fov ].y,
                c = color[ c_ix ]
            )

        X, Y = meshgrid( 
            linspace( 
                coord_df.stream_ix.min(), 
                coord_df.stream_ix.max(),
                coord_df.stream_ix.max() + 1
            ), 
            linspace( 
                coord_df.height_ix.min(), 
                coord_df.height_ix.max(), 
                coord_df.height_ix.max() + 1 
            ), 
        )

        U = zeros( X.shape )
        V = zeros( X.shape )

        cnt = 0
        for row in coord_df.iterrows():
            U[ row[1].height_ix , row[1].stream_ix ] = row[1].u
            V[ row[1].height_ix , row[1].stream_ix ] = row[1].v

        #axes[2].tricontourf( 
        #    coord_df.x, 
        #    coord_df.y, 
        #    coord_df.u, 
        #    #levels = arange(0,1.1,0.1),
        #    cmap   = plt.cm.Blues
        #)

        color=cm.Blues(linspace(0,cnt,10))
        stride = 3

        axes[2].contourf(
            X,Y,U/U.max(),
            cmap = cm.Reds,
            levels = arange(0,1,0.1)
        )
        axes[2].quiver(
            X[::stride][::stride],
            Y[::stride][::stride],
            U[::stride][::stride],
            V[::stride][::stride],
            scale = 100
        )
        axes[0].set_aspect( 'equal' )
        axes[1].set_aspect( 'equal' )
        axes[2].set_aspect( 'equal' )
        plt.savefig( plot , bbox_inches = 'tight' )

    return coord_df

def plot_phi( root = '' ):
    import matplotlib.pyplot as plt
    from matplotlib import rc
    from numpy      import meshgrid, arange, linspace
    from math       import pi
    from pandas     import read_pickle
    from os.path    import join
    from os         import listdir
        
    if not root:
        root = '/home/carlos/Documents/PhD/Articles/Article_3/' + \
                'Scripts/time_resolved/Results_v2'

    rc('text',usetex=True)
    rc('font',weight='normal')

    fontsize = 30

    rc('font',family='serif', serif='Linux Libertine')

    phi_pivots = [
        f for f in listdir( root ) if f.startswith( 'StreamwisePhi' ) and \
        f.endswith('.p')
    ]
    
    for pivot in phi_pivots:

        coh_df_xi_pivot = read_pickle( join( root, pivot ) )

        case = pivot.replace('.p','')

        X   = coh_df_xi_pivot.columns.values/1000.
        Y   = coh_df_xi_pivot.index.values
        x,y = meshgrid(X, Y)
        Z   = coh_df_xi_pivot.values

        fig_coh, ax_coh = plt.subplots( 1, 2, figsize = ( 13, 4 ) )

        xis = coh_df_xi_pivot.index.values
        line_color = plt.cm.Oranges( linspace(0, 1, len(xis )) )
        for xi, c in zip( xis, line_color ):
            ax_coh[1].plot( 
                coh_df_xi_pivot.T[ xi ].index.values[:-2]/1000.,
                coh_df_xi_pivot.T[ xi ].values[:-2]/pi,
                c = c
            )
        ax_coh[1].annotate(
            ' ',
            xy                  = (2.000, 4.),
            xytext              = (3.000, 1/2.),
            horizontalalignment = 'left',
            verticalalignment   = 'top',
            arrowprops          = dict(facecolor='black', shrink=0.025),
        )
        ax_coh[1].text(
            2.400,
            3.4,
            r'$\xi_x$',
            fontsize = fontsize
        )
        ax_coh[1].set_xlabel( r'$f$ [kHz]' , fontsize = fontsize)
        #ax_coh[1].set_ylabel( r'$\phi$ [rad]' , fontsize = fontsize)
        ax_coh[1].set_yticks( arange(0, 5.5, 1.0)  )
        ax_coh[1].set_yticklabels( [ "{0:.0f}".format(d) for d in \
            arange(0, 5.5, 1.0) ]
        )
        ax_coh[1].set_ylim( bottom = 0 )

        CS = ax_coh[0].contourf( X, Y, Z/pi, levels = arange(0, 5.5, 0.5) , 
                        cmap = plt.cm.Blues )

        cbar = plt.colorbar(CS, ax = ax_coh[0], shrink = 1.0)

        cbar.ax.set_ylabel(r'$\phi$ [rad$/ \pi$]', fontsize = fontsize)
        cbar.ax.tick_params( labelsize = fontsize )

        ax_coh[0].set_xlabel( r'$f$ [kHz]' , fontsize = fontsize)
        ax_coh[0].set_ylabel( r'$\xi_x$ [m]$\times10^{-3}$' , 
                             fontsize = fontsize)

        ax_coh[0].tick_params(axis='both', which='major', labelsize=fontsize)
        ax_coh[1].tick_params(axis='both', which='major', labelsize=fontsize)
        fig_coh.subplots_adjust( wspace = 0.1 )
        fig_coh.savefig( 
            join( output_root, '{0}.png'\
                 .format( case ) ),
            bbox_inches = 'tight' )

def do_the_streamwise_coherence_analysis( pickled_coh_data , root = '', 
                                         overwrite = False):
    from pandas                          import read_pickle, DataFrame
    from os.path                         import split, join, isfile
    from raw_data_processing_routines    import find_nearest
    from progressbar                     import ProgressBar,Percentage
    from progressbar                     import Bar,ETA,SimpleProgress

    if not root:
        root = '/home/carlos/Documents/PhD/Articles/Article_3/' + \
                'Scripts/time_resolved/Results_v2'

    case = split( pickled_coh_data )[1]\
            .replace('.p','').replace("StreamwiseCoherence_","") 

    if isfile( 
        join( root, 'Uc_data_{0}.p'.format( 
            case.replace("StreamwiseCoherence_","") 
        ) )
    ) and not overwrite:
        return 0

    coh_df        = read_pickle( join( root, pickled_coh_data ) )

    y_bl_locs_to_plot = [ 0.3, 0.6, 0.9 ]

    y_bl_available   = coh_df.y1_bl.unique()
    y_real_available = coh_df.y1.unique()

    nearest_y_bl_available = [ 
        find_nearest( n, y_bl_available ) for n in y_bl_locs_to_plot 
    ]

    y_ix_to_plot = coh_df[ 
        coh_df.y1_bl.isin( nearest_y_bl_available )
    ].height_ix.unique()

    Uc_df = DataFrame()

    progress = ProgressBar(
         widgets=[
             Bar(),' ',
             Percentage(),' ',
             ETA(), ' ( streamwise location ',
             SimpleProgress(),' )'], 
         maxval= len( y_bl_available ) 
         ).start()

    cnt = 0
    for y_ix, y_bl, y_re in zip( 
        coh_df.height_ix.unique(), y_bl_available, y_real_available
    ):
        # Do a plot of the coherence as a function of xi #######################
        for xi in coh_df[ coh_df.height_ix ==  y_ix  ]\
                  .xi.unique():

            coh_df.loc[ 
                ( coh_df.height_ix == y_ix ) & \
                ( coh_df.xi == xi ), 
                'phi' ] = remove_angle_jumps( 
                    coh_df[ 
                        ( coh_df.height_ix == y_ix ) & \
                        ( coh_df.xi == xi )
                    ]
                )

            Uc, pct_of_series = calculate_Uc( 
                coh_df[ 
                    ( coh_df.height_ix == y_ix ) & \
                    ( coh_df.xi == xi )
                ]
            )

            Uc_df = Uc_df.append(
                DataFrame( data = {
                    'Uc':            Uc,
                    'pct_of_series': pct_of_series,
                    'xi':            xi,
                    'y_real':        y_re,
                    'y_bl':          y_bl,
                    'case':          case
                }, index = [0]), ignore_index = True
            )

        if y_ix in y_ix_to_plot:

            coh_df_xi_pivot = coh_df[ 
                coh_df.height_ix == y_ix 
            ].pivot( 'xi', 'f', 'phi' )

            coh_df_xi_pivot.to_pickle(
                join( 
                    root, 
                    '{0}_ybl_{1:.1f}.p'.format( 
                        case , 
                        y_bl 
                    ) 
                )
            )

        cnt += 1
        progress.update( cnt )

    Uc_df.to_pickle( 
        join( root, 'Uc_data_{0}.p'.format( 
            case.replace("StreamwiseCoherence_","") 
        )
        ) 
    )
    progress.finish()

def get_streamwise_coherence_and_correlation( hdf_name, overwrite = False ,
                                            root = ''):
    from scipy.signal import csd
    from matplotlib   import rc
    from numpy        import corrcoef, sqrt, angle, round
    from pandas       import read_hdf, DataFrame
    from os.path      import split, join, isfile
    from progressbar  import ProgressBar,Percentage
    from progressbar  import Bar,ETA,SimpleProgress
    from scipy.signal import coherence

    rc('text',usetex=True)
    rc('font',weight='normal')

    rc('font',family='serif', serif='Linux Libertine')

    if not root:
        root = '/home/carlos/Documents/PhD/Articles/Article_3/' + \
                'Scripts/time_resolved/Results_v2'

    case = split( hdf_name )[1].replace('.hdf5','')

    case_coords = read_hdf( 
        hdf_name , where = ( 't=0' ),
        columns = [ 'x', 'y' , 'near_x_2h', 'near_y_delta', 'u', 'v'],
                          )
    
    case_coords = case_coords.sort_values( 
        by = [ 'y', 'x' ] 
    )

    if isfile( join( root, 'StreamwiseCorrelation_Values_{0}.p'.format( case ) )
             ) and not overwrite:
        return 0

    case_coords = sort_block_coordinates_by_height_and_streamwise( 
        case_coords,
        case     = case,
        plot     = "CoordExtraction_"+case+".png" ,
        truncate = True
    )

    print "   Performing the streamwise coherence and "
    print "   correlation analysis of {0}".format(case)

    cnt = 0
    progress = ProgressBar(
         widgets=[
             Bar(),' ',
             Percentage(),' ',
             ETA(), ' ( streamwise location ',
             SimpleProgress(),' )'], 
         maxval= len( case_coords ) 
         ).start()

    corr_df = DataFrame()
    coh_df  = DataFrame()

    for height_ix in case_coords.height_ix.unique():

        avaiable_xy = case_coords[ 
            case_coords.height_ix == height_ix 
        ][ ['x','y', 'near_x_2h', 'near_y_delta' , 'stream_ix'] ]

        ref_coord = avaiable_xy[ avaiable_xy.x == avaiable_xy.x.max() ]

        case_x_1 = read_hdf(
            hdf_name, 
            where = ( 
                'near_x_2h = {0} & near_y_delta = {1}'.format( 
                    ref_coord.near_x_2h.values[0],
                    ref_coord.near_y_delta.values[0],
                )
            ), 
            columns = [
                'near_x_2h',
                'near_y_delta',
                't',
                'u',
                'x',
                'y'
            ]
        )

        case_x_1.near_x_2h    = round( case_x_1.near_x_2h, 2 )
        case_x_1.near_y_delta = round( case_x_1.near_y_delta, 2 )

        case_x_1[ 'case' ] = case

        case_x_1.loc[ 
            is_outlier( case_x_1, 'u' , thresh = outlier_thresh) , 'u'
        ] = case_x_1[ ~is_outlier( case_x_1, 'u' , thresh = outlier_thresh) ]\
                .u.mean()

        case_x_1 = case_x_1.fillna( 0 )
        case_x_1 = case_x_1.reset_index( drop = True )

        case_x_1[ 'u_prime' ] = case_x_1.u - case_x_1.u.mean()

        for xcorr_coord in avaiable_xy.iterrows():

            case_x_2 = read_hdf(
                hdf_name, 
                where = ( 
                    'near_x_2h = {0} & near_y_delta = {1}'.format( 
                        xcorr_coord[1].near_x_2h,
                        xcorr_coord[1].near_y_delta,
                    )
                ), 
                columns = [
                    'near_x_2h',
                    'near_y_delta',
                    't',
                    'u',
                    'x',
                    'y'
                ]
            ).sort_values( by = 't' ).reset_index( drop = True )

            case_x_2.near_x_2h    = round( case_x_2.near_x_2h, 2 )
            case_x_2.near_y_delta = round( case_x_2.near_y_delta, 2 )

            case_x_2[ 'case' ] = case

            if case_x_2.empty:
                continue

            case_x_2.loc[ 
                is_outlier( case_x_2, 'u' , thresh = outlier_thresh) , 'u'
            ] = case_x_2[ 
                ~is_outlier( case_x_2, 'u' , thresh = outlier_thresh) 
            ].u.mean()

            case_x_2 = case_x_2.fillna( 0 )
            case_x_2 = case_x_2.reset_index( drop = True )

            case_x_2[ 'u_prime' ] = case_x_2.u - case_x_2.u.mean()

            len_diff = len( case_x_2 ) - len( case_x_1 )
            if len_diff != 0:
                if len_diff < 3 and len_diff > 0:
                    case_x_2 = case_x_2.iloc[ 0 : len( case_x_1 ) ]
                elif len_diff > -3 and len_diff < 0:
                    case_x_1 = case_x_1.iloc[ 0 : len( case_x_2 ) ]

                else:
                    print "   Warning, uneven case lengths:"
                    print "     {0} vs {1}".format(
                        len( case_x_1 ), len( case_x_2 )
                    ) 
                    continue

            xcorr = corrcoef( case_x_1.u_prime, case_x_2.u_prime )[1,0]

            f, Pxy = csd(
                case_x_1.u,
                case_x_2.u,
                fs      = fs,
                nperseg = nperseg,
            )

            f, gamma = coherence(
                case_x_2.u,
                case_x_1.u,
                fs      = fs,
                nperseg = nperseg,
            )

            phi = angle( Pxy )

            x1_real = round( case_x_1.x.unique()[0], 3 )
            x2_real = round( case_x_2.x.unique()[0], 3 )
            y1_real = round( case_x_1.y.unique()[0], 3 )
            y2_real = round( case_x_2.y.unique()[0], 3 )

            xi = sqrt( (x2_real - x1_real)**2 + (y2_real - y1_real)**2 )

            coh_df_loc = DataFrame( data = {
                'gamma' : gamma,
                'phi'   : phi,
                'f'     : f,
            } )

            coh_df_loc[ 'y1_bl'     ] = case_x_1.near_y_delta.unique()[0]
            coh_df_loc[ 'y2_bl'     ] = case_x_2.near_y_delta.unique()[0]
            coh_df_loc[ 'y1'        ] = y1_real
            coh_df_loc[ 'y2'        ] = y2_real
            coh_df_loc[ 'x1'        ] = x1_real
            coh_df_loc[ 'x2'        ] = x2_real
            coh_df_loc[ 'x1_2h'     ] = case_x_1.near_x_2h.unique()[0]
            coh_df_loc[ 'x2_2h'     ] = case_x_2.near_x_2h.unique()[0]
            coh_df_loc[ 'xi'        ] = xi
            coh_df_loc[ 'Um1'       ] = case_x_1.u.mean()
            coh_df_loc[ 'Um2'       ] = case_x_2.u.mean()
            coh_df_loc[ 'height_ix' ] = height_ix
            coh_df_loc[ 'stream_ix' ] = xcorr_coord[1].stream_ix

            coh_df = coh_df.append( coh_df_loc, ignore_index = True )

            corr_df = corr_df.append(
                DataFrame( data = {
                    'xcorr':     xcorr,
                    'y1_bl':     case_x_1.near_y_delta.unique()[0],
                    'y2_bl':     case_x_2.near_y_delta.unique()[0],
                    'y1':        y1_real,
                    'y2':        y2_real,
                    'x1':        x1_real,
                    'x2':        x2_real,
                    'x1_2h':     case_x_1.near_x_2h.unique()[0],
                    'x2_2h':     case_x_2.near_x_2h.unique()[0],
                    'xi':        xi,
                    'Um1':       case_x_1.u.mean(),
                    'Um2':       case_x_2.u.mean(),
                    'height_ix': height_ix,
                    'stream_ix': xcorr_coord[1].stream_ix
                } , index = [0] ), ignore_index = True
            )

            cnt += 1
            progress.update( cnt )

    progress.finish()

    corr_df.to_pickle( join( root, 'StreamwiseCorrelation_Values_{0}.p'.\
                      format( case ) )
                     )
    coh_df.to_pickle( join( root, 'StreamwiseCoherence_Values_{0}.p'.\
                      format( case ) )
                     )

def get_vertical_mean_values( hdf_names , root = '', plot_individual = True,
                            overwrite = False ):
    from numpy                           import nan
    from pandas                          import DataFrame, read_hdf
    from os.path                         import split, join, isfile
    from raw_data_processing_routines    import find_nearest
    from progressbar                     import ProgressBar,Percentage
    from progressbar                     import Bar,ETA,SimpleProgress

    if not root:
        root = '/home/carlos/Documents/PhD/Articles/Article_3/' + \
                'Scripts/time_resolved/Results_v2'

    for hdf_name in [ join(root, f) for f in hdf_names]:

        case = split( hdf_name )[1].replace('.hdf5','')

        if isfile(
            join( root, 'WallNormal_MeanValues_{0}.p'.format( case ) ) ):
            if not overwrite:
                continue

        print "  Processing vertical mean values for {0}".format( case )

        case_coords = read_hdf( 
            hdf_name , 
            where = ( 't=0 & near_y_delta > 0.0' ),
            columns = [ 'x', 'y' , 'near_x_2h', 'near_y_delta', 'u', 'v'],
        ).sort_values( by = [ 'x', 'y' ] )
        
        print "  from y/delta = {0:.2f} to {1:.2f}".format( 
            case_coords.near_y_delta.min(),
            case_coords.near_y_delta.max(),
        )

        #if 'z00' in hdf_name and not 'STE' in hdf_name:
        #    case_coords = case_coords[ case_coords.near_y_delta > 0.15 ]

        case_coords = sort_block_coordinates_by_height_and_streamwise( 
            case_coords , case = hdf_name, 
            plot = "CoordExtraction_VertCorr_"+case+".png" , truncate = False
        )

        color, marker, cmap = get_color_and_marker( 
            case
        )

        available_x = sorted( 
            case_coords[ case_coords.height_ix == 0 ].near_x_2h.unique() 
        )

        # Only evaluate the location over the TE # #########################
        if 'z10' in hdf_name :
            x = find_nearest( -2 / 40. , available_x )
        elif  'STE' in hdf_name:
            x = find_nearest( 2 / 40. , available_x )
        elif 'z05' in hdf_name:
            x = find_nearest( 18 / 40. , available_x )
        elif 'z00' in hdf_name:
            x = find_nearest( 38 / 40. , available_x )

        central_stream_ix = case_coords[ 
            ( case_coords.near_x_2h == x ) & \
            ( case_coords.height_ix == 0 )
        ].stream_ix.unique()[0]

        mean_df = DataFrame()

        height_ixs = case_coords[ 
            case_coords.stream_ix == central_stream_ix 
        ].height_ix.unique()

        progress = ProgressBar(
             widgets=[
                 Bar(),' ',
                 Percentage(),' ',
                 ETA(), ' ( streamwise location ',
                 SimpleProgress(),' )'], 
             maxval= len( height_ixs )
             ).start()

        cnt = 0

        for y1_ix in height_ixs:
            x1_loc = case_coords[
                ( case_coords.stream_ix == central_stream_ix ) & \
                ( case_coords.height_ix == y1_ix )
            ].near_x_2h.unique()[0]

            y1_loc = case_coords[
                ( case_coords.stream_ix == central_stream_ix ) & \
                ( case_coords.height_ix == y1_ix )
            ].near_y_delta.unique()[0]

            case_y_1 = read_hdf(
                hdf_name, 
                where = ( 
                    'near_x_2h = {0} & near_y_delta = {1}'.format( 
                        x1_loc ,
                        y1_loc
                    ) 
                )
            ).sort_values( by = 't' ).reset_index( drop = True )

            case_y_1[ 'case' ] = case

            case_y_1.loc[ 
                is_outlier( case_y_1, 'v' , thresh = outlier_thresh) , 'v'
            ] = nan 

            y1_real = case_y_1.y.unique()[0]

            mean_df = mean_df.append(
                DataFrame( data = {
                    'y1_bl':      y1_loc,
                    'y1_real':    y1_real,
                    'case':       case,
                    'height_ix1': y1_ix,
                    'Um':         case_y_1.u.mean(),
                    'w':          case_y_1.w.mean(),
                    'v':          case_y_1.v.mean(),
                    'urms':       case_y_1.u.std(),
                    'vrms':       case_y_1.v.std(),
                } , index = [0] ), 
                ignore_index = True
            )

            cnt += 1
            progress.update( cnt )


        progress.finish()

        mean_df.to_pickle( 
            join( root, 'WallNormal_MeanValues_{0}.p'.format( case ) )
        )

def get_vertical_correlation( hdf_names , root = '', plot_individual = True,
                            overwrite = False ):
    import matplotlib.pyplot as plt
    #import seaborn as sns
    from matplotlib                      import rc
    from scipy.signal                    import coherence
    from numpy                           import meshgrid, array, nan
    from numpy                           import sqrt,isnan
    from scipy.integrate                 import simps
    from pandas                          import DataFrame, read_hdf, concat
    from os.path                         import split, join, isfile
    from raw_data_processing_routines    import find_nearest
    from progressbar                     import ProgressBar,Percentage
    from progressbar                     import Bar,ETA,SimpleProgress

    rc('text',usetex=True)
    rc('font',weight='normal')

    #sns.set_context('paper')
    #sns.set(font='serif',font_scale=3.0,style='ticks')
    rc('font',family='serif', serif='Linux Libertine')

    if not root:
        root = '/home/carlos/Documents/PhD/Articles/Article_3/' + \
                'Scripts/time_resolved/Results_v2'

    for hdf_name in [ join(root, f) for f in hdf_names]:

        case = split( hdf_name )[1].replace('.hdf5','')

        if isfile(
            join( root, 'WallNormalCorrelation_Values_{0}.p'.format( case ) ) ):
            if not overwrite:
                continue

        print "  Processing vertical coh and corr for {0}".format( case )

        case_coords = read_hdf( 
            hdf_name , 
            where = ( 't=0 & near_y_delta > 0.2' ),
            columns = [ 'x', 'y' , 'near_x_2h', 'near_y_delta', 'u', 'v'],
        ).sort_values( by = [ 'x', 'y' ] )
        
        print "  from y/delta = {0:.2f} to {1:.2f}".format( 
            case_coords.near_y_delta.min(),
            case_coords.near_y_delta.max(),
        )

        #if 'z00' in hdf_name and not 'STE' in hdf_name:
        #    case_coords = case_coords[ case_coords.near_y_delta > 0.15 ]

        case_coords = sort_block_coordinates_by_height_and_streamwise( 
            case_coords , case = hdf_name, 
            plot = "CoordExtraction_VertCorr_"+case+".png" , truncate = False
        )

        color, marker, cmap = get_color_and_marker( 
            case
        )

        available_x = sorted( 
            case_coords[ case_coords.height_ix == 0 ].near_x_2h.unique() 
        )

        # Only evaluate the location over the TE # #########################
        if 'z10' in hdf_name :
            x = find_nearest( -2 / 40. , available_x )
        elif  'STE' in hdf_name:
            x = find_nearest( 2 / 40. , available_x )
        elif 'z05' in hdf_name:
            x = find_nearest( 15 / 40. , available_x )
        elif 'z00' in hdf_name:
            x = find_nearest( 38 / 40. , available_x )

        central_stream_ix = case_coords[ 
            ( case_coords.near_x_2h == x ) & \
            ( case_coords.height_ix == 0 )
        ].stream_ix.unique()[0]

        corr_df = DataFrame()

        height_ixs = case_coords[ 
            case_coords.stream_ix == central_stream_ix 
        ].height_ix.unique()

        progress = ProgressBar(
             widgets=[
                 Bar(),' ',
                 Percentage(),' ',
                 ETA(), ' ( streamwise location ',
                 SimpleProgress(),' )'], 
             maxval= len( height_ixs )**2 
             ).start()

        cnt = 0

        for y1_ix in height_ixs:

            x1_loc = case_coords[
                ( case_coords.stream_ix == central_stream_ix ) & \
                ( case_coords.height_ix == y1_ix )
            ].near_x_2h.unique()[0]

            y1_loc = case_coords[
                ( case_coords.stream_ix == central_stream_ix ) & \
                ( case_coords.height_ix == y1_ix )
            ].near_y_delta.unique()[0]

            case_y_1 = read_hdf(
                hdf_name, 
                where = ( 
                    'near_x_2h = {0} & near_y_delta = {1}'.format( 
                        x1_loc ,
                        y1_loc
                    ) 
                )
            ).sort_values( by = 't' ).reset_index( drop = True )

            case_y_1[ 'case' ] = case

            case_y_1.loc[ 
                is_outlier( case_y_1, 'v' , thresh = outlier_thresh) , 'v'
            ] = nan 

            gamma              = []
            delta_y            = []
            delta_y_bl         = []
            int_gamma_per_freq = []

            if plot_individual:
                fig, ax = plt.subplots(1, 2, sharex=True, figsize = (15,5))

            for y2_ix in case_coords[ 
                case_coords.stream_ix == central_stream_ix 
            ].height_ix.unique():

                x2_loc = case_coords[
                    ( case_coords.stream_ix == central_stream_ix ) & \
                    ( case_coords.height_ix == y2_ix )
                ].near_x_2h.unique()[0]

                y2_loc = case_coords[
                    ( case_coords.stream_ix == central_stream_ix ) & \
                    ( case_coords.height_ix == y2_ix )
                ].near_y_delta.unique()[0]

                case_y_2 = read_hdf(
                    hdf_name, 
                    where = ( 
                        'near_x_2h = {0} & near_y_delta = {1}'.format( 
                            x2_loc ,
                            y2_loc
                        ) 
                    )
                ).sort_values( by = 't' ).reset_index( drop = True )

                case_y_2[ 'case' ] = case

                case_y_2.loc[ 
                    is_outlier( case_y_2, 'v', thresh = outlier_thresh) , 'v'
                ] = nan

                case_y_1.v = case_y_1.v.interpolate( 
                            method='spline', order=3, s=0.
                        )
                case_y_2.v = case_y_2.v.interpolate( 
                            method='spline', order=3, s=0.
                        )

                case_y_1[ 'v_prime1' ] = case_y_1.v - case_y_1.v.mean()

                case_y_2[ 'v_prime2' ] = case_y_2.v - case_y_2.v.mean()

                tmp_corr_df = concat( 
                    [ case_y_1, case_y_2 ], axis = 1 
                ).dropna()

                #xcorr = corrcoef( 
                #    tmp_corr_df.v_prime1, 
                #    tmp_corr_df.v_prime2,
                #)[1,0]

                xcorr = tmp_corr_df.v_prime1.dot( tmp_corr_df.v_prime2 ) / \
                        sqrt( 
                            tmp_corr_df.v_prime1.dot( tmp_corr_df.v_prime1 ) * \
                            tmp_corr_df.v_prime2.dot( tmp_corr_df.v_prime2 ) 
                        )

                if isnan(xcorr):
                    print tmp_corr_df[ isnan(tmp_corr_df.v_prime1) ].T

                delta_y.append( 
                    case_y_2.y.unique()[0] - case_y_1.y.unique()[0] 
                )

                f, gamma_loc = coherence(
                    tmp_corr_df.v_prime1,
                    tmp_corr_df.v_prime2, 
                    fs      = 10000,
                    nperseg = 2**7,
                )

                gamma.append( gamma_loc )

                delta_y_bl = y2_loc - y1_loc

                y1_real = case_y_1.y.unique()[0]
                y2_real = case_y_2.y.unique()[0]

                #xi = sqrt( ( x1_real - x2_real )**2 + \
                #          ( y1_real - y2_real )**2 )
                xi = y1_real - y2_real 

                corr_df = corr_df.append(
                    DataFrame( data = {
                        'xcorr':      xcorr,
                        'y2_bl':      y2_loc,
                        'y1_bl':      y1_loc,
                        'y2_real':    y2_real,
                        'y1_real':    y1_real,
                        'xi':         xi,
                        'case':       case,
                        'height_ix1': y1_ix,
                        'height_ix2': y2_ix,
                        'Um':         case_y_1.u.mean(),
                        'w':          case_y_1.w.mean(),
                        'urms':       case_y_1.u.std(),
                        'vrms':       case_y_1.v.std(),
                    } , index = [0] ), 
                    ignore_index = True
                )

                cnt += 1
                progress.update( cnt )

            gamma   = array(gamma)
            delta_y = array(delta_y)

            for ix in range(len(f)):
                int_gamma_per_freq.append( 
                    simps( gamma.T[ ix ], delta_y )
                )

            if plot_individual:
                X, Y = meshgrid( f, delta_y_bl )

                ax[0].contourf( X/1000., Y, gamma )

                ax[1].plot( array(f)/1000., int_gamma_per_freq , 'o')

                ax[0].set_xlabel('$f$ [kHz]')
                ax[0].set_ylabel(r'$\Delta y$')
                ax[1].set_xlabel('$f$ [kHz]')
                ax[1].set_ylim( 0, 4 )
                ax[1].set_ylabel( r'$\Lambda_{2|22} = \int_0^\delta\sqrt{' + \
                                 '\gamma^2\left('+\
                                 r'\omega,\xi_y\right)}\,\textrm{d}'+\
                                 r'\left( \xi_y \right)$' 
                                )

                fig.savefig( join( output_root, 'VCoh_{0}_y{1:.2f}.png'.format( 
                    case, y1_loc 
                ) ),
                    bbox_inches = 'tight')
                plt.close( fig )

        progress.finish()

        corr_df.to_pickle( 
            join( root, 'WallNormalCorrelation_Values_{0}.p'.format( case ) )
        )

def return_correct_schematic( case ):
    from os      import listdir
    from os.path import join
    from re      import findall

    root = '/home/carlos/Documents/PhD/Articles/Article_2/Figures/'

    available_schematics = [ 
        f for f in listdir( root ) \
        if f.startswith('measurement_locations') \
        if f.endswith('.png')
    ]

    z_loc  = findall( 'z[0-9][0-9]', case )[0]
    device = findall( 'S[0-9rTER]+', case )[0]
    for s in available_schematics:
        if device in s and s.endswith('_STE.png'):
            return join( root, s )
        elif z_loc in s and not device == "STE":
            return join( root, s )

    return 0


def plot_vertical_correlation_from_pickle( pickled_correlation, 
                                          plot_schematic = True , root = ''):
    import matplotlib.pyplot as plt
    from pandas           import read_pickle
    from os.path          import split, join
    from numpy            import meshgrid, arange
    from matplotlib       import rc
    #from matplotlib.cbook import get_sample_data

    rc('text',usetex=True)

    rc('font',
       family = 'serif',
       serif  = 'Linux Libertine',
       size   = 40,
       weight = 'normal'
      )

    fontsize = 25

    if not root:
        root = '/home/carlos/Documents/PhD/Articles/Article_3/' + \
                'Scripts/time_resolved/Results_v2'

    print "   Plotting correlation map for {0}".format( 
        split(pickled_correlation)[1]
    )
    corr_df = read_pickle( join( root, pickled_correlation ) )

    case = split( pickled_correlation )[1].replace('.p','').\
            replace("WallNormalCorrelation_Values_",'')

    corr_df = corr_df.sort_values( by = [ 'height_ix1', 'height_ix2' ] )

    # Reshape the correlation plot to always have a range of xi = 0 -> delta
    # in the x axis
    
    corr_pivot = corr_df.pivot( 'y1_bl', 'y2_bl', 'xcorr' )

    X   = corr_pivot.columns.values
    Y   = corr_pivot.index.values
    x,y = meshgrid(X, Y)
    Z   = corr_pivot.values

    fig,ax = plt.subplots(1,1, figsize = (7,7) )

    CS = ax.contourf(x, y, Z, levels = arange(0,1.1,0.1) , cmap = plt.cm.Blues)
    cbar = plt.colorbar(CS, shrink = 0.4)
    cbar.ax.set_ylabel(r'$R_{vv}$', fontsize = fontsize)
    cbar.ax.tick_params( labelsize = fontsize )

    #if plot_schematic:
    #    schematic = return_correct_schematic( pickled_correlation )

    #    im = plt.imread( get_sample_data( schematic ) )
    #    newax = fig.add_axes([0.5, 0.34, 0.22, 0.22], anchor = 'SW', 
    #                                 zorder=100)
    #    newax.imshow(im)
    #    newax.axis('off')

    ax.set_xlabel(r'$y/\delta$', fontsize = fontsize)
    ax.set_ylabel(r'$y/\delta$', fontsize = fontsize)
    ax.set_aspect( 'equal' )
    ax.set_ylim( 0.2, 1.2 )
    ax.set_xticks(  [0.25, 0.5, 0.75, 1, 2] )
    ax.set_xticklabels( [0.25, 0.5, 0.75, 1, 2], rotation = 30 )
                 
    ax.set_yticks( [0.25, 0.5, 0.75, 1] )
    ax.set_yticklabels( [0.25, 0.5, 0.75, 1] )
    ax.tick_params(axis='both', which='major', labelsize=fontsize)
    ax.set_xlim( left = 0.2 )

    fig.savefig( 
        join( output_root, '{0}.png'.format( case.replace("_Values","") ) ),
        bbox_inches = 'tight'
    )

    corr_pivot.to_pickle(
        join( root, 'WallNormalCorrelation_Pivot_{0}.p'.format( case ) )
    )

    plt.close( fig )

def plot_streamwise_correlation_from_pickle( pickled_correlation ,
                                           plot_schematic = True, 
                                            root = ''
                                           ):
    import matplotlib.pyplot as plt
    from pandas           import read_pickle
    from os.path          import split, join
    from numpy            import meshgrid, arange, round
    from matplotlib       import rc
    #from matplotlib.cbook import get_sample_data

    if not root:
        root = '/home/carlos/Documents/PhD/Articles/Article_3/' + \
                'Scripts/time_resolved/Results_v2'

    import matplotlib as mpl

    rc('text',usetex=True)
    rc('font',family='sans-serif', serif='sans-serif')

    mpl.rcParams['text.latex.preamble'] = [
        r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
        r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
        r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
        r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
    ]



    df   = read_pickle( join( root, pickled_correlation ) )
    case = split( pickled_correlation )[1].replace('.p','')

    df[ 'near_x_2h' ] = 0
    for stream_ix in df.stream_ix.unique():
        df.loc[ df.stream_ix == stream_ix, 'near_x_2h' ] =\
                round( df[ df.stream_ix == stream_ix ].x2_2h.mean() - \
                df[ df.stream_ix == stream_ix ].x1_2h.mean(), 1 )

    fig,ax = plt.subplots(1,1, figsize = (5,5) )
    try:
        corr_pivot = df.pivot( 'y1_bl', 'near_x_2h', 'xcorr' )

        X   = corr_pivot.columns.values
        Y   = corr_pivot.index.values
        x,y = meshgrid(X, Y)
        Z   = corr_pivot.values

        ax.contourf(x, y, Z, levels = arange(0,1.1,0.1), cmap = plt.cm.Blues)

    except ValueError:
        ax.tricontourf( 
            -df.xi, 
            df.y1, 
            df.xcorr, 
            levels = arange(0,1.1,0.1),
            cmap   = plt.cm.Blues
        )

    #if plot_schematic:
    #    schematic = return_correct_schematic( pickled_correlation )

    #    im = plt.imread( get_sample_data( schematic ) )
    #    newax = fig.add_axes([0.15, 0.13, 0.35, 0.35], anchor = 'SW', 
    #                                 zorder=100)
    #    newax.imshow(im)
    #    newax.axis('off')

    ax.set_xlabel(r'$\xi/2h$')
    ax.set_ylabel(r'$y/\delta$')
    #ax.set_ylim( 0.1, 1.2 )
    ax.set_ylim( 0., 20 )

    fig.savefig( join( output_root, '{0}.png'\
                .format( case.replace('_Values','') ) ),
                bbox_inches = 'tight'
               )

    plt.close( fig )

def get_streamwise_length_scale_and_ke( root = 0 ):
    from pandas          import read_pickle, DataFrame
    from os              import listdir
    from os.path         import join
    from scipy.integrate import simps
    from scipy.special   import gamma as gamma_func
    from math            import pi, sqrt
    import matplotlib.pyplot as plt
    from matplotlib.cbook import get_sample_data
    from matplotlib      import rc

    if not root:
        root = '/home/carlos/Documents/PhD/Articles/Article_3/' + \
                'Scripts/time_resolved/Results_v2'

    corr_pickles = [ f for f in listdir( root )\
                    if f.startswith( 'StreamwiseCorrelation_Values' )
                    and f.endswith( '.p' )
                    #and not "STE" in f
                   ]

    import matplotlib as mpl

    rc('text',usetex=True)
    rc('font',family='sans-serif', serif='sans-serif')

    mpl.rcParams['text.latex.preamble'] = [
        r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
        r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
        r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
        r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
    ]



    fontsize = 30

    fig, ax             = plt.subplots(1, 2, figsize=(10,4), sharey = True )
    fig_real, ax_real   = plt.subplots(1, 2, figsize=(10,4), sharey = True )
    fig_joint, ax_joint = plt.subplots(1, 2, figsize=(10,4), sharex = True )

    length_scale_df = DataFrame()

    for c in corr_pickles:
        corr_df = read_pickle( join( root, c ) )

        case_name = c\
                .replace('StreamwiseCorrelation_Values','')\
                .replace('.p','')
        color, marker, cmap = get_color_and_marker( 
            case_name
        )
        plot_config = {
            'marker'          : marker,
            'markeredgewidth' : markeredgewidth,
            'markerfacecolor' : markerfacecolor,
            'markeredgecolor' : color,
            'markersize'      : markersize,
            'mew'             : mew,
            'color'           : color,
        }

        corr_df = corr_df.sort_values( by = ['height_ix'] )\
                .reset_index( drop = True )

        corr_df = corr_df[ corr_df.y1_bl <= 1 ]

        for y in corr_df.height_ix.unique():

            corr_df_y = corr_df[ corr_df.height_ix == y ]
            corr_df_y = corr_df_y.sort_values( by = 'xi' )\
                    .reset_index( drop = True )
            corr_df_y = corr_df_y.dropna()

            if corr_df_y.empty:
                continue
            length_scale = simps(
                corr_df_y.xcorr,
                corr_df_y.xi
            )
            ke = sqrt( pi ) * gamma_func( 5. / 6. ) / \
                    ( length_scale / 1000. * gamma_func( 1. / 3. ) )

            length_scale_df = length_scale_df.append(
                DataFrame( data = {
                    'case': case_name,
                    'y':    y,
                    'ls':   length_scale,
                    'ke':   ke
                }, index = [0] ), ignore_index = True
            )

            ax[0].plot(
                length_scale,
                corr_df_y.y1_bl.unique()[0],
                **plot_config
            )

            ax[1].plot(
                ke,
                corr_df_y.y1_bl.unique()[0],
                **plot_config
            )

            ax_real[0].plot(
                length_scale,
                corr_df_y.y1.unique()[0],
                **plot_config
            )

            ax_real[1].plot(
                ke,
                corr_df_y.y1.unique()[0],
                **plot_config
            )

            ax_joint[0].plot(
                ke,
                corr_df_y.y1_bl.unique()[0],
                **plot_config
            )

            ax_joint[1].plot(
                ke,
                corr_df_y.y1.unique()[0],
                **plot_config
            )

    ax[0].set_ylim(0, 1)
    ax[0].set_xlim(0, 8)
    ax[0].set_xlabel( r'$\Lambda_{1|11}$ [m]$\times10^{-3}$', 
                     fontsize = fontsize )
    ax[0].set_ylabel( r'$y/\delta$', fontsize = fontsize )
    ax[1].set_xlabel( r'$k_e$ [1/m]', fontsize = fontsize )
    ax[1].set_xlim(100, 450)

    ax_real[0].set_ylim(0, 15)
    ax_real[0].set_xlim(1, 9)
    ax_real[0].set_xlabel( r'$\Lambda_{x|uu}$ [m]$\times10^{-3}$', 
                          fontsize = fontsize )
    ax_real[0].set_ylabel( r'$y$ [mm]', fontsize = fontsize )
    ax_real[1].set_xlabel( r'$k_e$ [1/m]', fontsize = fontsize )
    ax_real[1].set_xlim(50, 450)
    #ax_real[0].set_xticklabels( [100, '', 200, '', 300, '', 400, '' ] )
    ax_real[1].set_xticklabels( ['', 100, '', 200, '', 300, '', 400, '' ] )

    ax_joint[0].set_xlabel( r'$k_e$ [1/m]', fontsize = fontsize )
    ax_joint[1].set_xlabel( r'$k_e$ [1/m]', fontsize = fontsize )
    ax_joint[0].set_xlim(100, 450)
    ax_joint[1].set_xlim(100, 450)
    ax_joint[0].set_ylim(0, 1)
    ax_joint[1].set_ylim(0, 15)
    ax_joint[1].set_ylabel( r'$y/\delta$', fontsize = fontsize )
    ax_joint[1].set_ylabel( r'$y$ [mm]', fontsize = fontsize )

    schematic = '/home/carlos/Documents/PhD/Articles/Article_3/Scripts/'\
            +'time_resolved/LegendAll.png'

    im = plt.imread( get_sample_data( schematic  ) )
    newax = fig_real.add_axes([0.29, 0.12, 0.20, 0.25], 
                             anchor = 'SW', zorder=100)
    newax.imshow(im)
    newax.axis('off')

    ax_real[0].tick_params(labelsize = fontsize)
    ax_real[1].tick_params(labelsize = fontsize)
    ax_joint[0].tick_params(labelsize = fontsize)
    ax_joint[1].tick_params(labelsize = fontsize)
    ax[0].tick_params(labelsize = fontsize)
    ax[1].tick_params(labelsize = fontsize)

    fig.savefig(
        join( output_root, "StreamwiseLengthScales.png" ),
        bbox_inches = 'tight' 
    )
    fig_real.savefig(
        join( output_root, "StreamwiseLengthScales_real_yloc.png" ),
        bbox_inches = 'tight' 
    )
    fig_joint.savefig(
        join( output_root, "StreamwiseLengthScales_joint.png" ),
        bbox_inches = 'tight' 
    )

    length_scale_df.to_pickle( 
        join( root, 'StreamwiseLengthScales.p' )
    )

def plot_vertical_mean_values( root = 0 ):
    from pandas          import read_pickle
    from os              import listdir
    from os.path         import join
    from matplotlib      import rc
    #from matplotlib.cbook import get_sample_data
    import matplotlib.pyplot as plt
    from matplotlib.cbook                import get_sample_data

    import matplotlib as mpl
    if not root:
        root = '/home/carlos/Documents/PhD/Articles/Article_3/' + \
                'Scripts/time_resolved/Results_v2'

    corr_pickles = [ f for f in listdir( root )\
                    if f.startswith( 'WallNormal_MeanValues' )
                    and f.endswith( '.p' )\
                   ]
    #corr_pickles = [ f for f in corr_pickles\
    #                if 'STE' in f\
    #                or '_mfr' in f
    #               ]

    #fig_Um,   ax_Um   = plt.subplots( 1, 2, figsize = (10, 4) )
    fig_Um,   ax_Um   = plt.subplots( 1, 1, figsize = (7, 6) )
    fig_w,    ax_w    = plt.subplots( 1, 1, figsize = (7, 6) )
    fig_v,    ax_v    = plt.subplots( 1, 1, figsize = (7, 6) )
    fig_urms, ax_urms = plt.subplots( 1, 2, figsize = (12, 6) )
    fig_vrms, ax_vrms = plt.subplots( 1, 2, figsize = (12, 6) )

    fontsize = 30

    for c in corr_pickles:
        corr_df = read_pickle( join( root, c ) ).sort_values( 
            by = ['height_ix1'] )

        color, marker, cmap = get_color_and_marker( 
            c
        )
        #if c in corrections.keys():
        #    print c
        #    corr_df.y1_real = corr_df.y1_real + corrections[ c ]
        
        plot_config = {
            'marker'          : marker,
            'markeredgewidth' : markeredgewidth,
            'markerfacecolor' : markerfacecolor,
            'markeredgecolor' : color,
            'markersize'      : markersize,
            'mew'             : mew,
            'color'           : color,
        }


        for y in corr_df.height_ix1.unique():

            ue = corr_df.Um.max()

            ax_Um.plot(
                corr_df[ corr_df.height_ix1 == y ].Um.unique()[0] / ue,
                corr_df[ corr_df.height_ix1 == y ].y1_bl.unique()[0],
                **plot_config
            )
            #ax_Um[1].plot(
            #    corr_df[ corr_df.height_ix1 == y ].Um.unique()[0] / ue,
            #    corr_df[ corr_df.height_ix1 == y ].y1_real.unique()[0],
            #    **plot_config
            #)
            ax_v.plot(
                corr_df[ corr_df.height_ix1 == y ].v.unique()[0] / ue * 10.,
                corr_df[ corr_df.height_ix1 == y ].y1_bl.unique()[0],
                **plot_config
            )
            #ax_v[1].plot(
            #    corr_df[ corr_df.height_ix1 == y ].v.unique()[0] / ue * 10.,
            #    corr_df[ corr_df.height_ix1 == y ].y1_real.unique()[0],
            #    **plot_config
            #)
            ax_urms[0].plot(
                corr_df[ corr_df.height_ix1 == y ].urms.unique()[0] / ue * 10.,
                corr_df[ corr_df.height_ix1 == y ].y1_bl.unique()[0],
                **plot_config
            )
            ax_urms[1].plot(
                corr_df[ corr_df.height_ix1 == y ].urms.unique()[0] / ue * 10.,
                corr_df[ corr_df.height_ix1 == y ].y1_real.unique()[0],
                **plot_config
            )
            ax_vrms[0].plot(
                corr_df[ corr_df.height_ix1 == y ].vrms.unique()[0] / ue * 100.,
                corr_df[ corr_df.height_ix1 == y ].y1_bl.unique()[0],
                **plot_config
            )
            ax_vrms[1].plot(
                corr_df[ corr_df.height_ix1 == y ].vrms.unique()[0] / ue * 100.,
                corr_df[ corr_df.height_ix1 == y ].y1_real.unique()[0],
                **plot_config
            )
            ax_w.plot(
                corr_df[ corr_df.height_ix1 == y ].w.unique()[0] / ue,
                corr_df[ corr_df.height_ix1 == y ].y1_bl.unique()[0],
                **plot_config
            )
            #ax_w[1].plot(
            #    corr_df[ corr_df.height_ix1 == y ].w.unique()[0] / ue,
            #    corr_df[ corr_df.height_ix1 == y ].y1_real.unique()[0],
            #    **plot_config
            #)

    #ax_Um[0].set_xlim(0, 1)
    ax_Um.set_ylim(0, 1)
    ax_Um.set_xlabel( r'$\overline{u} / u_e$', fontsize = fontsize )
    ax_Um.set_ylabel( r'$y/\delta$', fontsize = fontsize )
    ax_Um.set_xlim(0.2, 1)
    ax_Um.set_xticklabels( [ 0.2, '', 0.4, '', 0.6, '', 0.8, '', 1.0 ] )
    #ax_Um[1].set_xticklabels( [ 0.2, '', 0.4, '', 0.6, '', 0.8, '', 1.0 ] )
    #ax_Um[1].set_xlim(0.2, 1)
    #ax_Um[1].set_ylim(0, 15)
    #ax_Um[1].set_xlabel( r'$\overline{u} / u_e$', fontsize = fontsize )
    #ax_Um[1].set_ylabel( r'$y$ [mm]', fontsize = fontsize )

    ax_urms[0].set_ylim(0, 1)
    ax_urms[0].set_xlabel( r'$u_{\mathrm{rms}} / u_e\times 10^{-1}$', 
                          fontsize = fontsize )
    ax_urms[0].set_ylabel( r'$y/\delta$', fontsize = fontsize )

    ax_urms[1].set_ylim(0, 15)
    ax_urms[1].set_ylabel( r'$y$ [mm]', fontsize = fontsize )
    ax_urms[1].set_xlabel( r'$u_{\mathrm{rms}} / u_e\times 10^{-1}$', 
                          fontsize = fontsize )

    ax_vrms[0].set_ylim(0, 1)
    ax_vrms[0].set_xlabel( r'$v_{\mathrm{rms}} / u_e\times 10^{-2}$', 
                          fontsize = fontsize )
    ax_vrms[0].set_ylabel( r'$y/\delta$', fontsize = fontsize )
    ax_vrms[0].set_xticklabels( [ 0, '', 2, '', 4, '', 6, '', 8 ] )
    ax_vrms[1].set_xticklabels( [ 0, '', 2, '', 4, '', 6, '', 8 ] )
    
    ax_vrms[1].set_ylim(0, 15)
    ax_vrms[1].set_xlabel( r'$v_{\mathrm{rms}} / u_e\times 10^{-2}$', 
                          fontsize = fontsize )
    ax_vrms[1].set_ylabel( r'$y$ [mm]', fontsize = fontsize )

    ax_w.set_ylim(0, 1)
    ax_w.set_xlabel( r'$\overline{w} / u_e$', fontsize = fontsize )
    ax_w.set_ylabel( r'$y/\delta$', fontsize = fontsize )
    ax_w.set_xticklabels( [ -0.1, '', 0, '', 0.1, '', 0.2 ] )

    #ax_w[1].set_ylim(0, 15)
    #ax_w[1].set_xlabel( r'$\overline{w} / u_e$', fontsize = fontsize )
    #ax_w[1].set_ylabel( r'$y$ [mm]', fontsize = fontsize )
    #ax_w[1].set_xticklabels( [ -0.1, '', 0, '', 0.1, '', 0.2 ] )

    ax_v.set_ylim(0, 1)
    ax_v.set_xlim(-2.25, 0.75)
    ax_v.set_xticklabels( [ '', -2.0, '', -1.0, '', 0, 0.5, 1.0] )
    ax_v.set_xlabel( r'$\overline{v} / u_e\times10^{-1}$', fontsize = fontsize )
    ax_v.set_ylabel( r'$y/\delta$', fontsize = fontsize )

    #ax_v[1].set_ylim(0, 15)
    #ax_v[1].set_xlim(-2.25, 0.75)
    #ax_v[1].set_xticklabels( [ '', -2.0, '', -1.0, '', 0, 0.5, 1.0] )
    #ax_v[1].set_xlabel( r'$\overline{v} / u_e\times10^{-1}$', fontsize = fontsize )
    #ax_v[1].set_ylabel( r'$y$ [mm]', fontsize = fontsize )
    #ax_v[0].set_xticklabels( [ -0.1, '', 0, '', 0.1, '', 0.2 ] )
    #ax_v[1].set_xticklabels( [ -0.1, '', 0, '', 0.1, '', 0.2 ] )



    schematic = '/home/carlos/Documents/PhD/Articles/Article_3/Scripts/'\
            +'time_resolved/LegendAll.png'

    im = plt.imread( get_sample_data( schematic  ) )
    newax = fig_w.add_axes([0.50, 0.60, 0.40, 0.45], 
                             anchor = 'SW', zorder=100)
    newax.imshow(im)
    newax.axis('off')

    im = plt.imread( get_sample_data( schematic  ) )
    newax = fig_vrms.add_axes([0.70, 0.60, 0.20, 0.25], 
                             anchor = 'SW', zorder=100)
    newax.imshow(im)
    newax.axis('off')

    newax = fig_Um.add_axes([0.15, 0.60, 0.40, 0.45], 
                             anchor = 'SW', zorder=100)
    newax.imshow(im)
    newax.axis('off')

    newax = fig_v.add_axes([0.14, 0.11, 0.40, 0.45], 
                             anchor = 'SW', zorder=100)
    newax.imshow(im)
    newax.axis('off')

    ax_Um.tick_params(labelsize   = fontsize)
    #ax_Um[1].tick_params(labelsize   = fontsize)
    ax_w.tick_params(labelsize    = fontsize)
    #ax_w[1].tick_params(labelsize    = fontsize)
    ax_v.tick_params(labelsize    = fontsize)
    #ax_v[1].tick_params(labelsize    = fontsize)
    ax_urms[0].tick_params(labelsize = fontsize)
    ax_urms[1].tick_params(labelsize = fontsize)
    ax_vrms[0].tick_params(labelsize = fontsize)
    ax_vrms[1].tick_params(labelsize = fontsize)

    #fig_Um.subplots_adjust(wspace   = 0.3)
    fig_w.subplots_adjust(wspace    = 0.3)
    fig_v.subplots_adjust(wspace    = 0.3)
    fig_urms.subplots_adjust(wspace = 0.3)
    fig_vrms.subplots_adjust(wspace = 0.3)


    rc('text',usetex=True)
    rc('font',family='sans-serif', serif='sans-serif')

    mpl.rcParams['text.latex.preamble'] = [
        r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
        r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
        r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
        r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
    ]

    fig_Um.savefig( join( output_root, "Um.png" ), bbox_inches = 'tight' )
    fig_w.savefig( join( output_root, "w.png" ), bbox_inches = 'tight' )
    fig_v.savefig( join( output_root, "v.png" ), bbox_inches = 'tight' )
    fig_urms.savefig( join( output_root, "urms.png" ), bbox_inches = 'tight' )
    fig_vrms.savefig( join( output_root, "vrms.png" ), bbox_inches = 'tight' )

def get_vertical_length_scale( root = 0 ):
    from pandas          import read_pickle, DataFrame, concat
    from os              import listdir
    from os.path         import join
    from scipy.integrate import simps
    from matplotlib      import rc
    #from matplotlib.cbook import get_sample_data
    import matplotlib.pyplot as plt

    fontsize = 40

    if not root:
        root = '/home/carlos/Documents/PhD/Articles/Article_3/' + \
                'Scripts/time_resolved/Results_v2'

    corr_pickles = [ f for f in listdir( root )\
                    if f.startswith( 'WallNormalCorrelation_Values' )
                    and f.endswith( '.p' )\
                   ]
    corr_pickles = [ f for f in corr_pickles\
                    if 'STE' in f\
                    or '_mfr' in f
                   ]

    fig,    ax        = plt.subplots( 1, 2, figsize = (10,5) )

    length_scale_df = DataFrame()

    for c in corr_pickles:
        corr_df = read_pickle( join( root, c ) ).sort_values( 
            by = ['height_ix1', 'height_ix2'] )

        color, marker, cmap = get_color_and_marker( 
            c
        )

        #if c in corrections.keys():
        #    print c
        #    corr_df.y1_real = corr_df.y1_real + corrections[ c ]
        
        plot_config = {
            'marker'          : marker,
            'markeredgewidth' : markeredgewidth,
            'markerfacecolor' : markerfacecolor,
            'markeredgecolor' : color,
            'markersize'      : markersize,
            'mew'             : mew,
            'color'           : color,
        }


        ls         = []
        y_bl_loc   = []
        y_real_loc = []
        for y in corr_df.height_ix1.unique():

            corr_df_y = corr_df[ corr_df.height_ix1 == y ].dropna()
            corr_df_y.xi = -corr_df_y.xi

            corr_df_y = corr_df_y.sort_values( by = 'y2_bl', ascending = True )
            
            bl_thickness = ( corr_df_y.y1_real / corr_df_y.y1_bl ).mean()

            corr_df_y = corr_df_y[ corr_df_y.xi >= -0.01 ]
            
            # Ignore y locations for which the xi does not extend to the length
            # of the boundary layer thickness (less a bit of buffer)
            if corr_df_y.xi.max() < bl_thickness * 0.8 or corr_df_y.empty:
                continue
            #print corr_df_y[ ['xi','xcorr','y2_bl', 'y1_bl'] ]

            corr_df_y = corr_df_y[ corr_df_y.xi <= bl_thickness ]

            length_scale = simps(
                corr_df_y.xcorr,
                corr_df_y.xi
            )

            ls.append( length_scale )
            y_real_loc.append( 
                corr_df[ corr_df.height_ix1 == y ].y1_real.unique()[0] 
            )
            y_bl_loc.append( 
                corr_df[ corr_df.height_ix1 == y ].y1_bl.unique()[0] 
            )

            ax[0].plot(
                length_scale,
                corr_df_y.y1_bl.unique()[0],
                **plot_config
            )
            ax[1].plot(
                length_scale,
                corr_df_y.y1_real.unique()[0],
                **plot_config
            )

        length_scale_df = concat([ length_scale_df, 
            DataFrame( data = {
                ( corr_df_y.case.unique()[0], 'ls' ) : ls,
                ( corr_df_y.case.unique()[0], 'y' ) : \
                y_real_loc,
                ( corr_df_y.case.unique()[0], 'ybl' ) : \
                y_bl_loc,
            }, index = range( len( y_bl_loc ) ) ) ], axis = 1
        )

    ax[0].set_ylim(0, 1)
    ax[0].set_ylabel( r'$y/\delta$' , fontsize = fontsize)
    ax[1].set_ylabel( r'$y$ [mm]' , fontsize = fontsize)
    for axi in ax:
        axi.set_xlim(0, 5)
        axi.set_xticks( [1, 2, 3, 4, 5] )
        axi.set_xlabel( r'$\Lambda_{y|vv}$ [m]$\times10^{-3}$' , 
                      fontsize = fontsize)
        axi.tick_params(axis='both', which='major', labelsize=fontsize)

    #schematic = '/home/carlos/Documents/PhD/Articles/Article_2/'+\
            #'Figures/measurement_locations_TE_m2.png'

    #im = plt.imread( get_sample_data( schematic  ) )
    #newax = fig.add_axes([0.60, 0.5, 0.35, 0.35], anchor = 'SW', 
    #                             zorder=100)
    #newax.imshow(im)
    #newax.axis('off')

    import matplotlib as mpl

    rc('text',usetex=True)
    rc('font',family='sans-serif', serif='sans-serif')

    mpl.rcParams['text.latex.preamble'] = [
        r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
        r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
        r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
        r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
    ]



    fig.subplots_adjust(wspace = 0.3)

    fig.savefig( 
        join( output_root, "VerticalLengthScales.png"), 
        bbox_inches = 'tight' 
    )

    length_scale_df.to_pickle( join( root, 'VerticalLengthScales.p' ) )

def print_vertical_length_scale_bl_integration( root = 0 ):
    from pandas          import read_pickle
    from scipy.integrate import simps
    from os.path         import join

    if not root:
        root = '/home/carlos/Documents/PhD/Articles/Article_3//' + \
                'Scripts/time_resolved/Results_v2'

    length_scale_df = read_pickle( join( root, 'VerticalLengthScales.p' ) )

    cases = length_scale_df.columns.levels[0]

    for c in cases:
        df_case = length_scale_df[ c ].dropna().sort_values( by = 'y' )

        ls = simps(
            df_case[ df_case.ybl < 1 ].ls.values,
            df_case[ df_case.ybl < 1 ].y.values,
        )
        print c, round( ls, 1 )


def get_color_and_marker(case_name):
    import matplotlib.pyplot as plt
    import seaborn.apionly as sns

    fig, ax = plt.subplots( 1, 1 , figsize = ( 6, 4 ) )

    strings = [
        'Sr20R21_a0_p0',
        #'Sr20R21_a12_p0',
        'Sr20R21_a12_p6',
        #'Sr20R21_a-12_p0',
        'Sr20R21_a-12_p6',
        'STE_a12_p0',
        'STE_a-12_p0',
    ]

    labels = {
        'Sr20R21_a0_p0':   r'$\alpha = 0^\circ,\,\varphi = 0^\circ $',
        #'Sr20R21_a12_p0':  r'$\alpha = 12^\circ,\,\varphi = 0^\circ $, SS',
        'Sr20R21_a12_p6':  r'$\alpha = 12^\circ,\,\varphi = 6^\circ $, SS',
        #'Sr20R21_a-12_p0': r'$\alpha = 12^\circ,\,\varphi = 0^\circ $, PS',
        'Sr20R21_a-12_p6': r'$\alpha = 12^\circ,\,\varphi = 6^\circ $, PS',
        'STE_a12_p0':      r'Straight, $\alpha = 12^\circ$, SS',
        'STE_a-12_p0':     r'Straight, $\alpha = 12^\circ$, PS',
    }

    cmaps = [
        'Blues',
        #'Greens',
        'Oranges',
        #'Oranges',
        'Oranges',
        'Oranges',
        'Oranges',
    ]

    markers = [
        '.',
        #'1',
        '2',
        #'3',
        '4',
        '+',
        'x',
    ]

    current_palette = sns.color_palette( 
        'deep', 
        n_colors = len( strings ) 
    )
                                       

    color = 0

    for ix in range( len( strings ) ):
        if strings[ ix ] in case_name:
            color  = current_palette[ ix ]
            marker = markers [ ix ]
            cmap   = cmaps[ ix ]

        ax.scatter(
            [ ix ], [ 1 ],
            c      = current_palette[ ix ],
            marker = markers[ ix ],
            s      = 100,
            label  = labels[strings[ ix ]]
        )

    ax.legend( loc = 'best' )
    fig.savefig( 'Legend.png' )
    if not color:
        print case_name

    return color, marker, cmap


                        
def do_the_reynolds_stress_quadrant_analysis(cases_df,y_delta, plot_name = ''):

    from matplotlib import rc
    from seaborn import kdeplot
    import matplotlib.pyplot as plt
    from numpy import linspace, mean

    rc('text',usetex=True)
    rc('font',weight='normal')

    #sns.set_context('paper')
    #sns.set(font='serif',font_scale=5.0,style='whitegrid')
    rc('font',family='serif', serif='Linux Libertine')

    cases = sorted(cases_df.file.unique(),reverse=True)

    fig_uv,axes_uv = plt.subplots(
        1,len(cases), 
        figsize = (figsize[0]*len(cases), figsize[1]), 
        sharex=True,
        sharey=True, 
    )

    if not len(cases)>2:
        return 0
        
    for case_file, case_i, ax_uv in zip(cases, range(len(cases)), 
                                     axes_uv):

        case = cases_df[cases_df.file == case_file]\
                .sort_values( by = ['time_step'] )

        bl_data = return_bl_parameters( case , [ 20 ] )

        color, marker, cmap = get_color_and_marker( 
            case_file
        )

        if not color and not marker:
            break

        case["uprime"] = ( case.u - case.u.mean() ) / bl_data.Ue.values[0]
        case["vprime"] = ( case.v - case.v.mean() ) / bl_data.Ue.values[0]
        case["wprime"] = ( case.w - case.w.mean() ) / bl_data.Ue.values[0]
        
        
        kde_uv = kdeplot(case.uprime, case.vprime,
                             cmap         = cmap,
                             ax           = ax_uv,
                             shade        = True,
                             shade_lowers = False,
                             gridsize     = 80,
                             zorder = 120, 
                             alpha = 0.8,
                   )
        
        ax_uv.plot(
            case['uprime'].values[::100],
            case['vprime'].values[::100],
            ls              = '',
            marker          = marker,
            markeredgewidth = markeredgewidth,
            markerfacecolor = markerfacecolor,
            markeredgecolor = 'k',
            markersize      = markersize*1.8,
            mew             = mew,
            alpha           = 0.4,
            zorder = 100
        )

        t = linspace( -0.3, 0.3, 100 )

        up = 6 * -mean(case['uprime']*case['vprime'])/t
        ax_uv.plot( t, up , ls = '-.', c = 'k', lw = 3 , zorder = 500)
        ax_uv.plot( t, -up , ls = '-.', c = 'k', lw = 3 , zorder = 500)

        ax_uv.fill_between( t, 
                           -up,
                           up,
                           facecolor='w', 
                           alpha=1.0,
                           zorder = 80
                          )


        ax_uv.plot(
            case['uprime'].values[::10],
            case['vprime'].values[::10],
            ls              = '',
            marker          = marker,
            markeredgewidth = markeredgewidth,
            markerfacecolor = markerfacecolor,
            markeredgecolor = 'k',
            markersize      = markersize*1.8,
            mew             = mew,
            alpha           = 0.3,
            zorder = 1
        )

        ax_uv.set_aspect('equal')
        kde_uv.collections[0].set_alpha(0)

        ax_uv.set_xlabel('')
        ax_uv.set_ylabel('')

        ax_uv.set_xlim( -0.25 , 0.25 )
        ax_uv.set_ylim( -0.25 , 0.25 )
        ax_uv.set_xticks([-0.2,0,0.2])
        ax_uv.set_yticks([-0.2,0,0.2])

        ax_uv.axhline( 0, ls = '--', lw=3 , c = 'k', zorder = 150)
        ax_uv.axvline( 0, ls = '--', lw=3 , c = 'k', zorder = 150)

        if y_delta == 0.1:
            ax_uv.set_xlabel(r"$u'/u_e$")
            ax_uv.grid(False)

    axes_uv[0].set_ylabel(r"$v'/u_e$")

    t = axes_uv[0].text( -0.22, 0.15, 
                        r'$y/\delta_{{99}} = {0}$'.format(y_delta),
                       zorder = 300)
    t.set_bbox(dict(color='white', alpha=0.9, edgecolor='white', zorder = 300))

    t = axes_uv[-1].text( -0.22, 0.15, 
                        r'II'.format(y_delta),
                       zorder = 300)
    t.set_bbox(dict(color='white', alpha=0.9, edgecolor='white', zorder = 300))

    t = axes_uv[-1].text( 0.22, -0.18, 
                        r'IV'.format(y_delta), ha='right',
                       zorder = 300)
    t.set_bbox(dict(color='white', alpha=0.9, edgecolor='white', zorder = 300))

    plot_name = 'Results/ReynoldsQuadrant_{0}_ydelta{1:.2f}.png'\
        .format( plot_name,y_delta )

    fig_uv.savefig( 
        plot_name.replace('.','_').replace('_png','_uv.png'),
        bbox_inches = 'tight'
    )
    plt.cla()


def get_kolmogorov_law_curve( x_lim = (0.5,1.5) ):
    from numpy import linspace, log10

    slope_x = linspace( x_lim[0], x_lim[1] ,100 )

    slope_y = 10*log10(slope_x**(-5/3.))

    return [slope_x,slope_y]

def plot_mean_and_std( df , plot_name = '' ):
    import matplotlib.pyplot as plt
    from matplotlib import rc
    #import seaborn as sns
    from numpy import array

    rc('text',usetex=True)
    rc('font',weight='normal')

    #sns.set_context('paper')
    #sns.set(font='serif',font_scale=4.0,style='whitegrid')
    rc('font',family='serif', serif='Linux Libertine')

    for var in ['u','v','w']:
        fig_std,   ax_std   = plt.subplots(1, 1, figsize = figsize)
        fig_mean,  ax_mean  = plt.subplots(1, 1, figsize = figsize)
        for case in df.case_name.unique():
            case_df = df[df.case_name == case]
            bl_data = return_bl_parameters( case_df , [ 20 ] )
            y_locs = sorted( case_df.near_y_delta.unique() )
            Um  = []
            std = []
            real_y_locs = []
            for y_loc in y_locs:
                y_df = case_df[ 
                    (case_df.near_y_delta == y_loc) &\
                    (case_df.case_name == case) 
                ]
                Um.append(
                    y_df[var].mean()
                )
                std.append(
                    y_df[var].std()
                )
                real_y_locs.append(
                    case_df[ case_df.near_y_delta == y_loc ]\
                    .y.unique()[0]
                )

            color, marker, cmap = get_color_and_marker( case )

            plot_config = {
                'marker'          : marker,
                'markeredgewidth' : markeredgewidth,
                'markerfacecolor' : markerfacecolor,
                'markeredgecolor' : color,
                'markersize'      : markersize,
                'mew'             : mew,
                'color'           : color,
            }

            ax_std.plot(
                array(std) / bl_data.Ue.values[0], 
                real_y_locs / bl_data.delta_99.values[0],
                ls='',
                **plot_config
            )

            ax_mean.plot(
                array(Um) / bl_data.Ue.values[0], 
                real_y_locs / bl_data.delta_99.values[0],
                ls='',
                **plot_config
            )


        fig_std.savefig(
            "Results/std_{0}_{1}.png".format(
                plot_name.replace('.','_').replace('_png','.png'),
                var,
            ), bbox_inches = 'tight'
        )
        fig_mean.savefig(
            "Results/mean_{0}_{1}.png".format(
                plot_name.replace('.','_').replace('_png','.png'),
                var,
            ), bbox_inches = 'tight'
        )
        
def do_the_time_resolved_analysis():
    from os.path import join

    root = '/home/carlos/Documents/PhD/Articles/Article_3/' + \
            'Scripts/time_resolved/Results_v2'

    source_root = '/home/carlos/Documents/PhD/Articles/Article_3/' + \
            'Scripts/time_resolved/ReservedData_OrigProc'

    def do_the_vertical_coherence_analysis( 
        hdfs , 
        plot_individual = False,
        overwrite = False
    ):

        #get_vertical_correlation( 
        #    [ join( source_root, h ) for h in hdfs ], 
        #    root            = root,
        #    plot_individual = plot_individual ,
        #    overwrite       = overwrite
        #)
        get_vertical_mean_values( 
            [ join( source_root, h.replace( '_mfr', '' ) ) for h in hdfs ], 
            root            = root,
            overwrite       = overwrite
        )

        #for hdf in hdfs:
        #    plot_vertical_correlation_from_pickle( 
        #        join( 
        #            root, 
        #            'WallNormalCorrelation_Values_' + hdf.replace( 
        #                '.hdf5', '.p' 
        #            )
        #        ),
        #        root = root
        #    )

        get_vertical_length_scale( root = root )
        plot_vertical_mean_values( root = root )
        #print_vertical_length_scale_bl_integration( )

    def do_the_streamwise_analysis( hdfs , overwrite = False ):

        #for hdf in hdfs:
        #    #get_streamwise_coherence_and_correlation( 
        #    #    join( source_root, hdf ), 
        #    #    overwrite = overwrite 
        #    #)

        #    plot_streamwise_correlation_from_pickle( 
        #        'StreamwiseCorrelation_Values_' + hdf.replace( '.hdf5', '.p' )
        #    )

        #    do_the_streamwise_coherence_analysis(
        #        'StreamwiseCoherence_Values_' + hdf.replace( '.hdf5', '.p' ),
        #        overwrite = overwrite
        #    )

        #get_streamwise_length_scale_and_ke()

       # plot_pickled_Uc( 
       #     [ 'Uc_data_Values_' + f.replace('.hdf5','.p') for f in hdfs ]
       # )
        
        plot_frequency_spectra( 
            [ join( source_root, h.replace('_mfr','' ) ) for h in hdfs ], 
            var = 'u'
        )
        #plot_wavenumber_spectra( 
        #    [ join( source_root, h.replace('_mfr','' ) ) for h in hdfs ], 
        #    var = 'u'
        #)
        #plot_phi()

    hdf_list_to_process = [
        #'Sr20R21_a0_p0_U20_z05_tr.hdf5', 
        #'Sr20R21_a12_p0_U20_z05_tr_mfr.hdf5', 
        #'Sr20R21_a-12_p0_U20_z05_tr_mfr.hdf5', 
        'Sr20R21_a12_p6_U20_z05_tr_mfr.hdf5', 
        'Sr20R21_a-12_p6_U20_z05_tr_mfr.hdf5', 
        'STE_a12_p0_U20_z00_tr.hdf5', 
        'STE_a-12_p0_U20_z00_tr.hdf5', 
    ]

    #do_the_vertical_coherence_analysis(
    #    hdf_list_to_process,
    #    plot_individual = False,
    #    overwrite       = False
    #    )

    do_the_streamwise_analysis( 
        hdf_list_to_process,
        overwrite = False,
    )
    

# Constants ####################################################################

nperseg         = 2**6
fs              = 10000

St_min          = 0.225
St_max          = 2.4

markeredgewidth = 3
markerfacecolor = 'none'
markersize      = 17
mew             = 4 # Marker edge width

figsize         = (8,7)

serration_angle = 76.

outlier_thresh  = 3.5

# ##############################################################################

