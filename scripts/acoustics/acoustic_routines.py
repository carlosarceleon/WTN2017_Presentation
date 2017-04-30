def read_mat_file_to_df( mat_file ):
    from scipy.io import loadmat
    from pandas import DataFrame

    data = loadmat(mat_file)
    
    print data['psd'].T[0]
    df = DataFrame(
        data = {
            "f":   data['f'].T[0],
            'spl': data['psd'].T[0]
            }, 
        index=range(len(data['f'].T[0]))
    )

    return df

def return_comparisons( root = 'new_data' ):
    from os.path import join
    from pandas import DataFrame, concat

    prefix = 'psd_'
    suffix = '_new'

    alphas = [ 0,  6,  12 ]
    phis   = [ 0,  6      ]
    us     = [ 30, 35, 40 ]
    
    # Compare same alphas, for different Us, phis are separate #################
    cnt = 0
    for u in us:
        for alpha in alphas:
            # straight edge case ###############################################
            ste_case_name = prefix + 'STE_a{0:02d}_U{1}'.format( alpha, u ) \
                    + suffix + '.mat'

            ste_case = read_mat_file_to_df( join( root, ste_case_name ) )
            
            for phi in phis:

                sr_case_name_base = 'Sr20R21_a{0:02d}_p{1}_U{2}'.\
                        format( alpha, phi, u )
                sr_case_name = prefix \
                        +  sr_case_name_base \
                        + suffix + '.mat'

                try:
                    sr_case = read_mat_file_to_df( join( root, sr_case_name ) )
                except:
                    continue

                delta = ( ste_case.spl - sr_case.spl ).values

                result = DataFrame( data = { sr_case_name_base : delta },
                                              index = ste_case.f
                                             )
                if not cnt:
                    comparisons = result.copy()

                else:
                    comparisons = concat(
                        [comparisons, result], axis = 1
                    )

                cnt += 1
    return comparisons

def plot_comparisons_by_alpha_phi( root = 'new_data', output_folder = 'results' ):
    import matplotlib.pyplot as plt
    import seaborn.apionly as sns
    from os.path           import join
    from matplotlib.ticker import ScalarFormatter
    from matplotlib        import rc

    import matplotlib as mpl

    source_map_root = '/home/carlos/Documents/PhD/Articles/Article_3/'\
            + 'Scripts/new_acoustics/'
    root = join(source_map_root,root)

    sns.set(font_scale=3.6,style='ticks',
            rc={"axes.axisbelow": False,'image.cmap': 'OrRd'})

    rc('text',usetex=True)
    rc('font',family='sans-serif', serif='sans-serif')

    mpl.rcParams['text.latex.preamble'] = [
        r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
        r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
        r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
        r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
    ]
    comparisons = return_comparisons( root = root )

    alphas = [ 0,  6,  12 ]
    phis   = [ 0,  6      ]
    us     = [ 30, 35, 40 ]

    for phi in phis:
        fig, axes = plt.subplots( 1, len( alphas ), figsize = ( 20, 5 ), 
                                sharex = True, sharey = True )
        for alpha, ax in zip( alphas, axes ):

            for u in us:

                case_name = 'Sr20R21_a{0:02d}_p{1}_U{2}'.\
                        format( alpha, phi, u )

                if not case_name in comparisons.columns:
                    continue
                ax.plot(
                    comparisons[ case_name ].index / 1000.,
                    comparisons[ case_name ].values,
                    label = "{0} m/s".format( u ),
                    lw = 2
                )

            ax.set_title( r'$\alpha = {0}^\circ$'.format( alpha ) )
            ax.set_ylim( -4, 7 )
            ax.set_xlim( 1, 5 )
            ax.set_xlabel(r'$f$ [kHz]')
            xticks = [ 1, 2, 3, 4, 5 ]
            ax.set_xticks( xticks )
            ax.set_xticklabels( map( str, xticks ) )
            ax.set_xscale( 'log' )
            ax.get_xaxis().set_major_formatter(ScalarFormatter())
            plt.setp( ax.xaxis.get_majorticklabels(), rotation=0 )
            #ax.set_grid(True,which='both',ls="-", color='0.65')
            ax.grid(True,which='both')
            ax.axhline( y = 0 , c = 'k', lw = 2 )
            if ax == axes[0]:
                ax.legend( loc = 'best' )

        axes[0].set_ylabel("$\\Delta\\textrm{SPL}$ [dB]")
        fig.savefig( join( 
            output_folder, 
            'Relative_p{0}'.format( phi )
        ), bbox_inches = 'tight' )


def return_source_map_case_descriptions( cases ):
    from re import findall
    from pandas import DataFrame
    from os.path import split
    import scipy.io
    from numpy import meshgrid

    cases_df = DataFrame()
    for c in cases:

        velocity = findall( 'U[0-9][0-9]', split( c )[1] )[0]
        #device   = findall( '[TESrR0-9]+', split( c )[1] )[0]
        device   = split(c)[1].\
                replace( '-'+velocity, '' ).\
                replace( ' - 4000 Hz.mat', '' )
        
        mat = scipy.io.loadmat( c )

        plane      = mat.get('scan_plane')
        source_map = mat.get('amp_log')

        x = plane[0][0][0][0]
        y = plane[0][0][1][0]

        X, Y = meshgrid( x, y )

        case_df = DataFrame(
            data = {
                'x' : X.ravel(),
                'y' : Y.ravel(),
                'A' : source_map.T.ravel()
            }
        )

        
        case_df['device']   = device
        case_df['velocity'] = velocity
        case_df['file']     = c

        cases_df = cases_df.append( 
            case_df,
            ignore_index = True
        )

    return cases_df

def do_source_plots():
    import matplotlib.pyplot  as plt
    import seaborn.apionly            as sns
    import matplotlib.patches as patches

    from matplotlib import rc
    from numpy      import arange,meshgrid
    from os         import listdir
    from os.path    import join
    from math       import floor, ceil
    from collections import OrderedDict

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


    source_map_root = '/home/carlos/Documents/PhD/Articles/Article_3/'\
            + 'Scripts/new_acoustics/source_data'

    source_maps = [join( source_map_root, f ) for f in \
                   listdir( source_map_root ) if f.endswith('.mat') ]

    cases_df = return_source_map_case_descriptions( source_maps )

    # Get the total levels #####################################################
    maxA = ceil( cases_df.A.max() )+1
    minA = floor( cases_df.A.min() )

    minA = maxA - 12

    levels = arange( minA , maxA, 2 )

    fig, axes = plt.subplots( 
        1, len( cases_df.device.unique() ), 
        figsize = (18*0.7, 6.8*0.7),
        sharex  = True,
        sharey  = True,
    )


    cases_to_run = OrderedDict()
    cases_to_run['STE-a00'] = 'Straight edge'
    cases_to_run['Sr20R21-a00-p0'] = r'Serrated, $\varphi=0^\circ$'
    cases_to_run['Sr20R21-a00-p6'] = r'Serrated, $\varphi=6^\circ$'

    for dev, ax in zip( cases_to_run.keys(), axes ):

        dev_case_df = cases_df[ cases_df.device == dev ]

        x = dev_case_df.x.unique()
        y = dev_case_df.y.unique()

        X, Y = meshgrid( x, y )

        A = dev_case_df.A.values.reshape( X.shape )

        contf = ax.contourf( X, Y, A , 
                            vmin   = minA,
                            vmax   = maxA ,
                            levels = levels
                           )

        CS4 = ax.contour(X, Y, A, levels,
              colors=('k',),
              linewidths=(1,),
              )

        ax.clabel(CS4, fmt='%2.0f', colors='k', fontsize = 15)

        ax.text(
            0, 1.02,
            cases_to_run[ dev ],
            horizontalalignment = 'left',
            verticalalignment   = 'bottom',
            transform           = ax.transAxes,
            fontsize            = 20
        )

        ax.add_patch(
            patches.Rectangle(
                (-0.2, -0.2),
                0.4,
                0.2,
                fill=False,# remove background
                linewidth = 1
            )
        )

        serration_coords_x = arange(-0.2, 0.202, 0.01)
        serration_coords_y = [0,0.04]*(len(serration_coords_x)/2)+[0]

        if 'Sr20R21' in dev:
            ax.plot( serration_coords_x, serration_coords_y , 
                    c = 'k', lw = 1)

        #t = ax.text(
        #    0.2, -0.66,
        #    r'${0}$ m/s'.format( vel.replace("U","") ),
        #    zorder = 100, ha='right'
        #)
        #t.set_bbox(dict(color='white', alpha=0.8, edgecolor='white'))


        ax.set_xticks( [ -.2, 0, .2] )
        ax.set_ylim( -0.3, 0.3 )
        ax.set_xlim( dev_case_df.x.min() , dev_case_df.x.max() )
        ax.set_aspect( 'equal' )
        if dev == 'Sr20R21':
            ax.set_xlabel( "$z$ [m]" )
        ax.set_aspect( 'equal' )

    axes[0].set_ylabel( "$x$ [m]" )
    axes[1].set_xlabel( "$z$ [m]" )

    fig.subplots_adjust(right = 0.87)
    cbar_ax = fig.add_axes([0.9, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(contf, cax=cbar_ax)
    cbar.ax.set_ylabel('SPL [dB]')

    plt.savefig(
        'results/SourceMap.png',
        bbox_inches='tight'
    )


                

