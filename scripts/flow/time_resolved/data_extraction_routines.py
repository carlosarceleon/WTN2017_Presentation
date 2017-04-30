
def extract_relevant_data( case_list = [], exceptions = [], y_delta_locs = [],
                         x_2h_locs = [] , plot = False, output_root = '',
                         overwrite = False, just_plot = False ):
    """ This will extract the wall normal data at the spanwise location
    TE at a certain y density
    """

    from os                           import listdir
    from os.path                      import join,split,isfile
    from pandas                       import DataFrame, HDFStore
    from boundary_layer_routines      import return_bl_parameters
    from progressbar                  import ProgressBar,Percentage
    from progressbar                  import Bar,ETA,SimpleProgress
    from numpy                        import array, round, linspace, nan, sign
    from data_cleaning_routines       import show_surface_from_df

    x_2h_locs    = round( array( x_2h_locs ),    2 )
    y_delta_locs = round( array( y_delta_locs ), 2 )

    if not output_root:
        output_root = '/media/carlos/6E34D2CD34D29783/' +\
                '2015-02_SerrationPIV/TR_Data_Location_Calibrated_Article3'

    # Get the available HDF5 files #############################################
    hdf5_root = '/media/carlos/6E34D2CD34D29783/' +\
                '2015-02_SerrationPIV/TR_Data_Location_Calibrated_Article3'

    if not len(case_list):
        hdf5_files = [f for f in listdir( hdf5_root ) \
                      if f.endswith('.hdf5') \
                      and not f in exceptions ]
    else:
        hdf5_files = [f for f in listdir( hdf5_root ) \
                      if f.endswith('.hdf5') \
                      and f in case_list ]
    # ##########################################################################

    for hf in [join( hdf5_root, f ) for f in hdf5_files]:

        suffix = '_AirfoilNormal'

        mean_flow_rotated = False
        if hf.endswith("MeanFlowRotated.hdf5"):
            mean_flow_rotated = True
            suffix = ('_MeanFlowRotated')

        f = split( hf )[1].replace(suffix,'')\
                .replace('_Aligned.hdf5','').replace('.hdf5','')
        if mean_flow_rotated:
            f = f + '_mfr'

        # Get the normalization locations depending on the case ################
        if 'z00' in f and not 'STE' in f:
            x_bl_loc = 38
        elif 'z05' in f:
            x_bl_loc = 18
        elif 'z10' in f or 'STE' in f:
            x_bl_loc = -2

        #if mean_flow_rotated:
        #    mean_flow_df = read_pickle( mean_flow_pickle )
        #    mean_flow_angle = mean_flow_df.T.ix[ 
        #        f.replace('.hdf5','')
        #    ].ix[0].values[0]

        #    yshift = - x_bl_loc * sin( radians( mean_flow_angle ) )
        #else:
        #    yshift = 0

        print "   Processing {0}".format( f )
        if not overwrite and isfile( 'ReservedData/' + f + '.hdf5' ):
            print "   File exists; skipping"
            continue


        hdf_t = HDFStore( hf, 'r' )

        # Get the available coordinates ########################################
        hf_coords = hdf_t.select('data', where = [ 't = 0' ], 
                                 columns = [ 'x', 'y' , 'u' ] )
        # ######################################################################

        # Turn the non-dim requested locations into physical coords ############
        requested_locations = []
        requested_normalized_locations = []
        #for x,x_norm in zip(x_2h_locs * tooth_length, x_2h_locs):
        #    for y_d in y_delta_locs:
        #        bl_params = return_bl_parameters( f , [x] )
        #        d_99 = bl_params.delta_99.values[0]
        #        #if "STE" in f:
        #        #    d_99 = 9.4
        #        y = y_d * d_99
        #        requested_locations.append( (x,y) )
        #        requested_normalized_locations.append( ( x_norm, y_d ) )

        bl_params = return_bl_parameters( 
            f.replace( '_mfr', '' ) , 
            [x_bl_loc] )
        d_99 = bl_params.delta_99.values[0]

        for x,x_norm in zip(x_2h_locs * tooth_length, x_2h_locs):
            for y_d in y_delta_locs:
                y = y_d * d_99
                requested_locations.append( (x,y) )
                requested_normalized_locations.append( ( x_norm, y_d ) )
        print "    Normalizing to a BL thickness of {0:.2f} mm".\
                format(d_99)
        # ######################################################################

        hf_coords = hf_coords.replace( 0, nan )
        hf_coords = hf_coords.dropna()

        min_x_real_loc = min( x_2h_locs ) * tooth_length
        max_x_real_loc = max( x_2h_locs ) * tooth_length

        # This is the max y location for a full square in the rotated FoV ######
        max_y_field_of_view = min(
            hf_coords[ 
                ( hf_coords.x < max_x_real_loc * 1.1 ) & \
                ( hf_coords.x > max_x_real_loc * 0.9 ) 
            ].y.max(),

            hf_coords[ 
                ( hf_coords.x < \
                 min_x_real_loc * ( 1 + sign( min_x_real_loc ) * 0.1 ) ) & \
                ( hf_coords.x > \
                 min_x_real_loc * ( 1 - sign( min_x_real_loc ) * 0.1 ) ) 
            ].y.max()
        )
        # ######################################################################

        # This is the minimum (and x_loc aligned) y location for a full square
        # in the rotated FoV ###################################################
        min_y_field_of_view = max(
            hf_coords[ 
                ( hf_coords.x < \
                 min_x_real_loc * ( 1 + sign( min_x_real_loc ) * 0.1 ) ) & \
                ( hf_coords.x > \
                 min_x_real_loc * ( 1 - sign( min_x_real_loc ) * 0.1 ) ) 
            ].y.min(), 
            hf_coords[ 
                ( hf_coords.x < max_x_real_loc * 1.1 ) & \
                ( hf_coords.x > max_x_real_loc * 0.9 ) 
            ].y.min()
        )

        # ######################################################################

        if max( y_delta_locs ) * d_99 < max_y_field_of_view:
            max_y_loc = max( y_delta_locs ) * d_99
        else:
            max_y_loc = max_y_field_of_view

        if min( y_delta_locs ) * d_99 > min_y_field_of_view:
            min_y_loc = min( y_delta_locs ) * d_99
        else:
            min_y_loc = min_y_field_of_view

        available_xy_locs = hf_coords[
            ( hf_coords.x > min( x_2h_locs ) * 40. ) & \
            ( hf_coords.x < max( x_2h_locs ) * 40. ) & \
            ( hf_coords.y > min_y_loc ) & \
            ( hf_coords.y < max_y_loc )
        ][ ['x','y'] ]

        if available_xy_locs.empty:
            return 0
              
        available_xy_locs = [ tuple(x) for x in available_xy_locs.values ]

        #if mean_flow_rotated:
        #    y_correction = - x_bl_loc * sin( radians( mean_flow_angle ) )

        #if plot:

        #    trailing_edge,phi,alpha,U,z = decript_case_name( f )

        #    if trailing_edge   == 'serrated': device = 'Sr20R21'
        #    elif trailing_edge == 'straight': device = 'STE'
        #    elif trailing_edge == 'slitted':  device = 'Slit20R21'

        #    case_name = "{0}_phi{1}_alpha{2}_U{3}_loc{4}_tr.dat".format(
        #        device, phi, alpha, U, z
        #    )

        #    plot_name = ''
        #    if mean_flow_rotated:
        #        average_pickle_name = 'averaged_data/' + case_name + \
        #                '_mean_flow_rotated.p'
        #    else:
        #        average_pickle_name = 'averaged_data/' + case_name + '.p'
        #    df_av = read_pickle( average_pickle_name )

        #    if mean_flow_rotated:
        #        df_av.y = df_av.y - yshift

        #    plot_name = 'ReservedData/' + f + '.png'
        #    show_surface_from_df( 
        #        df_av , 
        #        points           = available_xy_locs ,
        #        plot_name        = plot_name,
        #        point_correction = ( 0 , 0 )
        #    )
        #    if not plot_name:
        #        continue

        # Don't try to get it all at once; split the vertical in 4 pieces
        vertical_split_blocks = 10

        y_ranges = linspace( 
            min_y_loc,
            max_y_loc,
            vertical_split_blocks
        ) 

        xmin = min(x_2h_locs) * tooth_length
        xmax = max(x_2h_locs) * tooth_length

        print "   Extracting data from {0}".format(f)
        print "     at the streamwise locations:",
        print "     {0:.2f} to {1:.2f}".format( xmin, xmax )
        print "     at the wall-normal locations:",
        print "     {0:.2f} to {1:.2f}".format( min_y_loc, max_y_loc )

        if plot:
            query = " x>={0} & x<{1} & y>={2} & y<{3} & t=0".\
                    format( xmin, xmax, min_y_loc, max_y_loc )

            df0 = hdf_t.select(
                key   = 'data',
                where = [ query ],
            )

            plot_name = 'ReservedData/' + f + '.png'
            if mean_flow_rotated:
                plot_name = plot_name.replace( '.png', '_mfr.png' )

            show_surface_from_df( 
                df0,
                plot_name        = plot_name,
                point_correction = ( 0 , 0 )
            )
            if just_plot:
                return 0

        query   = ''
        cnt_all = 0

        cnt = 0
        time_series_hdf = HDFStore( 'ReservedData/' + f + '.hdf5', 'w' )

        progress = ProgressBar(
             widgets=[
                 Bar(),' ',
                 Percentage(),' ',
                 ETA(), ' (query bunch  ',
                 SimpleProgress(),')'], 
             maxval = vertical_split_blocks
             ).start()

        for ymin, ymax in zip( y_ranges[:-1], y_ranges[1:] ):

            query = " x>={0} & x<{1} & y>={2} & y<{3} ".\
                    format( xmin, xmax, ymin, ymax )

            df_t = hdf_t.select(
                key   = 'data',
                where = [ query ],
            )

            df_t['near_x_2h']    = round( df_t.x / 40.,  4 )
            df_t['near_y_delta'] = round( df_t.y / d_99, 4 )

            # If it has been flow rotated, the zero y direction needs to be
            # brought back to the location of the wall at the reference x
            # location #########################################################
            #if mean_flow_rotated:
            #    df_t.y = df_t.y - y_correction
            # ##################################################################
            # The translation now happens when pickling the raw data

            if not cnt:
                time_series_hdf.put( 'data', df_t , 
                                    data_columns = [
                                        'near_x_2h',
                                        'near_y_delta',
                                        't'
                                    ],
                                    format = 't')
            else:
                time_series_hdf.append( 'data', df_t , 
                                       data_columns = [
                                           'near_x_2h',
                                           'near_y_delta',
                                           't'
                                       ],
                               format = 't')

            cnt_all += 1
            cnt     += 1

            progress.update(cnt_all)

            df_t = DataFrame()


        progress.finish()
        hdf_t.close()
        time_series_hdf.close()


def run_data_extraction( output_root = '', overwrite = False, root = '' ):
    from os import listdir
    from numpy import arange

    if not root:
        hdf5_root = '/media/carlos/6E34D2CD34D29783/' + \
                    '2015-02_SerrationPIV/TR_Data_Location_Calibrated_Article3'
    else:
        hdf5_root = root

    if not output_root:
        output_root = '/home/carlos/Documents/PhD/Articles/Article_3/' + \
                'Scripts/time_resolved/ReservedData'

    print hdf5_root
    case_list = sorted( [
        f for f in listdir( hdf5_root ) \
        if f.endswith('.hdf5')
    ] )

    #case_list = sorted( [
    #    f for f in case_list \
    #    #if f.endswith('MeanFlowRotated.hdf5') and '12_p6' in f
    #    #if f.endswith('MeanFlowRotated.hdf5') and '12_p6' in f
    #    #if 'STE' in f
    #] )

    print "   Going to retrieve selected locations from the following files"
    for c in case_list:
        print "      "+c
    
    #y_delta_locs = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    y_delta_locs = arange( 0.1, 2.0, 0.02 )

    for cl in case_list:
        x_2h_locs = []
        if cl.endswith('z00_tr_Aligned.hdf5') and not cl.startswith('STE'):
            x_2h_locs = arange( 0.76, 1.02, 0.01 ) 
        elif cl.endswith('z05_tr_Aligned.hdf5'):
            x_2h_locs = arange( 0.26, 0.52, 0.01 ) 
        elif cl.endswith('_AirfoilNormal.hdf5'):
            x_2h_locs = arange( -0.1, 0.16, 0.01 ) 
        elif cl.endswith('_MeanFlowRotated.hdf5') and 'z05' in cl:
            x_2h_locs = arange( 0.26, 0.52, 0.01 ) 

        if len( x_2h_locs ):
            extract_relevant_data( 
                case_list    = [cl],
                exceptions   = [],
                y_delta_locs = y_delta_locs,
                x_2h_locs    = x_2h_locs,
                plot         = True,
                output_root  = output_root,
                overwrite    = overwrite,
                just_plot    = False
            )

# CONSTANTS ####################################################################

tooth_length = 40

# ##############################################################################

mean_flow_pickle = 'mean_flow_angle.p'
