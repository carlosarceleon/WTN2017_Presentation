import publish
from os.path import isfile
from os import remove
#import case_dict_overall_correction as case_dict
#cases_df = case_dict.return_case_df()

#markers = ['o', 'v',  's', 'p', '^', '<', '>', '8',
#                  '*', 'h', 'H', 'D', 'd']
markers = ['x','2','+','o','s']

def get_color_and_marker(case_name):

    if "STE" in case_name:
        color  = (0.0, 0.4470588235294118, 0.6980392156862745)
        marker = 'x'
        cmap = 'Blues'
    elif 'z00' in case_name or 'loc00' in case_name:
        color = (0.0, 0.6196078431372549, 0.45098039215686275)
        marker = '2'
        cmap = 'Greens'
    elif 'z05' in case_name or 'loc05' in case_name:
        color = (0.8352941176470589, 0.3686274509803922, 0.0)
        marker = '+'
        cmap = 'Oranges'
    elif 'z10' in case_name or 'loc10' in case_name:
        color = (0.8, 0.4745098039215686, 0.6549019607843137)
        marker = 'o'
        cmap = 'RdPu'

    else: print case_name; return 0,0
    return color,marker,cmap
                        

#####################################################
#####################################################
# Analyze the boundary layer at the very trailing
# edge
#####################################################
#####################################################

def get_tr_streamlined_surface(z_loc = 0):
    import article2_time_averaged_routines as tar
    import tr_case_dict_overall_correction as case_dict
    cases_df = case_dict.return_case_df()

            
    pickle_files = [
        'data/STE_phi0_alpha0_U20_loc00_tr.dat.p',
        'data/Sr20R21_phi0_alpha0_U20_loc00_tr.dat.p',
        'data/Sr20R21_phi0_alpha0_U20_loc05_tr.dat.p',
        'data/Sr20R21_phi0_alpha0_U20_loc10_tr.dat.p',
    ]
    plot_names = [
        'surfaceplot_u_straight_edge.png',
        'surfaceplot_u_serrated_z00.png',
        'surfaceplot_u_serrated_z05.png',
        'surfaceplot_u_serrated_z10.png',
    ]

    for pf,pn in zip(
        pickle_files,plot_names
    ):
        df = tar.load_df_from_pickle(pf)

        height_correction = cases_df[ cases_df.file == pf ]\
                .height_correction.values[0]

        angle_correction = cases_df[ cases_df.file == pf ]\
                .rotation.values[0]
        
        streamwise_correction = cases_df[ cases_df.file == pf ]\
                .x_corr.values[0]


        tar.show_streamlined_surface_from_df(
            df                    = df,
            variable              = 'u',
            points                = [],
            mask                  = [],
            height_correction     = height_correction,
            angle_correction      = angle_correction,
            x_max                 = 40,
            x_min                 = -2,
            y_max                 = 6,
            y_min                 = -1,
            streamwise_correction = streamwise_correction,
            plot_name             = pn
        )

def get_streamlined_surface(z_loc = 0):
    import article2_time_averaged_routines as tar
    import case_dict_overall_correction as case_dict
    cases_df = case_dict.return_case_df()

            
    pickle_files = [
        'data/ste_phi0_alpha0_u20_loc00.dat_rotated.p',
        'data/sr20r21_phi0_alpha0_u20_loc00.dat_rotated.p',
        'data/sr20r21_phi0_alpha0_u20_loc05.dat_rotated.p',
        'data/sr20r21_phi0_alpha0_u20_loc10_even_newer.dat.p',
    ]
    plot_names = [
        'surfaceplot_u_straight_edge.png',
        'surfaceplot_u_serrated_z00.png',
        'surfaceplot_u_serrated_z05.png',
        'surfaceplot_u_serrated_z10.png',
    ]

    for pf,pn in zip(
        pickle_files,plot_names
    ):
        df = tar.load_df_from_pickle(pf)

        height_correction = cases_df[ cases_df.file == pf ]\
                .height_correction.values[0]

        angle_correction = cases_df[ cases_df.file == pf ]\
                .rotation.values[0]
        
        streamwise_correction = cases_df[ cases_df.file == pf ]\
                .x_corr.values[0]


        tar.show_streamlined_surface_from_df(
            df                    = df,
            variable              = 'u',
            points                = [],
            mask                  = [],
            height_correction     = height_correction,
            angle_correction      = angle_correction,
            x_max                 = 40,
            x_min                 = -2,
            y_max                 = 6,
            y_min                 = -1,
            streamwise_correction = streamwise_correction,
            plot_name             = pn
        )

def write_tr_wall_normal_lines_to_csv():
    import pandas as pd
    import article2_time_averaged_routines as tar
    import case_dict_overall_correction as case_dict
    cases_df = case_dict.return_case_df()

    case_STE = cases_df[ 
        (cases_df.file == 'data/STE_phi0_alpha0_U20_loc00_tr.dat.p') &\
        (cases_df.x_loc == 0) ]

    case_z00_x40 = cases_df[ 
        (cases_df.file == 'data/Sr20R21_phi0_alpha0_U20_loc00_tr.dat.p')&\
        (cases_df.x_loc == 40) ]
    
    case_z00_x20 = cases_df[ 
        (cases_df.file == 'data/Sr20R21_phi0_alpha0_U20_loc00_tr.dat.p')&\
        (cases_df.x_loc == 20) ]
    
    case_z00_x00 = cases_df[ 
        (cases_df.file == 'data/Sr20R21_phi0_alpha0_U20_loc00_tr.dat.p')&\
        (cases_df.x_loc == 0) ]
    
    case_z05 = cases_df[ 
        (cases_df.file == 'data/Sr20R21_phi0_alpha0_U20_loc05_tr.dat.p')&\
        (cases_df.x_loc == 20) ]

    case_z10 = cases_df[ 
        (cases_df.file == \
         'data/Sr20R21_phi0_alpha0_U20_loc10.dat.p')&\
        (cases_df.x_loc == 0) ]

    cases = case_STE
    cases = cases.append(case_z00_x00)
    cases = cases.append(case_z00_x20)
    cases = cases.append(case_z00_x40)
    cases = cases.append(case_z05)
    cases = cases.append(case_z10)

    df_cases = pd.DataFrame()

    for case in cases.iterrows():

        df = tar.get_wall_normal_line(
            pickle_file               = case[1].file,
            x_loc                     = case[1].x_loc,
            variable                  = ['u','v',"u_rms","v_rms",],
            plot                      = True,
            rotation_angle            = 0,
            trust_y_min               = 0,
            height_correction         = 0,
            streamwise_correction     = 0,
        )
        tar.write_to_csv_with_units(
            df, 
            '{0}_x{1:.2f}_z{2}_tr'.format(
                df.trailing_edge.unique()[0],
                df.x.unique()[0],
                df.z.unique()[0],
            )
        )

        df_cases = df_cases.append(df, ignore_index = True)

    tar.write_to_csv_with_units(
        df_cases, 
        'Serrations_BoundaryLayerData_StereoPIV_tr.csv'
    )

def write_wall_normal_lines_to_csv():
    import pandas as pd
    import article2_time_averaged_routines as tar
    import case_dict_overall_correction as case_dict
    cases_df = case_dict.return_case_df()

    case_STE = cases_df[ 
        (cases_df.file == 'data/STE_phi0_alpha0_U20_loc00.dat_rotated.p') &\
        (cases_df.x_loc == 0) ]

    case_z00_x40 = cases_df[ 
        (cases_df.file == 'data/Sr20R21_phi0_alpha0_U20_loc00.dat_rotated.p')&\
        (cases_df.x_loc == 40) ]
    
    case_z00_x20 = cases_df[ 
        (cases_df.file == 'data/Sr20R21_phi0_alpha0_U20_loc00.dat_rotated.p')&\
        (cases_df.x_loc == 20) ]
    
    case_z00_x00 = cases_df[ 
        (cases_df.file == 'data/Sr20R21_phi0_alpha0_U20_loc00.dat_rotated.p')&\
        (cases_df.x_loc == 0) ]
    
    case_z05 = cases_df[ 
        (cases_df.file == 'data/Sr20R21_phi0_alpha0_U20_loc05.dat_rotated.p')&\
        (cases_df.x_loc == 20) ]

    case_z10 = cases_df[ 
        (cases_df.file == \
         'data/Sr20R21_phi0_alpha0_U20_loc10_even_newer.dat.p')&\
        (cases_df.x_loc == 0) ]

    cases = case_STE
    cases = cases.append(case_z00_x00)
    cases = cases.append(case_z00_x20)
    cases = cases.append(case_z00_x40)
    cases = cases.append(case_z05)
    cases = cases.append(case_z10)

    df_cases = pd.DataFrame()

    for case in cases.iterrows():

        df = tar.get_wall_normal_line(
            pickle_file               = case[1].file,
            x_loc                     = case[1].x_loc,
            variable                  = ['u','v',"u_rms","v_rms",],
            plot                      = False,
            rotation_angle            = case[1].rotation,
            trust_y_min               = case[1].y_trust_min,
            height_correction         = case[1].height_correction,
            streamwise_correction     = case[1].x_corr,
        )
        tar.write_to_csv_with_units(
            df, 
            '{0}_x{1:.2f}_z{2}'.format(
                df.trailing_edge.unique()[0],
                df.x.unique()[0],
                df.z.unique()[0],
            )
        )

        df_cases = df_cases.append(df, ignore_index = True)

    tar.write_to_csv_with_units(
        df_cases, 
        'Serrations_BoundaryLayerData_StereoPIV_v2.csv'
    )

def get_trailing_edge_for_all_TR_cases_at_TE_m1():
    import tr_case_dict_overall_correction as case_dict
    cases_df = case_dict.return_case_df()

    case_STE = cases_df[ 
        (cases_df.file == 'data/STE_phi0_alpha0_U20_loc00_tr.dat.p') &\
        (cases_df.x_loc == -1) ]
    case_z00 = cases_df[ 
        (cases_df.file == 'data/Sr20R21_phi0_alpha0_U20_loc00_tr.dat.p')&\
        (cases_df.x_loc == 39) ]
    case_z05 = cases_df[ 
        (cases_df.file == 'data/Sr20R21_phi0_alpha0_U20_loc05_tr.dat.p')&\
        (cases_df.x_loc == 19) ]
    case_z10 = cases_df[ 
        (cases_df.file == \
         'data/Sr20R21_phi0_alpha0_U20_loc10_tr.dat.p')&\
        (cases_df.x_loc == -1) ]
    case_z00_2 = cases_df[ 
        (cases_df.file == 'data/Sr20R21_phi0_alpha0_U20_loc00_tr.dat.p')&\
        (cases_df.x_loc == 34) ]
    case_z05_2 = cases_df[ 
        (cases_df.file == 'data/Sr20R21_phi0_alpha0_U20_loc05_tr.dat.p')&\
        (cases_df.x_loc == 14) ]
    case_z10_2 = cases_df[ 
        (cases_df.file == \
         'data/Sr20R21_phi0_alpha0_U20_loc10_tr.dat.p')&\
        (cases_df.x_loc == 3) ]
    case_z10_2 = cases_df[ 
        (cases_df.file == \
         'data/Sr20R21_phi0_alpha0_U20_loc10_tr.dat.p')&\
        (cases_df.x_loc == -5) ]

    cases = case_STE
    cases = cases.append(case_z00, ignore_index = True)
    cases = cases.append(case_z05, ignore_index = True)
    cases = cases.append(case_z10, ignore_index = True)
    cases = cases.append(case_z00_2, ignore_index = True)
    cases = cases.append(case_z05_2, ignore_index = True)
    cases = cases.append(case_z10_2, ignore_index = True)

    schematic = '/home/carlos/Documents/PhD/Articles/Article_2/'+\
            'Figures/measurement_locations_TE_m2.png'

    write_boundary_layers(cases)
    plot_cases( cases , plot_name = "TR_At_trailing_edge", 
               schematic = schematic, time_resolved = True)

def get_trailing_edge_for_all_TR_cases_at_x_m1():
    import tr_case_dict_overall_correction as case_dict
    cases_df = case_dict.return_case_df()

    case_STE = cases_df[ 
        (cases_df.file == 'data/STE_phi0_alpha0_U20_loc00_tr.dat.p') &\
        (cases_df.x_loc == -1) ]
    case_z00 = cases_df[ 
        (cases_df.file == 'data/Sr20R21_phi0_alpha0_U20_loc00_tr.dat.p')&\
        (cases_df.x_loc == -1) ]
    case_z05 = cases_df[ 
        (cases_df.file == 'data/Sr20R21_phi0_alpha0_U20_loc05_tr.dat.p')&\
        (cases_df.x_loc == -1) ]
    case_z10 = cases_df[ 
        (cases_df.file == \
         'data/Sr20R21_phi0_alpha0_U20_loc10_tr.dat.p')&\
        (cases_df.x_loc == -1) ]

    cases = case_STE
    cases = cases.append(case_z00, ignore_index = True)
    cases = cases.append(case_z05, ignore_index = True)
    cases = cases.append(case_z10, ignore_index = True)

    schematic = '/home/carlos/Documents/PhD/Articles/Article_2/'+\
            'Figures/measurement_locations_x0_m2.png'

    write_boundary_layers(cases)
    plot_cases( cases , plot_name = "TR_At_x_m1", 
               schematic = schematic, time_resolved = True)

def get_trailing_edge_for_all_cases_at_TE_m1():
    import case_dict_overall_correction as case_dict
    cases_df = case_dict.return_case_df()

    case_STE = cases_df[ 
        (cases_df.file == 'data/STE_phi0_alpha0_U20_loc00.dat_rotated.p') &\
        (cases_df.x_loc == -1) ]
    case_z00 = cases_df[ 
        (cases_df.file == 'data/Sr20R21_phi0_alpha0_U20_loc00.dat_rotated.p')&\
        (cases_df.x_loc == 39) ]
    case_z05 = cases_df[ 
        (cases_df.file == 'data/Sr20R21_phi0_alpha0_U20_loc05.dat_rotated.p')&\
        (cases_df.x_loc == 19) ]
    case_z10 = cases_df[ 
        (cases_df.file == \
         'data/Sr20R21_phi0_alpha0_U20_loc10_even_newer.dat.p')&\
        (cases_df.x_loc == -1) ]

    cases = case_STE
    cases = cases.append(case_z00, ignore_index = True)
    cases = cases.append(case_z05, ignore_index = True)
    cases = cases.append(case_z10, ignore_index = True)

    schematic = '/home/carlos/Documents/PhD/Articles/Article_2/'+\
            'Figures/measurement_locations_TE_m2.png'

    write_boundary_layers(cases)
    plot_cases( cases , plot_name = "At_trailing_edge", schematic = schematic)

def get_trailing_edge_for_all_cases_at_x_m1():
    import case_dict_overall_correction as case_dict
    cases_df = case_dict.return_case_df()

    case_STE = cases_df[ 
        (cases_df.file == 'data/STE_phi0_alpha0_U20_loc00.dat_rotated.p') &\
        (cases_df.x_loc == -1) ]
    case_z00 = cases_df[ 
        (cases_df.file == 'data/Sr20R21_phi0_alpha0_U20_loc00.dat_rotated.p')&\
        (cases_df.x_loc == -1) ]
    case_z05 = cases_df[ 
        (cases_df.file == 'data/Sr20R21_phi0_alpha0_U20_loc05.dat_rotated.p')&\
        (cases_df.x_loc == -1) ]
    case_z10 = cases_df[ 
        (cases_df.file == \
         'data/Sr20R21_phi0_alpha0_U20_loc10_even_newer.dat.p')&\
        (cases_df.x_loc == -1) ]

    cases = case_STE
    cases = cases.append(case_z00, ignore_index = True)
    cases = cases.append(case_z05, ignore_index = True)
    cases = cases.append(case_z10, ignore_index = True)

    schematic = '/home/carlos/Documents/PhD/Articles/Article_2/'+\
            'Figures/measurement_locations_x0_m2.png'

    write_boundary_layers(cases)
    plot_cases( cases , plot_name = "At_x_m1", schematic = schematic)

def plot_cases(cases,plot_name = '', schematic = '', time_resolved = False):
    import article2_time_averaged_routines as tar
    import matplotlib.pyplot as plt
    from matplotlib import rc
    import seaborn as sns
    from numpy import linspace,log,diff
    import pandas as pd
    from matplotlib.cbook import get_sample_data
    from os.path import isfile,split,splitext

    #if time_resolved: 
    #    import tr_case_dict_overall_correction as case_dict
    #    cases = case_dict.return_case_df()

    markeredgewidth = 2
    markerfacecolor = 'none'
    markersize      = 12
    mew             = 4 # Marker edge width

    """ Cite
    Nagano1998
    Lee2008
    """

    rc('text',usetex=True)
    rc('font',weight='normal')

    sns.set_context('paper')
    sns.set(font='serif',font_scale=3.0,style='whitegrid')
    rc('font',family='serif', serif='Linux Libertine')

    figsize = (7,6)
    fig_nondim   , axes_nondim   = plt.subplots(1,1,sharex=True,sharey=True,
                                                figsize=figsize )
    fig_u        , axes_u        = plt.subplots(1,1,sharex=True,sharey=True,
                                                figsize=figsize )
    fig_v        , axes_v        = plt.subplots(1,1,sharex=True,sharey=True,
                                                figsize=figsize )
    fig_diff     , axes_diff     = plt.subplots(1,1,sharex=True,sharey=True,
                                                figsize=figsize )
    fig_diff_rms , axes_diff_rms = plt.subplots(1,1,sharex=True,sharey=True,
                                                figsize=figsize )
    fig_u_rms    , axes_u_rms    = plt.subplots(1,1,sharex=True,sharey=True,
                                                figsize=figsize )
    fig_v_rms    , axes_v_rms    = plt.subplots(1,1,sharex=True,sharey=True,
                                                figsize=figsize )
    fig_w        , axes_w        = plt.subplots(1,1,sharex=True,sharey=True,
                                                figsize=figsize )
    fig_reyn     , axes_reyn     = plt.subplots(1,1,sharex=True,sharey=True,
                                                figsize=figsize )

    df_all_cases = pd.DataFrame()

    for case in cases.iterrows():
        print "Processing {0}".format(case[1].file)

        df_pickle_name = "ReservedData/{0}_{1}_{2}_{3}_{4}.p".format(
            splitext(split(case[1].file)[1])[0],
            case[1].x_loc,
            case[1].rotation,
            case[1].height_correction,
            case[1].x_corr
        )


        if isfile(df_pickle_name):
            df = pd.read_pickle(df_pickle_name)
        else:
            df = tar.get_wall_normal_line(
                pickle_file               = case[1].file,
                x_loc                     = case[1].x_loc,
                variable                  = ['u','v','w',"u_rms","v_rms",
                                             'Reynold_stress_uv' ],
                plot                      = False,
                rotation_angle            = case[1].rotation,
                trust_y_min               = case[1].y_trust_min,
                height_correction         = case[1].height_correction,
                streamwise_correction     = case[1].x_corr,
            )

        U_e, delta_99, delta_displacement, delta_momentum = \
                tar.get_boundary_layer_values( df )

        color, marker, cmap = get_color_and_marker( case[1].file )

        df = tar.get_dimensionless_inner_variables(
            df,
            correction = 0,
            Cf = case[1].Cf,
        )

        axes_nondim.plot(
            df.y_plus,
            df.u_plus,
            label = case[1].case_name,
            marker = case[1].marker,
            markeredgewidth = markeredgewidth,
            markerfacecolor = markerfacecolor,
            markeredgecolor = color,
            markersize      = markersize,
            mew             = mew,
            color = [float(f) for f in case[1].color.split(',')]
        )
        axes_u.plot(
            #df.u/20.,
            df.u/U_e,
            df.y/delta_99,
            label = case[1].case_name,
            marker = case[1].marker,
            markeredgewidth=markeredgewidth,
            markerfacecolor = markerfacecolor,
            markeredgecolor = color,
            markersize      = markersize,
            mew             = mew,
            color = [float(f) for f in case[1].color.split(',')]
        )

        axes_u_rms.plot(
            df.u_rms/U_e,
            df.y/delta_99,
            label = case[1].case_name,
            marker = case[1].marker,
            markeredgewidth=markeredgewidth,
            markerfacecolor = markerfacecolor,
            markeredgecolor = color,
            markersize      = markersize,
            mew             = mew,
            color = [float(f) for f in case[1].color.split(',')]
        )

        axes_v_rms.plot(
            df.v_rms/U_e,
            df.y/delta_99,
            label = case[1].case_name,
            marker = case[1].marker,
            markeredgewidth=markeredgewidth,
            markerfacecolor = markerfacecolor,
            markeredgecolor = color,
            markersize      = markersize,
            mew             = mew,
            color = [float(f) for f in case[1].color.split(',')]
        )

        axes_diff.plot(
            df.y.values[0:-1],
            diff(df.u)/diff(df.y),
            label = case[1].case_name,
            marker = case[1].marker,
            markeredgewidth=markeredgewidth,
            markerfacecolor = markerfacecolor,
            markeredgecolor = color,
            markersize      = markersize,
            mew             = mew,
            color = [float(f) for f in case[1].color.split(',')]
        )
        
        axes_diff_rms.plot(
            df.y.values[0:-1],
            diff(df.u_rms)/diff(df.y),
            label = case[1].case_name,
            marker = case[1].marker,
            markeredgewidth=markeredgewidth,
            markerfacecolor = markerfacecolor,
            markeredgecolor = color,
            markersize      = markersize,
            mew             = mew,
            color = [float(f) for f in case[1].color.split(',')]
        )

        axes_v.plot(
            df.v/U_e,
            df.y/delta_99,
            label = case[1].case_name,
            marker = case[1].marker,
            markeredgewidth=markeredgewidth,
            markerfacecolor = markerfacecolor,
            markeredgecolor = color,
            markersize      = markersize,
            mew             = mew,
            color = [float(f) for f in case[1].color.split(',')]
        )

        axes_w.plot(
            df.w/U_e,
            df.y/delta_99,
            label = case[1].case_name,
            marker = case[1].marker,
            markeredgewidth=markeredgewidth,
            markerfacecolor = markerfacecolor,
            markeredgecolor = color,
            markersize      = markersize,
            mew             = mew,
            color = [float(f) for f in case[1].color.split(',')]
        )

        axes_reyn.plot(
            -df.Reynold_stress_uv/U_e**2 * 1000,
            df.y/delta_99,
            ls = '',
            label = case[1].case_name,
            marker = case[1].marker,
            markeredgewidth=markeredgewidth,
            markerfacecolor = markerfacecolor,
            markeredgecolor = color,
            markersize      = markersize,
            mew             = mew,
            color = [float(f) for f in case[1].color.split(',')]
        )

        df_all_cases = df_all_cases.append(
            df,
            ignore_index=True
        )
        if not isfile(df_pickle_name):
            df.to_pickle(df_pickle_name)

    # Plot u_plus = y_plus law
    line = axes_nondim.plot(
        linspace(1e0,1.2e1),
        linspace(1e0,1.2e1),
        '--'
    )
    # Plot u_plus = 1/0.4 * ln y_plus law
    axes_nondim.plot(
        linspace(8e0,1e3),
        log(linspace(8e0,1e3))/0.4+5,
        '--',
        color = line[0].get_color()
    )

    axes_nondim.set_xlabel("$y^+$")
    axes_nondim.set_ylabel("$u^+$")
    axes_nondim.set_xscale('log')
    axes_nondim.set_xlim(1e0, 1e3)
    axes_nondim.xaxis.grid(True, which='minor')

    list_of_axes = [
        axes_u,
        axes_v,
        axes_w,
        axes_v_rms,
        axes_u_rms,
        axes_diff,
        axes_diff_rms,
        axes_nondim,
        axes_reyn
    ]

    list_of_figures = [
        fig_u,
        fig_v,
        fig_w,
        fig_v_rms,
        fig_u_rms,
        fig_diff,
        fig_diff_rms,
        fig_nondim,
        fig_reyn
    ]

    if not schematic:
        for ax in list_of_axes:
            ax.legend(
                bbox_to_anchor = (0., 1.02, 1.0, .102),
                loc            = 3,
                ncol           = 2,
                mode           = "expand",
                borderaxespad  = 0.
            )

    if schematic:
        im = plt.imread( get_sample_data( schematic  ) )
        for fig in list_of_figures:
            if fig == fig_v_rms or fig == fig_u_rms or fig == fig_w\
               or fig == fig_reyn:
                newax = fig.add_axes([0.55, 0.55, 0.3, 0.3], anchor = 'SW', 
                                     zorder=15)
            else:
                newax = fig.add_axes([0.175, 0.55, 0.3, 0.3], anchor = 'SW', 
                                     zorder=15)
            newax.imshow(im)
            newax.axis('off')

    ylabel = r'$y / \delta_{99}$'
    axes_u.set_xlabel("$\overline u/u_e$")
    axes_u.set_ylabel(ylabel)
    axes_u.set_xlim(0, 1.1)
    axes_v.set_xlabel("$\overline v/u_e$")
    axes_v.set_ylabel(ylabel)
    axes_w.set_xlabel("$\overline w/u_e$")
    axes_w.set_ylabel(ylabel)
    axes_w.set_xlim(-0.04, 0.04)
    axes_w.set_xticks([-0.04, -0.02, 0, 0.02, 0.04])
    axes_u_rms.set_xlabel("$u'_\\textrm{rms}/u_e$")
    axes_u_rms.set_ylabel(ylabel)
    axes_u_rms.set_xlim(0,0.12)
    axes_v_rms.set_xlabel("$v'_\\textrm{rms}/u_e$")
    axes_v_rms.set_ylabel(ylabel)
    axes_v_rms.set_xlim(0,0.12)
    axes_diff.set_xlabel("$y$ [mm]")
    axes_diff.set_ylabel("$\\partial \overline u/ \\partial y$")
    axes_diff.set_yscale("log")
    axes_diff_rms.set_xlabel("$y$ [mm]")
    axes_diff_rms.set_ylabel("$\\partial u'_\\textrm{rms}/ \\partial y$")

    axes_reyn.set_xlabel(r"$-\overline{u'v'}/u_e^2 \times 10^{-3}$")
    axes_reyn.set_ylabel(ylabel)
    axes_reyn.set_xlim(-1.0, 6.0)
    axes_reyn.set_ylim(0, 1.5)

    axes_u.set_ylim(0, 1.5)
    axes_v.set_ylim(0, 1.5)
    axes_w.set_ylim(0, 1.5)
    axes_u_rms.set_ylim(0, 1.5)
    axes_v_rms.set_ylim(0, 1.5)
    axes_diff.set_xlim(left = 0)
    axes_diff_rms.set_xlim(left = 0)

    fig_u.savefig(
        'images/{0}_BoundaryLayerAnalysis_u.png'.format(plot_name),
        bbox_inches='tight'
    )
    fig_v.savefig(
        'images/{0}_BoundaryLayerAnalysis_v.png'.format(plot_name),
        bbox_inches='tight'
    )
    fig_w.savefig(
        'images/{0}_BoundaryLayerAnalysis_w.png'.format(plot_name),
        bbox_inches='tight'
    )
    fig_u_rms.savefig(
        'images/{0}_BoundaryLayerAnalysis_urms.png'.format(plot_name),
        bbox_inches='tight'
    )
    fig_v_rms.savefig(
        'images/{0}_BoundaryLayerAnalysis_vrms.png'.format(plot_name),
        bbox_inches='tight'
    )
    fig_diff.savefig(
        'images/{0}_BoundaryLayerAnalysis_dudy.png'.format(plot_name),
        bbox_inches='tight'
    )
    fig_diff_rms.savefig(
        'images/{0}_BoundaryLayerAnalysis_durmsdy.png'.format(plot_name),
        bbox_inches='tight'
    )
    fig_reyn.savefig(
        'images/{0}_BoundaryLayerAnalysis_ReynoldsStressUV.png'\
        .format(plot_name),
        bbox_inches='tight'
    )
    fig_nondim.savefig(
        'images/{0}_NonDimensionalBoundaryLayerAnalysis_u.png'.format(plot_name),
        bbox_inches='tight'
    )
    tar.write_to_csv_with_units(
        df_all_cases,
        "{0}_BoundaryLayerData_StereoPIV.csv".format(plot_name),
    )

def analyze_edge_velocities(df):
    from numpy import arange
    import article2_time_averaged_routines as tar
    import pandas as pd
    import matplotlib.pyplot as plt
    reload(tar)

    # This function was once used to evaluate the selection of 
    # U_e based on different parameters

    conditions = [
        'vorticity_integration_rate_of_change',
        #'u_rms_rate_of_change',
        #'v_rms_rate_of_change',
        #'u_rate_of_change',
    ]

    thresholds = arange(1,100,0.1)[::-1]

    U_e_df = pd.DataFrame()
    for c in conditions:
        U_e = []
        for t in thresholds:
            U_e.append(tar.get_edge_velocity(df, condition = c, threshold = t))
        U_e_df[c] = U_e

    U_e_df.plot()
    plt.axhline( df.u.max() )
    plt.show()

    return U_e_df

def write_boundary_layers(cases, time_resolved = False):
    import article2_time_averaged_routines as tar
    from os.path import isfile,splitext,split
    import pandas as pd

    for case in cases.iterrows():

        print case[1].x_loc, case[1].file
        df_pickle_name = "ReservedData/{0}_{1}_{2}_{3}_{4}.p".format(
            splitext(split(case[1].file)[1])[0],
            case[1].x_loc,
            case[1].rotation,
            case[1].height_correction,
            case[1].x_corr
        )

        if isfile(df_pickle_name):
            df = pd.read_pickle(df_pickle_name)
        else:
            df = tar.get_wall_normal_line(
                pickle_file               = case[1].file,
                x_loc                     = case[1].x_loc,
                variable                  = ['u','v',"u_rms","v_rms",],
                plot                      = False,
                rotation_angle            = case[1].rotation,
                trust_y_min               = case[1].y_trust_min,
                height_correction         = case[1].height_correction,
                streamwise_correction     = case[1].x_corr,
            )

        df = tar.get_dimensionless_inner_variables(
            df,
            correction = 0,
            Cf = case[1].Cf,
        )

        df['x_loc'] = case[1].x_loc


        if time_resolved:
            bl_file.replace('.csv','_TR.csv') 
        tar.write_boundary_layers(df, boundary_layers_file = bl_file)

bl_file = "Boundary_layer_information.csv"
if isfile(bl_file):
    remove(bl_file)

#get_tr_streamlined_surface(z_loc = 0)

#get_trailing_edge_for_all_cases_at_TE_m1()
#get_trailing_edge_for_all_cases_at_x_m1()

get_trailing_edge_for_all_TR_cases_at_TE_m1()
#get_trailing_edge_for_all_TR_cases_at_x_m1()

#write_wall_normal_lines_to_csv()
#publish.publish()
