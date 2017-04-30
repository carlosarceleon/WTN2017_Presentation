def return_case_df():
    import pandas as pd

    ste_height_correction   = -0.4
    ste_rotation_correction = -1.0
    ste_x_correction        = -2.5

    z00_height_correction   = 0.1
    z00_rotation_correction = -1.0
    z00_x_correction        = -1.0

    z05_height_correction   = -0.2
    z05_rotation_correction = -2.0
    z05_x_correction        = -0.5

    z10_height_correction   =  0.1
    z10_rotation_correction =  0.0
    z10_x_correction        =  0.0

    case_dicts = [ 

        # x = -1.0 ############################################################
        {
            'file':              'data/STE_phi0_alpha0_U20_loc00.dat_rotated.p',
            'case_name':         'Straight edge',
            'type':              'surface',
            'x_loc':             -1.0,
            'y_trust_min':       1.0,
            'Cf':                0.0015,
            'height_correction': ste_height_correction,
            'rotation':          ste_rotation_correction,
            'x_corr':            ste_x_correction,
            'marker':            'x',
            'color':             '0.0, 0.4470588235294118, 0.6980392156862745'
        },
        {
            'file':              'data/Sr20R21_phi0_alpha0_U20_loc00.dat_rotated.p',
            'case_name':         'Serrated, $z/\\lambda = 0$',
            'type':              'surface',
            'x_loc':             -1.0,
            'y_trust_min':       0.4,
            'Cf':                0.0015,
            'height_correction': z00_height_correction,
            'rotation':          z00_rotation_correction,
            'x_corr':            z00_x_correction,
            'marker':            '2',
            'color':             '0.0, 0.6196078431372549, 0.45098039215686275'
        },
        {
            'file':              'data/Sr20R21_phi0_alpha0_U20_loc05.dat_rotated.p',
            'case_name':         'Serrated, $z/\\lambda = 0.25$',
            'type':              'surface',
            'x_loc':             -1.0,
            'y_trust_min':       0.3,
            'Cf':                0.0015,
            'height_correction': z05_height_correction,
            'rotation':          z05_rotation_correction,
            'x_corr':            z05_x_correction,
            'marker':            '+',
            'color':             '0.8352941176470589, 0.3686274509803922, 0.0'
        },
        #{
        #    'file':              \
        #    'data/Sr20R21_phi0_alpha0_U20_loc10.dat.p',
        #    'case_name':         'Serrated, $z/\\lambda = 0.5$',
        #    'type':              'surface',
        #    'x_loc':             -1.0,
        #    'y_trust_min':       1.0,
        #    'Cf':                0.0015,
        #    'height_correction': z10_height_correction,
        #    'rotation':          z10_rotation_correction,
        #    'x_corr':            z10_x_correction,
        #    'marker':            'o',
        #    'color':             '0.8, 0.4745098039215686, 0.6549019607843137'
        #},
        {
            'file':              \
            'data/Sr20R21_phi0_alpha0_U20_loc10_even_newer.dat.p',
            'case_name':         'Serrated, $z/\\lambda = 0.5$',
            'type':              'surface',
            'x_loc':             -1.0,
            'y_trust_min':       1.0,
            'Cf':                0.0015,
            'height_correction': 0.5,
            'rotation':          0,
            'x_corr':            0,
            'marker':            'o',
            'color':             '0.8, 0.4745098039215686, 0.6549019607843137'
        },

        # x = 0 ###############################################################
        {
            'file':              'data/STE_phi0_alpha0_U20_loc00.dat_rotated.p',
            'case_name':         'Straight edge',
            'type':              'surface',
            'x_loc':             0,
            'y_trust_min':       1.0,
            'Cf':                0.00113,
            'height_correction': ste_height_correction,
            'rotation':          ste_rotation_correction,
            'x_corr':            ste_x_correction,
            'marker':            'x',
            'color':             '0.0, 0.4470588235294118, 0.6980392156862745'
        },
        {
            'file':              'data/Sr20R21_phi0_alpha0_U20_loc00.dat_rotated.p',
            'case_name':         'Serrated, $z/\\lambda = 0$',
            'type':              'surface',
            'x_loc':             0,
            'y_trust_min':       0,
            'Cf':                0.00113,
            'height_correction': z00_height_correction,
            'rotation':          z00_rotation_correction,
            'x_corr':            z00_x_correction,
            'marker':            '2',
            'color':             '0.0, 0.6196078431372549, 0.45098039215686275'
        },
        {
            'file':              'data/Sr20R21_phi0_alpha0_U20_loc05.dat_rotated.p',
            'case_name':         'Serrated, $z/\\lambda = 0.25$',
            'type':              'surface',
            'x_loc':             0,
            'y_trust_min':       0,
            'Cf':                0.00113,
            'height_correction': z05_height_correction,
            'rotation':          z05_rotation_correction,
            'x_corr':            z05_x_correction,
            'marker':            '+',
            'color':             '0.8352941176470589, 0.3686274509803922, 0.0'
        },
        {
            'file':              \
            'data/Sr20R21_phi0_alpha0_U20_loc10.dat.p',
            'case_name':         'Serrated, $z/\\lambda = 0.5$',
            'type':              'surface',
            'x_loc':             0,
            'y_trust_min':       1.0,
            'Cf':                0.00113,
            'height_correction': z10_height_correction,
            'rotation':          z10_rotation_correction,
            'x_corr':            z10_x_correction,
            'marker':            'o',
            'color':             '0.8, 0.4745098039215686, 0.6549019607843137'
        },

        # x = 9.0 #############################################################
        {
            'file':              'data/Sr20R21_phi0_alpha0_U20_loc00.dat_rotated.p',
            'case_name':         'Serrated, $z/\\lambda = 0$',
            'type':              'surface',
            'x_loc':             9.0,
            'y_trust_min':       0.3,
            'Cf':                0.0020,
            'height_correction': z00_height_correction,
            'rotation':          z00_rotation_correction,
            'x_corr':            z00_x_correction,
            'marker':            '2',
            'color':             '0.0, 0.6196078431372549, 0.45098039215686275'
        },
        {
            'file':              'data/Sr20R21_phi0_alpha0_U20_loc05.dat_rotated.p',
            'case_name':         'Serrated, $z/\\lambda = 0.25$',
            'type':              'surface',
            'x_loc':             9.0,
            'y_trust_min':       0.0,
            'Cf':                0.0020,
            'height_correction': z05_height_correction,
            'rotation':          z05_rotation_correction,
            'x_corr':            z05_x_correction,
            'marker':            '+',
            'color':             '0.8352941176470589, 0.3686274509803922, 0.0'
        },

        # x = 10.0 ############################################################
        {
            'file':              'data/Sr20R21_phi0_alpha0_U20_loc00.dat_rotated.p',
            'case_name':         'Serrated, $z/\\lambda = 0$',
            'type':              'surface',
            'x_loc':             10.0,
            'y_trust_min':       0.3,
            'Cf':                0.0020,
            'height_correction': z00_height_correction,
            'rotation':          z00_rotation_correction,
            'x_corr':            z00_x_correction,
            'marker':            '2',
            'color':             '0.0, 0.6196078431372549, 0.45098039215686275'
        },
        {
            'file':              'data/Sr20R21_phi0_alpha0_U20_loc05.dat_rotated.p',
            'case_name':         'Serrated, $z/\\lambda = 0.25$',
            'type':              'surface',
            'x_loc':             10.0,
            'y_trust_min':       0.0,
            'Cf':                0.0020,
            'height_correction': z05_height_correction,
            'rotation':          z05_rotation_correction,
            'x_corr':            z05_x_correction,
            'marker':            '+',
            'color':             '0.8352941176470589, 0.3686274509803922, 0.0'
        },


        # x = 19.0 ############################################################
        {
            'file':              'data/Sr20R21_phi0_alpha0_U20_loc00.dat_rotated.p',
            'case_name':         'Serrated, $z/\\lambda = 0$',
            'type':              'surface',
            'x_loc':             19.0,
            'y_trust_min':       0.0,
            'Cf':                0.0026,
            'height_correction': z00_height_correction,
            'rotation':          z00_rotation_correction,
            'x_corr':            z00_x_correction,
            'marker':            '2',
            'color':             '0.0, 0.6196078431372549, 0.45098039215686275'
        },
        {
            'file':              'data/Sr20R21_phi0_alpha0_U20_loc05.dat_rotated.p',
            'case_name':         'Serrated, $z/\\lambda = 0.25$',
            'type':              'surface',
            'x_loc':             19.0,
            'y_trust_min':       0.2,
            'Cf':                0.0026,
            'height_correction': z05_height_correction,
            'rotation':          z05_rotation_correction,
            'x_corr':            z05_x_correction,
            'marker':            '+',
            'color':             '0.8352941176470589, 0.3686274509803922, 0.0'
        },

        # x = 20.0 ############################################################
        {
            'file':              'data/Sr20R21_phi0_alpha0_U20_loc00.dat_rotated.p',
            'case_name':         'Serrated, $z/\\lambda = 0$',
            'type':              'surface',
            'x_loc':             20.0,
            'y_trust_min':       0.0,
            'Cf':                0.0026,
            'height_correction': z00_height_correction,
            'rotation':          z00_rotation_correction,
            'x_corr':            z00_x_correction,
            'marker':            '2',
            'color':             '0.0, 0.6196078431372549, 0.45098039215686275'
        },
        {
            'file':              'data/Sr20R21_phi0_alpha0_U20_loc05.dat_rotated.p',
            'case_name':         'Serrated, $z/\\lambda = 0.25$',
            'type':              'surface',
            'x_loc':             20.0,
            'y_trust_min':       0.0,
            'Cf':                0.0026,
            'height_correction': z05_height_correction,
            'rotation':          z05_rotation_correction,
            'x_corr':            z05_x_correction,
            'marker':            '+',
            'color':             '0.8352941176470589, 0.3686274509803922, 0.0'
        },


        # x = 39.0 ############################################################
        {
            'file':              'data/Sr20R21_phi0_alpha0_U20_loc00.dat_rotated.p',
            'case_name':         'Serrated, $z/\\lambda = 0$',
            'type':              'surface',
            'x_loc':             39.0,
            'y_trust_min':       0.4,
            'Cf':                0.0033,
            'height_correction': z00_height_correction,
            'rotation':          z00_rotation_correction,
            'x_corr':            z00_x_correction,
            'marker':            '2',
            'color':             '0.0, 0.6196078431372549, 0.45098039215686275'
        },

        # x = 40.0 ############################################################
        {
            'file':              'data/Sr20R21_phi0_alpha0_U20_loc00.dat_rotated.p',
            'case_name':         'Serrated, $z/\\lambda = 0$',
            'type':              'surface',
            'x_loc':             40.0,
            'y_trust_min':       0.4,
            'Cf':                0.0033,
            'height_correction': z00_height_correction,
            'rotation':          z00_rotation_correction,
            'x_corr':            z00_x_correction,
            'marker':            '2',
            'color':             '0.0, 0.6196078431372549, 0.45098039215686275'
        },

    ]

    for c in case_dicts:
        if c == case_dicts[0]:
            df = pd.DataFrame( data = c , index = [0])
        else:
            df = df.append( pd.DataFrame( data = c , index = [0]), 
                           ignore_index=True )
    return df
