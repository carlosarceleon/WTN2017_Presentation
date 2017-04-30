# ( x, y, rot )
import collections

correction_dict =  collections.OrderedDict()

correction_dict[ 'STE_a0_p0_U20_z00_tr'       ] = ( -1  , 0.6  , 180 + 1   , 0.5 )
correction_dict[ 'STE_a12_p0_U20_z00_tr'      ] = ( -1.5   , 1.3 + 1.2  , 2.0         , 0   )
correction_dict[ 'STE_a-12_p0_U20_z00_tr'     ] = ( 0   , -0.7 + 1  , 90        , 0   )

correction_dict[ 'Sr20R21_a0_p0_U20_z00_tr'   ] = ( 0   , -0.2 , 180 - 1.0 , 0.7 )
correction_dict[ 'Sr20R21_a0_p0_U20_z05_tr'   ] = ( 0   , 0.0  , 180 - 1.0 , 0.8 )
correction_dict[ 'Sr20R21_a0_p0_U20_z10_tr'   ] = ( 0   , -0.5 , 2.0       , 1.0 )

correction_dict[ 'Sr20R21_a12_p0_U20_z00_tr'  ] = ( 0   , -0.5 , -0.5      , 0   )
correction_dict[ 'Sr20R21_a12_p0_U20_z05_tr'  ] = ( 0   , 0    , 0.5       , 0   )
correction_dict[ 'Sr20R21_a12_p0_U20_z10_tr'  ] = ( 0   , 0    , 0         , 0   )

correction_dict[ 'Sr20R21_a12_p6_U20_z00_tr'  ] = ( 0   , -0.4 , 1         , 0   )
correction_dict[ 'Sr20R21_a12_p6_U20_z05_tr'  ] = ( -3  , 1 - 0.3   , -0.1      , 0.5 )
correction_dict[ 'Sr20R21_a12_p6_U20_z10_tr'  ] = ( -3  , -1.3 , 4.5       , 0   )

correction_dict[ 'Sr20R21_a-12_p0_U20_z00_tr' ] = ( 1.5 , -0.1 , 0         , 0   )
correction_dict[ 'Sr20R21_a-12_p0_U20_z05_tr' ] = ( -2  , -.5  , 180 + 0.8 , 0   )
correction_dict[ 'Sr20R21_a-12_p0_U20_z10_tr' ] = ( 0   , 0    , 1         , 0   )

correction_dict[ 'Sr20R21_a-12_p6_U20_z00_tr' ] = ( -1  , 0.6  , 180 - 0.8 , 0   )
correction_dict[ 'Sr20R21_a-12_p6_U20_z05_tr' ] = ( -10 , -2.4 - 1 , 180 - 1   , 0   )
correction_dict[ 'Sr20R21_a-12_p6_U20_z10_tr' ] = ( -10 , -2.5 , 180 + 1.  , 0   )

#correction_dict[ 'Slit20R21_a0_p0_U20_z00_tr'   ] = ( 0   , 0    , 180 - 1   , 0 )
#correction_dict[ 'Slit20R21_a0_p0_U20_z10_tr'   ] = ( 0   , -0.1 , 1.5       , 0 )
#correction_dict[ 'Slit20R21_a0_p0_U20_z05_tr'   ] = ( 0   , 0.6  , 180 - 1   , 0 )
#
#correction_dict[ 'Slit20R21_a12_p0_U20_z00_tr'  ] = ( -33 , 7.4  , 1.5       , 0 )
#correction_dict[ 'Slit20R21_a-12_p0_U20_z00_tr' ] = ( 0   , -0.3 , 0         , 0 )
#correction_dict[ 'Slit20R21_a-12_p0_U20_z10_tr' ] = ( 2.1 , -4.2 , 0         , 0 )
#correction_dict[ 'Slit20R21_a-12_p0_U20_z05_tr' ] = ( 2   , -4.6 , 180 - 1   , 0 )
