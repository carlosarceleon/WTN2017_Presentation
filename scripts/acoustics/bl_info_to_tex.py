def bl_info_to_tex():
    from pandas import read_pickle

    bl_pickle = 'outputs/BL_Data_Info.p'

    bldf = read_pickle( bl_pickle )
    bldf = bldf[ bldf.Test_section == 'open' ]

    delta_pivot = bldf.to

