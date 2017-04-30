def plot_gruber_strouhals():
    from matplotlib import pyplot as plt
    import pandas as pd
    import seaborn as sns
    from matplotlib import rc
    from numpy import array,argmin,unique

    rc('text',usetex=True)

    sns.set_context('paper')
    sns.set_style("whitegrid")
    sns.set(font='serif',font_scale=2.5,style='whitegrid')
    rc('font',family='serif', serif='cm10')

    strouhal_key = "$\\mathrm{St}_0=f_0\\delta/U_\\infty$"
    ratio_key    = "$\\lambda/h$"

    markers = [
            u'o', u'v', u'^', u'<', u'>', u'8', u's', u'p', u'*', 
        u'h', u'H', u'D', u'd'
    ]

    gruber_strouhals = [
        1.18,
        1.,
        1.35,
        1.35,
        1.4,
        1.45,
    ]
    gruber_ratios = [
        0.1,
        0.15,
        0.2,
        0.3,
        0.5,
        0.6,
    ]

    gruber_lengths = array([
        30,20,30,20,20,30,
    ])

    gruber_lambdas = array([
        1.5,1.5,3,3,5,9
    ])

    palette = sns.color_palette("cubehelix", 
                                n_colors=len(unique(gruber_lambdas))
                               )

    grouber_results_df = pd.DataFrame(
        columns=[ strouhal_key, ratio_key ]
    )

    for st,r,leng,lambd in zip(gruber_strouhals,gruber_ratios,
                         gruber_lengths,gruber_lambdas):
        grouber_results_df = grouber_results_df.append(
            {
                strouhal_key : st,
                ratio_key    : r,
                'length'     : leng,
                'lambda'     : lambd,
            },
            ignore_index=True
        )

    fig,ax = plt.subplots(1,1)
    plot_options = {
        's'      : 100,
        #'marker' : 's'
    }
    done_length = []
    done_lambda = []
    for index,row in grouber_results_df.iterrows():
        length_identifying_key = argmin( 
                abs(unique(gruber_lengths) - row.length) 
            ) 
        lambda_identifying_key = argmin( 
                abs(unique(gruber_lambdas) - row['lambda']) 
            ) 
        if row.length not in done_length:
            length_label = "$2h = {{{0}}}$ mm".format(row.length)
            done_length.append(row.length)
        else:
            length_label = '' 
        if row['lambda'] not in done_lambda:
            lambda_label = "$\\lambda = {{{0}}}$ mm".format(
                row['lambda']
            )
            done_lambda.append(row['lambda'])
        else:
            lambda_label = '' 
        ax.scatter(
            y      = row[strouhal_key],
            x      = row[ratio_key],
            c      = palette[lambda_identifying_key],
            #marker = markers[length_identifying_key],
            marker = markers[0],
            label  = lambda_label,
            **plot_options
        )
    print done_length
    ax.set_ylabel(strouhal_key)
    ax.set_xlabel(ratio_key)
    ax.legend(loc='lower right')
    plt.savefig('Gruber_Results.png',bbox_inches='tight')


def move_all_acoustic_data_to_local():
    import os
    import acoustic_functions as afunc

    destination = './AcousticData'
    acoustic_root = \
            '/media/carlos/6E34D2CD34D29783/2015-03_SerrationAcoustics/'
    acoustic_campaign = 'MarchData'

    acoustic_data_path = os.path.join(acoustic_root,acoustic_campaign)

    acoustic_cases = [f for f\
                      in os.listdir(acoustic_data_path)\
                      if os.path.isdir(
                          os.path.join(acoustic_data_path,f)
                      )]

    for ac in acoustic_cases:
        afunc.move_data_to_local(
            os.path.join(acoustic_data_path,ac),
            os.path.join(destination,"psd_"+ac+'.mat')
        )

def plot_interesting_cases( phi = 6, U = 35 ):
    from os import listdir
    from acoustic_functions import plot_spectra

    root = acoustic_data_folder

    cases_to_plot = [f for f in listdir(root)\
                     if 'STE' in f\
                     or 'Sr20R21' in f]
    cases_to_plot = [f for f in cases_to_plot\
                     if not "repitability" in f\
                     and not "Redo" in f\
                     and '{0}'.format( U ) in f\
                    ]
    cases_to_plot = [f for f in cases_to_plot\
                     if 'p{0}'.format( phi ) in f\
                     or "STE" in f
                    ]
    plot_spectra(
        root,
        cases_to_plot,
        third_octave = True,
        output       = './article_images/case{0}_spectra_p{1}.png'.\
        format( U, phi ),
        phi          = phi
    )

def plot_article_relative_cases(alpha = 0, phi = 0, article=True,
                               draw_crossover_points=True):
    import acoustic_functions as afunc
    from collections import OrderedDict

    cases = OrderedDict([
            ("psd_Sr20R21_a{0:02d}_p{1}_U30".format(alpha,phi), 
             "$U_\\infty = 30$ m/s".format(alpha)),
            ("psd_Sr20R21_a{0:02d}_p{1}_U35".format(alpha,phi), 
             "$U_\\infty = 35$ m/s".format(alpha)),
            ("psd_Sr20R21_a{0:02d}_p{1}_U40".format(alpha,phi), 
             "$U_\\infty = 40$ m/s".format(alpha)),
            ])

    relative_to = OrderedDict([
            ("psd_STE_a{0:02d}_U30".format(alpha), 
             "straight trailing edge, $\\alpha_g = 12^\circ$"),
            ("psd_STE_a{0:02d}_U35".format(alpha), 
             "straight trailing edge, $\\alpha_g = 12^\circ$"),
            ("psd_STE_a{0:02d}_U40".format(alpha), 
             "straight trailing edge, $\\alpha_g = 12^\circ$"),
            ])
    title = ""
    afunc.compare_cases_relative(
        acoustic_data_folder,
        cases=cases,
        relative_to=relative_to,
        plot_name="article_images/Relative_a{0}_p{1}.png"\
        .format(alpha,phi),
        title=title,
        article=article,
        draw_crossover_points=draw_crossover_points
    )


acoustic_data_folder = '/home/carlos/Documents/PhD/Articles/Article_3/'\
            + 'Scripts/acoustics_and_bl/AcousticDataRoberto/'
