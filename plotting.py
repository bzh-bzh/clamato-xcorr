from matplotlib import pyplot as plt


def plot_preamble():
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='small')
    plt.rc('ytick', labelsize='small')
    # plt.rc('legend', fontsize=7.5)
    plt.rc('legend', fontsize=10)
    plt.rc('axes', titlesize='large')     # fontsize of the axes title
    plt.rc('axes', labelsize='large')    # fontsize of the x and y labels
    # plt.rc('text', usetex=True)
    # plt.rc('axes', unicode_minus=False)
    #matplotlib.rcParams['axes.unicode_minus'] = False
    
# plot_preamble()