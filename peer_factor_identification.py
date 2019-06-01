import peer
import scipy as SP
import pylab as PL
import pandas as pd

def compute_peer_factors(tissue_rpkm_file='../Expression_by_Tissue/Blood Vessel.rpkm.bck'):
#    y = SP.loadtxt(tissue_rpkm_file, delimiter="\t")
    df = pd.read_csv(tissue_rpkm_file, delimiter='\t', header=1)
    df = df.drop(columns=['Name'])
#    print(y)
    df.set_index('Description', inplace=True)
    df = df.transpose()
#    print(df)
#    print(df.values)
    K = 15
    Nmax_iterations = 100
    model = peer.PEER()

    model.setNk(K) #number of factor for learning
    model.setPhenoMean(df.values) # data for inference
    # set priors (these are the default settings of PEER)
    model.setPriorAlpha(0.001,0.1);
    model.setPriorEps(0.1,10.);
    model.setNmax_iterations(Nmax_iterations)
    # perform inference
    model.update()

    #investigate results
    #factors:
    X = model.getX()
    #weights:
    W = model.getW()
    #ARD parameters
    Alpha = model.getAlpha()
    #get corrected dataset:
    Yc = model.getResiduals()

    # plot variance of factors - in this case, we expect a natural elbow where there are 5 active factors, as 5 were simulated
    plot_Alpha(Alpha)
    PL.savefig("demo_simple.pdf")
    print "Plotted factor relevance"


compute_peer_factors()

