import matplotlib
matplotlib.use('Agg')

import os

import peer
import scipy as SP
import pylab as PL
import pandas as pd


result_dir = 'PEER_results/'


def plot_Alpha(Alpha, color="blue"):
    PL.plot(1.0 / Alpha,lw=4, color=color)
    min_a,max_a = (1.0/Alpha).min(), (1.0/Alpha).max()
    PL.ylim(min_a - 0.1*(max_a - min_a), max_a + 0.1*(max_a - min_a))
    PL.xlabel("Factors")
    PL.ylabel("Factor relevance")


def compute_peer_factors(tissue_name='Blood Vessel'):
    tissue_rpkm_file = '../Expression_by_Subissue/%s.rpkm' % tissue_name
    K = 15
    peer_result_file = result_dir + 'peer_factors_%s_%s' % (tissue_name, K)
    if os.path.exists(peer_result_file):
        return
#    y = SP.loadtxt(tissue_rpkm_file, delimiter="\t")
#    print(y)
    df = pd.read_csv(tissue_rpkm_file, delimiter='\t', header=1)
    df = df.drop(columns=['Name'])
    if len(df.columns) < K + 1:
        return
    df.set_index('Description', inplace=True)
    df = df.transpose()
    return compute_and_store_peer_factors(df, peer_result_file, K)


def compute_and_store_peer_factors(df, peer_result_file, K=15):
#    print(df)
#    print(df.values)
    # df:
    #       Gene1   Gene2
    # Ind1  0.4     6
    # Ind2  14      0.5
    Nmax_iterations = 100
    Nmax_iterations = 3
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
#    plot_Alpha(Alpha)
#    PL.savefig("demo_simple.pdf")
    print "Plotted factor relevance"

    ids = list(df.index)
    with open(peer_result_file, 'w') as outfile:
        for i in range(len(X)):
            individual_id = ids[i]
            individual_factors = '\t'.join([str(e) for e in X[i]])
            outfile.write('%s\t%s\n' % (individual_id, individual_factors))
    return X


def get_simple_model_object(K=10, Nmax_iterations=100, expr_file="data/expression.csv"):
    #y = SP.loadtxt(expr_file,delimiter=",")
    df = pd.read_csv(expr_file, delimiter='\t', header=1)
    df = df.drop(columns=['Name'])
    df.set_index('Description', inplace=True)
    df = df.transpose()
    y = df.values
    model = peer.PEER()
    # set data and parameters
    model.setNk(K) #number of factor for learning
    model.setPhenoMean(y) # data for inference
    model.setNmax_iterations(Nmax_iterations)
    return model


def n_factors_demo():
    print "Comparing different numbers of factors in inference (n=2,4,6,10)."
    colors = ("yellow","red","green","blue","black")

    # plot factor weight variances for a large set of maximum number K of factors, and see if K has an effect
    alpha_10 = None
    n_factors = [10,6,4,2]
    for i,k in enumerate(n_factors):
        model = get_simple_model_object(K=k) # see simple_unsupervised_demo for how it is constructed
        model.update()
        plot_Alpha(model.getAlpha(),colors[i])
        if k == 10: alpha_10 = model.getAlpha()
    min_a,max_a = (1./alpha_10).min(), (1./alpha_10).max()
    PL.ylim(min_a - 0.1*(max_a - min_a), max_a + 0.1*(max_a - min_a))
    PL.legend(['K=%d'%f for f in n_factors])
    # expected result - as soon as there are at least 5 factors, the factor inference will not change any more
    PL.savefig("demo_factors.pdf")
    print "Plotting factor relevances"
    PL.show()


def eps_prior_demo():
    print "Comparing different noise priors to see the effect on how aggressively PEER explains variability."

    # plot factor weight variances for a large set of maximum number K of factors, and see if K has an effect
    for pa in (0.0001, 0.1, 1000):
        for pb in (0.1,10,1000):
            model = get_simple_model_object() # simple object using default simulated dataset; see simple_unsupervised_demo for how it is constructed
            model.setPriorEps(pa,pb);
            model.update()
            print "Eps pa=%.4f pb=%.4f mean(residuals^2)=%.4f"%(pa, pb, SP.mean(model.getResiduals()**2))


if __name__ == '__main__':
    import glob
    expression_files = glob.glob('../Expression_By_Subtissue/*')
    tissue_names = [ef.split('/')[-1][:-5] for ef in expression_files]
    for tissue in tissue_names:
        compute_peer_factors(tissue_name=tissue)

