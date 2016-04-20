import os
import data_sources as ds
import fit_grn as fg
import fused_reg as fr
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

#pct_fused is a list of what fraction of the genes are fused
#data_path is the directory where the synthetic data will live
#cv_folds is the number of cross validation folds
#returns array containing beta mean squared error relative to the unfused case
def fused_synthetic(pct_fused,data_path,cv_folds):    
    N_TF = 20
    N_G = 200
    pct_fused = list(np.linspace(0.3,1.0,10))
    reps = 1
    
    lamR = 2
    lamSs = list(np.linspace(0,5,10)) 
    lamP = 1.0
    aupr_array = np.zeros((len(pct_fused),len(lamSs)))

    out1 = data_path
    if not os.path.exists(out1):
        os.mkdir(out1)

    k = cv_folds

    for p in range(reps):
        for i, N in enumerate(pct_fused):
            out = os.path.join(out1,'dat_'+str(N))
            ds.write_fake_data1(N1=10, N2=200, out_dir = out, tfg_count1=(N_TF, N_G), tfg_count2 = (N_TF, N_G), measure_noise1 = 0.1, measure_noise2 = 0.1, sparse=0.75, fuse_std = 0, pct_fused=N, orth_falsepos = 0)        
            lamP = 1.0 
            seed = 10
            for j, lamS in enumerate(lamSs):
                (errd1, errd2) = fg.cv_model_m(out, lamP, lamR, lamS, k, solver='solve_ortho_direct', reverse=True, cv_both=(True,True), exclude_tfs=True, pct_priors=0, seed=seed, verbose=False)
                print errd1['B_mse'].mean()
                aupr_array[i,j]+= errd1['B_mse'].mean()

    aupr_array /= reps
    ar2 = np.zeros((aupr_array.shape[0],aupr_array.shape[1]))
    for i in range(len(aupr_array)):
        ar2[i,:] = aupr_array[i,:]/aupr_array[i,0]

    return ar2

#savef is file to save results in
#this function takes list of lamS values, cv_folds, and savef, and plots the precision recall curve 
def fused_bacterial(lamSs, cv_folds, savef):
    lamP = 1.0
    lamR = 0.5
    k = cv_folds
    metric = 'prc'
    normed = False
    cv_both = (True, False)
    roc_species = 0
    orgs = ['B_subtilis','B_anthracis']
    unfused = False
    orth_file = ['orth']
    pct_priors = 0
    test_all = 'part'
    lamS_opt = None
    out = os.path.join('data','bacteria_standard')
    errdls = plot_roc(out, lamP, lamR, lamSs, k, metric, savef,normed,scad, cv_both, roc_species, orgs, lamS_opt, unfused, orth_file, pct_priors,test_all)

    def plot_roc(out, lamP=1.0, lamR=5, lamSs=[0,1], k=20, metric='roc',savef=None,normed=False,scad=False, cv_both=(True,False), roc_species=0, orgs=None, lamS_opt = None, unfused = False, orth_file=['orth'], pct_priors=0, test_all='part'):

        seed = np.random.randn()
        scad = False
        settings = fr.get_settings()
        solver = 'solve_ortho_direct'
        xl = 'recall'
        yl = 'precision'
        summary_measure = 'aupr'
            chancek = 'chance'

        all_roc_curves = []
        import pickle
        if savef != None:
            savef = os.path.join('saves',savef)
        loaded = False
        if savef != None and os.path.exists(savef):
            with file(savef) as f:
                loaded = True
                errdls = pickle.load(f)
        else:
            errdls = []

        else:
            for i, lamS in enumerate(lamSs):
                errd = fg.cv_model_m(out, lamP=lamP, lamR=lamR, lamS=lamS, k=k, solver = solver, settings = settings, reverse = True, cv_both = cv_both, exclude_tfs=False, pct_priors=pct_priors,seed = seed, orgs = orgs, lamS_opt = None, orth_file = orth_file, test_all = test_all)
                errdf = errd[roc_species]
                errdls.append(errdf)
            if lamS_opt != None:
                errdf = fg.cv_model_m(out, lamP=lamP, lamR=lamR, lamS=lamS, k=k, solver = solver, settings = settings, reverse = True, cv_both = cv_both, exclude_tfs=False, pct_priors=pct_priors,seed = seed, orgs = orgs, lamS_opt = lamS_opt, orth_file = orth_file, test_all = test_all)[roc_species]
                errdls.append(errdf)
                    
        if savef != None and not loaded:
            with file(savef, 'w') as f:
                pickle.dump(errdls, f)

        chance_rates = []

        toplots = []
        toplots.extend(map(str, lamSs))
        for i, lamS in enumerate(lamSs):
            #print lamS
            #print i
            errd = errdls[i]      
            rocs = errd[metric]
            print (lamS, errd[summary_measure].mean())
            all_roc_curves.append(rocs)
            chance_rates.append(errd[chancek])
        if lamS_opt != None:
            print (lamS_opt, errd_opt[summary_measure].mean())
            errd_opt = errdls[-1]
            rocs_opt = errd_opt[metric]        
            all_roc_curves.append(rocs_opt)
            chance_rates.append(errd_opt[chancek])
            toplots.append(str(lamS_opt))
        linedesc = pd.Series(toplots, name='method')

        pss = []
        for roc_curve_group in all_roc_curves:
            (rs, ps, ts) = fg.pool_roc(roc_curve_group, all_roc_curves, max_x = 1000)
            pss.append(ps)
        
        if normed:
            to_plot = np.dstack(pss) / np.repeat(pss[0][:,:,None], len(pss), axis=2)
        else:
           to_plot = np.dstack(pss)
        xs = pd.Series(rs, name=xl)
        
        sns.tsplot(to_plot, time=xs, condition=linedesc, value=yl)
        plt.hold(True)

        chance_x = [min(xs), max(xs)]
        chance_y = [np.mean(chance_rates), np.mean(chance_rates)]
        #plt.plot(chance_x, chance_y,'--k')
        plt.show()
        return errdls