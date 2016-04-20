import numpy as np
import fused_reg as fl
import data_sources as ds
import random
from sklearn.metrics import roc_curve, auc, precision_recall_curve
import os
import collections
import matplotlib.pyplot as plt
#SECTION: ------------------UTILITY FUNCTIONS-----------------
#returns a mean-squared-error evaluator function, which has signature
# f(Bs) -> mean mse across betas, evaluated on Xs and Ys
#exclude_tfs: don't evaluate on transcription factors
def get_mse_evaluator(Xs_te, Ys_te, exclude_tfs=False):
    def inner(Bs):
        mse = 0
        for si in range(len(Bs)):
            mse += prediction_error(Xs_te[si], Bs[si], Ys_te[si], 'mse', exclude_tfs=exclude_tfs)
        mse /= len(Bs)
        return mse
    return inner

#returns a R2 evaluator function, which has signature
# f(Bs) -> mean R2 across betas, evaluated on Xs and Ys
#exclude_tfs: don't evaluate on transcription factors
def get_R2_evaluator(Xs_te, Ys_te, exclude_tfs=False):
    def inner(Bs):
        mse = 0
        for si in range(len(Bs)):
            mse += prediction_error(Xs_te[si], Bs[si], Ys_te[si], 'R2', exclude_tfs=exclude_tfs)
        mse /= len(Bs)
        return mse
    return inner

#returns a R2 evaluator function, which has signature
# f(Bs) -> mean R2 across betas, evaluated on Xs and Ys
#exclude_tfs: don't evaluate on transcription factors
def get_aupr_evaluator(Xs_te, Ys_te, genes, tfs, priors, tr_priors, exclude_tfs=False, test_all='part'):
    
    def inner(Bs):
        aupr_acc = 0
        for si in range(len(Bs)):
            Xsi = Xs_te[si]
            Ysi = Ys_te[si]
            Bsi = Bs[si]
            S = rescale_betas(Xsi, Ysi, Bsi)
        
            (aupr, prc) = eval_network_pr(S, genes[si], tfs[si], priors[si], tr_priors=tr_priors[si], exclude_tfs=exclude_tfs, constraints = None, test_all=test_all)        
            
            if not np.isnan(aupr):
                aupr_acc += aupr
        return aupr_acc
    return inner

#hey! what does this do?
def get_rank(scores, labels, coords):
    score_array = np.array(scores)
    label_array = np.array(labels)
    coord_array = np.array(coords)

    temp = abs(score_array).argsort()
    sorted_scores = score_array[temp]
    sorted_labels = label_array[temp].tolist()
    sorted_coords = coord_array[temp].tolist()

    rank = 0
    ranks = []
    tied_ranks = []
    for i in range(len(sorted_scores)):
        rank += 1
        if i < len(sorted_scores)-1:
            if sorted_scores[i] == sorted_scores[i+1]:
                tied_ranks.append(rank)
                for j in range(i+1, len(sorted_scores)):
                    while sorted_scores[i] == sorted_scores[j]:
                        rank += 1
                        tied_ranks.append(rank)
                    to_append = (j-i+1) * [np.mean(tied_ranks)]
                    ranks += to_append
            else:
                ranks.append(rank)
        if i == len(sorted_scores)-1:
            ranks.append(rank)

    return (ranks, sorted_labels, sorted_coords)

def to_rank(arr):
    arr = np.abs(arr)
    shape = arr.shape
    arr = arr.ravel()
    idx = np.argsort(arr)
    ranks = np.zeros(arr.shape[0])

    ranks[idx] = np.arange(arr.shape[0])
    return ranks.reshape(shape)

#take lists and averages ranks across coords, where each coords is a list of (tf, g) pairs corresponding to ranks
#labels_combined should be a list of 0s and 1s, unless there is disagreement among gold standards
def rank_combine(rankslist, labelslist, coordslist):
    rankslist = map(to_rank, rankslist)
    rank_dict = {}
    for i in range(len(coordslist)):
        for count, [tf,g] in enumerate(coordslist[i]):
            if (tf,g) in rank_dict:
                rank_dict[(tf,g)][0] += rankslist[i][count]
                rank_dict[(tf,g)][1] += labelslist[i][count]
                rank_dict[(tf,g)][2] += 1
            else:
                rank_dict[(tf,g)] = [rankslist[i][count], labelslist[i][count], 1]

    ranks_combined = []
    labels_combined = []
    coords_combined = []

    for (tf, g) in rank_dict:
        denom = rank_dict[(tf, g)][2]
        (rank, label, one) = map(lambda x: float(x)/denom, rank_dict[(tf, g)])
        ranks_combined.append(rank)
        labels_combined.append(label)
        coords_combined.append((tf,g))
    return (ranks_combined, labels_combined, coords_combined)

def rank_combine2(Ss, tfs, genes):
    def inds(x):
        idx = abs(x).argsort()
        y = np.empty(x.shape)
        y[idx] = np.arange(x.shape[0])
        return y

    Ss_ranks = map(lambda S: inds(S.ravel()).reshape(S.shape), Ss)
    coeff_to_ranks = collections.defaultdict(lambda: [])
    for i, S in enumerate(Ss_ranks):
        for row, tf in enumerate(tfs[i]):
            for col, gene in enumerate(genes[i]):
                coeff_to_ranks[(tf, gene)].append(S[row, col])
    
    tfs_u = list(reduce(lambda x,y: x.intersection(y), map(lambda z: set(z), tfs)))
    genes_u = list(reduce(lambda x,y: x.intersection(y), map(lambda z: set(z), genes)))
    rc_S = np.zeros((len(tfs_u), len(genes_u)))
    for row, tf in enumerate(tfs_u):
        for col, gene in enumerate(genes_u):
            rc_S[row, col] = np.mean(coeff_to_ranks[(tf, gene)])
    return (rc_S, tfs_u, genes_u)

def rank_combine2(Ss, tfs, genes):
    def inds(x):
        idx = x.argsort()
        y = np.empty(x.shape)
        y[idx] = np.arange(x.shape[0])
        return y

    Ss_ranks = map(lambda S: inds(S.ravel()).reshape(S.shape), Ss)
    coeff_to_ranks = collections.defaultdict(lambda: [])
    for i, S in enumerate(Ss_ranks):
        for row, tf in enumerate(tfs[i]):
            for col, gene in enumerate(genes[i]):
                coeff_to_ranks[(tf, gene)].append(S[row, col])
    
    tfs_u = list(reduce(lambda x,y: x.intersection(y), map(lambda z: set(z), tfs)))
    genes_u = list(reduce(lambda x,y: x.intersection(y), map(lambda z: set(z), genes)))
    rc_S = np.zeros((len(tfs_u), len(genes_u)))
    for row, tf in enumerate(tfs_u):
        for col, gene in enumerate(genes_u):
            rc_S[row, col] = np.mean(coeff_to_ranks[(tf, gene)])
    return (rc_S, tfs_u, genes_u)
#SECTION: ------------------FOR RUNNING BACTERIAL DATA

def fit_model(data_fn, lamP, lamR, lamS, solver='solve_ortho_direct', settings = None, orth_file ='orth'):
    ds1 = ds.standard_source(data_fn,0)
    ds2 = ds.standard_source(data_fn,1)
    orth_fn = os.path.join(data_fn, orth_file)

    organisms = [ds1.name, ds2.name]
    orth = ds.load_orth(orth_fn, organisms)

    
    (priors1, signs1) = ds1.get_priors()
    (priors2, signs2) = ds2.get_priors()

    (e1_tr, t1_tr, genes1, tfs1) = ds1.load_data()
    (e2_tr, t2_tr, genes2, tfs2) = ds2.load_data()

        # jam things together
    Xs = [t1_tr, t2_tr]
    Ys = [e1_tr, e2_tr]
    
    genes = [genes1, genes2]
    tfs = [tfs1, tfs2]
    priors = priors1 + priors2

    if solver == 'solve_ortho_direct':
        Bs = fl.solve_ortho_direct(organisms, genes, tfs, Xs, Ys, orth, priors, lamP, lamR, lamS, settings = settings)
    if solver == 'solve_ortho_direct_scad':
        Bs = fl.solve_ortho_direct_scad(organisms, genes, tfs, Xs, Ys, orth, priors, lamP, lamR, lamS, settings = settings)
    if solver == 'solve_ortho_ref':
        Bs = fl.solve_ortho_ref(organisms, genes, tfs, Xs, Ys, orth, priors, lamP, lamR, lamS, settings = settings)
    if solver == 'iter_solve':
        #solve solution paths then return the last value
        Bs = fl.solve_ortho_iter(organisms, genes, tfs, Xs, Ys, orth, priors, lamP, lamR, lamS, settings = settings)
        
    return Bs

#solves networks separately, then rank combines to get network
def cv_unfused(data_fn, lamP, lamR, k, solver='solve_ortho_direct', settings=None, reverse=False, cv_both=(True,True), exclude_tfs=True, pct_priors=0, seed=None, verbose=False, orth_file=['orth'], orgs=None, lamS_opt=None, test_all='part', use_TFA=False):
    lamS = 0
    if seed != None:
        random.seed(seed)

    metrics1 = ['aupr','auc','chance']
    num_species = 0
    dss = []
    err_dict = {m : np.zeros((k, 1)) for m in metrics1}
    err_dict['params'] = (lamP, lamR)
    all_orgs = []
    organisms = []

    if orgs == None:
        while os.path.isfile(os.path.join(data_fn, 'expression%d' % (num_species+1))):
            dsi = ds.standard_source(data_fn,num_species, use_TFA=use_TFA)
            dss.append(dsi)
            organisms.append(dsi.name)
            all_orgs.append(dsi.name)
            err_dict['params'] = (lamP, lamR)
            num_species += 1

    else:
        while os.path.isfile(os.path.join(data_fn, 'expression%d' % (num_species+1))):
            dsi = ds.standard_source(data_fn,num_species, use_TFA=use_TFA)
            all_orgs.append(dsi.name)
            if dsi.name in orgs:
                dss.append(dsi)
                organisms.append(dsi.name)
                err_dict['params'] = (lamP, lamR)
            num_species +=1
        num_species-=1

    metrics2 = ['prc','roc']
    for metric in metrics2:
        err_dict[metric] = map(lambda x: [], range(k))
    metrics = metrics1 + metrics2
    (constraints, marks, orth) = ds.load_constraints(data_fn, orgs=organisms)

    if len(orth_file) == 1:
        orth_fn = os.path.join(data_fn, orth_file[0])
        orth = ds.load_orth(orth_fn, all_orgs, organisms)
    else:
        orth = []
        for i in range(len(orth_file)):
            orth_fn = os.path.join(data_fn, orth_file[i])
            orth += ds.load_orth(orth_fn, all_orgs, organisms)
    
    folds = map((lambda x: x.partition_data(k)), dss)
    all_priors = map((lambda x: x.get_priors()[0]), dss)

    #helper to return all but ith entry of list x    
    excl = lambda x,i: x[0:i]+x[(i+1):] 
    
    #helper function for dividing priors
    def r_partition(x, t):
        inds = np.arange(len(x))
        random.shuffle(inds)
        p1 = map(lambda i: x[i], inds[0:t])
        p2 = map(lambda i: x[i], inds[t:])
        return (p1, p2)

    allpriors = map(lambda x: r_partition(x, int(pct_priors*len(x))), all_priors)
    priors_tr = map(lambda x: x[0], allpriors)
    priorstr = []
    for priorslist in priors_tr:
        priorstr.extend(priorslist)
    priors_te = map(lambda x: x[1], allpriors)
    priorste = []
    for priorslist in priors_te:
        priorste.extend(priorslist)
    priorstest = list(set(priorste))

    f_te = [None]*num_species
    f_tr = [None]*num_species
    genes = [None]*num_species
    tfs = [None]*num_species
    Xs = [None]*num_species
    Ys = [None]*num_species
    Xs_te = [None]*num_species
    Ys_te = [None]*num_species
    Bs = [None]*num_species

        #downsamples the ROC, or PRC, curve
    def downsample_roc(roc):
        if roc[0] == None:
            return roc
        (roc_f, roc_h, roc_t) = pool_roc([roc], max_x=10000)
            
            #hrm?
        return (roc_f[:], roc_h[0,:], roc_t[0,:])

    for fold in range(k):
        if verbose:
            print 'working on %d'  % fold
        #get conditions for current cross-validation fold
        for si in range(num_species):
            if cv_both[si] and k > 1:
                f_te[si] = folds[si][fold]
                f_tr[si] = np.hstack(excl(folds[si], fold))    
            if reverse:
                tmp = f_tr[si]
                f_tr[si] = f_te[si]
                f_te[si] = tmp
            (e_tr, t_tr, genes_si, tfs_si) = dss[si].load_data(f_tr[si])
            (e_te, t_te, genes_si, tfs_si) = dss[si].load_data(f_te[si])

            Xs[si] = t_tr
            Xs_te[si] = t_te

            Ys[si] = e_tr
            Ys_te[si] = e_te

            genes[si] = genes_si
            tfs[si] = tfs_si

        priors_tr_fl = reduce(lambda x,y: x+y, priors_tr)

        #solve the model
        for si in range(num_species):
            if solver == 'solve_ortho_direct':
                B = fl.solve_ortho_direct(organisms, genes, tfs, Xs, Ys, orth, priors_tr_fl, lamP, lamR, lamS, lamS_opt, settings = settings)[si]
            if solver == 'solve_ortho_direct_scad':
                B = fl.solve_ortho_direct_scad(organisms[si], genes[si], tfs[si], Xs[si], Ys[si], orth, priors_tr_fl[si], lamP, lamR, lamS, settings = settings)
            if solver == 'solve_ortho_direct_scad_plot':
                B = fl.solve_ortho_direct_scad_plot(data_fn, organisms[si], genes[si], tfs[si], Xs[si], Ys[si], orth, priors_tr_fl[si], lamP, lamR, lamS, settings = settings)
            if solver == 'solve_ortho_direct_mcp':
                B = fl.solve_ortho_direct_mcp(organisms[si], genes[si], tfs[si], Xs[si], Ys[si], orth, priors_tr_fl[si], lamP, lamR, lamS, settings = settings)
            if solver == 'solve_ortho_direct_em':
                B = fl.solve_ortho_direct_em(organisms[si], genes[si], tfs[si], Xs[si], Ys[si], orth, priors_tr_fl[si], lamP, lamR, lamS, settings = settings)
            Bs[si] = (B)

        allranks = []
        alllabels = []
        allcoords = []
        Ss = []
        #return (Xs, Ys, Bs, genes, tfs, all_priors, priors_tr)
        for si in range(num_species):
            Xsi = Xs[si]    
            Ysi = Ys[si]
            Bsi = Bs[si]

            if True: #always rescale the beta matrix now
                S = rescale_betas(Xsi, Ysi, Bsi)
            else:
                S = Bsi
            Ss.append(S)

        (rc_S, tfs_u, genes_u) = rank_combine2(Ss, tfs, genes)
                  
        (aupr, prc) = eval_network_pr(rc_S, genes_u, tfs_u, priorstest, tr_priors=priorstr, exclude_tfs=exclude_tfs, constraints = None, test_all=test_all)

        err_dict['aupr'][fold,0] = aupr
        err_dict['prc'][fold] = prc    
        
        (auc, roc) = eval_network_roc(rc_S, genes_u, tfs_u, priorstest, exclude_tfs=exclude_tfs, constraints = None, test_all=test_all)

        err_dict['auc'][fold,0] = auc
        err_dict['roc'][fold] = downsample_roc(roc)
        
        #chance precision and chance constrained interactions precision
        chance = compute_chance_precision(rc_S, genes_u, tfs_u, priorstest, tr_priors=priorstr, exclude_tfs=exclude_tfs, constraints = None, sub = si)
        err_dict['chance'][fold, 0] = chance
        
        #top 100 interactions
        err_dict['top_100'] = top_k_interactions(rc_S, genes_u, tfs_u, priorstest, org=organisms[si], k=np.inf)

    return err_dict



#this is the master cross-validator!
#if orgs=[], then use all the organisms in folder; otherwise use the organisms listed in orgs
#orth_file can be a list of orth file names. if it is a list then it uses all the orth files. 
def cv_model_m(data_fn, lamP, lamR, lamS, k, solver='solve_ortho_direct',settings = None, reverse=False, cv_both=(True,True), exclude_tfs=True, pct_priors=0, seed=None, verbose=False, orth_file=['orth'], orgs=None, lamS_opt = None, test_all = 'part', use_TFA=False):
    if seed != None:
        random.seed(seed)

    metrics1 = ['mse','R2','aupr','auc','corr', 'auc_con','aupr_con', 'auc_noncon', 'aupr_noncon', 'chance', 'chance_con', 'B_mse','top_100']
    num_species = 0
    dss = []
    all_orgs = []
    err_dicts = []
    organisms = []

    if orgs == None:
        while os.path.isfile(os.path.join(data_fn, 'expression%d' % (num_species+1))):
            dsi = ds.standard_source(data_fn,num_species, use_TFA=use_TFA)
            dss.append(dsi)
            organisms.append(dsi.name)
            all_orgs.append(dsi.name)
            err_dicts.append({m : np.zeros((k, 1)) for m in metrics1})
            err_dicts[num_species]['params'] = (lamP, lamR, lamS, settings)
            num_species += 1

    else:
        cand_species = 1
        while os.path.isfile(os.path.join(data_fn, 'expression%d' % (cand_species))):
            dsi = ds.standard_source(data_fn,cand_species - 1, use_TFA=use_TFA)
            all_orgs.append(dsi.name)            
            if dsi.name in orgs:
                dss.append(dsi)
                organisms.append(dsi.name)
                err_dicts.append({m : np.zeros((k, 1)) for m in metrics1})
                err_dicts[len(err_dicts)-1]['params'] = (lamP, lamR, lamS, settings)
                num_species +=1
            cand_species += 1

    #set up containers for results
    #prc and roc are special (not individual numbers)
    metrics2 = ['prc','roc', 'prc_con','roc_con', 'prc_noncon', 'roc_noncon']
    for err_dict in err_dicts:
        for metric in metrics2:
            err_dict[metric] = map(lambda x: [], range(k))
    metrics = metrics1 + metrics2
    (constraints, marks, orth) = ds.load_constraints(data_fn, orth_f='orth', orgs=organisms)

    if len(orth_file) == 1:
        orth_fn = os.path.join(data_fn, orth_file[0])
        orth = ds.load_orth(orth_fn, all_orgs, organisms)
    else:
        orth = []
        for i in range(len(orth_file)):
            orth_fn = os.path.join(data_fn, orth_file[i])
            orth += ds.load_orth(orth_fn, all_orgs, organisms)

    folds = map((lambda x: x.partition_data(k)), dss)
    all_priors = map((lambda x: x.get_priors()[0]), dss)
    
    #helper to return all but ith entry of list x    
    excl = lambda x,i: x[0:i]+x[(i+1):] 
    
    #helper function for dividing priors
    def r_partition(x, t):
        inds = np.arange(len(x))
        random.shuffle(inds)
        p1 = map(lambda i: x[i], inds[0:t])
        p2 = map(lambda i: x[i], inds[t:])
        return (p1, p2)

    allpriors = map(lambda x: r_partition(x, int(pct_priors*len(x))), all_priors)
    #print (len(allpriors))
    priors_tr = map(lambda x: x[0], allpriors)

    if test_all == 'all':
        priors_te = map(lambda x: x[0]+x[1], allpriors)

    elif test_all == 'part':
        priors_te = map(lambda x: x[1], allpriors)

    #right now only used for th17, which has a different gold standard than priors
    elif test_all == 'gold':    
        all_gold = map((lambda x: x.get_gold()[0]), dss)
        priors_te = all_gold

    else:
        print 'test_all must be all, part or gold'
    #print (len(priors_tr[0]), len(priors_tr[1]))
    #print (len(priors_te[0]), len(priors_te[1]))
    
    f_te = [None]*num_species
    f_tr = [None]*num_species
    genes = [None]*num_species
    tfs = [None]*num_species
    Xs = [None]*num_species
    Ys = [None]*num_species
    Xs_te = [None]*num_species
    Ys_te = [None]*num_species

        #downsamples the ROC, or PRC, curve
    def downsample_roc(roc):
        if roc[0] == None:
            return roc
        (roc_f, roc_h, roc_t) = pool_roc([roc], max_x=10000)
            
            #hrm?
        return (roc_f[:], roc_h[0,:], roc_t[0,:])

    for fold in range(k):
        if verbose:
            print 'working on %d' % fold
        #get conditions for current cross-validation fold
        for si in range(num_species):
            if cv_both[si] and k > 1:
                f_te[si] = folds[si][fold]
                f_tr[si] = np.hstack(excl(folds[si], fold))    
            if reverse:
                tmp = f_tr[si]
                f_tr[si] = f_te[si]
                f_te[si] = tmp
            (e_tr, t_tr, genes_si, tfs_si) = dss[si].load_data(f_tr[si])
            (e_te, t_te, genes_si, tfs_si) = dss[si].load_data(f_te[si])

            Xs[si] = t_tr
            Xs_te[si] = t_te

            Ys[si] = e_tr
            Ys_te[si] = e_te

            genes[si] = genes_si
            tfs[si] = tfs_si

        
        priors_tr_fl = reduce(lambda x,y: x+y, priors_tr)

        #solve the model
        if solver == 'solve_ortho_direct':
            Bs = fl.solve_ortho_direct(organisms, genes, tfs, Xs, Ys, orth, priors_tr_fl, lamP, lamR, lamS, lamS_opt, settings = settings)
        if solver == 'solve_ortho_direct_scad':
            Bs = fl.solve_ortho_direct_scad(organisms, genes, tfs, Xs, Ys, orth, priors_tr_fl, lamP, lamR, lamS, settings = settings)
        if solver == 'solve_ortho_direct_scad_plot':
            Bs = fl.solve_ortho_direct_scad_plot(data_fn, organisms, genes, tfs, Xs, Ys, orth, priors_tr_fl, lamP, lamR, lamS, settings = settings)
        if solver == 'solve_ortho_direct_mcp':
            Bs = fl.solve_ortho_direct_mcp(organisms, genes, tfs, Xs, Ys, orth, priors_tr_fl, lamP, lamR, lamS, settings = settings)
        if solver == 'solve_ortho_direct_em':
            Bs = fl.solve_ortho_direct_em(organisms, genes, tfs, Xs, Ys, orth, priors_tr_fl, lamP, lamR, lamS, settings = settings)
        if solver == 'iter_solve':
        #solve solution paths then return the last value
            if settings['iter_eval']:
                #iter_eval = get_R2_evaluator(Xs_te, Ys_te)
                #iter_eval = get_mse_evaluator(Xs_te, Ys_te)
                iter_eval = get_aupr_evaluator(Xs_te, Ys_te, genes, tfs, priors_te, priors_tr)
                settings['iter_eval'] = iter_eval
                Bs = fl.solve_ortho_iter(organisms, genes, tfs, Xs, Ys, orth, priors_tr_fl, lamP, lamR, lamS, settings = settings)
            
         #evaluate a bunch of metrics
        (corr, fused_coeffs) = fused_coeff_corr(organisms, genes, tfs, orth, Bs)
            
        for si in range(num_species):
            print si
            #correlation of fused coefficients
            err_dicts[si]['corr'][fold,0] = corr

            #mse
            mse = prediction_error(Xs_te[si], Bs[si], Ys_te[si], 'mse', exclude_tfs=exclude_tfs)
            err_dicts[si]['mse'][fold, 0] = mse

            #R2
            R2 = prediction_error(Xs_te[si], Bs[si], Ys_te[si], 'R2', exclude_tfs=exclude_tfs)
            err_dicts[si]['R2'][fold, 0] = R2
            
            Xsi = Xs[si]
            Ysi = Ys[si]
            Bsi = Bs[si]

            if True: #always rescale the beta matrix now
                S = rescale_betas(Xsi, Ysi, Bsi)
            else:
                S = Bsi
            #aupr and prc curves
            (aupr, prc) = eval_network_pr(S, genes[si], tfs[si], priors_te[si], tr_priors=priors_tr[si], exclude_tfs=exclude_tfs, constraints = None, test_all=test_all)
            err_dicts[si]['aupr'][fold,0] = aupr
            err_dicts[si]['prc'][fold] = prc            
            
            (aupr_noncon, prc_noncon) = eval_network_pr(S, genes[si], tfs[si], priors_te[si], tr_priors=priors_tr[si], exclude_tfs=exclude_tfs, constraints = constraints, non_con=True, sub = si, test_all=test_all)
            err_dicts[si]['aupr_noncon'][fold,0] = aupr_noncon                    
            err_dicts[si]['prc_noncon'][fold] = downsample_roc(prc_noncon)

            #constrained aupr and prc curves
            (aupr_con, prc_con) = eval_network_pr(S, genes[si], tfs[si], priors_te[si], tr_priors=priors_tr[si], exclude_tfs=exclude_tfs, constraints = constraints, sub = si, test_all=test_all)
            err_dicts[si]['aupr_con'][fold,0] = aupr_con                
            err_dicts[si]['prc_con'][fold] = downsample_roc(prc_con)

            #auc and roc curves
            (auc, roc) = eval_network_roc(S, genes[si], tfs[si], priors_te[si], tr_priors=priors_tr[si], exclude_tfs=exclude_tfs, constraints = None, test_all=test_all)
            err_dicts[si]['auc'][fold,0] = auc
            err_dicts[si]['roc'][fold] = downsample_roc(roc)

            #constrained auc and roc curves
            (auc_con, roc_con) = eval_network_roc(S, genes[si], tfs[si], priors_te[si], tr_priors=priors_tr[si], exclude_tfs=exclude_tfs, constraints = constraints, sub = si, test_all=test_all)
            err_dicts[si]['auc_con'][fold,0] = auc_con
            err_dicts[si]['roc_con'][fold] = downsample_roc(roc_con)

            (auc_noncon, roc_noncon) = eval_network_roc(S, genes[si], tfs[si], priors_te[si], tr_priors=priors_tr[si], exclude_tfs=exclude_tfs, constraints = constraints, non_con=True, sub = si,test_all=test_all)
            err_dicts[si]['auc_noncon'][fold,0] = auc_noncon
            err_dicts[si]['roc_noncon'][fold] = downsample_roc(roc_noncon)
            
            #beta error if data is simulated
            betafile = os.path.join(data_fn, 'beta%d' % (si+1))
            if os.path.exists(betafile):
                B_mse = mse_B(Bs[si], betafile)
                err_dicts[si]['B_mse'][fold,0] = B_mse
            
            #chance precision and chance constrained interactions precision
            chance = compute_chance_precision(S, genes[si], tfs[si], priors_te[si], tr_priors=priors_tr[si], exclude_tfs=exclude_tfs, constraints = None, sub = si)
            err_dicts[si]['chance'][fold, 0] = chance
            chance_con = compute_chance_precision(S, genes[si], tfs[si], priors_te[si], tr_priors=priors_tr[si], exclude_tfs=exclude_tfs, constraints = constraints, sub = si)
            err_dicts[si]['chance_con'][fold, 0] = chance
            
            #top 100 interactions
            err_dicts[si]['top_100'] = top_k_interactions(S, genes[si], tfs[si], priors_te[si], k=np.inf)
    return err_dicts


#runs the basic model with specified parameters under k-fold cross-validation, and stores a number of metrics
#k: the number of cv folds
#reverse: train on the little dude (reverse train and test)
#cv_both: if false, always use all the data for the corresponding species
#exclude_tfs: don't evaluate on transcription factors. this is useful for generated data, where you can't hope to get them right
def cv_model1(data_fn, lamP, lamR, lamS, k, solver='solve_ortho_direct',special_args=None, reverse=False, cv_both=(True,True), exclude_tfs=True, eval_con=False):
    print 'DEPRACATED cv_model1'


#runs the basic model with specified parameters under k-fold cross-validation, and stores a number of metrics
#returns array for plotting in seaborn
#k: the number of cv folds
#reverse: train on the little dude (reverse train and test)
#cv_both: if false, always use all the data for the corresponding species
#exclude_tfs: don't evaluate on transcription factors. this is useful for generated data, where you can't hope to get them right
def cv_model2(data_fn, lamP, lamR, lamS, k, solver='solve_ortho_direct',special_args=None, reverse=False, cv_both=(True,True), exclude_tfs=True):
    print 'DEPRACATED cv_model2'
    


#runs the basic model with specified parameters under k-fold cross-validation
#stores a bunch of metrics, applied to each CV-fold
#k: the number of cv folds
#reverse: train on the little dude (reverse train and test)
#cv_both: if false, always use all the data for the corresponding species
#exclude_tfs: don't evaluate on transcription factors. this is useful for generated data, where you can't hope to get them right
#doesn't output any files
def cv_model3(data_fn, lamP, lamR, lamS, k, solver='solve_ortho_direct',special_args=None, reverse=False, cv_both=(True,True), exclude_tfs=True):
    print 'DEPRACATED cv_model3'


#cv_model3, but with pct_priors
def cv_model4(data_fn, lamP, lamR, lamS, k, solver='solve_ortho_direct',special_args=None, reverse=False, cv_both=(True,True), exclude_tfs=True, pct_priors=0):
    print 'DEPRACATED cv_model4'

#cv_model1, but with percent_priors
#runs the basic model with specified parameters under k-fold cross-validation, and stores a number of metrics
#percent_priors is the percent of priors to use. these priors are removed from the test set
#k: the number of cv folds
#reverse: train on the little dude (reverse train and test)
#cv_both: if false, always use all the data for the corresponding species
#exclude_tfs: don't evaluate on transcription factors. this is useful for generated data, where you can't hope to get them right
def cv_model5(data_fn, lamP, lamR, lamS, k, solver='solve_ortho_direct',special_args=None, reverse=False, cv_both=(True,True), exclude_tfs=True, eval_con=False, pct_priors=0):
    print 'DEPRACATED cv_model5'

#SECTION: -------------------------CODE FOR EVALUATING THE OUTPUT

def mse_B(Bpred, net_fn, exclude_tfs = True):
    (B, genes, tfs) = ds.load_network(net_fn)
    num_tfs = len(tfs)
    if exclude_tfs:
        start_ind = num_tfs
    else:
        start_ind = 0
    msea = 0.0
    for c in range(start_ind, B.shape[1]):
        bp = Bpred[:,c]
        b = B[:,c]
        mse = ((bp - b)**2).mean()
        msea += mse
    return msea/(B.shape[1] - start_ind)

#model prediction error, using one of several metrics
#exclude_tfs doesn't evaluate predictions of the tfs, and assumes that TFs come before all other genes.
def prediction_error(X, B, Y, metric, exclude_tfs = True):
    Ypred = np.dot(X, B)
    y = Y[:,0]
    yp = Ypred[:,0]
    num_tfs = B.shape[0]
    if exclude_tfs:
        start_ind = num_tfs
    else:
        start_ind = 0
    
    if metric == 'R2':
        r2a = 0.0
        #from matplotlib import pyplot as plt
        for c in range(start_ind, Ypred.shape[1]):
            y = Y[:, c]
            yp = Ypred[:, c]
            r2 = 1 - ((y-yp)**2).sum()/ ((y-y.mean())**2).sum()
            r2a += r2
            
#            if c == start_ind:
#                plt.plot(y)
#                plt.plot(yp)
#                plt.show()
        return r2a/(Ypred.shape[1]-start_ind)
    if metric == 'mse':
        msea = 0.0
        for c in range(start_ind, Ypred.shape[1]):
            y = Y[:, c]
            yp = Ypred[:, c]
            mse = ((y-yp)**2).mean()
            msea += mse

        return msea / (Ypred.shape[1]-start_ind)
    if metric == 'corr':
        corra = 0.0
        for c in range(start_ind, Ypred.shape[1]):
            y = Y[:, c]
            yp = Ypred[:, c]
            corr = np.corrcoef(y, yp)[0,1]
            corra += corr
        return corra / (Ypred.shape[1] - start_ind)

#converts beta matrix B into scale matrix S by an approximation to the rule used in BBSR. The score i, j is assigned to be 1 - residual[j]/(B[i,j]^2*var(X[i]) + residual[j])
def rescale_betas(X, Y, B):
    residual = Y - np.dot(X,B)
    mse = (residual ** 2).mean(axis=0)
    
    vars_tf = np.var(X, axis=0)
    S = np.zeros(B.shape)

    for c in range(S.shape[1]):

        Sc = 1 - mse[c] / (B[:, c]**2 * vars_tf + mse[c])
        S[:, c] = Sc
    return S

#returns scores/labels for aupr or auc
def get_scores_labels(net, genes, tfs, priors, tr_priors=[], exclude_tfs = False, constraints = None, non_con = False, sub=None, test_all='part'):
    #from matplotlib import pyplot as plt
    if len(priors)==0:
        
        return ([], [], [])
    org = priors[0][0].organism
    priors_set = set(priors)
    if test_all == 'all':
        tr_priors_set = set()
    elif test_all == 'part':
        tr_priors_set = set(tr_priors)
    elif test_all == 'gold':
        tr_priors_set = set()
    else:
        print 'test_all must be all, part, or gold'
    gene_to_ind = {genes[x] : x for x in range(len(genes))}
    
    tf_to_ind = {tfs[x] : x for x in range(len(tfs))}
    gene_marked = np.zeros(len(genes)) != 0
    tf_marked = np.zeros(len(tfs)) != 0
    if constraints != None:
        con_set = set()
        for con in constraints:    
            if con.c1.sub == sub:
                con_set.add(con.c1)
            if con.c2.sub == sub:
                con_set.add(con.c2)
        
    #we only evaluate on interactions when the gene/tf is mentioned in a prior
    for prior in priors:
        
        gene_marked[gene_to_ind[prior[0].name]] = True
        gene_marked[gene_to_ind[prior[1].name]] = True
        if prior[0].name in tf_to_ind:
            tf_marked[tf_to_ind[prior[0].name]] = True
        if prior[1].name in tf_to_ind:
            tf_marked[tf_to_ind[prior[1].name]] = True

    genes = np.array(genes)[gene_marked]
    tfs = np.array(tfs)[tf_marked]
    net = net[:, gene_marked]
    net = net[tf_marked, :]
    scores = []#np.zeros(len(genes)*len(tfs))
    labels = []#np.zeros(len(genes)*len(tfs))
    coords = []
    i=0
    for tfi in range(len(tfs)):
        for gi in range(len(genes)):
            if exclude_tfs and gi < len(tfs):
                continue
            coeff = fl.coefficient(sub=sub, r=tfi, c=gi) #potential coefficient
            if constraints != None: 
                if non_con==True:
                    if coeff in con_set:
                        continue      
                else:
                    if not coeff in con_set:
                        continue
            tf = tfs[tfi]
            g = genes[gi]
            coord = (tf, g)
            score = np.abs(net[tfi, gi])
            label = 0
            if (fl.one_gene(tf, org), fl.one_gene(g, org)) in tr_priors_set:
                continue            
            if (fl.one_gene(g, org), fl.one_gene(tf, org)) in tr_priors_set:
                continue
            if (fl.one_gene(tf, org), fl.one_gene(g, org)) in priors_set:
                label = 1
            if (fl.one_gene(g, org), fl.one_gene(tf, org)) in priors_set:
                label = 1
            
            
            scores.append(score)#scores[i] = score
            labels.append(label)#labels[i] = label
            coords.append(coord)
            #print 'wat'
            #print score
            #print label
    return (scores, labels, coords)

#evaluates the area under the precision recall curve, with respect to some given priors
# exclude_tfs: do not evaluate on tf x tf interactions
# constraints: if not None, evaluates only on interactions which have fusion constraints
#sub: name of subproblem. used if constraints != None
#tr_priors are the training set priors
def eval_network_pr(net, genes, tfs, priors, tr_priors=[], exclude_tfs = False, constraints = None, non_con = False, sub=None, test_all='part'):
    (scores, labels, coords) = get_scores_labels(net, genes, tfs, priors, tr_priors, exclude_tfs, constraints, non_con, sub, test_all)
    
    if len(scores) and np.sum(labels)>0:
        #print scores[0:10]
        #print labels[0:10]
        
        (precision, recall,t) = precision_recall_curve(labels, scores)#prc(scores, labels)
        #precision recall, unlike roc_curve, has one fewer t than precision/recall value. I'm just going to copy the last element.
        t = np.concatenate((t, t[[-1]]))
        aupr = auc(recall, precision)
        return (aupr, (recall, precision, t))
    else:
        
        aupr = np.nan
    
        return (aupr, (None, None, None))
        

#evaluates the area under the roc, with respect to some given priors

def eval_network_roc(net, genes, tfs, priors, tr_priors=[], exclude_tfs = True, constraints = None, non_con = False, sub=None, test_all='part'):
    (scores, labels, coords) = get_scores_labels(net, genes, tfs, priors, tr_priors, exclude_tfs, constraints, non_con, sub, test_all)
    if len(scores) and np.sum(scores)>0:
        (fpr, tpr, t) = roc_curve(labels, scores)
        if any(np.isnan(fpr)) or any(np.isnan(tpr)):
            return (0.0, (None, None, None)) #no false positives        
        auroc = auc(fpr, tpr)
        return (auroc, (fpr, tpr, t))
    else:
        auroc = np.nan
        return (auroc, (None, None, None))

    

def eval_network_beta(net1, net2):
    return ((net1 - net2)**2).mean()
            
#generates fusion constraints, then computes the correlation between fused coefficients
def fused_coeff_corr(organisms, genes_l, tfs_l, orth, B_l):
    constraints = fl.orth_to_constraints(organisms, genes_l, tfs_l, orth, 1.0)
    fused_vals = [[],[]]
    
    if len(constraints) == 0:
        return (np.nan, np.zeros((2,0)))
    for con in constraints:
        s1 = con.c1.sub
        b1 = B_l[s1][con.c1.r, con.c1.c]
        s2 = con.c2.sub
        b2 = B_l[s2][con.c2.r, con.c2.c]
        #changed these to 0/1, instead of s1/s2. We don't care whether fusion constraints are within (or between) species
        fused_vals[0].append(b1)
        fused_vals[1].append(b2)
    fused_vals = np.array(fused_vals)
    return (np.corrcoef(fused_vals)[0,1], fused_vals)

#take list of lamP, lamR, lamS values and finds the optimal parameters using cv_model1
def grid_search_params(data_fn, lamP, lamR, lamS, k, solver='solve_ortho_direct',settings=None, reverse=False, cv_both=(True,True), exclude_tfs=True, eval_metric='mse'):

    seed = random.random()
    grid = dict()
    best_mse = 1000
    best_R2 = 0
    best_aupr = 0
    best_auroc = 0
    best_lamP = 1.0
    best_lamR = 0
    best_lamS = 0
    for r in range(len(lamR)):
        for s in range(len(lamS)):
            for p in range(len(lamP)):
                (errd1, errd2) = cv_model_m(data_fn, lamP[p], lamR[r], lamS[s], k, solver='solve_ortho_direct',settings=settings, reverse=False, cv_both=(True,True), exclude_tfs=True, seed=seed, verbose=True)

                if eval_metric == 'mse':
                    grid[str(lamR[r])+'_'+str(lamS[s])+'_'+str(lamP[p])] = errd1['mse']
                    score = errd1['mse'].mean()
                    if score < best_mse:
                        best_mse = score
                        best_lamP = lamP[p]
                        best_lamR = lamR[r]
                        best_lamS = lamS[s]
                if eval_metric == 'R2':
                    grid[str(lamR[r])+'_'+str(lamS[s])+'_'+str(lamP[p])] = errd1['R2']   
                    score = errd1['R2'].mean()
                    if score > best_R2:
                        best_R2 = score
                        best_lamP = lamP[p]
                        best_lamR = lamR[r]
                        best_lamS = lamS[s]
                if eval_metric == 'aupr':
                    grid[str(lamR[r])+'_'+str(lamS[s])+'_'+str(lamP[p])] = errd1['aupr']
                    score = errd1['aupr'].mean()
                    if score > best_aupr:
                        best_aupr = score
                        best_lamP = lamP[p]
                        best_lamR = lamR[r]
                        best_lamS = lamS[s]
                if eval_metric == 'auroc':
                    grid[str(lamR[r])+'_'+str(lamS[s])+'_'+str(lamP[p])] = errd1['auroc']
                    score = errd1['auroc'].mean()
                    if score > best_auroc:
                        best_auroc = score
                        best_lamP = lamP[p]
                        best_lamR = lamR[r]
                        best_lamS = lamS[s]
    if eval_metric == 'mse':
        return (best_mse, best_lamP, best_lamR, best_lamS, grid)
    if eval_metric == 'R2':
        return (best_R2, best_lamP, best_lamR, best_lamS, grid)
    if eval_metric == 'aupr':
        return (best_aupr, best_lamP, best_lamR, best_lamS, grid)
    if eval_metric == 'auroc':
        return (best_auroc, best_lamP, best_lamR, best_lamS, grid)


#given list of roc/prc curves of the form [(precision, recall, t), ...] produces a new set of curves interpolated at every unique value of recall/fpr present on the list all_roc_curves

def pool_roc(roc_curves, all_roc_curves = None, max_x = np.inf):
    all_f = set()
    if all_roc_curves == None:
        all_roc_curves = [roc_curves]

    for roc_curves_group in all_roc_curves:
        for (fs, hs, ts) in roc_curves_group:
            for f in fs:
                all_f.add(f)
    if len(all_f) > max_x:
        fs_interp = np.linspace(min(all_f), max(all_f), max_x)
    else:
        fs_interp = np.array(sorted(list(all_f)))        
    ts_interp = np.zeros((len(roc_curves), len(fs_interp)))
    hs_interp = np.zeros((len(roc_curves), len(fs_interp)))
    def interp_fl(x, xp, fp):
        reverse_interp = np.interp(x[::-1],xp[::-1],fp[::-1])
        return np.fliplr([reverse_interp])[0]
    for i, roc in enumerate(roc_curves):
        
        (fs, hs, ts) = roc
        #from matplotlib import pyplot as plt
        
        #for interp to work, the function must be increasing
        # if it isn't, reverse it, then reverse again
        if fs[0] < fs[-1]:            
            ts_interp[i, :] = np.interp(fs_interp, fs, ts)
            hs_interp[i, :] = np.interp(fs_interp, fs, hs)
        else:
            ts_interp[i, :] = interp_fl(fs_interp, fs, ts)
            hs_interp[i, :] = interp_fl(fs_interp, fs, hs)

        
    return (np.array(fs_interp),np.array(hs_interp), np.array(ts_interp))

#given list of roc/prc curves of the form [(false-alarm, hit, t), ...] produces a new set of Fbeta curves interpolated at every unique value of t

def pool_roc_f(roc_curves, all_roc_curves = None, beta=1.0):
    all_t = set()
    if all_roc_curves == None:
        all_roc_curves = [roc_curves]

    for roc_curves_group in all_roc_curves:
        for (fs, hs, ts) in roc_curves_group:
            for t in ts:
                all_t.add(t)
    ts_interp = sorted(list(all_t))
    fs_interp = np.zeros((len(roc_curves), len(fs_interp)))
    hs_interp = np.zeros((len(roc_curves), len(fs_interp)))
    
    def interp_fl(x, xp, fp):
        reverse_interp = np.interp(x[::-1],xp[::-1],fp[::-1])
        return np.fliplr([reverse_interp])[0]
                                   
    for i, roc in enumerate(roc_curves):
        
        (fs, hs, ts) = roc
        
        #for interp to work, the function must be increasing
        # if it isn't, reverse it, then reverse again
        if fs[0] < fs[-1]:
            fs_interp[i, :] = np.interp(fs_interp, fs, fs)
            hs_interp[i, :] = np.interp(fs_interp, fs, hs)
        else:
            fs_interp[i, :] = interp_fl(fs_interp, fs, fs)
            hs_interp[i, :] = interp_fl(fs_interp, fs, hs)
            
    return (fs_interp,hs_interp, ts_interp)

#returns the fraction of true labels for the given constraints and priors
def compute_chance_precision(net, genes, tfs, priors, tr_priors=[], exclude_tfs = False, constraints = None, sub=None):
    
    (scores, labels, coords) = get_scores_labels(net, genes, tfs, priors, tr_priors, exclude_tfs, constraints, sub)
    if len(scores) == 0:
        return np.nan
    return np.mean(labels)
    
#returns a list of the top k interactions (as tf gene name tuples) along with labels specifying whether they are in the prior set
#restricts to top k interactions for which some prior exists
def top_k_interactions(net, genes, tfs, priors, org=0, k=np.inf, restrict_genes_w_priors = True):
    (r, c) = np.indices(net.shape)
    rc_inds = zip(r.ravel(), c.ravel())
    
    rc_inds.sort(key = lambda rc: net[rc[0], rc[1]])
    
    k = min(k, len(rc_inds))
    priors_set = set(priors)
    rc_inds_k = rc_inds[0:k]
    interaction_names = []
    labels = []
    
    genes_w_priors = set()
    for (tf_og, gene_og) in priors:
        genes_w_priors.add(tf_og)
        genes_w_priors.add(gene_og)

    assigned = 0
    for ki, (r, c) in enumerate(rc_inds):
        if assigned >= k:
            break
        tf = tfs[r]
        gene = genes[c]
        label = 0
        tf_og = fl.one_gene(tf, org)
        gene_og = fl.one_gene(gene, org)
        if not restrict_genes_w_priors or tf_og in genes_w_priors or gene_og in genes_w_priors:
            if (tf_og, gene_og) in priors_set:
                label = 1
            if (gene_og, tf_og) in priors_set:
                label = 1
            labels.append(label)
            interaction_names.append((tf, gene))
            assigned += 1
    
    return (interaction_names, labels)
