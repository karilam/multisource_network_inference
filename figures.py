import pandas as pd
import seaborn as sns
import collections
import numpy as np
import os
import fit_grn as fg
import data_sources as ds
import fused_reg as fr
import matplotlib.pyplot as plt
import experiments as e
import random

def xspecies_perf(lamP=(0.03,0.05), lamR=(0.368,0.5), lamSs=[0,0.5,1], k=10,cv_both=(True,False), orgs=['B_subtilis','B_anthracis'], orth_file=['orth']):

    out = os.path.join('data','bacteria_standard')
    seed = np.random.randn()
    settings = fr.get_settings()
    solver = 'solve_ortho_direct'
    errdf = fg.cv_model_m(out, lamP=lamP, lamR=lamR, lamS=lamS, k=k, solver=solver, settings=settings, reverse=True, cv_both=cv_both, exclude_tfs=False, seed=seed, orgs=orgs, lamS_opt=None, orth_file=orth_file)
    #errdu = fg.cv_model_m(out, lamP=lamP, lamR=lamR, lamS=0, k=k, solver=solver, settings=settings, reverse=True, cv_both=cv_both, exclude_tfs=False, seed=seed, orgs=orgs, lamS_opt=None, orth_file=orth_file)
    #sanity check - errdu gives same answer as errdf

    a=collections.defaultdict(lambda:[])
    for i in range(len(errdf)):
        for j in range(len(errdf[i]['aupr'])):
            a['lamS'].append(lamSs[i])
            a['AUPR'].append(errdf[i]['aupr'][j][0])

    df = pd.DataFrame(a)
    sns.barplot(x="lamS", y="AUPR", data=df)
    plt.show()

    #return (errdf, errdu)

def xspecies_perf2(lamP=(0.03,0.05), lamR=(0.368,0.5), lamSs=[0,0.5,1], k=10,cv_both=(True,False), orgs=['B_subtilis','B_anthracis'], orth_file=['orth']):

    out = os.path.join('data','bacteria_standard')
    seed = np.random.randn()
    settings = fr.get_settings()
    solver = 'solve_ortho_direct'
    errd = []
    for i in range(len(lamSs)):
        errdf = fg.cv_model_m(out, lamP=lamP, lamR=lamR, lamS=lamSs[i], k=k, solver=solver, settings=settings, reverse=True, cv_both=cv_both, exclude_tfs=False, seed=seed, orgs=orgs, lamS_opt=None, orth_file=orth_file)
        errd.append(errdf)
    #errdu = fg.cv_model_m(out, lamP=lamP, lamR=lamR, lamS=0, k=k, solver=solver, settings=settings, reverse=True, cv_both=cv_both, exclude_tfs=False, seed=seed, orgs=orgs, lamS_opt=None, orth_file=orth_file)
    #sanity check - errdu gives same answer as errdf
    return errd

    a=collections.defaultdict(lambda:[])
    for i in range(len(errdf)):
        for j in range(len(errdf[i]['aupr'])):
            a['lamS'].append(lamSs[i])
            a['AUPR'].append(errdf[i]['aupr'][j][0])

def xspecies_perf3(lamP=(1.0,1.0), lamR=(0.368,0.5), lamSs=[0,0.5,5], k=20,cv_both=(True,False), orgs=['B_subtilis','B_anthracis'], orth_file=['orth']):

    out = os.path.join('data','bacteria_standard')
    seed = np.random.randn()
    settings = fr.get_settings()
    solver = 'solve_ortho_direct'
    errd = []
    for i in range(len(lamSs)):
        errdf = fg.cv_model_m(out, lamP=lamP, lamR=lamR, lamS=lamSs[i], k=k, solver=solver, settings=settings, reverse=True, cv_both=cv_both, exclude_tfs=False, seed=seed, orgs=orgs, lamS_opt=None, orth_file=orth_file)
        errd.append(errdf)
    #errdu = fg.cv_model_m(out, lamP=lamP, lamR=lamR, lamS=0, k=k, solver=solver, settings=settings, reverse=True, cv_both=cv_both, exclude_tfs=False, seed=seed, orgs=orgs, lamS_opt=None, orth_file=orth_file)
    #sanity check - errdu gives same answer as errdf
    return errd

    a=collections.defaultdict(lambda:[])
    for i in range(len(errdf)):
        for j in range(len(errdf[i]['aupr'])):
            a['lamS'].append(lamSs[i])
            a['AUPR'].append(errdf[i]['aupr'][j][0])

def xspecies_perf4(lamP=(0.03,0.05), lamR=(0.368,0.5), lamSs=[0,0.5,1], k=5,cv_both=(True,False), orgs=['B_subtilis','B_anthracis'], orth_file=['orth']):

    out = os.path.join('data','bacteria_standard')
    seed = np.random.randn()
    settings = fr.get_settings()
    solver = 'solve_ortho_direct'
    errd = []
    for i in range(len(lamSs)):
        errdf = fg.cv_model_m(out, lamP=lamP, lamR=lamR, lamS=lamSs[i], k=k, solver=solver, settings=settings, reverse=True, cv_both=cv_both, exclude_tfs=False, seed=seed, orgs=orgs, lamS_opt=None, orth_file=orth_file, use_TFA=0.5)
        errd.append(errdf)
    #errdu = fg.cv_model_m(out, lamP=lamP, lamR=lamR, lamS=0, k=k, solver=solver, settings=settings, reverse=True, cv_both=cv_both, exclude_t   fs=False, seed=seed, orgs=orgs, lamS_opt=None, orth_file=orth_file)
    #sanity check - errdu gives same answer as errdf
    return errd

def xspecies_perf5(lamP=(1.0,1.0), lamR=(5,5), lamSs=[0,0.5,5], k=10,cv_both=(True,False), orgs=['B_subtilis','B_anthracis'], orth_file=['orth']):

    out = os.path.join('data','bacteria_standard')
    seed = np.random.randn()
    settings = fr.get_settings()
    solver = 'solve_ortho_direct'
    errd = []
    for i in range(len(lamSs)):
        errdf = fg.cv_model_m(out, lamP=lamP, lamR=lamR, lamS=lamSs[i], k=k, solver=solver, settings=settings, reverse=True, cv_both=cv_both, exclude_tfs=False, seed=seed, orgs=orgs, lamS_opt=None, orth_file=orth_file)
        errd.append(errdf)
    #errdu = fg.cv_model_m(out, lamP=lamP, lamR=lamR, lamS=0, k=k, solver=solver, settings=settings, reverse=True, cv_both=cv_both, exclude_tfs=False, seed=seed, orgs=orgs, lamS_opt=None, orth_file=orth_file)
    #sanity check - errdu gives same answer as errdf
    return errd


def datasources_perf(lamP=(0.03,0.007), lamR=(0.368,0.0789), lamSs=[0,0.5], k=10,cv_both=(True,True), orgs=['B_subtilis','B_subtilis_eu'], orth_file=['orth']):
    out = os.path.join('data','bacteria_standard')

    seed = np.random.randn()

    settings = fr.get_settings()
    solver = 'solve_ortho_direct'

    errdr = fg.cv_unfused(out, lamP=lamP, lamR=lamR, k=k, solver = solver, settings = settings, reverse = True, cv_both = cv_both, exclude_tfs=False, seed = seed, orgs = orgs, lamS_opt = None)

    errdf=[]
    for i, lamS in enumerate(lamSs):
        errd = fg.cv_model_m(out, lamP=lamP, lamR=lamR, lamS=lamS, k=k, solver=solver, settings=settings, reverse=True, cv_both=cv_both, exclude_tfs=False, seed=seed, orgs=orgs, lamS_opt=None, orth_file=orth_file)
        errdf.append(errd)

    a = collections.defaultdict(lambda:[])

    for i in range(len(errdr['aupr'])):
        a['method'].append('rank combine')
        a['AUPR'].append(errdr['aupr'][i][0])

    for i in range(len(errdf[0])):
        for j in range(len(errdf[0][i]['aupr'])):
            a['method'].append('unfused')
            a['AUPR'].append(errdf[0][i]['aupr'][j][0])

    for i in range(len(errdf[1])):
        for j in range(len(errdf[1][i]['aupr'])):
            a['method'].append('fused-0.5')
            a['AUPR'].append(errdf[1][i]['aupr'][j][0])

    df = pd.DataFrame(a)
    sns.barplot(x='method',y='AUPR',data=df)
    plt.show()

    metric='prc'
    savef='datasources_perf'
    normed=False
    scad=False
    roc_species=0
    unfused=True
    lamS_opt=None
    
    e.plot_bacteria_roc(lamP, lamR, lamSs, k, metric, savef,normed, scad, cv_both, roc_species, orgs, unfused, lamS_opt, orth_file)

    return (errdr, errdf)

def datasources_perf3(lamP=(1,1), lamR=(0.5,0.5), lamSs=[0,1], k=10,cv_both=(True,True), orgs=['B_subtilis','B_subtilis_eu'], orth_file=['orth']):
    out = os.path.join('data','bacteria_standard')

    seed = np.random.randn()

    settings = fr.get_settings()
    solver = 'solve_ortho_direct'

    errdr = fg.cv_unfused(out, lamP=lamP, lamR=lamR, k=k, solver = solver, settings = settings, reverse = False, cv_both = cv_both, exclude_tfs=False, seed = seed, orgs = orgs, lamS_opt = None)

    errdf=[]
    for i, lamS in enumerate(lamSs):
        errd = fg.cv_model_m(out, lamP=lamP, lamR=lamR, lamS=lamS, k=k, solver=solver, settings=settings, reverse=False, cv_both=cv_both, exclude_tfs=False, seed=seed, orgs=orgs, lamS_opt=None, orth_file=orth_file)
        errdf.append(errd)

    return (errdr, errdf)

def datasources_perf_NOREVERSE(lamP=(0.03,0.007), lamR=(0.368,0.0789), lamSs=[0,0.5], k=10,cv_both=(True,True), orgs=['B_subtilis','B_subtilis_eu'], orth_file=['orth']):
    out = os.path.join('data','bacteria_standard')

    seed = np.random.randn()

    settings = fr.get_settings()
    solver = 'solve_ortho_direct'

    errdr = fg.cv_unfused(out, lamP=lamP, lamR=lamR, k=k, solver = solver, settings = settings, reverse = False, cv_both = cv_both, exclude_tfs=False, seed = seed, orgs = orgs, lamS_opt = None)

    errdf=[]
    for i, lamS in enumerate(lamSs):
        errd = fg.cv_model_m(out, lamP=lamP, lamR=lamR, lamS=lamS, k=k, solver=solver, settings=settings, reverse=False, cv_both=cv_both, exclude_tfs=False, seed=seed, orgs=orgs, lamS_opt=None, orth_file=orth_file)
        errdf.append(errd)

    a = collections.defaultdict(lambda:[])

    for i in range(len(errdr['aupr'])):
        a['method'].append('rank combine')
        a['AUPR'].append(errdr['aupr'][i][0])

    for i in range(len(errdf[0])):
        for j in range(len(errdf[0][i]['aupr'])):
            a['method'].append('unfused')
            a['AUPR'].append(errdf[0][i]['aupr'][j][0])

    for i in range(len(errdf[1])):
        for j in range(len(errdf[1][i]['aupr'])):
            a['method'].append('fused-0.5')
            a['AUPR'].append(errdf[1][i]['aupr'][j][0])

    df = pd.DataFrame(a)
    sns.barplot(x='method',y='AUPR',data=df)
    plt.show()

    metric='prc'
    savef='datasources_perf_noreverse2'
    normed=False
    scad=False
    roc_species=0
    unfused=True
    lamS_opt=None
    
    e.plot_bacteria_roc(lamP, lamR, lamSs, k, metric, savef,normed, scad, cv_both, roc_species, orgs, unfused, lamS_opt, orth_file)

    return (errdr, errdf)



#makes scatter plot of Bs for true/false constraints under fused l2, unfused, scad
#makes histogram of real B - fusedl2/scad/unfused B for true/false constraints
#makes histogram of deltabetas for true/false constraints
#plots histogram of lamS for true/false constraints
def scad_tf():
    if not os.path.exists(os.path.join('data','fake_data','fused_w_err')):
        os.mkdir(os.path.join('data','fake_data','fused_w_err'))

    out = os.path.join('data','fake_data','fused_w_err')
    #out = os.path.join('data','fake_data','deltb_lams')

    N_TF = 15
    N_G = 200
    amt_fused = 1.0
    orth_err = 0.5
    lamS = 0.5
    ds.write_fake_data1(N1 = 10, N2 = 15, out_dir = out, tfg_count1=(N_TF, N_G), tfg_count2 = (N_TF, N_G), pct_fused = amt_fused, orth_falsepos = orth_err, orth_falseneg = orth_err, measure_noise1 = 0.1, measure_noise2 = 0.1, sparse=0.0, fuse_std = 0.1)

    scad_settings = fr.get_settings({'s_it':5, 'per':40, 'return_cons':True})
    #'per':((1/(1+float(orth_err)))*100), 'return_cons' : True})

    lamP = 1
    lamR = 2
    (errd1s, errd2s) = fg.cv_model_m(out, lamP=lamP, lamR=lamR, lamS=lamS, k=5, solver='solve_ortho_direct_scad', reverse = False, settings = scad_settings, cv_both = (True, True))

    #(errd1u, errd2u) = fg.cv_model_m(out, lamP=lamP, lamR=lamR, lamS=0, k=5, solver='solve_ortho_direct', reverse = False, settings = None, cv_both = (True, True))

    #(errd1f, errd2f) = fg.cv_model_m(out, lamP=lamP, lamR=lamR,lamS=lamS, k=5, solver='solve_ortho_direct', reverse=False, settings=None, cv_both=(True,True))

    ds1 = ds.standard_source(out,0)
    ds2 = ds.standard_source(out,1)
    orth_fn=os.path.join(out, 'orth')
    organisms=[ds1.name, ds2.name]
    orth = ds.load_orth(orth_fn, organisms)
    (e1, t1, genes1, tfs1) = ds1.load_data()
    (e2, t2, genes2, tfs2) = ds2.load_data()
    Xs = [t1, t2]
    Ys = [e1, e2]
    (priors1, signs1) = ds1.get_priors()
    (priors2, signs2) = ds2.get_priors()
    priors = priors1 + priors2
    genes = [genes1, genes2]
    tfs = [tfs1, tfs2]
    (constraints, marks) = fr.orth_to_constraints_marked(organisms, genes, tfs, orth, lamS)
    betafile1 = os.path.join(out, 'beta1')
    betafile2 = os.path.join(out, 'beta2')
    
    lamS=0
    Bs_uf = fr.solve_ortho_direct(organisms, genes, tfs, Xs, Ys, orth, priors, lamP, lamR, lamS)
    #lamS=8
    lamS=4
    Bs_fr = fr.solve_ortho_direct(organisms, genes, tfs, Xs, Ys, orth, priors, lamP, lamR, lamS)
    lamS=4
    Bs_fs = fr.solve_ortho_direct_scad(organisms, genes, tfs, Xs, Ys, orth, priors, lamP, lamR, lamS, scad_settings)


    r1 = Bs_uf[0].shape[0]
    c1 = Bs_uf[0].shape[1]
    r2 = Bs_uf[1].shape[0]
    c2 = Bs_uf[1].shape[1]
    Buf1 = []
    Buf2 = []
    Bfr1 = []
    Bfr2 = []
    Bfs1 = []
    Bfs2 = []
    Bufd = []
    Bfrd = []
    Bfsd = []

    colors = []
    con_inds = np.random.permutation(range(len(constraints)))
    area_fr = []
    area_fs = []
    area_fm = []
    cons = []
    cons_fs = []

    subs = 2000
    for i in con_inds:#[0:subs]:
        con = constraints[i]
        cons.append(con.lam)
        mark = marks[i]
        Buf1.append(Bs_uf[con.c1.sub][con.c1.r, con.c1.c])
        Buf2.append(Bs_uf[con.c2.sub][con.c2.r, con.c2.c])
        Bufd.append(Bs_uf[con.c1.sub][con.c1.r, con.c1.c]-Bs_uf[con.c2.sub][con.c2.r, con.c2.c])
        Bfr1.append(Bs_fr[con.c1.sub][con.c1.r, con.c1.c])
        Bfr2.append(Bs_fr[con.c2.sub][con.c2.r, con.c2.c])
        Bfrd.append(Bs_fr[con.c1.sub][con.c1.r, con.c1.c]-Bs_fr[con.c2.sub][con.c2.r, con.c2.c])

        con_fs = scad_settings['cons'][i]
        cons_fs.append(con_fs.lam)
        Bfs1.append(Bs_fs[con_fs.c1.sub][con_fs.c1.r, con_fs.c1.c])
        Bfs2.append(Bs_fs[con_fs.c2.sub][con_fs.c2.r, con_fs.c2.c])
        Bfsd.append(Bs_fs[con_fs.c1.sub][con_fs.c1.r, con_fs.c1.c]-Bs_fs[con.c2.sub][con.c2.r, con.c2.c])

        if mark == 1:
            colors.append('g')
        else:
            colors.append('r')

    Buf1s = np.array(Buf1)
    Buf2s = np.array(Buf2)
    Bufds = np.array(Bufd)
    Bfr1s = np.array(Bfr1)
    Bfr2s = np.array(Bfr2)
    Bfrds = np.array(Bfrd)
    Bfs1s = np.array(Bfs1)
    Bfs2s = np.array(Bfs2)
    Bfsds = np.array(Bfsd)

    plt.close()
    plt.hist(cons_fs)
    plt.title('scad lamS dist')
    plt.show()

    plt.figure()
    plt.subplot(331)
    plt.scatter(Buf1s, Buf2s, c=colors, alpha=0.5)
    plt.xlabel('beta network 1')
    plt.ylabel('beta network 2')
    plt.title('Unfused')

    plt.subplot(332)
    plt.scatter(Bfr1s, Bfr2s, c=colors, alpha=0.5)
    plt.xlabel('beta network 1')
    plt.ylabel('beta network 2')
    plt.title('Fused L2')

    plt.subplot(333)
    plt.scatter(Bfs1s, Bfs2s, c=colors, alpha=0.5)
    plt.xlabel('beta network 1')
    plt.ylabel('beta network 2')
    plt.title('Adaptive fusion') 

    (B1, genes_1, tfs_1) = ds.load_network(betafile1)
    (B2, genes_2, tfs_2) = ds.load_network(betafile1)

    mse_t_uf = []
    mse_f_uf = []
    mse_t_fr = []
    mse_f_fr = []
    mse_t_sc = []
    mse_f_sc = []

    t_uf = []
    f_uf = []
    t_fr = []
    f_fr = []
    t_sc = []
    f_sc = []

    for i in range(len(constraints)):
        con = constraints[i]
        mark = marks[i]
        realdiff = B1[con.c1.r, con.c1.c]-B2[con.c2.r, con.c2.c]
        ufdiff = Bs_uf[con.c1.sub][con.c1.r, con.c1.c] - Bs_uf[con.c2.sub][con.c2.r, con.c2.c]
        frdiff = Bs_fr[con.c1.sub][con.c1.r, con.c1.c] - Bs_fr[con.c2.sub][con.c2.r, con.c2.c]
        scdiff = Bs_fs[con.c1.sub][con.c1.r, con.c1.c] - Bs_fs[con.c2.sub][con.c2.r, con.c2.c]
        Buf = abs(realdiff-ufdiff)
        Bfr = abs(realdiff-ufdiff)
        Bsc = abs(realdiff-scdiff)

        if mark == True:
            mse_t_uf.append(Buf)
            t_uf.append(ufdiff)
            mse_t_fr.append(Bfr)
            t_fr.append(frdiff)
            mse_t_sc.append(Bsc)
            t_sc.append(scdiff)
        if mark == False:
            mse_f_uf.append(Buf)
            f_uf.append(ufdiff)
            mse_f_fr.append(Bfr)
            f_fr.append(frdiff)
            mse_f_sc.append(Bsc)
            f_sc.append(scdiff)

    plt.subplot(334)
    bins = np.linspace(0,5,100)
    #bins = range(int(min(min(mse_t_uf), min(mse_f_uf))), int(max(max(mse_t_uf), max(mse_f_uf))), 100)
    plt.hist(mse_t_uf, bins, alpha=0.5, label='real')
    plt.hist(mse_f_uf, bins, alpha=0.5, label='fake')
    plt.ylim(0,1400)
    plt.xlim(0,4)
    plt.legend()

    #bins = range(int(min(min(mse_t_fr), min(mse_f_uf))), int(max(max(mse_t_fr), max(mse_f_fr))), 100)
    plt.subplot(335)
    plt.hist(mse_t_fr, bins, alpha=0.5, label='real')
    plt.hist(mse_f_fr, bins, alpha=0.5, label='fake')
    plt.ylim(0,1400)
    plt.xlim(0,4) 
    plt.legend()

    #bins = range(int(min(min(mse_t_sc), min(mse_f_sc))), int(max(max(mse_t_sc), max(mse_f_sc))), 100)
    plt.subplot(336)
    plt.hist(mse_t_sc, bins, alpha=0.5, label='real')
    plt.hist(mse_f_sc, bins, alpha=0.5, label='fake')
    plt.ylim(0,1400)
    plt.xlim(0,4)
    plt.legend()

    bins = np.linspace(min(f_fr), max(f_fr), 100)

    plt.subplot(337)
    plt.hist(t_uf, bins, alpha=0.5, label='real')
    plt.hist(f_uf, bins, alpha=0.5, label='fake')
    plt.legend()

    plt.subplot(338)
    plt.hist(t_fr, bins, alpha=0.5, label='real')
    plt.hist(f_fr, bins, alpha=0.5, label='fake')
    plt.legend()

    plt.subplot(339)
    plt.hist(t_sc, bins, alpha=0.5, label='real')
    plt.hist(f_sc, bins, alpha=0.5, label='fake')
    plt.legend()
    plt.show()

    e.plot_fuse_lams(out, scad_settings['cons'])
    #e.plot_betas_scad(out)

    return (errd1s, errd2s, errd1u, errd2u, errd1f, errd2f)




def scad_tf2():
    if not os.path.exists(os.path.join('data','fake_data','fused_w_err')):
        os.mkdir(os.path.join('data','fake_data','fused_w_err'))

    out = os.path.join('data','fake_data','fused_w_err')
    #out = os.path.join('data','fake_data','deltb_lams')

    N_TF = 35
    N_G = 200
    amt_fused = 1.0
    orth_err = 0.4
    lamS =10
    ds.write_fake_data1(N1 = 20, N2 = 20, out_dir = out, tfg_count1=(N_TF, N_G), tfg_count2 = (N_TF, N_G), pct_fused = amt_fused, orth_falsepos = orth_err, orth_falseneg = orth_err, measure_noise1 = 0.1, measure_noise2 = 0.1, sparse=0.0, fuse_std = 0.1)

    scad_settings = fr.get_settings({'s_it':5, 'per':95, 'return_cons':True})
    #scad_settings = fr.get_settings({'s_it':5, 'a':5000, 'return_cons':True})
    #'per':((1/(1+float(orth_err)))*100), 'return_cons' : True})

    lamP = 1
    lamR = 1
    #(errd1s, errd2s) = fg.cv_model_m(out, lamP=lamP, lamR=lamR, lamS=lamS, k=5, solver='solve_ortho_direct_scad', reverse = False, settings = scad_settings, cv_both = (True, True))

    #errd1u, errd2u) = fg.cv_model_m(out, lamP=lamP, lamR=lamR, lamS=0, k=5, solver='solve_ortho_direct', reverse = False, settings = None, cv_both = (True, True))

    #(errd1f, errd2f) = fg.cv_model_m(out, lamP=lamP, lamR=lamR,lamS=lamS, k=5, solver='solve_ortho_direct', reverse=False, settings=None, cv_both=(True,True))

    ds1 = ds.standard_source(out,0)
    ds2 = ds.standard_source(out,1)
    orth_fn=os.path.join(out, 'orth')
    organisms=[ds1.name, ds2.name]
    orth = ds.load_orth(orth_fn, organisms)
    (e1, t1, genes1, tfs1) = ds1.load_data()
    (e2, t2, genes2, tfs2) = ds2.load_data()
    Xs = [t1, t2]
    Ys = [e1, e2]
    (priors1, signs1) = ds1.get_priors()
    (priors2, signs2) = ds2.get_priors()
    priors = priors1 + priors2
    genes = [genes1, genes2]
    tfs = [tfs1, tfs2]
    (constraints, marks) = fr.orth_to_constraints_marked(organisms, genes, tfs, orth, lamS)
    betafile1 = os.path.join(out, 'beta1')
    betafile2 = os.path.join(out, 'beta2')

    lamS=0
    Bs_uf = fr.solve_ortho_direct(organisms, genes, tfs, Xs, Ys, orth, priors, lamP, lamR, lamS)
    #lamS=8
    lamS=10
    Bs_fr = fr.solve_ortho_direct(organisms, genes, tfs, Xs, Ys, orth, priors, lamP, lamR, lamS)

    Bs_fs = fr.solve_ortho_direct_scad(organisms, genes, tfs, Xs, Ys, orth, priors, lamP, lamR, lamS, scad_settings)

    #Bs_fr2 = fr.direct_solve_factor(Xs, Ys, fuse_con, ridge_con, lamR, adjust = scad_settings['adjust'])

    #r1 = Bs_uf[0].shape[0]
    #c1 = Bs_uf[0].shape[1]
    #r2 = Bs_uf[1].shape[0]
    #c2 = Bs_uf[1].shape[1]
    Buf1 = []
    Buf2 = []
    Bfr1 = []
    Bfr2 = []
    Bfs1 = []
    Bfs2 = []

    Bufd = []
    Bfrd = []
    Bfsd = []

    colors = []
    con_inds = np.random.permutation(range(len(constraints)))
    area_fr = []
    area_fs = []
    area_fm = []
    cons = []
    cons_fs = []

    subs = 2000
    for i in con_inds:#[0:subs]:
        con = constraints[i]
        cons.append(con.lam)
        mark = marks[i]
        Buf1.append(Bs_uf[con.c1.sub][con.c1.r, con.c1.c])
        Buf2.append(Bs_uf[con.c2.sub][con.c2.r, con.c2.c])
        Bufd.append(Bs_uf[con.c1.sub][con.c1.r, con.c1.c]-Bs_uf[con.c2.sub][con.c2.r, con.c2.c])
        Bfr1.append(Bs_fr[con.c1.sub][con.c1.r, con.c1.c])
        Bfr2.append(Bs_fr[con.c2.sub][con.c2.r, con.c2.c])
        Bfrd.append(Bs_fr[con.c1.sub][con.c1.r, con.c1.c]-Bs_fr[con.c2.sub][con.c2.r, con.c2.c])

        con_fs = scad_settings['cons'][i]
        cons_fs.append(con_fs.lam)
        Bfs1.append(Bs_fs[con_fs.c1.sub][con_fs.c1.r, con_fs.c1.c])
        Bfs2.append(Bs_fs[con_fs.c2.sub][con_fs.c2.r, con_fs.c2.c])
        Bfsd.append(Bs_fs[con_fs.c1.sub][con_fs.c1.r, con_fs.c1.c]-Bs_fs[con.c2.sub][con.c2.r, con.c2.c])

        if mark == 1:
            colors.append('g')
        else:
            colors.append('r')

    plt.figure()
    plt.subplot(131)
    plt.scatter(Buf1, Buf2, c=colors, alpha=0.5)
    plt.xlabel('beta network 1')
    plt.ylabel('beta network 2')
    plt.axis('equal') 
    plt.title('Unfused')

    plt.subplot(132)
    plt.scatter(Bfr1, Bfr2, c=colors, alpha=0.5)
    plt.xlabel('beta network 1')
    plt.ylabel('beta network 2')
    plt.axis('equal') 
    plt.title('Fused L2')

    plt.subplot(133)
    plt.scatter(Bfs1, Bfs2, c=colors, alpha=0.5)
    plt.xlabel('beta network 1')
    plt.ylabel('beta network 2')
    plt.title('Adaptive fusion')
    plt.axis('equal') 
    plt.show()

    (B1, genes_1, tfs_1) = ds.load_network(betafile1)
    (B2, genes_2, tfs_2) = ds.load_network(betafile1)

    a=collections.defaultdict(lambda:[])

    for i in range(len(constraints)):
        con = constraints[i]
        mark = marks[i]
        realdiff = B1[con.c1.r, con.c1.c]-B2[con.c2.r, con.c2.c]
        ufdiff1 = abs(Bs_uf[con.c1.sub][con.c1.r,con.c1.c] - B1[con.c1.r,con.c1.c])
        ufdiff2 = abs(Bs_uf[con.c2.sub][con.c2.r,con.c2.c] - B2[con.c2.r,con.c1.c])
        frdiff1 = abs(Bs_fr[con.c1.sub][con.c1.r,con.c1.c] - B1[con.c1.r,con.c1.c])
        frdiff2 = abs(Bs_fr[con.c2.sub][con.c2.r,con.c2.c] - B2[con.c2.r,con.c1.c])
        scdiff1 = abs(Bs_fs[con.c1.sub][con.c1.r,con.c1.c] - B1[con.c1.r,con.c1.c])
        scdiff2 = abs(Bs_fs[con.c2.sub][con.c2.r,con.c2.c] - B2[con.c2.r,con.c1.c])

        if mark == True:
            a['MSE'].append(ufdiff1)
            a['solver'].append('unfused')
            a['ortholog'].append('True')
            a['MSE'].append(ufdiff2)
            a['solver'].append('unfused')
            a['ortholog'].append('True')
            a['MSE'].append(frdiff1)
            a['solver'].append('fused L2')
            a['ortholog'].append('True')
            a['MSE'].append(frdiff2)
            a['solver'].append('fused L2')
            a['ortholog'].append('True')
            a['MSE'].append(scdiff1)
            a['solver'].append('adaptive fusion')
            a['ortholog'].append('True')
            a['MSE'].append(scdiff2)
            a['solver'].append('adaptive fusion')
            a['ortholog'].append('True')
        if mark == False:
            a['MSE'].append(ufdiff1)
            a['solver'].append('unfused')
            a['ortholog'].append('False')
            a['MSE'].append(ufdiff2)
            a['solver'].append('unfused')
            a['ortholog'].append('False')
            a['MSE'].append(frdiff1)
            a['solver'].append('fused L2')
            a['ortholog'].append('False')
            a['MSE'].append(frdiff2)
            a['solver'].append('fused L2')
            a['ortholog'].append('False')
            a['MSE'].append(scdiff1)
            a['solver'].append('adaptive fusion')
            a['ortholog'].append('False')
            a['MSE'].append(scdiff2)
            a['solver'].append('adaptive fusion')
            a['ortholog'].append('False')

    df = pd.DataFrame(a)
    sns.barplot(x='solver',y='MSE',hue='ortholog',data=df)
    plt.show()    
    return (scad_settings['cons'])



def scad_tf3():
    if not os.path.exists(os.path.join('data','fake_data','fused_w_err')):
        os.mkdir(os.path.join('data','fake_data','fused_w_err'))

    out = os.path.join('data','fake_data','fused_w_err')
    #out = os.path.join('data','fake_data','deltb_lams')

    N_TF = 35
    N_G = 200
    amt_fused = 1.0
    orth_err = 0.4
    lamS =15
    ds.write_fake_data1(N1 = 20, N2 = 20, out_dir = out, tfg_count1=(N_TF, N_G), tfg_count2 = (N_TF, N_G), pct_fused = amt_fused, orth_falsepos = orth_err, orth_falseneg = orth_err, measure_noise1 = 0.1, measure_noise2 = 0.1, sparse=0.0, fuse_std = 0.1)

    scad_settings = fr.get_settings({'s_it':5, 'per':90, 'return_cons':True})
    #scad_settings = fr.get_settings({'s_it':5, 'a':5000, 'return_cons':True})
    #'per':((1/(1+float(orth_err)))*100), 'return_cons' : True})

    lamP = 1
    lamR = 2
    #(errd1s, errd2s) = fg.cv_model_m(out, lamP=lamP, lamR=lamR, lamS=lamS, k=5, solver='solve_ortho_direct_scad', reverse = False, settings = scad_settings, cv_both = (True, True))

    #errd1u, errd2u) = fg.cv_model_m(out, lamP=lamP, lamR=lamR, lamS=0, k=5, solver='solve_ortho_direct', reverse = False, settings = None, cv_both = (True, True))

    #(errd1f, errd2f) = fg.cv_model_m(out, lamP=lamP, lamR=lamR,lamS=lamS, k=5, solver='solve_ortho_direct', reverse=False, settings=None, cv_both=(True,True))

    ds1 = ds.standard_source(out,0)
    ds2 = ds.standard_source(out,1)
    orth_fn=os.path.join(out, 'orth')
    organisms=[ds1.name, ds2.name]
    orth = ds.load_orth(orth_fn, organisms)
    (e1, t1, genes1, tfs1) = ds1.load_data()
    (e2, t2, genes2, tfs2) = ds2.load_data()
    Xs = [t1, t2]
    Ys = [e1, e2]
    (priors1, signs1) = ds1.get_priors()
    (priors2, signs2) = ds2.get_priors()
    priors = priors1 + priors2
    genes = [genes1, genes2]
    tfs = [tfs1, tfs2]
    (constraints, marks) = fr.orth_to_constraints_marked(organisms, genes, tfs, orth, lamS)
    betafile1 = os.path.join(out, 'beta1')
    betafile2 = os.path.join(out, 'beta2')

    lamS=0
    Bs_uf = fr.solve_ortho_direct(organisms, genes, tfs, Xs, Ys, orth, priors, lamP, lamR, lamS)
    #lamS=8
    lamS=15
    Bs_fr = fr.solve_ortho_direct(organisms, genes, tfs, Xs, Ys, orth, priors, lamP, lamR, lamS)

    Bs_fs = fr.solve_ortho_direct_scad(organisms, genes, tfs, Xs, Ys, orth, priors, lamP, lamR, lamS, scad_settings)

    r1 = Bs_uf[0].shape[0]
    c1 = Bs_uf[0].shape[1]
    r2 = Bs_uf[1].shape[0]
    c2 = Bs_uf[1].shape[1]
    Buf1 = []
    Buf2 = []
    Bfr1 = []
    Bfr2 = []
    Bfs1 = []
    Bfs2 = []
    Bufd = []
    Bfrd = []
    Bfsd = []

    colors = []
    con_inds = np.random.permutation(range(len(constraints)))
    area_fr = []
    area_fs = []
    area_fm = []
    cons = []
    cons_fs = []

    subs = 2000
    for i in con_inds:#[0:subs]:
        con = constraints[i]
        cons.append(con.lam)
        mark = marks[i]
        Buf1.append(Bs_uf[con.c1.sub][con.c1.r, con.c1.c])
        Buf2.append(Bs_uf[con.c2.sub][con.c2.r, con.c2.c])
        Bufd.append(Bs_uf[con.c1.sub][con.c1.r, con.c1.c]-Bs_uf[con.c2.sub][con.c2.r, con.c2.c])
        Bfr1.append(Bs_fr[con.c1.sub][con.c1.r, con.c1.c])
        Bfr2.append(Bs_fr[con.c2.sub][con.c2.r, con.c2.c])
        Bfrd.append(Bs_fr[con.c1.sub][con.c1.r, con.c1.c]-Bs_fr[con.c2.sub][con.c2.r, con.c2.c])

        con_fs = scad_settings['cons'][i]
        cons_fs.append(con_fs.lam)
        Bfs1.append(Bs_fs[con_fs.c1.sub][con_fs.c1.r, con_fs.c1.c])
        Bfs2.append(Bs_fs[con_fs.c2.sub][con_fs.c2.r, con_fs.c2.c])
        Bfsd.append(Bs_fs[con_fs.c1.sub][con_fs.c1.r, con_fs.c1.c]-Bs_fs[con.c2.sub][con.c2.r, con.c2.c])

        if mark == 1:
            colors.append('g')
        else:
            colors.append('r')

    Buf1s = np.array(Buf1)
    Buf2s = np.array(Buf2)
    Bufds = np.array(Bufd)
    Bfr1s = np.array(Bfr1)
    Bfr2s = np.array(Bfr2)
    Bfrds = np.array(Bfrd)
    Bfs1s = np.array(Bfs1)
    Bfs2s = np.array(Bfs2)
    Bfsds = np.array(Bfsd)

    plt.figure()
    plt.subplot(231)
    plt.scatter(Buf1s, Buf2s, c=colors, alpha=0.5)
    plt.xlabel('beta network 1')
    plt.ylabel('beta network 2')
    plt.title('Unfused')

    plt.subplot(232)
    plt.scatter(Bfr1s, Bfr2s, c=colors, alpha=0.5)
    plt.xlabel('beta network 1')
    plt.ylabel('beta network 2')
    plt.title('Fused L2')

    plt.subplot(233)
    plt.scatter(Bfs1s, Bfs2s, c=colors, alpha=0.5)
    plt.xlabel('beta network 1')
    plt.ylabel('beta network 2')
    plt.title('Adaptive fusion') 

    (B1, genes_1, tfs_1) = ds.load_network(betafile1)
    (B2, genes_2, tfs_2) = ds.load_network(betafile1)

    mse_t_uf = []
    mse_f_uf = []
    mse_t_fr = []
    mse_f_fr = []
    mse_t_sc = []
    mse_f_sc = []

    t_uf = []
    f_uf = []
    t_fr = []
    f_fr = []
    t_sc = []
    f_sc = []

    for i in range(len(constraints)):
        con = constraints[i]
        mark = marks[i]
        realdiff = B1[con.c1.r, con.c1.c]-B2[con.c2.r, con.c2.c]
        ufdiff = Bs_uf[con.c1.sub][con.c1.r, con.c1.c] - Bs_uf[con.c2.sub][con.c2.r, con.c2.c]
        frdiff = Bs_fr[con.c1.sub][con.c1.r, con.c1.c] - Bs_fr[con.c2.sub][con.c2.r, con.c2.c]
        scdiff = Bs_fs[con.c1.sub][con.c1.r, con.c1.c] - Bs_fs[con.c2.sub][con.c2.r, con.c2.c]
        Buf = abs(realdiff-ufdiff)
        Bfr = abs(realdiff-ufdiff)
        Bsc = abs(realdiff-scdiff)

        if mark == True:
            mse_t_uf.append(Buf)
            t_uf.append(ufdiff)
            mse_t_fr.append(Bfr)
            t_fr.append(frdiff)
            mse_t_sc.append(Bsc)
            t_sc.append(scdiff)
        if mark == False:
            mse_f_uf.append(Buf)
            f_uf.append(ufdiff)
            mse_f_fr.append(Bfr)
            f_fr.append(frdiff)
            mse_f_sc.append(Bsc)
            f_sc.append(scdiff)

    plt.subplot(234)
    bins = np.linspace(0,5,100)
    #bins = range(int(min(min(mse_t_uf), min(mse_f_uf))), int(max(max(mse_t_uf), max(mse_f_uf))), 100)
    plt.hist(mse_t_uf, bins, alpha=0.5, label='real')
    plt.hist(mse_f_uf, bins, alpha=0.5, label='fake')
    plt.ylim(0,1400)
    plt.xlim(0,4)
    plt.legend()

    #bins = range(int(min(min(mse_t_fr), min(mse_f_uf))), int(max(max(mse_t_fr), max(mse_f_fr))), 100)
    plt.subplot(235)
    plt.hist(mse_t_fr, bins, alpha=0.5, label='real')
    plt.hist(mse_f_fr, bins, alpha=0.5, label='fake')
    plt.ylim(0,1400)
    plt.xlim(0,4) 
    plt.legend()

    #bins = range(int(min(min(mse_t_sc), min(mse_f_sc))), int(max(max(mse_t_sc), max(mse_f_sc))), 100)
    plt.subplot(236)
    plt.hist(mse_t_sc, bins, alpha=0.5, label='real')
    plt.hist(mse_f_sc, bins, alpha=0.5, label='fake')
    plt.ylim(0,1400)
    plt.xlim(0,4)
    plt.legend()

    plt.show()
#    bins = np.linspace(min(f_fr), max(f_fr), 100)

#    plt.subplot(337)
#    plt.hist(t_uf, bins, alpha=0.5, label='real')
#    plt.hist(f_uf, bins, alpha=0.5, label='fake')
#    plt.legend()

#    plt.subplot(338)
#    plt.hist(t_fr, bins, alpha=0.5, label='real')
#    plt.hist(f_fr, bins, alpha=0.5, label='fake')
#    plt.legend()

#    plt.subplot(339)
#    plt.hist(t_sc, bins, alpha=0.5, label='real')
#    plt.hist(f_sc, bins, alpha=0.5, label='fake')
#    plt.legend()
#    plt.show()
#    (B1, genes_1, tfs_1) = ds.load_network(betafile1)
#    (B2, genes_2, tfs_2) = ds.load_network(betafile1)

    a=collections.defaultdict(lambda:[])

    for i in range(len(constraints)):
        con = constraints[i]
        mark = marks[i]
        realdiff = B1[con.c1.r, con.c1.c]-B2[con.c2.r, con.c2.c]
        ufdiff1 = abs(Bs_uf[con.c1.sub][con.c1.r,con.c1.c] - B1[con.c1.r,con.c1.c])
        ufdiff2 = abs(Bs_uf[con.c2.sub][con.c2.r,con.c2.c] - B2[con.c2.r,con.c1.c])
        frdiff1 = abs(Bs_fr[con.c1.sub][con.c1.r,con.c1.c] - B1[con.c1.r,con.c1.c])
        frdiff2 = abs(Bs_fr[con.c2.sub][con.c2.r,con.c2.c] - B2[con.c2.r,con.c1.c])
        scdiff1 = abs(Bs_fs[con.c1.sub][con.c1.r,con.c1.c] - B1[con.c1.r,con.c1.c])
        scdiff2 = abs(Bs_fs[con.c2.sub][con.c2.r,con.c2.c] - B2[con.c2.r,con.c1.c])

        if mark == True:
            a['MSE'].append(ufdiff1)
            a['solver'].append('unfused')
            a['ortholog'].append('True')
            a['MSE'].append(ufdiff2)
            a['solver'].append('unfused')
            a['ortholog'].append('True')
            a['MSE'].append(frdiff1)
            a['solver'].append('fused L2')
            a['ortholog'].append('True')
            a['MSE'].append(frdiff2)
            a['solver'].append('fused L2')
            a['ortholog'].append('True')
            a['MSE'].append(scdiff1)
            a['solver'].append('adaptive fusion')
            a['ortholog'].append('True')
            a['MSE'].append(scdiff2)
            a['solver'].append('adaptive fusion')
            a['ortholog'].append('True')
        if mark == False:
            a['MSE'].append(ufdiff1)
            a['solver'].append('unfused')
            a['ortholog'].append('False')
            a['MSE'].append(ufdiff2)
            a['solver'].append('unfused')
            a['ortholog'].append('False')
            a['MSE'].append(frdiff1)
            a['solver'].append('fused L2')
            a['ortholog'].append('False')
            a['MSE'].append(frdiff2)
            a['solver'].append('fused L2')
            a['ortholog'].append('False')
            a['MSE'].append(scdiff1)
            a['solver'].append('adaptive fusion')
            a['ortholog'].append('False')
            a['MSE'].append(scdiff2)
            a['solver'].append('adaptive fusion')
            a['ortholog'].append('False')

    df = pd.DataFrame(a)
    sns.barplot(x='solver',y='MSE',hue='ortholog',data=df)
    plt.show()    

#plots penalty function
def show_penalty():
    lamS = 1
    a=1.0
    def miniscad_dx(theta):
        atheta = np.abs(theta)
        if atheta < a/2:
            return atheta * lamS
        if atheta >= a/2:
            return lamS * max(0, (a - atheta))
    def miniL2dx(theta):
        return theta*lamS

    xs = np.linspace(0,1.5,100)
    dms = map(miniscad_dx, xs)
    dl2 = map(miniL2dx, xs)
    l2=[np.trapz(dl2[:i], xs[:i]) for i in range(len(xs))]
    sc=[np.trapz(dms[:i], xs[:i]) for i in range(len(xs))]
    plt.plot(xs, dms,'--')
    plt.plot(xs, dl2)
    plt.legend(('adaptive fusion', 'L2'))
    plt.xlabel('|B0 - B1|')
    plt.ylabel('derivative of penalty')
    plt.ylim(-0.1,1.6)
    plt.plot([a/2,a/2],[-0.1,0.3],color=(0.5,0.5,0.5))
    plt.show()


    #plt.plot(xs, (xs[1]-xs[0])*np.cumsum(dms))
    #plt.show()
    #plt.plot(xs, (xs[1]-xs[0])*np.cumsum(dl2),'--')
    plt.plot(xs, sc, '--')
    plt.plot(xs, l2)
    plt.legend(('adaptive fusion','L2'))
    plt.xlabel('|B0 - B1|')
    plt.ylabel('penalty')
    plt.show()

#makes heatmap showing effect of lamS on networks in L2 fusion; first is ideal case where networks are fused together perfectly; second is case where networks are imperfectly fused
def L2fusion_quick():    
    N_TF = 20
    N_G = 200
    pct_fused = list(np.linspace(0.3,1.0,10))
    reps = 1
    
    lamR = 2
    lamSs = list(np.linspace(0,5,10)) 
    lamP = 1.0
    aupr_array = np.zeros((len(pct_fused),len(lamSs)))

    out1 = os.path.join('data','fake_data','L2fusion_quick')
    k = 2#cv folds
    if not os.path.exists(out1):
        os.mkdir(out1)

    for p in range(reps):
        for i, N in enumerate(pct_fused):
            out = os.path.join(out1,'dat_'+str(N))
            ds.write_fake_data1(N1=10, N2=200, out_dir = out, tfg_count1=(N_TF, N_G), tfg_count2 = (N_TF, N_G), measure_noise1 = 0.1, measure_noise2 = 0.1, sparse=0.75, fuse_std = 0, pct_fused=N, orth_falsepos = 0)        
            lamP = 1.0 #priors don't matter
            seed = 10
            for j, lamS in enumerate(lamSs):
                (errd1, errd2) = fg.cv_model_m(out, lamP, lamR, lamS, k, solver='solve_ortho_direct', reverse=True, cv_both=(True,True), exclude_tfs=True, pct_priors=0, seed=seed, verbose=False)
                #print errd1['corr'].mean()
                print errd1['B_mse'].mean()
                aupr_array[i,j]+= errd1['B_mse'].mean()

    aupr_array /= reps
    ar2 = np.zeros((aupr_array.shape[0],aupr_array.shape[1]))
    for i in range(len(aupr_array)):
        ar2[i,:] = aupr_array[i,:]/aupr_array[i,0]
        #ar2[i,:] = aupr_array[i,:]
    df = pd.DataFrame(ar2[1:,:],index=pct_fused[1:],columns=lamSs)
    df.index.name = 'percent fused'
    df.columns.name = 'lamS'
    plt.figure()
    sns.heatmap(df,cmap="Blues", square=True)
    plt.title('1')
    #lt.show()

    lamSs = list(np.linspace(0,20,10)) 

    out3 = os.path.join('data','fake_data','L2fusion_quick3')
    if not os.path.isdir(out3):
        os.mkdir(out3)
    k = 2#cv folds

    for p in range(reps):
        for i, N in enumerate(pct_fused):
            out = os.path.join(out3,'dat_'+str(N))
            if not os.path.isdir(out):
                os.mkdir(out)
            ds.write_fake_data1(N1=10, N2=200, out_dir = out, tfg_count1=(N_TF, N_G), tfg_count2 = (N_TF, N_G), measure_noise1 = 0.1, measure_noise2 = 0.1, sparse=0.75, fuse_std = 0.75, pct_fused=N, orth_falsepos = 0)        
            lamP = 1.0 #priors don't matter
            seed = 10
            for j, lamS in enumerate(lamSs):
                (errd1, errd2) = fg.cv_model_m(out, lamP, lamR, lamS, k, solver='solve_ortho_direct', reverse=True, cv_both=(True,True), exclude_tfs=True, pct_priors=0, seed=seed, verbose=False)
                #print errd1['corr'].mean()
                print errd1['B_mse'].mean()
                aupr_array[i,j]+= errd1['B_mse'].mean()

    aupr_array /= reps
    ar2 = np.zeros((aupr_array.shape[0],aupr_array.shape[1]))
    for i in range(len(aupr_array)):
        ar2[i,:] = aupr_array[i,:]/aupr_array[i,0]
        #ar2[i,:] = aupr_array[i,:]
    df = pd.DataFrame(ar2[1:,:],index=pct_fused[1:],columns=lamSs)
    df.index.name = 'percent fused'
    df.columns.name = 'lamS'
    plt.figure()
    sns.heatmap(df,cmap="Blues", square=True)
    plt.title('3')
    plt.show()


def L2fusion_quick2():    
    N_TF = 20
    N_G = 200
    pct_fused = list(np.linspace(0.3,1.0,10))
    reps = 1
    
    lamR = 2
    lamSs = list(np.linspace(0,5,10)) 
    lamP = 1.0
    aupr_array = np.zeros((len(pct_fused),len(lamSs)))

    out1 = os.path.join('data','fake_data','L2fusion_quick')
    k = 2#cv folds
    if not os.path.exists(out1):
        os.mkdir(out1)

    for p in range(reps):
        for i, N in enumerate(pct_fused):
            out = os.path.join(out1,'dat_'+str(N))
            ds.write_fake_data1(N1=10, N2=200, out_dir = out, tfg_count1=(N_TF, N_G), tfg_count2 = (N_TF, N_G), measure_noise1 = 0.1, measure_noise2 = 0.1, sparse=0.75, fuse_std = 0, pct_fused=N, orth_falsepos = 0)        
            lamP = 1.0 #priors don't matter
            seed = 10
            for j, lamS in enumerate(lamSs):
                (errd1, errd2) = fg.cv_model_m(out, lamP, lamR, lamS, k, solver='solve_ortho_direct', reverse=True, cv_both=(True,True), exclude_tfs=True, pct_priors=0, seed=seed, verbose=False)
                #print errd1['corr'].mean()
                print errd1['B_mse'].mean()
                aupr_array[i,j]+= errd1['B_mse'].mean()

    aupr_array /= reps
    ar2 = np.zeros((aupr_array.shape[0],aupr_array.shape[1]))
    for i in range(len(aupr_array)):
        ar2[i,:] = aupr_array[i,:]/aupr_array[i,0]
        #ar2[i,:] = aupr_array[i,:]
    df = pd.DataFrame(ar2[1:,:],index=pct_fused[1:],columns=lamSs)
    df.index.name = 'percent fused'
    df.columns.name = 'lamS'
    plt.figure()
    sns.heatmap(df,cmap="Blues", square=True)
    plt.title('1')
    #lt.show()

    lamSs = list(np.linspace(0,10,10)) 

    out3 = os.path.join('data','fake_data','L2fusion_quick3')
    if not os.path.isdir(out3):
        os.mkdir(out3)
    k = 2#cv folds

    for p in range(reps):
        for i, N in enumerate(pct_fused):
            out = os.path.join(out3,'dat_'+str(N))
            if not os.path.isdir(out):
                os.mkdir(out)
            ds.write_fake_data1(N1=10, N2=200, out_dir = out, tfg_count1=(N_TF, N_G), tfg_count2 = (N_TF, N_G), measure_noise1 = 0.1, measure_noise2 = 0.1, sparse=0.75, fuse_std = 0.75, pct_fused=N, orth_falsepos = 0)        
            lamP = 1.0 #priors don't matter
            seed = 10
            for j, lamS in enumerate(lamSs):
                (errd1, errd2) = fg.cv_model_m(out, lamP, lamR, lamS, k, solver='solve_ortho_direct', reverse=True, cv_both=(True,True), exclude_tfs=True, pct_priors=0, seed=seed, verbose=False)
                #print errd1['corr'].mean()
                print errd1['B_mse'].mean()
                aupr_array[i,j]+= errd1['B_mse'].mean()

    aupr_array /= reps
    ar2 = np.zeros((aupr_array.shape[0],aupr_array.shape[1]))
    for i in range(len(aupr_array)):
        ar2[i,:] = aupr_array[i,:]/aupr_array[i,0]
        #ar2[i,:] = aupr_array[i,:]
    df = pd.DataFrame(ar2[1:,:],index=pct_fused[1:],columns=lamSs)
    df.index.name = 'percent fused'
    df.columns.name = 'lamS'
    plt.figure()
    sns.heatmap(df,cmap="Blues", square=True)
    plt.title('3')
    plt.show()


def bac_priors():
    e.plot_bacteria_roc(lamP=(0.0368,0.00789), lamR=(0.368,0.0789), lamSs=[0], k=10, metric='prc',savef='bac_priors',normed=False,scad=False, cv_both=(True,False), roc_species=0, orgs=['B_subtilis','B_subtilis_eu'], unfused=False, lamS_opt=None, orth_file=['orth'], pct_priors=0)
    e.plot_bacteria_roc(lamP=(0.0368,0.00789), lamR=(0.368,0.0789), lamSs=[0], k=10, metric='prc',savef='bac_priors1',normed=False,scad=False, cv_both=(True,False), roc_species=0, orgs=['B_subtilis','B_subtilis_eu'], unfused=False, lamS_opt=None, orth_file=['orth'], pct_priors=.99)

def L2fusion_dataex():    
    N_TF = 20
    N_G = 200
    pct_fused = 0.75
    reps = 1

    N_1 = list(np.linspace(5,50,10))
    N_2 = list(np.linspace(5,50,10))
    
    lamR = 2
    lamS = 10 
    lamP = 1.0
    aupr_array = np.zeros((len(N_1),len(N_2)))

    out1 = os.path.join('data','fake_data','L2fusion_dataex')
    k = 2#cv folds
    if not os.path.exists(out1):
        os.mkdir(out1)

    for p in range(reps):
        for i, n1 in enumerate(N_1):
            for j, n2 in enumerate(N_2):
                out = os.path.join(out1, 'dat_' + str(n1) + str(n2))
                ds.write_fake_data1(N1=n1, N2=n2, out_dir=out, tfg_count1=(N_TF,N_G), tfg_count2=(N_TF,N_G), measure_noise1=0.1, measure_noise2=0.1, sparse=0.75, fuse_std=0.3, pct_fused=pct_fused, orth_falsepos=0)
                lamP=1.0
                seed=10
                (errd1, errd2) = fg.cv_model_m(out, lamP, lamR, lamS, k, solver='solve_ortho_direct', reverse=True, cv_both=(True,True), exclude_tfs=True, pct_priors=0, seed=seed, verbose=False)
                aupr_array[i,j]+= errd1['B_mse'].mean()

    aupr_array /= reps
    ar2 = np.zeros((aupr_array.shape[0],aupr_array.shape[1]))
    for i in range(len(aupr_array)):
    #    ar2[i,:] = aupr_array[i,:]/aupr_array[i,0]
        ar2[i,:] = aupr_array[i,:]
    #return ar2
    df = pd.DataFrame(ar2,index=N_2,columns=N_1)
    df.index.name = 'conditions, network 2'
    df.columns.name = 'conditions, network 1'
    plt.figure()
    sns.heatmap(df,cmap="Blues", square=True)
    #plt.title('')
    plt.show()



def adapt_perc():
    if not os.path.exists(os.path.join('data','fake_data','a')):
        os.mkdir(os.path.join('data','fake_data','a'))

    out = os.path.join('data','fake_data','a')

    N_TF = 35
    N_G = 200
    amt_fused = 1.0
    orth_err = 0.3
    lamS =10
    ds.write_fake_data1(N1 = 20, N2 = 20, out_dir = out, tfg_count1=(N_TF, N_G), tfg_count2 = (N_TF, N_G), pct_fused = amt_fused, orth_falsepos = orth_err, orth_falseneg = orth_err, measure_noise1 = 0.1, measure_noise2 = 0.1, sparse=0.5, fuse_std = 0.1)
    ds1 = ds.standard_source(out,0)
    ds2 = ds.standard_source(out,1)
    orth_fn=os.path.join(out, 'orth')
    organisms=[ds1.name, ds2.name]
    orth = ds.load_orth(orth_fn, organisms)
    (e1, t1, genes1, tfs1) = ds1.load_data()
    (e2, t2, genes2, tfs2) = ds2.load_data()
    Xs = [t1, t2]
    Ys = [e1, e2]
    (priors1, signs1) = ds1.get_priors()
    (priors2, signs2) = ds2.get_priors()
    priors = priors1 + priors2
    genes = [genes1, genes2]
    tfs = [tfs1, tfs2]
    (constraints, marks) = fr.orth_to_constraints_marked(organisms, genes, tfs, orth, lamS)
    betafile1 = os.path.join(out, 'beta1')
    betafile2 = os.path.join(out, 'beta2')
    percs = np.linspace(50,95,5)

    lamP = 1
    lamR = 1
    lamS = 0
    Bs_uf = fr.solve_ortho_direct(organisms, genes, tfs, Xs, Ys, orth, priors, lamP, lamR, lamS)
    lamS = 10
    #plt.hold(True)
    for i, N in enumerate(percs):

        scad_settings = fr.get_settings({'s_it':5, 'per':N, 'return_cons':True})

        lamP = 1
        lamR = 1

        Bs_fs = fr.solve_ortho_direct_scad(organisms, genes, tfs, Xs, Ys, orth, priors, lamP, lamR, lamS, scad_settings)

        false = []
        true = []
        cons = []

        for i in range(len(constraints)):#[0:subs]:
            con = constraints[i]
            cons.append(con.lam)
            mark = marks[i]

            con_fs = scad_settings['cons'][i]
            lams = con_fs.lam**2

            if mark == True:
                true.append(lams)
            if mark == False:
                false.append(lams)

        plt.figure(i)
        plt.hist(true, label='true', alpha = 0.5)
        plt.hist(false, label='false',alpha=0.5)
        plt.title('percentile: %s' % N)
        plt.legend()
        plt.show()
        #plt.show(block=False)

    (B1, genes_1, tfs_1) = ds.load_network(betafile1)
    (B2, genes_2, tfs_2) = ds.load_network(betafile1)
    B_true = []
    B_false = []
    for i in range(len(constraints)):
        con = constraints[i]
        mark = marks[i]
        if mark == True:
            B_true.append(abs(B1[con.c1.r,con.c1.c]-B2[con.c2.r,con.c2.c]))
        if mark == False:
            B_false.append(abs(B1[con.c1.r,con.c1.c]-B2[con.c2.r,con.c2.c]))
    bins = np.linspace(0,5,20)
    plt.hist(B_true, bins, histtype='stepfilled',label='true', alpha=0.5)
    plt.hist(B_false, bins, histtype = 'stepfilled',label='false',alpha=0.5)
    plt.legend()
    plt.show()

    bactf = out
    ds1 = ds.standard_source(bactf,0)
    ds2 = ds.standard_source(bactf,1)
    (priors1, signs1) = ds1.get_priors()
    
    (priors2, signs2) = ds2.get_priors()
    (constraints, marks, orths) = ds.load_constraints(bactf)
    (e1_tr, t1_tr, genes1, tfs1) = ds1.load_data()
    (e2_tr, t2_tr, genes2, tfs2) = ds2.load_data()
    
    tfs_set_subt = set(tfs1)
    tfs_set_anth = set(tfs2)

    orth_set_subt = set(map(lambda orth: orth.genes[0].name, orths))
    orth_set_anth = set(map(lambda orth: orth.genes[1].name, orths))
    #map from one to the other, OK if orthology 1-1
    subt_to_anth = {orths[x].genes[0].name : orths[x].genes[1].name for x in range(len(orths))}
    anth_to_subt = {orths[x].genes[1].name : orths[x].genes[0].name for x in range(len(orths))}


    tfs_orth_subt = filter(lambda tf: tf in orth_set_subt and subt_to_anth[tf] in tfs_set_anth, tfs1)
    tfs_orth_anth = filter(lambda tf: tf in orth_set_anth and anth_to_subt[tf] in tfs_set_subt, tfs2)

    gen_orth_subt = filter(lambda g: g in orth_set_subt, genes1)
    gen_orth_anth = filter(lambda g: g in orth_set_anth, genes2)

    gene_inds_subt = {genes1[i] : i for i in range(len(genes1))}
    gene_inds_anth = {genes2[i] : i for i in range(len(genes2))}

    subt_corr_mat = np.zeros((len(tfs_orth_subt), len(gen_orth_subt)))
    anth_corr_mat = np.zeros((len(tfs_orth_anth), len(gen_orth_anth)))
    subt_B_mat = np.zeros((len(tfs_orth_subt), len(gen_orth_subt)))
    anth_B_mat = np.zeros((len(tfs_orth_anth), len(gen_orth_anth)))
    subt_Breal_mat = np.zeros((len(tfs_orth_subt), len(gen_orth_subt)))
    anth_Breal_mat = np.zeros((len(tfs_orth_anth), len(gen_orth_anth)))
   #subt_corr_mat = np.random.randn(50,50)#len(tfs_orth_subt), len(gen_orth_subt))
    #anth_corr_mat = np.random.randn(50,50)#len(tfs_orth_anth), len(gen_orth_anth))
    
    
    for r, tf in enumerate(tfs_orth_subt):
        for c, g in enumerate(gen_orth_subt):
            
            tfi = gene_inds_subt[tf]
            gi = gene_inds_subt[g]
            corr = np.corrcoef( np.hstack((e1_tr[:, [tfi]], e1_tr[:, [gi]])).T )[0,1]            
            subt_corr_mat[r, c] = corr
            subt_B_mat[r, c] = Bs_uf[0][tfi, gi]
            subt_Breal_mat[r, c] = B1[tfi, gi]

    anth_to_subt = {orths[x].genes[1].name : orths[x].genes[0].name for x in range(len(orths))}
    #for the anthracis correlations, we map the tf or gene onto the corresponding subtilis gene, then compute that index
    tfs_orth_subt_ind = {tfs_orth_subt[i] : i for i in range(len(tfs_orth_subt))}
    gen_orth_subt_ind = {gen_orth_subt[i] : i for i in range(len(gen_orth_subt))}

    for r, tf in enumerate(tfs_orth_anth):
        for c, g in enumerate(gen_orth_anth):
            tfi = gene_inds_anth[tf]
            gi = gene_inds_anth[g]
            corr = np.corrcoef( np.hstack((e2_tr[:, [tfi]], e2_tr[:, [gi]])).T )[0, 1]

            tf_subt = anth_to_subt[tf]
            g_subt = anth_to_subt[g]
       
            r_subt = tfs_orth_subt_ind[tf_subt]
            c_subt = gen_orth_subt_ind[g_subt]
            
            anth_corr_mat[r_subt, c_subt] = corr
            anth_B_mat[r_subt, c_subt] = Bs_uf[1][tfi,gi]
            anth_Breal_mat[r_subt, c_subt] = B2[1][tfi,gi] 
            
    subt_corrs = subt_corr_mat.ravel()
    anth_corrs = anth_corr_mat.ravel()
    subt_B = subt_B_mat.ravel()
    anth_B = anth_B_mat.ravel()
    subt_B_real = subt_Breal_mat.ravel()
    anth_B_real = anth_Breal_mat.ravel()
    
    plt.scatter(np.abs(subt_B - anth_B), np.abs(subt_corrs - anth_corrs), alpha=0.5)
    plt.xlabel('delta beta, unfused')
    plt.ylabel('delta correlation')
    plt.show()
    plt.scatter(np.abs(subt_B_real - anth_B_real), np.abs(subt_corrs - anth_corrs), alpha=0.5)
    plt.xlabel('delta beta, true')
    plt.ylabel('delta correlation')
    plt.show()

    sns.kdeplot(np.abs(subt_corrs - anth_corrs), shade=True, label = 'orthologs, %f'% np.mean(np.abs(subt_corrs - anth_corrs)))
    plt.hold(True)
    np.random.shuffle(subt_corrs)
    sns.kdeplot(np.abs(subt_corrs - anth_corrs), shade=True, label = 'non-orths, %f'% np.mean(np.abs(subt_corrs - anth_corrs)))
    plt.xlabel('correlation difference')
    plt.ylabel('frequency')
    plt.show()


def bsub_fakeorth(falsepos=0.3,use_TFA=True):
    out = os.path.join('data','bacteria_standard')
    fil = file(os.path.join(out, 'orth'))
    fs = fil.read()
    fsn = filter(len, fs.split('\n'))
    fsnt = map(lambda x: x.split('\t'), fsn)    
    orths = []

    organisms = ['B_subtilis','B_anthracis']
    for o in fsnt:
        real = True
        if len(o[0]) > 0:
            if len(o[1]) > 1:
                orth = fr.orthology(genes = (fr.one_gene(name=o[1],organism=organisms[1]), fr.one_gene(name=o[0], organism=organisms[0])), real = real)
        orths.append(orth)

    falseneg=0

    ds_sub = ds.standard_source(out,0, use_TFA=use_TFA)
    ds_ant = ds.standard_source(out,1, use_TFA=use_TFA)
    (exp_sub, tfexp_sub, genes_sub, tfs_sub) = ds_sub.load_data()
    (exp_ant, tfexp_ant, genes_ant, tfs_ant) = ds_ant.load_data()
    priors_sub = ds_sub.get_priors()[0]
    priors_ant = ds_ant.get_priors()[0]

    tfs_sub2 = map(lambda x: x, tfs_sub)
    genes_sub2 = map(lambda x: x, genes_sub)
    tfs_ant2 = map(lambda x: x, tfs_ant)
    genes_ant2 = map(lambda x: x, genes_ant)

    new_orths = ds.generate_faulty_orth3(orths, genes_sub2, tfs_sub2, genes_ant2, tfs_ant2, organisms, falsepos, falseneg)
    ds.write_orth(out+os.sep+'fake_orth', new_orths, organisms)

    seed = np.random.randn()
    settings = fr.get_settings()
    solver = 'solve_ortho_direct_scad'
    N = 20
    scad_settings = fr.get_settings({'s_it':5, 'per':N, 'return_cons':True})
    genes = [genes_sub, genes_ant]
    tfs = [tfs_sub, tfs_ant]
    Xs = [tfexp_sub, tfexp_ant]
    Ys = [exp_sub, exp_ant]
    priors = []
    (constraints, marks, orths) = ds.load_constraints(out, 'fake_orth', ['B_subtilis','B_anthracis'])

    constraint_dict=dict()
    lamS=0
    lamP=1
    lamR=0.5
    Bs_uf = fr.solve_ortho_direct(organisms, genes, tfs, Xs, Ys, new_orths, priors, lamP, lamR, lamS)

    false_constraints_dbs = []
    true_constraints_dbs = []
    for i in range(len(constraints)):
        con = constraints[i]
        mark = marks[i]
        if mark == False:
            #constraint_dict[con] = abs(Bs_uf[con.c1.sub][con.c1.r, con.c1.c] - Bs_uf[con.c2.sub][con.c2.r,con.c2.c])
            false_constraints_dbs.append(abs(Bs_uf[con.c1.sub][con.c1.r, con.c1.c] - Bs_uf[con.c2.sub][con.c2.r,con.c2.c]))
        if mark == True:
            true_constraints_dbs. append(abs(Bs_uf[con.c1.sub][con.c1.r, con.c1.c] - Bs_uf[con.c2.sub][con.c2.r,con.c2.c]))

    lamS = 1
    Bs_fs = fr.solve_ortho_direct_scad(organisms, genes, tfs, Xs, Ys, new_orths, priors, lamP, lamR, lamS, scad_settings)

    false = []
    true = []
    true_count = 0
    true_unfused = 0
    false_count = 0
    false_unfused = 0
    cons = []

    for i in range(len(constraints)):#[0:subs]:
        con = constraints[i]
        cons.append(con.lam)
        mark = marks[i]

        con_fs = scad_settings['cons'][i]
        lams = con_fs.lam**2

        if mark == True:
            true.append(lams)
            true_count +=1
            if lams == 0:
                true_unfused +=1
        if mark == False:
            false.append(lams)
            false_count +=1
            if lams == 0:
                false_unfused +=1

    fun = false_unfused / float(false_count)
    tun = true_unfused / float(true_count)

    plt.bar([0,1], [fun,tun])

    bins = np.linspace(0,8,13)
    plt.hist(false_constraints_dbs, bins, alpha=0.5)
    plt.hist(true_constraints_dbs, bins)

    (a1,b1,c1)=plt.hist(true_constraints_dbs, bins, alpha=0.5)
    (a2,b2,c2)=plt.hist(false_constraints_dbs, bins, alpha=0.5)
    plt.plot((b1[0:-1] + b1[1:])/2, a1/a1.sum())
    plt.plot((b2[0:-1] + b2[1:])/2, a2/a2.sum())