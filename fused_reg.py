import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt
import pandas as pd
import random
import collections
import time
import scipy.sparse
import scipy.sparse.linalg
from sklearn import mixture
try:
    from glmnetpython import ElasticNet
    from glmnetpython import CVGlmNet
except ImportError:
    print 'glmnetpython error'

try:
    import rpy2
    import rpy2.robjects.packages as rpackages
    import rpy2.robjects as robjects
    import rpy2.robjects.numpy2ri
except ImportError:
    print 'rpy2 not installed'



from sklearn import linear_model

def get_settings(override = None):
    d = dict()
    d['return_cons'] = False
    d['cons'] = []
    d['adjust'] = False
    d['solver'] = 'basic'#use this later
    d['a'] = 0
    d['per'] = 0
    d['s_it'] = 2
    d['m_it'] = 5
    d['it'] = 100 #number of iterative solver iterations (lamS steps)
    d['lamR_steps'] = 100
    d['sol_changes'] = [] #for iterative solver, how fast the solution changes
    d['iter_eval'] = False #for iterative solver, evaluation function for scoring lambdaS paths. If true, will get filled in as appropriate
    d['lamS_path'] = [] #gets filled in if iter_eval is supplied
    d['lamS_scores'] = []
    if override != None:
        for k in override.keys():
            if k in d:
                d[k] = override[k]
            else:
                print 'attempting to override nonexistent key'
                print d[k]
    return d    

#SECTION: -------------------DATA STRUCTURES--------------
coefficient = collections.namedtuple('coefficient',['sub','r','c'])

#c2 == None is a constant constraint
constraint = collections.namedtuple('constraint',['c1','c2','lam'])

one_gene = collections.namedtuple('gene',['name','organism'])

orthology = collections.namedtuple('orthology',['genes','real'])

column = collections.namedtuple('column',['sub','c'])
class dictl(dict):
    def __getitem__(self, i):
        if not i in self:
            self[i] = []
        return super(dictl, self).__getitem__(i)

#SECTION: ------------------------UTILITY FUNCTIONS-------
#diagonally concatenates two matrices
#TODO: sparse
def diag_concat(mats):

    tot_rows=sum(map(lambda x: x.shape[0], mats))
    tot_cols=sum(map(lambda x: x.shape[1], mats))
    A=np.zeros((tot_rows, tot_cols))
    row=0
    col=0
    for i in range(0,len(mats)):
        mi = mats[i]
        A[row:(row+mi.shape[0]), col:(col+mi.shape[1])] = mi
        row += mi.shape[0]
        col += mi.shape[1]
    return A

def diag_concat_sparse(mats):
    tot_rows=sum(map(lambda x: x.shape[0], mats))
    tot_cols=sum(map(lambda x: x.shape[1], mats))
    #csr_matrix((data, (row_ind, col_ind)), [shape=(M, N)])
    #where data, row_ind and col_ind satisfy the relationship a[row_ind[k], col_ind[k]] = data[k].
    r_inds_l = []
    c_inds_l = []
    data_l = []
    r_offs = 0
    c_offs = 0
    
    for m in mats:
        data = m.ravel()
        inds = np.arange(m.shape[0]*m.shape[1])
        (r_inds, c_inds) = np.unravel_index(inds, m.shape)
        r_inds = np.array(r_inds) #make writeable
        c_inds = np.array(c_inds)
        
        r_inds += r_offs
        c_inds += c_offs
        r_offs += m.shape[0]
        c_offs += m.shape[1]
        data_l.append(data)
        r_inds_l.append(r_inds)
        c_inds_l.append(c_inds)
    data_all = np.hstack(data_l)
    r_inds_all = np.hstack(r_inds_l)
    c_inds_all = np.hstack(c_inds_l)
    
    M = scipy.sparse.csr_matrix((data_all, (r_inds_all, c_inds_all)), shape=(tot_rows, tot_cols))
    return M

#returns a list of lists of all pairs of entries from l, where the first entry in the pair occurs before the second in l
def all_pairs(l):
    pl = []
    for li1 in range(len(l)):
        for li2 in range(li1+1, len(l)):
           pl.append([l[li1], l[li2]])
    return pl

#solves a quadratic equation, returning both solutions
def quad(a, b, c):
    x1 = 0.5*(-b + (b**2 - 4*a*c)**0.5)/a
    x2 = 0.5*(-b - (b**2 - 4*a*c)**0.5)/a
    
    return (x1, x2)



#SECTION: --------CODE FOR GENERATING CONSTRAINTS---------

#given a list of orthology objects, returns a (larger) list of one-to-one orthologies. that is, orthologies with only two genes, consisting of all pairings within each orth
def orth_to_one_one_orth(orths):
    orths2 = []
    for orth in orths:
        gene_pairs = all_pairs(orth.genes)
        for gene_pair in gene_pairs:
            new_orth = orthology(genes=(gene_pair[0], gene_pair[1]), real=orth.real)
            orths2.append(new_orth)
    return orths2
#enumerates fusion constraints from orthology
#args as in solve_ortho_direct
#returns list of constraints. individual constraints are pairs of coefficients, with an associated weight lamS
#also returns a list of True/False for real or fake orthologies
#lamS_opt: specify a dictionary mapping tuples (pairs) of organisms to lamS for that organism
def orth_to_constraints_marked(organisms, gene_ls, tf_ls, orth, lamS, lamS_opt=None):
    #build some maps
    
    org_to_ind = {organisms[x] : x for x in range(len(organisms))}
    gene_to_inds = map(lambda o: {gene_ls[o][x] : x for x in range(len(gene_ls[o]))}, range(len(organisms)))
    tf_to_inds = map(lambda o: {tf_ls[o][x] : x for x in range(len(tf_ls[o]))}, range(len(organisms)))
    
    #turn orth into a list of all pairs of coefficients
    #orth_pairs = reduce(lambda x,y: x+y, map(all_pairs, orth),[])
    orth_pairs = orth_to_one_one_orth(orth)

    #build a set to use for checking on the realness of an orthology
    
    real_orths = set()
    for orth in orth_pairs:
        if orth.real:
            real_orths.add((orth.genes[0], orth.genes[1]))
            real_orths.add((orth.genes[1], orth.genes[0]))

    #build a dictionary mapping each gene to all of its orthologs
    #the orthology file can be messy, and specify the same orthology multiple times, so changing from a list of orthologs to a set of orthologs
    
    ortho_dict = collections.defaultdict(lambda: set())#dictl()
    for op in orth_pairs:
        (gene1, gene2) = op.genes
        
        ortho_dict[gene1].add(gene2)
        ortho_dict[gene2].add(gene1)
    #for the purposes of operons, it's nice if every gene is an ortholog of itself
    for org_i in range(len(organisms)):
        for g_i in range(len(gene_ls[org_i])):
            g = one_gene(gene_ls[org_i][g_i], organisms[org_i])
            #ortho_dict[g].add(g)

    constraints = []
    marks = []
    for org_i in range(len(organisms)):
        for tf_i in range(len(tf_ls[org_i])):
            for g_i in range(len(gene_ls[org_i])):
                tf = one_gene(tf_ls[org_i][tf_i], organisms[org_i])
                g = one_gene(gene_ls[org_i][g_i], organisms[org_i])

                for tf_orth in ortho_dict[tf]:
                    lamS_pair = lamS
                    if lamS_opt != None:
                        if (g.organism, tf_orth.organism) in lamS_opt:
                            lamS_pair = lamS_opt[(g.organism, tf_orth.organism)]
                        if (tf_orth.organism, g.organism) in lamS_opt:
                            lamS_pair = lamS_opt[(g.organism, tf_orth.organism)]

                    for g_orth in ortho_dict[g]:
                        
                        sub1 = org_i
                        sub2 = org_to_ind[tf_orth.organism]
                        sub3 = org_to_ind[g_orth.organism]
                        
                        if not sub2 == sub3:
                            continue
                        if not tf_orth.name in tf_to_inds[sub2]:
                            continue #if a tf is orthologous to a non-tf
                        if not tf_orth.name in tf_to_inds[sub2] or not g_orth.name in gene_to_inds[sub2]:
                            continue #orth file might specify genes that we don't have

                        if tf == tf_orth and g == g_orth:
                            continue #no point in self constraints

                        #now check if it's real
                        real = (g, g_orth) in real_orths and (tf, tf_orth) in real_orths

                        coeff1 = coefficient(sub1, tf_to_inds[sub1][tf.name], gene_to_inds[sub1][g.name])
                        coeff2 = coefficient(sub2, tf_to_inds[sub2][tf_orth.name], gene_to_inds[sub2][g_orth.name])
                        
                        constr = constraint(coeff1, coeff2, lamS_pair)
                        constraints.append(constr)
                        marks.append(real)
    return (constraints, marks)

def orth_to_constraints(organisms, gene_ls, tf_ls, orth, lamS, lamS_opt=None):
    (con, marks) = orth_to_constraints_marked(organisms, gene_ls, tf_ls, orth, lamS, lamS_opt=None)
    return con

#enumerates fusion constraints from orthology
#args as in solve_ortho_direct
#returns list of constraints. individual constraints are pairs of coefficients, with an associated weight lamS
def orth_to_constraints_old(organisms, gene_ls, tf_ls, orth, lamS, lamS_opt=None):
    #build some maps
    
    org_to_ind = {organisms[x] : x for x in range(len(organisms))}
    gene_to_inds = map(lambda o: {gene_ls[o][x] : x for x in range(len(gene_ls[o]))}, range(len(organisms)))
    tf_to_inds = map(lambda o: {tf_ls[o][x] : x for x in range(len(tf_ls[o]))}, range(len(organisms)))
    
    #turn orth into a list of all pairs of coefficients
    #orth_pairs = reduce(lambda x,y: x+y, map(all_pairs, orth),[])
    orth_pairs = reduce(lambda x,y: x+y, map(all_pairs, map(lambda z: z.genes, orth)),[]) #we just want the genes that are in the orthology groups
    
    ortho_dict = dictl()
    for op in orth_pairs:
    
        (gene1, gene2) = op
        ortho_dict[gene1].append(gene2)
        ortho_dict[gene2].append(gene1)
    
    constraints = []
    for org_i in range(len(organisms)):
        for tf_i in range(len(tf_ls[org_i])):
            for g_i in range(len(gene_ls[org_i])):
                tf = one_gene(tf_ls[org_i][tf_i], organisms[org_i])
                g = one_gene(gene_ls[org_i][g_i], organisms[org_i])
                
                for tf_orth in ortho_dict[tf]:
                    for g_orth in ortho_dict[g]:
                        sub1 = org_i
                        sub2 = org_to_ind[tf_orth.organism]
                        sub3 = org_to_ind[g_orth.organism]
            
                        if not sub2 == sub3:
                            continue
                        if not tf_orth.name in tf_to_inds[sub2]:
                            continue #if a tf is orthologous to a non-tf
                        if not tf_orth.name in tf_to_inds[sub2] or not g_orth.name in gene_to_inds[sub2]:
                            continue #orth file might specify genes that we don't have                                     
                        coeff1 = coefficient(sub1, tf_to_inds[sub1][tf.name], gene_to_inds[sub1][g.name])
                        coeff2 = coefficient(sub2, tf_to_inds[sub2][tf_orth.name], gene_to_inds[sub2][g_orth.name])
                        

                        if lamS_opt != None:
                            nlamS = 0
                            if (tf_orth.organism, g_orth.organism) in lamS_opt:
                                nlamS = lamS_opt[(tf_orth.organism, g_orth.organism)]
                                print 'new lamS'
                                print nlamS
                            if (g_orth.organism, tf_orth.organism) in lamS_opt:
                                nlamS = lamS_opt[(g_orth.organism, tf_orth.organism)]
                                print 'new lamS'
                                print nlamS
                            else:
                                print 'whats going on'
                            constr = constraint(coeff1, coeff2, nlamS)
                        else:
                            constr = constraint(coeff1, coeff2, lamS)
                        constraints.append(constr)
    
    return constraints

#given priors, returns a list of constraints. constraints representing priors have their second coefficient equal to None.
#NOTE: this function does not allow priors to have different weights.
def priors_to_constraints(organisms, gene_ls, tf_ls, priors, lam):
    
    org_to_ind = {organisms[x] : x for x in range(len(organisms))}
    gene_to_inds = map(lambda o: {gene_ls[o][x] : x for x in range(len(gene_ls[o]))}, range(len(organisms)))
    tf_to_inds = map(lambda o: {tf_ls[o][x] : x for x in range(len(tf_ls[o]))}, range(len(organisms)))
    constraints = []
    
    for prior in priors:
        (gene1, gene2) = prior
        sub1 = org_to_ind[gene1.organism]
        sub2 = org_to_ind[gene2.organism]
        if not sub1 == sub2:
            print '!?!?!?!?!?'
            continue

        tfi = tf_to_inds[sub1][gene1.name]
        gi = gene_to_inds[sub2][gene2.name]
        constr = constraint(coefficient(sub1, tfi, gi), None, lam[sub1])
        constraints.append(constr)
    return constraints

#adds a "prior" on transcription factors not self-regulating. This is an attempt to avoid degenerate learning of an identity matrix for that part of B
def no_self_reg_constraints(organisms, gene_ls, tf_ls, lam):
    constraints = []
    gene_to_inds = map(lambda o: {gene_ls[o][x] : x for x in range(len(gene_ls[o]))}, range(len(organisms)))
    tf_to_inds = map(lambda o: {tf_ls[o][x] : x for x in range(len(tf_ls[o]))}, range(len(organisms)))

    for subi, org in enumerate(organisms):
        for tfi, tf in enumerate(tf_ls[subi]):
            gi = gene_to_inds[subi][tf]
            constr = constraint(coefficient(subi, tfi, gi), None, lam)
            constraints.append(constr)
    return constraints

#SECTION: --------CODE FOR ADJUSTING CONSTRAINTS---------

#for the given columns, figures out the inverse covariance matrix of the equivalent prior for fused/unfused, takes the trace of the inverses as the variance, and figures out a constant to multiply fused interactions by that will equalize them
def adjust_var(Xs, cols, fuse_cons, ridge_cons, lamR):
#returns the inverse covariance matrix of the prior equivalent to the given constraint structure
    counter = 0
    b_to_n = dict()
    n_to_b = dict()    
    for co in cols:
        sub = co.sub
        for r in range(Xs[sub].shape[1]):
            b_to_n[(sub, co.c, r)] = counter
            n_to_b[counter] = (sub, co.c, r)
            counter += 1
    N = counter

    def inv_cov(Xs, cols, fuse_cons, ridge_cons, lamR):
        Einv = np.zeros((N, N))
        for n in range(N):
            Einv[n, n] = lamR
        for ridge_con in ridge_cons:
            n = b_to_n[(ridge_con.c1.sub, ridge_con.c1.c, ridge_con.c1.r)]
            Einv[n, n] = ridge_con.lam
        for fuse_con in fuse_cons:
            n = b_to_n[(fuse_con.c1.sub, fuse_con.c1.c, fuse_con.c1.r)]
            m = b_to_n[(fuse_con.c2.sub, fuse_con.c2.c, fuse_con.c2.r)]
            Einv[n, m] += -fuse_con.lam
            Einv[m, n] += -fuse_con.lam
            Einv[n, n] += fuse_con.lam
            Einv[m, m] += fuse_con.lam
        return Einv

    Einv_fused = inv_cov(Xs, cols, fuse_cons, ridge_cons, lamR)
    Einv_unfused = inv_cov(Xs, cols, [], ridge_cons, lamR)

    
    v_f = np.sum(np.trace(np.linalg.inv(Einv_fused)))
    v_u = np.sum(np.trace(np.linalg.inv(Einv_unfused)))

    #to fix the cov mat ratio we want to multiply by v_u/v_f
    #because we're multiplying by the inverse, v_f/v_u
    c = v_f/v_u
    fuse_cons_adj = []
    ridge_cons_adj = []
    lamR = lamR*c
    #print 'living with %d constraints' % len(fuse_cons)
    #print 'adjusting by (%f / %f) = %f'%(v_f, v_u, c)
    for fcon in fuse_cons:
        fuse_cons_adj.append(constraint(fcon.c1, fcon.c2, fcon.lam*c))
    for rcon in ridge_cons:
        ridge_cons_adj.append(constraint(rcon.c1, rcon.c2, rcon.lam*c))
        
    return (fuse_cons_adj, ridge_cons_adj, lamR)

#SECTION: ----------------CODE FOR SOLVING THE MODEL-------------------

#solves for the optimal lamR value for Xs, Ys. maxlamR is the max lamR value tried; lamR steps is the number of lamR values between 0 and maxlamR tried. glmnet fits models for each value of lambda for each cv fold, then chooses the best model by minimizing deviance, which is minus twice the log-likelihood on the left-out data.
#returns list of optimal lamRs for each iteration; pandas dataframe with lamR and deviance; list of lambdaRs used
def opt_lamR(Xs, Ys, folds, maxlamR, lamR_steps, it):
    lambdaRs = np.linspace(0, maxlamR, lamR_steps)
    optlamR = []

    #tot_devs is a list of dicts, one for each species, mapping lamR to the deviances
    tot_devs = map(lambda y: collections.defaultdict(list), range(len(Xs)))
    #tot_dev=collections.defaultdict(list)
    #dev_unit=collections.defaultdict(list)

    for i in range(it):
        for j in range(len(Xs)):
            for k in range(Ys[j].shape[1]):
                print k
                enet=ElasticNet(alpha=0)
                enet_cv = CVGlmNet(enet, n_folds=folds)
                enet_cv.fit(Xs[j],Ys[j][:,k], lambdas=lambdaRs)
                #tot_dev.append(enet_cv._oof_deviances)
                dev = enet_cv._oof_deviances
                for m in range(len(lambdaRs)):
                    unit = len(tot_devs[j][lambdaRs[m]])
                    tot_devs[j][lambdaRs[m]].append((dev[m],unit))


    all_devs = []
    for i in range(len(Xs)):
        sorted_devitems = sorted(tot_devs[i].items(), key= lambda x:x[0])
        sorted_devvals = map(lambda x: x[1], sorted_devitems)
        summed_devs = map(sum, map(lambda y: map(lambda x: x[0], y), sorted_devvals))
        all_devs.append(summed_devs)

    for i in range(len(all_devs)):
        opt_ind = np.where(all_devs[i] == min(all_devs[i]))
        print opt_ind
        optlamR.append(lambdaRs[opt_ind])


    devs = map(lambda tot_dev: pd.DataFrame([[lamR, deviance, iteration] for lamR, devlist in tot_dev.items() for deviance, iteration in devlist], columns=['lamR', 'deviance','unit']), tot_devs)  

    return (optlamR, devs, lambdaRs)


def lstsq_dumb(A, b):
    return np.dot(np.linalg.pinv(A), b)

#iterative solver
def iter_solve(Xs, Ys, fuse_constraints, ridge_constraints, lambdaR, settings):
    it = settings['it']
    iter_eval = settings['iter_eval']
    #from glmnetpython import ElasticNet
    Bs = []
    #set the initial guess (all zero)
    for j in range(len(Xs)):
        Bs.append(np.zeros((Xs[j].shape[1], Ys[j].shape[1])))
    cB = Bs
    pB = map(lambda x: x.copy(), Bs)
    lamS_path = []
    lam_ramp = np.linspace(0, 1.0, it)
    sol_changes = []
    scores = []
    #build a map from subproblem/column to constraint
    sc_to_con = collections.defaultdict(lambda: [])
    for con in fuse_constraints:
        sc_to_con[(con.c1.sub, con.c1.c)].append(con)
            
    #builds a penalty matrix associated with column c of subproblem s
    def build_pen_targ(s, c):
        P = np.zeros((len(sc_to_con[s,c]), Bs[s].shape[0]))
        targ = []
        for i, con in enumerate(sc_to_con[s, c]):
            r_fr = con.c1.r
            r_to = con.c2.r
            P[i, r_fr] = con.lam
            targ.append((con.c2.sub, con.c2.r, con.c2.c, con.lam))
        return (P, targ)

    #given a list of targets as returned by build_pen_targ, returns vector of targeted values
    def targ_to_yp(targ, lamsc):
        yp = np.zeros((len(targ), 1))
        for i, (s, r, c, lam) in enumerate(targ):
            yp[i, 0] = (lam*lamsc)**0.5 * Bs[s][r,c]
        return yp

    if False:
        P_targ_saved = []
        for s in range(len(Xs)):
            P_targ_saved.append([])
            for c in range(Ys[s].shape[1]):
                P_targ_saved[-1].append(build_pen_targ(s, c))
        
    for i in range(it):
        #first we swap cB and pB
        tmp = pB
        pB = cB
        cB = tmp
        for s in range(len(Xs)):
            X = Xs[s] #X, Y for current subproblem s
            Y = Ys[s]
            #print 'here1'
            lamsc = lam_ramp[i]
            
            for c in range(Y.shape[1]): #column of B we are solving
                print 'iteration %d, lamS=%f, sub %d, col %d'%(i,lam_ramp[i],s,c),
                print '\r',

                I = np.eye(X.shape[1])*lambdaR[s]
                (P, targ) = build_pen_targ(s,c)
                #(P, targ) = P_targ_saved[s][c]
                P_cat = (P*lamsc)**0.5

                yp = targ_to_yp(targ, lamsc)
                #print P*lamsc
                F = np.vstack((X, I, P_cat))
                
                yI = np.zeros((I.shape[0], 1))
                y = np.vstack((Y[:, [c]], yI, yp))
                                
                (b, resid, rank, sing) = np.linalg.lstsq(F, y)
                #b=lstsq_dumb(F, y)
                
                cB[s][:, [c]] = b
            sol_change = ((pB[s] - cB[s])**2).sum()
            sol_changes.append(sol_change)
        if iter_eval:
            scores.append(iter_eval(cB))
            print 'iteration %d: score: %f'%(i, scores[-1])
    #backchannel return sol_changes
    settings['sol_changes'].append(sol_changes)
    settings['lamS_path'].append(lam_ramp)
    settings['lamS_scores'].append(scores)
    return cB


#no cleverness at all
#this is here as a sanity check
def direct_solve(Xs, Ys, fuse_constraints, ridge_constraints, lambdaR):
    ncols = sum(map(lambda y: y.shape[1], Ys))
    nrowsX = sum(map(lambda x: x.shape[0], Xs))
    ncolsX = sum(map(lambda x: x.shape[1], Xs))
    X = diag_concat(Xs * ncols)
    y = diag_concat(Ys).ravel(order='F')[:, np.newaxis]
    #y = np.vstack(map(lambda x: x.ravel()[:, np.newaxis], Ys)*ncols)
    
    #finds the flattened index in b of subproblem sub, row r, col c
    def indb(sub, r, c):
        ind = ncols * c + r
    xpad_l = [X]

    ypad_l = [y]
    I = np.eye(X.shape[1])*lambdaR
    for con in ridge_constraints:
        ind = indb(con.c1.sub, con.c1.r, con.c1.c)
        I[ind, ind] = con.lam
    for con in fuse_constraints:
        xpad_l.append(np.zeros((1, X.shape[1])))
        ypad_l.append(0)
        ind1 = indb(con.c1.sub, con.c1.r, con.c1.c)
        ind2 = indb(con.c2.sub, con.c2.r, con.c2.c)
        xpad_l[-1][0, ind1] = con.lam
        xpad_l[-1][0, ind2] = -con.lam
    xpad_l.append(I)
    ypad_l.append(np.zeros((X.shape[1], 1)))
    
    F = np.vstack(xpad_l)
    p = np.vstack(ypad_l)
    
    (b, resid, rank, sing) = np.linalg.lstsq(F, p)    

    B = np.reshape(b, (ncolsX, ncols),order='F')

        #now collect Bs
    Bs = []
    crow = 0
    ccol = 0
    
    for i in range(len(Xs)):
        extr = Xs[i].shape[1]
        extc = Ys[i].shape[1]
        Bs.append(B[crow:(crow + extr), ccol:(ccol+extc)])
    return Bs


#reference solver, that does not involve breaking up constraints
#gene_ls/tf_ls: lists of gene names and tf names for each problem
#Xs: list of TF expression matrices
#YS: list of gene expression matrices
#Orth: list of lists of one_genes
#priors: list of lists of one_gene pairs
def solve_ortho_ref(organisms, gene_ls, tf_ls, Xs, Ys, orth, priors,lamP, lamR, lamS, lamS_opt, settings=None):
    if settings == None:
        settings = get_settings()
    (lamR, lamP) = pad_lamRP(lamR, lamP, len(organisms))    
    ridge_con = priors_to_constraints(organisms, gene_ls, tf_ls, priors, lamP)
    fuse_con = orth_to_constraints(organisms, gene_ls, tf_ls, orth, lamS, lamS_opt)
    Bs = direct_solve(Xs, Ys, fuse_con, ridge_con, lamR, adjust = settings['adjust'])
    return Bs

#most basic solver
#gene_ls/tf_ls: lists of gene names and tf names for each problem
#Xs: list of TF expression matrices
#YS: list of gene expression matrices
#Orth: list of lists of one_genes
#priors: list of lists of one_gene pairs
#lamR: list of lambdaR for each species
def solve_ortho_direct(organisms, gene_ls, tf_ls, Xs, Ys, orth, priors,lamP, lamR, lamS, lamS_opt=None, adjust=False, self_reg_pen = 0, settings=None):   
    if settings == None:
        settings = get_settings()    
    (lamR, lamP) = pad_lamRP(lamR, lamP, len(organisms))    
    ridge_con = priors_to_constraints(organisms, gene_ls, tf_ls, priors, lamP)
    fuse_con = orth_to_constraints(organisms, gene_ls, tf_ls, orth, lamS**0.5, lamS_opt)
    print 'there are %d fuse cons'%len(fuse_con)
    Bs = direct_solve_factor(Xs, Ys, fuse_con, ridge_con, lamR, adjust = settings['adjust'])    
    
    return Bs

#makes sure that lamR/lamS are lists of length N.
#if lamP is supplied as a scalar, then lamP is assumed to represent a constant multiple of lamR. Otherwise, lamP is the actual lamP value
def pad_lamRP(lamR, lamP, N):
    try:
        lamR[0]
    except:
        lamR = [lamR] * N
    try:
        lamP[0]
    except:
        lamP = [lamP * lamR[i] for i in range(N)]
    return (lamR, lamP)

#iterative solver
#gene_ls/tf_ls: lists of gene names and tf names for each problem
#Xs: list of TF expression matrices
#YS: list of gene expression matrices
#Orth: list of lists of one_genes
#priors: list of lists of one_gene pairs
#CURRENTLY IGNORING lamS_opt
def solve_ortho_iter(organisms, gene_ls, tf_ls, Xs, Ys, orth, priors,lamP, lamR, lamS, lamS_opt=None, adjust=False, self_reg_pen = 0, settings=None):   
    if settings == None:
        settings = get_settings()    
    
    (lamR, lamP) = pad_lamRP(lamR, lamP, len(organisms))    
    

    ridge_con = priors_to_constraints(organisms, gene_ls, tf_ls, priors, lamP)
    fuse_con = orth_to_constraints(organisms, gene_ls, tf_ls, orth, lamS, lamS_opt)

    #actually pass settings to this one
    Bs = iter_solve(Xs, Ys, fuse_con, ridge_con, lamR, settings)
    
    return Bs

#parameters as solve_ortho_direct. 
#s_it defines the number of scad-like iterations to do
def solve_ortho_direct_scad(organisms, gene_ls, tf_ls, Xs, Ys, orth, priors, lamP, lamR, lamS, settings = None):
    if settings == None:
        settings = get_settings()  
    (lamR, lamP) = pad_lamRP(lamR, lamP, len(organisms))    
    ridge_con = priors_to_constraints(organisms, gene_ls, tf_ls, priors, lamP)    
    fuse_con = orth_to_constraints(organisms, gene_ls, tf_ls, orth, lamS**0.5)
    Bs = solve_scad(Xs, Ys, fuse_con, ridge_con, lamR, lamS**0.5, settings['s_it'], settings = settings)
    return Bs

#same as solve_ortho_direct_scad, but also plots lams at each iteration
def solve_ortho_direct_scad_plot(out, organisms, gene_ls, tf_ls, Xs, Ys, orth, priors, lamP, lamR, lamS, settings = None):
    if settings == None:
        settings = get_settings()  
    (lamR, lamP) = pad_lamRP(lamR, lamP, len(organisms))    
    ridge_con = priors_to_constraints(organisms, gene_ls, tf_ls, priors, lamP)    
    fuse_con = orth_to_constraints(organisms, gene_ls, tf_ls, orth, lamS)
    Bs = solve_scad_plot(out, Xs, Ys, fuse_con, ridge_con, lamR, lamS, settings['s_it'], settings = settings)
    return Bs

#parameters as solve_ortho_direct. 
#m_it defines the number of mcp-like iterations to do
def solve_ortho_direct_mcp(organisms, gene_ls, tf_ls, Xs, Ys, orth, priors, lamP, lamR, lamS, settings = None):
    if settings == None:
        settings = get_settings()  
    (lamR, lamP) = pad_lamRP(lamR, lamP, len(organisms))    
    ridge_con = priors_to_constraints(organisms, gene_ls, tf_ls, priors, lamP)    
    fuse_con = orth_to_constraints(organisms, gene_ls, tf_ls, orth, lamS)
    Bs = solve_mcp(Xs, Ys, fuse_con, ridge_con, lamR, lamS, settings['s_it'], settings = settings)
    return Bs

#this is like direct solve, but it breaks up unrelated columns
#solves W = argmin_W ((XW - Y)**2).sum() + constraint related terms
#Xs, Ys: X and Y for each subproblem
#fuse_constraints: fusion constraints
#ridge_constraints: ridge regression constraints. constraints not mentioned are assumed to exist with lam=lambdaR
#it: number of iterations to run
def direct_solve_factor(Xs, Ys, fuse_constraints, ridge_constraints, lambdaRs, adjust=False):
    
    (cols_l, cons_l, ridg_l) = factor_constraints_columns(Xs, Ys, fuse_constraints, ridge_constraints)
    
    Bs = []
    #Initialize matrices to hold solutions
    for i in range(len(Xs)):
        Bs.append(np.zeros((Xs[i].shape[1], Ys[i].shape[1])))
        
    #iterate over subproblems
    for f, cols in enumerate(cols_l):
        if adjust:
            (fcons, rcons, lamR) = adjust_var(Xs, cols_l[f], cons_l[f], ridg_l[f] ,lambdaR)
        else:
            fcons = cons_l[f]
            rcons = ridg_l[f]
            lamR = lambdaRs                  
            
        (F, y, b_inds) = get_Design(Xs, Ys, cols_l[f], fcons, rcons, lamR)

        bsp = scipy.sparse.linalg.lsqr(F, y)#the first result is b, and it needs to be turned into a column vector
        b = bsp[0][:, None] #god this is annoying
        
        #we've now solved a b vector that contains potentially several columns of B. figure out what indices go where, and put them into the right B
        for co_i, co in enumerate(cols):
            (start_ind, end_ind) = b_inds[co_i]
            Bs[co.sub][:, [co.c]] = b[start_ind:end_ind]
    return Bs

#lambdaRs is a list of lambdaRs, one for each species
def get_Design(Xs, Ys, columns, cons_f, cons_r, lambdaRs):
    col_order = dict()
    X_l = []
    Y_l = []
    R_l = []
    counter = 0
    for co in columns:
        sub = co.sub
        col_order[co] = counter
        counter += 1
    for co in columns:
        X = Xs[co.sub]
        X_l.append(X)
        for i in range(X.shape[1]):
            R_l.append(lambdaRs[co.sub])
    for co in columns:
        Y = Ys[co.sub]
        Y_l.append(Y[:, [co.c]])
    
   #compute a cumulative sum over the number of columns in each sub-block, for use as offsets when computing coefficient indices for penalties
    #what this means is that column c in subproblem s has start index cums[s] + r
    cums = [0]+list(np.cumsum(map(lambda x: x.shape[1], X_l)))
    
    #I_vals = np.ones(cums[-1]) * lambdaR
    I_vals = np.ones(cums[-1])*R_l
        #now we go through the ridge constraints and set entries of I
    for con in cons_r:
        ind1 = col_order[(con.c1.sub, con.c1.c)]
        pc1 = cums[ind1] + con.c1.r
        I_vals[pc1] = con.lam
            
    I = scipy.sparse.csr_matrix((I_vals, np.vstack((np.arange(cums[-1]), np.arange(cums[-1])))), shape=(cums[-1],cums[-1]))
    
    #initialize P. we know it has as many rows as constraints
    if len(cons_f):
        P = scipy.sparse.dok_matrix((len(cons_f), cums[-1]))#np.zeros((len(constraints), cums[-1]))
    else:
        P = None
    
    for P_r, con in enumerate(cons_f):
        #get the indices of the diagonal blocks in X corresponding to the coefficients in this constraint
        ind1 = col_order[(con.c1.sub, con.c1.c)]
        ind2 = col_order[(con.c2.sub, con.c2.c)]
        #the locations corresponding to the coefficients we want to penalize are the start index of that block plus the row
        pc1 = cums[ind1] + con.c1.r
        pc2 = cums[ind2] + con.c2.r
            
        P[P_r, pc1] = con.lam
        P[P_r, pc2] = -con.lam
        
    X = diag_concat_sparse(X_l)        
    
    if len(cons_f) and len(P.keys()):
        P = P.asformat('csr')
        F = scipy.sparse.vstack((X, P, I))#F is the penalized design matrix
        Y_l.append(np.zeros((cums[-1] + len(cons_f), 1)))
    else:
        Y_l.append(np.zeros((cums[-1], 1)))
        F = scipy.sparse.vstack((X, I), format='csr')#F is the penalized design matrix    
    y = np.vstack(Y_l)

    #now we get the indices for b. specifically, for each column in columns, we return a start and end index [ ) for b
    b_inds = []
    for co_i in range(len(columns)):
        co = columns[co_i]
        start_ind = cums[co_i]
        end_ind = cums[co_i+1]
        b_inds.append((start_ind, end_ind))
    
    return (F, y, b_inds)

#computes unbiased weighted sample covariance.
#assumes that samples are row? vectors
def weighted_covar(weights, samples):
    u = samples.mean(axis = 1)
    acc = 0
    w_acc = 0
    for i in range(len(weights)):
        xi_u = samples[i, :] - u
        acc += np.dot(xi_u.T, xi_u)
        w_acc += weights[i]
    return (1/(w_acc-1.0)) * acc

#returns differences between betas which have fusion constraints 
def beta_diff(Bs, fuse_cons):
    beta_diffs = np.zeros((len(fuse_cons)))
    for i, con in enumerate(fuse_cons):
        b1 = Bs[con.c1.sub][con.c1.r, con.c1.c]
        b2 = Bs[con.c2.sub][con.c2.r, con.c2.c]
        beta_diffs[i] = np.abs(b1 - b2)
    return beta_diffs

def beta_diff_rand(Bs, k):
    rs0 = np.random.randint(0, Bs[0].shape[0], k)
    rs1 = np.random.randint(0, Bs[1].shape[0], k)
    cs0 = np.random.randint(0, Bs[0].shape[1], k)
    cs1 = np.random.randint(0, Bs[1].shape[1], k)
    return Bs[0][rs0, cs0] - Bs[1][rs1, cs1]
    
        
def gaussian_density(u, s2, x):
    nrm = 1.0 / (2*np.pi*s2)**0.5
    return nrm * np.exp(-1 * (x-u)**2 / (2*s2))


def solve_mcp(Xs, Ys, fuse_con, ridge_con, lamR, lamS, m_it, settings):  
    
    a = settings['a']

    Bs = direct_solve_factor(Xs, Ys, fuse_con, ridge_con, lamR)

    for i in range(m_it-1):
        fuse_con = mcp(Bs, fuse_con, lamS, a=a)
        Bs = direct_solve_factor(Xs, Ys, fuse_con, ridge_con, lamR)
    if settings['return_cons']:
        settings['cons'] = fuse_con
    return Bs

#returns a new set of fusion constraints corresponding to a saturating penalty
def mcp(Bs_init, fuse_constraints, lamS, lamW=None, a=2):
    new_fuse_constraints = []
    
    for i in range(len(fuse_constraints)):
        con = fuse_constraints[i]
        b_init_1 = Bs_init[con.c1.sub][con.c1.r, con.c1.c]
        b_init_2 = Bs_init[con.c2.sub][con.c2.r, con.c2.c]
        theta_init = np.abs(b_init_1 - b_init_2)

        if np.abs(theta_init) <= lamS*a:
            nlamS = lamS*theta_init - (theta_init**2)/(2*a)
            #nlamS = np.abs(lamS-1/a)
            #nlamS = (lamS*theta_init - theta_init/a) / theta_init
            #nlamS = lamS*(theta_init**2) - (theta_init**2)/(2*a) / theta_init

        else:
            nlamS = (a*lamS**2)/2
            #nlamS = ((a*lamS**2)/(2*theta_init**2)) #/ theta_init

        new_con = constraint(con.c1, con.c2, nlamS)
        new_fuse_constraints.append(new_con)
    return new_fuse_constraints



#iteratively adjusts fusion constraint weight to approximate saturating penalty
def solve_scad(Xs, Ys, fuse_con, ridge_con, lamR, lamS, s_it, settings): 
    Bs = direct_solve_factor(Xs, Ys, fuse_con, ridge_con, lamR)
    import copy
    Bs0 = copy.copy(Bs)
    if settings['a'] > 0:
        a = settings['a']
    else:
        a = pick_a(Bs, fuse_con, settings['per'])
        settings['a'] = a     
    for i in range(s_it-1):
        fuse_con = scad(Bs, fuse_con, lamS, a=a)
        Bs = direct_solve_factor(Xs, Ys, fuse_con, ridge_con, lamR)
    deltabetas=[]
    fuse_pen=[]
    for i in range(len(fuse_con)):
        con = fuse_con[i]
        db = Bs0[con.c1.sub][con.c1.r,con.c1.c]-Bs0[con.c2.sub][con.c2.r,con.c2.c]
        fp = con.lam*db**2
        deltabetas.append(db)
        fuse_pen.append(fp)
    plt.scatter(deltabetas,fuse_pen)
    plt.title('solvescad bs0')
    plt.show()

    deltabetas=[]
    fuse_pen=[]
    for i in range(len(fuse_con)):
        con = fuse_con[i]
        db = Bs[con.c1.sub][con.c1.r,con.c1.c]-Bs[con.c2.sub][con.c2.r,con.c2.c]
        fp = con.lam*db**2
        deltabetas.append(db)
        fuse_pen.append(fp)
    plt.scatter(deltabetas,fuse_pen)
    plt.title('solvescad bs_s')
    plt.show()
    if settings['return_cons']:
        settings['cons'] = fuse_con
    return Bs

#same as solve_scad but includes plot_fuse_lams
def solve_scad_plot(out, Xs, Ys, fuse_con, ridge_con, lamR, lamS, s_it, settings): 
    import experiments as e
    import data_sources as ds
    import matplotlib.pyplot as plt
    Bs = direct_solve_factor(Xs, Ys, fuse_con, ridge_con, lamR)
    if settings['a'] > 0:
        a = settings['a']
    else:
        a = pick_a(Bs, fuse_con, settings['per'])
    (constraints, marks, orths) = ds.load_constraints(out)
    marks = np.array(marks)
    true_cons = []
    false_cons = []
    
    for i in range(s_it-1):
        fuse_con = scad(Bs, fuse_con, lamS, a=a)
        Bs = direct_solve_factor(Xs, Ys, fuse_con, ridge_con, lamR)
        lams = np.array(map(lambda x: x.lam, fuse_con))
        true_cons.append(lams[marks==True])
        false_cons.append(lams[marks==False])

    #print false_cons 
    true_arr = np.array(true_cons).T
    false_arr = np.array(false_cons).T
    T = range(true_arr.shape[1])
    for i in range(true_arr.shape[0]):
        plt.plot(T, true_arr[i,:])
    plt.title('true orth')
    plt.show()
    F = range(false_arr.shape[1])
    for i in range(false_arr.shape[0]):
        plt.plot(F, false_arr[i,:])
    plt.title('false orth')
    plt.show()

    if settings['return_cons']:
        settings['cons'] = fuse_con
    return Bs

def scad_prime(theta, lamS, a):
    if theta < lamS:
        return theta * lamS
    else:
        return max(0, (a*lamS - theta) / (a-1))
    
def scad2_prime(theta, lamS, a):

    if theta < a:
        return theta * lamS
    if theta >= a:
        return lamS * max(0, (2*a - theta))

#returns a new set of fusion constraints corresponding to a saturating penalty
def scad(Bs_init, fuse_constraints, lamS, a):
    import matplotlib.pyplot as plt
    count=0
    new_fuse_constraints = []
    import math
    deltabeta=[]
    newlams=[]
    for i in range(len(fuse_constraints)):
        con = fuse_constraints[i]
        b_init_1 = Bs_init[con.c1.sub][con.c1.r, con.c1.c]
        b_init_2 = Bs_init[con.c2.sub][con.c2.r, con.c2.c]
        theta_init = np.abs(b_init_1 - b_init_2)
        #if theta_init is too small, don't want to cause numerical problems
        if theta_init <= 0.001 * (np.abs(b_init_1) + np.abs(b_init_2)):
            nlamS = lamS
            deltabeta.append(theta_init)
            newlams.append(nlamS)
        else:
            nlamS = scad2_prime(theta_init, lamS, a) / theta_init
            deltabeta.append(theta_init)
            newlams.append(nlamS)
        new_con = constraint(con.c1, con.c2, nlamS)
        new_fuse_constraints.append(new_con)
    return new_fuse_constraints

def pick_a(Bs_init, fuse_constraints, perc):
    newpercentile = 100 - 0.5*(100 - perc)
    deltabetas = np.abs(beta_diff(Bs_init, fuse_constraints))
    a = percentile(deltabetas, newpercentile)
    count1 = 0
    count2 = 0
    for i in range(len(deltabetas)):
        if deltabetas[i] >= a:
            count1+=1
        else:
            count2+=1
    return a
    #return a*2

def percentile(scores, perc):
    import matplotlib.pyplot as plt
    #plt.hist(scores)
    #plt.show()
    n = int(round((float(perc)/100) * len(scores), 1))
    scores_sorted = np.sort(scores)
    return scores_sorted[n]

#this code cuts up columns by depth first search
#returns a list of lists of columns associated with each subproblem
#and a list of lists of constraints associated with each subproblem
#and a list of ridge constraints associated with each subproblem
def factor_constraints_columns(Xs, Ys, constraints, ridge):
    columns = []
    
    for sub in range(len(Xs)):
        for c in range(Ys[sub].shape[1]):
            columns.append(column(sub=sub, c=c))
            
    #first we construct an adjacency map on columns
    col_adj = collections.defaultdict(lambda: [])
    
    for con in constraints:
        col_adj[column(sub = con.c1.sub, c=con.c1.c)].append(column(sub = con.c2.sub, c=con.c2.c))
        col_adj[column(sub = con.c2.sub, c=con.c2.c)].append(column(sub = con.c1.sub, c=con.c1.c))

    not_visited = set(columns)
    
    def enumerate_linked_columns(start_col):
        linked_cols = []
        
        def elc_inner(col):
            linked_cols.append(col)
            for adj in col_adj[col]:        
                if adj in not_visited:
                    not_visited.remove(adj)
                    elc_inner(adj)
        
        elc_inner(start_col)
        return linked_cols
    
    col_subls = []
    while len(not_visited):
        #print('\r factoring: %d remaining' % len(not_visited))
        col_fr = not_visited.pop()
        linked_cols = enumerate_linked_columns(col_fr)
        col_subls.append(linked_cols)
        
    #now we look at the linked columns and break up constraints by subproblems, and coefficients
    cols_l = col_subls
    cons_l = []
    ridg_l = []
    #we know exactly how many subproblems there are, so initialize coeffs_l and cons_l
    for i in range(len(col_subls)):
        cons_l.append([])
        ridg_l.append([])
    #build a map from column to subproblem number
    col_to_sub = dict()
    for sub, cols in enumerate(col_subls):
        for col in cols:
            col_to_sub[col] = sub
    
    #now add columns and constraints to the right sub lists
    for con in constraints:
        col = column(sub = con.c1.sub, c = con.c1.c)
        sub = col_to_sub[col]
        cons_l[sub].append(con)
        
    for con in ridge:
        col = column(sub = con.c1.sub, c = con.c1.c)
        sub = col_to_sub[col]
        ridg_l[sub].append(con)
    
    return (cols_l, cons_l, ridg_l)


#this code cuts up columns by depth first search
#returns a list of lists of coefficients associated with each subproblem
#and a list of lists of constraints associated with each subproblem
#NOTE: this is pretty horribly written and should be redone at some point
def factor_constraints_columns_o(Xs, Ys, constraints):
    columns = []
    for sub in range(len(Xs)):
        for c in range(Ys[sub].shape[1]):
            columns.append((sub, c))
    
    constraints_c = set()
    for con in constraints:
        constraints_c.add(constraint(coefficient(con.c1.sub, None, con.c1.c), coefficient(con.c2.sub, None, con.c2.c), None))
        constraints_c.add(constraint( coefficient(con.c2.sub, None, con.c2.c), coefficient(con.c1.sub, None, con.c1.c), None))

    
    not_visited = set(columns)
    
    def factor_helper(col_fr, col_l):
        for col_to in columns:
            if not col_to in not_visited:
                continue
            potential_con = constraint(coefficient(col_fr[0], None, col_fr[1]), coefficient(col_to[0], None, col_to[1]), None)
            
            
            if potential_con in constraints_c:
                col_l.append(col_to)
                not_visited.remove(col_to)
                factor_helper(col_to, col_l)
    col_subls = []
    while len(not_visited):
        #print('factoring: \r %d remaining' % len(not_visited))
        col_fr = not_visited.pop()
        
        col_l = [col_fr]
        factor_helper(col_fr, col_l)
        col_subls.append(col_l)
        
        
    coeffs_l = []
    cons_l = []

    coeff_to_constraints = dict()
    for con in constraints:
        if con.c1 in coeff_to_constraints:
            coeff_to_constraints[con.c1].append(con)
        else:
            coeff_to_constraints[con.c1] = [con]
        if con.c2 in coeff_to_constraints:
            coeff_to_constraints[con.c2].append(con)
        else:
            coeff_to_constraints[con.c2] = [con]
            


    for col_subl in col_subls:
        #now we need to get all the coefficients that really go here
        coeffs = []
        cons = []
        for col in col_subl:
            (sub, c) = col
            for r in range(Xs[sub].shape[1]):
                coeff = coefficient(sub, r, c)
                coeffs.append(coeff)
                if coeff in coeff_to_constraints:
                    for con in coeff_to_constraints[coeff]:    
                        cons.append(con)
        coeffs_l.append(coeffs)
        cons_l.append(cons)
    #print '\ndone: %d subproblems, %d columns'%(len(coeffs_l), len(columns))
    return (coeffs_l, cons_l)
