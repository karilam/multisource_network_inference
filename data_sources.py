import numpy as np
import fused_reg as fr
import random
import os
import collections
from itertools import combinations
from itertools import izip
import csv
import pandas as pd

#SECTION: -------------UTILITY FUNCTIONS------------------------------
#quantile normalizes conditions AND scales to mean zero/unit variance
#AND DOES NOT divides expression matrices by the square root of the sample size
#CHANGED! mean subtract after quantile normalizing
def normalize(exp_a, mean_zero = False):
    
    canonical_dist = np.sort(exp_a, axis=1).mean(axis=0)
    #if mean_zero:
    #    canonical_dist = canonical_dist - canonical_dist.mean()
    canonical_mean = canonical_dist.mean()
    canonical_dist = (canonical_dist - canonical_mean) / canonical_dist.std() + canonical_mean
    
    exp_n_a = np.zeros(exp_a.shape)
    for r in range(exp_a.shape[0]):
        order = np.argsort(exp_a[r, :])
        exp_n_a[r, order] = canonical_dist
    #exp_n_a / np.sqrt(exp_n_a.shape[0])
    if mean_zero:
        exp_n_a = exp_n_a - exp_n_a.mean(axis=0)

    return exp_n_a

#quantile normalizes multiple expression matrices; takes 
def quant_normalize_mats(exp_mats):
    canonical_dist = np.sort(exp_mats[0], axis=1).mean(axis=0)
    canonical_mean = canonical_dist.mean()
    canonical_dist = (canonical_dist - canonical_mean)/canonical_dist.std() + canonical_mean

    new_exp_mats = []
    for i in range(len(exp_mats)):
        exp_norm = np.zeros(exp_mat[i].shape)
        new_colnum = exp_mat[i].shape[1]
        ints = 100/float(new_colnum)
        percentiles = np.arange(0,100,ints)
        new_canonical_dist = map(lambda per: scipy.stats.scoreatpercentile(canonical_dist, per), percentiles)

        for r in range(exp_mat[i].shape[0]):
            order = np.argsort(exp_mat[i][r,:])
            exp_norm[r, order] = new_canonical_dist
        new_exp_mats.append(exp_norm)

    return new_exp_mats
    


def normalize_zscore(exp_a):
    gene_means = exp_a.mean(axis=0)
    gene_stds = exp_a.std(axis=0)

    exp_n_a = np.zeros(exp_a.shape)
    for c in range(exp_a.shape[1]):
        exp_n_a[:,c] = (exp_a[:,c] - gene_means[c]) / gene_stds[c]

    return exp_n_a



#creates a new expression matrix by joining two exp mats, with arbitrary gene coverage
#keeps only the genes that are shared in common
def join_expr_data(names1, names2, exp_a1, exp_a2):
    names = list(set(names1).intersection(set(names2)))
    name_to_ind = {names[x] : x for x in range(len(names))}
    
    def n_to_i(n):
        if n in name_to_ind:
            return name_to_ind[n]
        return -1
    
    exp_a = np.zeros((exp_a1.shape[0] + exp_a2.shape[0], len(names)))
    
    i=0
    #name_ind_to_name1_ind[i] maps an index in names to an index in names1
    name_ind_to_name1_ind = {n_to_i(names1[x]) : x for x in range(len(names1))}

    #name1_fr_inds[i] is the index in names1 that corresponds ti names[i]
    name1_fr_inds = map(lambda x: name_ind_to_name1_ind[x], range(len(names)))
    #array of the previous
    name1_fr_inds_a = np.array(name1_fr_inds)
    #copy into the array
    for r1 in range(exp_a1.shape[0]):
        exp_a[r1, :] = exp_a1[r1, name1_fr_inds_a]
    name_ind_to_name2_ind = {n_to_i(names2[x]) : x for x in range(len(names2))}
    name2_fr_inds = map(lambda x: name_ind_to_name2_ind[x], range(len(names)))
    name2_fr_inds_a = np.array(name2_fr_inds)
    for r2 in range(exp_a2.shape[0]):
        exp_a[r2+exp_a1.shape[0], :] = exp_a2[r2, name2_fr_inds_a]

    return (exp_a, names)

def TFA(X, priors, genes, tfs):
    gene_to_inds = {genes[x] : x for x in range(len(genes))}
    tf_to_inds = {tfs[x] : x for x in range(len(tfs))}

    Xt = X.T
    P = np.zeros((Xt.shape[0],len(tfs)))

    for prior in priors:
        (gene1, gene2) = prior
        if gene1.name in tf_to_inds:
            if gene2.name in gene_to_inds:
                tfi = tf_to_inds[gene1.name]
                gi = gene_to_inds[gene2.name]
                P[gi, tfi] = 1
            elif gene2.name not in gene_to_inds:
                print 'prior gene not in genes'
                print gene2.name
        elif gene1.name not in tf_to_inds:
            print 'prior tf not in tfs'
            print gene1.name

    for tf in tfs:
        tfi = tf_to_inds[tf]
        P[tfi, tfi] = 1

    #for i in range(P.shape[1]):

    #    P[i,i] = 1

    Pi = np.linalg.pinv(P)
    TFA = np.dot(Pi, Xt)

    #TFA = Pi * Xt
    TFA_norm = normalize_zscore(TFA.T)
    return TFA_norm


#SECTION: -----------------------------DATA LOADERS------------------

class data_source():
    def __init__(self):
        print 'NOT IMPLEMENTED!'
        self.N = 0
        self.name = 'cookie monster junk DNA'

    #returns a list of k (approximately) equally sized random folds.
    #folds are lists of integers referring to condition numbers 
    def partition_data(self,k, seed=None):
        state = random.getstate() #don't want to set the seed for everything
        if not seed == None:            
            random.seed(seed)
        conds = np.arange(self.N)
        random.shuffle(conds)
        random.setstate(state)

        incr = float(self.N)/k
        upper_bounds = np.arange(k+1) * self.N / float(k)
        partitions = []
        for i in range(k):
            low = int(upper_bounds[i])
            high = int(upper_bounds[i+1])
            partitions.append(conds[low:high])
        return partitions
    
    #loads data associated with a list of conditions
    def load_data(self,conditions=None):
        if conditions is None:
            conditions = np.arange(self.N)
        conditions = np.array(conditions)
        exp_mat = self.exp_mat[conditions, :]
        tf_mat = self.tf_mat[conditions, :]
        return (exp_mat, tf_mat, self.genes, self.tfs)
    
    #returns a list of priors associated with this data source
    def get_priors(self):
        print 'NOT IMPLEMENTED!'
        return None

#rearranges genes (and corresponding expression matrices) so that tfs appear first, in the order that they appear within genes
def arrange_tfs_first(tfs, genes, tf_mat, exp_mat):
    tfs_set = set(tfs)
            #compute the (gene) indices of tfs and non-tfs
    tf_indices = filter(lambda x: genes[x] in tfs_set, range(len(genes)))
    ntf_indices = filter(lambda x: not genes[x] in tfs_set, range(len(genes)))

    new_order = tf_indices + ntf_indices
            
    tfs = map(lambda x: genes[x], tf_indices)            
    genes = map(lambda x: genes[x], new_order)
    
    tf_mat = exp_mat[:, tf_indices]            
    exp_mat = exp_mat[:, new_order]
    return (tfs, genes, tf_mat, exp_mat)

#turns expression matrices along with time-difference list of lists of strings that specify
#cond1, cond2, gap into time difference approximations of the derivative according to
# [X(t_2) - X(t_1)] / (t_2 - t_1) + x(t_1)/ tc
def generate_timeseries(tf_mat, exp_mat, td, tc):
    N = len(td)
    M_tf = tf_mat.shape[1]
    M_exp = exp_mat.shape[1]
    td_tf_mat = np.zeros((N, M_tf))
    td_exp_mat = np.zeros((N, M_exp))
    
    for td_row in range(N):
        (fr_s, to_s, s_s) = td[td_row]
        fr = int(fr_s)
        to = map(int, to_s.split(','))
        s = float(s_s)

        #print 'current: %d, preceding: %d' % (to, fr)
        T1 = tf_mat[fr, :]
        X1 = exp_mat[fr, :]
    
        X2 = exp_mat[to, :].mean(axis=0)
    
        Xdt = (X2 - X1)/s
        Xtc = X1 / tc
        #print Xdt
        #print Xtc
        td_tf_mat[td_row, :] = T1
        td_exp_mat[td_row, :] = Xdt - Xtc
        #print 'using for row %d' % td_row
        #print Xdt + Xtc
    return (td_tf_mat, td_exp_mat)

#this class assumes directory layout we use for generated data
#if use_TFA != False, 0 < use_TFA <=1, where 1 means TFA uses all the priors. then get_priors()[0] returns the priors that TFA is NOT using. get_priors()[1] returns the priors TFA is using.
class standard_source(data_source):
    def __init__(self, datadir, org_ind, use_TFA=False):
        #strip last separator
        if datadir[-1] == os.sep:
            datadir = datadir[0:-1]
        
        expr_fn = os.path.join(datadir, 'expression%d' % (org_ind+1))
        prior_fn = os.path.join(datadir,'priors%d') % (org_ind+1)
        gold_fn = os.path.join(datadir, 'gold%d') % (org_ind+1)
        tfs_fn = os.path.join(datadir, 'tfnames%d') % (org_ind+1)
        time_diff_fn = os.path.join(datadir,'time_diffs%d' % (org_ind+1))
        name_key = 'organism%d' % (org_ind+1)
        tc_key = 'tc%d' % (org_ind+1)
        
        name = 'dunno'
        
        def spl_ln_t(s):
            return map(lambda x: x.split('\t'), filter(len, s.split('\n')))
        #scan description file for keys
        with file(os.path.join(datadir, 'description')) as descr_f:
            fslt = spl_ln_t(descr_f.read())
            for entry in fslt:
                if len(entry) == 2:
                    if entry[0] == name_key:
                        name = entry[1]
                    if entry[0] == tc_key:
                        tc = float(entry[1])
        #now get the expression data
        with file(expr_fn) as expr_f:
            expr_fsnt = spl_ln_t(expr_f.read())
            genes = expr_fsnt[0]
            nconds = len(expr_fsnt)-1
            exp_mat = np.zeros((nconds, len(genes)))
            for row in range(nconds):
                expr = expr_fsnt[row + 1]
                expr_arr = np.array(expr)
                exp_mat[row,:] = expr_arr
        #now get the tf names
        with file(tfs_fn) as tfs_f:
            tfs_fs = tfs_f.read()
            tfs = filter(len, tfs_fs.split('\n'))
        
        #read the time difference file
        with file(time_diff_fn) as time_diff_f:
            tdlt = spl_ln_t(time_diff_f.read())
            #td = np.array(tdlt).astype(float)
            td = tdlt
        
        #now we extract the transcription factor expression matrix
        gene_to_ind = {genes[i] : i for i in range(len(genes))}
        tf_inds = map(lambda x: gene_to_ind[x], tfs)
        tf_mat = exp_mat[:, np.array(tf_inds)]
        
        #rearrange order so that tfs appear first
        (tfs, genes, tf_mat, exp_mat) = arrange_tfs_first(tfs, genes, tf_mat, exp_mat)

        self.prior_fn = prior_fn
        self.gold_fn = gold_fn
        self.name = name
        self.use_TFA = use_TFA

        if use_TFA != False:
            tf_mat = TFA(exp_mat, self.get_priors()[1], genes, tfs)        
            (tf_mat, exp_mat) = generate_timeseries(tf_mat, exp_mat, td, tc)

        else:
            (tf_mat, exp_mat) = generate_timeseries(tf_mat, exp_mat, td, tc)
        #DO NOT normalize here
        #exp_mat = normalize(exp_mat, True)
        #tf_mat = normalize(tf_mat, False)
    
        self.exp_mat = exp_mat
        self.tf_mat = tf_mat
        self.datadir = datadir

        self.genes = genes
        self.tfs = tfs
        self.N = exp_mat.shape[0]         
        
    def get_priors(self):
        def r_partition(x, t):
            inds = np.arange(len(x))
            random.shuffle(inds)
            p1 = map(lambda i: x[i], inds[0:t])
            p2 = map(lambda i: x[i], inds[t:])
            return (p1, p2)

        p = file(self.prior_fn)
        ps = p.read()
        psn = filter(len, ps.split('\n'))
        psnt = map(lambda x: x.split('\t'), psn)
        priors = map(lambda x: (fr.one_gene(x[0], self.name), fr.one_gene(x[1], self.name)), psnt)
        
        signs = []
        for x in psnt:
            sign = x[2]
            
            if sign == '1':
                signs.append(1)
            elif sign == '-1':
                signs.append(-1)
            else:
                signs.append(0)
        p.close()

        if self.use_TFA != False:
            pct_priors = self.use_TFA
            t = int(pct_priors*len(priors))
            (te_priors, tfa_priors) = r_partition(priors, t)
            return (te_priors, tfa_priors, signs)    

        else:
            return (priors, signs)


    def get_gold(self):
        p = file(self.gold_fn)
        ps = p.read()
        psn = filter(len, ps.split('\n'))
        psnt = map(lambda x: x.split('\t'), psn)
        priors = map(lambda x: (fr.one_gene(x[0], self.name), fr.one_gene(x[1], self.name)), psnt)
        
        signs = []
        for x in psnt:
            sign = x[2]
            
            if sign == '1':
                signs.append(1)
            elif sign == '-1':
                signs.append(-1)
            else:
                signs.append(0)
        p.close()
        


        return (priors, signs)

#SUBSECTION: ----------------------REAL DATA----------------------



class ba_timeseries(data_source):
    def __init__(self):
        f = file('data/bacteria1/Normalized_data_RMA._txt')
        fs = f.read()
        fsn = filter(len, fs.split('\n'))
        fsnt = map(lambda x: x.split('\t'), fsn)
        conds = fsnt[0][1:]
    #first line is SCAN REF
    #second line is composite element REF
    #lines 3-end are data
        exp_mat_t = np.zeros((len(fsnt)-2, len(conds)))#first col gene
        genes = []
        f_tf = file('data/bacteria1/tfNamesAnthracis')
        f_tfs = f_tf.read()
        tfs = filter(len, f_tfs.split('\n'))


        for r in range(exp_mat_t.shape[0]):
            gene_str_full = fsnt[r+2][0]
        #gene name is 4th element, separated by ':'
            gene_str = gene_str_full.split(':')[3]
            gene_str = gene_str.replace('_pXO1_','').replace('_pXO2','')#what is this? dunno!
            expr = np.array(fsnt[r+2][1:])
            exp_mat_t[r, :] = expr
            genes.append(gene_str)

    #require that tfs be genes that we have data for!
        tfs = filter(lambda x: x in genes, tfs)

        tf_mat_t = np.zeros((len(tfs), len(conds)))

        gene_to_ind = {genes[x] : x for x in range(len(genes))}
        for ti in range(len(tfs)):
            gi = gene_to_ind[tfs[ti]]
            tf_mat_t[ti, :] = exp_mat_t[gi, :]
        exp_mat = exp_mat_t.T
    
        tf_mat = tf_mat_t.T

        #self.exp_mat = normalize(exp_mat, True)
        #self.tf_mat = normalize(tf_mat, False)
        self.exp_mat = normalize_zscore(exp_mat)
        self.tf_mat = normalize_zscore(tf_mat)        
        self.genes = genes
        self.tfs = tfs
        self.N = exp_mat.shape[0]
        
        self.name = 'ba_iron'

    #returns a list of priors associated with this data source
    def get_priors(self):
        return ([],[])
        


#B anthracis data relating to iron starvation conditions
class ba_iron(data_source):
    def __init__(self):
        f = file('data/bacteria1/normalizedgeneexpressionvalues.txt')
        fs = f.read()
        fsn = filter(len, fs.split('\n'))
        fsnt = map(lambda x: x.split('\t'), fsn)
        conds = fsnt[0][1:]
    #first line is SCAN REF
    #second line is composite element REF
    #lines 3-end are data
        f_tf = file('data/bacteria1/tfNamesAnthracis')
        f_tfs = f_tf.read()
        tfs = filter(len, f_tfs.split('\n'))
        
        exp_mat_t = np.zeros((len(fsnt)-2, len(conds)))#first col gene
        genes = []

        for r in range(exp_mat_t.shape[0]):
            gene_str = fsnt[r+2][0]
        
        #gene name is 4th element, separated by ':'
        
            expr = np.array(fsnt[r+2][1:])
            exp_mat_t[r, :] = expr
            genes.append(gene_str)

    #require that tfs be genes that we have data for!
        tfs = filter(lambda x: x in genes, tfs)

        tf_mat_t = np.zeros((len(tfs), len(conds)))

        gene_to_ind = {genes[x] : x for x in range(len(genes))}
        for ti in range(len(tfs)):
            gi = gene_to_ind[tfs[ti]]
            tf_mat_t[ti, :] = exp_mat_t[gi, :]
        exp_mat = exp_mat_t.T

        tf_mat = tf_mat_t.T
        #NOTE removing normalization here
        #NOTE: adding back in, removed from standard source
        #self.exp_mat = normalize(exp_mat, True)
        #self.tf_mat = normalize(tf_mat, False)
        self.exp_mat = normalize_zscore(exp_mat)
        self.tf_mat = normalize_zscore(tf_mat)
        #self.exp_mat = exp_mat
        #self.tf_mat = tf_mat
        self.genes = genes
        self.tfs = tfs
        self.N = exp_mat.shape[0]
        self.name = 'ba_time'
        
    #returns a list of priors associated with this data source
    def get_priors(self):
        return ([],[])

class subt(data_source):
    def __init__(self):
    
        expr_fn = 'data/bacteria1/B_subtilis.csv'
        tfs_fn = 'data/bacteria1/tfNames_subtilis.txt'
        f = file(expr_fn)
        fs = f.read()
        fsl = filter(len, fs.split('\n'))
        fslc = map(lambda x: x.split(','), fsl)
        f.close()
    
        t = file(tfs_fn)
        ts = t.read()
        tfs = filter(len, ts.split('\n'))
        t.close()

        tfs_set = set(tfs)

    
        conds = fslc[0]
        genes = map(lambda x: x[0], fslc[1:])
        exp_mat_t = np.zeros((len(genes), len(conds)))
        for r in range(len(genes)):
            conds_f = map(float, fslc[1+r][1:])
            conds_a = np.array(conds_f)
            exp_mat_t[r, :] = conds_a

        tf_mat_t = np.zeros((len(tfs), len(conds)))
        gene_to_ind = {genes[x] : x for x in range(len(genes))}
    
        for ti in range(len(tfs)):
            gi = gene_to_ind[tfs[ti]]
            tf_mat_t[ti, :] = exp_mat_t[gi, :]

        exp_mat = exp_mat_t.T
        tf_mat = tf_mat_t.T
        #NOTE: removed normalization
        #NOTE: adding back in, removed from standard source
        #self.exp_mat = normalize(exp_mat, True)
        #self.tf_mat = normalize(tf_mat, False)
        self.exp_mat = normalize_zscore(exp_mat)
        self.tf_mat = normalize_zscore(tf_mat)
        #self.exp_mat = exp_mat
        #self.tf_mat = tf_mat
        self.genes = genes
        self.tfs = tfs
        self.N = exp_mat.shape[0]
        self.offset = 0
        self.name = 'B_subtilis'

    def get_priors(self):
        #priors_fn = 'data/bacteria1/bsubt_priors_011916'
        priors_fn = 'data/bacteria1/gsSDnamesWithActivitySign082213'
        p = file(priors_fn)
        ps = p.read()
        psn = filter(len, ps.split('\n'))
        psnt = map(lambda x: x.split('\t'), psn)
        priors = map(lambda x: (fr.one_gene(x[0], self.name), fr.one_gene(x[1], self.name)), psnt)
        signs = map(lambda x: [-1,1][x[2]=='activation'], psnt)
        p.close()
        return (priors, signs)

#just copied subt
class subt_eu(data_source):
    def __init__(self):
    
        expr_fn = 'data/bacteria1/B_subtilis_eu.csv'
        tfs_fn = 'data/bacteria1/tfNames_subtilis_eu.txt'
        f = file(expr_fn)
        fs = f.read()
        fsl = filter(len, fs.split('\n'))
        fslc = map(lambda x: x.split(','), fsl)
        f.close()
    
        t = file(tfs_fn)
        ts = t.read()
        tfs = filter(len, ts.split('\n'))
        t.close()

        tfs_set = set(tfs)


        conds = fslc[0]
        genes = map(lambda x: x[0], fslc[1:])
        exp_mat_t = np.zeros((len(genes), len(conds)))
        for r in range(len(genes)):
            conds_f = map(float, fslc[1+r][1:])
            conds_a = np.array(conds_f)
            exp_mat_t[r, :] = conds_a

        tf_mat_t = np.zeros((len(tfs), len(conds)))
        gene_to_ind = {genes[x] : x for x in range(len(genes))}
    
        for ti in range(len(tfs)):
            gi = gene_to_ind[tfs[ti]]
            tf_mat_t[ti, :] = exp_mat_t[gi, :]

        exp_mat = exp_mat_t.T
        tf_mat = tf_mat_t.T
        #NOTE: removed normalization
        #NOTE: adding back in, removed from standard source
        #self.exp_mat = normalize(exp_mat, True)
        #self.tf_mat = normalize(tf_mat, False)
        self.exp_mat = normalize_zscore(exp_mat)
        self.tf_mat = normalize_zscore(tf_mat)
        #self.exp_mat = exp_mat
        #self.tf_mat = tf_mat
        self.genes = genes
        self.tfs = tfs
        self.N = exp_mat.shape[0]
        self.offset = 0
        self.name = 'B_subtilis'

    def get_priors(self):
        #priors_fn = 'data/bacteria1/bsubt_priors_011916'
        priors_fn = 'data/bacteria1/gsSDnamesWithActivitySign082213'
        p = file(priors_fn)
        ps = p.read()
        psn = filter(len, ps.split('\n'))
        psnt = map(lambda x: x.split('\t'), psn)
        priors = []
        signs = []
        for line in psnt:
            if line[0] in self.genes:
                if line[1] in self.genes:
                    priors.append((fr.one_gene(line[0], self.name),fr.one_gene(line[1], self.name)))
                    signs.append([-1,1][line[2]=='activation'])
        p.close()
        return (priors, signs)


#NOTE: I just changed anthracis to quantile normalize separately before combining. I can't see a downside to doing it this way, or the other way
class anthr(data_source):
    def __init__(self):
        self.ba_i = ba_iron()
        self.ba_t = ba_timeseries()
        self.name = 'B_anthracis'
        (e, genes) = join_expr_data(self.ba_i.genes, self.ba_t.genes, self.ba_i.exp_mat, self.ba_t.exp_mat)
        (t, tfs) = join_expr_data(self.ba_i.tfs, self.ba_t.tfs, self.ba_i.tf_mat, self.ba_t.tf_mat)
        genes = map(lambda x: x.replace('_at',''), genes)
        tfs = map(lambda x: x.replace('_at',''), tfs)


        self.exp_mat = e
        self.genes = genes
        self.tf_mat = t
        self.tfs = tfs
        self.N = self.exp_mat.shape[0]

    #returns a list of k (approximately) equally sized random folds.
    #folds are lists of integers referring to condition numbers 
    def partition_data(self, k):
        par1 = self.ba_i.partition_data(k)
        par2 = self.ba_t.partition_data(k)
        partitions = []
        off = self.ba_i.N #add this offset to the timeseries conds
        for i in range(k):
            partitions.append(np.hstack((par1[i], off + par2[i])))
        return partitions
    
    def get_priors(self):
        return ([], [])

#NOTE: I am missing the (largely unsuccessful) loader for subtilis that takes into account timeseries data



#SECTION: ----------------------------------ORTHOLOGY LOADERS----------

#orgs is a list of organisms of interest, if we are interested in just a subset of the organisms in a folder
def load_orth(orth_fn, organisms, orgs=None):
    f = file(orth_fn)
    fs = f.read()
    fsn = filter(len, fs.split('\n'))
    fsnt = map(lambda x: x.split('\t'), fsn)
    f.close()    

    cols = []
    if orgs != None:
        cols = map(lambda x: organisms.index(x), orgs)
    else:
        cols = np.arange(len(organisms))

    orths = []
    for o in fsnt:
        real = o[-1] == 'True'
        
        genes = []
        for org in range(len(o)-1):
            if org in cols:
                gs = filter(len, o[org].split(','))
                for g in gs:
                    genes.append(fr.one_gene(name = g, organism = organisms[org]))

        orths_row = map(lambda g1g2: fr.orthology(genes = g1g2, real = real), combinations(genes, 2))
        for orth in orths_row:
            orths.append(orth)
    return orths

#returns the MARKED constraints associated with a particular data directory, in standard format

def load_constraints(data_fn, orth_f='orth', orgs=None):

    organisms = []
    all_orgs = []
    gene_ls = []
    tf_ls = []
    num_species = 0
    if orgs != None:
        while os.path.isfile(os.path.join(data_fn, 'expression%d' % (num_species+1))):
            dsi = standard_source(data_fn,num_species)
            all_orgs.append(dsi.name)
            if dsi.name in orgs:
                organisms.append(dsi.name)
                gene_ls.append(dsi.genes)
                tf_ls.append(dsi.tfs)
            num_species +=1
        num_species-=1

    else:
        while os.path.isfile(os.path.join(data_fn, 'expression%d' % (num_species+1))):
            dsi = standard_source(data_fn,num_species)
            organisms.append(dsi.name)
            all_orgs.append(dsi.name)
            gene_ls.append(dsi.genes)
            tf_ls.append(dsi.tfs)
            num_species +=1
        num_species-=1

    orth_fn = os.path.join(data_fn, orth_f)
    orth = load_orth(orth_fn, all_orgs, organisms)

    (constraints, marks) = fr.orth_to_constraints_marked(organisms, gene_ls, tf_ls, orth, 1.0)
    return (constraints, marks, orth)

ba_bs_orth = lambda: load_orth('data/bacteria1/bs_ba_ortho_804',['B_anthracis','B_subtilis'])

#Okuda, S. and Yoshizawa, A.C.; ODB: a database for operon organizations, 2011 update. Nucleic Acids Res.39(Database issue):D552-555 (2011). [pubmed]
#Okuda S, Katayama T, Kawashima S, Goto S, and Kanehisa M.; ODB: a database of operons accumulating known operons across multiple genomes. Nucleic Acids Res. 34(Database issue):D358-362 (2006).[pubmed]

#takes downloaded operon file and writes new file with just bsu operons
#makes tfs orthologs of themselves
def bsu_operon_to_orth(op_file, new_orth_file):
    f = file(op_file)
    w = file(new_orth_file, 'w')
    f1 = f.read().split('\n')
    f1s = filter(len, f1)
    f.close()
    bsu = filter((lambda x: x.split('\t')[1] == 'bsu'), f1s)
    for line in range(len(bsu)):
        genes = bsu[line].split('\t')[3].split(',')
        w.write(','.join(genes) + '\t' + '\t' + 'True' + '\n')
    sub = subt()
    (bs_e, bs_t, bs_genes, bs_tfs) = sub.load_data()
    for tf in bs_tfs:
        w.write(tf + ',' + tf + '\t' + '\t' + 'True' + '\n')
    w.close()


#SECTION: ----------------------------DATA GENERATORS------------------
#generates a matrix Y from the linear model specified by B.
#x is sampled uniformly from -1 to 1
#input N: the number of samples
#input B: matrix specifying linear model
#input noise_std: std of noise!
#returns: (X, Y)
def generate_from_linear(N, B, noise_std):
    X_pre = np.random.randn(N, B.shape[0])#
    #X_pre = 1-2*np.random.random((N, B.shape[0]))
    Y_pre = np.hstack((X_pre, np.zeros((N, B.shape[1]-B.shape[0]))))

    Y = Y_pre + np.dot(X_pre,B) + noise_std*np.random.randn(N, B.shape[1])
    X = Y[:, 0:B.shape[0]]
    return (X, Y, X_pre, Y_pre)    

#turns a dictionary mapping genes to lists of orthologs into a list of orthology objects, marked as real
 
def omap_to_orths(omap, real):
    orths = []
    for gene1 in omap.keys():       
        for gene2 in omap[gene1]:
            orths.append( fr.orthology(genes=(gene1, gene2),real=real))
    return orths

#this functions somehow generates orthology mappings, and adds them to omap
#this strange design was chosen over simply returning an orthology mapping because i am a sadist
#this works by generating orthology groups of random size from 2 to max_grp_size (inclusive). Each organism must contribute at least one gene to the orthology group (this could change later). This process continues until the total number of genes in orthology groups is approximately pct_fused percent of the total genes. 
#what goes into the orthology mapping! well.... it seems to be a dictionary mapping organism name/gene tuples to lists of organism name/gene tuples. This could be rewritten to use one_gene objects, because that's what they're for, but really the whole thing should be burned to the ground
def build_orth(genes1, genes2, max_grp_size, pct_fused, omap, organisms, shuffle=True):
    shuffle=False
    if shuffle:
        random.shuffle(genes1)
        random.shuffle(genes2)
    
    amt_fused = np.floor((len(genes1)+len(genes2))*pct_fused) #round down
    ind1 = 0
    ind2 = 0
    while ind1 + ind2 < amt_fused:
        grp_size = random.randrange(2, max_grp_size+1)

        grp1_size = random.randrange(1, grp_size)
        grp2_size = grp_size - grp1_size
        #modify the group sizes to deal with not enough of one sub
        grp1_size = min(grp1_size, len(genes1)-ind1)
        grp2_size = min(grp2_size, len(genes2)-ind2) 
        
        
        #enumerate ALL of the 1-1 orthologies that span different organisms within this group
        #fixed a bug here wherein if the two species shared gene names, things fucked up. this can be fixed by only considering one direction of orthology.
        for i in range(ind1, ind1 + grp1_size):
            for j in range(ind2, ind2 + grp2_size):
                
                sub1c = fr.one_gene(genes1[i], organism = organisms[0])
                sub2c = fr.one_gene(genes2[j], organism = organisms[1])
                omap[sub1c].append(sub2c)
                
                
        ind1 += grp1_size
        ind2 += grp2_size
    

#TO THE BEST OF MY KNOWLEDGE THIS IS WHAT HAPPENS. this function builds two bata matrices, with dims specified by b1/tfg_count2. pct_fused percent of genes and TFs are (separately) assigned to orthology groups of size 2-max_grp_size (sampled uniformly). Coefficientsa are assigned value from standard normal distributions. Coefficients that are fused (that is, both gene and tf are fused) are assigned to the same value, perturbed by fuse_std. Sparse specifies the proportion of fused groups of coefficients (including groups of size 1) that are set to 0
#NOTE: tfg_count1 is the number of tfs and non-tf genes in a tuple
#shuffle - do you shuffle the order of genes??
def fuse_bs_orth(tfg_count1, tfg_count2, max_grp_size, pct_fused, fuse_std, sparse, organisms,shuffle=False):
    #create empty expression matrices
    
    b1 = np.nan * np.ones((tfg_count1[0], tfg_count1[0] + tfg_count1[1]))
    b2 = np.nan * np.ones((tfg_count2[0], tfg_count2[0] + tfg_count2[1]))
    
#b1 = np.nan * np.ones(tfg_count1)
    #b2 = np.nan * np.ones(tfg_count2)
    
    tfs1 = map(lambda x: organisms[0]+'t'+str(x), range(b1.shape[0]))
    tfs2 = map(lambda x: organisms[1]+'t'+str(x), range(b2.shape[0]))   
    #THESE ARE GENES THAT ARE NOT ALSO TFS
    #This is done so that we can generate orthologies between genes and tfs separately
    non_tf_genes1 = map(lambda x: organisms[0]+'g'+str(x), range(tfg_count1[1]))
    non_tf_genes2 = map(lambda x: organisms[1]+'g'+str(x), range(tfg_count2[1]))
    
    
    
    omap = collections.defaultdict(lambda: [])#dict()
    build_orth(tfs1, tfs2, max_grp_size, pct_fused, omap, organisms,shuffle=shuffle)
    build_orth(non_tf_genes1, non_tf_genes2, max_grp_size, pct_fused, omap, organisms,shuffle=shuffle)
    
    #now we append the tfs
    genes1 = tfs1 + non_tf_genes1
    genes2 = tfs2 + non_tf_genes2
    orths = omap_to_orths(omap, True)
    
    bs = [b1, b2]
    genes = [genes1, genes2]
    tfs = [tfs1, tfs2]

    gene_inds = [{genes1[i] : i for i in range(len(genes1))}, {genes2[i] : i for i in range(len(genes2))}]
    tf_inds = [{tfs1[i] : i for i in range(len(tfs1))}, {tfs2[i] : i for i in range(len(tfs2))}]
    orths_query = set(orths)
    
    #recursively fills everything connected to r,c,organism
    #NOTE: we are relying on the fact that transcription factors are also genes
    def fill(r, c, val, std, organism):
        fusion_noise = np.random.randn()*std
        fill_val = val + fusion_noise
        organism_ind = organisms.index(organism) #blegh
        if not np.isnan(bs[organism_ind][r,c]): #already filled, return
            return           
        bs[organism_ind][r,c] = fill_val
        
        tf_fill = fr.one_gene(tfs[organism_ind][r], organism)
        gene_fill = fr.one_gene(genes[organism_ind][c], organism)

        for tf_orth in omap[tf_fill]:
            for g_orth in omap[gene_fill]:
                #we need to verify that both orths are in the same organism as each other
                if tf_orth.organism == g_orth.organism:
                    orth_org_ind = organisms.index(tf_orth.organism)
                    
                    orth_r = tf_inds[orth_org_ind][tf_orth.name]
                    orth_c = gene_inds[orth_org_ind][g_orth.name]
                    fill(orth_r, orth_c, val, std, organisms[orth_org_ind])

    #go through each value in each beta matrix and, if not yet filled, fill it in along with all of its linked values
    for organism_ind, b in enumerate(bs):
        for r in range(b.shape[0]):
            for c in range(b.shape[1]):
                if np.isnan(b[r,c]):
                    if np.random.random() < sparse:
                        val = 0
                        std = 0
                    else:
                        val = np.random.randn()
                        std = fuse_std
                    fill(r, c, val, std, organisms[organism_ind])

    return (b1, b2, orths, genes1, tfs1, genes2, tfs2)


# generates priors from a given beta matrix B, some of which may be wrong or missing
def generate_faulty_priors(B, genes, tfs, falsepos, falseneg):
    priors = []
    fakepriors = []
    for r in range(B.shape[0]):
        for c in range(B.shape[1]):
            if B[r, c] != 0:
                priors.append((tfs[r], genes[c]))
                
            if B[r, c] == 0:
                fakepriors.append((tfs[r], genes[c]))
    num_to_remove = int(falseneg * len(priors))
    num_to_add = int(falsepos*(len(priors) - num_to_remove)/(1-falsepos))    
    
    random.shuffle(priors)
    random.shuffle(fakepriors)
    to_keep = priors[0:(len(priors) - num_to_remove)]
    to_add = fakepriors[0:num_to_add]
    
    return  to_keep + to_add 


#now also marks orths as real or fake
#orth_falsepos: number of fake orthologies to add, in units of the original number of orthologies
#orth_falseneg: fraction of original orthologies to remove
def generate_faulty_orth(orths, genes1, tfs1, genes2, tfs2, organisms, falsepos, falseneg):
    #make a list of sets containing gene fusion groups to prevent from adding false orths that result in unduly large fusion groups
        
    num_to_remove = int(falseneg * len(orths))
    num_to_add = int(falsepos * len(orths))
    
    random.shuffle(orths)
    final_real_ind = max(0, len(orths)-num_to_remove)
    orths_retain = orths[0:final_real_ind]
    orth_genes = set()
    for orth in orths_retain:
        orth_genes.add(orth.genes[0])
        orth_genes.add(orth.genes[1])
    
    all_possible_orths = []

    
    for gene1 in genes1:
        for gene2 in genes2:
            possible_orth = fr.orthology(genes = (fr.one_gene(name=gene1, organism = organisms[0]), fr.one_gene(name=gene2, organism = organisms[1])), real = False)
            all_possible_orths.append(possible_orth)
            
    for tf1 in tfs1:
        for tf2 in tfs2:
            possible_orth = fr.orthology(genes = (fr.one_gene(name=tf1, organism = organisms[0]), fr.one_gene(name=tf2, organism = organisms[1])), real = False)
            all_possible_orths.append(possible_orth)

    random.shuffle(all_possible_orths)
    
    to_add = []
    #add potential orthologies until enough have been added
    #don't add orthologies for orths that we already have
    for candidate_orth in all_possible_orths:
    
        if len(to_add) == num_to_add:
            break
        if candidate_orth.genes[0] in orth_genes or candidate_orth.genes[1] in orth_genes:
            continue
        #add these, to make sure no double dipping
        orth_genes.add(candidate_orth.genes[0])
        orth_genes.add(candidate_orth.genes[1])
        to_add.append(candidate_orth)

    #print '%d real, %d fake' % (len(orths_retain), len(to_add))
    real_fake_orths = orths_retain + to_add
    
    return real_fake_orths

#adds false positives/negatives without generating all possible orthologies.
#may run forever with large values of false positives
#now also marks orths as real or fake
#orth_falsepos: number of fake orthologies to add, in units of the original number of orthologies
#orth_falseneg: fraction of original orthologies to remove
def generate_faulty_orth2(orths, genes1, tfs1, genes2, tfs2, organisms, falsepos, falseneg):
    #make a list of sets containing gene fusion groups to prevent from adding false orths that result in unduly large fusion groups
        
    num_to_remove = int(falseneg * len(orths))
    num_to_add = int(falsepos * len(orths))
    
    random.shuffle(orths)
    final_real_ind = max(0, len(orths)-num_to_remove)
    orths_retain = orths[0:final_real_ind]
    orth_genes = set()
    for orth in orths_retain:
        orth_genes.add(orth.genes[0])
        orth_genes.add(orth.genes[1])
    
    #take tfs out of genes
    tfs1_s = set(tfs1)
    tfs2_s = set(tfs2)
    genes1 = filter(lambda x: not x in tfs1_s, genes1)
    genes2 = filter(lambda x: not x in tfs1_s, genes2)
            
    to_add = []
    #add potential orthologies until enough have been added
    #don't add orthologies for orths that we already have
    while len(to_add) < num_to_add:
        #first decide if a tf or gene orth. 

        if random.random() < float(len(tfs1)*len(tfs2))/(len(genes1)*len(genes2)):
            g1 = random.sample(tfs1,1)[0]
            g2 = random.sample(tfs2,1)[0]
        else:
            g1 = random.sample(genes1,1)[0]
            g2 = random.sample(genes2,1)[0]
        candidate_orth = fr.orthology(genes = (fr.one_gene(name=g1, organism = organisms[0]), fr.one_gene(name=g2, organism = organisms[1])), real = False)
        if candidate_orth.genes[0] in orth_genes or candidate_orth.genes[1] in orth_genes:
            continue
        #add these, to make sure no double dipping
        orth_genes.add(candidate_orth.genes[0])
        orth_genes.add(candidate_orth.genes[1])
        to_add.append(candidate_orth)

    #print '%d real, %d fake' % (len(orths_retain), len(to_add))
    real_fake_orths = orths_retain + to_add
    
    return real_fake_orths

#creates sets of tfs, genes which are not already in orth, and creates faulty orths from them
def generate_faulty_orth3(orths, genes1, tfs1, genes2, tfs2, organisms, falsepos, falseneg):

    num_to_remove = int(falseneg * len(orths))
    num_to_add = int(falsepos * len(orths))
    
    random.shuffle(orths)
    final_real_ind = max(0, len(orths)-num_to_remove)
    orths_retain = orths[0:final_real_ind]
    orth_genes = set()
    for orth in orths_retain:
        orth_genes.add(orth.genes[0])
        orth_genes.add(orth.genes[1])
    
    #take tfs out of genes
    tfs1_s = set(tfs1)
    tfs2_s = set(tfs2)
    genes1 = set(filter(lambda x: not x in tfs1_s, genes1))
    genes2 = set(filter(lambda x: not x in tfs1_s, genes2))
            
    to_add = []
    #add potential orthologies until enough have been added
    #don't add orthologies for orths that we already have
    while len(to_add) < num_to_add:
        #first decide if a tf or gene orth. 
        if random.random() < float(len(tfs1_s)*len(tfs2_s))/(len(genes1)*len(genes2)):
            g1 = random.sample(tfs1_s,1)[0]
            g2 = random.sample(tfs2_s,1)[0]
            tfs1_s.remove(g1)
            tfs2_s.remove(g2)
        else:
            g1 = random.sample(genes1,1)[0]
            g2 = random.sample(genes2,1)[0]
            genes1.remove(g1)
            genes2.remove(g2)
        candidate_orth = fr.orthology(genes = (fr.one_gene(name=g1, organism = organisms[0]), fr.one_gene(name=g2, organism = organisms[1])), real = False)
        #add these, to make sure no double dipping
        orth_genes.add(candidate_orth.genes[0])
        orth_genes.add(candidate_orth.genes[1])
        to_add.append(candidate_orth)

    #print '%d real, %d fake' % (len(orths_retain), len(to_add))
    real_fake_orths = orths_retain + to_add
    
    return real_fake_orths



#writes fake data, assumes some reasonable defaults
def write_fake_data1(out_dir=None, tfg_count1=(5,10), tfg_count2=(5,10), N1=10, N2=10, max_grp_size=2, pct_fused=1.0, fuse_std=0.5, sparse=0.0, organisms = ['uno','dos'], prior_falsepos=0.0, prior_falseneg=0.0, orth_falsepos=0.0, orth_falseneg=0.0, measure_noise1=0.1, measure_noise2=0.1):

    (B1, B2, orths, genes1, tfs1, genes2, tfs2) = fuse_bs_orth(tfg_count1, tfg_count2, max_grp_size, pct_fused, fuse_std, sparse, organisms)
    
    (x1, y1, x1_pre, y1_pre) = generate_from_linear(N1, B1, measure_noise1)
    (x2, y2, x2_pre, y2_pre) = generate_from_linear(N2, B2, measure_noise2)
    
    orths = generate_faulty_orth(orths, genes1, tfs1, genes2, tfs2, organisms, orth_falsepos, orth_falseneg)
    
    expr1 = np.vstack( (y1_pre, y1) )
    expr2 = np.vstack( (y2_pre, y2) )
    
    priors1 = generate_faulty_priors(B1, genes1, tfs1, prior_falsepos, prior_falseneg)
    priors2 = generate_faulty_priors(B2, genes2, tfs2, prior_falsepos, prior_falseneg)

    
    if out_dir == None:
        out_dir = os.tempnam('data/fake_data','gdat')
        print 'putting data in ' + out_dir
    #now we do writing to a file stuff
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    
    write_fake_td(os.path.join(out_dir, 'time_diffs1'), expr1)
    write_fake_td(os.path.join(out_dir, 'time_diffs2'), expr2)
    write_expr_mat(out_dir+os.sep+'expression1', expr1, genes1)
    write_expr_mat(out_dir+os.sep+'expression2', expr2, genes2)
    write_priors(out_dir+os.sep+'priors1',priors1)
    write_priors(out_dir+os.sep+'priors2',priors2)
    write_tfnames(out_dir+os.sep+'tfnames1',tfs1)
    write_tfnames(out_dir+os.sep+'tfnames2',tfs2)
    #since this is simulated data the first part of B is an identity matrix. That is, the TF expressions appear out of nowhere.
    #does this make me happy? no it does not.
    
    write_network(genes1, tfs1, B1, os.path.join(out_dir, 'beta1'))
    write_network(genes2, tfs2, B2, os.path.join(out_dir, 'beta2'))



    write_orth(out_dir+os.sep+'orth', orths, organisms)
    descr_dict = {'tfg_count' : str(tfg_count1), 'tfg_count2':str(tfg_count2), 'N1':N1,'N2':N2, 'max_grp_size':max_grp_size,'pct_fused':pct_fused,'fuse_std':fuse_std,'sparse':sparse,'organism1':organisms[0],'organism2':organisms[1],'prior_falsepos':prior_falsepos,'prior_falseneg':prior_falseneg,'measure_noise1':measure_noise1,'measure_noise2':measure_noise2,'tc1':'inf','tc2':'inf'}

    with file(out_dir+os.sep+'description', 'w') as f:
        for k in sorted(descr_dict.keys()):
            f.write('%s\t%s\n' % (k, str(descr_dict[k])))
    return out_dir

#adds tf expressions to the gene expression matrix, and tfs to the list of genes. This better matches the real data, in which tfs are genes.
#returns new genes, and gene expression mat
#NOTE: this requires that genes and tfs have disjoint names.
def concat_tfs_genes(genes, tfs, x, y):
    
    expr_mat = np.hstack((x,y))
    #check if they have disjoint names
    if len(set(genes).intersection(set(tfs))):
        print 'WAKA WAKA WAKA ALERT YOU HAVE GENES AND TFS WITH THE SAME NAMES'
        raise Warning #what does this do?
    genes = tfs + genes
    return (genes, expr_mat)

#SECTION: code for writing things that may or may not be fake

#writes a td file with gap time of 1 appropriate for steady state experiments
def write_ss_td(outf, expr):
    N = expr.shape[0]
    with file(outf,'w') as f:
        for i in range(N):
            f.write('%d\t%d\t%f\n' % (i, i, 1.0))

#writes a td file with a gap time of 1, and each of N/2 conditions i preceded by condition i+N/2
def write_fake_td(outf, expr):
    N = expr.shape[0]
    with file(outf,'w') as f:
        for i in range(N/2):
            f.write('%d\t%d\t%f\n' % (i, i+N/2, 1.0))

#writes a td file with a gap time of 1 for steady state, where each condition is preceded by itself
def write_ss_td(outf, expr):
    N = expr.shape[0]
    with file(outf,'w') as f:
        for i in range(N):
            f.write('%d\t%d\t%f\n' % (i, i, 1.0))

#write an expression matrix
#format is: line 1, gene names. line 2-(N+1) expressions.
#tab delimited
def write_expr_mat(outf, expr, genes):
    f = file(outf,'w')
    head = '\t'.join(genes) + '\n' #os independent linebreaks? fuck that.
    f.write(head)
    
    for cond in range(expr.shape[0]):
        gene_expressions = map(str, expr[cond, :])
        f.write('\t'.join(gene_expressions) + '\n')
    f.close()

#writes priors
#format is: lines 1-N: tf gene sign
#signs are 'activation' 'repression' 'dunno'
#priors as input is a list of tuples containing gene/tf names
def write_priors(outf, priors, signs=None):
    f = file(outf, 'w')
    for i, prior in enumerate(priors):
        (tf, gene) = prior
        if signs == None:
            sign = 'dunno'
        else:
            sign = signs[i]
        f.write('%s\t%s\t%s\n' % (tf, gene, sign))
        #f.write('%s\t%s\t%s\n' % (tf.name, gene.name, sign))
    f.close()

#write orth
#format is: each column contains comma delimited lists of genes, each column corresponds to a species, last column is mark
#orth is pairs of one_genes
def write_orth(outf, orth, organisms):
    f = file(outf, 'w')

    def which_org(org, organisms):
        arr = map(lambda x: x == org, organisms)
        ind = [i for i, x in enumerate(arr) if x]
        return ind[0]

    for o in orth:
        if len(o.genes) > 2:
            
            print 'WARNING TRYING TO WRITE NON 1-1 ORTHOLOGY'
        gene1 = o.genes[0].name
        org1 = which_org(o.genes[0].organism, organisms)  
        gene2 = o.genes[1].name
        org2 = which_org(o.genes[1].organism, organisms)
        real = o.real

        row_str = map(lambda x: [], range(len(organisms)))
        row_str[org1].append(gene1)
        row_str[org2].append(gene2)
        f.write('\t'.join(map(lambda x: ','.join(x), row_str))+ '\t'
 + str(real) + '\n')

    f.close()

#writes tf names separated by newlines
def write_tfnames(outf, tfnames):
    f = file(outf, 'w')
    f.write('\n'.join(tfnames)+'\n')
    f.close()

#writes a network
def write_network(genes, tfs, B, outf):
        
    Bs_str_l = []
    Bs_str_l.append('\t'.join(tfs))
    for gi in range(len(genes)):
        gene = genes[gi]
        regulators = B[:, gi]
        
        Bs_str_l.append(gene +'\t'+ '\t'.join(map(str, regulators)))
    f = file(outf, 'w')
    f.write('\n'.join(Bs_str_l))
    f.close()

#loads a network
def load_network(net_fn):
    f = file(net_fn)
    fs = f.read()
    fsl = filter(len, fs.split('\n'))
    fslt = map(lambda x: x.split('\t'), fsl)
    tfs = fslt[0]
    genes = map(lambda x: x[0], fslt[1:])
    #the network is written as the transpose of the matrix we want
    net = np.zeros((len(genes), len(tfs)))
    for g in range(len(genes)):
        targets = np.array(map(float, fslt[g+1][1:]))
        net[g, :] = targets
    return (net.T, genes, tfs)

#writes real data to a standard formatted folder. ONLY RUN ONCE!
#if you need to run again, change f.write('%s\t%s\t%s\n' % (tf, gene, sign)) to f.write('%s\t%s\t%s\n' % (tf.name, gene.name, sign))

def voodoo():
#writes priors
#format is: lines 1-N: tf gene sign
#signs are 'activation' 'repression' 'dunno'
#priors as input is a list of tuples containing gene/tf names
    def write_priors_voodoo(outf, priors, signs=None):
        f = file(outf, 'w')
        for i, prior in enumerate(priors):
            (tf, gene) = prior
            if signs == None:
                sign = 'dunno'
            else:
                sign = signs[i]
            f.write('%s\t%s\t%s\n' % (tf.name, gene.name, sign))
        f.close()
    #we need to reverse the order of the orthology to be consistent with 0=subtilis 1=anthracis
    def load_orth_voodoo(orth_fn, organisms):
        f = file(orth_fn)
        fs = f.read()
        fsn = filter(len, fs.split('\n'))
        fsnt = map(lambda x: x.split('\t'), fsn)
    
        orths = []
        for o in fsnt:
            real = True
            #reversing here
            orth = fr.orthology(genes = (fr.one_gene(name=o[1],organism=organisms[1]), fr.one_gene(name=o[0], organism=organisms[0])), real = real)
            orths.append(orth)
        return orths

    sub = subt()
    sub_eu = subt_eu()
    anth = anthr()
    
    (bs_e, bs_t, bs_genes, bs_tfs) = sub.load_data()
    (bseu_e, bseu_t, bseu_genes, bseu_tfs) = sub_eu.load_data()
    (ba_e, ba_t, ba_genes, ba_tfs) = anth.load_data()

    #[bs_e, bseu_e, ba_e] = quant_normalize_mats([bs_e, bseu_e, ba_e])
    #[bs_t, bseu_t, ba_t] = quant_normalize_mats([bs_t, bseu_t, ba_t])

    def orth_sub(bs_genes, bseu_genes, bs_tfs, bseu_tfs):
        organisms = ['B_subtilis', 'B_subtilis_eu']
        orths = []
        for g in bs_genes:
            if g in bseu_genes:	
                orth = fr.orthology(genes = (fr.one_gene(name=g,organism=organisms[0]), fr.one_gene(name=g, organism=organisms[1])), real = True)
                orths.append(orth)
        return orths

    (bs_priors, bs_sign) = sub.get_priors()
    (bseu_priors, bseu_sign) = sub_eu.get_priors()
    (ba_priors, ba_sign) = anth.get_priors()
    orths = load_orth_voodoo('data/bacteria1/bs_ba_ortho_804',['B_anthracis','B_subtilis'])
    orths_sub = orth_sub(bs_genes, bseu_genes, bs_tfs, bseu_tfs)
    allorths = orths + orths_sub
    bsu_operon_to_orth(os.path.join('data','bacteria1','known_operon.download.txt'), os.path.join('data','bacteria_standard','operon'))
    #operons = load_orth('data/bacteria1/bsu_operon_orth',['B_subtilis','B_anthracis'])
    out_dir = os.path.join('data','bacteria_standard')
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    write_expr_mat(out_dir+os.sep+'expression1', bs_e, bs_genes)
    write_expr_mat(out_dir+os.sep+'expression2', ba_e, ba_genes)
    write_expr_mat(out_dir+os.sep+'expression3', bseu_e, bseu_genes)
    #ss td files
    #write_ss_td(os.path.join(out_dir, 'time_diffs1'), bs_e)
    #write_ss_td(os.path.join(out_dir, 'time_diffs2'), ba_e)
    #'real' td files
    subtilis_td()
    anthracis_td()
    subtilis_eu_td()
    
    write_priors_voodoo(out_dir+os.sep+'priors1',bs_priors, bs_sign)
    write_priors_voodoo(out_dir+os.sep+'priors2',ba_priors, ba_sign)
    write_priors_voodoo(out_dir+os.sep+'priors3',bseu_priors, bseu_sign)
    write_tfnames(out_dir+os.sep+'tfnames1',bs_tfs)
    write_tfnames(out_dir+os.sep+'tfnames2',ba_tfs)
    write_tfnames(out_dir+os.sep+'tfnames3',bseu_tfs)
    write_orth(out_dir+os.sep+'orth', allorths, ['B_subtilis', 'B_anthracis', 'B_subtilis_eu'])
    #write_orth(out_dir+os.sep+'operon',operons)
    #tc: 1/tc * log(2) = 10 minutes 
    tc = np.log(2)/10
    with file(out_dir+os.sep+'description', 'w') as f:
        f.write('organism1\t%s\norganism2\t%s\norganism3\t%s\n' % ('B_subtilis','B_anthracis','B_subtilis_eu'))
        f.write('tc1\t%f\ntc2\t%f\ntc3\t%f' %(tc,tc,tc))
    


#time series loader
#tu (time unit) scales the time series. basically arbitrary. howtoset? related to decay rate
def load_subtilis(tu=1.0):
    (e, t, gene_n, tf_n) = load_B_subtilis()
    #already normalized data
    ftd = file('data/subtilis_td')
    fss = file('data/subtilis_ss')
    ftd_sl = filter(len, ftd.read().split('\n'))
    fss_sl = filter(len, fss.read().split('\n'))
    #append rows we want as a list (or differences of rows)
    gene_exps = []
    tf_exps = []
    if use_time_series:
        for l in ftd_sl:
            (from_cond_s, to_cond_s, td_s) = l.split('\t')
            (from_cond, to_cond, td) = (int(from_cond_s), int(to_cond_s), float(td_s))
        #try and predict time difference by the mean TF expr at each end
            td_expr = e[from_cond,:]+tu*(e[to_cond,:] - e[from_cond,:])/td
            tf_expr = t[from_cond,:]#(t[to_cond,:]+t[from_cond,:])/2
            
            gene_exps.append(td_expr)
            tf_exps.append(tf_expr)
    if use_steady_state:
        for l in fss_sl:
            cond = int(l)
            tf_exps.append(t[cond, :])
            gene_exps.append(e[cond, :])
    exp_mat = np.vstack(gene_exps)
    tf_mat = np.vstack(tf_exps)
    return (exp_mat, tf_mat, gene_n, tf_n)


#makes sure that XB = Y for generated data

def verify_data_integrity(N=100, N_TF=10, N_G=10):
    out = os.path.join('data','fake_data','verify_data_integrity')
    ds.write_fake_data1(N1=N, N2=N, out_dir = 'data/fake_data/simpletest', tfg_count1=(N_TF, N_G), tfg_count2 = (N_TF, N_G), measure_noise1 = 0.0, measure_noise2 = 0.0, sparse=0.0, fuse_std = 0.0, pct_fused=1.0)
    
    (n1, g1, t1) = ds.load_network('data/fake_data/simpletest/beta1')
    (n2, g2, t2) = ds.load_network('data/fake_data/simpletest/beta2')

    d1 = ds.standard_source('data/fake_data/simpletest/', 0)
    d2 = ds.standard_source('data/fake_data/simpletest/', 1)

    err1 = ((np.dot(d1.tf_mat, n1) - d1.exp_mat)**2).mean()
    err2 = ((np.dot(d2.tf_mat, n2) - d2.exp_mat)**2).mean()
    print err1
    print err2

#returns the b subtilis time-difference file associated with the metadata file subtilis_metadata. currently there are some (38) conditions that don't appear in the metadata, and are excluded
def subtilis_td():
    import pandas
    mdata = pandas.read_table(os.path.join('data','bacteria1','subtilis_metadata'), sep=' ')
    exp = pandas.read_table(os.path.join('data','bacteria1','B_subtilis.csv'), sep=',')
    cond_names = exp.columns
    cond_names_ind = {cond_names[i] : i for i in range(len(cond_names))}
    
    next_cond = dict() #maps condition number to (next_cond, time elapsed) tuples
    for r in range(len(mdata.index)):
        #figure out what  condition number this is
        cond_ind = cond_names_ind[mdata['condName'][r]]
        if not mdata['isTs'][r]:
            next_cond[cond_ind] = (cond_ind, 1.0)
        else:
            prev = mdata['prevCol'][r]
            if not str(prev) == 'nan': #this is insane
                prev_ind = cond_names_ind[prev]
                elapsed = mdata['del.t'][r]
                next_cond[prev_ind] = (cond_ind, elapsed)
    
    with file('data/bacteria_standard/time_diffs1', 'w') as td_file:
        #let's make the conditions come out in basically the same order as they appera in the data
        
        for cond1 in range(len(cond_names)):
            
            if cond1 in next_cond:
            
                (cond2, td) = next_cond[cond1]
                td_file.write('%d\t%d\t%f\n' % (cond1, cond2, td))

def subtilis_eu_td():
    import pandas
    mdata = pandas.read_table(os.path.join('data','bacteria1','subtilis_eu_metadata'), sep=' ')
    exp = pandas.read_table(os.path.join('data','bacteria1','B_subtilis_eu.csv'), sep=',')
    cond_names = exp.columns
    cond_names_ind = {cond_names[i] : i for i in range(len(cond_names))}
    
    next_cond = dict() #maps condition number to (next_cond, time elapsed) tuples
    for r in range(len(mdata.index)):
        #figure out what  condition number this is
        cond_ind = cond_names_ind[mdata['condName'][r]]
        if not mdata['isTs'][r]:
            next_cond[cond_ind] = (cond_ind, 1.0)
        else:
            prev = mdata['prevCol'][r]
            if not str(prev) == 'nan': #this is insane
                prev_ind = cond_names_ind[prev]
                elapsed = mdata['del.t'][r]
                next_cond[prev_ind] = (cond_ind, elapsed)
    
    with file('data/bacteria_standard/time_diffs3', 'w') as td_file:
        #let's make the conditions come out in basically the same order as they appera in the data
        
        for cond1 in range(len(cond_names)):
            
            if cond1 in next_cond:
            
                (cond2, td) = next_cond[cond1]
                td_file.write('%d\t%d\t%f\n' % (cond1, cond2, td))

#loads the metadata file obtained (by hand) from http://www.google.com/url?q=http%3A%2F%2Fwww.ebi.ac.uk%2Farrayexpress%2Fexperiments%2FE-MEXP-788%2Fsamples%2F&sa=D&sntz=1&usg=AFQjCNGy0TRm9r6L5MSwwl3ogfOeQZLI7Q and converts it into condition numbers. Treats all conditions from the iron starvation experiment as steady-state
def anthracis_td():
    import pandas
    #timeseries metadata
    mdata_ts = pandas.read_table('data/bacteria1/time_diffs_ant', names = ['from','to','td'],sep=' ')
    #iron is steady state!
    
    ts_data =  pandas.read_table('data/bacteria1/Normalized_data_RMA._txt')
    ts_conds = ts_data.columns[1:]
    iron_data = pandas.read_table('data/bacteria1/normalizedgeneexpressionvalues.txt')
    iron_conds = iron_data.columns[1:]
    cond_index = dict()
    #joining expression matrices doesn't change the order of conditions
    for i, ts_cond in enumerate(ts_conds):
        cond_index[ts_cond] = i
    for i, ir_cond in enumerate(iron_conds):
        if ir_cond in cond_index:
            print('overlapping condition names!!!!!! disaster!')
        cond_index[ir_cond] = i + len(ts_conds)

    with file('data/bacteria_standard/time_diffs2','w') as td_file:
        #fill in timeseries
        for r in range(len(mdata_ts.index)):
            fr = mdata_ts['from'][r]
            to = mdata_ts['to'][r].split(',')
            td = mdata_ts['td'][r]

            fr_ind = cond_index[fr]
            to_inds = ','.join(map(lambda x: str(cond_index[x]), to))
            td_file.write('%d\t%s\t%f\n' % (fr_ind, to_inds , td))
                          
        for r in range(len(iron_conds)):
            fr = cond_index[iron_conds[r]]
            to = cond_index[iron_conds[r]]
            td = 1.0
            td_file.write('%d\t%d\t%f\n' % (fr, to , td))
            
#SECTION -----------------TH17----------------------------------------

def th17_ss(use_TFA = True):
    out_dir = 'data/Th17_standard'

    rna_df = pd.DataFrame.from_csv('data/TH17/Th17Paper_rnaseq_ratios_np.tsv',sep='\t',header=0)
    rna_genes = list(rna_df.index)

    rna_tffile_old = file('data/TH17/Th17Paper_rnaseq_tfNames.tsv')
    rna_tflist = rna_tffile_old.read().split('\n')
    rna_tffile = file('data/Th17_standard/tfnames2','w')
    rna_tfs = []
    for i in range(len(rna_tflist)):
        if rna_tflist[i] in rna_genes:
            rna_tfs.append(rna_tflist[i])
            rna_tffile.write(rna_tflist[i]+'\n')

    rna_tffile_old.close()
    rna_tffile.close()

    ma_df = pd.DataFrame.from_csv('data/TH17/Tr1_Th17_noBatch_Th17PapCut.tsv',sep='\t',header=0)
    ma_genes = list(ma_df.index)

    ma_tffile_old = file('data/TH17/Tr1_Th17_noBatch_Th17PapCut_TFs.tsv')
    ma_tflist = ma_tffile_old.read().split('\n')
    ma_tffile = file('data/Th17_standard/tfnames1','w')
    ma_tfs = []
    for i in range(len(ma_tflist)):
        if ma_tflist[i] in ma_genes:
            ma_tfs.append(ma_tflist[i])
            ma_tffile.write(ma_tflist[i]+'\n')

    ma_tffile_old.close()
    ma_tffile.close()

    priors1 = file('data/Th17_standard/priors1','w')
    priors2 = file('data/Th17_standard/priors2','w')
    priorsr = []
    priorsm = []
    dfc = pd.DataFrame.from_csv('data/TH17/th17_whole_C_cut_prcnt_0_num_tfs_28_sam_0_deseq_cut_1_Aug_8_2012_priorCut0p75.tsv',sep='\t',header=0)
    nonzerosy = dfc.apply(np.nonzero, axis=1)
    for i in range(len(nonzerosy)):
        for j in range(len(nonzerosy[i][0])):
            k = nonzerosy[i][0][j]
            if dfc.columns[k] in ma_tfs:
                if dfc.index[i] in ma_genes:
                    priors1.write(dfc.columns[k] + '\t' + dfc.index[i] + '\t' + str(abs(dfc.loc[dfc.index[i], dfc.columns[k]])) + '\n')
                    priorsm.append((fr.one_gene(dfc.columns[k], 'microrarray'), fr.one_gene(dfc.index[i], 'microarray')))
            if dfc.columns[k] in ma_genes:
                if dfc.index[i] in ma_tfs:
                    priors1.write(dfc.index[i] + '\t' + dfc.columns[k] + '\t' + str(abs(dfc.loc[dfc.index[i], dfc.columns[k]])) + '\n')
                    priorsm.append((fr.one_gene(dfc.index[i], 'microrarray'), fr.one_gene(dfc.columns[k], 'microarray')))
            if dfc.columns[k] in rna_tfs:
                if dfc.index[i] in rna_genes:
                    priors2.write(dfc.columns[k] + '\t' + dfc.index[i] + '\t' + str(abs(dfc.loc[dfc.index[i], dfc.columns[k]])) + '\n')
                    priorsr.append((fr.one_gene(dfc.columns[k], 'RNAseq'), fr.one_gene(dfc.index[i], 'RNAseq')))
            if dfc.columns[k] in rna_genes:
                if dfc.index[i] in rna_tfs:
                    priors2.write(dfc.index[i] + '\t' + dfc.columns[k] + '\t' + str(abs(dfc.loc[dfc.index[i], dfc.columns[k]])) + '\n')
                    priorsr.append((fr.one_gene(dfc.index[i], 'RNAseq'), fr.one_gene(dfc.columns[k], 'RNAseq')))
    priors1.close()
    priors2.close()

    gold1 = file('data/Th17_standard/gold1','w')
    gold2 = file('data/Th17_standard/gold2','w')
    dfk = pd.DataFrame.from_csv('data/TH17/th17_whole_K_cut_prcnt_20_num_tfs_28_sam_0_deseq_cut_0.25_Aug_8_2012_priorCut0p75.tsv',sep='\t',header=0)
    nonzerosk = dfk.apply(np.nonzero, axis=1)
    for i in range(len(nonzerosk)):
        for j in range(len(nonzerosk[i][0])):
            k = nonzerosk[i][0][j]
            if dfk.columns[k] in ma_tfs:
                if dfk.index[i] in ma_genes:                    
                    gold1.write(dfk.columns[k] + '\t' + dfk.index[i] + '\t' + str(abs(dfk.loc[dfk.index[i], dfk.columns[k]])) + '\n')
            if dfk.columns[k] in ma_genes:
                if dfk.index[i] in ma_tfs:
                    gold1.write(dfk.index[i] + '\t' + dfk.columns[k] + '\t' + str(abs(dfk.loc[dfk.index[i], dfk.columns[k]])) + '\n')
            if dfk.columns[k] in rna_tfs:
                if dfk.index[i] in rna_genes:                    
                    gold2.write(dfk.columns[k] + '\t' + dfk.index[i] + '\t' + str(abs(dfk.loc[dfk.index[i], dfk.columns[k]])) + '\n')
            if dfk.columns[k] in rna_genes:
                if dfk.index[i] in rna_tfs:
                    gold2.write(dfk.index[i] + '\t' + dfk.columns[k] + '\t' + str(abs(dfk.loc[dfk.index[i], dfk.columns[k]])) + '\n')
    gold1.close()
    gold2.close()


    rna_df2 = rna_df.T
    rna_n = np.array(rna_df2)
    rna_norm = normalize_zscore(rna_n)
    rna_df3 = pd.DataFrame(rna_norm)
    rna_df3.columns = rna_df2.columns
    rna_df3.to_csv('data/Th17_standard/expression2', sep='\t',header=True, index=False)

    rna_condfile = file('data/Th17_standard/conditions2','w')
    rna_conditions = list(rna_df2.index)
    for i in range(len(rna_conditions)):
        rna_condfile.write(rna_conditions[i]+'\n')

    ma_df2 = ma_df.T
    ma_n = np.array(ma_df2)
    ma_norm = normalize_zscore(ma_n)
    ma_df3 = pd.DataFrame(ma_norm)
    ma_df3.columns = ma_df2.columns
    ma_df3.to_csv('data/Th17_standard/expression1', sep='\t', header=True, index=False)

    ma_condfile = file('data/Th17_standard/conditions1','w')
    ma_conditions = list(ma_df2.index)
    for i in range(len(ma_conditions)):
        ma_condfile.write(ma_conditions[i]+'\n')



    tc = np.log(2)/10
    with file(out_dir+os.sep+'description', 'w') as f:
        f.write('organism1\t%s\norganism2\t%s\n' % ('microarray','RNAseq'))
        f.write('tc1\t%f\ntc2\t%f' %(tc,tc))


    orthfile = file('data/Th17_standard/orth','w')
    for i in range(len(ma_genes)):
        if ma_genes[i] in rna_genes:
            orthfile.write(ma_genes[i] + '\t' + ma_genes[i] + '\t' + 'True' + '\n')

    #write_fake_td(os.path.join(out_dir, 'time_diffs1'), np.array(ma_df2))
    #write_fake_td(os.path.join(out_dir, 'time_diffs2'), np.array(rna_df2))

    write_ss_td(os.path.join(out_dir, 'time_diffs1'), np.array(ma_df2))
    write_ss_td(os.path.join(out_dir, 'time_diffs2'), np.array(rna_df2))
