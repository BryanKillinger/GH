import numpy as np
import pandas as pd
#For multiple t-test comparisons
import scipy.stats
#For Louvain algorithm clustering
#import louvain
#import igraph as ig
from os import listdir
from os.path import isfile, join
#For multiple hypothesis testing
from statsmodels.sandbox.stats.multicomp import multipletests
import copy #For generating datasets that are modified differently for different analyses
#pvals = [random.random() for _ in range(10)]
#is_reject, corrected_pvals, _, _ = multipletests(pvals, alpha=0.1, method='fdr_bh')
#from sklearn.metrics import jaccard_similarity_score
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri

pandas2ri.activate()
dunn = importr('ClassComparison') #Need to have the ClassComparison
obase = importr('oompaBase') #For obtaining the oompa summary method for the BUM class
base = importr('base') #For accessing individual values in the BUM class

#%%
def import_data():
    peptide_path = './/Data//Full_peptides_locus.csv'
    metadata_path = './/Data//metadata.csv'
    experiment_names_path = './/Data//Experiment_names.csv'
    annotations_path = './/Data/protein_annotations.csv'
    comparison_path = './/Data/comparison_table.csv'
    #Import the crazy large file
    full_peptides = pd.read_csv(peptide_path)
    metadata = pd.read_csv(metadata_path)
    annotations = pd.read_csv(annotations_path)
    comparison_table = pd.read_csv(comparison_path)
    experiment_names = pd.read_csv(experiment_names_path, index_col=0)
    peptide_table = pd.pivot_table(full_peptides, values = 'MT_Abundance_Unscaled', index = ['Locus', 'Sequence'], columns = 'Sample_Name')

    return peptide_table, metadata, experiment_names, annotations, comparison_table

#%%
def sub_peptide_data(peptide_table, experiment_names):
    data_key = {}
    for row in experiment_names.itertuples():
        pd_row = pd.Series(row)
        filtered = pd_row.dropna()
        data_key[filtered[0]] = peptide_table[list(filtered[1:])]
    
    return data_key

#Multiply all values by the ratio of the mean of means to the column's mean to normalize b/w reps
def normalize_reps(peptides):
    for key in peptides:
        rep_mean = np.mean(peptides[key])
        mean_of_reps = np.mean(rep_mean)
        ratios = mean_of_reps/rep_mean
        peptides[key] = peptides[key]*ratios
    return peptides

def generate_data(num_samples, desired_mean, desired_std, samples, protein_name):
    np.random.seed(42)
    random_data = np.random.normal(loc=0.0, scale=desired_std, size=num_samples)
    new_data = (random_data - np.mean(random_data)) *        \
                (desired_std/(np.std(random_data-np.mean(random_data)))) + desired_mean
    formatted_data = pd.Series(new_data, index=samples)
    formatted_data.name = protein_name
    return formatted_data

def fill_data(probe_data, control_data, n_min_probe, n_min_control):

    '''
    1. Fill in missing data.
    2. Calculate q-values.
    3. Calculate fold-changes.
    4. Isolate protein targets.
    5. Place targets into dictionary with the name of the comparison.
    
    1. Fill in missing data.
        1. Remove any proteins where the probe has less than 3 datapoints. Need at least 3 datapoints in the probe sample to move forward.
        2. If a control sample is completely missing data, populate with 50% the minimum value for that protein across all datasets with a standard deviation of the probe protein for that sample.
        3. If a sample has only 1 datapoint, fill the missing values in with random datapoints with the mean of the present datapoint and a standard deviation of the probe sample for that same protein.
        4. If a sample has 2 datapoints, fill in the missing values with datapoints with a mean of the present datapoints and their standard deviation.

    '''
    #Errors for incorrect data.
    '''
    if len(probe_data.index) != len(control_data.index):
        raise Exception('ERORR! Number of rows in probe data does not equal rows in control data.')
    if len(control_data.columns) < n_min_control:
         raise Exception('Need at least {} control samples for {}'.format(n_min_control, comp[1]['Control']))
    if len(probe_data.columns) < n_min_probe:
         raise Exception('Need at least {} control samples for {}'.format(n_min_control, comp[1]['Probe']))
'''
    probe_rep_mins = np.min(probe_data)
    np_rep_mins = np.min(probe_data)
    
    if n_min_control >= n_min_probe:
        raise Exception('n_min_control must be smaller than n_min_probe')

    new_probe_protein = []
    new_control_protein = []
    probe_names = []
    control_names = []
    protein_list = probe_data.index
    extra_probes = []
    not_in_np_probes = []

#Actual analysis
    
    #Assume data coming in is of the correct form. No single replicate samples.
    for protein in protein_list:
        probe_protein = probe_data.loc[protein]
        control_protein = control_data.loc[protein]
        
        #Number of datapoints for this protein
        n_probe_protein = len(probe_protein) - np.isnan(probe_protein).sum()
        n_control_protein = len(control_protein) - np.isnan(control_protein).sum()
        
        #If all control protein replicates are present
        if (n_control_protein == len(control_protein)):

            
            #if n_min_control  is 2, then the probe will be added to probe list. Otherwise, leave it off the standard list and add it to the not_in_np_probes list.
            #n_min_control is 2
            #n_min_probe is 3
            #So here, if n_probe_protein = 2, generate data for it.
            if (n_probe_protein >= n_min_control) & (n_probe_protein < n_min_probe): 
                probe_protein = generate_data(len(probe_protein), np.mean(probe_protein), np.std(probe_protein), probe_data.columns, protein)
                new_probe_protein.append(probe_protein)
                probe_names.append(protein)
                extra_probes.append(protein)
                new_control_protein.append(control_protein)
                control_names.append(protein)
            elif (n_probe_protein == 1):
                probe_protein = generate_data(len(probe_protein), np.mean(probe_protein), np.std(control_protein), probe_data.columns, protein)
                new_probe_protein.append(probe_protein)
                probe_names.append(protein)
                new_control_protein.append(control_protein)
                control_names.append(protein)
            elif (n_probe_protein == 0):
                probe_protein = generate_data(len(probe_protein), np.mean(probe_rep_mins), np.std(control_protein), probe_data.columns, protein)
                new_probe_protein.append(probe_protein)
                probe_names.append(protein)
                new_control_protein.append(control_protein)
                control_names.append(protein)

            '''   
            elif (n_probe_protein < n_min_control):
                not_in_np_probes.append(protein) #These are probe peptides that have 0 or 1 data points, but 3 data points for the corresponding no_probe peptides.
            ''' 

        #n_min_probe is 3, so if all 3 data points are present for the probe:
        if (n_probe_protein >= n_min_probe):
            
            #Generate data for control_protein if missing data points
            if n_control_protein == 0:
                control_protein = generate_data(len(control_protein), np.mean(np_rep_mins), np.std(probe_protein), control_data.columns, protein)
                new_control_protein.append(control_protein)
                control_names.append(protein)
            
            elif n_control_protein == 1:
                control_protein = generate_data(len(control_protein), np.mean(control_protein), np.std(probe_protein), control_data.columns, protein)
                new_control_protein.append(control_protein)
                control_names.append(protein)
            
            elif (n_control_protein > 1) & (n_control_protein < len(control_protein)):
                control_protein = generate_data(len(control_protein), np.mean(control_protein), np.std(control_protein), control_data.columns, protein)
                new_control_protein.append(control_protein)
                control_names.append(protein)
                
            elif (n_control_protein == len(control_protein)):
                new_control_protein.append(control_protein)
                control_names.append(protein)

            #If there are missing data points for the probe_protein, generate data.
            #This shouldn't ever occur for the current data set.
            if n_probe_protein != len(probe_protein):
                probe_protein = generate_data(len(probe_protein), np.mean(probe_protein), np.std(probe_protein), probe_data.columns, protein)
                
            
            #Append new protein data to list
            new_probe_protein.append(probe_protein)
            probe_names.append(protein)
        
        
    probe_index = pd.MultiIndex.from_tuples(probe_names, names=['RefID', 'Peptide'])
    control_index = pd.MultiIndex.from_tuples(control_names, names=['RefID', 'Peptide'])
    return pd.DataFrame(new_probe_protein, index=probe_index), pd.DataFrame(new_control_protein, index=control_index), probe_names, extra_probes, not_in_np_probes

#%%
#If t-statistic is negative, then the np is larger than the probe. This isn't a target. The real 1-side p-value is 1 - p for these values.
def check_for_negative(t_statistic, p_value):
    if t_statistic < 0:
        return 1 - p_value
    else:
        return p_value

#%%
#From the two samples, generate a t-test with q-values for each and fold changes
    
def target_statistics(probe_data, control_data, protein_list, one_sided):
    t_test = scipy.stats.ttest_ind(probe_data, control_data, axis=1, equal_var=False)
    if (one_sided == 1):
        p_vals = list(map(check_for_negative, t_test[0], t_test[1]))
    else:
        p_vals = list(t_test[1])
    is_reject, corrected_pvals, _, _ = multipletests(p_vals, alpha=0.05, method='fdr_bh')
    q = pd.Series(corrected_pvals, index = protein_list)
    p = pd.Series(p_vals, index = protein_list)
    #Check if t-statistic is -. If it is, then replace the p-value with 1 - p.
    
    
    return q, p
    '''
    Will return t-statistic that is positive or negative. Need the real p-value statistic. Will need to look at the t-statistic to compute the p-value. IF t is -, then 1 - p. If t is +, then p = p (or vice versa)
    '''


def group_peptides(peptides):
    proteins = {}
    for ref, mynewdef in peptides.groupby(level=0):
        proteins[ref] = mynewdef
    return proteins

# %%

def filter_proteins(probe, no_probe, protein_list, min_count, fc):
    proteins_filt = []
    for key in protein_list:
        ind_fc = (np.mean(probe[key], axis=1) - np.mean(no_probe[key], axis=1))
        if (len(ind_fc) > min_count) & (np.max(ind_fc) > fc): #Check # peptides and max fc > 2
            proteins_filt.append(key)
    return proteins_filt

#%%
def test_for_target(q_vals):
    #For a set of q_values, go through each protein and model as a mixture of uniform and beta distributions
    proteins = {}
    for ref, mynewdef in q_vals.groupby(level=0):
        #print(mynewdef)
        proteins[ref] = mynewdef 
        #AND/OR can run direct test here for the beta uniform model  
    return proteins

'''
#Select a specific protein and amino acid:
b = mydf.loc[308, 'VLDHIAILR']

#Select all peptides for a specific protein:
c = mydf.loc[308]

'''

# %%


def bumfit(p_vals, tau):
        b = pandas2ri.py2ri(pd.Series(p_vals))
        c = dunn.Bum(b)
        d = obase.summary(c, tau)
        at_sym = base.__dict__["@"]
        estimates = at_sym(d, "estimates")
        bum = at_sym(d, "bum")
        #pvals = at_sym(bum, "pvals")
        ahat = pandas2ri.ri2py(at_sym(bum, "ahat"))
        lhat = pandas2ri.ri2py(at_sym(bum, "lhat"))
        pihat = pandas2ri.ri2py(at_sym(bum, "pihat"))
        q = pandas2ri.ri2py(estimates)
        p = pd.DataFrame([ahat, lhat, pihat], index=['ahat', 'lhat', 'pihat'])
        return q, ahat

#p_vals = [0.01, .02, .05, .95, .23, .35, .235, .32, .35, .1, .05, .001]
#p_vals = [.01, .05, .03, .1, .5, .01, .001, .6]


def calc_sig(q_vals, alpha):
    bum_estimates, a = bumfit(q_vals, alpha) 
    sig = np.float((bum_estimates['TP']/(bum_estimates['FN']))/(bum_estimates['FP'] + bum_estimates['TN'])) 
    tn_tp = np.float(bum_estimates['TN']/bum_estimates['TP'])
    fdr = np.float(bum_estimates['FP']/(bum_estimates['FP'] + bum_estimates['TP']))
    tp = np.float(bum_estimates['TP'])
    avg_q = np.average(q_vals)
    beta_uniform = np.float((bum_estimates['TP'] + bum_estimates['FN'])/(bum_estimates['FP'] + bum_estimates['TN']))

    return sig, tn_tp, fdr, tp, avg_q, a, beta_uniform

def protein_looper(protein_dict, tau, false_discovery, true_positive, tn_tp_thresh, a_thresh):
    sig_proteins = {}
    nonsig_proteins = {}
    for protein in protein_dict:
        if len(protein_dict[protein]) > 1: #Another test for how many peptides are present
            proteins, tn_tp, fdr, tp, q, a, bu = calc_sig(protein_dict[protein], tau)
            if (a < a_thresh) & (tn_tp < tn_tp_thresh):
                sig_proteins[protein] = [proteins, tn_tp, fdr, tp, q, a, bu]
            else:
                nonsig_proteins[protein] = [proteins, tn_tp, fdr, tp, q, a, bu]
    return sig_proteins, nonsig_proteins


# %%
def find_targets(peptide_probe, peptide_np, peptide_list, protein_annotations):
    grouped_peptide_probe = group_peptides(peptide_probe)
    grouped_peptide_np = group_peptides(peptide_np)
    protein_list = list(grouped_peptide_probe.keys())
    proteins = filter_proteins(grouped_peptide_probe, grouped_peptide_np, protein_list, 3, 2)
    probe_pass_peptides = peptide_probe.loc[proteins]
    np_pass_peptides = peptide_np.loc[proteins]
    peptide_list = list(probe_pass_peptides.index)
    added_q_vals, p_vals = target_statistics(probe_pass_peptides, np_pass_peptides, peptide_list, 1)
    added_q_vals.index = pd.MultiIndex.from_tuples(added_q_vals.index, names = ['RefID', 'Peptide'])
    
     
    proteins = test_for_target(added_q_vals)
    all_vals, nonsig_vals = protein_looper(proteins, .05, .1, .1, 1, .4)
    #tau, fdr, TP, TP/FN, ahat
    qq = pd.DataFrame.from_dict(all_vals, orient='index')
    qq.columns = ['Sig', 'tn_tp', 'fdr', 'tp', 'q', 'ahat', 'beta_uniform']
    qq_cutoff = qq#[qq['Sig'] > 1]
    qq_annotations = protein_annotations[protein_annotations['LOCUS'].isin(qq_cutoff.index)]
    qq_annotations.set_index('LOCUS', inplace=True)
    qq_together = pd.concat([qq_cutoff, qq_annotations], axis=1)
    qq_together.sort_values('Sig', ascending=False, inplace=True)
    protein_targets = list(qq_cutoff.index)
    qq_np = peptide_np.loc[protein_targets]
    qq_p = peptide_probe.loc[protein_targets]
    #np_proteins = {}
    """
    for ref, mynewdef in qq_np.groupby(level=0):
            np_proteins[ref] = mynewdef 
    p_proteins = {}
    for ref, mynewdef in qq_p.groupby(level=0):
            p_proteins[ref] = mynewdef 
    """
    
    #Needs to be only probe targets (qq_p), and only the probe targets that have alpha < .05.
    
    significant_peptides = added_q_vals[added_q_vals < .05].index #All q's less than .05
    probe_peptides = qq_p.index #Indices of all peptides from proteins that made cutoff
    short_peptide_list = significant_peptides.intersection(probe_peptides) #Intersection
    short_peptide_q_values = added_q_vals.loc[short_peptide_list] #q-vals of intersection
    
    lst = {k: proteins[k] for k in protein_targets}
    return qq_together, lst, qq_p, qq_np, short_peptide_q_values, proteins, nonsig_vals

#%%
def find_shared_proteins(list_1, list_2):
    return list_1.intersection(list_2)

def find_shared_proteins_dict(protein_dict, annotations):
    shared_list = pd.DataFrame()
    unique_set = set()
    for keys in protein_dict:
        shared_list[keys] = pd.Series(protein_dict[keys].index)
        if len(unique_set) == 0:
            unique_set = shared_list[keys]
        else:
            unique_set = set(shared_list[keys]).intersection(unique_set)
    shared_proteins = list(unique_set)
    shared_annotations = annotations[annotations['LOCUS'].isin(list(set(shared_proteins)))]
    
    return shared_annotations, shared_list



# %%
#a = find_shared_proteins(sample_targets)

def list_all_targets(target, protein_annotations):
    targets_list = []
    for keys in target:
        targets_list.extend(list(target[keys].index))
    target_annotations = protein_annotations[protein_annotations['LOCUS'].isin(list(set(targets_list)))]
    return target_annotations

# %%
peptide_table, metadata, experiment_names, protein_annotations, comparison_table = import_data() #import data
peptide_table = np.log2(peptide_table) #log transform
unscaled_peptides = sub_peptide_data(peptide_table, experiment_names)
peptides = normalize_reps(unscaled_peptides)

def target_loop():
    targets = {}
    target_peptides = {}
    sig_peptides = {}
    q_values_of_probe_targets = {}
    peptides_of_probe_targets = {}
    prot_stats = {}
    nonsig_prot_stats = {}
    for comparison in comparison_table.iterrows():
        probe_comp = comparison[1].loc['Probe']
        np_comp = comparison[1].loc['NoProbe']
        
        peptide_probe, peptide_np, peptide_list, ex, no_np = fill_data(peptides[probe_comp], peptides[np_comp], 3, 2)
        
        target_list, prot_list, probe_peps, noprobe_peps, sig_peps, ind_stats, nonsig_prots = find_targets(peptide_probe, peptide_np, peptide_list, protein_annotations)
        
        nonsig_prot_stats[probe_comp] = nonsig_prots
        prot_stats[probe_comp] = ind_stats
        targets[probe_comp] = target_list
        target_peptides[probe_comp] = peptide_probe
        sig_peptides[probe_comp] = sig_peps
        peptides_of_probe_targets[probe_comp] = probe_peps
  
    all_targets = list_all_targets(targets, protein_annotations)
    shared_targets, targets_df = find_shared_proteins_dict(targets, protein_annotations)
      
    return targets, all_targets, peptides, comparison_table, unscaled_peptides, target_peptides, shared_targets, sig_peptides, q_values_of_probe_targets, prot_stats, nonsig_prot_stats

#USING PEPTIDES_OF_PROBE_TARGETS for actual target peptides

sample_targets, all_targets, peptides, comparison_table, unscaled_peptides, target_peptides, shared_targets, sig_peptides, q_values_of_probe_targets, prot_stats, nonsig_prot_stats = target_loop()


#%% 
#All targets that are shared across every condition
shared_targets_dict, targets_df = find_shared_proteins_dict(sample_targets, protein_annotations)

#%%Targets
#Compare activity between probe targets in probe samples
target_comparison_path = './/Data/target_comparisons.csv'
target_comparisons = pd.read_csv(target_comparison_path)


avg_q_vals = {}
avg_fc_vals = {}
comp1 = {}
comp2 = {}
target_comp = {}
target_comp_neg = {}
proteins_vs = {}
all_together = {}
all_mine = {}
all_sig = {}
protein_q_vals = {}
activity_all_vals = {}
activity_nonsig_vals = {}
bum_activity = {}
all_bum = {}
bum_activity_fc = {}
probe_comp_table = pd.DataFrame(columns=['Comp', 'ID'])
fold_change = {}
fold_change_sign = {}
sign_avg = pd.Series()


proteins_vs_neg = {}
activity_all_vals_neg = {}
activity_nonsig_vals_neg = {}
bum_activity_pos = {}
bum_activity_neg = {}

for row in target_comparisons.iterrows():

    avg_q = pd.Series()
    avg_fc = pd.Series()
    together = pd.DataFrame(columns=['Activity_Q_value', 'Activity_Fold_Change'])
    
    comp = row[1][0] + '_vs_' + row[1][1]
    new_comp = pd.Series([comp, row[1]['ID']], index=['Comp', 'ID'])
    probe_comp_table = probe_comp_table.append(new_comp, ignore_index=True)
    
    
    
    
    
    #target_peptides is all peptides that are filled in (COULD HAVE BEEN IMPUTED). These peptides are then subsetted with the sig_peptides, which are peptides used in the p vs. np comparison, that also had a q-value of less than .05. 
    
    #Might want to only subset peptides that are NOT filled in, raw data. Each should have all reps present. But might want to keep since these peptides are only with z-values less than .05, might be good to keep filled in data.
    comp1_peptides = peptides[row[1][0]].loc[peptides[row[1][0]].index]
    comp2_peptides = peptides[row[1][1]].loc[peptides[row[1][1]].index]
    
    #But we want only peptides that are deemed targets across ANY condition, and get ALL of them
    comp1_peptides = comp1_peptides.loc[all_targets["LOCUS"]]
    comp2_peptides = comp2_peptides.loc[all_targets["LOCUS"]]
    

    #Now fill in missing data between the two samples the same way as for probe vs np:
    comp1_peptides, comp2_peptides, peptide_list, ex, no_np = fill_data(comp1_peptides, comp2_peptides, 3, 2)
    
    
    
    """
    a = comp1_peptides.index.intersection(comp2_peptides.index)
    b = list(a)
    comp1[comp] = comp1_peptides.loc[b]
    comp2[comp] = comp2_peptides.loc[b]
   """
    comp1[comp] = comp1_peptides
    comp2[comp] = comp2_peptides
    
    #Q-test between comp1 and comp2
    target_comp[comp] = target_statistics(comp1[comp], comp2[comp], comp1[comp].index, 1)[0]
    target_comp_neg[comp] = target_statistics(comp2[comp], comp1[comp], comp1[comp].index, 1)[0]
    
    #Roll up to a dictionary with proteins as keys with q-value series inside each protein.
    proteins_vs[comp] = test_for_target(target_comp[comp])
    proteins_vs_neg[comp] = test_for_target(target_comp_neg[comp])
    
    #Next need to obtain average q-values or model with BUM. Perhaps want to model the q-values similarily to the probe vs np to get an accurate comparison, but keep it two tailed. Use BUM fit thresholds as well as average fold-change. Need to include average fold change because we must use a two-tailed test since either sample could have a a higher value, requiring two-tailed test. So we need to make sure we aren't saying something is different based just on q-value distribution. Want absolute value of fold change to be greater than threshold (2).
    
    #Will need to run "find_targets", perhaps want to include average FC next to the average q-value in find_targets, but not absolutely necessary since we calculate it below.

    activity_all_vals[comp], activity_nonsig_vals[comp] = protein_looper(proteins_vs[comp], .05, .1, .1, .5, .4)
    activity_all_vals_neg[comp], activity_nonsig_vals_neg[comp] = protein_looper(proteins_vs_neg[comp], .05, .1, .1, 1, .4)
    #tau, fdr, TP, TP/FN, ahat
    bum_activity_pos[comp] = pd.DataFrame.from_dict(activity_all_vals[comp], orient='index')
    bum_activity_neg[comp] = pd.DataFrame.from_dict(activity_all_vals_neg[comp], orient='index')
    
    
    bum_activity[comp] = bum_activity_pos[comp].append(bum_activity_neg[comp])
    
    
    if len(bum_activity[comp] > 0): #If there is at least 1 protein that meets the BUM criteria
        bum_activity[comp].columns = ['Sig', 'tn_tp', 'fdr', 'tp', 'q', 'ahat', 'beta_uniform']
        fold_change[comp] = pd.Series()      
        for protein in bum_activity[comp].iterrows():
            print(protein[0])
            
            fold_change[comp][protein[0]] = np.average(np.array(comp1[comp].loc[protein[0]]) - np.array(comp2[comp].loc[protein[0]]))
            fc_sign = np.average(np.array(comp1[comp].loc[protein[0]]) - np.array(comp2[comp].loc[protein[0]]), axis=1)
            fc_sign[fc_sign > 0] = 1
            fc_sign[fc_sign < 0] = -1
            fc_sign[fc_sign == 0] = 0
            sign_avg[protein[0]] = np.average(fc_sign)
            
        bum_activity[comp]['Fold_Change_Sign_Average'] = sign_avg
        bum_activity[comp]['Fold_Change'] = fold_change[comp]
        bum_annotations = protein_annotations[protein_annotations['LOCUS'].isin(bum_activity[comp].index)]
        bum_annotations.set_index('LOCUS', inplace=True)
        all_bum[comp] = pd.concat([bum_activity[comp], bum_annotations], axis=1)     
        bum_activity_fc[comp] = all_bum[comp][np.abs(all_bum[comp]['Fold_Change_Sign_Average']) > .5]
        
#%%Globals
for row in target_comparisons.iterrows():

    avg_q = pd.Series()
    avg_fc = pd.Series()
    together = pd.DataFrame(columns=['Activity_Q_value', 'Activity_Fold_Change'])
    
    comp = row[1][0] + '_vs_' + row[1][1]
    new_comp = pd.Series([comp, row[1]['ID']], index=['Comp', 'ID'])
    probe_comp_table = probe_comp_table.append(new_comp, ignore_index=True)
    
    
    
    
    
    #target_peptides is all peptides that are filled in (COULD HAVE BEEN IMPUTED). These peptides are then subsetted with the sig_peptides, which are peptides used in the p vs. np comparison, that also had a q-value of less than .05. 
    
    #Might want to only subset peptides that are NOT filled in, raw data. Each should have all reps present. But might want to keep since these peptides are only with z-values less than .05, might be good to keep filled in data.
    comp1_peptides = peptides[row[1][0]].loc[peptides[row[1][0]].index]
    comp2_peptides = peptides[row[1][1]].loc[peptides[row[1][1]].index]
    
    #But we want only peptides that are deemed targets across ANY condition, and get ALL of them
    comp1_peptides = comp1_peptides.loc[all_targets["LOCUS"]]
    comp2_peptides = comp2_peptides.loc[all_targets["LOCUS"]]
    

    #Now fill in missing data between the two samples the same way as for probe vs np:
    comp1_peptides, comp2_peptides, peptide_list, ex, no_np = fill_data(comp1_peptides, comp2_peptides, 3, 2)
    
    
    
    """
    a = comp1_peptides.index.intersection(comp2_peptides.index)
    b = list(a)
    comp1[comp] = comp1_peptides.loc[b]
    comp2[comp] = comp2_peptides.loc[b]
   """
    comp1[comp] = comp1_peptides
    comp2[comp] = comp2_peptides
    
    #Q-test between comp1 and comp2
    target_comp[comp] = target_statistics(comp1[comp], comp2[comp], comp1[comp].index, 1)[0]
    target_comp_neg[comp] = target_statistics(comp2[comp], comp1[comp], comp1[comp].index, 1)[0]
    
    #Roll up to a dictionary with proteins as keys with q-value series inside each protein.
    proteins_vs[comp] = test_for_target(target_comp[comp])
    proteins_vs_neg[comp] = test_for_target(target_comp_neg[comp])
    
    #Next need to obtain average q-values or model with BUM. Perhaps want to model the q-values similarily to the probe vs np to get an accurate comparison, but keep it two tailed. Use BUM fit thresholds as well as average fold-change. Need to include average fold change because we must use a two-tailed test since either sample could have a a higher value, requiring two-tailed test. So we need to make sure we aren't saying something is different based just on q-value distribution. Want absolute value of fold change to be greater than threshold (2).
    
    #Will need to run "find_targets", perhaps want to include average FC next to the average q-value in find_targets, but not absolutely necessary since we calculate it below.

    activity_all_vals[comp], activity_nonsig_vals[comp] = protein_looper(proteins_vs[comp], .05, .1, .1, .5, .4)
    activity_all_vals_neg[comp], activity_nonsig_vals_neg[comp] = protein_looper(proteins_vs_neg[comp], .05, .1, .1, 1, .4)
    #tau, fdr, TP, TP/FN, ahat
    bum_activity_pos[comp] = pd.DataFrame.from_dict(activity_all_vals[comp], orient='index')
    bum_activity_neg[comp] = pd.DataFrame.from_dict(activity_all_vals_neg[comp], orient='index')
    
    
    bum_activity[comp] = bum_activity_pos[comp].append(bum_activity_neg[comp])
    
    
    if len(bum_activity[comp] > 0): #If there is at least 1 protein that meets the BUM criteria
        bum_activity[comp].columns = ['Sig', 'tn_tp', 'fdr', 'tp', 'q', 'ahat', 'beta_uniform']
        fold_change[comp] = pd.Series()      
        for protein in bum_activity[comp].iterrows():
            print(protein[0])
            
            fold_change[comp][protein[0]] = np.average(np.array(comp1[comp].loc[protein[0]]) - np.array(comp2[comp].loc[protein[0]]))
            fc_sign = np.average(np.array(comp1[comp].loc[protein[0]]) - np.array(comp2[comp].loc[protein[0]]), axis=1)
            fc_sign[fc_sign > 0] = 1
            fc_sign[fc_sign < 0] = -1
            fc_sign[fc_sign == 0] = 0
            sign_avg[protein[0]] = np.average(fc_sign)
            
        bum_activity[comp]['Fold_Change_Sign_Average'] = sign_avg
        bum_activity[comp]['Fold_Change'] = fold_change[comp]
        bum_annotations = protein_annotations[protein_annotations['LOCUS'].isin(bum_activity[comp].index)]
        bum_annotations.set_index('LOCUS', inplace=True)
        all_bum[comp] = pd.concat([bum_activity[comp], bum_annotations], axis=1)     
        bum_activity_fc[comp] = all_bum[comp][np.abs(all_bum[comp]['Fold_Change_Sign_Average']) > .5]

#%%
    """
    #Use average q-value and average fold-change.
    for key in proteins_vs[comp]:
        avg_q[key] = np.average(proteins_vs[comp][key], axis=0)
        avg_fc[key] = np.average(np.array(comp1[comp].loc[key]) - np.array(comp2[comp].loc[key]))
        
    #Also need average fold change for each protein
    together['Activity_Q_value'] = avg_q
    together['Activity_Fold_Change'] = avg_fc
    annotations = protein_annotations[protein_annotations['LOCUS'].isin(together.index)]
    annotations.set_index('LOCUS', inplace=True)
    avg_q_vals[comp]= avg_q
    avg_fc_vals[comp] = avg_fc
    all_together[comp] = together
    all_mine[comp] = pd.concat([all_together[comp], annotations], axis=1)
    
    #Isolate Q-values less than .05
    all_sig[comp] = all_mine[comp][all_mine[comp]['Activity_Q_value'] < .05]
    
    #Can isolate proteins that are from shared organisms, since of course any organism present in one sample vs another it isn't in will have more protein targets. Could just keep, but must be aware.
    """

    
#%% Make the above into a function and use it to compare probe vs np (select single peptide for fold_change), activities, and globals. Reason for a single peptide fold-change is to reduce sample complexity. For comparing activities and globals, we are using a two-tailed test since we want to see if either is larger. A fold-change allows us to discriminate random noise from a single protein whose peptides lie on both sides of the distribution.
'''
Comparison of global values
For each comparison, fill in data using normal method.
Should model using BUM distribution as well
'''
global_comparison_table = pd.read_csv('.//Data//global_comparison_table_probe_only.csv')
global_comps = {}
global_fc = {}
global_comps_fc_q = {}

def compare_globals(global_comparison_table):
    comp_table = pd.DataFrame(columns=['Comp', 'ID'])
    for global_comparison in global_comparison_table.iterrows():
        comp = global_comparison[1][0] + '_vs_' + global_comparison[1][1]
        
        #pd.concat([comp_table, pd.Series(comp, global_comparison[1]['Index'])])
        new_data = pd.Series([comp, global_comparison[1]['ID']], index=['Comp', 'ID'])
        comp_table = comp_table.append(new_data, ignore_index=True)

        #comp_table['Comp'].append(pd.Series(comp), inplace=True)
        #comp_table['ID'].append(pd.Series(global_comparison[1]['Index']), inplace=True)
        
        print(global_comparison)
        global_comp_1 = global_comparison[1].loc['Global_Comp_1']
        global_comp_2 = global_comparison[1].loc['Global_Comp_2']
        
        
        peptide_comp_1, peptide_comp_2, peptide_list, ex, no_np = fill_data(peptides[global_comp_1].loc[all_targets['LOCUS']], peptides[global_comp_2].loc[all_targets['LOCUS']], 3, 2)
        
        
        #Add in BUM comparisons here.
        
        
        
        all_fc_1 = peptide_comp_1.apply(np.average, axis=1)
        all_fc_2 = peptide_comp_2.apply(np.average, axis=1)
        avg_fc = all_fc_1 - all_fc_2
        global_fc = test_for_target(avg_fc)
        
        peptide_q, peptide_p = target_statistics(peptide_comp_1, peptide_comp_2, peptide_list, 0)
        peptide_q.index = pd.MultiIndex.from_tuples(peptide_q.index, names = ['RefID', 'Peptide'])
        proteins = test_for_target(peptide_q)
        
        all_vals = {}
        for spec_protein in proteins:
            prots, tn_tp, fdr, tp, q, a, bu = calc_sig(proteins[spec_protein], .05)
            fc = np.average(global_fc[spec_protein])
            all_vals[spec_protein] = [prots, tn_tp, fdr, tp, q, a, bu, fc]
        

        qq = pd.DataFrame.from_dict(all_vals, orient='index')
        qq.columns = ['Sig', 'tn_tp', 'fdr', 'tp', 'q', 'ahat', 'beta_uniform', 'global_fc']
        
        global_comps[comp] = qq
        global_comps_fc_q[comp] = pd.DataFrame(columns=['global_q', 'global_fc'])
        global_comps_fc_q[comp]['global_q'] = qq['q']
        global_comps_fc_q[comp]['global_fc'] = qq['global_fc']
        
    return global_comps, global_fc, comp_table, global_comps_fc_q
        


global_comps, global_fc, global_comp_table, global_comps_fc_q = compare_globals(global_comparison_table)

#%%
new_df = {}
concat_activity_global = {}
def targets_probe_global(probe_comps, global_comps, probe_dict, global_dict):
    for comp in probe_comps.iterrows():
        numerical_id = comp[1].loc['ID'] #Extract numerical id from probe comp
        probe_comp = comp[1].loc['Comp']
        global_comp = list(global_comps[global_comps['ID'] == numerical_id]['Comp']) #Select comparison from globals
        target_index = probe_dict[probe_comp].index
        #global_index = global_dict[global_comp[0]].index.isin(target_index)
        #new_df = global_index
        #print(target_index)
        #print(type(target_index))
        globals_of_targets = global_dict[global_comp[0]].loc[global_dict[global_comp[0]].index.isin(target_index)]
        
        #new_df[numerical_id] = pd.concat([probe_dict[probe_comp], globals_of_targets], ignore_index=True)
        z = probe_dict[probe_comp]
        #concat_globals_activity = probe_dict[comp[1]]
        global_activity_annotations = all_targets[all_targets['LOCUS'].isin(z.index)]
        global_activity_annotations.index = global_activity_annotations['LOCUS']
        concat_activity_global[probe_comp] = pd.concat([z, globals_of_targets, global_activity_annotations], axis=1)
        
        #Select correct comparisons matched between probe_dict and global_dict based on index
        #Return new dataframe that concatenates both dataframe while matching indices
    return concat_activity_global
        

target_activity_global = targets_probe_global(probe_comp_table, global_comp_table, all_together, global_comps_fc_q)


'''
Verify all data is correct and what is expected manually
'''
        
    #Do essentially the exact same thing as with the probe vs np data, but keep the first t-test (as well as the second) two-sided so that it doesn't invert any p-values. This is essentially all that is required.
    
    #Next, look at biological findings. Can cluster all targets based on activity. Cluster all proteomics data based on abundance. How to actually quantify? Find shared peptide between all samples OR scale detected peptides to match those of undetected. If nothing is shared between samples, can't quantify.
    
    
#%%
#target_activity_globa['.to_csv('out5_1_hits.csv')
import scipy.stats
def correlate_peptides(peptide_matrix, method, min_periods):
    return peptide_matrix.corr(method=method, min_periods=min_periods)

#Select only target peptides, not everything
binary_conditions = pd.read_csv('.//Data/binary_sample_data.csv')
activity_samples = list(binary_conditions[binary_conditions['PROBE'] == 1]['NEW NAME'])
activity_peptides = peptide_table[activity_samples]
activity_peptides = activity_peptides.loc[all_targets['LOCUS']]
activity_peptides_t = activity_peptides.transpose()
activity_peptides_t = activity_peptides.loc[activity_peptides_t.count(axis=0) == len(activity_peptides_t)]
activity_peptides = activity_peptides_t.transpose()



#cor_act = scipy.stats.spearmanr(activity_peptides_t)
#correlation_activity_peptides = correlate_peptides(activity_peptides_t, 'spearman', 0)
#corr_act_peps = test_for_target(correlation_activity_peptides)

#OR just cluster ABP target peptides
import sklearn.cluster as sc
#%%For troubleshooting and isolating individual peptides from a specific protein
pep = "Stan_51_0054"
samp1 = "GH2B_HL69_HL91"
samp2 = "GH2B_HL69_NoB12"
t1 = peptides[samp1].loc[pep]
t2 = peptides[samp2].loc[pep]
#t3 = target_peptides[samp1].loc[pep]
#t4 = target_peptides[samp2].loc[pep]


#%%
#Targets that are not present in the comparison, with the corresponding number of peptides detected in each, along with the average q-values and fold change for the probe vs np comparison.




