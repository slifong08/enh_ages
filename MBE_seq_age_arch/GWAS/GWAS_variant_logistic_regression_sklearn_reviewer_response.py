
#%% In[1]:

from collections import Counter

import glob

import matplotlib.pyplot as plt

import numpy as np
from numpy import mean

import os, sys
import pandas as pd
from scipy import stats
import seaborn as sns
from sklearn.model_selection import train_test_split
import statsmodels.api as sm
import statsmodels.formula.api as smf
import subprocess
from scipy.stats import norm
import seaborn as sns
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.model_selection import RepeatedStratifiedKFold, cross_val_score
from sklearn import metrics, preprocessing
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix, plot_roc_curve
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, precision_recall_curve, auc, plot_precision_recall_curve
from scipy import stats
import subprocess


colors = [ "amber", "faded green"]
palette = sns.xkcd_palette(colors)
sns.palplot(palette)
plt.rcParams.update({'font.size': 15})
sns.set_style("white")

import datetime
LAST_RUN = datetime.datetime.now()
TODAY = (datetime.date.today())
RE = "/dors/capra_lab/projects/enhancer_ages/gwas_catalog/results/for_publication/"

print("last run", LAST_RUN)


#%% pleiotropy data


multipath = "/dors/capra_lab/projects/enhancer_ages/fantom/data/multiintersect/trimmed/"
multifile = "%strimmed_all_fantom_enh_112_tissues_multiintersect_0.5_count.bed"%multipath


#%% Import species data


syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t')
syn_gen_bkgd["mrca"] = syn_gen_bkgd["mrca"].round(3)
syn_gen_bkgd["mrca_2"] = syn_gen_bkgd["mrca_2"].round(3)
syn_gen_bkgd["mya"] = syn_gen_bkgd["taxon2"].apply(lambda x: x.split(" ")[-1])

syn_gen_bkgd.head()


#%% intersect with GWAS overlaps only


gwas_path = "/dors/capra_lab/projects/enhancer_ages/gwas_catalog/data/"
gwasF = "%sgwasCatalog_2019-09-24_hg19_unique_cleaned_LDEx_p5e-8.bed" % gwas_path
outF= "%senh_x_var/FANTOM_x_gwas19_ldex.bed" % gwas_path

cmd = "bedtools intersect -a %s -b %s -wao > %s" % (multifile, gwasF, outF)
os.system(cmd)

print(outF)

#%% open pleiotropy


multi = pd.read_csv(outF, sep ='\t', header = None, usecols=[0,1,2,3,4,5,6,7,15])
multi.head()


multi.columns = ["chr_enh", "start_enh", "end_enh", "old_len", "core_remodeling",
"mrca_2", "datatype", "count_overlap", "gwas_overlap"] # rename columns
multi = multi.loc[multi.chr_enh != "chrX"]

multi.head()


#%%

multi.groupby("gwas_overlap")["datatype"].count()
#%%
"""
gwas_overlap
0    30218
1     1248
"""
#%%
multi["gwas_overlap"].unique()

#%% LOGISTIC REGRESSION FUNCTION


def logit_pvalue(model, x):
    """ Calculate z-scores for scikit-learn LogisticRegression.
    parameters:
        model: fitted sklearn.linear_model.LogisticRegression with intercept and large C
        x:     matrix on which the model was fit
    This function uses asymtptics for maximum likelihood estimates.

    https://stackoverflow.com/questions/25122999/scikit-learn-how-to-check-coefficients-significance
    """
    p = model.predict_proba(x)
    n = len(p)
    m = len(model.coef_[0]) + 1
    coefs = np.concatenate([model.intercept_, model.coef_[0]])
    x_full = np.matrix(np.insert(np.array(x), 0, 1, axis = 1))
    ans = np.zeros((m, m))
    for i in range(n):
        ans = ans + np.dot(np.transpose(x_full[i, :]), x_full[i, :]) * p[i,1] * p[i, 0]
    vcov = np.linalg.inv(np.matrix(ans))
    se = np.sqrt(np.diag(vcov))
    t =  coefs/se
    p = (1 - norm.cdf(abs(t))) * 2
    return p



def logistic_regression(Xvars, yvars, df, ):


    sid = "+".join(Xvars) # create a unique sample id
    sid = "logit"
    X = df[Xvars].to_numpy()#.reshape(-1, 1)
    y = df[yvars].to_numpy().ravel()
    print(X.shape, y.shape)


    # if training/test sets are large enough
    model = LogisticRegression(max_iter = 5000, class_weight = "balanced")
    # 5-fold cross validation
    #model = LogisticRegressionCV(cv=5, random_state=0, max_iter = 1000)
    cv = RepeatedStratifiedKFold(n_splits=10, n_repeats=1, random_state=1)
    # define evaluation procedure
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    # evaaluate model
    scores = cross_val_score(model, X, y, scoring='roc_auc', cv=cv, n_jobs=-1)
    # summarize performance
    print('Mean ROC AUC: %.3f' % mean(scores))

    tprs = []
    aucs = []
    prs = []
    pr_aucs = []
    mean_recall = np.linspace(0, 1, 100)
    mean_baseline = []
    mean_fpr = np.linspace(0, 1, 100)

    fig, (ax, ax2) = plt.subplots(ncols = 2, figsize = (12,6)) # for roc_auc

    for i, (train, test) in enumerate(cv.split(X, y)):
        model.fit(X[train], y[train])
        viz = plot_roc_curve(model, X[test], y[test],
                             name='ROC fold {}'.format(i),
                             alpha=0.3, lw=1, ax=ax)
        interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr) # interpolate mean_fpr
        interp_tpr[0] = 0.0

        tprs.append(interp_tpr)
        aucs.append(viz.roc_auc)

        y_pred = model.predict(X[test]) # make predictions on test
        y_probs = model.predict_proba(X[test])[:,1]
        coef = model.coef_[0]

        print (coef)
        #print("pvalues", logit_pvalue(model, X_train))
        print("pvalues", logit_pvalue(model, X))

        cm = metrics.confusion_matrix(y[test], y_pred) # confusion matrix

        precision, recall, _ = precision_recall_curve(y[test],y_probs)
        prs.append(np.interp(mean_recall, precision, recall)) # add predictions
        pr_auc = auc(recall, precision) # calculate aucs
        pr_aucs.append(pr_auc) # append to pr_aucs
        prcurve = plot_precision_recall_curve(model, X[test], y[test],\
        name='PR fold {}'.format(i),
        alpha=0.3, lw=1, ax=ax2)

        baseline = cm[1].sum()/cm.sum() # calculate baseline as true /total pos.
        mean_baseline.append(baseline)

        ax2.set(xlabel = "Recall", ylabel= "Precision",\
        title = "PR - %s" % (sid))
        ax2.legend(loc="upper right")
        ax2.axhline(baseline, ls = '--')#, color ='k')

    #summarize roc_auc
    ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
            label='Chance', alpha=.8)

    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    ax.plot(mean_fpr, mean_tpr, color='b',
            label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
            lw=2, alpha=.8)

    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                    label=r'$\pm$ 1 std. dev.')

    ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05],
           title="auROC %s" % (sid))
    ax.legend(bbox_to_anchor =(1,-1))

    # summarize pr

    mean_precision = np.mean(prs, axis=0)
    mean_auc = auc(mean_recall, mean_precision)
    std_auc = np.std(pr_aucs)
    ax2.plot(mean_precision, mean_recall, color='navy',
             label=r'Mean (AUCPR = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
             lw=2)
    ax2.axhline(np.mean(mean_baseline), ls = "--", color ="r",\
      label='mean_baseline')

    ax2.set(xlim = [-0.05, 1.05], ylim = [-0.05, 1.05],
    xlabel = 'Recall', ylabel = 'Precision', title="PR %s" % (sid))
    ax2.legend(bbox_to_anchor =(1,-1))

    plt.show()

    return fig


#%%


multi["constant"] = 1 # add intercept

multi[["old_len", "core_remodeling", "count_overlap" ]] = multi[["old_len", "core_remodeling", "count_overlap" ]].astype(int)
multi["mrca_2"] = multi["mrca_2"].astype(float)
Xvars = ["old_len", "core_remodeling", "count_overlap", "mrca_2"]
yvars = ['gwas_overlap']
fig = logistic_regression(Xvars, yvars, multi)
