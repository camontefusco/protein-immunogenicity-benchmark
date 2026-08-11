from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import GroupShuffleSplit

CANON=Path(__file__).resolve().parents[1];OUT=CANON/'results/forensic_sensitivities';OUT.mkdir(parents=True,exist_ok=True)
SEED=42

class UF:
    def __init__(self,n):self.p=list(range(n));self.rank=[0]*n
    def f(self,x):
        while self.p[x]!=x:self.p[x]=self.p[self.p[x]];x=self.p[x]
        return x
    def u(self,a,b):
        a,b=self.f(a),self.f(b)
        if a==b:return
        if self.rank[a]<self.rank[b]:a,b=b,a
        self.p[b]=a
        if self.rank[a]==self.rank[b]:self.rank[a]+=1

d=pd.read_csv(CANON/'data/splits/peptide_split_assignments_v1.csv');d=d[d.valid_sequence.astype(bool)].reset_index(drop=True)

def groups(combine_edit):
    u=UF(len(d));seen={};cluster_seen={}
    for i,row in d.iterrows():
        if pd.notna(row.pmids):
            for pmid in str(row.pmids).split(';'):
                if pmid in seen:u.u(i,seen[pmid])
                else:seen[pmid]=i
        if combine_edit:
            c=int(row.edit_distance_2_cluster)
            if c in cluster_seen:u.u(i,cluster_seen[c])
            else:cluster_seen[c]=i
    roots=[u.f(i) for i in range(len(d))];return pd.factorize(roots)[0]

params=json.load(open(CANON/'results/corrected_sequence/selected_parameters.json'))['cluster']
rows=[];assign=[]
for design,combine in [('pmid_disjoint',False),('pmid_and_edit2_disjoint',True)]:
    g=groups(combine);target=d.label_hiconf.mean();best=None
    for seed in range(4096):
        a,b=next(GroupShuffleSplit(n_splits=1,test_size=.2,random_state=seed).split(d,d.label_hiconf,g))
        score=abs(len(b)/len(d)-.2)+abs(d.iloc[b].label_hiconf.mean()-target)
        cand=(score,seed,a,b)
        if best is None or cand[:2]<best[:2]:best=cand
    score,seed,a,b=best;tr=d.iloc[a].copy();te=d.iloc[b].copy();ytr=tr.label_hiconf.to_numpy(int);yte=te.label_hiconf.to_numpy(int)
    v=TfidfVectorizer(analyzer='char',ngram_range=(2,4),lowercase=False,min_df=2,sublinear_tf=True,dtype=np.float32);xtr=v.fit_transform(tr.peptide);xte=v.transform(te.peptide)
    models={
      'elastic_logistic':SGDClassifier(loss='log_loss',penalty='elasticnet',alpha=params['elastic_logistic']['alpha'],l1_ratio=params['elastic_logistic']['l1_ratio'],class_weight='balanced',max_iter=3000,tol=1e-4,random_state=SEED),
      'random_forest':RandomForestClassifier(n_estimators=500,max_features=params['random_forest']['max_features'],min_samples_leaf=params['random_forest']['min_samples_leaf'],class_weight='balanced_subsample',random_state=SEED,n_jobs=-1),
      'xgboost':xgboost.XGBClassifier(n_estimators=500,max_depth=params['xgboost']['max_depth'],learning_rate=params['xgboost']['learning_rate'],min_child_weight=params['xgboost']['min_child_weight'],subsample=.8,colsample_bytree=.8,reg_lambda=1.,objective='binary:logistic',eval_metric='aucpr',tree_method='hist',scale_pos_weight=float((ytr==0).sum()/(ytr==1).sum()),random_state=SEED,n_jobs=8),
    }
    train_pmids=set().union(*[set(str(x).split(';')) for x in tr.pmids.dropna()]);test_pmids=set().union(*[set(str(x).split(';')) for x in te.pmids.dropna()])
    for name,m in models.items():
        m.fit(xtr,ytr);p=m.predict_proba(xte)[:,1];order=np.argsort(-p,kind='stable')
        rows.append({'design':design,'model':name,'selected_seed':seed,'selection_score':score,'n_groups':len(np.unique(g)),'largest_group':int(pd.Series(g).value_counts().max()),'n_train':len(tr),'n_test':len(te),'test_prevalence':float(yte.mean()),'shared_pmids':len(train_pmids&test_pmids),'edit2_group_overlap':len(set(tr.edit_distance_2_cluster)&set(te.edit_distance_2_cluster)),'roc_auc':roc_auc_score(yte,p),'pr_auc':average_precision_score(yte,p),'normalized_pr_auc':average_precision_score(yte,p)/yte.mean(),'brier_score':brier_score_loss(yte,p),'precision_at_20':float(yte[order[:20]].mean()),'precision_at_100':float(yte[order[:100]].mean())})
    for i in a:assign.append({'design':design,'peptide':d.iloc[i].peptide,'partition':'train','group':int(g[i])})
    for i in b:assign.append({'design':design,'peptide':d.iloc[i].peptide,'partition':'test','group':int(g[i])})
pd.DataFrame(rows).to_csv(OUT/'study_disjoint_sequence_metrics.csv',index=False);pd.DataFrame(assign).to_csv(OUT/'study_disjoint_assignments.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False))
