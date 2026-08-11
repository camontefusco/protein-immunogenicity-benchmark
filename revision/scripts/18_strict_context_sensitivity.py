from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import FunctionTransformer, OneHotEncoder, StandardScaler

CANON=Path(__file__).resolve().parents[1]
OUT=CANON/'results/forensic_sensitivities';OUT.mkdir(parents=True,exist_ok=True)
SEED=42;N_BOOT=2000
BIO=['virus_species','virus_strain','antigen','protein'];ASSAY=['assay_method','readout','host'];SAFE=['peptide',*BIO,*ASSAY]

def pipe(with_context):
    ts=[('sequence',TfidfVectorizer(analyzer='char',ngram_range=(2,4),lowercase=False,min_df=2,sublinear_tf=True,dtype=np.float32),'peptide')]
    if with_context:
        ts += [('categorical',OneHotEncoder(handle_unknown='ignore',min_frequency=2),BIO+ASSAY),
               ('numeric',Pipeline([('log1p',FunctionTransformer(np.log1p,feature_names_out='one-to-one')),('scale',StandardScaler())]),['length'])]
    return Pipeline([('features',ColumnTransformer(ts,remainder='drop')),('classifier',LogisticRegression(C=1.,class_weight='balanced',max_iter=3000,random_state=SEED,solver='liblinear'))])

raw=pd.read_csv(CANON/'data/curated/assay_level_with_year.csv')
raw.peptide=raw.peptide.astype(str).str.upper().str.strip();alphabet=set('ACDEFGHIKLMNPQRSTVWY');raw=raw[raw.peptide.map(lambda x:set(x)<=alphabet)].copy()
for c in BIO+ASSAY:raw[c]=raw[c].fillna('__MISSING__').astype(str)
raw=raw.drop_duplicates()
latest=raw.groupby(SAFE,as_index=False,dropna=False).agg(latest_pub_year_reconstructed=('latest_pub_year','max'))
ctx=pd.read_csv(CANON/'data/curated/context_dataset_and_splits.csv').merge(latest,on=SAFE,how='left',validate='one_to_one')
constructions={
 'majority_all':np.ones(len(ctx),bool),
 'majority_excluding_ties':~ctx.majority_tie.astype(bool).to_numpy(),
 'consistent_labels_only':ctx.n_unique_labels.eq(1).to_numpy(),
}
rows=[];preds=[]
for cname,mask in constructions.items():
    sub=ctx.loc[mask].copy();train=sub[sub.temporal_split.eq('train')&sub.latest_pub_year_reconstructed.lt(2021)].copy();test=sub[sub.temporal_split.eq('test')].copy()
    for name,withctx in [('sequence_only',False),('sequence_plus_all_no_evidence',True)]:
        cols=['peptide']+(BIO+ASSAY+['length'] if withctx else []);m=pipe(withctx);m.fit(train[cols],train.y);p=m.predict_proba(test[cols])[:,1]
        rows.append({'construction':cname,'model':name,'n_train':len(train),'n_test':len(test),'prevalence':test.y.mean(),'roc_auc':roc_auc_score(test.y,p),'pr_auc':average_precision_score(test.y,p),'brier_score':brier_score_loss(test.y,p)})
        q=test[['peptide','y']].copy();q['row_id']=q.index;q['construction']=cname;q['model']=name;q['probability']=p;preds.append(q)
pred=pd.concat(preds,ignore_index=True);pred.to_csv(OUT/'context_strict_temporal_predictions.csv',index=False)
diffs=[]
for cname,frame in pred.groupby('construction'):
    wide=frame.pivot(index='row_id',columns='model',values='probability').sort_index();meta=frame.drop_duplicates('row_id').set_index('row_id').loc[wide.index,['peptide','y']]
    y=meta.y.to_numpy(int);units=meta.peptide.astype(str).to_numpy();uniq=np.unique(units);idx={u:np.flatnonzero(units==u) for u in uniq};rng=np.random.default_rng(SEED);vals=[]
    while len(vals)<N_BOOT:
        take=np.concatenate([idx[u] for u in rng.choice(uniq,size=len(uniq),replace=True)])
        if np.unique(y[take]).size<2:continue
        vals.append(average_precision_score(y[take],wide.iloc[take]['sequence_plus_all_no_evidence'])-average_precision_score(y[take],wide.iloc[take]['sequence_only']))
    point=average_precision_score(y,wide.sequence_plus_all_no_evidence)-average_precision_score(y,wide.sequence_only)
    diffs.append({'construction':cname,'difference_context_minus_sequence':point,'ci_low':np.quantile(vals,.025),'ci_high':np.quantile(vals,.975),'replicates':N_BOOT})
pd.DataFrame(rows).to_csv(OUT/'context_strict_temporal_construction_metrics.csv',index=False)
pd.DataFrame(diffs).to_csv(OUT/'context_strict_temporal_paired_uplift.csv',index=False)
with open(OUT/'context_strict_temporal_provenance.json','w') as f:json.dump({'training_rule':'temporal_split=train and reconstructed latest publication year < 2021','test_rule':'temporal_split=test (earliest publication year >= 2021)','models':'fixed balanced logistic C=1.0','bootstrap':'2000 peptide-grouped paired replicates','seed':SEED},f,indent=2)
print(pd.DataFrame(rows).to_string(index=False));print('\n',pd.DataFrame(diffs).to_string(index=False))
