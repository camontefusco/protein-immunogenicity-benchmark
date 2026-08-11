from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.linear_model import LogisticRegression, SGDClassifier
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import FunctionTransformer, OneHotEncoder, StandardScaler

CANON=Path(__file__).resolve().parents[1]
OUT=CANON/'results/forensic_sensitivities';OUT.mkdir(parents=True,exist_ok=True)
SEED=42

SEQ_PARAMS=json.load(open(CANON/'results/corrected_sequence/selected_parameters.json'))

def vectorizer():
    return TfidfVectorizer(analyzer='char',ngram_range=(2,4),lowercase=False,min_df=2,sublinear_tf=True,dtype=np.float32)

def estimator(name,params,y):
    weight=float((y==0).sum()/(y==1).sum())
    if name=='elastic_logistic': return SGDClassifier(loss='log_loss',penalty='elasticnet',alpha=params['alpha'],l1_ratio=params['l1_ratio'],class_weight='balanced',max_iter=3000,tol=1e-4,random_state=SEED)
    if name=='random_forest': return RandomForestClassifier(n_estimators=500,max_features=params['max_features'],min_samples_leaf=params['min_samples_leaf'],class_weight='balanced_subsample',random_state=SEED,n_jobs=-1)
    if name=='xgboost': return xgboost.XGBClassifier(n_estimators=500,max_depth=params['max_depth'],learning_rate=params['learning_rate'],min_child_weight=params['min_child_weight'],subsample=.8,colsample_bytree=.8,reg_lambda=1.,objective='binary:logistic',eval_metric='aucpr',tree_method='hist',scale_pos_weight=weight,random_state=SEED,n_jobs=8)
    raise KeyError(name)

def metrics(y,p):
    order=np.argsort(-p,kind='stable');r={'roc_auc':roc_auc_score(y,p),'pr_auc':average_precision_score(y,p),'brier_score':brier_score_loss(y,p)}
    for k in [20,50,100,200,500]:
        kk=min(k,len(y));prec=float(y[order[:kk]].mean());r[f'precision_at_{k}']=prec;r[f'enrichment_at_{k}']=prec/float(y.mean())
    return r

def run_sequence(design,train,test,label):
    ytr=train.label_hiconf.to_numpy(int);yte=test.label_hiconf.to_numpy(int)
    probs={}
    inner=StratifiedGroupKFold(n_splits=3,shuffle=True,random_state=SEED)
    folds=list(inner.split(train.peptide,ytr,train.edit_distance_2_cluster))
    oof=np.zeros((len(train),3))
    for j,name in enumerate(['elastic_logistic','random_forest','xgboost']):
        params=SEQ_PARAMS[design][name]
        for a,b in folds:
            v=vectorizer();xa=v.fit_transform(train.iloc[a].peptide);xb=v.transform(train.iloc[b].peptide)
            m=estimator(name,params,ytr[a]);m.fit(xa,ytr[a]);oof[b,j]=m.predict_proba(xb)[:,1]
        v=vectorizer();xtr=v.fit_transform(train.peptide);xte=v.transform(test.peptide)
        m=estimator(name,params,ytr);m.fit(xtr,ytr);probs[name]=m.predict_proba(xte)[:,1]
    meta=LogisticRegression(C=1.,class_weight='balanced',max_iter=2000,random_state=SEED).fit(oof,ytr)
    probs['stacked_ensemble']=meta.predict_proba(np.column_stack([probs[x] for x in ['elastic_logistic','random_forest','xgboost']]))[:,1]
    rows=[]
    for name,p in probs.items(): rows.append({'analysis':label,'design':design,'model':name,'n_train':len(train),'n_test':len(test),'test_prevalence':float(yte.mean()),**metrics(yte,p)})
    return rows

seq=pd.read_csv(CANON/'data/splits/peptide_split_assignments_v1.csv')
seq=seq[seq.valid_sequence.astype(bool)].copy()
rows=[]
# Strict temporal: exclude any training peptide whose aggregated evidence reaches 2021 or later.
tr=seq[seq.temporal_split.eq('train') & seq.latest_pub_year.lt(2021)].copy();te=seq[seq.temporal_split.eq('test')].copy()
rows+=run_sequence('temporal',tr,te,'strict_temporal_latest_before_2021')
# Evidence sensitivity: retain at least 10 tested participants/units in train and test.
for design in ['exact','cluster','temporal']:
    tr=seq[seq[f'{design}_split'].eq('train') & seq.total_tested.ge(10)].copy();te=seq[seq[f'{design}_split'].eq('test') & seq.total_tested.ge(10)].copy()
    rows+=run_sequence(design,tr,te,'minimum_total_tested_10')
pd.DataFrame(rows).to_csv(OUT/'sequence_sensitivity_metrics.csv',index=False)

# Reconstruct latest year per context from the recovered assay-level table.
raw=pd.read_csv(CANON/'data/curated/assay_level_with_year.csv')
bio=['virus_species','virus_strain','antigen','protein'];assay=['assay_method','readout','host'];safe=['peptide',*bio,*assay]
raw.peptide=raw.peptide.astype(str).str.upper().str.strip();alphabet=set('ACDEFGHIKLMNPQRSTVWY');raw=raw[raw.peptide.map(lambda x:set(x)<=alphabet)].copy()
for c in bio+assay:raw[c]=raw[c].fillna('__MISSING__').astype(str)
raw=raw.drop_duplicates()
latest=raw.groupby(safe,as_index=False,dropna=False).agg(latest_pub_year_reconstructed=('latest_pub_year','max'))
context=pd.read_csv(CANON/'data/curated/context_dataset_and_splits.csv').merge(latest,on=safe,how='left',validate='one_to_one')

def context_pipe(seq,cats,nums):
    ts=[]
    if seq:ts.append(('sequence',vectorizer(),'peptide'))
    if cats:ts.append(('categorical',OneHotEncoder(handle_unknown='ignore',min_frequency=2),cats))
    if nums:ts.append(('numeric',Pipeline([('log1p',FunctionTransformer(np.log1p,feature_names_out='one-to-one')),('scale',StandardScaler())]),nums))
    return Pipeline([('features',ColumnTransformer(ts,remainder='drop')),('classifier',LogisticRegression(C=1.,class_weight='balanced',max_iter=3000,random_state=SEED,solver='liblinear'))])

specs={
 'sequence_only':(True,[],[]),'biological_only':(False,bio,['length']),'assay_only':(False,assay,[]),
 'sequence_plus_biological':(True,bio,['length']),'sequence_plus_assay':(True,assay,[]),
 'sequence_plus_all_no_evidence':(True,bio+assay,['length']),'all_context_no_sequence_no_evidence':(False,bio+assay,['length'])}
tr=context[context.temporal_split.eq('train') & context.latest_pub_year_reconstructed.lt(2021)].copy();te=context[context.temporal_split.eq('test')].copy()
crows=[]
for name,(withseq,cats,nums) in specs.items():
    cols=(['peptide'] if withseq else [])+cats+nums;m=context_pipe(withseq,cats,nums);m.fit(tr[cols],tr.y);p=m.predict_proba(te[cols])[:,1]
    crows.append({'analysis':'strict_temporal_latest_before_2021','model':name,'n_train':len(tr),'n_test':len(te),'test_prevalence':float(te.y.mean()),**metrics(te.y.to_numpy(int),p)})
pd.DataFrame(crows).to_csv(OUT/'context_strict_temporal_metrics.csv',index=False)
print('\nSEQUENCE SENSITIVITIES')
print(pd.DataFrame(rows)[['analysis','design','model','n_train','n_test','test_prevalence','pr_auc','roc_auc','precision_at_20']].to_string(index=False))
print('\nCONTEXT STRICT TEMPORAL')
print(pd.DataFrame(crows)[['model','n_train','n_test','test_prevalence','pr_auc','roc_auc','brier_score']].to_string(index=False))
