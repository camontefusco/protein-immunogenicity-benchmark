from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score,brier_score_loss,roc_auc_score
from sklearn.model_selection import GroupShuffleSplit
from sklearn.pipeline import Pipeline

CANON=Path(__file__).resolve().parents[1];OUT=CANON/'results/forensic_sensitivities';OUT.mkdir(parents=True,exist_ok=True);SEED=42

def model():return Pipeline([('tfidf',TfidfVectorizer(analyzer='char',ngram_range=(2,4),lowercase=False,min_df=2,sublinear_tf=True,dtype=np.float32)),('classifier',RandomForestClassifier(n_estimators=500,min_samples_leaf=2,max_features='sqrt',class_weight='balanced_subsample',random_state=SEED,n_jobs=-1))])
def ece(y,p,kind):
    if kind=='equal_width':edges=np.linspace(0,1,11);membership=np.clip(np.digitize(p,edges[1:-1],right=True),0,9)
    else:
        edges=np.quantile(p,np.linspace(0,1,11));membership=np.clip(np.digitize(p,edges[1:-1],right=True),0,9)
    return sum((membership==i).mean()*abs(y[membership==i].mean()-p[membership==i].mean()) for i in range(10) if (membership==i).any())
d=pd.read_csv(CANON/'data/splits/peptide_split_assignments_v1.csv');rows=[];preds=[]
for split in ['exact','cluster','temporal']:
    dev=d[d[f'{split}_split'].eq('train')].copy();test=d[d[f'{split}_split'].eq('test')].copy();cands=[]
    for cid,(a,b) in enumerate(GroupShuffleSplit(n_splits=256,test_size=.2,random_state=SEED).split(dev,dev.label_hiconf,dev.edit_distance_2_cluster)):
        cal=dev.iloc[b];score=abs(len(cal)/len(dev)-.2)+abs(cal.label_hiconf.mean()-dev.label_hiconf.mean());cands.append((score,cid,a,b))
    _,cid,a,b=min(cands,key=lambda x:x[:2]);fit=dev.iloc[a];cal=dev.iloc[b];m=model().fit(fit.peptide,fit.label_hiconf);pc=m.predict_proba(cal.peptide)[:,1];pt=m.predict_proba(test.peptide)[:,1];eps=1e-6
    logc=np.log(np.clip(pc,eps,1-eps)/np.clip(1-pc,eps,1-eps));logt=np.log(np.clip(pt,eps,1-eps)/np.clip(1-pt,eps,1-eps));pl=LogisticRegression(C=1e6,max_iter=2000,random_state=SEED).fit(logc.reshape(-1,1),cal.label_hiconf);pp=pl.predict_proba(logt.reshape(-1,1))[:,1]
    y=test.label_hiconf.to_numpy(int);null=float(y.mean()*(1-y.mean()))
    for status,p in [('uncalibrated',pt),('platt',pp)]:
        bs=brier_score_loss(y,p);rows.append({'split':split,'status':status,'n_fit':len(fit),'n_calibration':len(cal),'n_test':len(test),'prevalence':y.mean(),'roc_auc':roc_auc_score(y,p),'pr_auc':average_precision_score(y,p),'brier_score':bs,'null_brier_score':null,'brier_skill_score':1-bs/null,'ece_equal_width_10':ece(y,p,'equal_width'),'ece_equal_frequency_10':ece(y,p,'equal_frequency'),'mean_probability':p.mean()})
        q=test[['peptide','label_hiconf']].copy();q['split']=split;q['status']=status;q['probability']=p;preds.append(q)
pd.DataFrame(rows).to_csv(OUT/'calibration_extended_metrics.csv',index=False);pd.concat(preds).to_csv(OUT/'calibration_predictions_recomputed.csv',index=False);print(pd.DataFrame(rows).to_string(index=False))
