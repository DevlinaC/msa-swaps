import numpy as np, pandas as pd, csv
from scipy.stats import ks_2samp

AF2_SW='/mnt/data/AF2_TMs.csv'; AF2_SINGLE='/mnt/data/AF2_single_TMs.csv'
AF3_SW='/mnt/data/tm_align_full_results.csv'; AF3_SINGLE='/mnt/data/singleseq_af3_tmscores.csv'

def read_col(p,c):
    return pd.read_csv(p)[c].dropna().to_numpy(float)

full3=pd.read_csv(AF3_SW)
junk3=full3[full3['Job_Type'].str.lower().eq('junk')].copy()
junk3['query']=junk3['Job_Label'].str.split('_').str[0]
af3q=junk3[junk3['Compared_Protein'].astype(str)==junk3['query'].astype(str)].copy()
s3=pd.read_csv(AF3_SINGLE)
queries=s3['Target'].drop_duplicates().tolist()

# AF3 dicts
A3={q:af3q.loc[af3q['query']==q,'TM_Score'].to_numpy(float) for q in queries}
B3={q:s3.loc[s3['Target']==q,'TM_Score'].to_numpy(float) for q in queries}

# AF2 reconstruct query labels
sw2=read_col(AF2_SW,'AF2_A'); si2=read_col(AF2_SINGLE,'AF2_A')
job_labels=junk3['Job_Label'].drop_duplicates().tolist()
A2={q:[] for q in queries}
for i,label in enumerate(job_labels):
    q=label.split('_')[0]
    A2[q].extend(sw2[i*25:(i+1)*25].tolist())
A2={q:np.asarray(v,float) for q,v in A2.items()}
B2={q:si2[i*25:(i+1)*25] for i,q in enumerate(queries)}

def weighted_cluster_ks(A,B,nperm=20000,seed=20261029):
    qs=list(queries)
    pooled=np.unique(np.concatenate([np.concatenate([A[q] for q in qs]), np.concatenate([B[q] for q in qs])]))
    # Difference in query-specific ECDFs at each pooled point.
    Dq=np.empty((len(qs), len(pooled)),dtype=np.float32)
    for i,q in enumerate(qs):
        a=np.sort(A[q]); b=np.sort(B[q])
        Fa=np.searchsorted(a,pooled,side='right')/len(a)
        Fb=np.searchsorted(b,pooled,side='right')/len(b)
        Dq[i]=(Fa-Fb).astype(np.float32)
    obs=float(np.max(np.abs(Dq.mean(axis=0))))
    rng=np.random.default_rng(seed)
    ge=0
    batch=250
    done=0
    permDs=[]
    while done<nperm:
        m=min(batch,nperm-done)
        signs=rng.choice(np.array([-1.0,1.0],dtype=np.float32),size=(m,len(qs)))
        # paired label swap means +/- query-specific difference
        diffs=(signs @ Dq)/len(qs)
        vals=np.max(np.abs(diffs),axis=1)
        ge += int(np.sum(vals >= obs-1e-12))
        permDs.extend(vals.tolist())
        done += m
    p=(ge+1)/(nperm+1)
    return obs,p,np.quantile(permDs,[0.5,0.95,0.99]),len(pooled)

for name,A,B in [('AF2',A2,B2),('AF3',A3,B3)]:
    rawA=np.concatenate([A[q] for q in queries]); rawB=np.concatenate([B[q] for q in queries])
    naive=ks_2samp(rawA,rawB,alternative='two-sided',method='asymp')
    obs,p,qsupp,nuniq=weighted_cluster_ks(A,B)
    print(name)
    print('raw n',len(rawA),len(rawB),'naive KS D',naive.statistic,'p',naive.pvalue)
    print('query-weighted KS D',obs,'cluster paired permutation p',p,'perm quantiles',qsupp,'support',nuniq)
    # Also average per-query KS D descriptive
    per=[ks_2samp(A[q],B[q],method='asymp').statistic for q in queries]
    print('per-query KS D median',np.median(per),'mean',np.mean(per),'range',min(per),max(per))
