from pathlib import Path
import pandas as pd
R=Path(__file__).resolve().parents[2]; D=R/"data"/"synthetic"
ref=pd.read_csv(D/"ref2_metadata.csv")
b=pd.read_csv(D/"ref2_candidate_beta.csv",index_col=0)
M=pd.read_csv(D/"ref2_candidate_mvalues.csv",index_col=0)
k=pd.read_csv(D/"known_oscc_annotation.csv")
c=pd.read_csv(D/"cosmx_annotation.csv")
assert len(ref)==28
assert (ref.Progression=="Progressor").sum()==25
assert (ref.Progression=="NonProgressor").sum()==3
assert len(k)==342 and len(c)==355
assert list(b.columns)==list(M.columns)
assert set(b.columns)==set(ref.matrix_sample_id)
assert ((b>=0)&(b<=1)).all().all()
print("Validated expanded synthetic OED/OSCC layer: Ref2 25/3, 342 known CpGs, 355 CosMx CpGs.")
