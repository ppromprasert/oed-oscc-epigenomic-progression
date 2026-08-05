from pathlib import Path
import pandas as pd
R=Path(__file__).resolve().parents[2]; D=R/"data"/"synthetic"
m=pd.read_csv(D/"ref2_metadata.csv"); b=pd.read_csv(D/"ref2_beta_matrix.csv",index_col=0); M=pd.read_csv(D/"ref2_m_matrix.csv",index_col=0)
assert m.sample_id.is_unique
assert list(b.columns)==list(M.columns)
assert set(b.columns)==set(m.sample_id)
assert ((b>=0)&(b<=1)).all().all()
print(f"Validated {len(m)} synthetic Ref2-like samples and {b.shape[0]} synthetic candidate CpGs.")
