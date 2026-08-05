from pathlib import Path
R=Path(__file__).resolve().parents[2]
tokens=[
 "/phi_"+"data_01/","/Users/"+"pp761/","Box-"+"Box/Projects",
 "RECORD_"+"ID",
]
for p in R.rglob("*"):
    if p.name=="privacy_check.py" or not p.is_file(): continue
    if p.suffix.lower() in {".md",".r",".py",".yml",".yaml",".json",".csv",".txt"}:
        t=p.read_text(errors="ignore")
        for x in tokens: assert x not in t, f"Potential restricted token {x} in {p}"
print("Privacy scan passed.")
