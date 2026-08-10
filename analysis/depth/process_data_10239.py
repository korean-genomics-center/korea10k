# %%
from collections import Counter

import pandas as pd

# %%# %%
df = pd.read_excel("/BiO/Access/kyungwhan1998/genome/depthCoverage/10243sample_list.xlsx")
dict_id_cohort = dict(zip(df["ID"], df["모집군"]))
Counter(dict_id_cohort.values())

# %%
list_samples_exclude = ["KU10K-10433", "KU10K-10689", "KU10K-04846", "KU10K-10007"]

dict_id_cohort_filt = dict()
for sample, cohort in dict_id_cohort.items():
    if sample not in list_samples_exclude:
        dict_id_cohort_filt[sample] = cohort
    
sum(Counter(dict_id_cohort_filt.values()).values())