# First need to merge subset of data we have salivary flow rate + pH for to get subset of metadata

```bash
cd /Users/mann/OneDrive - University of Wyoming/01-Publications/2024-Long_oral_microbiome/07-Revisions-10.7.24
# these are only from visit two and three, these data were not collected for visit one
awk -F"\t" '$13 != 1' metadata-fixed.10.11.24.txt > temp
# get new column in both that has study ID and visit date/sample date
awk -F"\t" '{print $6 "_" $15}' temp > temp2 # opened in excel and merged together as temp_metadata
# same for saliva data
awk -F"\t" '{print $1 "_" $4}' SFR_pH_data_Paul.txt > temp
# should have at minimum 267 matches if this works
awk -F"\t" '{print $7}' SFR_pH_data_Paul.txt | while read line; do grep $line temp_metadata.txt; done > to_merge
```

Now need to hash together to match multiple plaque samples from the same children sampled on the same day

```python
import pandas as pd 
df1 = pd.read_csv("SFR_pH_data_Paul.txt", sep="\t")
df2 = pd.read_csv("temp_metadata.txt", sep="\t")
merged = pd.merge(df1, df2, on="studyID_date")
# write to file
merged.to_csv("merged_saliva_mapping.txt", sep="\t", index=False)
```