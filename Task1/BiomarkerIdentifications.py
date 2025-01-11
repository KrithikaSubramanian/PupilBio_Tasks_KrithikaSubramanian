import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact
import matplotlib.pyplot as plt
import sys

file_path = 'PupilBioTest_PMP_revA.csv'  
data = pd.read_csv(file_path,sep=",")
methylation_columns =['`000', '`001', '`010', '`011', '`100', '`101', '`110', '`111']
data['Total_Reads'] = data[methylation_columns].sum(axis=1)

#Froup by CpG_Coordinates and Tissue to calculate total reads
grouped_data = data.groupby(['CpG_Coordinates', 'Tissue']).agg(
    Total_Reads=('Total_Reads', 'sum')
).reset_index()

#Pivot the data to create a contingency table for each PMP
pivot_table = grouped_data.pivot(index='CpG_Coordinates', columns='Tissue', values='Total_Reads').fillna(0)

#Perform Chi-Square for each PMP
results = []
for pmp in pivot_table.index:
    observed = pivot_table.loc[pmp].values  # Observed counts for cfDNA and Islet
    if sum(observed) > 0:  # Only test if there are reads for the PMP
        try:
            # Use Chi-Square Test
            chi2, p, _, _ = chi2_contingency([observed, observed])
            results.append({'CpG_Coordinates': pmp, 'p-value': p})
        except ValueError:  # If Chi-Square fails due to low counts, use Fisher's Exact Test
            odds_ratio, p = fisher_exact([[observed[0], observed[1]], [observed[1], observed[0]]])
            results.append({'CpG_Coordinates': pmp, 'p-value': p})

#Convert results to a DataFrame
pmp_pvalues = pd.DataFrame(results)

# Merge significant PMPs with the original dataset and save
significant_pmps = pmp_pvalues[pmp_pvalues['p-value'] < 0.05] 
significant_full_data = data.merge(significant_pmps, on='CpG_Coordinates', how='inner')  
significant_output_file_path = 'Significant_PMPs_Data.csv'
significant_full_data.to_csv(significant_output_file_path, index=False)
print(f"Significant PMPs saved to: {significant_output_file_path}")

# Summarize significant PMPs by tissue
summary = significant_full_data.groupby('Tissue').agg(
    Total_Significant_PMPs=('CpG_Coordinates', 'nunique'),
    Average_Reads=('Total_Reads', 'mean')
).reset_index()
print(summary)

significant_full_data['Total_Reads'].hist(bins=20, color='blue', alpha=0.7)
plt.title('Distribution of Total Reads for Significant PMPs')
plt.xlabel('Total Reads')
plt.ylabel('Frequency')
plt.savefig('PMP_distributions.png')
plt.close()
