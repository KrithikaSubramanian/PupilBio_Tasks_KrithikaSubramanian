import pandas as pd

file_path = 'PupilBioTest_PMP_revA.csv'  
data = pd.read_csv(file_path,sep=",")
print(data.head())

methylation_columns =['`000', '`001', '`010', '`011', '`100', '`101', '`110', '`111']
data['Total_Reads'] = data[methylation_columns].sum(axis=1)

#Calculate total reads for each tissue
total_reads_by_tissue = data.groupby('Tissue')['Total_Reads'].sum().to_dict()

#Calculate VRF for each PMP
data['VRF'] = data.apply(
    lambda row: row['Total_Reads'] / total_reads_by_tissue[row['Tissue']], axis=1
)

#Filter significant PMPs
significant_pmps_path = 'Significant_PMPs_Full_Data.csv'  # Replace with the path to significant PMPs
significant_pmps = pd.read_csv(significant_pmps_path)
significant_pmps['VRF'] = significant_pmps.apply(
    lambda row: row['Total_Reads'] / total_reads_by_tissue[row['Tissue']], axis=1
)

#Group by PMP and Tissue to calculate mean VRF
mean_vrf_all = data.groupby(['CpG_Coordinates', 'Tissue'])['VRF'].mean().reset_index()
mean_vrf_significant = significant_pmps.groupby(['CpG_Coordinates', 'Tissue'])['VRF'].mean().reset_index()

#Save results
all_output_file = 'All_PMP_Mean_VRF.csv'
significant_output_file = 'Significant_PMP_Mean_VRF.csv'

mean_vrf_all.to_csv(all_output_file, index=False)
mean_vrf_significant.to_csv(significant_output_file, index=False)

print(f"Mean VRF for all PMPs saved to: {all_output_file}")
print(f"Mean VRF for significant PMPs saved to: {significant_output_file}")
