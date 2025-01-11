import pandas as pd
import matplotlib.pyplot as plt
import sys
import seaborn as sns

# Load the dataset
file_path = 'PupilBioTest_PMP_revA.csv'  
data = pd.read_csv(file_path,sep=",")
num_samples = data['Sample_ID'].nunique()
print(f"Number of Samples: {num_samples}")
# Count unique values in the 'Tissue' column (Number of Tissues)
num_tissues = data['Tissue'].nunique()
print(f"Number of Tissues: {num_tissues}")
# Get a list of unique values in the 'Tissue' column (Names of Tissues)
unique_tissues = data['Tissue'].unique()
print(f"Tissues: {unique_tissues}")


# Step 1: Calculate total coverage for each row
coverage_columns = ['`000', '`001', '`010', '`011', '`100', '`101', '`110', '`111']
data['Total_Coverage'] = data[coverage_columns].sum(axis=1)

# Step 2: Group by Tissue and calculate median and CV
tissue_stats = data.groupby('Tissue')['Total_Coverage'].agg(
    Median='median',
    Mean='mean',
    StdDev='std'
).reset_index()

# Add CV to the table
tissue_stats['CV'] = (tissue_stats['StdDev'] / tissue_stats['Mean']) * 100

# Print results
print("Median and CV for each tissue:")
print(tissue_stats)


# Boxplot for Total Coverage by Tissue
plt.figure(figsize=(8, 6))
sns.boxplot(x='Tissue', y='Total_Coverage', data=data)
plt.title('Coverage Distribution by Tissue')
plt.xlabel('Tissue')
plt.ylabel('Total Coverage')
plt.savefig('Coverage_Distribution_Boxplot.png')
plt.close()

# Violin Plot for Total Coverage by Tissue (Optional)
plt.figure(figsize=(8, 6))
sns.violinplot(x='Tissue', y='Total_Coverage', data=data, inner="quartile")
plt.title('Coverage Density by Tissue')
plt.xlabel('Tissue')
plt.ylabel('Total Coverage')
plt.savefig('Coverage_Distribution_ViolinPlot.png')
plt.close()