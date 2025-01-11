import pandas as pd

# Function to extract allele frequency (AF)
def extract_af(ad_field):
    ad_values = ad_field.split(",")
    if len(ad_values) < 2:
        return 0.0 
    ref, alt = int(ad_values[0]), int(ad_values[1]) 
    total_dp = ref + alt
    af = alt / total_dp if total_dp > 0 else 0.0
    return af

#Calculate statistics 
def calculate_statistics(df, label):
    tumor_af_median = df["Tumor_AF"].median()
    normal_af_median = df["Normal_AF"].median()
    tumor_af_threshold = tumor_af_median + 3 * df["Tumor_AF"].std()
    normal_af_threshold = normal_af_median + 3 * df["Normal_AF"].std()
    tumor_reads_per_million = tumor_af_threshold * 1e6
    normal_reads_per_million = normal_af_threshold * 1e6
    
    return {
        "Category": label,
        "Tumor Median Background AF": tumor_af_median,
        "Normal Median Background AF": normal_af_median,
        "Tumor Confidence Threshold": tumor_af_threshold,
        "Normal Confidence Threshold": normal_af_threshold,
        "Tumor Reads per Million": tumor_reads_per_million,
        "Normal Reads per Million": normal_reads_per_million,
    }

#Load the VCF data into a pandas DataFrame
df = pd.read_csv("Varaint_bcftools_Query.csv", sep=r'\s+', names=[
    "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "DP", "FS", "MBQ", "POPAF", "TLOD",
    "PA220KH-lib09-P19-Tumor_S2_GT", "PA221MH-lib09-P19-Norm_S1_GT",
    "PA220KH-lib09-P19-Tumor_S2_AD", "PA221MH-lib09-P19-Norm_S1_AD",
    "PA220KH-lib09-P19-Tumor_S2_DP", "PA221MH-lib09-P19-Norm_S1_DP",
    "PA220KH-lib09-P19-Tumor_S2_FAD", "PA221MH-lib09-P19-Norm_S1_FAD",
    "PA220KH-lib09-P19-Tumor_S2_SB", "PA221MH-lib09-P19-Norm_S1_SB"
])

#Split CHROM field for additional columns
df[['Gene', 'CHR', 'POS_start', 'POS_end']] = df['CHROM'].str.split('_|:|-', expand=True)
df['GRChr30_POS'] = df['POS_start'].astype('int') + df['POS'].astype('int')

#Calculate Normal_AF and Tumor_AF
df["Normal_AF"] = df["PA221MH-lib09-P19-Norm_S1_AD"].apply(extract_af)
df["Tumor_AF"] = df["PA220KH-lib09-P19-Tumor_S2_AD"].apply(extract_af)

#Classify rows into multiallelic, indel, and snv
multiallelic_df = df[df['ALT'].str.contains(",")]  
remaining_df = df[~df.index.isin(multiallelic_df.index)]  
indel_df = remaining_df[(remaining_df['REF'].str.len() > 1) | (remaining_df['ALT'].str.len() > 1)]  
remaining_df = remaining_df[~remaining_df.index.isin(indel_df.index)]  
snv_df = remaining_df  

#Calculate statistics for each category
stats = []
stats.append(calculate_statistics(df, "Whole Data"))
stats.append(calculate_statistics(snv_df, "SNV"))
stats.append(calculate_statistics(indel_df, "Indel"))
stats.append(calculate_statistics(multiallelic_df, "Multiallelic"))

#Create a statistics table
stats_df = pd.DataFrame(stats)

#Write the data into an Excel file with separate sheets
output_excel = "variant_analysis.xlsx"

with pd.ExcelWriter(output_excel, engine='xlsxwriter') as writer:
    stats_df.to_excel(writer, sheet_name="Statistics", index=False)
    df.to_excel(writer, sheet_name="Whole Data", index=False)
    snv_df.to_excel(writer, sheet_name="SNV", index=False)
    indel_df.to_excel(writer, sheet_name="Indel", index=False)
    multiallelic_df.to_excel(writer, sheet_name="Multiallelic", index=False)

print(f"Results written to {output_excel}")
