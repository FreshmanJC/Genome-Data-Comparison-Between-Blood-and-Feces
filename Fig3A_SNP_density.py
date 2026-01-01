import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Set font and output format as vector
rcParams['pdf.fonttype'] = 42  # Use Type 3 font, better compatibility with AI software
rcParams['ps.fonttype'] = 42

# Read Excel file
df = pd.read_excel("step6.blood.xlsx")

# Calculate SNP density
window_size_kb = 100  # 100 kb
df['SNP_Density'] = df['SNP_Count'] / window_size_kb  # SNP density per kb

# Assign a continuous starting position index for each chromosome, displaying chromosome numbers on X-axis
chromosome_order = [f'Chrwds{i}' for i in range(1, 21)]  # Chrwds1 to Chrwds20
df['Chromosome'] = pd.Categorical(df['Chromosome'], categories=chromosome_order, ordered=True)

# Add new column: position index for each chromosome segment
chromosome_offset = {}
offset = 0

for chrom in chromosome_order:
    chrom_length = df[df['Chromosome'] == chrom]['Start'].max()  # Maximum starting position for each chromosome
    chromosome_offset[chrom] = offset
    offset += chrom_length + 1e6  # 1 Mb interval

df['Relative_Position'] = df.apply(lambda row: row['Start'] + chromosome_offset[row['Chromosome']], axis=1)

# Data smoothing: rolling average, window size adjusted to 50
df['SNP_Density_Smoothed'] = df.groupby('Chromosome')['SNP_Density'].transform(lambda x: x.rolling(50, min_periods=1).mean())

# Data sampling: take one point every 5 rows
df = df.iloc[::5, :]

# Custom color scheme (single population color)
custom_color = '#e66f51'

# Plotting: use `Relative_Position` instead of `Start`
plt.figure(figsize=(14, 4))  # Increase figure size
sns.lineplot(
    data=df, 
    x='Relative_Position', 
    y='SNP_Density_Smoothed', 
    color=custom_color, 
    ci=None, 
    linewidth=1,  # Reduce line width
    alpha=1       # Lower transparency
)

# Set chromosome labels as X-axis ticks
chrom_labels = {v: k for k, v in chromosome_offset.items()}
plt.xticks(list(chrom_labels.keys()), list(chrom_labels.values()), rotation=45)
plt.xlabel('Chromosome', fontsize=14)
plt.ylabel('SNP Density (per kb)', fontsize=14)
plt.title('SNP Density Distribution by Chromosome (Blood)', fontsize=16)
plt.tight_layout()

# Save as vector PDF format
plt.savefig("step6_snp_density_blood_adjusted.pdf", format='pdf', bbox_inches='tight')
plt.show()