import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
import numpy as np
from matplotlib import rcParams
from statannotations.Annotator import Annotator

# Set font and vector output
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42

# Read data
input_file = "step1_Fraction of Mapped Reads.xlsx"
df = pd.read_excel(input_file)

# Data grouping
blood_data = df[df['Style'] == 'blood']['Comparison Rate']
feces_data = df[df['Style'] == 'feces']['Comparison Rate']

# t-test
t_stat, p_value = stats.ttest_ind(blood_data, feces_data)

# Set figure size and style
plt.figure(figsize=(3, 4))
sns.set_style("whitegrid")

# Draw box plot and data points
custom_palette = sns.color_palette(["#244552", "#dc6c50"])
sns.boxplot(x='Style', y='Comparison Rate', data=df, palette=custom_palette, width=0.4, showfliers=False)

# Modify swarmplot - add black edge and fill with sample type color
sns.swarmplot(x='Style', y='Comparison Rate', data=df, 
              palette=custom_palette,  # Use sample type color for fill
              size=6, 
              alpha=0.8,
              edgecolor='black',  # Black edge
              linewidth=1)        # Edge line width

# Add significance annotation
pairs = [("blood", "feces")]
annotator = Annotator(plt.gca(), pairs, data=df, x='Style', y='Comparison Rate')
annotator.configure(test=None, text_format='star', loc='outside', verbose=2)
annotator.set_pvalues([p_value])
annotator.annotate()

# Label median
for style, data in df.groupby('Style'):
    median_val = data['Comparison Rate'].median()
    plt.text(x=0 if style == 'blood' else 1, 
             y=median_val, s=f"{median_val:.3f}", 
             ha='center', va='center', fontsize=10, color='white', weight='bold', 
             bbox=dict(facecolor='#244552' if style == 'blood' else '#dc6c50', alpha=0.8))

# Add sample size information
group_sizes = df.groupby('Style').size()
labels = [f"{style}\n(n={size})" for style, size in group_sizes.items()]
plt.gca().set_xticklabels(labels)

# Set title and labels
#plt.title("Mapping Rate Comparison (Blood vs Feces)", fontsize=14)
plt.ylabel("Mapping Rate", fontsize=12)
plt.xlabel("Sample Type", fontsize=12)

# Save and display
plt.tight_layout()
output_file = "step1_boxplot_updated.pdf"
plt.savefig(output_file, format="pdf", bbox_inches="tight")
plt.show()

print(f"Chart has been successfully saved as PDF file: {output_file}")
print(f"T-test p-value: {p_value:.3e}")