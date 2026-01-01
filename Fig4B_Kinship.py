import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Set font and output format as vector
rcParams['pdf.fonttype'] = 42  # Use Type 3 font, better compatibility with AI software
rcParams['ps.fonttype'] = 42

# Read Excel file
df = pd.read_excel("step8.feces_relationship.xlsx", sheet_name="Sheet1")

# Create symmetric matrix (handle asymmetric relationships)
all_individuals = sorted(set(df['Ind1']).union(set(df['Ind2'])))
matrix = pd.DataFrame(0.0, index=all_individuals, columns=all_individuals)

# Fill matrix data
for _, row in df.iterrows():
    i, j, val = row['Ind1'], row['Ind2'], row['Relationship']
    matrix.at[i, j] = val
    matrix.at[j, i] = val  # Assume relationship is symmetric

# Draw heatmap
plt.figure(figsize=(8, 6))
sns.heatmap(
    matrix,
    annot=False,
    cmap="YlGnBu",
    linewidths=0.5,
    square=True,
    cbar_kws={'label': 'Kinship Coefficient'}
)

# Adjust graphic format
plt.title("Kinship Relationship Heatmap", fontsize=14, pad=20)
plt.xticks(rotation=45, ha='right', fontsize=8)
plt.yticks(rotation=0, fontsize=8)
plt.xlabel("Individual 2", fontsize=10)
plt.ylabel("Individual 1", fontsize=10)
plt.tight_layout()

# Save image
# Save as PDF file (vectorized text)
output_file = "step8_2feces.pdf"
plt.savefig(output_file, format="pdf", bbox_inches="tight")
plt.show()

print(f"Chart has been successfully saved as PDF file: {output_file}")