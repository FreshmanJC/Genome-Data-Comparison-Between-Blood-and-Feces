import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rcParams

# Set font and output format as vector
rcParams['pdf.fonttype'] = 42  # Use Type 3 font, better compatibility with AI software
rcParams['ps.fonttype'] = 42

# Read Excel file
file_path = "step7.pi.xlsx"
df = pd.read_excel(file_path, sheet_name="Sheet1")

# Data cleaning (ensure no non-numeric content)
df = df.apply(pd.to_numeric, errors='coerce').dropna()

# Set plot style
plt.figure(figsize=(8, 7))
sns.set(style="whitegrid", font_scale=1.2)

# Convert data to long format (suitable for seaborn plotting)
df_long = df.melt(var_name="Sample Type", value_name="PI Value")

# Draw box plot
boxplot = sns.boxplot(
    x="Sample Type",
    y="PI Value",
    data=df_long,
    palette=["#e66f51", "#299d92"],
    width=0.5,
    showfliers=True,  # Ensure outliers are displayed
    flierprops={
        'marker': 'o',          # Set as circle
        'markerfacecolor': 'black',  # Fill color
        'markeredgecolor': 'black',  # Edge color
        'markersize': 5          # Size
    }
)

# Add title and labels
plt.ylabel("Nucelotide diversity across windows", labelpad=10)

# Optimize axis display
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# Adjust layout and display chart
plt.tight_layout()
plt.grid(axis='y', linestyle='--', alpha=0)  # Add grid lines for easier reading

# Save as PDF file (vectorized text)
output_file = "step7.pi.pdf"
plt.savefig(output_file, format="pdf", bbox_inches="tight")
plt.show()