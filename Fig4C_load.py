import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind

# Read data
file_path = 'step10_3.load.xlsx'  # Modify to your file path
df = pd.read_excel(file_path)

# Ensure 'group' and 'Style' columns are categorical data
df['group'] = df['group'].astype('category')
df['Style'] = df['Style'].astype('category')

# Set figure size
plt.figure(figsize=(6, 8))

# Create violin plot, x-axis as group, grouped by Style, modify colors
palette = {'blood': '#c9664c', 'feces': '#2d9187'}

# Create violin plot, split=False to not separate display of sample types
sns.violinplot(x='group', y='2*Hom/(2*Hom+1*Het)', hue='Style', data=df, split=False, inner="box", dodge=True, palette=palette, width=0.4)

# Perform t-test and annotate p-values
gs_blood = df[(df['group'] == 'GS') & (df['Style'] == 'blood')]['2*Hom/(2*Hom+1*Het)']
gs_feces = df[(df['group'] == 'GS') & (df['Style'] == 'feces')]['2*Hom/(2*Hom+1*Het)']
lof_blood = df[(df['group'] == 'LOF') & (df['Style'] == 'blood')]['2*Hom/(2*Hom+1*Het)']
lof_feces = df[(df['group'] == 'LOF') & (df['Style'] == 'feces')]['2*Hom/(2*Hom+1*Het)']

# Calculate p-values
p_value_gs = ttest_ind(gs_blood, gs_feces)[1]  # p-value for GS group
p_value_lof = ttest_ind(lof_blood, lof_feces)[1]  # p-value for LOF group

# Add p-values
plt.text(0, 0.78, f'{p_value_gs:.2E}', horizontalalignment='center', verticalalignment='center', fontsize=12, color='black')
plt.text(1, 0.78, f'{p_value_lof:.2E}', horizontalalignment='center', verticalalignment='center', fontsize=12, color='black')

# Set title and labels
plt.xlabel("Group", fontsize=14)
plt.ylabel("2*Hom/(2*Hom+1*Het)", fontsize=14)

# Set y-axis range
plt.ylim(0.1, 0.8)  # Set y-axis range, from 0.1 to 0.7

# Set outer border color of box plot as hollow
for ax in plt.gcf().get_axes():
    for artist in ax.artists:
        artist.set_edgecolor('black')  # Set outer border color
        artist.set_facecolor('none')  # Remove fill color

# Save figure as PDF
plt.tight_layout()
plt.savefig('step10_3.load.pdf', format='pdf')

# Display figure
plt.show()