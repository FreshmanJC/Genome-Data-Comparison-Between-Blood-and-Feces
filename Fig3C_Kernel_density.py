import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams

# File paths
faece_file_path = "faece_allele_freq.frq"  # Fecal data file path
blood_file_path = "blood_allele_freq.frq"  # Blood data file path

# Ensure output as vectorized font for easy font modification in AI software
rcParams['pdf.fonttype'] = 42  # Embed vector font
rcParams['ps.fonttype'] = 42  # Embed vector font

# Function: Extract minor allele frequency (< 0.5)
def extract_minor_allele_freq(file_path):
    df = pd.read_csv(file_path, sep="\t")  # Read .frq file
    frequencies = []  # Store minor allele frequencies
    for freq_str in df['{ALLELE:FREQ}']:
        alleles = freq_str.split()  # Split string, e.g., "A:0.8 T:0.2"
        freqs = [float(a.split(":")[1]) for a in alleles]  # Extract frequency values
        minor_freq = min(freqs)  # Take the smaller frequency value, minor allele frequency
        if minor_freq < 0.5:  # Keep minor allele frequency values < 0.5
            frequencies.append(minor_freq)
    return frequencies

# Extract minor allele frequencies for fecal and blood samples
faece_frequencies = extract_minor_allele_freq(faece_file_path)
blood_frequencies = extract_minor_allele_freq(blood_file_path)

# Create DataFrame for plotting
combined_df = pd.DataFrame({
    'Frequency': faece_frequencies + blood_frequencies,
    'Sample Type': ['Faeces'] * len(faece_frequencies) + ['Blood'] * len(blood_frequencies)
})

# Draw stacked density curve plot
plt.figure(figsize=(4, 4))

# Set custom color scheme
custom_palette = {'Faeces': '#7e9d31', 'Blood': '#ce8862'}  # Can adjust color codes as needed

sns.kdeplot(data=combined_df, x='Frequency', hue='Sample Type', fill=True, alpha=0.6, bw_adjust=1.2, palette=custom_palette)
plt.xlabel('Minor Allele Frequency', fontsize=14)
plt.ylabel('Density', fontsize=14)
plt.title('Density Plot of Minor Allele Frequency (Faeces vs Blood)', fontsize=16)

# Remove legend outer border, optimize spacing
plt.legend(title='Sample Type', frameon=False, loc='upper right')
plt.tight_layout()  # Automatically adjust spacing

# Save chart as PDF format
plt.savefig("allele_frequency_density_plot_customized.pdf", format='pdf', bbox_inches='tight')
plt.show()