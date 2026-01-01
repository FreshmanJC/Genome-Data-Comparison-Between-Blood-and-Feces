import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
from matplotlib import rcParams

# Set font and graphic export format
rcParams['pdf.fonttype'] = 42  # Type 3 font (compatible with AI software)
rcParams['ps.fonttype'] = 42

# Read data
input_file = "step2.coverage_tendencies.xlsx"
df = pd.read_excel(input_file)

# Set plot style
sns.set(style="whitegrid")
colors = {"blood": "#244552", "feces": "#dc6c50"}
coverage_depths = ["1X Coverage", "5X Coverage", "10X Coverage"]

# Create horizontally arranged subplots
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(6, 4), sharey=True)

for i, depth in enumerate(coverage_depths):
    ax = axes[i]
    
    # Prepare paired data
    blood_df = df[df['Style'] == 'blood'][["ID", depth]].set_index("ID")
    feces_df = df[df['Style'] == 'feces'][["ID", depth]].set_index("ID")
    merged = pd.merge(blood_df, feces_df, left_index=True, right_index=True, suffixes=('_blood', '_feces')).reset_index()
    
    # Long format for plotting
    plot_df = pd.melt(merged, id_vars="ID", value_vars=[f"{depth}_blood", f"{depth}_feces"],
                      var_name="Style", value_name="Coverage Rate")
    plot_df["Style"] = plot_df["Style"].apply(lambda x: "blood" if "blood" in x else "feces")

    # Draw box plot
    sns.boxplot(x="Style", y="Coverage Rate", data=plot_df, ax=ax, palette=colors, width=0.5, showfliers=False, linewidth=0.6)

    # Paired connection lines and circles
    for j in range(len(merged)):
        x_vals = [0, 1]
        y_vals = [merged.iloc[j][f"{depth}_blood"], merged.iloc[j][f"{depth}_feces"]]
        
        # Draw connecting line
        ax.plot(x_vals, y_vals, color="gray", alpha=0.3, linewidth=1)
        
        # Draw circles - blood samples use dark blue fill, feces samples use orange-red fill
        ax.scatter([0], [y_vals[0]], facecolor=colors["blood"], edgecolor="black", 
                  s=10, alpha=0.6, linewidth=0.8, zorder=3)  # blood sample circle
        ax.scatter([1], [y_vals[1]], facecolor=colors["feces"], edgecolor="black", 
                  s=10, alpha=0.6, linewidth=0.8, zorder=3)  # feces sample circle

    # t-test
    t_stat, p_value = stats.ttest_rel(merged[f"{depth}_blood"], merged[f"{depth}_feces"])
    ax.text(0.5, max(plot_df["Coverage Rate"]) + 0.01, f"p = {p_value:.2e}", ha='center', fontsize=10)

    # Set title and labels
    ax.set_title(f"{depth}", fontsize=13)
    ax.set_xlabel("Sample Type", fontsize=11)
    if i == 0:
        ax.set_ylabel("Coverage Rate", fontsize=11)
    else:
        ax.set_ylabel("")

# Adjust layout
plt.tight_layout(w_pad=3.0)  # Increase horizontal spacing to avoid crowding

# Save as PDF
output_file = "step2_paired_boxplots_with_circles.pdf"
plt.savefig(output_file, format="pdf", bbox_inches="tight")
plt.show()

print(f"Image has been saved as: {output_file}")