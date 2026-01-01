import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib import rcParams
from scipy.stats import t, f
import matplotlib.colors as mcolors

# Set font and vector output
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42

# Step 1: Read data
df = pd.read_excel("step13_recall.xlsx")
df['Sequence depth'] = df['Sequence depth'].astype(str).str.replace('X', '').astype(int)
df_sorted = df.sort_values(by='EEN', ascending=True)

x = df_sorted['Sequence depth'].values
y = df_sorted['Recall'].values
c = df_sorted['EEN'].values

# Step 2: Define fitting function (exponential saturation model)
def exp_saturation(x, a, b, k):
    return a - b * np.exp(-k * x)

# Step 3: Fit model and get covariance matrix
p0 = [1.0, 0.5, 0.1]
popt, pcov = curve_fit(exp_saturation, x, y, p0=p0, maxfev=10000)
a, b, k = popt

# Step 4: Construct fitting curve
x_fit = np.linspace(min(x), max(x), 300)
y_fit = exp_saturation(x_fit, *popt)

# Step 5: Calculate 95% confidence interval
def calculate_confidence_interval(x_vals, func, popt, pcov, alpha=0.05):
    dof = max(0, len(x_vals) - len(popt))
    tval = t.ppf(1.0 - alpha / 2., dof)
    delta = 1e-6
    J = np.array([
        [ (func(xi, *(popt + delta*np.eye(len(popt))[i])) - func(xi, *popt)) / delta
          for i in range(len(popt)) ]
        for xi in x_vals
    ])
    sigma = np.sqrt(np.sum(J @ pcov * J, axis=1))
    y_upper = func(x_vals, *popt) + tval * sigma
    y_lower = func(x_vals, *popt) - tval * sigma
    return y_upper, y_lower

y_upper, y_lower = calculate_confidence_interval(x_fit, exp_saturation, popt, pcov)

# Step 6: Calculate R² and p-value
# Calculate R²
residuals = y - exp_saturation(x, *popt)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((y - np.mean(y))**2)
r2 = 1 - (ss_res / ss_tot)

# Calculate F statistic and p-value
n = len(x)
k = len(popt)
f_stat = (r2 / k) / ((1 - r2) / (n - k - 1))
p_value = 1 - f.cdf(f_stat, k, n - k - 1)

# Step 7: Plotting
vmin, vmax = 1e-5, 9e-3
boundaries = [vmin, 5e-5, 9e-5, 4e-4, 9e-4, 4e-3, 9e-3]
norm = mcolors.BoundaryNorm(boundaries=boundaries, ncolors=256)

plt.figure(figsize=(6, 4))

# Sample scatter plot
sc = plt.scatter(x, y, c=c, cmap='viridis', norm=norm,
                 s=100, edgecolors='k', alpha=0.8, label='Samples')

# Fitting curve
plt.plot(x_fit, y_fit, 'r-', linewidth=2, label='Exponential Saturation Fit')

# Add confidence interval shading
plt.fill_between(x_fit, y_lower, y_upper, color='red', alpha=0.1,
                 label='95% Confidence Interval')

# Labels and title
plt.title('Recall vs Sequencing Depth (Saturation Model Fit)')
plt.xlabel('Sequencing Depth (X)')
plt.ylabel('Recall')
plt.xticks(sorted(df['Sequence depth'].unique()))
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend()

# Add colorbar
cbar = plt.colorbar(sc)
cbar.set_label('EEN (Enrichment Efficiency)', rotation=270, labelpad=15)

# Add R² and p-value text
plt.text(0.05, 0.95,
         f'R² = {r2:.4f}\np = {p_value:.2e}',
         transform=plt.gca().transAxes,
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
         verticalalignment='top', fontsize=10)

# Save image
plt.tight_layout()
plt.savefig("step13_Recall_with_CI.pdf", format="pdf", bbox_inches="tight")
plt.show()

# Output fitting function and statistical information
print(f"Fitting function form: Recall = {a:.4f} - {b:.4f} * exp(-{k:.4f} * x)")
print(f"R² = {r2:.4f}, p-value = {p_value:.2e}")