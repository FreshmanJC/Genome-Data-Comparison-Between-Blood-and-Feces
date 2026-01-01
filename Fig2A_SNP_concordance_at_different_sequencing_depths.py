import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib import rcParams
import matplotlib.colors as mcolors
from scipy.stats import t, f

# Set font and vector output
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42

# Step 1: Read data and preprocess
df = pd.read_excel("step13_1_ERROR_rate.xlsx")
df['Sequence depth'] = df['Sequence depth'].astype(str).str.replace('X', '').astype(int)
df_sorted = df.sort_values(by='Sequence depth')
x = df_sorted['Sequence depth'].values
y = df_sorted['Genotyping Error Rate'].values  # Actually F1-score

# Step 2: Define multiple fitting models
def exp_decay(x, a, b, k):
    return a * np.exp(-b * x) + k

def power_decay(x, a, b, k):
    return a * (x ** (-b)) + k

def rational_function(x, a, b, c, d):
    return (a * x + b) / (c * x + d)

def sigmoid_decay(x, a, b, c, d):
    return a / (1 + b * np.exp(c * x)) + d

models = {
    'Exponential Decay': (exp_decay, [max(y), 0.1, min(y)]),
    'Power Decay': (power_decay, [max(y), 0.5, min(y)]),
    'Rational Function': (rational_function, [1, max(y), 0.1, 1]),
    'Sigmoid Decay': (sigmoid_decay, [max(y)-min(y), 1, -0.1, min(y)])
}

# Step 3: Fit and select best model
best_model, best_r2, best_params, best_name = None, -np.inf, None, ""
for name, (func, p0) in models.items():
    try:
        popt, pcov = curve_fit(func, x, y, p0=p0, maxfev=10000)
        y_pred = func(x, *popt)
        
        # Calculate R²
        ss_res = np.sum((y - y_pred)**2)
        ss_tot = np.sum((y - np.mean(y))**2)
        r2 = 1 - (ss_res / ss_tot)
        
        if r2 > best_r2:
            best_r2, best_model, best_params, best_name = r2, func, popt, name
        print(f"{name}: R² = {r2:.4f}")
    except:
        print(f"{name}: Fitting failed")

# Step 4: Prepare fitting curve and confidence interval function
def calculate_confidence_interval(x_vals, func, popt, pcov, alpha=0.05):
    dof = max(0, len(x_vals) - len(popt))
    tval = t.ppf(1.0 - alpha / 2., dof)
    delta = 1e-6
    J = np.array([
        [ (func(x, *(popt + delta*np.eye(len(popt))[i])) - func(x, *popt)) / delta
          for i in range(len(popt)) ] for x in x_vals
    ])
    y_sigma = np.sqrt(np.sum(J @ pcov * J, axis=1))
    y_fit = func(x_vals, *popt)
    return y_fit + tval * y_sigma, y_fit - tval * y_sigma

# Step 5: Calculate F-value and p-value
# Calculate F-value and p-value
n = len(x)  # Number of samples
k = len(best_params)  # Degrees of freedom
ss_res = np.sum((y - best_model(x, *best_params))**2)
ss_tot = np.sum((y - np.mean(y))**2)
f_stat = (ss_tot - ss_res) / k / (ss_res / (n - k - 1))
p_value = 1 - f.cdf(f_stat, k, n - k - 1)

# Get fitting parameter covariance matrix, calculate confidence interval
x_fit = np.linspace(min(x), max(x), 300)
y_fit = best_model(x_fit, *best_params)

try:
    _, pcov = curve_fit(best_model, x, y, p0=models[best_name][1], maxfev=10000)
    y_upper, y_lower = calculate_confidence_interval(x_fit, best_model, best_params, pcov)
except Exception as e:
    print(f"Failed to calculate confidence interval: {e}")
    y_upper, y_lower = None, None

# Plotting
plt.figure(figsize=(5, 4))
plt.scatter(x, y, color='#8d8d8e', s=100,  linewidth=0.5, alpha=0.6, label='Samples')

plt.plot(x_fit, y_fit, 'r-', linewidth=2.5, label=f'{best_name} Fit')

if y_upper is not None and y_lower is not None:
    plt.fill_between(x_fit, y_lower, y_upper, color='#3b66bc', alpha=0.25, label='95% Confidence Interval')

plt.title('F1-Score vs Sequencing Depth (Best Fit Model)', fontsize=14, fontweight='bold')
plt.xlabel('Sequencing Depth (X)', fontsize=12)
plt.ylabel('Genotyping Error Rate', fontsize=12)
unique_depths = sorted(df['Sequence depth'].unique())
plt.xticks(unique_depths, rotation=45 if len(unique_depths) > 8 else 0)
plt.grid(True, linestyle='--', alpha=0.3)
plt.legend(frameon=True, fancybox=True, shadow=True)

plt.text(0.05, 0.95, f'R² = {best_r2:.4f}', transform=plt.gca().transAxes,
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
         verticalalignment='top', fontsize=10)

# Display p-value
plt.text(0.05, 0.90, f'p = {p_value:.2e}', transform=plt.gca().transAxes,
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
         verticalalignment='top', fontsize=10)

plt.tight_layout()
output_file = "step13_1_errorcurve.pdf"
plt.savefig(output_file, format="pdf", bbox_inches="tight", dpi=300)
plt.show()

# Step 6: Print results
if best_model is not None:
    print(f"\nBest fitting model: {best_name} (R² = {best_r2:.4f}, p = {p_value:.2e})")
    print(f"Fitting parameters: {best_params}")
    if best_name == 'Exponential Decay':
        a, b, k = best_params
        print(f"F1-score = {a:.4f} * exp(-{b:.4f} * x) + {k:.4f}")
    elif best_name == 'Power Decay':
        a, b, k = best_params
        print(f"F1-score = {a:.4f} * x^(-{b:.4f}) + {k:.4f}")
    elif best_name == 'Rational Function':
        a, b, c, d = best_params
        print(f"F1-score = ({a:.4f} * x + {b:.4f}) / ({c:.4f} * x + {d:.4f})")
    elif best_name == 'Sigmoid Decay':
        a, b, c, d = best_params
        print(f"F1-score = {a:.4f} / (1 + {b:.4f} * exp({c:.4f} * x)) + {d:.4f}")