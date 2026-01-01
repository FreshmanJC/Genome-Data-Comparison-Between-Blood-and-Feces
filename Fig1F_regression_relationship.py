import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
from matplotlib import rcParams

# Set font vector output
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42

# Read Excel file
df = pd.read_excel('step14.CTvalue.xlsx')

# Data transformation
x = df['EEN']
y = df['10X Coverage %'] * 100  # Convert to percentage form

# Define different fitting functions
def log_model(x, a, b, c):
    """Logarithmic model: y = a * log(x) + b + c"""
    return a * np.log10(x) + b

def saturation_model(x, a, b, c):
    """Saturation model (similar to Michaelis-Menten): y = a * x / (b + x) + c"""
    return a * x / (b + x) + c

def exp_saturation_model(x, a, b, c):
    """Exponential saturation model: y = a * (1 - exp(-b * x)) + c"""
    return a * (1 - np.exp(-b * x)) + c

def power_model(x, a, b, c):
    """Power function model: y = a * x^b + c"""
    return a * np.power(x, b) + c

# Fit different models and select the best one
models = {
    'Log Model': log_model,
    'Saturation Model': saturation_model,
    'Power Model': power_model
}

best_model = None
best_r2 = -1
best_params = None
best_name = ""

results = {}

for name, model_func in models.items():
    try:
        # Fit model
        if name == 'Log Model':
            # Initial parameters for logarithmic model
            popt, pcov = curve_fit(model_func, x, y, 
                                 p0=[10, 90], 
                                 maxfev=10000)
        elif name == 'Saturation Model':
            # Initial parameters for saturation model
            popt, pcov = curve_fit(model_func, x, y, 
                                 p0=[100, 0.001, 80], 
                                 maxfev=10000)
        elif name == 'Power Model':
            # Initial parameters for power function model
            popt, pcov = curve_fit(model_func, x, y, 
                                 p0=[50, 0.1, 80], 
                                 maxfev=10000)
        
        # Predicted values
        y_pred = model_func(x, *popt)
        
        # Calculate R² value
        r2 = r2_score(y, y_pred)
        
        # Calculate p-value
        n = len(y)
        k = len(popt)
        mse = np.mean((y - y_pred) ** 2)
        
        # F statistic
        if r2 < 1:
            f_stat = (r2 / (k-1)) / ((1 - r2) / (n - k))
            p_value = 1 - stats.f.cdf(f_stat, k-1, n - k)
        else:
            p_value = 0.0
        
        results[name] = {
            'params': popt,
            'r2': r2,
            'p_value': p_value,
            'model_func': model_func
        }
        
        # Select best model
        if r2 > best_r2:
            best_r2 = r2
            best_model = model_func
            best_params = popt
            best_name = name
            
        print(f"{name}: R² = {r2:.4f}, P-value = {p_value:.2E}")
        
    except Exception as e:
        print(f"{name} fitting failed: {e}")
        continue

# Use best model for prediction and confidence interval calculation
if best_model is not None:
    print(f"\nBest model: {best_name}")
    print(f"R² = {best_r2:.4f}")
    
    # Generate dense points for fitting curve
    x_fit = np.logspace(np.log10(x.min()), np.log10(x.max()), 300)
    y_fit = best_model(x_fit, *best_params)
    
    # Bootstrap confidence interval calculation
    def bootstrap_confidence_interval(x_data, y_data, model_func, params, x_pred, n_bootstrap=1000, confidence=0.95):
        n_samples = len(x_data)
        predictions = []
        
        for _ in range(n_bootstrap):
            # Bootstrap sampling
            indices = np.random.choice(n_samples, n_samples, replace=True)
            x_boot = x_data[indices]
            y_boot = y_data[indices]
            
            try:
                # Refit model
                if best_name == 'Log Model':
                    popt_boot, _ = curve_fit(model_func, x_boot, y_boot, 
                                           p0=params, maxfev=5000)
                elif best_name == 'Saturation Model':
                    popt_boot, _ = curve_fit(model_func, x_boot, y_boot, 
                                           p0=params, maxfev=5000)
                elif best_name == 'Power Model':
                    popt_boot, _ = curve_fit(model_func, x_boot, y_boot, 
                                           p0=params, maxfev=5000)
                
                # Predict
                y_pred_boot = model_func(x_pred, *popt_boot)
                predictions.append(y_pred_boot)
            except:
                # If fitting fails, use original parameters
                y_pred_boot = model_func(x_pred, *params)
                predictions.append(y_pred_boot)
        
        predictions = np.array(predictions)
        
        # Calculate confidence interval
        alpha = 1 - confidence
        lower_percentile = (alpha / 2) * 100
        upper_percentile = (1 - alpha / 2) * 100
        
        ci_lower = np.percentile(predictions, lower_percentile, axis=0)
        ci_upper = np.percentile(predictions, upper_percentile, axis=0)
        
        return ci_lower, ci_upper
    
    # Calculate confidence interval
    ci_lower, ci_upper = bootstrap_confidence_interval(x.values, y.values, best_model, best_params, x_fit)
    
    # Create figure
    plt.figure(figsize=(5, 4))
    
    # Draw confidence interval
    plt.fill_between(x_fit, ci_lower, ci_upper, alpha=0.2, color='#1f778d', label='95% CI')
    
    # Draw scatter plot
    plt.scatter(x, y, color='#1f778d', alpha=0.9, s=70, zorder=5)
    
    # Draw fitting curve
    plt.plot(x_fit, y_fit, linestyle='--', color='black', linewidth=4, zorder=4)
    
    # Set x-axis as logarithmic scale
    plt.xscale('log')
    plt.xticks([5e-5, 5e-4, 5e-3], ['5e-5', '5e-4', '5e-3'], fontsize=12)
    
    # Set y-axis ticks and range
    plt.ylim(75, 100)
    plt.yticks(np.arange(77.5, 100, 2.5), fontsize=12)
    
    # Add axis labels
    plt.xlabel('EEN', fontsize=14, fontweight='bold')
    plt.ylabel('10X Coverage %', fontsize=14, fontweight='bold')
    
    # Add axis tick lines
    plt.tick_params(axis='both', which='major', direction='in', length=6, width=1.5)
    plt.tick_params(axis='both', which='minor', direction='in', length=3, width=1)
    plt.minorticks_on()
    
    # Add grid
    plt.grid(True, linestyle='--', alpha=0.3)
    
    # Add statistical information text box
    plt.text(0.95, 0.95,
             f'$R^2$ = {best_r2:.3f}\n$P$ = {results[best_name]["p_value"]:.2E}\nModel: {best_name}',
             transform=plt.gca().transAxes,
             fontsize=10,
             verticalalignment='top',
             horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='white', edgecolor='black', linewidth=1.5))
    
    # Optimize layout and save
    plt.tight_layout()
    output_file = f"step14_EEN_{best_name.replace(' ', '_').lower()}.pdf"
    plt.savefig(output_file, format="pdf", bbox_inches="tight")
    plt.show()
    
    print(f"\nChart has been saved as: {output_file}")
    print(f"Model used: {best_name}")
    print(f"Model parameters: {best_params}")
    
    # Check trend at highest EEN value
    x_max = x.max()
    x_test_high = np.array([x_max * 0.9, x_max, x_max * 1.1])
    y_test_high = best_model(x_test_high, *best_params)
    
    print(f"\nHigh EEN value trend check:")
    for i, (x_val, y_val) in enumerate(zip(x_test_high, y_test_high)):
        print(f"EEN = {x_val:.2E}, Coverage = {y_val:.2f}%")
    
    # Calculate slope change
    slopes = np.diff(y_test_high) / np.diff(x_test_high)
    print(f"Slope: {slopes}")
    
    if np.all(slopes < 1e-6):  # Approximately stable
        print("✓ Curve tends to stabilize at high EEN values")
    elif slopes[-1] > 0:
        print("! Curve is still rising, model adjustment may be needed")
    else:
        print("✗ Curve is still declining, consider trying other models")

else:
    print("All model fittings failed, please check data or adjust parameters")