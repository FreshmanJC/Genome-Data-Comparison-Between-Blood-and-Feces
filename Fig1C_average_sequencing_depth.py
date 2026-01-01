# JupyterLab Sequencing Depth Difference Visualization Script
from matplotlib import rcParams
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter
import warnings
import matplotlib.font_manager as fm
warnings.filterwarnings('ignore')

# Set font and vector output
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
rcParams['font.size'] = 12

# Solve font issue, set font to system-common DejaVu Sans
rcParams['font.family'] = 'DejaVu Sans'

# Set matplotlib parameters to ensure high-quality graphics display in JupyterLab
plt.rcParams.update({
    'font.size': 12,
    'axes.linewidth': 1.5,
    'xtick.major.width': 1.5,
    'ytick.major.width': 1.5,
    'xtick.minor.width': 1.0,
    'ytick.minor.width': 1.0,
    'lines.linewidth': 1.5,
    'patch.linewidth': 1.5,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1
})

def load_data(file_path):
    """
    Load Excel data
    """
    print("Loading data...")
    df = pd.read_excel(file_path)
    
    # Clean column names
    df.columns = df.columns.str.strip()
    print(f"Data dimensions: {df.shape}")
    print(f"Column names: {df.columns.tolist()}")
    
    # Calculate midpoint positions for x-axis
    df['Position_Mid'] = (df.iloc[:, 0] + df.iloc[:, 1]) / 2
    df['Position_Kb'] = df['Position_Mid'] / 1000  # Convert to Kb
    
    return df

def smooth_data(data, window_size=101, polyorder=3):
    """
    Smooth data using Savitzky-Golay filter
    """
    if len(data) < window_size:
        window_size = len(data) if len(data) % 2 == 1 else len(data) - 1
    if window_size < polyorder + 2:
        polyorder = max(1, window_size - 2)
    
    return savgol_filter(data, window_size, polyorder)

def create_depth_line_plot(df, smooth=True, figsize=(8, 6), save_pdf=False, pdf_name="depth_analysis.pdf", fecal_color='#d3694e', blood_color='#234451'):
    """
    Create sequencing depth line plot
    
    Parameters:
    df: DataFrame
    smooth: Whether to apply data smoothing
    figsize: Figure size
    save_pdf: Whether to save as PDF
    pdf_name: PDF file name
    fecal_color: Color for fecal samples (default #dc6c50)
    blood_color: Color for blood samples (default #244552)
    """
    
    # Extract data
    x_pos = df['Position_Kb'].values
    fecal_data = df.iloc[:, 3].values  # Third column: fecal sample depth
    blood_data = df.iloc[:, 2].values  # Fourth column: blood sample depth
    
    print(f"Plotting {len(fecal_data)} data points")
    
    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    
    # Data smoothing (suitable for large datasets)
    if smooth and len(fecal_data) > 1000:
        print("Applying data smoothing...")
        # Adaptive window size
        window_size = min(501, len(fecal_data) // 100)
        if window_size % 2 == 0:
            window_size += 1
        
        fecal_smooth = smooth_data(fecal_data, window_size)
        blood_smooth = smooth_data(blood_data, window_size)
        
        # Draw main smoothed lines
        ax.plot(x_pos, fecal_smooth, color=fecal_color, linewidth=2.5, label='Fecal sample', alpha=1, zorder=3)
        ax.plot(x_pos, blood_smooth, color=blood_color, linewidth=2.5, label='Blood sample', alpha=1, zorder=3)
        
        # Add transparent background for original data (optional)
        ax.plot(x_pos, fecal_data, color=fecal_color, linewidth=0.3, alpha=0.2, zorder=1)
        ax.plot(x_pos, blood_data, color=blood_color, linewidth=0.3, alpha=0.2, zorder=1)
    else:
        # Directly plot original data
        ax.plot(x_pos, fecal_data, color=fecal_color, linewidth=2.0, label='Fecal sample', alpha=1)
        ax.plot(x_pos, blood_data, color=blood_color, linewidth=2.0, label='Blood sample', alpha=1)
    
    # Set axis labels and title
    ax.set_xlabel('Chromosomal Position (kb)', fontsize=16, fontweight='regular')
    ax.set_ylabel('Average Sequencing Depth', fontsize=16, fontweight='regular')
    
    # Axis beautification
    ax.tick_params(axis='both', which='major', labelsize=6)
    ax.spines['top'].set_visible(True)   # Show top axis
    ax.spines['right'].set_visible(True) # Show right axis
    ax.spines['top'].set_linewidth(1.5)
    ax.spines['right'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    
    # Hide top and right ticks
    ax.tick_params(axis='x', which='both', bottom=False, top=False)
    ax.tick_params(axis='y', which='both', left=False, right=False)
    
    # Smart x-axis tick setting
    n_points = len(x_pos)
    if n_points > 10000:
        n_ticks = 8
    elif n_points > 1000:
        n_ticks = 10
    else:
        n_ticks = min(10, n_points // 100)
    
    if n_ticks > 1:
        tick_positions = np.linspace(x_pos.min(), x_pos.max(), n_ticks)
        ax.set_xticks(tick_positions)
        ax.set_xticklabels([f'{int(pos)}' for pos in tick_positions])
    
    # Adjust layout
    plt.tight_layout()
    
    # Display basic statistical information
    fecal_mean = np.mean(fecal_data)
    blood_mean = np.mean(blood_data)
    fold_change = blood_mean / fecal_mean
    
    print(f"\n=== Basic Statistical Information ===")
    print(f"Fecal sample average depth: {fecal_mean:.3f}")
    print(f"Blood sample average depth: {blood_mean:.3f}")
    print(f"Fold change (Blood/Fecal): {fold_change:.2f}")
    
    # Save PDF (optional)
    if save_pdf:
        plt.savefig(pdf_name, format='pdf', bbox_inches='tight')
        print(f"Chart has been saved as: {pdf_name}")
    
    # Display figure in JupyterLab
    plt.show()
    
    return fig, ax

# ====================
# Main execution code
# ====================

# 1. Load data
df = load_data("AT_1367.xlsx")

# 2. Create line plot
fig, ax = create_depth_line_plot(
    df, 
    smooth=True,           # Whether to smooth data (recommended for large datasets)
    figsize=(8, 6),        # Figure size
    save_pdf=True,         # Whether to save PDF
    pdf_name="sequencing_depth_comparison1.pdf",  # PDF file name
    fecal_color='#d3694e', # Color for fecal samples
    blood_color='#234451'  # Color for blood samples
)