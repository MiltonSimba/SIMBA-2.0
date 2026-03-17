import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Set the aesthetic style for Q1 publication
sns.set_theme(style="white", palette="muted")
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10

# def add_spike_domains(ax, max_val):
#     """Adds coloured shading for functional domains of Spike"""
#     # Standard SARS-CoV-2 Spike Domain Coordinates
#     domains = [
#         (13, 305, 'NTD', '#e1f5fe'),
#         (319, 541, 'RBD', '#fff9c4'),
#         (1, 685, 'S1', '#f3e5f5'),
#         (686, 1273, 'S2', '#e8f5e9')
#     ]
    
#     for start, end, label, color in domains:
#         ax.axvspan(start, end, alpha=0.3, color=color, zorder=0)
#         ax.text((start + end) / 2, max_val * 0.95, label, 
#                 fontweight='bold', ha='center', alpha=0.6)
def add_spike_domains(ax, max_val):
    """Adds coloured shading for functional domains of Spike including RBM"""
    # Standard SARS-CoV-2 Spike Domain Coordinates
    # Format: (Start, End, Label, Colour, Y_Offset_Multiplier)
    domains = [
        (13, 305, 'NTD', '#e1f5fe', 0.95),
        (319, 541, 'RBD', '#fff9c4', 0.95),
        (437, 508, 'RBM', '#ffccbc', 0.88), # Nested inside RBD
        (1, 685, 'S1', '#f3e5f5', 0.05),   # Placed at the bottom
        (686, 1273, 'S2', '#e8f5e9', 0.05) # Placed at the bottom
    ]
    
    for start, end, label, color, y_mult in domains:
        ax.axvspan(start, end, alpha=0.4, color=color, zorder=0)
        ax.text((start + end) / 2, max_val * y_mult, label, 
                fontweight='bold', ha='center', alpha=0.7, fontsize=9)

def plot_selection_landscape(df, output_name):
    """Creates a high-impact Manhattan plot of selection probabilities"""
    plt.figure(figsize=(12, 5))
    ax = sns.scatterplot(
        data=df, 
        x='Codon_Pos', 
        y='Prob_Positive_Selection',
        hue='Prob_Positive_Selection',
        palette='magma',
        size='Beta_Mean',
        sizes=(20, 200),
        edgecolor='black',
        linewidth=0.5,
        legend=None
    )
    
    add_spike_domains(ax, 1.0)
    
    plt.axhline(y=0.95, color='red', linestyle='--', linewidth=1, label='Significance Threshold')
    plt.title("Evolutionary Selection Probability Across the Spike Gene")
    plt.xlabel("Codon Position")
    plt.ylabel("Probability of Positive Selection")
    plt.ylim(0, 1.05)
    plt.tight_layout()
    plt.savefig(f"{output_name}_landscape.png", dpi=300)
    print(f"Saved landscape plot as {output_name}_landscape.png")

def plot_rate_comparison(df, output_name):
    """Scatter plot comparing synonymous and non-synonymous rates"""
    plt.figure(figsize=(7, 7))
    
    # Use log scale for better visualisation of diverse rates
    df_plot = df.copy()
    df_plot['Alpha_Log'] = np.log10(df['Alpha_Mean'] + 1e-3)
    df_plot['Beta_Log'] = np.log10(df['Beta_Mean'] + 1e-3)

    sns.scatterplot(
        data=df_plot,
        x='Alpha_Log',
        y='Beta_Log',
        hue='Label',
        palette={'Positive': '#d32f2f', 'Neutral/Negative': '#455a64'},
        alpha=0.6
    )
    
    # The Neutrality Line (Alpha = Beta)
    lims = [min(df_plot['Alpha_Log'].min(), df_plot['Beta_Log'].min()),
            max(df_plot['Alpha_Log'].max(), df_plot['Beta_Log'].max())]
    plt.plot(lims, lims, color='black', linestyle='--', alpha=0.5, label='Neutral Evolution')
    
    plt.title("Comparison of Synonymous and Non-synonymous Substitution Rates")
    plt.xlabel("Log10 Synonymous Rate (Alpha)")
    plt.ylabel("Log10 Non-synonymous Rate (Beta)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{output_name}_rates.png", dpi=300)
    print(f"Saved rate comparison as {output_name}_rates.png")

def generate_visuals(csv_file):
    """Main function to load data and produce all visuals"""
    df = pd.read_csv(csv_file)
    base_name = csv_file.replace(".csv", "")
    
    plot_selection_landscape(df, base_name)
    plot_rate_comparison(df, base_name)

if __name__ == "__main__":
    # You can run this directly if you have the CSV
    generate_visuals("simba_results_zimbabwe.csv")