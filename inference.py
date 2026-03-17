### THE BRAIN   ###

import numpy as np
import pandas as pd
from scipy.stats import poisson

def create_rate_grid(resolution=20, lower=0.01, upper=100.0):
    """
    Creates a 20x20 log-spaced grid for Alpha (dS) and Beta (dN) rates.
    """
    # Log-spaced points ensure we capture low rates (0.01) and high rates (100) with equal fidelity
    grid_points = np.logspace(np.log10(lower), np.log10(upper), resolution)
    
    # Create meshgrid: alpha_grid (rows), beta_grid (cols)
    alpha_mesh, beta_mesh = np.meshgrid(grid_points, grid_points, indexing='ij')
    
    return alpha_mesh, beta_mesh

def calculate_site_likelihoods(df, alpha_mesh, beta_mesh):
    """
    Calculates P(Data | alpha, beta) for every site and every grid point.
    Returns a 3D tensor: (N_sites, 20, 20).
    """
    n_sites = len(df)
    n_grid = alpha_mesh.shape[0]
    
    # Extract data as numpy arrays (N_sites, 1, 1) for broadcasting
    S_obs = df['S_obs'].values[:, None, None]
    N_obs = df['N_obs'].values[:, None, None]
    S_exp = df['S_exp_hky'].values[:, None, None]
    N_exp = df['N_exp_hky'].values[:, None, None]
    
    # Broadcast grid to match sites: (1, 20, 20)
    alpha_broad = alpha_mesh[None, :, :]
    beta_broad = beta_mesh[None, :, :]
    
    # Expected counts for each grid point: Rate * Opportunity
    # Add epsilon to expected values to avoid log(0) errors if S_exp is 0
    epsilon = 1e-12
    lambda_S = (alpha_broad * S_exp) + epsilon
    lambda_N = (beta_broad * N_exp) + epsilon
    
    # Calculate Log-Likelihoods (Poisson)
    # ll_S = log P(S_obs | alpha * S_exp)
    ll_S = poisson.logpmf(S_obs, lambda_S)
    ll_N = poisson.logpmf(N_obs, lambda_N)
    
    # Total Log-Likelihood for the grid point is sum of Syn and Non-syn
    total_ll = ll_S + ll_N
    
    return total_ll

def run_inference(df_counts):
    """
    Main engine: 
    1. Grids the rates.
    2. Calculates Likelihoods.
    3. Learns the Global Prior (Empirical Bayes).
    4. Computes Posterior Probabilities.
    """
    print("Initialising Bayesian Inference Grid (20x20)...")
    
    # 1. Setup Grid
    alpha_mesh, beta_mesh = create_rate_grid()
    
    # 2. Calculate Raw Log-Likelihoods per site
    # Shape: (N_sites, 20, 20)
    log_likelihoods = calculate_site_likelihoods(df_counts, alpha_mesh, beta_mesh)
    
    # Convert Log-Likelihoods to Probabilities (numerically stable softmax)
    # We subtract the max per site to prevent overflow when doing exp()
    max_ll = np.max(log_likelihoods, axis=(1, 2), keepdims=True)
    likelihoods = np.exp(log_likelihoods - max_ll)
    
    # Normalize so each site sums to 1.0 (P(Grid | Data_i))
    site_posteriors_flat = likelihoods / np.sum(likelihoods, axis=(1, 2), keepdims=True)
    
    # 3. Learn Global Prior (Empirical Bayes)
    # Sum the normalized posteriors across ALL sites to see "Global" trends
    global_weights = np.sum(site_posteriors_flat, axis=0)
    global_weights /= np.sum(global_weights) # Normalize to sum to 1
    
    # 4. Calculate Final Posterior for each site
    # P(Grid | Data, Prior) propto P(Data | Grid) * P(Prior)
    # We multiply the site's likelihood surface by the Global Weights
    final_posteriors = likelihoods * global_weights[None, :, :]
    
    # Normalize again
    epsilon2 = 1e-15
    final_posteriors /= (np.sum(final_posteriors, axis=(1, 2), keepdims=True) + epsilon2)
    
    # 5. Extract Metrics
    results = []
    
    # Mask where Beta > Alpha (Positive Selection)
    pos_sel_mask = beta_mesh > alpha_mesh
    
    for i in range(len(df_counts)):
        grid_probs = final_posteriors[i]
        
        # Mean Alpha and Beta (Weighted average across grid)
        mean_alpha = np.sum(grid_probs * alpha_mesh)
        mean_beta = np.sum(grid_probs * beta_mesh)
        
        # Prob(Beta > Alpha) -> Sum of probabilities where beta > alpha
        prob_pos_sel = np.sum(grid_probs[pos_sel_mask])
        
        results.append({
            'Codon_Pos': df_counts.iloc[i]['Codon_Pos'],
            'Alpha_Mean': mean_alpha,
            'Beta_Mean': mean_beta,
            'Prob_Positive_Selection': prob_pos_sel,
            'Label': 'Positive' if prob_pos_sel > 0.95 else 'Neutral/Negative'
        })
        
    return pd.DataFrame(results)