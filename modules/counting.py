### THE ENGINE ###

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Data import CodonTable

def filter_sequences(fasta_file, threshold=0.05):
    """
    Filters sequences based on the proportion of 'N' or gap characters.
    """
    sequences = []
    # Parsing sequences (Generator is memory efficient)
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        if len(seq) == 0: continue
        
        # Calculate ambiguity
        ambiguity = (seq.count('N') + seq.count('-')) / len(seq)
        if ambiguity <= threshold:
            sequences.append(seq)
            
    print(f"Filtered Alignment: Kept {len(sequences)} sequences.")
    return sequences

def calculate_base_frequencies(sequences):
    """
    Calculates global base frequencies (pi) from the filtered alignment.
    """
    # Use a Counter or simple sum for memory efficiency with large lists
    total_len = 0
    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    
    for s in sequences:
        total_len += len(s)
        counts['A'] += s.count('A')
        counts['C'] += s.count('C')
        counts['G'] += s.count('G')
        counts['T'] += s.count('T')
    
    return {k: v / total_len for k, v in counts.items()}

def calculate_hky_opportunities(ref_seq, kappa, pi):
    """
    Calculates S_exp and N_exp using HKY85 weights.
    """
    table = CodonTable.unambiguous_dna_by_id[1]
    # Create a fast translation map
    forward_table = table.forward_table
    stop_codons = table.stop_codons
    
    n_codons = len(ref_seq) // 3
    S_exp_site = np.zeros(n_codons)
    N_exp_site = np.zeros(n_codons)

    bases = ['A', 'C', 'G', 'T']

    for i in range(n_codons):
        p = i * 3
        codon = ref_seq[p:p+3]
        
        if 'N' in codon or '-' in codon or len(codon) != 3:
            continue
            
        ref_aa = forward_table.get(codon, '*')
        
        # Iterate over 3 positions in the codon
        for pos in range(3):
            original_base = codon[pos]
            
            for target_base in bases:
                if target_base == original_base: continue
                
                # Check Transition (Ts) vs Transversion (Tv)
                is_ts = (original_base in 'AG' and target_base in 'AG') or \
                        (original_base in 'CT' and target_base in 'CT')
                
                weight = (kappa * pi[target_base]) if is_ts else pi[target_base]
                
                # Construct mutated codon
                mut_list = list(codon)
                mut_list[pos] = target_base
                mut_codon = "".join(mut_list)
                
                if mut_codon in stop_codons: continue # Skip stop codons
                
                mut_aa = forward_table.get(mut_codon, '*')
                
                if mut_aa == ref_aa:
                    S_exp_site[i] += weight
                else:
                    N_exp_site[i] += weight
                    
    return S_exp_site, N_exp_site

def observe_mutations_vectorised(aln_matrix, ref_seq):
    """
    Optimized counting using NumPy broadcasting.
    Input: 
      aln_matrix: (N_seq, L_bases) numpy array of characters
      ref_seq: string
    """
    n_seqs, n_bases = aln_matrix.shape
    n_codons = n_bases // 3
    
    # Pre-compute translation table for speed
    table = CodonTable.unambiguous_dna_by_id[1]
    forward_table = table.forward_table
    
    S_obs = np.zeros(n_codons)
    N_obs = np.zeros(n_codons)
    
    # Convert Reference to Array for broadcasting
    ref_arr = np.array(list(ref_seq))
    
    print("Starting vectorized counting...")
    
    # Iterate over CODONS (stride of 3), not nucleotides
    for i in range(n_codons):
        p = i * 3
        ref_codon_arr = ref_arr[p:p+3]
        ref_codon_str = ref_seq[p:p+3]
        
        if 'N' in ref_codon_str or '-' in ref_codon_str: continue
        ref_aa = forward_table.get(ref_codon_str, '*')

        # Extract the column for this codon (N_seqs x 3)
        site_codons = aln_matrix[:, p:p+3]
        
        # 1. Quick Filter: Find sequences that do NOT match the reference at this site
        # casting to string for comparison is safer than assuming simple broadcasting matches
        # (N, 3) != (3,) broadcasts correctly
        is_mutated = np.any(site_codons != ref_codon_arr, axis=1)
        
        if not np.any(is_mutated):
            continue
            
        # 2. Strict Masking: Ignore 'N' or '-' in the mutated sequences
        has_ambiguity = np.any((site_codons == 'N') | (site_codons == '-'), axis=1)
        
        # We only care about: Mutated AND Clean
        valid_mutants_idx = np.where(is_mutated & ~has_ambiguity)[0]
        
        if len(valid_mutants_idx) == 0:
            continue
            
        # 3. Process only the valid mutants (Sparse loop is fast)
        for idx in valid_mutants_idx:
            mut_codon_arr = site_codons[idx]
            # Fast join
            mut_codon_str = mut_codon_arr[0] + mut_codon_arr[1] + mut_codon_arr[2]
            
            mut_aa = forward_table.get(mut_codon_str, '*')
            
            if mut_aa == ref_aa:
                S_obs[i] += 1
            else:
                N_obs[i] += 1
                
    return S_obs, N_obs

def process_alignment(fasta_file, ref_seq, kappa=4.0):
    sequences = filter_sequences(fasta_file)
    if not sequences:
        raise ValueError("No sequences passed the quality filter.")
        
    pi = calculate_base_frequencies(sequences)
    print(f"Global Frequencies: {pi}")

    # Convert to 2D Char Array
    # Ensure all sequences are same length as reference (truncate or pad if needed)
    # For GISAID, we assume alignment is pre-aligned.
    L = len(ref_seq)
    valid_sequences = [list(s[:L]) for s in sequences if len(s) >= L]
    aln_matrix = np.array(valid_sequences)

    S_exp, N_exp = calculate_hky_opportunities(ref_seq, kappa, pi)
    S_obs, N_obs = observe_mutations_vectorised(aln_matrix, ref_seq)

    return pd.DataFrame({
        'Codon_Pos': np.arange(1, len(S_exp) + 1),
        'S_obs': S_obs,
        'N_obs': N_obs,
        'S_exp_hky': S_exp,
        'N_exp_hky': N_exp
    })
















# Copilot generation of counting.py
# import numpy as np
# import pandas as pd
# from Bio import SeqIO
# from Bio.Data import CodonTable

# def filter_sequences(fasta_file, threshold=0.05):
#     """
#     Filters sequences based on the proportion of 'N' or gap characters.

#     Parameters:
#         fasta_file (str): Path to the FASTA file.
#         threshold (float): Maximum allowed proportion of 'N' or gaps.

#     Returns:
#         list: Filtered list of sequences.
#     """
#     sequences = []
#     for record in SeqIO.parse(fasta_file, "fasta"):
#         seq = str(record.seq).upper()
#         if (seq.count('N') + seq.count('-')) / len(seq) <= threshold:
#             sequences.append(seq)
#     return sequences

# def calculate_base_frequencies(sequences):
#     """
#     Calculates global base frequencies (pi_A, pi_C, pi_G, pi_T).

#     Parameters:
#         sequences (list): List of sequences.

#     Returns:
#         dict: Base frequencies.
#     """
#     concatenated = ''.join(sequences)
#     total_bases = len(concatenated)
#     return {
#         'A': concatenated.count('A') / total_bases,
#         'C': concatenated.count('C') / total_bases,
#         'G': concatenated.count('G') / total_bases,
#         'T': concatenated.count('T') / total_bases
#     }

# def calculate_hky_opportunities(ref_seq, kappa, pi):
#     """
#     Calculates synonymous and non-synonymous opportunities using the HKY85 model.

#     Parameters:
#         ref_seq (str): Reference sequence.
#         kappa (float): Transition/Transversion ratio.
#         pi (dict): Base frequencies.

#     Returns:
#         tuple: Arrays of synonymous and non-synonymous opportunities.
#     """
#     codon_table = CodonTable.unambiguous_dna_by_id[1]
#     S_exp_site = np.zeros(len(ref_seq) // 3)
#     N_exp_site = np.zeros(len(ref_seq) // 3)

#     for i in range(0, len(ref_seq), 3):
#         codon = ref_seq[i:i+3]
#         if '-' in codon or 'N' in codon:
#             continue

#         for pos in range(3):
#             for base in 'ACGT':
#                 if base == codon[pos]:
#                     continue

#                 mutated_codon = list(codon)
#                 mutated_codon[pos] = base
#                 mutated_codon = ''.join(mutated_codon)

#                 if '-' in mutated_codon or 'N' in mutated_codon:
#                     continue

#                 is_transition = (codon[pos] in 'AG' and base in 'AG') or (codon[pos] in 'CT' and base in 'CT')
#                 weight = kappa * pi[base] if is_transition else pi[base]

#                 if mutated_codon in codon_table.synonyms[codon]:
#                     S_exp_site[i // 3] += weight
#                 else:
#                     N_exp_site[i // 3] += weight

#     return S_exp_site, N_exp_site

# def observe_mutations_numpy(aln_matrix, ref_vector):
#     """
#     Observes mutations in the alignment matrix compared to the reference vector.

#     Parameters:
#         aln_matrix (np.ndarray): Alignment matrix (N_seq, L_codon).
#         ref_vector (np.ndarray): Reference sequence vector.

#     Returns:
#         tuple: Observed synonymous and non-synonymous counts.
#     """
#     mask = (ref_vector != 'N') & (ref_vector != '-')
#     aln_matrix = aln_matrix[:, mask]
#     ref_vector = ref_vector[mask]

#     S_obs = np.zeros(ref_vector.shape[0])
#     N_obs = np.zeros(ref_vector.shape[0])

#     for i in range(ref_vector.shape[0]):
#         ref_codon = ref_vector[i]
#         for seq in aln_matrix:
#             if seq[i] == ref_codon:
#                 continue
#             if seq[i] == 'N' or seq[i] == '-':
#                 continue

#             if seq[i] in CodonTable.unambiguous_dna_by_id[1].synonyms[ref_codon]:
#                 S_obs[i] += 1
#             else:
#                 N_obs[i] += 1

#     return S_obs, N_obs

# def process_alignment(fasta_file, ref_seq, kappa=4.0):
#     """
#     Processes the alignment and calculates mutation statistics.

#     Parameters:
#         fasta_file (str): Path to the FASTA file.
#         ref_seq (str): Reference sequence.
#         kappa (float): Transition/Transversion ratio.

#     Returns:
#         pd.DataFrame: DataFrame containing mutation statistics.
#     """
#     sequences = filter_sequences(fasta_file)
#     pi = calculate_base_frequencies(sequences)

#     aln_matrix = np.array([list(seq) for seq in sequences])
#     ref_vector = np.array(list(ref_seq))

#     S_exp, N_exp = calculate_hky_opportunities(ref_seq, kappa, pi)
#     S_obs, N_obs = observe_mutations_numpy(aln_matrix, ref_vector)

#     return pd.DataFrame({
#         'Codon_Pos': np.arange(1, len(S_exp) + 1),
#         'S_obs': S_obs,
#         'N_obs': N_obs,
#         'S_exp_hky': S_exp,
#         'N_exp_hky': N_exp
#     })