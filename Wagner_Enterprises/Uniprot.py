import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

from Protein import Protein
from GlobalCounts import sort_dictionary, character_count_map, two_char, three_char, four_char, two_character_count_map, three_character_count_map, four_character_count_map

def read_protein_sequences_from_file(directory, file):
    try:
        filename = f"{directory}/{file}"
        with open(filename, 'r') as fh:
            content = fh.read().replace(" ", "").splitlines()
        
        sequences = []
        current_seq = []
        
        for i, line in enumerate(content):
            if line[0] == '>':
                if current_seq:  # Save previous sequence if exists
                    sequences.append(''.join(current_seq))
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_seq:  # Add the last sequence
            sequences.append(''.join(current_seq))
            
        return sequences
    except Exception as e:
        print(f"Error reading file {file}: {str(e)}")
        return []

def get_top_n_peptides(protein_sequences, protein_class_name, n_chars, top_n=20):
    """Generic function to get top N peptides of specified length"""
    combined_peptides = {}
    
    for seq in protein_sequences:
        protein_obj = Protein(seq, protein_class_name)
        peptide_dict = getattr(protein_obj, f"{n_chars}_character_dict")
        
        for peptide, count in peptide_dict.items():
            combined_peptides[peptide] = combined_peptides.get(peptide, 0) + count
    
    sorted_peptides = dict(sorted(combined_peptides.items(), key=lambda item: item[1], reverse=True))
    return dict(list(sorted_peptides.items())[:top_n])

def create_dipeptide_heatmap(protein_sequences, protein_class_name):
    """Create a heatmap of dipeptide frequencies"""
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'  # Standard amino acids
    matrix = np.zeros((20, 20))
    
    # Count occurrences
    total_dipeptides = 0
    single_aa_counts = {aa: 0 for aa in amino_acids}
    
    for seq in protein_sequences:
        protein_obj = Protein(seq, protein_class_name)
        # Count single amino acids
        for aa in seq:
            if aa in amino_acids:
                single_aa_counts[aa] += 1
        
        # Count dipeptides
        for dipeptide, count in protein_obj.two_character_dict.items():
            if len(dipeptide) == 2 and all(aa in amino_acids for aa in dipeptide):
                i = amino_acids.index(dipeptide[0])
                j = amino_acids.index(dipeptide[1])
                matrix[i][j] += count
                total_dipeptides += count
    
    # Convert to frequencies
    matrix = matrix / total_dipeptides if total_dipeptides > 0 else matrix
    
    # Calculate expected frequencies (assuming random distribution)
    individual_freqs = np.array([single_aa_counts[aa] / sum(single_aa_counts.values()) for aa in amino_acids])
    expected_matrix = np.outer(individual_freqs, individual_freqs)
    
    # Calculate log2 ratio of observed/expected
    enrichment_matrix = np.log2(matrix / expected_matrix)
    enrichment_matrix[np.isnan(enrichment_matrix)] = 0
    enrichment_matrix[np.isinf(enrichment_matrix)] = 0
    
    # Create heatmap
    plt.figure(figsize=(12, 10))
    plt.imshow(enrichment_matrix, cmap='RdBu_r', aspect='equal')
    plt.colorbar(label='log2(Observed/Expected)')
    
    # Add labels
    plt.xticks(range(20), list(amino_acids))
    plt.yticks(range(20), list(amino_acids))
    plt.xlabel('Second Amino Acid')
    plt.ylabel('First Amino Acid')
    plt.title(f'Dipeptide Enrichment in {protein_class_name}')
    
    return plt.gcf()

def analyze_enriched_pairs(protein_sequences, protein_class_name, enrichment_threshold=0.5):
    """
    Analyze enriched amino acid pairs and their chemical properties
    enrichment_threshold: log2 value above which pairs are considered enriched 
    (0.5 means 1.4x more frequent than expected, 1.0 means 2x more frequent)
    """
    # Chemical properties of amino acids with clearer naming
    properties = {
        'hydrophobic': 'AILMFWV',
        'polar': 'STNQ',
        'positive': 'RHK',
        'negative': 'DE',
        'structure_breaking': 'GP',  # Changed from 'special'
        'disulfide_bonding': 'C',    # Split from 'special'
        'aromatic': 'FWY'
    }
    
    # Get enrichment matrix
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    matrix = np.zeros((20, 20))
    single_aa_counts = {aa: 0 for aa in amino_acids}
    total_dipeptides = 0
    
    # Count occurrences (same as in heatmap function)
    for seq in protein_sequences:
        protein_obj = Protein(seq, protein_class_name)
        for aa in seq:
            if aa in amino_acids:
                single_aa_counts[aa] += 1
        
        for dipeptide, count in protein_obj.two_character_dict.items():
            if len(dipeptide) == 2 and all(aa in amino_acids for aa in dipeptide):
                i = amino_acids.index(dipeptide[0])
                j = amino_acids.index(dipeptide[1])
                matrix[i][j] += count
                total_dipeptides += count
    
    # Calculate enrichment
    matrix = matrix / total_dipeptides
    individual_freqs = np.array([single_aa_counts[aa] / sum(single_aa_counts.values()) for aa in amino_acids])
    expected_matrix = np.outer(individual_freqs, individual_freqs)
    enrichment_matrix = np.log2(matrix / expected_matrix)
    
    # Find enriched pairs
    enriched_pairs = []
    for i in range(20):
        for j in range(20):
            if enrichment_matrix[i][j] >= enrichment_threshold:
                aa1, aa2 = amino_acids[i], amino_acids[j]
                score = enrichment_matrix[i][j]
                
                # Get properties for each amino acid
                props1 = [prop for prop, aas in properties.items() if aa1 in aas]
                props2 = [prop for prop, aas in properties.items() if aa2 in aas]
                
                enriched_pairs.append({
                    'pair': f"{aa1}{aa2}",
                    'enrichment': score,
                    'properties': f"{props1} -> {props2}"
                })
    
    # Sort by enrichment score
    enriched_pairs.sort(key=lambda x: x['enrichment'], reverse=True)
    
    # Analyze patterns
    property_patterns = {}
    for pair in enriched_pairs:
        props = pair['properties']
        property_patterns[props] = property_patterns.get(props, 0) + 1
    
    # Modified print output to show threshold
    print(f"\nEnriched Dipeptide Analysis for {protein_class_name}")
    print(f"(showing pairs with >{2**enrichment_threshold:.1f}x enrichment)")
    print("=" * 50)
    
    if not enriched_pairs:
        print("\nNo significantly enriched pairs found at this threshold.")
    else:
        print("\nTop 10 Most Enriched Pairs:")
        for pair in enriched_pairs[:10]:
            print(f"{pair['pair']}: {2**pair['enrichment']:.2f}x enrichment ({pair['properties']})")
    
    print("\nCommon Property Patterns:")
    for pattern, count in sorted(property_patterns.items(), key=lambda x: x[1], reverse=True)[:5]:
        print(f"{pattern}: {count} pairs")
    
    return enriched_pairs, property_patterns

def process_protein_files(directory):
    try:
        files = [f for f in os.listdir(directory) if not f.startswith('.')]
    except Exception as e:
        print(f"Error accessing directory {directory}: {str(e)}")
        return
    
    results = {
        'di': {},
        'tri': {},
        'tetra': {}
    }
    
    for file in files:
        protein_class_name = file.split('.')[0]
        sequences = read_protein_sequences_from_file(directory, file)
        
        if not sequences:
            continue
            
        # Create and save heatmap
        fig = create_dipeptide_heatmap(sequences, protein_class_name)
        fig.savefig(f"{protein_class_name}_dipeptide_heatmap.png")
        plt.close(fig)
        
        # Analyze enriched pairs
        analyze_enriched_pairs(sequences, protein_class_name)
            
        results['di'][protein_class_name] = get_top_n_peptides(sequences, protein_class_name, 'two')
        results['tri'][protein_class_name] = get_top_n_peptides(sequences, protein_class_name, 'three')
        results['tetra'][protein_class_name] = get_top_n_peptides(sequences, protein_class_name, 'four')
    
    return results

def extract_sequence_features(sequence):
    """
    Extract numerical features from a protein sequence for ML
    Returns a dictionary of features
    """
    # Basic amino acid properties
    hydrophobicity = {
        'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8,
        'G': -0.4, 'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8,
        'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5,
        'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3
    }
    
    # Feature dictionary
    features = {
        'length': len(sequence),
        'avg_hydrophobicity': sum(hydrophobicity.get(aa, 0) for aa in sequence) / len(sequence),
        'charge_balance': sequence.count('R') + sequence.count('K') - sequence.count('D') - sequence.count('E'),
    }
    
    # Add amino acid frequencies
    for aa in 'ACDEFGHIKLMNPQRSTVWY':
        features[f'freq_{aa}'] = sequence.count(aa) / len(sequence)
    
    # Add dipeptide frequencies (top N most informative based on our previous analysis)
    dipeptides = ['LL', 'AA', 'GG', 'SS', 'EE', 'PP']  # Example - we can expand this
    for di in dipeptides:
        count = sum(1 for i in range(len(sequence)-1) if sequence[i:i+2] == di)
        features[f'dipep_{di}'] = count / (len(sequence)-1)
    
    return features

def prepare_ml_dataset(directory):
    """
    Prepare a dataset for machine learning
    Returns X (features) and y (labels) for ML models
    """
    X = []  # Features
    y = []  # Labels (protein classes)
    
    try:
        files = [f for f in os.listdir(directory) if not f.startswith('.')]
    except Exception as e:
        print(f"Error accessing directory {directory}: {str(e)}")
        return None, None
    
    for file in files:
        protein_class = file.split('.')[0]
        sequences = read_protein_sequences_from_file(directory, file)
        
        for seq in sequences:
            features = extract_sequence_features(seq)
            X.append(features)
            y.append(protein_class)
    
    return X, y

def train_protein_classifier():
    """
    Train a model to classify proteins based on sequence features
    """
    from sklearn.model_selection import train_test_split
    from sklearn.preprocessing import StandardScaler
    from sklearn.ensemble import RandomForestClassifier
    import pandas as pd
    
    # Prepare dataset
    directory = "d:\\Coding\\Wagner_Enterprises\\Uniprot"  # Update path as needed
    X, y = prepare_ml_dataset(directory)
    
    if X is None or y is None:
        return
    
    # Convert to DataFrame for easier handling
    X_df = pd.DataFrame(X)
    
    # Split dataset
    X_train, X_test, y_train, y_test = train_test_split(X_df, y, test_size=0.2, random_state=42)
    
    # Scale features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # Train model
    model = RandomForestClassifier(n_estimators=100, random_state=42)
    model.fit(X_train_scaled, y_train)
    
    # Evaluate
    train_score = model.score(X_train_scaled, y_train)
    test_score = model.score(X_test_scaled, y_test)
    
    print(f"Training accuracy: {train_score:.3f}")
    print(f"Testing accuracy: {test_score:.3f}")
    
    # Feature importance analysis
    feature_importance = pd.DataFrame({
        'feature': X_df.columns,
        'importance': model.feature_importances_
    }).sort_values('importance', ascending=False)
    
    print("\nTop 10 most important features:")
    print(feature_importance.head(10))
    
    return model, scaler, feature_importance

def main():
    # Existing analysis
    directory = "d:\\Coding\\Wagner_Enterprises\\Uniprot"
    results = process_protein_files(directory)
    
    # ML analysis
    model, scaler, feature_importance = train_protein_classifier()
    
    # Print results
    for peptide_type, class_data in results.items():
        print(f"\nTop 20 {peptide_type}peptides by protein class:")
        for protein_class, peptides in class_data.items():
            print(f"{protein_class}: {peptides}")

if __name__ == "__main__":
    main()