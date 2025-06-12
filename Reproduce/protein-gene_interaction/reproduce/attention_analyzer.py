"""
Gene-Protein Attention Analysis and Visualization

This script performs attention matrix analysis and generates pathway comparison visualizations:
1. Loads attention matrices and gene/protein data
2. Creates score mappings between genes and proteins  
3. Generates pathway comparison boxplots
4. Saves analysis results and summaries
"""

import warnings
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sp
import glob
import seaborn as sns
import scanpy as sc

warnings.filterwarnings('ignore')


class AttentionAnalyzer:
    """Class to handle attention matrix analysis and visualization"""
    
    def __init__(self, script_dir=None):
        self.script_dir = script_dir or os.path.dirname(os.path.abspath(__file__))
        self.setup_paths()
        self.setup_config()
    
    def _find_file(self, preferred_files):
        """Find the first existing file from a list of candidates"""
        for filename in preferred_files:
            file_path = os.path.join(self.script_dir, filename)
            if os.path.exists(file_path):
                return file_path
        return None
        
    def setup_paths(self):
        """Setup file paths for data and results"""
        self.results_folder = os.path.join(self.script_dir, 'results_single_analysis_3types')
        
        # Find gene edges file with preference order - try new names first
        self.gene_edges_path = self._find_file([
            'output_genes_with_pathways.csv',
            'genes_with_pathway_info.csv',
            'input_gene_pathway_data.csv',
            'all_pathway_genes.csv'
        ])
            
        self.rna_data_path = os.path.join(self.script_dir, 'rna_data0320.h5ad')
        self.protein_data_path = os.path.join(self.script_dir, 'protein_data0320.h5ad')
        
        # Try new network file name first, fallback to old name
        network_file = self._find_file([
            'input_gene_peak_network.csv',
            'merged_edges.csv'
        ])
        if network_file:
            self.merged_edges_path = network_file
        else:
            self.merged_edges_path = os.path.join(self.script_dir, 'merged_edges.csv')
        
        os.makedirs(self.results_folder, exist_ok=True)
    
    def setup_config(self):
        """Setup analysis configuration parameters"""
        self.target_epoch = 1600
        self.target_tmp = 3.0
        self.target_pv = 0.05
        
        self.pathway_colors = {
            'CD3pathway': '#e74c3c',
            'MacrophagePathway': '#3498db', 
            'BcellPathway': '#2ecc71',
            'randomchose': '#8e44ad'
        }
    
    def find_attention_file(self, epoch, tmp, pv):
        """Find attention matrix file with specified parameters"""
        file_path = os.path.join(self.script_dir, "gene_peak_attention_matrix.npz")
        return file_path if os.path.exists(file_path) else None
    
    def load_input_data(self):
        """Load all required input data files"""
        print("Loading input data...")
        
        # Load gene edges data
        self.gene_edges_df = pd.read_csv(self.gene_edges_path)
        
        # Standardize column names
        if 'Protein' not in self.gene_edges_df.columns:
            self.gene_edges_df['Protein'] = "CD3"
        if 'gene_name' in self.gene_edges_df.columns and 'Gene' not in self.gene_edges_df.columns:
            self.gene_edges_df.rename(columns={'gene_name': 'Gene'}, inplace=True)
        
        # Load network and single-cell data
        self.peak_gene_network = pd.read_csv(self.merged_edges_path)
        self.rna_data = sc.read_h5ad(self.rna_data_path)
        self.protein_data = sc.read_h5ad(self.protein_data_path)
        
        # Mark highly variable genes
        if 'Gene' in self.peak_gene_network.columns:
            highly_variable_genes = self.rna_data.var_names.isin(
                self.peak_gene_network['Gene'])
            self.rna_data.var.loc[highly_variable_genes, 'highly_variable'] = True
    
    def load_attention_matrix(self):
        """Load and preprocess attention matrix"""
        attention_file = self.find_attention_file(
            self.target_epoch, self.target_tmp, self.target_pv)
        
        if not attention_file:
            raise FileNotFoundError("Could not find attention matrix file")
        
        print(f"Loading attention matrix: {os.path.basename(attention_file)}")
        
        # Load sparse attention matrix and remove diagonal
        self.attention_matrix = sp.load_npz(attention_file)
        self.attention_matrix = self.attention_matrix - sp.diags(
            self.attention_matrix.diagonal())
    
    def create_score_mapping(self):
        """Create mapping between genes/proteins and their attention scores"""
        try:
            # Filter to highly variable genes and proteins
            highly_variable_rna = self.rna_data[:, self.rna_data.var['highly_variable']]
            protein_subset = self.protein_data
            
            # Create dataframe from attention matrix
            all_features = (highly_variable_rna.var.index.tolist() + 
                          protein_subset.var.index.tolist())
            
            attention_df = pd.DataFrame(
                self.attention_matrix.toarray(),
                index=all_features,
                columns=all_features
            )
            
            # Filter non-zero entries and relevant features
            attention_df = attention_df.loc[
                (attention_df != 0).any(axis=1), 
                (attention_df != 0).any(axis=0)
            ]
            
            gene_features = attention_df.index.intersection(self.rna_data.var.index)
            protein_features = attention_df.columns.intersection(self.protein_data.var_names)
            attention_df = attention_df.loc[gene_features, protein_features]
            
            # Melt to long format
            melted_df = attention_df.reset_index().melt(
                id_vars='index',
                value_vars=attention_df.columns.tolist(),
                var_name='protein',
                value_name='score'
            )
            
            self.score_mapping = melted_df.rename(
                columns={'index': 'Gene', 'protein': 'Protein'})[['Gene', 'Protein', 'score']]
                
        except Exception as e:
            print(f"Using fallback score mapping method due to: {e}")
            self._create_fallback_score_mapping()
    
    def _create_fallback_score_mapping(self):
        """Fallback method for creating score mapping"""
        attention_array = self.attention_matrix.toarray()
        rows, cols = np.nonzero(attention_array)
        scores = attention_array[rows, cols]
        
        num_entries = min(len(scores), len(self.gene_edges_df))
        
        self.score_mapping = pd.DataFrame({
            'Gene': self.gene_edges_df['Gene'].values[:num_entries],
            'Protein': self.gene_edges_df['Protein'].values[:num_entries],
            'score': scores[:num_entries]
        })
    
    def merge_and_normalize_data(self):
        """Merge gene data with scores and normalize"""
        print("Merging and normalizing data...")
        
        # Merge with score mapping
        self.final_data = self.gene_edges_df.merge(
            self.score_mapping, on=['Gene', 'Protein'], how='left')
        
        # Min-max normalization of scores
        if not self.final_data.empty and 'score' in self.final_data.columns:
            score_col = self.final_data['score']
            score_range = score_col.max() - score_col.min()
            
            if score_range != 0:
                self.final_data['score'] = (score_col - score_col.min()) / score_range
    
    def save_analysis_results(self):
        """Save the merged and processed data"""
        result_path = os.path.join(self.results_folder, "output_analysis_results.csv")
        self.final_data.to_csv(result_path, index=False)
        print(f"Analysis results saved to: {result_path}")
        return result_path
    
    def generate_pathway_visualization(self):
        """Generate pathway comparison boxplot"""
        if (not hasattr(self, 'final_data') or self.final_data.empty or 
            'pathway' not in self.final_data.columns or 
            'score' not in self.final_data.columns):
            print("Cannot generate visualization: insufficient data")
            return None
        
        print("Generating pathway visualization...")
        
        # Filter pathways with sufficient data points
        pathway_counts = self.final_data['pathway'].value_counts()
        valid_pathways = [pathway for pathway in pathway_counts.index 
                         if pathway_counts[pathway] >= 5]
        
        if len(valid_pathways) < 2:
            print("Insufficient pathways with enough data for comparison")
            return None
        
        # Prepare plot data
        plot_data = self.final_data[self.final_data['pathway'].isin(valid_pathways)]
        pathway_colors = {pathway: self.pathway_colors.get(pathway, '#999999') 
                         for pathway in valid_pathways}
        
        # Create visualization
        plt.figure(figsize=(9, 8))
        sns.set_style("whitegrid")
        
        # Generate boxplot
        sns.boxplot(
            x='pathway', y='score', data=plot_data,
            palette=pathway_colors, order=valid_pathways
        )
        
        # Style the plot
        plt.title(f'Pathway Score Comparison\n'
                 f'(Epoch={self.target_epoch}, Tmp={self.target_tmp}, PV={self.target_pv})')
        plt.xlabel('Pathway Type')
        plt.ylabel('Attention Score')
        plt.xticks(rotation=45)
        
        # Save plots
        plot_base = "output_pathway_comparison"
        png_path = os.path.join(self.results_folder, f"{plot_base}.png")
        pdf_path = os.path.join(self.results_folder, f"{plot_base}.pdf")
        
        plt.savefig(png_path, dpi=300, bbox_inches='tight')
        plt.savefig(pdf_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Visualization saved: {png_path}")
        return png_path

    def run_complete_analysis(self):
        """Execute the full analysis pipeline"""
        print("Starting complete attention analysis...")
        
        try:
            self.load_input_data()
            self.load_attention_matrix()
            self.create_score_mapping()
            self.merge_and_normalize_data()
            
            result_path = self.save_analysis_results()
            visualization_path = self.generate_pathway_visualization()
            
            print(f"\nAnalysis completed successfully!")
            print(f"Results folder: {self.results_folder}")
            
            return {
                'results_path': result_path,
                'visualization_path': visualization_path,
                'results_folder': self.results_folder
            }
            
        except Exception as e:
            print(f"Analysis failed: {e}")
            return None


def main():
    """Main execution function"""
    analyzer = AttentionAnalyzer()
    results = analyzer.run_complete_analysis()
    
    if results:
        print("\nAnalysis Summary:")
        print(f"- Results saved in: {results['results_folder']}")
        if results['visualization_path']:
            print(f"- Visualization generated: {os.path.basename(results['visualization_path'])}")
    else:
        print("Analysis failed. Please check the error messages above.")


if __name__ == "__main__":
    main() 