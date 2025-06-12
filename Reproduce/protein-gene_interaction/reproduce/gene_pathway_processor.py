"""
Gene-Pathway Analysis and Visualization Data Preparation

This script processes gene pathway data by:
1. Reading pathway reference files and gene data
2. Creating gene-to-pathway mappings
3. Adding pathway names and attention scores to genes
4. Generating final output for visualization
"""

import pandas as pd
import os


class GenePathwayProcessor:
    """Class to handle gene-pathway mapping and score assignment"""
    
    def __init__(self, base_directory="."):
        self.base_dir = base_directory
        self.pathway_mappings = {}
        self.gene_to_score = {}
    
    def _find_file(self, preferred_files):
        """
        Find the first existing file from a list of candidates
        
        Args:
            preferred_files: List of filenames to check in order of preference
            
        Returns:
            str: Path to first existing file or None if none exist
        """
        for filename in preferred_files:
            file_path = os.path.join(self.base_dir, filename)
            if os.path.exists(file_path):
                return file_path
        return None
    
    def load_data_files(self):
        """Load all required CSV files"""
        print("Loading data files...")
        
        # Main gene data - try new name first, fallback to old name
        gene_data_file = self._find_file([
            'input_gene_pathway_data.csv',
            'all_pathway_genes.csv'
        ])
        if not gene_data_file:
            raise FileNotFoundError("No gene pathway data file found")
        self.all_genes_df = pd.read_csv(gene_data_file)
        
        # Pathway reference files - try new names first, fallback to old names
        pathway_files = {
            'cd3': self._find_file([
                'input_cd3_pathway_genes.csv',
                'pathway_gene_summary.csv'
            ]),
            'bcell': self._find_file([
                'input_bcell_pathway_genes.csv',
                'bcr_pathway_summary.csv'
            ]),
            'macrophage': self._find_file([
                'input_macrophage_pathway_genes.csv',
                'macrophage_pathway_gene_summary.csv'
            ])
        }
        
        self.pathway_dfs = {}
        for pathway_type, file_path in pathway_files.items():
            if not file_path:
                raise FileNotFoundError(f"No {pathway_type} pathway file found")
            self.pathway_dfs[pathway_type] = pd.read_csv(file_path)
        
        # Gene scores - try new name first, fallback to old names
        score_file = self._find_file([
            'input_attention_scores.csv',
            'gene_protein_attention_scores.csv',
            'result_df_epoch1600_tmp3.0_pv0.05.csv'
        ])
        
        if not score_file:
            raise FileNotFoundError("No score file found")
            
        score_df = pd.read_csv(score_file)
        self.gene_to_score = dict(zip(score_df['Gene'], score_df['score']))
    
    def create_pathway_mapping(self, pathway_df):
        """
        Create gene-to-pathway mapping from pathway dataframe
        
        Args:
            pathway_df: DataFrame with 'Pathway' and 'Genes' columns
            
        Returns:
            dict: Mapping from gene names (lowercase) to pathway names
        """
        gene_to_pathway = {}
        
        for _, row in pathway_df.iterrows():
            pathway_name = row['Pathway']
            genes_string = row['Genes']
            
            # Parse comma-separated gene list
            gene_list = [gene.strip() for gene in genes_string.split(',')]
            
            for gene in gene_list:
                if gene:  # Skip empty strings
                    gene_to_pathway[gene.lower()] = pathway_name
        
        return gene_to_pathway
    
    def build_all_pathway_mappings(self):
        """Build pathway mappings for all pathway types"""
        print("Building pathway mappings...")
        
        pathway_mapping_configs = {
            'CD3pathway': 'cd3',
            'BcellPathway': 'bcell', 
            'MacrophagePathway': 'macrophage'
        }
        
        for pathway_key, df_key in pathway_mapping_configs.items():
            self.pathway_mappings[pathway_key] = self.create_pathway_mapping(
                self.pathway_dfs[df_key])
    
    def find_pathway_for_gene(self, gene_name, pathway_type):
        """
        Find pathway name for a gene with flexible matching
        
        Args:
            gene_name: Name of the gene to search
            pathway_type: Type of pathway to search in
            
        Returns:
            str: Pathway name or 'Unknown' if not found
        """
        if pathway_type not in self.pathway_mappings:
            return "Unknown"
        
        mapping = self.pathway_mappings[pathway_type]
        gene_lower = gene_name.lower()
        
        # Direct match
        if gene_lower in mapping:
            return mapping[gene_lower]
        
        # Substring matching
        for known_gene, pathway in mapping.items():
            if gene_lower in known_gene or known_gene in gene_lower:
                return pathway
        
        return "Unknown"
    
    def assign_pathway_names(self):
        """Assign pathway names to all genes in the dataset"""
        print("Assigning pathway names...")
        
        # Initialize pathway name column
        self.all_genes_df['pathwayname'] = ''
        
        for idx, row in self.all_genes_df.iterrows():
            gene_name = row['gene_name']
            pathway_type = row['pathway']
            
            if pathway_type in self.pathway_mappings:
                pathway_name = self.find_pathway_for_gene(gene_name, pathway_type)
            elif pathway_type == 'randomchose':
                pathway_name = 'Randomchoose'
            else:
                pathway_name = 'Unknown'
            
            self.all_genes_df.at[idx, 'pathwayname'] = pathway_name
    
    def add_metadata_columns(self):
        """Add protein and attention score columns"""
        print("Adding metadata columns...")
        
        # Add protein column (all CD3 for this dataset)
        self.all_genes_df['Protein'] = 'CD3'
        
        # Add attention scores
        self.all_genes_df['Attentionscore'] = self.all_genes_df['gene_name'].map(
            self.gene_to_score)
    
    def save_results(self, output_filename='output_genes_with_pathways.csv'):
        """Save processed data to CSV file"""
        output_path = os.path.join(self.base_dir, output_filename)
        self.all_genes_df.to_csv(output_path, index=False)
        print(f"Results saved to: {output_path}")
        return output_path
    
    def process_all(self):
        """Run the complete processing pipeline"""
        self.load_data_files()
        self.build_all_pathway_mappings()
        self.assign_pathway_names()
        self.add_metadata_columns()
        return self.save_results()


def main():
    """Main execution function"""
    processor = GenePathwayProcessor()
    output_file = processor.process_all()
    print(f"Processing complete. Updated file saved as '{output_file}'")


if __name__ == "__main__":
    main()
