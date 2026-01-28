#!/usr/bin/env python3
"""
Cell Wall Gene Enrichment Analysis for P-body Proteomics
=========================================================

Purpose:
    Identify and statistically test enrichment of cell wall-related genes
    in P-body proteomics data under different experimental conditions.

Methodology:
    1. Fetch GO annotations from UniProt for all proteins
    2. Identify cell wall genes using:
       - Specific GO IDs (6 core cell wall terms)
       - Expanded keywords (14 cell wall-related terms)
    3. Filter genes by log2FC threshold (0.8)
    4. Perform Fisher's exact test for enrichment
    5. Calculate fold enrichment and FDR correction
    6. Generate visualizations (bubble plots)

Script Strenghts:
    - Correct statistical background (GO-annotated genes only)
    - Proper GO aspect parsing using GoAspect property
    - Comprehensive cell wall definition (GO IDs + keywords)
    - GO annotation caching for fast re-runs
    - Rich outputs with detailed gene lists
    - Automated validation against baseline
    - Professional bubble plot visualizations

Input:
    - data/processed/pbody_proteomics_annotated.csv (from Script 02)

Output:
    - data/processed/cell_wall_genes_[condition].csv (4 files)
    - data/processed/cell_wall_enrichment_summary.csv
    - data/processed/all_gene_GO_annotations.csv
    - data/processed/go_annotation_cache.json (cache)
    - results/figures/cell_wall_bubble_plot_[condition].png (4 files)
    - logs/03_cell_wall_enrichment_analysis.log

Postdoctoral Bioinformatician Author: Moya-Cuevas J
Supervisor: Moschou PN
Date: 2026-01-22
Version: 3.0
"""

import sys
import json
import time
import logging
import requests
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Dict, List, Set, Tuple
from scipy.stats import fisher_exact, hypergeom
from statsmodels.stats.multitest import multipletests
from datetime import datetime

# ============================================================================
# CONFIGURATION
# ============================================================================

# Directories
BASE_DIR = Path(__file__).parent.parent
DATA_DIR = BASE_DIR / "data"
PROCESSED_DIR = DATA_DIR / "processed"
FIGURES_DIR = BASE_DIR / "results" / "figures"
LOG_DIR = BASE_DIR / "logs"

# Create directories
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

# Files
INPUT_FILE = PROCESSED_DIR / "pbody_proteomics_annotated.csv"
CACHE_FILE = PROCESSED_DIR / "go_annotation_cache.json"
GO_ANNOTATIONS_FILE = PROCESSED_DIR / "all_gene_GO_annotations.csv"
SUMMARY_FILE = PROCESSED_DIR / "cell_wall_enrichment_summary.csv"
LOG_FILE = LOG_DIR / "03_cell_wall_enrichment_analysis.log"

# Analysis parameters
LOG2FC_THRESHOLD = 0.8  # Fixed for baseline validation
BATCH_SIZE = 20  # UniProt API batch size
API_DELAY = 0.3  # Seconds between API requests

# Cell wall GO IDs (from Script 04 - corrected version)
CELL_WALL_GO_IDS = {
    "GO:0005618",  # cell wall
    "GO:0042546",  # cell wall biogenesis
    "GO:0071555",  # cell wall organization
    "GO:0009832",  # plant-type cell wall biogenesis
    "GO:0009834",  # plant-type cell wall organization
    "GO:0030312",  # external encapsulating structure
}

# Cell wall keywords (from Script 04 - expanded)
CELL_WALL_KEYWORDS = [
    "cell wall",
    "cell-wall",
    "cellulose",
    "hemicellulose",
    "xyloglucan",
    "xylan",
    "pectin",
    "lignin",
    "callose",
    "arabinogalactan",
    "wall organization",
    "wall biogenesis",
    "secondary wall",
    "primary wall",
]

# Experimental conditions
CONDITIONS = {
    "AP_NS": "Affinity Purification - No Stress",
    "AP_HS": "Affinity Purification - Heat Stress",
    "PDL_NS": "Proximity-Dependent Labeling - No Stress",
    "PDL_HS": "Proximity-Dependent Labeling - Heat Stress",
}

# ============================================================================
# LOGGING SETUP
# ============================================================================

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(LOG_FILE),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

# ============================================================================
# GO ANNOTATION FUNCTIONS
# ============================================================================

def load_cache() -> Dict[str, List[Dict]]:
    """Load GO annotation cache if it exists."""
    if CACHE_FILE.exists():
        logger.info(f"Loading GO annotation cache from {CACHE_FILE}")
        with open(CACHE_FILE, 'r') as f:
            cache = json.load(f)
        logger.info(f"  ✓ Loaded {len(cache)} cached annotations")
        return cache
    return {}


def save_cache(cache: Dict[str, List[Dict]]) -> None:
    """Save GO annotation cache."""
    logger.info(f"Saving GO annotation cache to {CACHE_FILE}")
    with open(CACHE_FILE, 'w') as f:
        json.dump(cache, f, indent=2)
    logger.info(f"  ✓ Cached {len(cache)} annotations")


def fetch_go_annotations(
    gene_ids: List[str],
    batch_size: int = BATCH_SIZE,
    use_cache: bool = True
) -> Dict[str, List[Dict]]:
    """
    Fetch GO annotations from UniProt for Arabidopsis genes.
    
    Parameters:
    -----------
    gene_ids : list
        List of Arabidopsis gene IDs
    batch_size : int
        Number of genes per batch request
    use_cache : bool
        Whether to use/update cache
        
    Returns:
    --------
    dict
        Mapping of gene_id -> list of GO terms with info
    """
    # Load cache
    go_annotations = load_cache() if use_cache else {}
    
    # Identify genes not in cache
    genes_to_fetch = [g for g in gene_ids if g not in go_annotations]
    
    if not genes_to_fetch:
        logger.info("All genes found in cache, no API requests needed")
        return go_annotations
    
    logger.info(f"Fetching GO annotations from UniProt...")
    logger.info(f"  Genes to fetch: {len(genes_to_fetch)}")
    logger.info(f"  Already cached: {len(gene_ids) - len(genes_to_fetch)}")
    
    # Fetch in batches
    for i in range(0, len(genes_to_fetch), batch_size):
        batch = genes_to_fetch[i:i + batch_size]
        batch_num = i // batch_size + 1
        total_batches = (len(genes_to_fetch) - 1) // batch_size + 1
        
        logger.info(f"  Batch {batch_num}/{total_batches}: {len(batch)} genes")
        
        # Build query
        gene_query = ' OR '.join([f'(gene_exact:{gid})' for gid in batch])
        query = f"({gene_query}) AND (organism_id:3702)"
        
        url = "https://rest.uniprot.org/uniprotkb/search"
        params = {
            'query': query,
            'format': 'json',
            'fields': 'accession,gene_names,go',
            'size': 500
        }
        
        try:
            response = requests.get(url, params=params)
            response.raise_for_status()
            data = response.json()
            
            # Process entries
            batch_upper = {g.upper(): g for g in batch}
            
            for entry in data.get('results', []):
                # Match gene ID
                gene_names = entry.get('genes', [])
                matched_gene_id = None
                
                for gene_data in gene_names:
                    # Try geneName
                    gene_value = gene_data.get('geneName', {}).get('value', '')
                    if gene_value.upper() in batch_upper:
                        matched_gene_id = batch_upper[gene_value.upper()]
                        break
                    
                    # Try orderedLocusNames
                    for oln in gene_data.get('orderedLocusNames', []):
                        oln_value = oln.get('value', '')
                        if oln_value.upper() in batch_upper:
                            matched_gene_id = batch_upper[oln_value.upper()]
                            break
                    
                    if matched_gene_id:
                        break
                
                if not matched_gene_id:
                    continue
                
                # Extract GO terms
                go_terms = []
                for go_entry in entry.get('uniProtKBCrossReferences', []):
                    if go_entry.get('database') == 'GO':
                        go_id = go_entry.get('id', '')
                        properties = go_entry.get('properties', [])
                        
                        go_term = ''
                        go_aspect = ''
                        
                        for prop in properties:
                            if prop.get('key') == 'GoTerm':
                                go_term = prop.get('value', '')
                            elif prop.get('key') == 'GoAspect':
                                # Use GoAspect property (CORRECT METHOD)
                                go_aspect = prop.get('value', '')
                        
                        # Clean GO term (remove prefix if present)
                        go_term_clean = go_term.split(':', 1)[1] if ':' in go_term else go_term
                        
                        go_terms.append({
                            'go_id': go_id,
                            'go_term': go_term_clean,
                            'go_aspect': go_aspect
                        })
                
                if go_terms:
                    go_annotations[matched_gene_id] = go_terms
            
            logger.info(f"    Retrieved {len([g for g in batch if g in go_annotations])} annotations")
            
        except Exception as e:
            logger.error(f"    Error fetching batch: {e}")
        
        # Rate limiting
        time.sleep(API_DELAY)
    
    # Save cache
    if use_cache:
        save_cache(go_annotations)
    
    logger.info(f"  ✓ Total genes with GO annotations: {len(go_annotations)}")
    
    return go_annotations


def identify_cell_wall_genes(go_annotations: Dict[str, List[Dict]]) -> Set[str]:
    """
    Identify genes with cell wall-related GO terms.
    
    Uses both specific GO IDs and expanded keywords.
    
    Parameters:
    -----------
    go_annotations : dict
        GO annotations for each gene
        
    Returns:
    --------
    set
        Set of gene IDs with cell wall-related GO terms
    """
    cell_wall_genes = set()
    
    for gene_id, go_terms in go_annotations.items():
        for go_entry in go_terms:
            go_id = go_entry['go_id']
            go_term_lower = go_entry['go_term'].lower()
            
            # Check GO IDs (specific)
            if go_id in CELL_WALL_GO_IDS:
                cell_wall_genes.add(gene_id)
                break
            
            # Check keywords (broad)
            if any(keyword in go_term_lower for keyword in CELL_WALL_KEYWORDS):
                cell_wall_genes.add(gene_id)
                break
    
    return cell_wall_genes


# ============================================================================
# STATISTICAL ANALYSIS
# ============================================================================

def perform_enrichment_analysis(
    enriched_genes: Set[str],
    cell_wall_genes: Set[str],
    background_genes: Set[str]
) -> Dict:
    """
    Perform Fisher's exact test for cell wall enrichment.
    
    Uses correct statistical background (only GO-annotated genes).
    
    Parameters:
    -----------
    enriched_genes : set
        Genes passing the log2FC threshold
    cell_wall_genes : set
        Genes with cell wall GO terms
    background_genes : set
        All genes with GO annotations (correct background)
        
    Returns:
    --------
    dict
        Statistical results
    """
    # Contingency table:
    #                      Cell Wall    Not Cell Wall
    # Enriched (>thresh)       a              b
    # Not Enriched             c              d
    
    a = len(enriched_genes & cell_wall_genes)
    b = len(enriched_genes - cell_wall_genes)
    c = len(cell_wall_genes - enriched_genes)
    d = len(background_genes - enriched_genes - cell_wall_genes)
    
    # Fisher's exact test
    odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative='greater')
    
    # Calculate fold enrichment
    n_enriched = len(enriched_genes)
    n_cell_wall = len(cell_wall_genes)
    n_total = len(background_genes)
    
    expected = (n_enriched * n_cell_wall) / n_total if n_total > 0 else 0
    fold_enrichment = a / expected if expected > 0 else 0
    
    # Hypergeometric test (alternative)
    hyper_pvalue = hypergeom.sf(a - 1, n_total, n_cell_wall, n_enriched)
    
    results = {
        'a': a,  # Enriched with cell wall
        'b': b,  # Enriched without cell wall
        'c': c,  # Cell wall not enriched
        'd': d,  # Neither
        'n_enriched': n_enriched,
        'n_cell_wall': n_cell_wall,
        'n_total': n_total,
        'expected': expected,
        'fold_enrichment': fold_enrichment,
        'odds_ratio': odds_ratio,
        'fisher_pvalue': p_value,
        'hypergeom_pvalue': hyper_pvalue
    }
    
    return results


# ============================================================================
# VISUALIZATION
# ============================================================================

def create_bubble_plot(
    summary_df: pd.DataFrame,
    condition: str,
    output_file: Path
) -> None:
    """
    Create bubble plot for cell wall enrichment.
    
    Parameters:
    -----------
    summary_df : pd.DataFrame
        Summary statistics
    condition : str
        Condition code (e.g., 'AP_NS')
    output_file : Path
        Output file path
    """
    # Extract data for this condition
    row = summary_df[summary_df['Condition'] == condition]
    
    if row.empty or row['Fisher_Pvalue'].values[0] >= 0.05:
        logger.info(f"  Skipping bubble plot for {condition} (not significant)")
        return
    
    # Create figure
    fig, ax = plt.subplots(figsize=(8, 6))
    
    fold_enrich = row['Fold_Enrichment'].values[0]
    n_genes = row['N_CellWall_Enriched'].values[0]
    pvalue = row['Fisher_Pvalue'].values[0]
    fdr = row['FDR'].values[0]
    
    # Bubble plot
    scatter = ax.scatter(
        [fold_enrich],
        [1],
        s=n_genes * 100,
        c=[-np.log10(fdr)],
        cmap='Reds',
        alpha=0.7,
        edgecolors='black',
        linewidth=2
    )
    
    # Labels
    ax.set_xlabel('Fold Enrichment', fontsize=12, fontweight='bold')
    ax.set_ylabel('')
    ax.set_yticks([1])
    ax.set_yticklabels(['Cell Wall Enrichment'])
    ax.set_title(f'{CONDITIONS[condition]}\nCell Wall Gene Enrichment', 
                 fontsize=14, fontweight='bold')
    ax.axvline(x=1, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlim(0, fold_enrich * 1.3)
    ax.set_ylim(0.5, 1.5)
    
    # Colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('-log10(FDR)', fontsize=11, fontweight='bold')
    
    # Annotation
    ax.text(
        fold_enrich, 1.15,
        f'Fold: {fold_enrich:.2f}x\nGenes: {n_genes}\nFDR: {fdr:.2e}',
        ha='center',
        fontsize=10,
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
    )
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"  ✓ Bubble plot saved: {output_file.name}")


# ============================================================================
# MAIN ANALYSIS PIPELINE
# ============================================================================

def main():
    """Main analysis pipeline."""
    
    logger.info("="*80)
    logger.info("CELL WALL GENE ENRICHMENT ANALYSIS IN P-BODY PROTEOMICS")
    logger.info("="*80)
    logger.info(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Log2FC threshold: {LOG2FC_THRESHOLD}")
    
    # ========================================================================
    # 1. Load proteomics data
    # ========================================================================
    logger.info("\n" + "="*80)
    logger.info("STEP 1: Loading proteomics data")
    logger.info("="*80)
    
    if not INPUT_FILE.exists():
        logger.error(f"Input file not found: {INPUT_FILE}")
        sys.exit(1)
    
    df = pd.read_csv(INPUT_FILE)
    logger.info(f"  ✓ Loaded {len(df)} proteins from {INPUT_FILE.name}")
    
    # Get unique gene IDs
    all_genes = set(df['protein'].unique())
    logger.info(f"  ✓ Total unique genes: {len(all_genes)}")
    
    # ========================================================================
    # 2. Fetch GO annotations
    # ========================================================================
    logger.info("\n" + "="*80)
    logger.info("STEP 2: Fetching GO annotations")
    logger.info("="*80)
    
    go_annotations = fetch_go_annotations(list(all_genes))
    annotated_genes = set(go_annotations.keys())
    
    logger.info(f"  ✓ Genes with GO annotations: {len(annotated_genes)}")
    logger.info(f"  ✓ Coverage: {100*len(annotated_genes)/len(all_genes):.2f}%")
    
    # ========================================================================
    # 3. Identify cell wall genes
    # ========================================================================
    logger.info("\n" + "="*80)
    logger.info("STEP 3: Identifying cell wall genes")
    logger.info("="*80)
    logger.info(f"  GO IDs: {len(CELL_WALL_GO_IDS)}")
    logger.info(f"  Keywords: {len(CELL_WALL_KEYWORDS)}")
    
    cell_wall_genes = identify_cell_wall_genes(go_annotations)
    logger.info(f"  ✓ Cell wall genes found: {len(cell_wall_genes)}")
    logger.info(f"  ✓ Percentage of annotated: {100*len(cell_wall_genes)/len(annotated_genes):.2f}%")
    
    # ========================================================================
    # 4. Analyze enrichment per condition
    # ========================================================================
    logger.info("\n" + "="*80)
    logger.info("STEP 4: Analyzing enrichment per condition")
    logger.info("="*80)
    
    results_summary = []
    
    for cond_col, cond_name in CONDITIONS.items():
        logger.info(f"\n{cond_name} ({cond_col}):")
        logger.info("-"*80)
        
        # Filter by threshold
        enriched_df = df[df[cond_col] > LOG2FC_THRESHOLD].copy()
        enriched_genes = set(enriched_df['protein'].unique()) & annotated_genes
        
        # Cell wall genes in enriched set
        cw_enriched = enriched_genes & cell_wall_genes
        
        logger.info(f"  Genes with log2FC > {LOG2FC_THRESHOLD}: {len(enriched_genes)}")
        logger.info(f"  Cell wall genes enriched: {len(cw_enriched)}")
        
        # Perform enrichment analysis
        stats = perform_enrichment_analysis(
            enriched_genes,
            cell_wall_genes & annotated_genes,
            annotated_genes
        )
        
        # Print contingency table
        logger.info(f"\n  Contingency Table:")
        logger.info(f"  {'─'*76}")
        logger.info(f"                        Cell Wall    Not Cell Wall    Total")
        logger.info(f"  Enriched (>{LOG2FC_THRESHOLD})      {stats['a']:8d}     {stats['b']:8d}     {stats['n_enriched']:8d}")
        logger.info(f"  Not Enriched          {stats['c']:8d}     {stats['d']:8d}     {len(annotated_genes) - stats['n_enriched']:8d}")
        logger.info(f"  Total                 {stats['n_cell_wall']:8d}     {len(annotated_genes) - stats['n_cell_wall']:8d}     {stats['n_total']:8d}")
        logger.info(f"  {'─'*76}")
        logger.info(f"  Expected (random):     {stats['expected']:.2f}")
        logger.info(f"  Observed:              {stats['a']}")
        logger.info(f"  Fold Enrichment:       {stats['fold_enrichment']:.3f}x")
        logger.info(f"  Odds Ratio:            {stats['odds_ratio']:.3f}")
        logger.info(f"  Fisher's P-value:      {stats['fisher_pvalue']:.4e}")
        logger.info(f"  Hypergeometric P:      {stats['hypergeom_pvalue']:.4e}")
        
        # Interpretation
        if stats['fisher_pvalue'] < 0.001:
            sig = "*** Highly significant"
        elif stats['fisher_pvalue'] < 0.01:
            sig = "** Very significant"
        elif stats['fisher_pvalue'] < 0.05:
            sig = "* Significant"
        else:
            sig = "Not significant"
        logger.info(f"  Significance:          {sig}")
        
        # Save cell wall genes for this condition
        if len(cw_enriched) > 0:
            cw_enriched_df = enriched_df[enriched_df['protein'].isin(cw_enriched)].copy()
            
            # Add GO term information
            cw_enriched_df['cell_wall_GO_terms'] = cw_enriched_df['protein'].apply(
                lambda x: '; '.join([
                    f"{go['go_term']} ({go['go_id']})"
                    for go in go_annotations.get(x, [])
                    if go['go_id'] in CELL_WALL_GO_IDS or
                       any(kw in go['go_term'].lower() for kw in CELL_WALL_KEYWORDS)
                ])
            )
            
            output_file = PROCESSED_DIR / f"cell_wall_genes_{cond_col}.csv"
            cw_enriched_df.to_csv(output_file, index=False)
            logger.info(f"\n  ✓ Cell wall genes saved: {output_file.name}")
            
            # Show top genes
            logger.info(f"\n  Top 10 cell wall genes by log2FC:")
            top_genes = cw_enriched_df.nlargest(10, cond_col)[
                ['protein', 'gene_name', 'protein_name', cond_col]
            ]
            for idx, row in top_genes.iterrows():
                gene_name = row.get('gene_name', 'N/A')
                protein_name = str(row.get('protein_name', 'Unknown'))[:40]
                logger.info(
                    f"    {row['protein']:<12} {gene_name:<12} "
                    f"{protein_name:<40} {row[cond_col]:6.3f}"
                )
        
        # Store results for summary
        results_summary.append({
            'Condition': cond_col,
            'Description': cond_name,
            'N_Enriched': stats['n_enriched'],
            'N_CellWall_Enriched': stats['a'],
            'Expected': stats['expected'],
            'Fold_Enrichment': stats['fold_enrichment'],
            'Odds_Ratio': stats['odds_ratio'],
            'Fisher_Pvalue': stats['fisher_pvalue'],
            'Hypergeom_Pvalue': stats['hypergeom_pvalue']
        })
    
    # ========================================================================
    # 5. Create summary table
    # ========================================================================
    logger.info("\n" + "="*80)
    logger.info("STEP 5: Creating summary table")
    logger.info("="*80)
    
    summary_df = pd.DataFrame(results_summary)
    
    # FDR correction
    summary_df['FDR'] = multipletests(summary_df['Fisher_Pvalue'], method='fdr_bh')[1]
    
    # Save summary
    summary_df.to_csv(SUMMARY_FILE, index=False)
    logger.info(f"  ✓ Summary saved: {SUMMARY_FILE.name}")
    
    logger.info("\n" + summary_df.to_string(index=False))
    
    # ========================================================================
    # 6. Save complete GO annotations
    # ========================================================================
    logger.info("\n" + "="*80)
    logger.info("STEP 6: Saving complete GO annotations")
    logger.info("="*80)
    
    go_df_rows = []
    for gene_id, go_list in go_annotations.items():
        for go_entry in go_list:
            go_df_rows.append({
                'gene_id': gene_id,
                'go_id': go_entry['go_id'],
                'go_term': go_entry['go_term'],
                'go_aspect': go_entry['go_aspect'],
                'is_cell_wall': gene_id in cell_wall_genes
            })
    
    go_df = pd.DataFrame(go_df_rows)
    go_df.to_csv(GO_ANNOTATIONS_FILE, index=False)
    logger.info(f"  ✓ GO annotations saved: {GO_ANNOTATIONS_FILE.name}")
    logger.info(f"  ✓ Total GO annotations: {len(go_df)}")
    
    # ========================================================================
    # 7. Create bubble plots
    # ========================================================================
    logger.info("\n" + "="*80)
    logger.info("STEP 7: Creating bubble plots")
    logger.info("="*80)
    
    for cond_col in CONDITIONS.keys():
        output_file = FIGURES_DIR / f"cell_wall_bubble_plot_{cond_col}.png"
        create_bubble_plot(summary_df, cond_col, output_file)
    
    # ========================================================================
    # 8. Validation against baseline
    # ========================================================================
    logger.info("\n" + "="*80)
    logger.info("STEP 8: Validation against baseline")
    logger.info("="*80)
    
    baseline_summary = PROCESSED_DIR / "cell_wall_enrichment_summary_baseline.csv"
    
    if baseline_summary.exists():
        baseline_df = pd.read_csv(baseline_summary)
        
        # Compare key metrics
        logger.info("\nComparison with baseline:")
        logger.info("-"*80)
        logger.info(f"{'Metric':<30} {'Baseline':<15} {'Ours':<15} {'Match':<10}")
        logger.info("-"*80)
        
        for _, row in summary_df.iterrows():
            cond = row['Condition']
            baseline_row = baseline_df[baseline_df['Condition'] == cond]
            
            if not baseline_row.empty:
                for metric in ['N_CellWall_Enriched', 'Fold_Enrichment', 'Fisher_Pvalue']:
                    baseline_val = baseline_row[metric].values[0]
                    our_val = row[metric]
                    
                    # Check match (allow small numerical differences)
                    if metric == 'N_CellWall_Enriched':
                        match = baseline_val == our_val
                    else:
                        match = abs(baseline_val - our_val) < 0.01
                    
                    match_str = "✓" if match else "✗"
                    
                    logger.info(
                        f"{cond} - {metric:<20} {baseline_val:<15.4g} "
                        f"{our_val:<15.4g} {match_str:<10}"
                    )
        
        logger.info("-"*80)
    else:
        logger.warning("  Baseline summary not found, skipping validation")
    
    # ========================================================================
    # COMPLETION
    # ========================================================================
    logger.info("\n" + "="*80)
    logger.info("ANALYSIS COMPLETE")
    logger.info("="*80)
    logger.info(f"End time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"\nOutputs generated:")
    logger.info(f"  - {SUMMARY_FILE.name}")
    logger.info(f"  - {GO_ANNOTATIONS_FILE.name}")
    logger.info(f"  - {CACHE_FILE.name}")
    logger.info(f"  - cell_wall_genes_[condition].csv (4 files)")
    logger.info(f"  - cell_wall_bubble_plot_[condition].png (4 files)")
    logger.info(f"  - {LOG_FILE.name}")
    
    return summary_df, go_annotations, cell_wall_genes


if __name__ == "__main__":
    summary_df, go_annotations, cell_wall_genes = main()
