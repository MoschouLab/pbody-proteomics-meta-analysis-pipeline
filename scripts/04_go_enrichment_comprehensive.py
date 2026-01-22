#!/usr/bin/env python3
"""
Comprehensive GO Term Enrichment Analysis for P-body Proteomics
================================================================

Purpose:
    Perform genome-wide GO enrichment analysis across all biological processes,
    molecular functions, and cellular components. Unlike Script 03 (cell wall
    specific), this script identifies ALL significantly enriched GO terms.

Features:
    - Generic analysis: ALL GO terms tested
    - Separation by GO aspect: BP, MF, CC
    - g:Profiler validation (gold-standard tool)
    - CLI parametrizable for flexibility
    - Bubble plot visualizations per aspect
    - Reuses GO annotation cache from Script 03

Input:
    - data/processed/pbody_proteomics_annotated.csv (from Script 02)
    - data/processed/go_annotation_cache.json (from Script 03)

Output:
    - data/processed/go_enrichment_[condition]_[aspect].csv (12 files)
    - data/processed/go_enrichment_summary_all_conditions.csv
    - results/figures/go_bubble_[condition]_[aspect].png (12 plots)
    - results/gprofiler/ (optional, if g:Profiler enabled)
    - logs/04_go_enrichment_comprehensive.log

Author: Moschou PN and Moya-Cuevas J
Date: 2026-01-22
Version: 4.0 (comprehensive, CLI, g:Profiler)
"""

import sys
import json
import logging
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Dict, List, Set
from collections import defaultdict, Counter
from scipy.stats import fisher_exact
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
GPROFILER_DIR = BASE_DIR / "results" / "gprofiler"
LOG_DIR = BASE_DIR / "logs"

# Create directories
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
GPROFILER_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

# Files
INPUT_FILE = PROCESSED_DIR / "pbody_proteomics_annotated.csv"
CACHE_FILE = PROCESSED_DIR / "go_annotation_cache.json"
LOG_FILE = LOG_DIR / "04_go_enrichment_comprehensive.log"

# GO Aspects
GO_ASPECTS = {
    'BP': 'Biological Process',
    'MF': 'Molecular Function',
    'CC': 'Cellular Component'
}

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

def setup_logging(verbose: bool = False):
    """Setup logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(LOG_FILE),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)

# ============================================================================
# GO TERM ANALYSIS FUNCTIONS
# ============================================================================

def load_go_cache() -> Dict[str, List[Dict]]:
    """Load GO annotation cache from Script 03."""
    if not CACHE_FILE.exists():
        raise FileNotFoundError(
            f"GO cache not found: {CACHE_FILE}\n"
            "Please run Script 03 first to generate GO annotations."
        )
    
    with open(CACHE_FILE, 'r') as f:
        cache = json.load(f)
    
    return cache


def build_go_term_gene_sets(
    go_annotations: Dict[str, List[Dict]],
    min_genes: int = 3
) -> Dict[str, Dict]:
    """
    Build gene sets for each GO term.
    
    Returns dict: go_id -> {genes, term, aspect}
    """
    go_term_genes = defaultdict(lambda: {
        'genes': set(),
        'term': '',
        'aspect': ''
    })
    
    for gene_id, go_list in go_annotations.items():
        for go_entry in go_list:
            go_id = go_entry['go_id']
            go_term_genes[go_id]['genes'].add(gene_id)
            go_term_genes[go_id]['term'] = go_entry['go_term']
            go_term_genes[go_id]['aspect'] = go_entry['go_aspect']
    
    # Filter by minimum genes
    filtered = {
        k: v for k, v in go_term_genes.items()
        if len(v['genes']) >= min_genes
    }
    
    return filtered


def perform_go_enrichment(
    enriched_genes: Set[str],
    go_term_genes: Dict[str, Dict],
    background_genes: Set[str],
    aspect: str = None,
    min_genes_in_term: int = 2
) -> pd.DataFrame:
    """
    Perform GO enrichment analysis using Fisher's exact test.
    
    Parameters:
    -----------
    enriched_genes : set
        Genes passing threshold
    go_term_genes : dict
        GO term -> gene sets
    background_genes : set
        All genes with GO annotations
    aspect : str
        Filter by GO aspect (BP/MF/CC) or None for all
    min_genes_in_term : int
        Minimum genes in GO term to test
        
    Returns:
    --------
    pd.DataFrame
        Enrichment results sorted by p-value
    """
    results = []
    
    # Filter by aspect if specified
    if aspect:
        filtered_terms = {
            k: v for k, v in go_term_genes.items()
            if v['aspect'] == aspect
        }
    else:
        filtered_terms = go_term_genes
    
    for go_id, go_info in filtered_terms.items():
        go_genes = go_info['genes'] & background_genes
        
        # Overlap with enriched genes
        a = len(enriched_genes & go_genes)
        
        if a < min_genes_in_term:
            continue
        
        b = len(enriched_genes - go_genes)
        c = len(go_genes - enriched_genes)
        d = len(background_genes - enriched_genes - go_genes)
        
        try:
            odds_ratio, p_value = fisher_exact(
                [[a, b], [c, d]],
                alternative='greater'
            )
        except:
            continue
        
        # Calculate metrics
        n_enriched = len(enriched_genes)
        n_go = len(go_genes)
        n_total = len(background_genes)
        
        expected = (n_enriched * n_go) / n_total if n_total > 0 else 0
        fold_enrichment = a / expected if expected > 0 else 0
        
        results.append({
            'GO_ID': go_id,
            'GO_Term': go_info['term'],
            'GO_Aspect': go_info['aspect'],
            'Genes_in_GO': n_go,
            'Genes_Enriched': a,
            'Expected': expected,
            'Fold_Enrichment': fold_enrichment,
            'Odds_Ratio': odds_ratio,
            'P_value': p_value
        })
    
    if not results:
        return pd.DataFrame()
    
    df = pd.DataFrame(results)
    df['FDR'] = multipletests(df['P_value'], method='fdr_bh')[1]
    
    return df.sort_values('P_value')


# ============================================================================
# g:PROFILER INTEGRATION
# ============================================================================

def run_gprofiler(
    gene_list: List[str],
    condition: str,
    logger
) -> pd.DataFrame:
    """
    Run g:Profiler enrichment analysis.
    
    Parameters:
    -----------
    gene_list : list
        List of gene names
    condition : str
        Condition identifier
    logger : logging.Logger
        Logger instance
        
    Returns:
    --------
    pd.DataFrame
        g:Profiler results
    """
    try:
        from gprofiler import GProfiler
    except ImportError:
        logger.warning("g:Profiler not installed, skipping")
        return pd.DataFrame()
    
    if not gene_list:
        logger.warning(f"Empty gene list for {condition}, skipping g:Profiler")
        return pd.DataFrame()
    
    logger.info(f"  Running g:Profiler for {condition} ({len(gene_list)} genes)")
    
    try:
        gp = GProfiler(return_dataframe=True)
        
        results = gp.profile(
            organism="athaliana",
            query=gene_list,
            sources=["GO:BP", "GO:MF", "GO:CC"],
            user_threshold=0.05,
            significance_threshold_method="fdr"
        )
        
        if results is not None and not results.empty:
            # Save results
            output_file = GPROFILER_DIR / f"gprofiler_{condition}.csv"
            results.to_csv(output_file, index=False)
            logger.info(f"    ✓ g:Profiler results saved: {output_file.name}")
            
            # Separate by aspect
            for aspect in ["BP", "MF", "CC"]:
                subset = results[results["source"] == f"GO:{aspect}"]
                if not subset.empty:
                    aspect_file = GPROFILER_DIR / f"gprofiler_{condition}_{aspect}.csv"
                    subset.to_csv(aspect_file, index=False)
            
            return results
        else:
            logger.warning(f"    No g:Profiler results for {condition}")
            return pd.DataFrame()
            
    except Exception as e:
        logger.error(f"    g:Profiler error: {e}")
        return pd.DataFrame()


# ============================================================================
# VISUALIZATION
# ============================================================================

def create_bubble_plot(
    enrichment_df: pd.DataFrame,
    condition: str,
    aspect: str,
    output_file: Path,
    top_n: int = 20
) -> None:
    """Create bubble plot for GO enrichment."""
    if enrichment_df.empty:
        return
    
    # Filter significant results
    sig_df = enrichment_df[enrichment_df['FDR'] < 0.05]
    
    if sig_df.empty:
        return
    
    plot_df = sig_df.head(top_n).copy()
    plot_df['-log10(FDR)'] = -np.log10(plot_df['FDR'].clip(lower=1e-10))
    plot_df['GO_Term_Short'] = plot_df['GO_Term'].apply(
        lambda x: x[:50] + '...' if len(x) > 50 else x
    )
    plot_df = plot_df.iloc[::-1]  # Reverse for plotting
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, max(6, top_n * 0.4)))
    
    scatter = ax.scatter(
        plot_df['Fold_Enrichment'],
        range(len(plot_df)),
        s=plot_df['Genes_Enriched'] * 50,
        c=plot_df['-log10(FDR)'],
        cmap='Reds',
        alpha=0.7,
        edgecolors='black',
        linewidth=1
    )
    
    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels(plot_df['GO_Term_Short'], fontsize=9)
    ax.set_xlabel('Fold Enrichment', fontsize=12, fontweight='bold')
    ax.set_ylabel('GO Term', fontsize=12, fontweight='bold')
    ax.set_title(
        f'{CONDITIONS[condition]} - {GO_ASPECTS[aspect]}',
        fontsize=14,
        fontweight='bold'
    )
    ax.axvline(x=1, color='gray', linestyle='--', alpha=0.5)
    
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('-log10(FDR)', fontsize=11, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()


# ============================================================================
# MAIN ANALYSIS PIPELINE
# ============================================================================

def main():
    """Main analysis pipeline."""
    
    # Parse arguments
    parser = argparse.ArgumentParser(
        description='Comprehensive GO enrichment analysis for P-body proteomics'
    )
    parser.add_argument(
        '--threshold',
        type=float,
        default=0.8,
        help='log2FC threshold (default: 0.8)'
    )
    parser.add_argument(
        '--no-gprofiler',
        action='store_true',
        help='Disable g:Profiler validation'
    )
    parser.add_argument(
        '--min-genes',
        type=int,
        default=3,
        help='Minimum genes in GO term (default: 3)'
    )
    parser.add_argument(
        '--top-n',
        type=int,
        default=20,
        help='Top N terms for plots (default: 20)'
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Verbose logging'
    )
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging(args.verbose)
    
    logger.info("="*80)
    logger.info("COMPREHENSIVE GO TERM ENRICHMENT ANALYSIS")
    logger.info("="*80)
    logger.info(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Threshold: log2FC > {args.threshold}")
    logger.info(f"g:Profiler: {'Disabled' if args.no_gprofiler else 'Enabled'}")
    logger.info(f"Min genes per GO term: {args.min_genes}")
    
    # ========================================================================
    # 1. Load data
    # ========================================================================
    logger.info("\n" + "="*80)
    logger.info("STEP 1: Loading data")
    logger.info("="*80)
    
    if not INPUT_FILE.exists():
        logger.error(f"Input file not found: {INPUT_FILE}")
        sys.exit(1)
    
    df = pd.read_csv(INPUT_FILE)
    logger.info(f"  ✓ Loaded {len(df)} proteins")
    
    # Load GO cache
    logger.info("  Loading GO annotation cache from Script 03...")
    go_annotations = load_go_cache()
    annotated_genes = set(go_annotations.keys())
    logger.info(f"  ✓ GO cache loaded: {len(annotated_genes)} genes")
    
    # ========================================================================
    # 2. Build GO term gene sets
    # ========================================================================
    logger.info("\n" + "="*80)
    logger.info("STEP 2: Building GO term gene sets")
    logger.info("="*80)
    
    go_term_genes = build_go_term_gene_sets(go_annotations, args.min_genes)
    
    aspect_counts = Counter([v['aspect'] for v in go_term_genes.values()])
    logger.info(f"  ✓ GO terms to test: {len(go_term_genes)}")
    logger.info(f"    BP: {aspect_counts.get('BP', 0)}")
    logger.info(f"    MF: {aspect_counts.get('MF', 0)}")
    logger.info(f"    CC: {aspect_counts.get('CC', 0)}")
    
    # ========================================================================
    # 3. Enrichment analysis per condition & aspect
    # ========================================================================
    logger.info("\n" + "="*80)
    logger.info("STEP 3: Performing enrichment analysis")
    logger.info("="*80)
    
    all_results = []
    
    for cond_col, cond_name in CONDITIONS.items():
        logger.info(f"\n{cond_name} ({cond_col}):")
        logger.info("-"*80)
        
        # Get enriched genes
        enriched_df = df[df[cond_col] > args.threshold]
        enriched_genes = set(enriched_df['protein'].unique()) & annotated_genes
        
        logger.info(f"  Enriched genes (>{args.threshold}): {len(enriched_genes)}")
        
        if len(enriched_genes) == 0:
            logger.warning(f"  No enriched genes for {cond_col}, skipping")
            continue
        
        # Analyze each GO aspect
        for aspect_code, aspect_name in GO_ASPECTS.items():
            logger.info(f"\n  {aspect_name} ({aspect_code}):")
            
            # Perform enrichment
            results = perform_go_enrichment(
                enriched_genes,
                go_term_genes,
                annotated_genes,
                aspect=aspect_code,
                min_genes_in_term=2
            )
            
            if results.empty:
                logger.info(f"    No enrichment")
                continue
            
            # Filter significant
            sig_results = results[results['FDR'] < 0.05]
            logger.info(f"    Total terms tested: {len(results)}")
            logger.info(f"    Significant (FDR<0.05): {len(sig_results)}")
            
            if not sig_results.empty:
                logger.info(f"    Top term: {sig_results.iloc[0]['GO_Term']}")
                logger.info(f"      Fold: {sig_results.iloc[0]['Fold_Enrichment']:.2f}x")
                logger.info(f"      FDR: {sig_results.iloc[0]['FDR']:.2e}")
            
            # Save results
            output_file = PROCESSED_DIR / f"go_enrichment_{cond_col}_{aspect_code}.csv"
            results.to_csv(output_file, index=False)
            logger.info(f"    ✓ Saved: {output_file.name}")
            
            # Create bubble plot
            plot_file = FIGURES_DIR / f"go_bubble_{cond_col}_{aspect_code}.png"
            create_bubble_plot(results, cond_col, aspect_code, plot_file, args.top_n)
            if plot_file.exists():
                logger.info(f"    ✓ Plot: {plot_file.name}")
            
            # Store for summary
            all_results.append({
                'Condition': cond_col,
                'Aspect': aspect_code,
                'N_Enriched': len(enriched_genes),
                'N_Terms_Tested': len(results),
                'N_Significant': len(sig_results)
            })
        
        # g:Profiler validation
        if not args.no_gprofiler:
            logger.info(f"\n  g:Profiler validation:")
            gene_list = enriched_df['gene_name'].dropna().unique().tolist()
            run_gprofiler(gene_list, cond_col, logger)
    
    # ========================================================================
    # 4. Summary
    # ========================================================================
    logger.info("\n" + "="*80)
    logger.info("STEP 4: Creating summary")
    logger.info("="*80)
    
    summary_df = pd.DataFrame(all_results)
    summary_file = PROCESSED_DIR / "go_enrichment_summary_all_conditions.csv"
    summary_df.to_csv(summary_file, index=False)
    logger.info(f"  ✓ Summary saved: {summary_file.name}")
    
    logger.info("\n" + summary_df.to_string(index=False))
    
    # ========================================================================
    # COMPLETION
    # ========================================================================
    logger.info("\n" + "="*80)
    logger.info("ANALYSIS COMPLETE")
    logger.info("="*80)
    logger.info(f"End time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"\nOutputs generated:")
    logger.info(f"  - go_enrichment_[condition]_[aspect].csv (up to 12 files)")
    logger.info(f"  - go_bubble_[condition]_[aspect].png (up to 12 plots)")
    logger.info(f"  - go_enrichment_summary_all_conditions.csv")
    if not args.no_gprofiler:
        logger.info(f"  - results/gprofiler/ (g:Profiler validation)")
    logger.info(f"  - {LOG_FILE.name}")


if __name__ == "__main__":
    main()
