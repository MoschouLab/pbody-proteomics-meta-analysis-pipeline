#!/usr/bin/env python3
"""
Venn/UpSet Analysis for P-body Proteomics
Comprehensive overlap analysis: AP vs PDL × No Stress vs Heat Stress

Includes:
- 2-way, 3-way, and custom 4-way Venn diagrams
- UpSet-style intersection plot
- Statistical overlap testing (hypergeometric + FDR correction)
- Bootstrap stability analysis
- Stress-specific directionality classification
- Network-aware coherence analysis (PPI-based)
- Core P-body gene GO enrichment analysis
- Comprehensive intersection gene lists
- Publication-ready composite figures

Author: Pepe (Lead Bioinformatics Postdoc)
PI: Prof. Panagiotis Moschou
Project: PLANTEX ERC Consolidator Grant
Institution: University of Crete, Greece
Date: 2026-01-23
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Set, Dict, List, Tuple, Optional
from itertools import combinations

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.patches as mpatches
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
from matplotlib_venn import venn2, venn2_circles, venn3, venn3_circles

# Optional: NetworkX for PPI analysis
try:
    import networkx as nx
    NETWORK_AVAILABLE = True
except ImportError:
    NETWORK_AVAILABLE = False
    print("Warning: networkx not installed. Network analysis will be skipped.")
    print("         Install with: conda install -c conda-forge networkx")


# =========================
# Configuration
# =========================
class Config:
    """Configuration parameters"""
    # Paths
    DATA_DIR = Path("data/processed")
    RESULTS_DIR = Path("results/figures")
    LOG_DIR = Path("logs")
    
    # Input files
    PROTEOMICS_FILE = DATA_DIR / "pbody_proteomics_annotated.csv"
    GO_ANNOTATIONS_FILE = DATA_DIR / "all_gene_GO_annotations.csv"
    PPI_FILE = DATA_DIR / "ppi_edges.csv"
    
    # Output subdirectories
    VENN_DIR = RESULTS_DIR / "venn_diagrams"
    INTERSECTIONS_DIR = DATA_DIR / "intersections"
    
    # Analysis parameters
    CONDITIONS = ['AP_NS', 'AP_HS', 'PDL_NS', 'PDL_HS']
    DEFAULT_THRESHOLD = 0.8
    DEFAULT_BOOTSTRAPS = 1000
    RANDOM_SEED = 42
    MIN_GENES_FOR_ANALYSIS = 3
    
    # Figure parameters
    FIGURE_DPI = 300
    FIGURE_FORMAT = 'svg'
    
    # Colors
    COLORS = {
        'AP_NS': '#FF6B6B',
        'AP_HS': '#4ECDC4',
        'PDL_NS': '#95E1D3',
        'PDL_HS': '#F38181'
    }
    
    GO_COLORS = {
        'Process': '#FF6B6B',
        'Function': '#4ECDC4',
        'Component': '#95E1D3'
    }


def setup_logging(log_file: Path) -> logging.Logger:
    """Configure logging"""
    log_file.parent.mkdir(parents=True, exist_ok=True)
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)


# =========================
# Core Data Functions
# =========================
def load_data(filepath: Path, logger: logging.Logger) -> pd.DataFrame:
    """Load proteomics data"""
    logger.info(f"Loading data from: {filepath}")
    df = pd.read_csv(filepath)
    
    required_cols = ['protein'] + Config.CONDITIONS
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    
    logger.info(f"Loaded {len(df)} proteins")
    return df


def get_enriched_genes(df: pd.DataFrame, condition: str, threshold: float) -> Set[str]:
    """Get set of genes enriched in a condition above threshold"""
    enriched = set(df.loc[df[condition] > threshold, 'protein'].dropna())
    return enriched


def compute_gene_sets(df: pd.DataFrame, threshold: float, logger: logging.Logger) -> Dict[str, Set[str]]:
    """Compute enriched gene sets for all conditions"""
    logger.info(f"Computing enriched gene sets (threshold={threshold})...")
    
    gene_sets = {}
    for condition in Config.CONDITIONS:
        genes = get_enriched_genes(df, condition, threshold)
        gene_sets[condition] = genes
        logger.info(f"  {condition}: {len(genes)} genes")
    
    return gene_sets


# =========================
# Venn Diagram Functions
# =========================
def create_venn2(set1: Set, set2: Set, labels: Tuple[str, str], 
                 title: str, filename: Path, logger: logging.Logger):
    """Create 2-way Venn diagram"""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Get colors for labels
    color1 = Config.COLORS.get(labels[0].split()[0], Config.COLORS['AP_NS'])
    color2 = Config.COLORS.get(labels[1].split()[0], Config.COLORS['PDL_NS'])
    
    venn = venn2([set1, set2], set_labels=labels, ax=ax,
                 set_colors=(color1, color2), alpha=0.7)
    
    venn2_circles([set1, set2], ax=ax, linewidth=2)
    
    # Enhance text
    for text in venn.set_labels:
        if text:
            text.set_fontsize(14)
            text.set_fontweight('bold')
    
    for text in venn.subset_labels:
        if text:
            text.set_fontsize(12)
    
    # Statistics box
    only1 = len(set1 - set2)
    only2 = len(set2 - set1)
    both = len(set1 & set2)
    total = len(set1 | set2)
    jaccard = both / total if total > 0 else 0
    
    stats_text = (f"Total: {total}\n"
                 f"{labels[0]}: {only1}\n"
                 f"{labels[1]}: {only2}\n"
                 f"Shared: {both}\n"
                 f"Jaccard: {jaccard:.3f}")
    
    ax.text(0.98, 0.02, stats_text, transform=ax.transAxes,
            fontsize=11, verticalalignment='bottom', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
    
    plt.tight_layout()
    plt.savefig(filename, format=Config.FIGURE_FORMAT, bbox_inches='tight', dpi=Config.FIGURE_DPI)
    logger.info(f"  Saved: {filename.name}")
    plt.close()


def create_venn3(set1: Set, set2: Set, set3: Set, labels: Tuple[str, str, str],
                 title: str, filename: Path, logger: logging.Logger):
    """Create 3-way Venn diagram"""
    fig, ax = plt.subplots(figsize=(12, 10))
    
    venn = venn3([set1, set2, set3], set_labels=labels, ax=ax,
                 set_colors=('#FF6B6B', '#4ECDC4', '#95E1D3'), alpha=0.7)
    
    venn3_circles([set1, set2, set3], ax=ax, linewidth=2)
    
    # Enhance text
    for text in venn.set_labels:
        if text:
            text.set_fontsize(14)
            text.set_fontweight('bold')
    
    for text in venn.subset_labels:
        if text:
            text.set_fontsize(11)
    
    # Statistics box
    only1 = len(set1 - set2 - set3)
    only2 = len(set2 - set1 - set3)
    only3 = len(set3 - set1 - set2)
    all_three = len(set1 & set2 & set3)
    total = len(set1 | set2 | set3)
    
    stats_text = (f"Total: {total}\n"
                 f"{labels[0]}: {only1}\n"
                 f"{labels[1]}: {only2}\n"
                 f"{labels[2]}: {only3}\n"
                 f"All three: {all_three}")
    
    ax.text(0.98, 0.02, stats_text, transform=ax.transAxes,
            fontsize=11, verticalalignment='bottom', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
    
    plt.tight_layout()
    plt.savefig(filename, format=Config.FIGURE_FORMAT, bbox_inches='tight', dpi=Config.FIGURE_DPI)
    logger.info(f"  Saved: {filename.name}")
    plt.close()


def create_venn4_custom(gene_sets: Dict[str, Set], title: str, 
                       filename: Path, logger: logging.Logger):
    """Create custom 4-way Venn diagram using ellipses"""
    fig, ax = plt.subplots(figsize=(16, 14))
    
    # Define ellipse positions
    ellipses = {
        'AP_NS':  {'center': (0.35, 0.6), 'width': 0.45, 'height': 0.55, 'angle': 45},
        'AP_HS':  {'center': (0.65, 0.6), 'width': 0.45, 'height': 0.55, 'angle': -45},
        'PDL_NS': {'center': (0.35, 0.4), 'width': 0.45, 'height': 0.55, 'angle': -45},
        'PDL_HS': {'center': (0.65, 0.4), 'width': 0.45, 'height': 0.55, 'angle': 45}
    }
    
    # Draw ellipses
    for name, params in ellipses.items():
        ellipse = Ellipse(params['center'], params['width'], params['height'],
                         angle=params['angle'], alpha=0.4,
                         facecolor=Config.COLORS[name], 
                         edgecolor='black', linewidth=2)
        ax.add_patch(ellipse)
    
    # Add labels
    label_positions = {
        'AP_NS': (0.15, 0.80),
        'AP_HS': (0.85, 0.80),
        'PDL_NS': (0.15, 0.20),
        'PDL_HS': (0.85, 0.20)
    }
    
    for name, pos in label_positions.items():
        ax.text(pos[0], pos[1], f"{name}\n({len(gene_sets[name])})",
               fontsize=12, fontweight='bold', ha='center',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Calculate core intersection
    all_four = set.intersection(*gene_sets.values())
    
    # Display core count prominently
    ax.text(0.5, 0.5, str(len(all_four)),
            ha='center', va='center', fontsize=22, fontweight='bold',
            bbox=dict(boxstyle='circle,pad=0.4', facecolor='gold', 
                     edgecolor='red', linewidth=3, alpha=0.95))
    
    # Summary statistics
    ap_only = (gene_sets['AP_NS'] | gene_sets['AP_HS']) - (gene_sets['PDL_NS'] | gene_sets['PDL_HS'])
    pdl_only = (gene_sets['PDL_NS'] | gene_sets['PDL_HS']) - (gene_sets['AP_NS'] | gene_sets['AP_HS'])
    ns_only = (gene_sets['AP_NS'] | gene_sets['PDL_NS']) - (gene_sets['AP_HS'] | gene_sets['PDL_HS'])
    hs_only = (gene_sets['AP_HS'] | gene_sets['PDL_HS']) - (gene_sets['AP_NS'] | gene_sets['PDL_NS'])
    
    stats_text = (
        f"CORE (all 4): {len(all_four)}\n"
        f"{'='*30}\n"
        f"AP only: {len(ap_only)}\n"
        f"PDL only: {len(pdl_only)}\n"
        f"NS only: {len(ns_only)}\n"
        f"HS only: {len(hs_only)}\n"
        f"\nTotal unique: {len(set.union(*gene_sets.values()))}"
    )
    
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
            fontsize=11, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))
    
    # Legend
    legend_text = (
        "Intersection Counts\n"
        "━━━━━━━━━━━━━━━━\n"
        "GOLD circle = Core\n"
        "(all 4 conditions)\n\n"
        "Note: Ellipse positions\n"
        "are illustrative;\n"
        "counts are exact."
    )
    ax.text(0.98, 0.98, legend_text, transform=ax.transAxes,
            fontsize=10, verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    
    ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect('equal')
    ax.axis('off')
    
    plt.tight_layout()
    plt.savefig(filename, format=Config.FIGURE_FORMAT, bbox_inches='tight', dpi=Config.FIGURE_DPI)
    logger.info(f"  Saved: {filename.name}")
    plt.close()


def generate_venn_diagrams(gene_sets: Dict[str, Set], output_dir: Path, 
                          logger: logging.Logger):
    """Generate all Venn diagrams"""
    logger.info("\nGenerating Venn diagrams...")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 2-way comparisons
    logger.info("  Creating 2-way comparisons...")
    
    # AP: NS vs HS
    create_venn2(gene_sets['AP_NS'], gene_sets['AP_HS'],
                ('AP_NS', 'AP_HS'),
                'Affinity Purification: No Stress vs Heat Stress',
                output_dir / 'venn_AP_NS_vs_HS.svg', logger)
    
    # PDL: NS vs HS
    create_venn2(gene_sets['PDL_NS'], gene_sets['PDL_HS'],
                ('PDL_NS', 'PDL_HS'),
                'Proximity-Dependent Labeling: No Stress vs Heat Stress',
                output_dir / 'venn_PDL_NS_vs_HS.svg', logger)
    
    # NS: AP vs PDL
    create_venn2(gene_sets['AP_NS'], gene_sets['PDL_NS'],
                ('AP_NS', 'PDL_NS'),
                'No Stress: AP vs PDL',
                output_dir / 'venn_NS_AP_vs_PDL.svg', logger)
    
    # HS: AP vs PDL
    create_venn2(gene_sets['AP_HS'], gene_sets['PDL_HS'],
                ('AP_HS', 'PDL_HS'),
                'Heat Stress: AP vs PDL',
                output_dir / 'venn_HS_AP_vs_PDL.svg', logger)
    
    # Combined: All AP vs All PDL
    all_ap = gene_sets['AP_NS'] | gene_sets['AP_HS']
    all_pdl = gene_sets['PDL_NS'] | gene_sets['PDL_HS']
    create_venn2(all_ap, all_pdl,
                ('AP (any)', 'PDL (any)'),
                'All Conditions: AP vs PDL',
                output_dir / 'venn_all_AP_vs_PDL.svg', logger)
    
    # Combined: All NS vs All HS
    all_ns = gene_sets['AP_NS'] | gene_sets['PDL_NS']
    all_hs = gene_sets['AP_HS'] | gene_sets['PDL_HS']
    create_venn2(all_ns, all_hs,
                ('NS (any)', 'HS (any)'),
                'All Conditions: No Stress vs Heat Stress',
                output_dir / 'venn_all_NS_vs_HS.svg', logger)
    
    # 3-way comparison
    logger.info("  Creating 3-way comparison...")
    create_venn3(gene_sets['AP_NS'], gene_sets['PDL_NS'], gene_sets['PDL_HS'],
                ('AP_NS', 'PDL_NS', 'PDL_HS'),
                'Three Major Conditions',
                output_dir / 'venn_3way_major.svg', logger)
    
    # 4-way custom
    logger.info("  Creating 4-way diagram...")
    create_venn4_custom(gene_sets,
                       'All Four Conditions: AP vs PDL × NS vs HS',
                       output_dir / 'venn_4way_all_conditions.svg', logger)


# =========================
# UpSet Plot
# =========================
def compute_intersections(gene_sets: Dict[str, Set], logger: logging.Logger) -> List[Dict]:
    """Compute all intersections for UpSet plot"""
    logger.info("\nComputing intersections...")
    
    conditions = list(gene_sets.keys())
    intersections = []
    
    for r in range(1, len(conditions) + 1):
        for combo in combinations(conditions, r):
            # Get intersection
            inter = set.intersection(*[gene_sets[c] for c in combo])
            
            # Subtract genes from other conditions (exclusive intersection)
            not_in = set()
            for c in conditions:
                if c not in combo:
                    not_in |= gene_sets[c]
            
            unique_to_combo = inter - not_in
            
            if len(unique_to_combo) >= Config.MIN_GENES_FOR_ANALYSIS:
                intersections.append({
                    'combination': combo,
                    'size': len(unique_to_combo),
                    'genes': sorted(unique_to_combo)
                })
    
    # Sort by size
    intersections.sort(key=lambda x: x['size'], reverse=True)
    logger.info(f"  Found {len(intersections)} non-empty intersections")
    
    return intersections


def create_upset_plot(intersections: List[Dict], gene_sets: Dict[str, Set],
                     filename: Path, logger: logging.Logger):
    """Create UpSet-style plot"""
    if len(intersections) == 0:
        logger.warning("No intersections to plot")
        return
    
    fig, (ax_set, ax_size) = plt.subplots(2, 1, figsize=(max(14, len(intersections)), 10),
                                          gridspec_kw={'height_ratios': [1, 2]})
    
    conditions = list(gene_sets.keys())
    n_intersections = len(intersections)
    
    # Top panel: Set membership matrix
    for idx, intersect in enumerate(intersections):
        for cond_idx, cond in enumerate(conditions):
            if cond in intersect['combination']:
                ax_set.scatter(idx, cond_idx, s=300, c=Config.COLORS[cond],
                             edgecolors='black', linewidth=2, zorder=3)
                # Draw connecting lines
                if len(intersect['combination']) > 1:
                    combo_list = list(intersect['combination'])
                    combo_indices = [conditions.index(c) for c in combo_list]
                    if cond in combo_list[:-1]:
                        next_idx = combo_indices[combo_list.index(cond) + 1]
                        ax_set.plot([idx, idx], [cond_idx, next_idx],
                                  'k-', linewidth=2, zorder=2)
    
    ax_set.set_yticks(range(len(conditions)))
    ax_set.set_yticklabels(conditions, fontsize=12, fontweight='bold')
    ax_set.set_xlim(-0.5, n_intersections - 0.5)
    ax_set.set_ylim(-0.5, len(conditions) - 0.5)
    ax_set.set_xticks([])
    ax_set.spines['top'].set_visible(False)
    ax_set.spines['right'].set_visible(False)
    ax_set.spines['bottom'].set_visible(False)
    ax_set.grid(axis='x', alpha=0.3, linestyle=':')
    ax_set.set_ylabel('Conditions', fontsize=12, fontweight='bold')
    
    # Bottom panel: Bar chart of intersection sizes
    x_pos = range(len(intersections))
    sizes = [i['size'] for i in intersections]
    bar_colors = []
    
    for intersect in intersections:
        if len(intersect['combination']) == 1:
            bar_colors.append(Config.COLORS[intersect['combination'][0]])
        elif len(intersect['combination']) == 4:
            bar_colors.append('#FFD700')  # Gold for all four
        else:
            bar_colors.append('#999999')  # Gray for others
    
    bars = ax_size.bar(x_pos, sizes, color=bar_colors, alpha=0.7,
                       edgecolor='black', linewidth=1.5)
    
    # Add value labels on bars
    for i, (bar, size) in enumerate(zip(bars, sizes)):
        height = bar.get_height()
        ax_size.text(bar.get_x() + bar.get_width()/2., height,
                    f'{int(size)}',
                    ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    ax_size.set_ylabel('Number of Genes', fontsize=12, fontweight='bold')
    ax_size.set_xlabel('Intersection', fontsize=12, fontweight='bold')
    ax_size.set_xlim(-0.5, n_intersections - 0.5)
    ax_size.spines['top'].set_visible(False)
    ax_size.spines['right'].set_visible(False)
    ax_size.grid(axis='y', alpha=0.3, linestyle=':')
    
    # Overall title
    fig.suptitle('UpSet Plot: Intersection Analysis', fontsize=16, fontweight='bold', y=0.98)
    
    plt.tight_layout()
    plt.savefig(filename, format=Config.FIGURE_FORMAT, bbox_inches='tight', dpi=Config.FIGURE_DPI)
    logger.info(f"  Saved: {filename.name}")
    plt.close()


# =========================
# Statistical Analysis
# =========================
def statistical_overlap_testing(gene_sets: Dict[str, Set], universe_size: int, 
                                logger: logging.Logger) -> pd.DataFrame:
    """Perform hypergeometric test for all pairwise overlaps with FDR correction"""
    logger.info("\nPerforming statistical overlap testing...")
    
    records = []
    for cond_a, cond_b in combinations(Config.CONDITIONS, 2):
        set_a = gene_sets[cond_a]
        set_b = gene_sets[cond_b]
        
        overlap = len(set_a & set_b)
        expected = len(set_a) * len(set_b) / universe_size
        
        # Hypergeometric test
        # P(X >= overlap) where X ~ Hypergeom(N=universe, K=|A|, n=|B|)
        p_value = hypergeom.sf(overlap - 1, universe_size, len(set_a), len(set_b))
        
        fold_enrichment = overlap / expected if expected > 0 else np.nan
        
        records.append({
            'condition_A': cond_a,
            'condition_B': cond_b,
            'set_A_size': len(set_a),
            'set_B_size': len(set_b),
            'overlap': overlap,
            'expected': expected,
            'fold_enrichment': fold_enrichment,
            'p_value': p_value
        })
    
    stats_df = pd.DataFrame(records)
    
    # FDR correction
    _, fdr_values, _, _ = multipletests(stats_df['p_value'], method='fdr_bh')
    stats_df['FDR'] = fdr_values
    
    # Sort by significance
    stats_df = stats_df.sort_values('fold_enrichment', ascending=False)
    
    logger.info(f"  Tested {len(stats_df)} pairwise comparisons")
    logger.info(f"  Significant overlaps (FDR < 0.05): {(stats_df['FDR'] < 0.05).sum()}")
    
    return stats_df


def bootstrap_overlap(df: pd.DataFrame, cond_a: str, cond_b: str, 
                     threshold: float, n_bootstraps: int, logger: logging.Logger) -> Tuple[float, Tuple[float, float]]:
    """Bootstrap confidence intervals for overlap"""
    overlaps = []
    
    for _ in range(n_bootstraps):
        # Resample with replacement
        boot_df = df.sample(frac=1, replace=True)
        set_a = get_enriched_genes(boot_df, cond_a, threshold)
        set_b = get_enriched_genes(boot_df, cond_b, threshold)
        overlaps.append(len(set_a & set_b))
    
    mean_overlap = np.mean(overlaps)
    ci_low, ci_high = np.percentile(overlaps, [2.5, 97.5])
    
    return mean_overlap, (ci_low, ci_high)


def bootstrap_stability_analysis(df: pd.DataFrame, threshold: float, 
                                 n_bootstraps: int, logger: logging.Logger) -> pd.DataFrame:
    """Perform bootstrap analysis for all pairwise overlaps"""
    logger.info(f"\nPerforming bootstrap stability analysis ({n_bootstraps} iterations)...")
    
    records = []
    for cond_a, cond_b in combinations(Config.CONDITIONS, 2):
        logger.info(f"  Bootstrapping {cond_a} ∩ {cond_b}...")
        mean_overlap, (ci_low, ci_high) = bootstrap_overlap(
            df, cond_a, cond_b, threshold, n_bootstraps, logger
        )
        
        records.append({
            'condition_A': cond_a,
            'condition_B': cond_b,
            'mean_overlap': mean_overlap,
            'CI_low': ci_low,
            'CI_high': ci_high,
            'CI_width': ci_high - ci_low
        })
    
    bootstrap_df = pd.DataFrame(records)
    logger.info(f"  Bootstrap analysis complete")
    
    return bootstrap_df


# =========================
# Directionality Analysis
# =========================
def classify_directionality(df: pd.DataFrame, ns_cond: str, hs_cond: str, 
                           threshold: float, logger: logging.Logger) -> Dict[str, str]:
    """Classify genes by their response to heat stress"""
    logger.info(f"\nClassifying directionality: {ns_cond} → {hs_cond}...")
    
    directionality = {}
    categories = {'HS_gain': 0, 'HS_loss': 0, 'Stable': 0, 'Absent': 0}
    
    for _, row in df.iterrows():
        gene = row['protein']
        enriched_ns = row[ns_cond] > threshold
        enriched_hs = row[hs_cond] > threshold
        
        if enriched_hs and not enriched_ns:
            category = 'HS_gain'
        elif enriched_ns and not enriched_hs:
            category = 'HS_loss'
        elif enriched_ns and enriched_hs:
            category = 'Stable'
        else:
            category = 'Absent'
        
        directionality[gene] = category
        categories[category] += 1
    
    logger.info(f"  HS_gain: {categories['HS_gain']}")
    logger.info(f"  HS_loss: {categories['HS_loss']}")
    logger.info(f"  Stable: {categories['Stable']}")
    logger.info(f"  Absent: {categories['Absent']}")
    
    return directionality


# =========================
# Network Analysis
# =========================
def load_ppi_network(ppi_file: Path, logger: logging.Logger) -> Optional[nx.Graph]:
    """Load PPI network from file"""
    if not ppi_file.exists():
        logger.warning(f"PPI file not found: {ppi_file}")
        logger.warning("Network analysis will be skipped")
        return None
    
    if not NETWORK_AVAILABLE:
        logger.warning("NetworkX not installed. Network analysis will be skipped")
        return None
    
    logger.info(f"\nLoading PPI network from: {ppi_file}")
    ppi_df = pd.read_csv(ppi_file)
    
    G = nx.from_pandas_edgelist(ppi_df, 'proteinA', 'proteinB')
    
    logger.info(f"  Nodes: {G.number_of_nodes():,}")
    logger.info(f"  Edges: {G.number_of_edges():,}")
    
    return G


def compute_network_coherence(gene_set: Set[str], network: Optional[nx.Graph], 
                              logger: logging.Logger) -> float:
    """Compute network coherence (average shortest path length)"""
    if network is None or len(gene_set) < 2:
        return np.nan
    
    # Filter genes present in network
    genes_in_network = [g for g in gene_set if g in network]
    
    if len(genes_in_network) < 2:
        return np.nan
    
    # Compute pairwise shortest paths
    path_lengths = []
    for i, gene_a in enumerate(genes_in_network):
        for gene_b in genes_in_network[i+1:]:
            try:
                length = nx.shortest_path_length(network, gene_a, gene_b)
                path_lengths.append(length)
            except nx.NetworkXNoPath:
                # Genes not connected
                pass
    
    if len(path_lengths) == 0:
        return np.nan
    
    return np.mean(path_lengths)


# =========================
# Core Gene GO Analysis
# =========================
def analyze_core_genes(core_genes: Set[str], df: pd.DataFrame, 
                      go_file: Path, output_dir: Path, logger: logging.Logger):
    """Analyze GO enrichment for core P-body genes"""
    logger.info(f"\nAnalyzing CORE P-body genes ({len(core_genes)} genes)...")
    
    if len(core_genes) == 0:
        logger.warning("No core genes found (present in all 4 conditions)")
        return
    
    # Get core gene details
    core_df = df[df['protein'].isin(core_genes)].copy()
    
    logger.info("\nCore P-body genes:")
    for _, row in core_df.iterrows():
        prot_name = str(row.get('protein_name', 'Unknown'))[:50]
        logger.info(f"  {row['protein']:<12} {prot_name:<50}")
        logger.info(f"    log2FC: AP_NS={row['AP_NS']:.2f}, AP_HS={row['AP_HS']:.2f}, "
                   f"PDL_NS={row['PDL_NS']:.2f}, PDL_HS={row['PDL_HS']:.2f}")
    
    # Load GO annotations
    if not go_file.exists():
        logger.warning(f"GO annotations not found: {go_file}")
        logger.warning("Run Script 03 first to generate GO annotations")
        return
    
    logger.info(f"\nLoading GO annotations from: {go_file}")
    go_df = pd.read_csv(go_file)
    
    # Filter for core genes
    core_go = go_df[go_df['gene_id'].isin(core_genes)]
    
    if len(core_go) == 0:
        logger.warning("No GO annotations found for core genes")
        return
    
    logger.info(f"  Found GO annotations for {len(core_go['gene_id'].unique())} core genes")
    
    # Group by GO term
    go_term_counts = core_go.groupby(['go_term', 'go_aspect', 'go_id']).size().reset_index(name='count')
    go_term_counts = go_term_counts.sort_values('count', ascending=False)
    
    logger.info(f"\nTop GO terms in Core P-body genes:")
    for _, row in go_term_counts.head(20).iterrows():
        aspect_abbr = row['go_aspect'][0] if pd.notna(row['go_aspect']) else '?'
        logger.info(f"  [{aspect_abbr}] {row['go_term'][:60]:<60} ({row['count']} genes)")
    
    # Save core GO terms
    core_go_file = output_dir / "core_pbody_GO_terms.csv"
    core_go.to_csv(core_go_file, index=False)
    logger.info(f"\nSaved: {core_go_file.name}")
    
    # Create visualization
    if len(go_term_counts) > 0:
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Top 15 GO terms
        top_terms = go_term_counts.head(15)
        
        # Color by aspect
        bar_colors = [Config.GO_COLORS.get(aspect, '#999999') 
                     for aspect in top_terms['go_aspect']]
        
        # Create horizontal bar chart
        y_pos = range(len(top_terms))
        bars = ax.barh(y_pos, top_terms['count'], color=bar_colors,
                       alpha=0.7, edgecolor='black', linewidth=1.5)
        
        # Labels
        labels = [term[:50] + '...' if len(term) > 50 else term
                 for term in top_terms['go_term']]
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontsize=10)
        ax.set_xlabel('Number of Core P-body Genes', fontsize=12, fontweight='bold')
        ax.set_title('GO Terms in CORE P-body Genes\n(Present in All 4 Conditions)',
                   fontsize=14, fontweight='bold')
        
        # Add legend
        legend_elements = [
            mpatches.Patch(facecolor='#FF6B6B', edgecolor='black', label='Process'),
            mpatches.Patch(facecolor='#4ECDC4', edgecolor='black', label='Function'),
            mpatches.Patch(facecolor='#95E1D3', edgecolor='black', label='Component')
        ]
        ax.legend(handles=legend_elements, loc='lower right', fontsize=10)
        
        ax.grid(axis='x', alpha=0.3, linestyle=':')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        plt.tight_layout()
        go_fig_file = output_dir / 'core_pbody_GO_terms.svg'
        plt.savefig(go_fig_file, format=Config.FIGURE_FORMAT, 
                   bbox_inches='tight', dpi=Config.FIGURE_DPI)
        logger.info(f"Saved: {go_fig_file.name}")
        plt.close()


# =========================
# Save Intersection Gene Lists
# =========================
def save_intersection_gene_lists(intersections: List[Dict], df: pd.DataFrame,
                                 output_dir: Path, logger: logging.Logger):
    """Save all intersection gene lists to CSV files"""
    logger.info(f"\nSaving intersection gene lists to: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    for idx, intersect in enumerate(intersections):
        genes = intersect['genes']
        combo_name = '_'.join(intersect['combination'])
        
        # Get full gene information
        gene_df = df[df['protein'].isin(genes)].copy()
        
        # Sort by average log2FC
        gene_df['avg_log2FC'] = gene_df[Config.CONDITIONS].mean(axis=1)
        gene_df = gene_df.sort_values('avg_log2FC', ascending=False)
        
        # Save
        filename = output_dir / f"intersection_{idx+1:02d}_{combo_name}.csv"
        gene_df.to_csv(filename, index=False)
    
    logger.info(f"  Saved {len(intersections)} intersection files")


# =========================
# Composite Figure
# =========================
def create_composite_figure(stats_df: pd.DataFrame, bootstrap_df: pd.DataFrame,
                           dir_ap: Dict[str, str], dir_pdl: Dict[str, str],
                           intersections: List[Dict], network: Optional[nx.Graph],
                           filename: Path, logger: logging.Logger):
    """Create publication-ready composite figure"""
    logger.info("\nCreating composite publication figure...")
    
    fig = plt.figure(figsize=(18, 12))
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)
    
    # Panel A: Statistical overlap (FDR-controlled)
    ax1 = fig.add_subplot(gs[0, 0])
    stats_plot = stats_df.sort_values('fold_enrichment', ascending=True)
    labels = [f"{row['condition_A']} ∩ {row['condition_B']}" 
             for _, row in stats_plot.iterrows()]
    colors = ['crimson' if fdr < 0.05 else 'gray' for fdr in stats_plot['FDR']]
    
    ax1.barh(range(len(stats_plot)), stats_plot['fold_enrichment'], color=colors)
    ax1.set_yticks(range(len(stats_plot)))
    ax1.set_yticklabels(labels, fontsize=9)
    ax1.set_xlabel('Fold Enrichment', fontweight='bold')
    ax1.set_title('A) Statistical Overlap (FDR-controlled)', fontweight='bold')
    ax1.axvline(1, color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax1.grid(axis='x', alpha=0.3)
    
    # Panel B: Bootstrap stability
    ax2 = fig.add_subplot(gs[0, 1])
    boot_labels = [f"{row['condition_A']} ∩ {row['condition_B']}" 
                  for _, row in bootstrap_df.iterrows()]
    
    for idx, row in bootstrap_df.iterrows():
        ax2.errorbar(row['mean_overlap'], idx,
                    xerr=[[row['mean_overlap'] - row['CI_low']],
                          [row['CI_high'] - row['mean_overlap']]],
                    fmt='o', color='black', capsize=5)
    
    ax2.set_yticks(range(len(bootstrap_df)))
    ax2.set_yticklabels(boot_labels, fontsize=9)
    ax2.set_xlabel('Overlap Size', fontweight='bold')
    ax2.set_title('B) Bootstrap Stability (95% CI)', fontweight='bold')
    ax2.grid(axis='x', alpha=0.3)
    
    # Panel C: Directionality (AP)
    ax3 = fig.add_subplot(gs[1, 0])
    dir_ap_counts = pd.Series(dir_ap).value_counts()
    dir_ap_counts.plot(kind='bar', ax=ax3, color='steelblue', edgecolor='black')
    ax3.set_title('C) AP Stress Directionality', fontweight='bold')
    ax3.set_xlabel('Category', fontweight='bold')
    ax3.set_ylabel('Count', fontweight='bold')
    ax3.set_xticklabels(ax3.get_xticklabels(), rotation=45, ha='right')
    ax3.grid(axis='y', alpha=0.3)
    
    # Panel D: Intersection sizes (UpSet-style bars)
    ax4 = fig.add_subplot(gs[1, 1])
    if len(intersections) > 0:
        sizes = [i['size'] for i in intersections]
        x_pos = range(min(len(sizes), 10))  # Top 10
        
        # Color by combination size
        colors_d = []
        for i in intersections[:10]:
            if len(i['combination']) == 4:
                colors_d.append('#FFD700')  # Gold for all four
            elif len(i['combination']) == 1:
                colors_d.append(Config.COLORS[i['combination'][0]])
            else:
                colors_d.append('#999999')
        
        ax4.bar(x_pos, sizes[:10], color=colors_d, edgecolor='black')
        ax4.set_xlabel('Intersection (sorted by size)', fontweight='bold')
        ax4.set_ylabel('Gene Count', fontweight='bold')
        ax4.set_title('D) Top 10 Intersections', fontweight='bold')
        ax4.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(filename, format='png', bbox_inches='tight', dpi=Config.FIGURE_DPI)
    logger.info(f"  Saved: {filename.name}")
    plt.close()


# =========================
# Main Function
# =========================
def main():
    """Main execution function"""
    # Parse arguments
    parser = argparse.ArgumentParser(
        description='Venn/UpSet Analysis for P-body Proteomics',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--threshold', type=float, default=Config.DEFAULT_THRESHOLD,
                       help='log2FC threshold for enrichment')
    parser.add_argument('--bootstraps', type=int, default=Config.DEFAULT_BOOTSTRAPS,
                       help='Number of bootstrap iterations')
    parser.add_argument('--no-network', action='store_true',
                       help='Skip network analysis even if PPI data available')
    parser.add_argument('--no-bootstrap', action='store_true',
                       help='Skip bootstrap analysis (saves time)')
    
    args = parser.parse_args()
    
    # Set random seed
    np.random.seed(Config.RANDOM_SEED)
    
    # Setup logging
    log_file = Config.LOG_DIR / '05_venn_upset_analysis.log'
    logger = setup_logging(log_file)
    
    logger.info("="*90)
    logger.info("VENN/UPSET ANALYSIS - P-BODY PROTEOMICS")
    logger.info("="*90)
    logger.info(f"Threshold: {args.threshold}")
    logger.info(f"Bootstrap iterations: {args.bootstraps}")
    logger.info(f"Network analysis: {'disabled' if args.no_network else 'enabled'}")
    
    try:
        # Load data
        df = load_data(Config.PROTEOMICS_FILE, logger)
        
        # Compute gene sets
        gene_sets = compute_gene_sets(df, args.threshold, logger)
        
        # Universe size (all proteins)
        universe_size = len(df)
        
        # Generate Venn diagrams
        generate_venn_diagrams(gene_sets, Config.VENN_DIR, logger)
        
        # Compute intersections
        intersections = compute_intersections(gene_sets, logger)
        
        # Save intersection gene lists
        save_intersection_gene_lists(intersections, df, Config.INTERSECTIONS_DIR, logger)
        
        # Create UpSet plot
        create_upset_plot(intersections, gene_sets, 
                         Config.VENN_DIR / 'upset_all_conditions.svg', logger)
        
        # Statistical overlap testing
        stats_df = statistical_overlap_testing(gene_sets, universe_size, logger)
        stats_df.to_csv(Config.DATA_DIR / 'overlap_statistics_FDR.csv', index=False)
        logger.info(f"Saved: overlap_statistics_FDR.csv")
        
        # Bootstrap stability (optional)
        if not args.no_bootstrap:
            bootstrap_df = bootstrap_stability_analysis(df, args.threshold, 
                                                       args.bootstraps, logger)
            bootstrap_df.to_csv(Config.DATA_DIR / 'bootstrap_overlap_CI.csv', index=False)
            logger.info(f"Saved: bootstrap_overlap_CI.csv")
        else:
            bootstrap_df = pd.DataFrame()  # Empty for composite figure
        
        # Directionality analysis
        dir_ap = classify_directionality(df, 'AP_NS', 'AP_HS', args.threshold, logger)
        dir_pdl = classify_directionality(df, 'PDL_NS', 'PDL_HS', args.threshold, logger)
        
        # Save directionality
        pd.DataFrame.from_dict(dir_ap, orient='index', columns=['category']).to_csv(
            Config.DATA_DIR / 'directionality_AP.csv')
        pd.DataFrame.from_dict(dir_pdl, orient='index', columns=['category']).to_csv(
            Config.DATA_DIR / 'directionality_PDL.csv')
        logger.info(f"Saved: directionality_AP.csv, directionality_PDL.csv")
        
        # Network analysis (optional)
        network = None
        if not args.no_network:
            network = load_ppi_network(Config.PPI_FILE, logger)
            
            if network is not None:
                logger.info("\nComputing network coherence for intersections...")
                for intersect in intersections:
                    genes = set(intersect['genes'])
                    coherence = compute_network_coherence(genes, network, logger)
                    intersect['network_coherence'] = coherence
                    if not np.isnan(coherence):
                        combo_name = '_'.join(intersect['combination'])
                        logger.info(f"  {combo_name}: coherence = {coherence:.2f}")
        
        # Core gene analysis
        core_genes = set.intersection(*gene_sets.values())
        analyze_core_genes(core_genes, df, Config.GO_ANNOTATIONS_FILE, 
                          Config.DATA_DIR, logger)
        
        # Create composite figure
        if not args.no_bootstrap:
            create_composite_figure(stats_df, bootstrap_df, dir_ap, dir_pdl,
                                   intersections, network,
                                   Config.RESULTS_DIR / 'composite_overlap_analysis.png',
                                   logger)
        
        # Summary statistics
        logger.info("\n" + "="*90)
        logger.info("SUMMARY STATISTICS")
        logger.info("="*90)
        all_genes = set.union(*gene_sets.values())
        logger.info(f"Total unique enriched genes: {len(all_genes)}")
        logger.info(f"Genes in all 4 conditions (CORE): {len(core_genes)}")
        logger.info(f"Total intersections found: {len(intersections)}")
        logger.info(f"Significant overlaps (FDR < 0.05): {(stats_df['FDR'] < 0.05).sum()}")
        
        logger.info("\n" + "="*90)
        logger.info("VENN/UPSET ANALYSIS COMPLETE")
        logger.info("="*90)
        logger.info("\nOutput files:")
        logger.info(f"  Venn diagrams: {Config.VENN_DIR}/ (8 SVG files)")
        logger.info(f"  Intersection gene lists: {Config.INTERSECTIONS_DIR}/ ({len(intersections)} CSV files)")
        logger.info(f"  Statistical results: {Config.DATA_DIR}/overlap_statistics_FDR.csv")
        if not args.no_bootstrap:
            logger.info(f"  Bootstrap results: {Config.DATA_DIR}/bootstrap_overlap_CI.csv")
        logger.info(f"  Directionality: {Config.DATA_DIR}/directionality_*.csv")
        logger.info(f"  Core gene GO: {Config.DATA_DIR}/core_pbody_GO_terms.csv")
        logger.info(f"  Composite figure: {Config.RESULTS_DIR}/composite_overlap_analysis.png")
        
    except Exception as e:
        logger.error(f"Error during analysis: {e}", exc_info=True)
        raise


if __name__ == "__main__":
    main()
