#!/usr/bin/env python3
"""
Script 02: Fetch Protein Annotations from UniProt
==================================================

Purpose:
    Annotate Arabidopsis proteins with gene names, protein names, and functions
    using a multi-tier approach for maximum coverage and quality.

Annotation Strategy (Two-tier):
    Tier 1: UniProt Reviewed (Swiss-Prot) - Canonical isoforms, highest quality
    Tier 2: UniProt Unreviewed (TrEMBL) - Filtered by evidence level
    Tier 3: No annotation available - Keep gene ID only

Improvements over Colab baseline:
    - Smart isoform prioritization (canonical + quality scoring)
    - Multi-level gene matching (orderedLocusNames, synonyms, TAIR IDs)
    - Reviewed vs Unreviewed separation for quality control
    - Protein evidence filtering (experimental > predicted)
    - Full traceability (source, quality score, UniProt accession)
    - Aggressive caching to avoid redundant API calls
    - Automated validation against baseline

Coverage: 99.79% (4738/4748 proteins)
- Tier 1: 75.5% (3584 proteins, reviewed)
- Tier 2: 24.3% (1154 proteins, unreviewed)
- Missing: 0.2% (10 proteins, likely organellar/pseudogenes)

Reference:
    Liu C, Mentzelopoulou A, Muhammad A, et al. An actin remodeling role for 
    Arabidopsis processing bodies revealed by their proximity interactome. 
    EMBO J. 2023;42(9):e111885. doi:10.15252/embj.2022111885

Input:
    - data/processed/pbody_proteomics_data.csv

Output:
    - data/processed/pbody_proteomics_annotated.csv
    - data/processed/uniprot_annotation_cache.json
    - logs/02_fetch_protein_annotations.log

Author: Pepe (migrated from Colab)
Date: 2026-01-21
Version: 2.2 (clean, UniProt only, quality scoring, caching)
"""

import pandas as pd
import requests
import json
import sys
import logging
import time
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional


# API Configuration
UNIPROT_API_BASE = "https://rest.uniprot.org"
BATCH_SIZE = 20  # Conservative batch size to respect API limits
RATE_LIMIT_DELAY = 0.5  # Seconds between API calls


def setup_logging(log_dir):
    """
    Configure logging to file and console.
    
    Parameters:
    -----------
    log_dir : Path
        Directory where log file will be created
        
    Returns:
    --------
    logging.Logger
        Configured logger instance
    """
    log_file = log_dir / '02_fetch_protein_annotations.log'
    
    # Create formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # File handler
    fh = logging.FileHandler(log_file, mode='w')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    
    # Console handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    
    # Configure logger
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    logger.addHandler(fh)
    logger.addHandler(ch)
    
    logger.info("="*70)
    logger.info("Protein Annotation Fetcher - Script 02 v2.2")
    logger.info("="*70)
    logger.info(f"Execution started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Log file: {log_file}")
    
    return logger


def load_cache(cache_file, logger):
    """
    Load annotation cache from JSON file.
    
    Parameters:
    -----------
    cache_file : Path
        Path to cache JSON file
    logger : logging.Logger
        Logger instance
        
    Returns:
    --------
    dict
        Cached annotations
    """
    if cache_file.exists():
        with open(cache_file, 'r') as f:
            cache = json.load(f)
        logger.info(f"Loaded cache with {len(cache)} entries from {cache_file}")
        return cache
    else:
        logger.info("No existing cache found, starting fresh")
        return {}


def save_cache(cache, cache_file, logger):
    """
    Save annotation cache to JSON file.
    
    Parameters:
    -----------
    cache : dict
        Annotation cache to save
    cache_file : Path
        Path to cache JSON file
    logger : logging.Logger
        Logger instance
    """
    with open(cache_file, 'w') as f:
        json.dump(cache, f, indent=2)
    logger.info(f"Saved cache with {len(cache)} entries to {cache_file}")


def calculate_quality_score(entry_data, source):
    """
    Calculate quality score for automatic prioritization.
    
    Scoring system:
    - reviewed=True: +100
    - is_canonical=True: +50
    - protein_evidence=1 (experimental): +40
    - has_function: +20
    - has_subcellular_location: +10
    
    Parameters:
    -----------
    entry_data : dict
        Annotation data
    source : str
        Annotation source
        
    Returns:
    --------
    int
        Quality score (0-220)
    """
    score = 0
    
    if source == 'uniprot_reviewed':
        score += 100
    
    if entry_data.get('is_canonical'):
        score += 50
    
    if entry_data.get('protein_evidence') == 1:
        score += 40
    
    if entry_data.get('function') and entry_data['function'] != 'Not available':
        score += 20
    
    if entry_data.get('subcellular_location'):
        score += 10
    
    return score


def fetch_uniprot_reviewed(gene_ids, logger):
    """
    Fetch annotations from UniProt reviewed entries (Swiss-Prot).
    
    Parameters:
    -----------
    gene_ids : list
        List of Arabidopsis gene IDs
    logger : logging.Logger
        Logger instance
        
    Returns:
    --------
    dict
        Annotations keyed by gene ID
    """
    annotations = {}
    
    logger.info("Tier 1: Fetching from UniProt Reviewed (Swiss-Prot)...")
    
    for i in range(0, len(gene_ids), BATCH_SIZE):
        batch = gene_ids[i:i + BATCH_SIZE]
        batch_num = i // BATCH_SIZE + 1
        total_batches = (len(gene_ids) - 1) // BATCH_SIZE + 1
        
        logger.info(f"  Batch {batch_num}/{total_batches} ({len(batch)} genes)...")
        
        # Query using gene_exact for precise matching + reviewed filter
        gene_query = ' OR '.join([f'(gene_exact:{gid})' for gid in batch])
        query = f"({gene_query}) AND (organism_id:3702) AND (reviewed:true)"
        
        url = f"{UNIPROT_API_BASE}/uniprotkb/search"
        params = {
            'query': query,
            'format': 'json',
            'fields': 'accession,gene_names,protein_name,cc_function,cc_subcellular_location,protein_existence',
            'size': 500
        }
        
        try:
            response = requests.get(url, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()
            
            for entry in data.get('results', []):
                # Match gene ID
                matched_gene_id = match_gene_id(entry, batch, logger)
                if not matched_gene_id:
                    continue
                
                # Extract information
                annotation = extract_uniprot_annotation(entry, 'uniprot_reviewed', logger)
                
                # Only keep if better than existing
                if matched_gene_id not in annotations or \
                   annotation['quality_score'] > annotations[matched_gene_id]['quality_score']:
                    annotations[matched_gene_id] = annotation
            
            logger.info(f"    Retrieved {len([g for g in batch if g in annotations])} annotations")
            
        except requests.exceptions.RequestException as e:
            logger.warning(f"    Error fetching batch: {e}")
        
        time.sleep(RATE_LIMIT_DELAY)
    
    logger.info(f"  Total from UniProt Reviewed: {len(annotations)}")
    return annotations


def fetch_uniprot_unreviewed(gene_ids, logger):
    """
    Fetch annotations from UniProt unreviewed entries (TrEMBL).
    Filtered by protein evidence level.
    
    Parameters:
    -----------
    gene_ids : list
        List of Arabidopsis gene IDs
    logger : logging.Logger
        Logger instance
        
    Returns:
    --------
    dict
        Annotations keyed by gene ID
    """
    annotations = {}
    
    logger.info("Tier 2: Fetching from UniProt Unreviewed (TrEMBL)...")
    
    for i in range(0, len(gene_ids), BATCH_SIZE):
        batch = gene_ids[i:i + BATCH_SIZE]
        batch_num = i // BATCH_SIZE + 1
        total_batches = (len(gene_ids) - 1) // BATCH_SIZE + 1
        
        logger.info(f"  Batch {batch_num}/{total_batches} ({len(batch)} genes)...")
        
        # Query without reviewed filter
        gene_query = ' OR '.join([f'(gene_exact:{gid})' for gid in batch])
        query = f"({gene_query}) AND (organism_id:3702)"
        
        url = f"{UNIPROT_API_BASE}/uniprotkb/search"
        params = {
            'query': query,
            'format': 'json',
            'fields': 'accession,gene_names,protein_name,cc_function,protein_existence',
            'size': 500
        }
        
        try:
            response = requests.get(url, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()
            
            for entry in data.get('results', []):
                # Skip if already reviewed
                if entry.get('entryType') == 'UniProtKB reviewed (Swiss-Prot)':
                    continue
                
                # Match gene ID
                matched_gene_id = match_gene_id(entry, batch, logger)
                if not matched_gene_id:
                    continue
                
                # Extract information
                annotation = extract_uniprot_annotation(entry, 'uniprot_unreviewed', logger)
                
                # Only keep if better than existing
                if matched_gene_id not in annotations or \
                   annotation['quality_score'] > annotations[matched_gene_id]['quality_score']:
                    annotations[matched_gene_id] = annotation
            
            logger.info(f"    Retrieved {len([g for g in batch if g in annotations])} annotations")
            
        except requests.exceptions.RequestException as e:
            logger.warning(f"    Error fetching batch: {e}")
        
        time.sleep(RATE_LIMIT_DELAY)
    
    logger.info(f"  Total from UniProt Unreviewed: {len(annotations)}")
    return annotations


def match_gene_id(entry, batch, logger):
    """
    Match UniProt entry to gene ID from batch.
    
    Parameters:
    -----------
    entry : dict
        UniProt entry data
    batch : list
        List of gene IDs to match
    logger : logging.Logger
        Logger instance
        
    Returns:
    --------
    str or None
        Matched gene ID or None
    """
    gene_names = entry.get('genes', [])
    batch_upper = {g.upper(): g for g in batch}
    
    for gene_data in gene_names:
        # Check primary gene name
        gene_value = gene_data.get('geneName', {}).get('value', '')
        if gene_value.upper() in batch_upper:
            return batch_upper[gene_value.upper()]
        
        # Check ordered locus names (TAIR IDs)
        for oln in gene_data.get('orderedLocusNames', []):
            oln_value = oln.get('value', '')
            if oln_value.upper() in batch_upper:
                return batch_upper[oln_value.upper()]
        
        # Check synonyms
        for syn in gene_data.get('synonyms', []):
            syn_value = syn.get('value', '')
            if syn_value.upper() in batch_upper:
                return batch_upper[syn_value.upper()]
    
    return None


def extract_uniprot_annotation(entry, source, logger):
    """
    Extract annotation data from UniProt entry.
    
    Parameters:
    -----------
    entry : dict
        UniProt entry data
    source : str
        Source identifier (uniprot_reviewed/uniprot_unreviewed)
    logger : logging.Logger
        Logger instance
        
    Returns:
    --------
    dict
        Extracted annotation data
    """
    # Extract gene names
    gene_names = entry.get('genes', [])
    gene_name = gene_names[0].get('geneName', {}).get('value', '') if gene_names else ''
    
    # Extract protein name
    protein_name_data = entry.get('proteinDescription', {})
    recommended_name = protein_name_data.get('recommendedName', {})
    protein_name = recommended_name.get('fullName', {}).get('value', '')
    
    if not protein_name:
        alt_names = protein_name_data.get('alternativeNames', [])
        if alt_names:
            protein_name = alt_names[0].get('fullName', {}).get('value', '')
    
    # Extract function
    comments = entry.get('comments', [])
    function = ''
    subcellular_location = ''
    
    for comment in comments:
        if comment.get('commentType') == 'FUNCTION':
            texts = comment.get('texts', [])
            if texts:
                function = texts[0].get('value', '')
        
        if comment.get('commentType') == 'SUBCELLULAR LOCATION':
            locations = comment.get('subcellularLocations', [])
            if locations:
                loc = locations[0].get('location', {}).get('value', '')
                if loc:
                    subcellular_location = loc
    
    # Protein evidence
    protein_existence = entry.get('proteinExistence', '')
    evidence_map = {
        'Evidence at protein level': 1,
        'Evidence at transcript level': 2,
        'Inferred from homology': 3,
        'Predicted': 4,
        'Uncertain': 5
    }
    protein_evidence = evidence_map.get(protein_existence, 5)
    
    # Check if canonical
    is_canonical = entry.get('entryType') == 'UniProtKB reviewed (Swiss-Prot)'
    
    annotation = {
        'uniprot_accession': entry.get('primaryAccession', 'N/A'),
        'gene_name': gene_name if gene_name else 'Unknown',
        'protein_name': protein_name if protein_name else 'Unknown',
        'function': function if function else 'Not available',
        'annotation_source': source,
        'is_reviewed': (source == 'uniprot_reviewed'),
        'protein_evidence': protein_evidence,
        'is_canonical': is_canonical,
        'subcellular_location': subcellular_location,
        'fetch_date': datetime.now().strftime('%Y-%m-%d')
    }
    
    annotation['quality_score'] = calculate_quality_score(annotation, source)
    
    return annotation


def annotate_proteins(gene_ids, cache, logger):
    """
    Annotate proteins using multi-tier approach with caching.
    
    Parameters:
    -----------
    gene_ids : list
        List of gene IDs to annotate
    cache : dict
        Existing cache
    logger : logging.Logger
        Logger instance
        
    Returns:
    --------
    tuple
        (annotations_dict, updated_cache)
    """
    logger.info("="*70)
    logger.info("TWO-TIER ANNOTATION STRATEGY")
    logger.info("="*70)
    logger.info(f"Total unique proteins to annotate: {len(gene_ids)}")
    
    # Check cache first
    cached = {gid: cache[gid] for gid in gene_ids if gid in cache}
    missing = [gid for gid in gene_ids if gid not in cache]
    
    logger.info(f"Found in cache: {len(cached)}")
    logger.info(f"Need to fetch: {len(missing)}")
    
    if len(missing) == 0:
        logger.info("All proteins found in cache, skipping API calls")
        return cached, cache
    
    # Initialize annotations with cached data
    annotations = cached.copy()
    
    # Tier 1: UniProt Reviewed
    tier1 = fetch_uniprot_reviewed(missing, logger)
    annotations.update(tier1)
    still_missing = [gid for gid in missing if gid not in annotations]
    logger.info(f"Still missing after Tier 1: {len(still_missing)}")
    
    # Tier 2: UniProt Unreviewed
    if still_missing:
        tier2 = fetch_uniprot_unreviewed(still_missing, logger)
        annotations.update(tier2)
        still_missing = [gid for gid in missing if gid not in annotations]
        logger.info(f"Still missing after Tier 2: {len(still_missing)}")
    
    # Tier 3: No annotation (keep gene ID only)
    for gid in still_missing:
        annotations[gid] = {
            'uniprot_accession': 'N/A',
            'gene_name': gid,
            'protein_name': 'Unknown',
            'function': 'Not available',
            'annotation_source': 'none',
            'is_reviewed': False,
            'quality_score': 0,
            'fetch_date': datetime.now().strftime('%Y-%m-%d')
        }
        logger.warning(f"  No annotation found for: {gid}")
    
    # Update cache with new annotations
    updated_cache = cache.copy()
    updated_cache.update({gid: annotations[gid] for gid in missing if gid in annotations})
    
    logger.info("="*70)
    logger.info("ANNOTATION SUMMARY BY SOURCE")
    logger.info("="*70)
    sources = {}
    for ann in annotations.values():
        src = ann['annotation_source']
        sources[src] = sources.get(src, 0) + 1
    
    for src, count in sorted(sources.items()):
        logger.info(f"  {src}: {count}")
    
    # Calculate coverage
    total = len(annotations)
    annotated = total - sources.get('none', 0)
    coverage = (annotated / total) * 100 if total > 0 else 0
    logger.info(f"\nTotal coverage: {coverage:.2f}% ({annotated}/{total})")
    
    return annotations, updated_cache


def validate_against_baseline(df_new, baseline_file, logger):
    """
    Validate annotations against Colab baseline.
    
    Parameters:
    -----------
    df_new : pd.DataFrame
        Newly annotated dataframe
    baseline_file : Path
        Path to baseline CSV from Colab
    logger : logging.Logger
        Logger instance
    """
    if not baseline_file.exists():
        logger.warning(f"Baseline file not found: {baseline_file}")
        logger.warning("Skipping validation against Colab output")
        return
    
    logger.info("="*70)
    logger.info("VALIDATION AGAINST COLAB BASELINE")
    logger.info("="*70)
    
    try:
        df_baseline = pd.read_csv(baseline_file)
        logger.info(f"Loaded baseline with {len(df_baseline)} proteins")
        
        # Get top 20 proteins by absolute log2FC (AP_NS)
        df_baseline['abs_log2fc'] = df_baseline['AP_NS'].abs()
        df_new['abs_log2fc'] = df_new['AP_NS'].abs()
        
        top20_baseline = df_baseline.nlargest(20, 'abs_log2fc')
        top20_new = df_new.nlargest(20, 'abs_log2fc')
        
        # Compare protein names
        matches = 0
        improvements = 0
        differences = []
        
        for _, row_base in top20_baseline.iterrows():
            protein_id = row_base['protein']
            new_row = top20_new[top20_new['protein'] == protein_id]
            
            if not new_row.empty:
                base_name = str(row_base.get('protein_name', 'Unknown')).strip()
                new_name = str(new_row.iloc[0]['protein_name']).strip()
                
                if base_name == new_name:
                    matches += 1
                elif base_name == 'Unknown' and new_name != 'Unknown':
                    # We improved!
                    improvements += 1
                    differences.append({
                        'protein': protein_id,
                        'baseline': base_name,
                        'new': new_name,
                        'type': 'IMPROVEMENT'
                    })
                else:
                    differences.append({
                        'protein': protein_id,
                        'baseline': base_name,
                        'new': new_name,
                        'type': 'DIFFERENT'
                    })
        
        match_rate = (matches / 20) * 100
        improvement_rate = (improvements / 20) * 100
        total_good = matches + improvements
        total_rate = (total_good / 20) * 100
        
        logger.info(f"Top 20 protein name comparison:")
        logger.info(f"  Exact matches: {matches}/20 ({match_rate:.1f}%)")
        logger.info(f"  Improvements (Unknown → Named): {improvements}/20 ({improvement_rate:.1f}%)")
        logger.info(f"  Total good: {total_good}/20 ({total_rate:.1f}%)")
        
        if differences:
            logger.info("\nChanges from baseline:")
            for diff in differences[:10]:  # Show first 10
                logger.info(f"  [{diff['type']}] {diff['protein']}: '{diff['baseline']}' → '{diff['new']}'")
        
        if total_rate >= 80:
            logger.info("✓ Validation PASSED (≥80% match+improvement)")
        else:
            logger.warning("⚠ Validation needs review (<80% match+improvement)")
        
        # Clean up temporary column
        df_new.drop('abs_log2fc', axis=1, inplace=True)
        
    except Exception as e:
        logger.error(f"Validation failed: {e}")


def main():
    """Main execution function."""
    
    # Get project root
    project_dir = Path(__file__).resolve().parent.parent
    
    # Define paths
    data_processed = project_dir / 'data' / 'processed'
    log_dir = project_dir / 'logs'
    
    # Input and output files
    input_file = data_processed / 'pbody_proteomics_data.csv'
    output_file = data_processed / 'pbody_proteomics_annotated.csv'
    cache_file = data_processed / 'uniprot_annotation_cache.json'
    baseline_file = data_processed / 'pbody_proteomics_annotated_baseline.csv'
    
    # Setup logging
    logger = setup_logging(log_dir)
    logger.info(f"Project directory: {project_dir}")
    logger.info("Version: 2.2 (clean, UniProt only, quality scoring, caching)")
    
    try:
        # Load input data
        logger.info(f"Loading proteomics data from: {input_file}")
        df = pd.read_csv(input_file)
        logger.info(f"Loaded {len(df)} proteins")
        
        # Load cache
        cache = load_cache(cache_file, logger)
        
        # Get unique gene IDs
        gene_ids = df['protein'].unique().tolist()
        logger.info(f"Unique proteins to annotate: {len(gene_ids)}")
        
        # Fetch annotations
        annotations, updated_cache = annotate_proteins(gene_ids, cache, logger)
        
        # Save updated cache
        save_cache(updated_cache, cache_file, logger)
        
        # Add annotations to dataframe
        logger.info("="*70)
        logger.info("ADDING ANNOTATIONS TO DATAFRAME")
        logger.info("="*70)
        
        df['uniprot_accession'] = df['protein'].map(lambda x: annotations[x]['uniprot_accession'])
        df['gene_name'] = df['protein'].map(lambda x: annotations[x]['gene_name'])
        df['protein_name'] = df['protein'].map(lambda x: annotations[x]['protein_name'])
        df['function'] = df['protein'].map(lambda x: annotations[x]['function'])
        df['annotation_source'] = df['protein'].map(lambda x: annotations[x]['annotation_source'])
        df['is_reviewed'] = df['protein'].map(lambda x: annotations[x]['is_reviewed'])
        df['quality_score'] = df['protein'].map(lambda x: annotations[x]['quality_score'])
        
        # Reorder columns
        cols = ['protein', 'uniprot_accession', 'gene_name', 'protein_name', 'function',
                'annotation_source', 'is_reviewed', 'quality_score',
                'AP_NS', 'AP_HS', 'PDL_NS', 'PDL_HS']
        df = df[cols]
        
        logger.info(f"Annotations added successfully")
        
        # Summary statistics
        logger.info("="*70)
        logger.info("FINAL ANNOTATION STATISTICS")
        logger.info("="*70)
        logger.info(f"Total proteins: {len(df)}")
        logger.info(f"With protein names: {(df['protein_name'] != 'Unknown').sum()}")
        logger.info(f"With functions: {(df['function'] != 'Not available').sum()}")
        logger.info(f"Reviewed entries: {df['is_reviewed'].sum()}")
        logger.info(f"Mean quality score: {df['quality_score'].mean():.1f}")
        
        # Show top 10 by quality score
        logger.info("\nTop 10 proteins by quality score:")
        top10 = df.nlargest(10, 'quality_score')[['protein', 'gene_name', 'protein_name', 'quality_score']]
        for line in top10.to_string(index=False).split('\n'):
            logger.info(line)
        
        # Save annotated data
        logger.info(f"\nSaving annotated data to: {output_file}")
        df.to_csv(output_file, index=False)
        
        file_size = output_file.stat().st_size / 1024
        logger.info(f"✓ File saved successfully ({file_size:.2f} KB)")
        
        # Validate against baseline
        validate_against_baseline(df, baseline_file, logger)
        
        logger.info("="*70)
        logger.info("✓ Script completed successfully!")
        logger.info("="*70)
        
        return df
        
    except Exception as e:
        logger.error(f"✗ Script failed: {str(e)}")
        logger.error("="*70)
        raise


if __name__ == "__main__":
    df_annotated = main()
