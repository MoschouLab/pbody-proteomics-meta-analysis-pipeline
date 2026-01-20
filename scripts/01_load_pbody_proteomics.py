#!/usr/bin/env python3
"""
Script 01: Load P-body Proteomics Data
========================================

Purpose:
    Load P-body proteomics data from EMBOJ-2022-111885 Source File 1.
    Extracts protein IDs and log2FC values for:
    - Affinity Purification (AP): NS and HS
    - Proximity-Dependent Labeling (PDL): NS and HS

Reference:
    Liu C, Mentzelopoulou A, Muhammad A, et al. An actin remodeling role for 
    Arabidopsis processing bodies revealed by their proximity interactome. 
    EMBO J. 2023;42(9):e111885. doi:10.15252/embj.2022111885

Input:
    - data/raw/EMBOJ-2022-111885_SourceFile1 .xlsx

Output:
    - data/processed/pbody_proteomics_data.csv
    - logs/01_load_pbody_proteomics.log

Author: Pepe (migrated from Colab)
Date: 2026-01-20
Version: 1.1 (column search by name)
"""

import pandas as pd
import sys
import logging
from pathlib import Path
from datetime import datetime


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
    log_file = log_dir / '01_load_pbody_proteomics.log'
    
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
    logger.info("P-Body Proteomics Data Loader - Script 01 v1.1")
    logger.info("="*70)
    logger.info(f"Execution started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Log file: {log_file}")
    
    return logger


def validate_input_file(input_file, logger):
    """
    Validate that input file exists and log file information.
    
    Parameters:
    -----------
    input_file : Path
        Path to input Excel file
    logger : logging.Logger
        Logger instance
    """
    if not input_file.exists():
        logger.error(f"Input file not found: {input_file}")
        raise FileNotFoundError(f"Input file not found: {input_file}")
    
    logger.info(f"Input file found: {input_file}")
    file_size = input_file.stat().st_size / 1024  # KB
    logger.info(f"File size: {file_size:.2f} KB")


def find_column_indices(df_raw, logger):
    """
    Find column indices by searching multi-level headers.
    
    The Excel has a 3-level header structure:
    - Row 0: "LOG2 RATIO DCP1/GFP" (spanning cols 32-35)
    - Row 1: Method - "AP" or "PDL" (may have nans for continuation)
    - Row 2: Condition - "Log2(FC) NS" or "Log2(FC) HS"
    
    Parameters:
    -----------
    df_raw : pd.DataFrame
        Raw data from Excel
    logger : logging.Logger
        Logger instance
        
    Returns:
    --------
    dict
        Dictionary mapping column names to indices
    """
    logger.info("Searching for column indices in multi-level headers...")
    
    # Extract header rows
    row1 = df_raw.iloc[1, :].copy()  # Method row (AP/PDL)
    row2 = df_raw.iloc[2, :].copy()  # Condition row (NS/HS)
    
    # Forward fill row1 to handle merged cells (nans)
    row1_filled = row1.ffill()
    
    logger.info("Header structure (after forward fill):")
    for idx in range(30, 36):
        logger.info(f"  Col {idx}: Method='{row1_filled.iloc[idx]}', Condition='{row2.iloc[idx]}'")
    
    # Search for target columns
    column_map = {}
    
    for idx in range(len(row1_filled)):
        method = str(row1_filled.iloc[idx]).strip()
        condition = str(row2.iloc[idx]).strip()
        
        # AP Log2(FC) NS
        if method == 'AP' and 'NS' in condition and 'Log2(FC)' in condition:
            column_map['AP_NS'] = idx
            logger.info(f"  ✓ Found AP_NS at column index {idx}")
        
        # AP Log2(FC) HS
        elif method == 'AP' and 'HS' in condition and 'Log2(FC)' in condition:
            column_map['AP_HS'] = idx
            logger.info(f"  ✓ Found AP_HS at column index {idx}")
        
        # PDL Log2(FC) NS
        elif method == 'PDL' and 'NS' in condition and 'Log2(FC)' in condition:
            column_map['PDL_NS'] = idx
            logger.info(f"  ✓ Found PDL_NS at column index {idx}")
        
        # PDL Log2(FC) HS
        elif method == 'PDL' and 'HS' in condition and 'Log2(FC)' in condition:
            column_map['PDL_HS'] = idx
            logger.info(f"  ✓ Found PDL_HS at column index {idx}")
    
    # Validate that all columns were found
    expected_cols = ['AP_NS', 'AP_HS', 'PDL_NS', 'PDL_HS']
    missing = [col for col in expected_cols if col not in column_map]
    
    if missing:
        logger.error(f"Failed to find columns: {missing}")
        raise ValueError(f"Could not locate columns: {missing}")
    
    logger.info(f"✓ All {len(expected_cols)} columns found successfully")
    return column_map


def load_pbody_data(input_file, logger):
    """
    Load P-body proteomics data with dynamic column search.
    
    Parameters:
    -----------
    input_file : Path
        Path to input Excel file
    logger : logging.Logger
        Logger instance
    
    Returns:
    --------
    pd.DataFrame
        DataFrame with columns: protein, AP_NS, AP_HS, PDL_NS, PDL_HS
    """
    logger.info("Reading Excel file...")
    
    try:
        # Read the file
        df_raw = pd.read_excel(input_file, sheet_name=0, header=None)
        logger.info(f"Raw data shape: {df_raw.shape}")
        
        # Find column indices dynamically
        col_map = find_column_indices(df_raw, logger)
        
        # Extract protein column (column A, index 0) - skip first 3 header rows
        proteins = df_raw.iloc[3:, 0].reset_index(drop=True)
        logger.info(f"Extracted {len(proteins)} protein entries from column A (index 0)")
        
        # Extract log2FC columns using found indices
        logger.info("Extracting log2FC data using found column indices...")
        result = pd.DataFrame({
            'protein': proteins,
            'AP_NS': df_raw.iloc[3:, col_map['AP_NS']].reset_index(drop=True),
            'AP_HS': df_raw.iloc[3:, col_map['AP_HS']].reset_index(drop=True),
            'PDL_NS': df_raw.iloc[3:, col_map['PDL_NS']].reset_index(drop=True),
            'PDL_HS': df_raw.iloc[3:, col_map['PDL_HS']].reset_index(drop=True)
        })
        
        # Convert log2FC columns to numeric (handle any non-numeric values)
        logger.info("Converting log2FC columns to numeric...")
        for col in ['AP_NS', 'AP_HS', 'PDL_NS', 'PDL_HS']:
            result[col] = pd.to_numeric(result[col], errors='coerce')
            n_valid = result[col].notna().sum()
            n_missing = result[col].isna().sum()
            logger.info(f"  {col}: {n_valid} valid values, {n_missing} missing")
        
        # Reset index
        result = result.reset_index(drop=True)
        
        logger.info(f"Final dataset shape: {result.shape}")
        return result
        
    except Exception as e:
        logger.error(f"Error loading data: {str(e)}")
        raise


def validate_data(df, logger):
    """
    Validate loaded data dimensions and value ranges.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Loaded proteomics data
    logger : logging.Logger
        Logger instance
    """
    logger.info("="*70)
    logger.info("DATA VALIDATION")
    logger.info("="*70)
    
    # Validate dimensions
    n_proteins = len(df)
    logger.info(f"Number of proteins: {n_proteins}")
    
    # Expected ~4776 proteins based on Excel info
    if n_proteins < 4000:
        logger.warning(f"⚠️  Low protein count: {n_proteins} (expected ~4776)")
    elif n_proteins > 5000:
        logger.warning(f"⚠️  High protein count: {n_proteins} (expected ~4776)")
    else:
        logger.info(f"✓ Protein count within expected range (~4776)")
    
    # Validate log2FC ranges (typically -10 to +10 for proteomics)
    logger.info("\nValidating log2FC value ranges:")
    for col in ['AP_NS', 'AP_HS', 'PDL_NS', 'PDL_HS']:
        valid_data = df[col].dropna()
        if len(valid_data) > 0:
            min_val = valid_data.min()
            max_val = valid_data.max()
            logger.info(f"  {col}: [{min_val:.2f}, {max_val:.2f}]")
            
            # Check for extreme values
            if abs(min_val) > 15 or abs(max_val) > 15:
                logger.warning(f"    ⚠️  Extreme values detected (|log2FC| > 15)")
            else:
                logger.info(f"    ✓ Values within typical range")
        else:
            logger.warning(f"  {col}: No valid data")
    
    # Check for completely missing proteins
    n_missing_protein_ids = df['protein'].isna().sum()
    if n_missing_protein_ids > 0:
        logger.warning(f"⚠️  {n_missing_protein_ids} rows with missing protein IDs")


def display_summary(df, logger):
    """
    Display and log summary statistics for the loaded data.
    
    Parameters:
    -----------
    df : pd.DataFrame
        The loaded proteomics data
    logger : logging.Logger
        Logger instance
    """
    logger.info("="*70)
    logger.info("DATA SUMMARY")
    logger.info("="*70)
    
    logger.info(f"Total proteins: {len(df)}")
    logger.info(f"Proteins with valid IDs: {df['protein'].notna().sum()}")
    
    logger.info("\nSummary statistics for log2FC values:")
    summary = df[['AP_NS', 'AP_HS', 'PDL_NS', 'PDL_HS']].describe()
    for line in summary.to_string().split('\n'):
        logger.info(line)
    
    logger.info("\nMissing values per column:")
    missing = df.isnull().sum()
    for col, count in missing.items():
        logger.info(f"  {col}: {count}")
    
    logger.info("\nFirst 10 entries:")
    for line in df.head(10).to_string(index=False).split('\n'):
        logger.info(line)
    
    logger.info("\nLast 10 entries:")
    for line in df.tail(10).to_string(index=False).split('\n'):
        logger.info(line)


def save_data(df, output_file, logger):
    """
    Save processed data to CSV.
    
    Parameters:
    -----------
    df : pd.DataFrame
        The processed proteomics data
    output_file : Path
        Path to output CSV file
    logger : logging.Logger
        Logger instance
    """
    logger.info(f"Saving processed data to: {output_file}")
    df.to_csv(output_file, index=False)
    
    # Verify save
    if output_file.exists():
        file_size = output_file.stat().st_size / 1024  # KB
        n_lines = len(df) + 1  # +1 for header
        logger.info(f"✓ File saved successfully")
        logger.info(f"  Size: {file_size:.2f} KB")
        logger.info(f"  Lines: {n_lines} (including header)")
    else:
        logger.error("✗ Failed to save output file")
        raise IOError("Failed to save output file")


def main():
    """Main execution function."""
    
    # Get project root (parent of scripts directory)
    project_dir = Path(__file__).resolve().parent.parent
    
    # Define paths
    data_raw = project_dir / 'data' / 'raw'
    data_processed = project_dir / 'data' / 'processed'
    log_dir = project_dir / 'logs'
    
    # Input and output files
    input_file = data_raw / 'EMBOJ-2022-111885_SourceFile1 .xlsx'
    output_file = data_processed / 'pbody_proteomics_data.csv'
    
    # Setup logging
    logger = setup_logging(log_dir)
    logger.info(f"Project directory: {project_dir}")
    logger.info("Version: 1.1 (dynamic column search with forward fill)")
    
    try:
        # Validate input
        validate_input_file(input_file, logger)
        
        # Load data
        df = load_pbody_data(input_file, logger)
        
        # Validate data
        validate_data(df, logger)
        
        # Display summary
        display_summary(df, logger)
        
        # Save processed data
        save_data(df, output_file, logger)
        
        logger.info("="*70)
        logger.info("✓ Script completed successfully!")
        logger.info("="*70)
        
        return df
        
    except Exception as e:
        logger.error(f"✗ Script failed: {str(e)}")
        logger.error("="*70)
        raise


if __name__ == "__main__":
    df = main()
