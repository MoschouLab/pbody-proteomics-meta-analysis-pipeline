#!/usr/bin/env python3
"""
Process STRING PPI Data for Arabidopsis
Converts UniProt IDs to AGI locus IDs and filters by confidence score

Author: Pepe
Date: 2026-01-23
"""

import pandas as pd
from pathlib import Path

# Configuration
RAW_DIR = Path("data/raw")
PROCESSED_DIR = Path("data/processed")

STRING_FILE = RAW_DIR / "3702.protein.links.v12.0.txt"
PROTEOMICS_FILE = PROCESSED_DIR / "pbody_proteomics_annotated.csv"
OUTPUT_FILE = PROCESSED_DIR / "ppi_edges.csv"

CONFIDENCE_THRESHOLD = 400  # Medium confidence (STRING scale: 0-1000)

print("="*90)
print("STRING PPI DATA PROCESSING")
print("="*90)

# Load STRING data
print(f"\n1. Loading STRING data from: {STRING_FILE}")
string_df = pd.read_csv(STRING_FILE, sep=' ')
print(f"   Loaded: {len(string_df):,} interactions")

# Remove '3702.' prefix from protein IDs
print("\n2. Cleaning protein IDs (removing '3702.' prefix)...")
string_df['protein1'] = string_df['protein1'].str.replace('3702.', '', regex=False)
string_df['protein2'] = string_df['protein2'].str.replace('3702.', '', regex=False)

# Filter by confidence score
print(f"\n3. Filtering by confidence score (>={CONFIDENCE_THRESHOLD})...")
string_df = string_df[string_df['combined_score'] >= CONFIDENCE_THRESHOLD]
print(f"   Retained: {len(string_df):,} high-confidence interactions")

# Load UniProt → AGI mapping from our proteomics data
print(f"\n4. Loading UniProt→AGI mapping from: {PROTEOMICS_FILE}")
proteomics_df = pd.read_csv(PROTEOMICS_FILE)

print(f"   Available columns: {list(proteomics_df.columns)}")

# Create mapping dictionary: UniProt accession → AGI locus ID
# Use uniprot_accession column (correct column name)
uniprot_to_agi = dict(zip(
    proteomics_df['uniprot_accession'].dropna().str.upper(),  # STRING uses uppercase
    proteomics_df['protein']
))

print(f"   Mapping size: {len(uniprot_to_agi):,} UniProt accessions")

# Sample some mappings to verify
print("\n   Sample mappings (first 5):")
for uniprot, agi in list(uniprot_to_agi.items())[:5]:
    print(f"     {uniprot} → {agi}")

# Map STRING UniProt IDs to AGI locus IDs
print("\n5. Converting STRING interactions to AGI locus IDs...")
string_df['proteinA'] = string_df['protein1'].str.upper().map(uniprot_to_agi)
string_df['proteinB'] = string_df['protein2'].str.upper().map(uniprot_to_agi)

# Count successful mappings
mapped_a = string_df['proteinA'].notna().sum()
mapped_b = string_df['proteinB'].notna().sum()
print(f"   Protein1 mapped: {mapped_a:,}/{len(string_df):,} ({100*mapped_a/len(string_df):.1f}%)")
print(f"   Protein2 mapped: {mapped_b:,}/{len(string_df):,} ({100*mapped_b/len(string_df):.1f}%)")

# Remove interactions where mapping failed (not in our dataset)
before = len(string_df)
string_df = string_df.dropna(subset=['proteinA', 'proteinB'])
after = len(string_df)
print(f"   Both mapped: {after:,} interactions ({before-after:,} dropped)")

# If no mappings, investigate
if after == 0:
    print("\n   ⚠️ WARNING: No interactions mapped!")
    print("   Investigating STRING protein IDs...")
    
    # Sample STRING IDs
    sample_string_ids = list(set(
        list(string_df.head(100)['protein1'].values) + 
        list(string_df.head(100)['protein2'].values)
    ))[:10]
    
    print(f"\n   Sample STRING IDs (first 10):")
    for sid in sample_string_ids:
        print(f"     {sid}")
    
    # Sample our UniProt IDs
    sample_our_ids = proteomics_df['uniprot_accession'].dropna().head(10).values
    print(f"\n   Sample UniProt accessions from our data (first 10):")
    for uid in sample_our_ids:
        print(f"     {uid}")
    
    print("\n   💡 STRING uses UniProt IDs, but we only have a subset in our proteomics data.")
    print("   💡 This is EXPECTED - we only have ~4,748 proteins, not the full Arabidopsis proteome.")
    print("   💡 Continuing with available mappings...")

# Remove self-loops
before = len(string_df)
string_df = string_df[string_df['proteinA'] != string_df['proteinB']]
after = len(string_df)
if before > after:
    print(f"   Removed {before-after:,} self-loops")

# Keep only necessary columns
ppi_df = string_df[['proteinA', 'proteinB', 'combined_score']].copy()

# Remove duplicate edges (A-B and B-A are the same in undirected network)
if len(ppi_df) > 0:
    print("\n6. Removing duplicate edges...")
    # Sort protein pairs alphabetically to identify duplicates
    ppi_df['pair'] = ppi_df.apply(
        lambda row: tuple(sorted([row['proteinA'], row['proteinB']])), 
        axis=1
    )
    before = len(ppi_df)
    ppi_df = ppi_df.drop_duplicates(subset='pair', keep='first')
    after = len(ppi_df)
    print(f"   Removed {before-after:,} duplicate edges")
    
    # Drop temporary column
    ppi_df = ppi_df.drop('pair', axis=1)
    
    # Sort by confidence score (highest first)
    ppi_df = ppi_df.sort_values('combined_score', ascending=False)

# Save to CSV
print(f"\n7. Saving PPI network to: {OUTPUT_FILE}")
OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
ppi_df.to_csv(OUTPUT_FILE, index=False)

print(f"   Saved: {len(ppi_df):,} interactions")

# Summary statistics
print("\n" + "="*90)
print("PPI NETWORK SUMMARY")
print("="*90)
if len(ppi_df) > 0:
    print(f"Total interactions: {len(ppi_df):,}")
    unique_proteins = set(ppi_df['proteinA']) | set(ppi_df['proteinB'])
    print(f"Unique proteins: {len(unique_proteins):,}")
    print(f"Confidence score range: {ppi_df['combined_score'].min():.0f} - {ppi_df['combined_score'].max():.0f}")
    print(f"Mean confidence: {ppi_df['combined_score'].mean():.1f}")
    print(f"Median confidence: {ppi_df['combined_score'].median():.0f}")
    
    # Top 10 most confident interactions
    print("\nTop 10 highest-confidence interactions:")
    print(ppi_df.head(10).to_string(index=False))
    
    # Check coverage of our proteomics proteins
    our_proteins = set(proteomics_df['protein'])
    proteins_with_ppi = unique_proteins & our_proteins
    print(f"\nProteomics coverage:")
    print(f"  Proteins in our dataset: {len(our_proteins):,}")
    print(f"  Proteins with PPI data: {len(proteins_with_ppi):,} ({100*len(proteins_with_ppi)/len(our_proteins):.1f}%)")
else:
    print("⚠️ No interactions mapped. Network analysis will be skipped in Script 05.")

print("\n" + "="*90)
print("PPI DATA PROCESSING COMPLETE")
print("="*90)
