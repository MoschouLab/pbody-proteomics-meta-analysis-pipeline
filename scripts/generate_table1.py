#!/usr/bin/env python3
"""
Generate Table 1 for Review Manuscript
Matches Box1 formatting style (blue header)

Postdoctoral Bioinformatician Author: Moya-Cuevas J
Date: 2026-01-23
"""

from docx import Document
from docx.shared import Pt, RGBColor, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

def shade_cell(cell, fill_color):
    """Add blue background to table cell"""
    shading_elm = OxmlElement('w:shd')
    shading_elm.set(qn('w:fill'), fill_color)
    cell._element.get_or_add_tcPr().append(shading_elm)

def create_table1():
    """Create Table 1 with Box1 blue header style"""
    
    doc = Document()
    
    # Title
    title = doc.add_paragraph()
    title_run = title.add_run('Table 1. Summary of P-body proteomics enrichment, Gene Ontology analysis, and overlap patterns across experimental conditions.')
    title_run.bold = True
    title_run.font.size = Pt(11)
    
    # Headers
    headers = ['Category', 'Condition', 'Enriched\nGenes (n)', 'Cell Wall\nGenes (n)', 
               'Fold\nEnrichment', 'P-value /\nFDR', 'Significant GO Terms\n(FDR<0.05)', 
               'Top GO Categories', 'Stress Response']
    
    # Data
    data = [
        ['Proteome-wide\nenrichment', 'AP_NS', '299', '35', '2.44× (CW)', '4.6×10⁻⁷', 
         '2 (2 BP)', 'Protein phosphorylation;\nResponse to cadmium', 
         '262 lost;\n37 stable;\n120 gained'],
        ['', 'AP_HS', '157', '12', '1.59× (CW)', '0.073 (NS)', '7 (4 BP, 3 CC)', 
         'Ribosome biogenesis;\nNucleolus', '—'],
        ['', 'PDL_NS', '582', '53', '1.88× (CW)', '2.1×10⁻⁶', '63 (38 BP, 10 MF, 15 CC)', 
         'Apoplast (3.80×);\nCytosol (1.76×);\nCopper binding (3.40×)', 
         '537 lost;\n46 stable;\n150 gained'],
        ['', 'PDL_HS', '196', '4', '0.44× (CW)', '0.98\n(depleted)', '32 (24 BP, 3 MF, 5 CC)', 
         'Lsm1-7-Pat1 complex (21.5×);\nmRNA binding (2.29×);\nRNA processing', '—'],
        ['Overlap\nanalysis', 'AP_NS ∩ AP_HS', '37', '—', '3.76×', '3.2×10⁻¹²', '—', 
         'Method-specific core\n(stable across stress)', '—'],
        ['', 'PDL_NS ∩ PDL_HS', '46', '—', '1.93×', '1.6×10⁻⁵', '—', 
         'Method-specific core\n(stable across stress)', '—'],
        ['', 'All 4 conditions\n(CORE)', '1', '—', '—', '—', '—', 
         'AT1G08370: mRNA-decapping enzyme\n(universal P-body marker)', '—'],
        ['', 'Total unique', '1,095', '88', '—', '—', '—', '—', 
         'Net P-body depletion\nunder heat stress']
    ]
    
    # Create table
    table = doc.add_table(rows=1 + len(data), cols=len(headers))
    table.style = 'Table Grid'
    
    # Header row - BLUE BACKGROUND like Box1
    header_cells = table.rows[0].cells
    for i, header_text in enumerate(headers):
        cell = header_cells[i]
        paragraph = cell.paragraphs[0]
        run = paragraph.add_run(header_text)
        
        # WHITE text on BLUE background (Box1 style)
        run.font.size = Pt(10)
        run.font.bold = True
        run.font.color.rgb = RGBColor(255, 255, 255)  # White
        paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
        
        shade_cell(cell, '4472C4')  # Blue background
    
    # Data rows
    for row_idx, row_data in enumerate(data, start=1):
        row_cells = table.rows[row_idx].cells
        
        for col_idx, cell_text in enumerate(row_data):
            cell = row_cells[col_idx]
            paragraph = cell.paragraphs[0]
            run = paragraph.add_run(cell_text)
            
            run.font.size = Pt(9)
            
            # Bold category names
            if col_idx == 0 and cell_text.strip():
                run.font.bold = True
            
            # Center numeric columns
            if col_idx in [2, 3, 4, 5, 6]:
                paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    
    # Caption (italic)
    doc.add_paragraph()
    caption = doc.add_paragraph()
    caption_text = (
        "P-body-enriched proteins were identified using affinity purification (AP) and "
        "proximity-dependent labeling (PDL) under no stress (NS) and heat stress (HS) "
        "conditions at log₂FC > 0.8. Cell wall-associated genes were identified using "
        "Fisher's exact test against GO-annotated background (n=4,593 genes). GO enrichment "
        "analysis was performed using hypergeometric test with FDR correction (Benjamini-Hochberg). "
        "Overlap significance was assessed using hypergeometric test with FDR correction. "
        "Fold enrichment values represent observed/expected ratios. Stress response indicates "
        "gene count changes between NS and HS conditions within each method. BP, Biological Process; "
        "MF, Molecular Function; CC, Cellular Component; CW, Cell Wall; NS, not significant (FDR ≥ 0.05). "
        "AT1G08370 encodes an mRNA-decapping enzyme-like protein, representing the only gene enriched "
        "across all four experimental conditions. The net P-body depletion under heat stress "
        "(799 genes lost vs 270 gained across both methods) suggests stress-induced disassembly "
        "as a cellular response mechanism."
    )
    
    caption_run = caption.add_run(caption_text)
    caption_run.font.size = Pt(9)
    caption_run.font.italic = True
    
    # Save
    doc.save('Table1_PBody_Summary.docx')
    print("✅ Table 1 created: Table1_PBody_Summary.docx")
    print("   Style: Blue header with white text (Box1 style)")

if __name__ == "__main__":
    create_table1()
