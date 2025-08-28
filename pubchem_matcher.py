import pandas as pd

def match_pubchem_data(input_file, output_file, pubchem_file):
    """
    Match MSP data with PubChem annotations based on InChIKey or Smiles,
    and insert PubMed_Count and Patent_Count before 'Num peaks' line.
    """
    pubchem_df = pd.read_csv(pubchem_file)
    # pubchem_df.set_index("InChIKey", inplace=True)
    pubchem_df.set_index("FirstBlock", inplace=True)

    # Initialize statistics
    stats = {
        "total_spectra": 0,
        "matched_by_inchikey": 0,
        "matched_by_smiles": 0,
        "total_matched": 0,
        "spectra_with_inchikey": 0,
        "spectra_with_smiles": 0,
        "spectra_with_both": 0,
        "spectra_with_neither": 0,
        "pubmed_counts": [],
        "patent_counts": []
    }

    with open(input_file, "r", encoding="utf-8") as fi, \
         open(output_file, "w", encoding="utf-8") as fo:

        spectrum_block = []
        inchi_key = None
        smiles = None

        for line in fi:
            spectrum_block.append(line)
            if line.startswith("Inchikey:") or line.startswith("InChIKey:"):
                inchi_key = line.split(":", 1)[1].strip()
            elif line.startswith("Smiles:"):
                smiles = line.split(":", 1)[1].strip()

            if line.strip() == "":
                stats["total_spectra"] += 1
                
                # Track identifier availability
                has_inchikey = bool(inchi_key)
                has_smiles = bool(smiles)
                
                if has_inchikey and has_smiles:
                    stats["spectra_with_both"] += 1
                elif has_inchikey:
                    stats["spectra_with_inchikey"] += 1
                elif has_smiles:
                    stats["spectra_with_smiles"] += 1
                else:
                    stats["spectra_with_neither"] += 1
                
                # Process entire block
                pubchem_row = None
                
                if inchi_key and inchi_key in pubchem_df.index:
                    pubchem_row = pubchem_df.loc[inchi_key]
                    stats["matched_by_inchikey"] += 1
                elif smiles:
                    match = pubchem_df[pubchem_df["SMILES"] == smiles]
                    if not match.empty:
                        pubchem_row = match.iloc[0]
                        stats["matched_by_smiles"] += 1

                if pubchem_row is not None:
                    stats["total_matched"] += 1
                    stats["pubmed_counts"].append(pubchem_row['PubMed_Count'])
                    stats["patent_counts"].append(pubchem_row['Patent_Count'])
                    
                    for idx, line_content in enumerate(spectrum_block):
                        if line_content.lower().startswith("num peaks"):
                            spectrum_block.insert(
                                idx, 
                                f"PubMed_Count: {pubchem_row['PubMed_Count']}\n"
                            )
                            spectrum_block.insert(
                                idx + 1, 
                                f"Patent_Count: {pubchem_row['Patent_Count']}\n"
                            )
                            break

                fo.writelines(spectrum_block)
                spectrum_block = []
                inchi_key = None
                smiles = None

    # Calculate and print detailed statistics
    print(f"\n=== PubChem Annotation Statistics ===")
    print(f"Input file: {input_file}")
    print(f"Output file: {output_file}")
    print(f"PubChem database: {pubchem_file}")
    print(f"PubChem database size: {len(pubchem_df):,} records")
    print()
    
    print(f"Total spectra processed: {stats['total_spectra']:,}")
    print(f"Total matched: {stats['total_matched']:,} ({stats['total_matched']/stats['total_spectra']*100:.1f}%)")
    print()
    
    print(f"Identifier availability:")
    print(f"  - InChIKey only: {stats['spectra_with_inchikey']:,}")
    print(f"  - SMILES only: {stats['spectra_with_smiles']:,}")
    print(f"  - Both InChIKey & SMILES: {stats['spectra_with_both']:,}")
    print(f"  - Neither: {stats['spectra_with_neither']:,}")
    print()
    
    print(f"Matching breakdown:")
    print(f"  - Matched by InChIKey: {stats['matched_by_inchikey']:,}")
    print(f"  - Matched by SMILES: {stats['matched_by_smiles']:,}")
    print(f"  - Unmatched: {stats['total_spectra'] - stats['total_matched']:,}")
    print()
    
    if stats['pubmed_counts']:
        pubmed_array = pd.Series(stats['pubmed_counts'])
        patent_array = pd.Series(stats['patent_counts'])
        
        print(f"PubChem Count statistics:")
        print(f"  - Min: {pubmed_array.min()}")
        print(f"  - Max: {pubmed_array.max()}")
        print(f"  - Mean: {pubmed_array.mean():.1f}")
        print(f"  - Median: {pubmed_array.median():.1f}")
        print()
        
        print(f"Patent Count statistics:")
        print(f"  - Min: {patent_array.min()}")
        print(f"  - Max: {patent_array.max()}")
        print(f"  - Mean: {patent_array.mean():.1f}")
        print(f"  - Median: {patent_array.median():.1f}")
        print()
    
    print(f"Annotation complete. Annotated MSP saved to {output_file}")
    
    return stats

if __name__ == "__main__":
    # Use the config-based paths for EFS storage
    base_path = "/reference_library/raw_msp"
    
    input_file = f"{base_path}/input_test.msp"
    output_file = f"{base_path}/output_test.msp" 
    pubchem_file = f"{base_path}/PubChemLite.csv"

    match_pubchem_data(input_file, output_file, pubchem_file)
