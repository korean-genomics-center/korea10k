#!/usr/bin/env python3
import os
import shutil
import subprocess
import sys
import tempfile


def run_cmd(cmd, check=True, capture_output=False):
    """Run a shell command safely and return output if requested."""
    result = subprocess.run(
        cmd,
        shell=True,
        text=True,
        check=check,
        capture_output=capture_output
    )
    if capture_output:
        return result.stdout.strip()
    return None


def fix_vcf_header(input_vcf, ref_fai, output_vcf):
    # Ensure bcftools is available
    if shutil.which("bcftools") is None:
        sys.exit("[ERROR] bcftools not found in PATH. Please install it first.")

    # --- 1. Extract current VCF header ---
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as tmp_header:
        header_file = tmp_header.name
        run_cmd(f"bcftools view -h {input_vcf} > {header_file}")

    # --- 2. Check if contig definitions already exist ---
    with open(header_file) as f:
        has_contigs = any(line.startswith("##contig=") for line in f)

    if has_contigs:
        print("[✓] Contig definitions already exist — copying file.")
        shutil.copy2(input_vcf, output_vcf)
        return

    print(f"[!] No contig definitions found. Generating from: {ref_fai}")

    # --- 3. Generate new contig header lines from .fai ---
    contig_file = tempfile.NamedTemporaryFile(mode='w+', delete=False)
    contig_path = contig_file.name
    with open(ref_fai) as fai:
        for line in fai:
            fields = line.strip().split('\t')
            if len(fields) < 2:
                continue
            chrom, length = fields[0], fields[1]
            contig_file.write(f"##contig=<ID={chrom},length={length}>\n")
    contig_file.close()

    # --- 4. Insert contigs before the #CHROM line in header ---
    fixed_header = tempfile.NamedTemporaryFile(mode='w+', delete=False)
    fixed_header_path = fixed_header.name

    with open(header_file) as orig, open(contig_path) as contigs, open(fixed_header_path, 'w') as out:
        contig_lines = contigs.read()
        for line in orig:
            if line.startswith("#CHROM"):
                out.write(contig_lines)
            out.write(line)

    # --- 5. Reheader the VCF using bcftools ---
    run_cmd(f"bcftools reheader -h {fixed_header_path} -o {output_vcf} {input_vcf}")

    print(f"[+] Fixed VCF saved to: {output_vcf}")

    # --- 6. Cleanup temporary files ---
    for fpath in [header_file, contig_path, fixed_header_path]:
        try:
            os.remove(fpath)
        except OSError:
            pass


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python fix_vcf_missing_contigs.py <input.vcf(.gz)> <reference.fai> <output_fixed.vcf(.gz)>")
        sys.exit(1)

    input_vcf, ref_fai, output_vcf = sys.argv[1:4]

    if not os.path.exists(input_vcf):
        sys.exit(f"[ERROR] Input VCF not found: {input_vcf}")
    if not os.path.exists(ref_fai):
        sys.exit(f"[ERROR] Reference FAI not found: {ref_fai}")

    fix_vcf_header(input_vcf, ref_fai, output_vcf)
