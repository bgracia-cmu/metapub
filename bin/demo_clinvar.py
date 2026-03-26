import tempfile
from metapub import ClinVarFetcher

# Demo with known valid ClinVar variant IDs
# Note: Most low-numbered IDs are invalid, so we use a curated list
VALID_IDS = [4, 8, 1013, 10000, 12000, 12003, 12004, 12005, 12006, 12007]

# Use a temporary cache directory to avoid conflicts
with tempfile.TemporaryDirectory() as tmpdir:
    cvfetch = ClinVarFetcher(cachedir=tmpdir)

    # Search directly by gene:
    # for variant in cvfetch.variants_by_gene("CFTR"):
    #     print("Variant:", variant)

    print("ClinVar Fetcher Demo - showing valid variants")
    print("=" * 50)
    hgvs_list = []

    for varid in VALID_IDS:
        print(f"Variant ID: {varid}")
        try:
            var = cvfetch.variant(varid)
            hgvs_list.append(var.hgvs_c[0])
        except Exception as error:
            print(f"  ERROR: {error}")

        print()

    print("ClinVar demo completed.")
    print("LIST:", hgvs_list)

    # Now pretend these are HGVS sequences not associated with a variant
    variant_list = cvfetch.annotate_variants(hgvs_list)
    print("Variants", variant_list)

