import bibtexparser

def load(path):
    with open(path) as f:
        db = bibtexparser.load(f)
    return {e['ID']: e for e in db.entries}

a = load("pubs/shared/references.bib")
b = load("pubs/shared/zotero-meem-refs.bib") # jfm.bib, vendored from OpenFLASH pubs/JFM

keys = set(a) | set(b)

num_dups = sum(1 for k in keys if k in a and k in b)
print(f"Found {num_dups} duplicate keys out of {len(keys)} total keys")

for k in sorted(keys):
    if k in a and k in b:
        print (f"\nDUPLICATE key: {k}")
        if a[k] != b[k]:
            print(f"\nDIFF in key: {k}")
            for field in sorted(set(a[k]) | set(b[k])):
                if a[k].get(field) != b[k].get(field):
                    print(f"  {field}:\n    file1: {a[k].get(field)}\n    file2: {b[k].get(field)}")
        else:
            print(f"  entries are identical for key: {k}")