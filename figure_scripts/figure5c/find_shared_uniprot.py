from bioservices import UniProt

def extract_uniprot_ids_and_queries(diamond_file):
    """Return:
    - set of UniProt IDs
    - dict mapping UniProt ID -> first query ID seen
    """
    uniprot_ids = set()
    query_map = {}
    with open(diamond_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 2:
                continue
            query_id = fields[0]
            subject_id = fields[1]
            if subject_id.startswith("sp|"):
                parts = subject_id.split('|')
                if len(parts) >= 2:
                    uid = parts[1]
                    uniprot_ids.add(uid)
                    if uid not in query_map:
                        query_map[uid] = query_id  # store only first match
    return uniprot_ids, query_map

def fetch_info(uid, service):
    try:
        entry_txt = service.retrieve(uid, frmt="txt")
        lines = entry_txt.splitlines()

        protein_name = "Unknown"
        organism = "Unknown"

        for line in lines:
            if line.startswith("DE   RecName: Full="):
                protein_name = line.split("=", 1)[1].rstrip(";")
            elif line.startswith("OS   "):
                organism = line[5:].strip()

        return protein_name, organism
    except Exception as e:
        return "Error", str(e)

def main():
    file1 = "scaffold_GCA_017434555.proteins.matches.tsv"
    file2 = "scaffold_GCA_017406775.proteins.matches.tsv"

    ids1, queries1 = extract_uniprot_ids_and_queries(file1)
    ids2, queries2 = extract_uniprot_ids_and_queries(file2)
    shared_ids = sorted(ids1 & ids2)

    output_file = "shared_uniprot_output.tsv"
    with open(output_file, "w") as out:
        out.write("UniProt ID\tProtein Name\tOrganism\tQuery ID File 1\tQuery ID File 2\n")
        
        print("\nShared UniProt IDs and Info:\n")
        print("{:<12}  {:<50}  {:<25}  {:<20}  {}".format("UniProt ID", "Protein Name", "Organism", "File 1", "File 2"))
        print("-" * 130)

        u = UniProt()
        for uid in shared_ids:
            protein, organism = fetch_info(uid, u)
            query1 = queries1.get(uid, "")
            query2 = queries2.get(uid, "")
            print("{:<12}  {:<50}  {:<25}  {:<20}  {}".format(uid, protein, organism, query1, query2))
            out.write(f"{uid}\t{protein}\t{organism}\t{query1}\t{query2}\n")

    print("\n✅ Results saved to '{}'".format(output_file))

if __name__ == "__main__":
    main()
