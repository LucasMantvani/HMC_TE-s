import pandas as pd
import numpy as np

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# Sua sequência de exemplo
sequence = "TCACCCCCAGCCCTCAGGGACTCTCCTGTCTGTCCCAGCTACTCTCCAACCACGCCCAGATTTCAGCNGGAGTCAGTTCCAGGCACCCAGGAATCACCNNCAAACTCACCGATTTCACTGAGTTACTCCCCAGTCTCTCTCACGTCT"

# Realizar uma busca de domínios conservados usando BLAST no NCBI
result_handle = NCBIWWW.qblast("blastn", "nt", sequence)

# Parsear os resultados do BLAST
blast_records = NCBIXML.read(result_handle)

# Exibir resultados
for alignment in blast_records.alignments:
    for hsp in alignment.hsps:
        print("****Alignment****")
        print("sequence:", alignment.title)
        print("length:", alignment.length)
        print("e value:", hsp.expect)
        print(hsp.query[0:75] + "...")
        print(hsp.match[0:75] + "...")
        print(hsp.sbjct[0:75] + "...")

