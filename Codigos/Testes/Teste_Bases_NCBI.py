from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

NCBIWWW.email = 'jbaioco@usp.br'

# Sua sequência de exemplo
sequence = "TANTTAAGGGATAATGTTCACGGCGGAGGAGTATACGAAGCAATAAACGGCTTTTGCGGGTGATTAAACGCCGAAGTGAAGCTGAGGCGTTTGATNAACCGCAAAAGCCGTTTATTGCGAGTATACTCCAATGCCACGAACATTATTCCTATTACACGACAAGANAGA"

# Realizar uma busca de domínios conservados usando BLAST no NCBI
result_handle = NCBIWWW.qblast("blastn", "nt", sequence)

# Parsear os resultados do BLAST
blast_records = NCBIXML.read(result_handle)

print(type(blast_records))

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

