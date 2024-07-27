import pandas as pd 
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def file_fasta(sequencias: pd.DataFrame, file_path: str) -> None:

    sequencias = sequencias.set_index('accession').to_dict()['consensus_sequence']

    with open(file_path, "w") as fasta_file:

        for seq_id, sequence in sequencias.items():

            fasta_file.write(f">{seq_id}\n")
            fasta_file.write(f"{sequence}\n")

    print(f"Arquivo {file_path} criado com sucesso!")

    return None


def new_file_fasta(sequencias: pd.DataFrame, file_path: str) -> None:
   
    a = [SeqRecord(Seq(row.consensus_sequence), id=row.accession) for row in sequencias.itertuples(index=False)]

    SeqIO.write(a, file_path, "fasta")

    return None


def main() -> None:

    #disposições gerais

    diretorio = "/home/mantovani/Documents/trabalhos/Em_produção/IC/HMC_TE-s/Temp_Data/Arquivos_Fasta/"

    #reduzindo o data_frame

    df = pd.read_csv('/home/mantovani/Documents/trabalhos/Em_produção/IC/HMC_TE-s/Temp_Data/Dados_Genes/dado.csv')

    df = df[['accession', 'consensus_sequence']]

    df = df.iloc[range(0,1)]

    
    # file_fasta(df, diretorio + "sequences.fasta")
    new_file_fasta(df, diretorio + "sequences2.fasta")

    return None


if __name__ == "__main__":

    main()