import requests 
import pandas as pd

import time

class Data:

    def __init__(self, taxa_dowload: int = 1000) -> None:

        self.url:          str  = 'https://dfam.org/api/families'
        self.taxa_dowload: int  = taxa_dowload
        self.params:       dict = {"format": "full", 
                                   "start" : 0, 
                                   "limit" : self.taxa_dowload}

        return None
    
    def get_data(self) -> pd.DataFrame:

        resposta: requests = requests.get(self.url, self.params)

        try:

            resultados: dict = resposta.json()['results']
            
        except:

            raise ValueError("O formato da entrada não é compativel com jason, veja se eles não mudaram a API.")
        

        saida: pd.DataFrame = pd.DataFrame(resultados)[['accession', 'classification', 'consensus_sequence']]

        self.params["start"] += self.taxa_dowload
        
        return saida


def main() -> None:

    a = time.time()

    data = Data(1000)
    # data = Data(5)

    data.get_data().to_csv("C:\\Users\\UsuarioOct\\Documentos\\trabalhos\\Em produção\\IC\\code\\temp_data\\dado.csv")

    b = time.time()

    print(b-a)

    return None

if __name__ == '__main__':

    main()