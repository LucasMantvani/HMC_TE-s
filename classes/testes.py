import requests

# Endpoint da API para detalhes de uma família específica
url = 'https://dfam.org/api/families/DF000001643/assemblies'
url = 'https://dfam.org/api/classes'
url = 'https://dfam.org/api/version'
url = 'https://dfam.org/api/families/DF000000012'




response = requests.get(url)

# Verificando se a requisição foi bem-sucedida
if response.status_code == 200:

    try:
        data = response.json()
        print(data)

    except:

        print(response.text)
else:
    print('Erro:', response.status_code)
