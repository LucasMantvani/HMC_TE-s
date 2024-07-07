import requests
import time

a = time.time()



url = "https://dfam.org/api/families"

params = {
    "format": "full",

    "start": 10,

    "limit": 1
    }

response = requests.get(url, params=params)
results = response.json()

# Prints "Vingi-2_CE" at the time of this writing

b = time.time()

print(b-a)
# print(results)

print(results)