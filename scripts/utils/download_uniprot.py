#download_uniprot.py
#script para descargar secuencias de proteinas de UniProt con biopython
import requests

# URL de la secuencia en formato FASTA
url = 'https://www.uniprot.org/uniprotkb/Q9Y238.fasta'

# Ruta donde se guardará el archivo
output_path = '/data/UniProt/Q9Y238.fasta'

# Hacer la solicitud HTTP y guardar el archivo
response = requests.get(url)

# Verificar que la descarga fue exitosa
if response.status_code == 200:
    with open(output_path, 'wb') as file:
        file.write(response.content)
    print(f"Archivo descargado y guardado en {output_path}")
else:
    print(f"Error en la descarga: {response.status_code}")
