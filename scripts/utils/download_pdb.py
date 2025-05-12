# download_pdb.py
# Este script descarga estructuras de proteínas desde el PDB utilizando Biopython.

import os
from Bio.PDB import PDBList

# Ruta donde se guardará el archivo PDB
output_path = '/Users/Olivermop/Documents/bioinformatics_portfolio/data/PDB'

# Crear la carpeta PDB si no existe
if not os.path.exists(output_path):
    os.makedirs(output_path)

# ID de la proteína en el PDB (puedes cambiar este ID a la proteína que desees)
pdb_id = '1TUP'

# Descargar la estructura PDB usando Biopython
pdbl = PDBList()
print(f"Descargando estructura PDB '{pdb_id}'...")
pdb_file = pdbl.retrieve_pdb_file(pdb_id, pdir=output_path, file_format='pdb')

# Cambiar el nombre del archivo a la extensión correcta (de .ent a .pdb)
pdb_file_corrected = pdb_file.replace(".ent", ".pdb")
os.rename(pdb_file, pdb_file_corrected)

print(f"Estructura PDB '{pdb_id}' descargada y guardada en {pdb_file_corrected}")
