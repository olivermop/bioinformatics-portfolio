# download_vcf.py
# Este script descarga un archivo VCF específico del cromosoma 22 desde el Proyecto 1000 Genomas.
import os
import subprocess

# URL del archivo VCF del cromosoma 22 (del Proyecto 1000 Genomas)
vcf_url = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz'

# Directorio donde se almacenará el archivo VCF descargado
output_dir = '/Users/Olivermop/Documents/bioinformatics_portfolio/data/VCF'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Ruta completa para guardar el archivo descargado
output_path = os.path.join(output_dir, 'ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz')

# Comando wget para descargar el archivo VCF
wget_command = f"wget -O {output_path} {vcf_url}"

# Ejecutar el comando wget
try:
    print(f"Descargando archivo desde {vcf_url} usando wget...")
    subprocess.run(wget_command, shell=True, check=True)
    print(f"Archivo descargado y guardado en {output_path}")
except subprocess.CalledProcessError as e:
    print(f"Error en la descarga: {e}")
