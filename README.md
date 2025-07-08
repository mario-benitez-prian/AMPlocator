# 🧬 AMPlocator

**AMPminer** es una herramienta de línea de comandos para la detección de **péptidos antimicrobianos (AMPs)** en archivos FASTA de secuencias proteicas completas, utilizando un modelo de inteligencia artificial previamente entrenado.

---

## 📦 Instalación

1. Clona el repositorio y navega a la carpeta del proyecto:

```bash
git clone https://github.com/usuario/AMP_locator.git
cd AMPlocator
```

2. Crea un entorno virtual e instala las dependencias:

```bash
# Usando pip
pip install -r requirements.txt

# O usando conda
conda create -n amplocator python=3.10
conda activate amplocator
pip install -r requirements.txt
```

3. Instala el paquete como herramienta CLI:

```bash
pip install .
```

Esto habilita el comando `amplocator` desde cualquier parte del sistema.

---

## 🧪 Uso

```bash
amplocator <input.fasta> <output.fasta> --mode [precursor, full, locator]
```

### Argumentos:

```
Posicionales:
  input.fasta           Archivo FASTA de entrada con el proteoma completo.
  output.fasta          Archivo FASTA de salida con las secuencias positivas.

Opcionales:
  --mode                [precursor, full, locator]
```

### Ejemplo:

```bash
amplocator proteoma.fasta AMPs_predichos.fasta --mode precursor
```

---

## 📄 Requisitos

Las dependencias necesarias están listadas en `requirements.txt`. Las principales son:

- `tensorflow >= 2.x`
- `numpy`
- `pandas`
- `biopython`
- `argparse`

Instálalas con:

```bash
pip install -r requirements.txt
```

---

## 👨‍🔬 Autor

Desarrollado por **Mario Benitez-Prián** como parte de un proyecto de detección de péptidos antimicrobianos con inteligencia artificial.

---
