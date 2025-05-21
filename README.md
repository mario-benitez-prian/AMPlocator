# 🧬 AMPminer

**AMPminer** es una herramienta de línea de comandos para la detección de **péptidos antimicrobianos (AMPs)** en archivos FASTA de secuencias proteicas completas, utilizando un modelo de inteligencia artificial previamente entrenado.

---

## 📦 Instalación

1. Clona el repositorio y navega a la carpeta del proyecto:

```bash
git clone https://github.com/usuario/AMPminer.git
cd AMPminer
```

2. Crea un entorno virtual e instala las dependencias:

```bash
# Usando pip
pip install -r requirements.txt

# O usando conda
conda create -n ampminer python=3.10
conda activate ampminer
pip install -r requirements.txt
```

3. Instala el paquete como herramienta CLI:

```bash
pip install .
```

Esto habilita el comando `ampminer` desde cualquier parte del sistema.

---

## 🧪 Uso

```bash
ampminer <input.fasta> <modelo.keras> <output.fasta> [--max_length MAX_LENGTH]
```

### Argumentos:

```
Posicionales:
  input.fasta           Archivo FASTA de entrada con el proteoma completo.
  modelo.keras          Ruta al modelo entrenado (.keras o .h5).
  output.fasta          Archivo FASTA de salida con las secuencias positivas.

Opcionales:
  --max_length          Longitud máxima para las secuencias (por defecto: 300).
```

### Ejemplo:

```bash
ampminer proteoma.fasta models/amp_locator_model.keras AMPs_predichos.fasta --max_length 300
```

También puede ejecutarse directamente sin instalar como paquete:

```bash
python -m ampminer_app proteoma.fasta models/amp_locator_model.keras salida.fasta --max_length 300
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

## ✅ Pruebas

Puedes crear pruebas unitarias dentro de la carpeta `tests/` para verificar el correcto funcionamiento de los módulos de la app.

---

## 👨‍🔬 Autor

Desarrollado por **Mario Benitez-Prián** como parte de un proyecto de detección de péptidos antimicrobianos con inteligencia artificial.

---
