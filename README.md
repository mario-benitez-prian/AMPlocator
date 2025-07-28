# 🧬 AMPlocator

**AMPlocator** is a command-line tool that employs two deep-learning models to predict **antimicrobial peptide** precursors and mature/active peptide location.

---

## 📦 Instalation

1. Clone the repository and navigate to the project folder:

```bash
git clone https://github.com/usuario/AMP_locator.git
cd AMPlocator
```

2. Create a virtual environment and install dependencies:

```bash
# Without creating a virtual environment
pip install -r requirements.txt

# Creating a virtual environment (recommended)
conda create -n amplocator python=3.10
conda activate amplocator
pip install -r requirements.txt
```

3. Install the package as a CLI tool (you can run the command 'amplocator' wherever you are in your system):

```bash
pip install .
```

---

## 🧪 How to use?

```bash
amplocator --input_file <input.fasta> --output_file <output file> --mode [precursor, full, locator]
```

### Arguments:

```
  --input_file          File in fasta format with protein data (amino acids)
  --output_file         Output file prefix or path to save the results (several files will be saved in .fasta .tsv and .excel format)
  --mode                [precursor, full, locator]
```

### Examples of usage:

#### Predict only AMP precursors

```bash
amplocator proteome.fasta precursors --mode precursor
```

#### Predict only AMP mature peptide regions in precursor candidates

```bash
amplocator precursor_candidates.fasta amp_regions --mode locator
```

#### Predict both precursors and AMP mature peptide regions

```bash
amplocator proteome.fasta precursors_and_amp_regions --mode full
```

---

## 📄 Requirements

The necessary dependencies are listed below:

- `tensorflow >= 2.x`
- `numpy`
- `pandas`
- `biopython`
- `argparse`

You can install them with:

```bash
pip install -r requirements.txt
```

---

## 👨‍🔬 Autor

Developed by Mario Benítez-Prián et al.

---
