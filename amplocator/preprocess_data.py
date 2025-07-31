import numpy as np

def preprocess_fasta_sequences(sequences):

    import tensorflow as tf

    np.random.seed(42)
    amino_acids = 'ACDEFGHIKLMNPQRSTVWYX'
    aa_to_index = {aa: i for i, aa in enumerate(amino_acids)}
    num_classes = len(amino_acids)

    max_length = 300

    #print(f"[DEBUG] Se cargaron en preprocess {len(sequences)} secuencias desde el FASTA.")

    processed = []

    for i, seq in enumerate(sequences):
        #print(f"Secuencia {i+1}/{len(sequences)} original: {seq}")

        # Clean seq 
        seq = seq.strip().upper()
        cleaned = ''.join([aa if aa in amino_acids else 'X' for aa in seq])
        #print(f"  Secuencia limpiada: {cleaned}")

        # Cut sequences to max_length
        cleaned = cleaned[:max_length]
        seq_len = len(cleaned)
        pad_total = max_length - seq_len
        #print(f"  Longitud tras truncamiento: {seq_len}")

        # Random padding at both ends
        pad_left = np.random.randint(0, pad_total + 1)
        pad_right = pad_total - pad_left
        #print(f"  Padding ? izquierda: {pad_left}, derecha: {pad_right}")

        # One-hot encoding 
        indices = [aa_to_index[aa] for aa in cleaned]
        one_hot = tf.keras.utils.to_categorical(indices, num_classes=num_classes)

        # Padding with zero vectors
        one_hot = np.pad(one_hot, ((pad_left, pad_right), (0, 0)), mode='constant')

        #print(f"  Shape final: {one_hot.shape}\n")
        processed.append(one_hot)

    #print("Preprocesamiento completo.\n")
    return np.array(processed)

"""secuencias = ["MKTIIALSYIFCLVFAD", "GVLKKLGX*QY"]
datos = preprocess_fasta_sequences(secuencias, max_length=30)
print(datos[-1])"""
