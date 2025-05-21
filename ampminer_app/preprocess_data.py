import numpy as np
import tensorflow as tf

def preprocess_data(data, max_length):
    sequences = []
    amino_acids = 'ACDEFGHIKLMNPQRSTVWYX'
    aa_to_index = {aa: i for i, aa in enumerate(amino_acids)}

    for _, row in data.iterrows():
        raw_seq = row['Full_Seq'].strip().upper()
        cleaned_seq = ''.join([aa if aa in amino_acids else 'X' for aa in raw_seq])
        indices = [aa_to_index.get(aa, aa_to_index['X']) for aa in cleaned_seq]
        one_hot_seq = tf.keras.utils.to_categorical(indices, num_classes=21)
        one_hot_seq = one_hot_seq[:max_length]
        one_hot_seq = np.pad(one_hot_seq, ((0, max_length - len(one_hot_seq)), (0, 0)), mode='constant')
        sequences.append(one_hot_seq)

    return np.array(sequences)
