import tensorflow as tf
from ampminer_app.preprocess_data import preprocess_data
from ampminer_app.fasta_parser import load_fasta_as_dataframe, write_fasta

def run_prediction(fasta_file, model_file, output_file, max_length):
    print("[INFO] Cargando modelo...")
    model = tf.keras.models.load_model(model_file)

    print("[INFO] Cargando datos del proteoma...")
    df = load_fasta_as_dataframe(fasta_file)

    print("[INFO] Preprocesando secuencias...")
    X = preprocess_data(df, max_length)

    print("[INFO] Realizando predicciones...")
    predictions = model.predict(X)
    predicted_labels = (predictions > 0.5).astype(int).flatten()

    print("[INFO] Filtrando secuencias predichas como positivas...")
    df["Predicted"] = predicted_labels
    df_positive = df[df["Predicted"] == 1]

    print(f"[INFO] Guardando {len(df_positive)} secuencias positivas en: {output_file}")
    write_fasta(df_positive, output_file)

    print("[INFO] Proceso completado.")
