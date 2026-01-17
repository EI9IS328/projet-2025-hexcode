#!/usr/bin/env python3

import os
import sys
import pandas as pd

# -----------------------
# Fonctions utilitaires
# -----------------------

def parse_value(token):
    """Convertit un token en int ou float si possible, sinon retourne la chaîne"""
    try:
        if "." in token or "e" in token.lower():
            return float(token)
        return int(token)
    except ValueError:
        return token

def expand_token(token):
    """Décompresse un token de type Nx(...) ou NxVAL ou simple VAL."""
    token = token.strip()

    # Nx(...)   (valeurs séparées par des virgules)
    if "x(" in token and token.endswith(")"):
        n_str, values_str = token.split("x(", 1)
        try:
            n = int(n_str)
        except ValueError:
            n = 1
        values_str = values_str[:-1]
        values = [parse_value(v.strip()) for v in values_str.split(",")]
        return values * n
    if "(" in token and token.endswith(")"):
        return [token[1:-1]]

    # NxVAL
    if "x" in token:
        n_str, value_str = token.split("x", 1)
        try:
            n = int(n_str)
        except ValueError:
            n = 1
        value = parse_value(value_str)
        return [value] * n

    # simple VAL
    return [parse_value(token)]

# -----------------------
# Décompression .rle
# -----------------------

def decompress_rle_text(text):
    """Décompression pour les fichiers .rle classiques"""
    result = {}
    for line in text.strip().splitlines():
        parts = line.split()
        if not parts:
            continue
        name = parts[0]
        data = []
        for token in parts[1:]:
            data.extend(expand_token(token))
        result[name] = data
    return result

def write_rle_csv(data, output_file):
    indices = data["index"]
    timesteps = data["timestep"]
    xs = data["x"]
    ys = data["y"]
    zs = data["z"]
    pressures = data["pressure"]

    num_points = len(indices)
    num_timesteps = len(timesteps) // num_points

    rows = []
    for t in range(num_timesteps):
        for i in range(num_points):
            rows.append({
                "index": indices[i],
                "timestep": timesteps[t * num_points + i],
                "x": xs[i],
                "y": ys[i],
                "z": zs[i],
                "pressure": pressures[t * num_points + i]
            })
    df = pd.DataFrame(rows)
    df.to_csv(output_file, sep=" ", index=False, float_format="%.6f")

# -----------------------
# Décompression .quanti.rle
# -----------------------

def decompress_quanti_rle_text(text):
    """
    Décompression d'un fichier .quanti.rle.
    Retourne :
      - first_line : la première ligne (metadata)
      - df : DataFrame avec toutes les colonnes
    """
    lines = text.strip().splitlines()
    first_line = lines[0]  # on garde la première ligne
    data = {}

    for line in lines[1:]:
        if not line.strip():
            continue
        parts = line.split(maxsplit=1)
        if len(parts) < 2:
            continue
        col_name, tokens_str = parts
        tokens = tokens_str.split()

        # décompression NxVAL / Nx(...)
        values = []
        for token in tokens:
            values.extend(expand_token(token))

        data[col_name] = values

    # Déterminer la longueur maximale
    max_len = max(len(v) for v in data.values())

    # Remplir les colonnes plus courtes avec NaN
    for k in data:
        if len(data[k]) < max_len:
            data[k].extend([float('nan')] * (max_len - len(data[k])))

    df = pd.DataFrame(data)
    return first_line, df


def write_quanti_rle_csv(first_line, df, output_file):
    """
    Écriture du DataFrame en .quanti avec la première ligne metadata.
    On réorganise dynamiquement selon les colonnes présentes.
    """
    # Détecter un ordre standard
    if all(c in df.columns for c in ["snap", "timestep", "x", "y", "z", "pressure"]):
        desired_order = ["snap", "timestep", "x", "y", "z", "pressure"]
    elif all(c in df.columns for c in ["plane", "timestep", "i", "j", "pressure"]):
        desired_order = ["plane", "timestep", "i", "j", "pressure"]
    else:
        desired_order = list(df.columns)  # sinon on garde l'ordre trouvé

    df = df[[c for c in desired_order if c in df.columns]]  # filtrer les colonnes existantes

    with open(output_file, "w") as f:
        f.write(first_line + "\n")
        df.to_csv(f, sep=" ", index=False, float_format="%.6f")

# -----------------------
# Décompression .quanti
# -----------------------

def decompress_quanti_file(file_path):
    """Décompression pour les fichiers .quanti"""
    with open(file_path, "r") as f:
        pmin, pmax, dcompress = map(float, f.readline().split())

    df = pd.read_csv(file_path, skiprows=1, sep=r"\s+", engine="python")
    
    if "pressure" in df.columns:
        df["pressure"] = df["pressure"] * dcompress + pmin

    return df

def write_quanti_csv(df, output_file):
    df.to_csv(output_file, sep=" ", index=False, float_format="%.6f")

# -----------------------
# Parcours des fichiers
# -----------------------

if len(sys.argv) != 2:
    print("Usage: python decompress_pandas.py <directory>")
    sys.exit(1)

root_dir = sys.argv[1]

for dirpath, dirnames, filenames in os.walk(root_dir):
    for filename in filenames:
        file_path = os.path.join(dirpath, filename)
        try:
            if filename.endswith(".quanti.rle"):
                print(f"Traitement .quanti.rle : {file_path}")
                with open(file_path, "r") as f:
                    text = f.read()
                first_values, data = decompress_quanti_rle_text(text)
                temp_file = file_path.replace(".quanti.rle", ".quanti")
                write_quanti_rle_csv(first_values, data, temp_file)
                df = decompress_quanti_file(temp_file)
                output_file = file_path.replace(".quanti.rle", ".csv")
                write_quanti_csv(df, output_file)

            elif filename.endswith(".rle") and not filename.endswith(".quanti.rle"):
                print(f"Traitement .rle : {file_path}")
                with open(file_path, "r") as f:
                    text = f.read()
                data = decompress_rle_text(text)
                output_file = file_path.replace(".rle", ".csv")
                write_rle_csv(data, output_file)

            elif filename.endswith(".quanti"):
                print(f"Traitement .quanti : {file_path}")
                df = decompress_quanti_file(file_path)
                output_file = file_path.replace(".quanti", ".csv")
                write_quanti_csv(df, output_file)

            print(f"Décompression terminée : {file_path}")

        except Exception as e:
            print(f"Erreur lors du traitement de {file_path} : {e}")
