import re
import pandas as pd


def to_snake(name: str) -> str:
    """
    Convert camelCase/PascalCase column names to snake_case,
    preserving dots introduced by pandas.json_normalize.
    """
    parts = name.split(".")
    snake_parts = []

    for part in parts:
        part = re.sub(r"(.)([A-Z][a-z]+)", r"\1_\2", part)
        part = re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", part)
        snake_parts.append(part.lower())

    return ".".join(snake_parts)


def normalize_metadata_cols(df: pd.DataFrame) -> pd.DataFrame:
    """
    Return a copy of the DataFrame with normalized snake_case column names.
    """
    df_norm = df.copy()
    df_norm.columns = [to_snake(str(col)) for col in df_norm.columns]

    return df_norm