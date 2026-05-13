import pandas as pd
import re


def parse_clade(name):
    if re.match(r'^\d+\.pilon\.fasta$', name):
        return 'MG'
    m = re.search(r'_cl([iv]+)_?\.fasta', name)
    if m:
        return m.group(1).upper()
    m = re.search(r'_([iv]{2,})\.fasta', name)
    if m:
        return m.group(1).upper()
    return None  # ex: ca17lbn.fasta


def get_clades(df_all, gt_df):
    vcf_fixed_cols = {'chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'format'}
    sample_cols = [c for c in df_all.columns if c not in vcf_fixed_cols]

    gt_list = list(gt_df.T.columns)
    clades = []
    for iso in gt_list:
        clades.append(parse_clade(iso))
    
    return clades


def parse_year(name):
    m = re.search(r'_((?:19|20)\d{2})(?:_|\.fasta)', name)
    return int(m.group(1)) if m else None


def parse_country(name):
    if re.match(r'^\d+\.pilon\.fasta$', name):
        return None
    body = re.sub(r'\.fasta.*$', '', name)                  # remove .fasta(.ref)?
    body = re.sub(r'_(cl[iv]+|[iv]{2,})_?$', '', body)      # remove clade suffix
    body = re.sub(r'_((?:19|20)\d{2})$', '', body)          # remove year
    # percorre tokens de trás pra frente, coleta os puramente alfabéticos (palavras do país)
    country_tokens = []
    for token in reversed(body.split('_')[2:]):             # skip accession + version
        if token.isalpha():
            country_tokens.insert(0, token)
        else:
            break
    return '_'.join(country_tokens) if country_tokens else None

def get_sample_meta(gt_df):
    sample_meta = pd.DataFrame({
        'clade':   gt_df['clade'],
        'country': [parse_country(s) for s in gt_df.index],
        'year':    [parse_year(s) for s in gt_df.index],
    }, index=gt_df.index)
    return sample_meta
