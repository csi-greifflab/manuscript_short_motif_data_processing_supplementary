import pandas as pd
from pathlib import Path
from argparse import ArgumentParser


def main():
    parser = ArgumentParser()
    parser.add_argument(
        '-i',
        '--input_dir',
        help='directory with csv files containing unique cdr3s from different samples',
        required=True,
    )
    parser.add_argument(
        '-o',
        '--output_dir',
        required=True,
    )
    args = parser.parse_args()
    
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    files = list(input_dir.iterdir())
    files = [f for f in files if f.suffix == '.csv']
    dfs = [pd.read_csv(f) for f in files]
    sets = [set(df['aa_cdr3_seq']) for df in dfs]
    seen_cdr3s = set()
    for f, df in zip(files, dfs):
        df['aa_cdr3_seq'] = df['aa_cdr3_seq'].map(lambda s: s[2:-1])
        df = df.loc[
            ~df['aa_cdr3_seq'].isin(seen_cdr3s) & \
            (df['count'] > 1) \
        ].sort_values(by='count', ascending=False)
        df.to_csv()
        df = df.loc[~df['aa_cdr3_seq'].isin(seen_cdr3s)]
        df.to_csv(output_dir / f.name, index=False)
        seen_cdr3s.update(df['aa_cdr3_seq'])


if __name__ == '__main__':
    main()