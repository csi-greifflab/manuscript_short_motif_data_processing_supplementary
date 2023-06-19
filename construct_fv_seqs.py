import pandas as pd
from argparse import ArgumentParser


def construct_fv_seqs(input_file, output_file):
    '''
    gets fv seqs from cdr3s by grafting them to the original tras
    '''
    lseq = 'DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIKR'
    hseq = 'EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSSA'
    h_istart = 98
    h_istop = 108
    df = pd.read_csv(input_file)
    outseqs = []
    for i,row in df.iloc[:,:].iterrows():
        hcdr3 = row.aa_cdr3_seq
        fvseq = hseq[:h_istart] + hcdr3 + hseq[h_istop:]
        outseqs.append([fvseq, lseq])
    outcols = ['hseq', 'lseq']
    outdf = pd.DataFrame(outseqs, columns=outcols)
    outdf['seq_index']= outdf.index.tolist()
    outdf.to_csv(output_file, index=False)


def main():
    parser = ArgumentParser()
    parser.add_argument(
        '-i', '--input', help='csv file containing unique cdr3s and their counts', required=True,
    )
    parser.add_argument('-o', '--output', required=True,)
    args = parser.parse_args()
    construct_fv_seqs(args.input, args.output)


if __name__ == '__main__':
    main()
