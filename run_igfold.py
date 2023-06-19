import os
import sys
from argparse import ArgumentParser

import pandas as pd

from igfold import IgFoldRunner
from igfold.refine.pyrosetta_ref import init_pyrosetta
init_pyrosetta()
from mpi4py import MPI


def run_igfold_mpi(input_file, output_dir, ncpu):
    '''
    runs IgFold
    '''
    df = pd.read_csv(input_file)
    nbins = int(ncpu)-1
    nsamples = int(len(df)/nbins)
    frames = [ df.iloc[i*nsamples:(i+1)*nsamples].copy() for i in range(nbins+1) ]
    outtag = input_file.split('/')[-1].split('.')[0]
    # add parallel computing with mpi (scatter)
    # MPI stuff
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    if rank == 0:
       data = frames
    else:
        data=None
    data = comm.scatter(data, root=0)
    igfold = IgFoldRunner()
    counter = 0
    for i, row in data.iloc[:,:].iterrows():
        ndata = data.shape[0]
        print('folding %s/%s' % (counter, ndata)) 
        counter += 1
        hseq = row['hseq']
        lseq = row['lseq']
        pred_pdb = os.path.join(output_dir, outtag + '_' + str(row['seq_index']) + '_igfold.pdb')
        sequences={
            'H': hseq,
            'L': lseq,
        }
        with open(pred_pdb, 'w') as otp:
            pass
        # igfold.fold(
        #     pred_pdb, # Output PDB file
        #     sequences=sequences, # Antibody sequences
        #     do_refine=True, # Refine the antibody structure with PyRosetta
        #     do_renum=True, # Send predicted structure to AbNum server for Chothia renumbering
        # ) 


def main():
    parser = ArgumentParser()
    parser.add_argument(
        '-i', '--input', required=True,
    )
    parser.add_argument('-o', '--output_dir', required=True)
    parser.add_argument('-n', '--n_cpu', required=True)
    args = parser.parse_args()
    run_igfold_mpi(args.input, args.output_dir, args.n_cpu)


if __name__ == '__main__':
    main()
