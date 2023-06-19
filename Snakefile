import os 


DATA_DIR = 'data'
WORK_DIR = 'intermediate_results'
UNIQUE_CDR3_DIR = os.path.join(WORK_DIR, 'unique_cdr3s')
REFINED_UNIQUE_CDR3_DIR = os.path.join(WORK_DIR, 'unique_cdr3s_refined')
IGFOLD_DIR = os.path.join(WORK_DIR, 'igfold')
N_CPU_IGFOLD = 4

samples = {
    'high': ('1_small.R1.fastq.gz', '1_small.R2.fastq.gz'),
}


def get_reads(sample, i):
    return os.path.join(DATA_DIR, samples[sample][i])


rule all:
    input:
        'count_unique_cdr3s.out',
        # expand(os.path.join(REFINED_UNIQUE_CDR3_DIR, '{sample}_unique_cdr3s_fv.csv'), sample=samples),
        expand(os.path.join(IGFOLD_DIR, '{sample}'), sample=samples),


rule compile_cpp:
    input:
        'count_unique_cdr3s.cpp'
    output:
        'count_unique_cdr3s.out'
    shell:
        'g++ {input} -o {output}'


rule merge_reads:
    input:
        r1=lambda wildcards: get_reads(wildcards.sample, 0),
        r2=lambda wildcards: get_reads(wildcards.sample, 1),
    output:
        os.path.join(WORK_DIR, '{sample}_merged.assembled.fastq')
    shell:
        f'mkdir -p {UNIQUE_CDR3_DIR} && ' 
        'pear -f {input.r1} -r {input.r2} -o ' + os.path.join(WORK_DIR, '{wildcards.sample}_merged')


rule extract_unique_cdr3s:
    input:
        p='count_unique_cdr3s.out',
        r=os.path.join(WORK_DIR, '{sample}_merged.assembled.fastq'),
    output:
        os.path.join(UNIQUE_CDR3_DIR, '{sample}_unique_cdr3s.csv'),
    shell:
        f'mkdir -p {UNIQUE_CDR3_DIR} && '
        './count_unique_cdr3s.out {input.r} {output}'


rule refine_unique_cdr3s:
    input:
        os.path.join(UNIQUE_CDR3_DIR, '{sample}_unique_cdr3s.csv'),
    output:
        os.path.join(REFINED_UNIQUE_CDR3_DIR, '{sample}_unique_cdr3s.csv'),
    shell:
        f'mkdir -p {REFINED_UNIQUE_CDR3_DIR} && '
        f'python refine_cdr3s.py -i {UNIQUE_CDR3_DIR} -o {REFINED_UNIQUE_CDR3_DIR}'


rule construct_fv_sequences:
    input:
        os.path.join(REFINED_UNIQUE_CDR3_DIR, '{sample}_unique_cdr3s.csv'),
    output:
        os.path.join(REFINED_UNIQUE_CDR3_DIR, '{sample}_unique_cdr3s_fv.csv'),
    shell:
        'python construct_fv_seqs.py -i {input} -o {output}'


rule run_igfold:
    input:
        os.path.join(REFINED_UNIQUE_CDR3_DIR, '{sample}_unique_cdr3s_fv.csv'),
    output:
        directory(os.path.join(IGFOLD_DIR, '{sample}')),
    shell:
        'mkdir -p {output}' + ' && '
        f'mpiexec -n {N_CPU_IGFOLD} python run_igfold.py -n {N_CPU_IGFOLD} -i ' +'{input} -o {output}'
