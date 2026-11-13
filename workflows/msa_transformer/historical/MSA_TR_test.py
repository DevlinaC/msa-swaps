import esm
import torch
from Bio import SeqIO
import itertools
from typing import List, Tuple
import string
import numpy as np
import pandas as pd
torch.set_grad_enabled(False)

deletekeys = dict.fromkeys(string.ascii_lowercase)
deletekeys["."] = None
deletekeys["*"] = None
translation = str.maketrans(deletekeys)

def remove_insertions(sequence: str) -> str:
    return sequence.translate(translation)

def read_msa(filename: str, nseq: int) -> List[Tuple[str, str]]:
    return [(record.description, remove_insertions(str(record.seq)))
            for record in itertools.islice(SeqIO.parse(filename, "fasta"), nseq)]

def Make_DataFrame(tensor_in, output, query_length):
    mrf = tensor_in
    temp = {'i': [], 'j': [], 'zscore': []}
    for i in range(mrf.shape[0]):
        for j in range(mrf.shape[1]):
            temp['i'].append(i)
            temp['j'].append(j)
            temp['zscore'].append(mrf[i][j].item())
    pd_mtx = pd.DataFrame.from_dict(temp)

    N = int(query_length)
    zero = np.zeros((N, N))
    for index, row in pd_mtx.iterrows():
        zero[int(row['i']), int(row['j'])] = row['zscore']
        zero[int(row['j']), int(row['i'])] = row['zscore']

    top = pd_mtx.loc[pd_mtx['j'] - pd_mtx['i'] > 3].sort_values("zscore", ascending=False)
    return top.head(int(output * mrf.shape[0])), zero

def MSA_TR(msa_dir, files, output, query_length, Beg=0, End=0, name="msa_00"):
    print('Start MSA Transformer...')
    msa_data = []
    for i in files:
        temp = read_msa(f"{msa_dir}/{i}", 64)
        if Beg:
            temp = [(desc, seq[Beg:]) for desc, seq in temp]
        if End:
            temp = [(desc, seq[:-End]) for desc, seq in temp]
        msa_data.append(temp)

    msa_transformer, msa_alphabet = esm.pretrained.esm_msa1b_t12_100M_UR50S()
    msa_transformer = msa_transformer.to(torch.device("cpu"))
    msa_batch_converter = msa_alphabet.get_batch_converter()
    msa_batch_labels, msa_batch_strs, msa_batch_tokens = msa_batch_converter(msa_data)

    print(msa_batch_tokens.size(), msa_batch_tokens.dtype)

    msa_contacts = msa_transformer.predict_contacts(msa_batch_tokens).cpu()

    for n in range(len(files)):
        top, full = Make_DataFrame(msa_contacts[n], output[n], query_length)
        top['i'] += Beg
        top['j'] += Beg
        np.savetxt(f'full_{name}.csv', full, delimiter=",")