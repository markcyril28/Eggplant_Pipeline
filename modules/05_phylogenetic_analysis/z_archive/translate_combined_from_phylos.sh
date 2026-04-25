#!/bin/bash

>Combined_Phylo_GRF_ProtSeq_complete_named.fasta

seqkit translate -f 1 \
    Combined_Phylo_GRF_NucSeq_complete_named.fasta > Combined_Phylo_GRF_ProtSeq_complete_named.fasta

seqkit translate -f 1 \
    Combined_Phylo_GIF_NucSeq_complete_named.fasta > Combined_Phylo_GIF_ProtSeq_complete_named.fasta
