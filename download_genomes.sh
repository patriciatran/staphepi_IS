#!/bin/bash

datasets download genome taxon 'Staphylococcus epidermidis' --assembly-level complete --include genome,cds,protein,gbff --mag exclude --assembly-source 'RefSeq'

ls -lhs
