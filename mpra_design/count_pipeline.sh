#!/bin/bash
. "/nfs/apps/test/conda_test/etc/profile.d/conda.sh"
conda activate /lab/solexa_page/hannah/mpraflow/MPRAflow/conda-env/MPRAflow
cd  /lab/solexa_page/hannah/mpraflow/MPRAflow
nextflow run count.nf -w /lab/solexa_page/hannah/220516_mpra/Count_Basic/work --experiment-file "/lab/solexa_page/hannah/220516_mpra/Count_Basic/data/experiment.csv" --dir "/lab/solexa_page/hannah/220516_mpra/FASTQ" --outdir "/lab/solexa_page/hannah/220516_mpra/Count_Basic/output" --design "/lab/solexa_page/hannah/lib_mpra/lib_mpra_ordered2.fa" --association "/lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021/concat_results/TCGTCAA/TCGTCAA_filtered_coords_to_barcodes.pickle"
