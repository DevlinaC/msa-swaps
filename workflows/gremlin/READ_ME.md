## READ ME ##

# GREMLIN Biowulf environment-variable workflow

This package runs the uploaded GREMLIN code once per .fas file under /data/AF_interpretation/msa_swaps/original_proper_msas and under on /data/AF_interpretation/msa_swaps/original_junk_msas.

The Python runner is gremlin_biowulf_env.py. It reads configuration from environment variables & doesn't call a command-line parser (local environment I tested on had issues w/ this). The shell scripts set those environment variables for each task.

# Code
-gremlin_biowulf_env.py: runs one .fas alignment & writes CSV/PNG/TXT outputs.
-gremlin_original_proper_msas_job.sh: starts a batch job to run gremlin_biowulf_env.py for each .fas file in /data/AF_interpretation/msa_swaps/original_proper_msas.
-gremlin_original_junk_msas_job.sh: starts a batch job to run gremlin_biowulf_env.py for each .fas file in /data/AF_interpretation/msa_swaps/original_junk_msas.

# Files
-full*.csv: full GREMLIN-scored residue pair table; APC scores included from correction; zscore normalized APC scores included; amino-acid labels included; first 2 numbers are the pandas row index from the OG full dataframe; i = the 1st residue position in the pair, j = the 2nd
-top*.csv: top GREMLIN-scored residue pair table; residue pairs w/ the top 10 direct coupling signals; candidate residue;residue contacts/coevolving residue pairs; APC scores included from correction; zscore normalized APC scores included; amino-acid labels included; pairs filtered to where j-i > 5 & then sorted by decreasing z-score; first 2 numbers are the pandas row index from the OG full dataframe; i = the 1st residue position in the pair, j = the 2nd
-*_heatmap.png: visual version of GREMLIN pair scores from the top-contact table; squares = 1 residue-position pair; darker blue = pair w/ a stronger GREMLIN coupling score
-*_summary.txt: list of environment variables & values from the run; includes associated files (& their path from the run)