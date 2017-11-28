# HMMcopy utilities

## Mappability Counter


Please use run.py for

1. simulate reads from fasta (fasta_to_read.py)
2. align simulated reads (bowtie_index.py)
3. create bigwig file with information from aligner (read_to_map.py)
4. convert bigwig to mappability wig file (map_counter.py)



```

usage: run.py [-h] [--chromosomes [CHROMOSOMES [CHROMOSOMES ...]]]
              [--window_size WINDOW_SIZE]
              [--mapcounter_window_size MAPCOUNTER_WINDOW_SIZE]
              [--aligner {bowtie}] [--maxhits MAXHITS] [--cleanup]
              reference temp_dir output

positional arguments:
  reference             path to the reference fasta file. The reference must
                        be indexed with bowtie-build and samtools index
  temp_dir              path to a temporary dir for storing intermediate
                        files. the directory will be removed after the run is
                        complete if the clean up flag is set.
  output                path to the output mappability wig file

optional arguments:
  -h, --help            show this help message and exit
  --chromosomes [CHROMOSOMES [CHROMOSOMES ...]]
                        specify target chromosomes
  --window_size WINDOW_SIZE
                        specify window size for simulating reads
  --mapcounter_window_size MAPCOUNTER_WINDOW_SIZE
                        specify the window size for mapcounter wig file
                        generation
  --aligner {bowtie}    specify the aligner to use for aligning the simulated
                        reads
  --maxhits MAXHITS     threshold for the num_hits value from bowtie
  --cleanup             if set, the tempdir will be deleted at the end of
                        execution

```


# chromosomes is optional
# default window_size is 50
# only bowtie is supported for --aligner option
# default maxhits is 4