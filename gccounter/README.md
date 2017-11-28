# HMMcopy utilities

## GC Counter



```

usage: gc_counter.py [-h] [--chromosomes [CHROMOSOMES [CHROMOSOMES ...]]]
                     [--window_size WINDOW_SIZE]
                     reference output

positional arguments:
  reference             path to the reference fasta file. The reference must
                        be indexed with bowtie-build and samtools index
  output                path to the output file

optional arguments:
  -h, --help            show this help message and exit
  --chromosomes [CHROMOSOMES [CHROMOSOMES ...]]
                        specify target chromosomes
  --window_size WINDOW_SIZE
                        specify window size.


```


# chromosomes is optional
# default window_size is 1000