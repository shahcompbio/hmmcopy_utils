# HMMcopy utilities

## read Counter



```

usage: read_counter.py [-h] [--chromosomes [CHROMOSOMES [CHROMOSOMES ...]]]
                       [-w WINDOW_SIZE] [-m MAPPING_QUALITY_THRESHOLD] [--seg]
                       bam output

positional arguments:
  bam                   specify the path to the input bam file
  output                specify path to the output file

optional arguments:
  -h, --help            show this help message and exit
  --chromosomes [CHROMOSOMES [CHROMOSOMES ...]]
                        specify target chromosomes
  -w WINDOW_SIZE, --window_size WINDOW_SIZE
                        specify bin size
  -m MAPPING_QUALITY_THRESHOLD, --mapping_quality_threshold MAPPING_QUALITY_THRESHOLD
                        threshold for the mapping quality, reads with quality
                        lower than threshold will be ignored
  --seg                 write the output in seg format

```

1. default window_size is 1000
2. seg flag is for debugging purposes only.