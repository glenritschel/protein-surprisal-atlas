# Smoke Test Results

- **Proteins Scored**: 10
- **Wall-clock Time**: 37 seconds
- **Proteins/Hour**: 972.97
- **Peak RAM Usage**: 0.13 GB (Note: this is the memory of the reporting script, not the peak of the scoring run itself, though CPU usage was minimal)
- **GPU Memory**: N/A (Run on CPU)

## Protein Table Head

```
  uniprot_id                       gene_symbol  total_surprisal_bits  mean_residue_surprise
0     P00533              EGFR ERBB ERBB1 HER1           3837.380806               3.171389
1     P01308                               INS            409.915901               3.726508
2     P02452                            COL1A1           3499.089060               2.390088
3     P02649                              APOE           1034.126968               3.262230
4     P04406  GAPDH GAPD CDABP0047 OK/SW-cl.12           1290.694076               3.852818
```

## Residue Table Head

```
  uniprot_id           gene_symbol  position reference_residue  surprisal_bits  probability  window_id  distance_from_window_edge
0     P00533  EGFR ERBB ERBB1 HER1         1                 M        0.013730     0.990528          0                          0
1     P00533  EGFR ERBB ERBB1 HER1         2                 R        3.138238     0.113579          0                          1
2     P00533  EGFR ERBB ERBB1 HER1         3                 P        3.226922     0.106807          0                          2
3     P00533  EGFR ERBB ERBB1 HER1         4                 S        4.117838     0.057598          0                          3
4     P00533  EGFR ERBB ERBB1 HER1         5                 G        3.580451     0.083594          0                          4
```
