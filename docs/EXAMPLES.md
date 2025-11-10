# Example Outputs

## Identity-bin table snippet (from chr_2 run)
```
group	n	median_len	median_id	top_motifs	notes
genome:all	4344	6957.0	0.991	TGCA (3984, 92%), TACA (156, 4%), TGTA (92, 2%), TATA (59, 1%), TGGA (14, 0%)	TSD consensus <40%
chr_2:0.900-0.940	6	5507.5	0.927	TGCA (5, 83%), TACA (1, 17%)	LOW N (n=6); TSD consensus <40%
chr_2:>=0.940	530	7762.0	0.991	TGCA (485, 92%), TACA (18, 3%), TGTA (12, 2%), TATA (10, 2%), TGAC (3, 1%)	TSD consensus <40%
```

## ASCII postcard (chr_2 ≥ 0.94)
```
> AVG chr_2:>=0.940  range:0.940–-
  n: 530
  median len: 7762.0 bp
  IQR: 5074–10784 bp
  identity mean:0.989 median:0.991 range:0.940–-
  LTR5 median: 1047.5 bp   LTR3 median: 1047.5 bp
  Internal median: 5667.0 bp
  strand +:200 (38%)/-:220 (42%)

=============-------------------------------------------------------------------------==============
LTR5:1048 | INT:6714 | LTR3:7762
Q25:5074bp  Q75:10784bp
```

## SVG postcard example

![](example_chr2.svg)

The SVG includes:
- Header + warning
- Metrics table (n, medians, IQR, identity, strand mix)
- LTR/internal bars with median-length labels
- Dashed `Q25` / `Q75` lines indicating the interquartile range
