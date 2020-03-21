# Lung Cancer


For Extracting required information for running CTPsingle from a vcf file, the following line of code can be used.

```python
python input.py <normal_id+tumor_id.vcf> <Male, Female>
```

The output will be a file with name tumor_id.frq with the format specified in https://github.com/nlgndnmz/CTPsingle for running CTPsingle tool.

The coverage distribution for the sample SC284808 is as shown in figure below.
![alt text](plots/combined_SC284808_cov.png?raw=true "Title")


The cellular prevalence plot for the sample SC284808 is as shown in the figure below.
The x-axis represents 101 intervals that a VAF value may fall into (0 <= cp <= 0.01, 0.01 < cp <= 0.02, ..., 0.99 < cp <= 1, cp > 1), and the y-axis is the counts for each interval.
![alt text](plots/VAF_count.png?raw=true "Title")

CONETT input graph for 69 samples (out of 82) is as follows:
![alt text](plots/CONETT_graph.png?raw=true "Title")
