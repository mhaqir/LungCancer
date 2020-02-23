# Lung Cancer


For Extracting required information for running CTPsingle from a vcf file, the following line of code can be used.

```python
python input.py <normal_id+tumor_id.vcf> <Male, Female>
```

The output will be a file with name tumor_id.frq with the format specified in https://github.com/nlgndnmz/CTPsingle for running CTPsingle tool.

The coverage distribution for the sample SC284808 is as shown in figure below.
![alt text](plots/combined_SC284808_cov.png?raw=true "Title")
