# modelestimator --- Infer sequence evolution rate matrices from a MSA


## Example usage

``` shell
    modelestimator fasta -t 0.001 file1.fa file2.fa file3.fa
```
Infer a rate matrix (written to stdout) from three alignment files in Fasta format.

``` shell
    modelestimator fasta -b 200 file.fa
```
Try the experimental bootstrapping feature (200 replicates) on a Fasta multialignment.

## Syntax

```
modelestimator <format> <options> infiles
```


`<format>` should be either FASTA, STOCKHOLM or PHYLIP format.

Output is a rate matrix and residue distribution vector.

### Options

```
    -threshold or -t <f>
	Stop when consecutive iterations do not change by more than <f>. Default is 0.001.

    -bootstrap or -b <r>
	Perform bootstrapping on multialignment with <r> resamplings. Only one infile should be given in this mode. Returns bootstrap norm.
```
