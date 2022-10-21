# exoclasma-note
Annotation of genomic variants

## Dependencies

Perl:

```bash
sudo apt install perl
```

Annovar:

```bash
GENOME="hg19" # hg38 if you prefer
curl --output annovar.latest.tar.gz www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz
tar -xzf annovar.latest.tar.gz
cd annovar
perl annotate_variation.pl -downdb -buildver ${GENOME} -webfrom annovar refGene humandb/
perl annotate_variation.pl -downdb -buildver ${GENOME} -webfrom annovar knownGene humandb/
perl annotate_variation.pl -downdb -buildver ${GENOME} -webfrom annovar ensGene humandb/
cd ..
rm annovar.latest.tar.gz
```
