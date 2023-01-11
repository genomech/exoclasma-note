# exoclasma-note

Annotation of genomic variants.

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

## Install

```bash
sudo python3 -m pip install exoclasma-note
```

## Usage

```bash
exoclasma-note -u ${unit_json} -a ${annovar_folder} -g ${genome}
```

* `unit_json`: Unit JSON file which was created by exoclasma-pipe
* `annovar_folder`: Path to ANNOVAR folder where perl scripts are located
* `genome`: Genome assembly which use ANNOVAR (i.e., hg19, hg38, etc.)
