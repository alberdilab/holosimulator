# holosimulator
Code repository of the HoloSimulator software


## Examples

### Simulation of two populations of lizards

Prepare the working directory
```sh
mkdir holosimulator_lizards
cd holosimulator_lizards
```

#### Prepare host genome

Download the genome and check contig headers.
```sh
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/329/235/GCF_004329235.1_PodMur_1.0/GCF_004329235.1_PodMur_1.0_genomic.fna.gz
gunzip GCF_004329235.1_PodMur_1.0_genomic.fna.gz
grep ">" GCF_004329235.1_PodMur_1.0_genomic.fna | head
```

Subset the genome to the first two chromosomes to avoid unnecessary data volumes. 
```sh
awk '/^>NC_041312.1/{flag=1;print;next} /^>/{flag=0} flag' GCF_004329235.1_PodMur_1.0_genomic.fna > host.fna
awk '/^>NC_041313.1/{flag=1;print;next} /^>/{flag=0} flag' GCF_004329235.1_PodMur_1.0_genomic.fna >> host.fna
```
