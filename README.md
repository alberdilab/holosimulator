# holosimulator
Code repository of the HoloSimulator software


## Examples

### Simulation of two populations of lizards

Prepare the working directory
```sh
mkdir holosimulator_lizards
cd holosimulator_lizards
```

#### Fetch and subset host genome

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

#### Simulate genomic divergence

Using 'holosimulate mutations' simulate two reference genomes that are ca. 1% (99% ANI) divergent from the original. 
Each of these genomes will represent the centroid of the populations A and B, respectivelly.
```sh
holosimulate mutations -i host.fna -o host1.fna -a 0.99
holosimulate mutations -i host.fna -o host5.fna -a 0.99
```

Using 'holosimulate mutations' simulate four more reference genomes that are ca. 0.1% (99.9% ANI) divergent from the reference population genomes. 
```sh
holosimulate mutations -i host1.fna -o host2.fna -a 0.999
holosimulate mutations -i host1.fna -o host3.fna -a 0.999
holosimulate mutations -i host1.fna -o host4.fna -a 0.999
holosimulate mutations -i host1.fna -o host5.fna -a 0.999
holosimulate mutations -i host5.fna -o host6.fna -a 0.999
holosimulate mutations -i host5.fna -o host7.fna -a 0.999
holosimulate mutations -i host5.fna -o host8.fna -a 0.999
holosimulate mutations -i host5.fna -o host9.fna -a 0.999
```

