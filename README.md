# holosimulator
Code repository of the HoloSimulator software

## Installation

The easiest way to install **Holosimulator** is to create a conda environment containing the software and all its dependencies.

```sh
conda env create -f https://raw.githubusercontent.com/alberdilab/holosimulator/refs/heads/main/environment.yaml
```

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

Using 'holosimulate mutations' simulate two population-reference genomes that are ca. 1% (99% ANI) divergent from the original. 
Each of these genomes will represent the centroid of the populations A and B, respectivelly.
```sh
holosimulator mutations -i host.fna -o host1.fna -a 0.99
holosimulator mutations -i host.fna -o host6.fna -a 0.99
```

Using 'holosimulate mutations' simulate four more reference genomes that are ca. 0.1% (99.9% ANI) divergent from the reference population genomes. 
```sh
holosimulator mutations -i host1.fna -o host2.fna -a 0.999
holosimulator mutations -i host1.fna -o host3.fna -a 0.999
holosimulator mutations -i host1.fna -o host4.fna -a 0.999
holosimulator mutations -i host1.fna -o host5.fna -a 0.999
holosimulator mutations -i host6.fna -o host7.fna -a 0.999
holosimulator mutations -i host6.fna -o host8.fna -a 0.999
holosimulator mutations -i host6.fna -o host9.fna -a 0.999
holosimulator mutations -i host6.fna -o host10.fna -a 0.999
```



Bacterial taxa

Fusobacterium_A ulcerans_A	
Buttiauxella ferragutiae	
Hafnia alvei
Bacteroides thetaiotaomicron		
Bacteroides fragilis	
Bacteroides sp902362375	
Spiroplasma phoeniceum	
Stenotrophomonas maltophilia_AA	
Clostridioides sp020641055	
Parabacteroides sp900540715	
Wolbachia sp001439985	
Buttiauxella davisae_A	
Morganella morganii_A	
Escherichia marmotae	
Robinsoniella peoriensis	
Lactococcus lactis	
Rahnella variigena	
Klebsiella variicola	
Phascolarctobacterium faecium	
UMGS946 sp900544965	

#### Validate genomic divergence

Using 'fastANI', confirm the simulated genomes exhibit the expected level of dissimilarity.


### Simulation of five species of corvids

Prepare the working directory
```sh
mkdir holosimulator_corvids
cd holosimulator_corvids
```

#### Fetch and subset host genome

Download the genomes and check contig headers.
```sh
# Corvus hawaiiensis
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/020/740/725/GCF_020740725.1_bCorHaw1.pri.cur/GCF_020740725.1_bCorHaw1.pri.cur_genomic.fna.gz
gunzip GCF_020740725.1_bCorHaw1.pri.cur_genomic.fna.gz
grep ">" GCF_020740725.1_bCorHaw1.pri.cur_genomic.fna | head

# Corvus moneduloides
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/650/955/GCF_009650955.1_bCorMon1.pri/GCF_009650955.1_bCorMon1.pri_genomic.fna.gz
gunzip GCF_009650955.1_bCorMon1.pri_genomic.fna.gz
grep ">" GCF_009650955.1_bCorMon1.pri_genomic.fna | head

# Corvus cornix
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/738/735/GCF_000738735.6_ASM73873v6/GCF_000738735.6_ASM73873v6_genomic.fna.gz
gunzip GCF_000738735.6_ASM73873v6_genomic.fna.gz
grep ">" GCF_000738735.6_ASM73873v6_genomic.fna | head

# Aphelocoma coerulescens
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/041/296/385/GCF_041296385.1_UR_Acoe_1.0/GCF_041296385.1_UR_Acoe_1.0_genomic.fna.gz
gunzip GCF_041296385.1_UR_Acoe_1.0_genomic.fna.gz
grep ">" GCF_041296385.1_UR_Acoe_1.0_genomic.fna | head

# Cyanocitta cristata
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/046/129/655/GCA_046129655.1_bCyaCrs1.hap1/GCA_046129655.1_bCyaCrs1.hap1_genomic.fna.gz
gunzip GCA_046129655.1_bCyaCrs1.hap1_genomic.fna.gz
grep ">" GCA_046129655.1_bCyaCrs1.hap1_genomic.fna | head
```

Subset genomes to the two first chromosomes to avoid unnecessary data volumes. 
```sh
# Corvus hawaiiensis
awk '/^>NC_063213.1/{flag=1;print;next} /^>/{flag=0} flag' GCF_020740725.1_bCorHaw1.pri.cur_genomic.fna > corvus_hawaiiensis.fna
awk '/^>NC_063214.1/{flag=1;print;next} /^>/{flag=0} flag' GCF_020740725.1_bCorHaw1.pri.cur_genomic.fna >> corvus_hawaiiensis.fna

# Corvus moneduloides
awk '/^>NC_045476.1/{flag=1;print;next} /^>/{flag=0} flag' GCF_009650955.1_bCorMon1.pri_genomic.fna > corvus_moneduloides.fna
awk '/^>NC_045477.1/{flag=1;print;next} /^>/{flag=0} flag' GCF_009650955.1_bCorMon1.pri_genomic.fna >> corvus_moneduloides.fna

# Corvus cornix
awk '/^>NC_046332.1/{flag=1;print;next} /^>/{flag=0} flag' GCF_000738735.6_ASM73873v6_genomic.fna > corvus_cornix.fna
awk '/^>NC_046333.1/{flag=1;print;next} /^>/{flag=0} flag' GCF_000738735.6_ASM73873v6_genomic.fna >> corvus_cornix.fna

# Aphelocoma coerulescens
awk '/^>NC_091013.1/{flag=1;print;next} /^>/{flag=0} flag' GCF_041296385.1_UR_Acoe_1.0_genomic.fna > aphelocoma_coerulescens.fna
awk '/^>NC_091014.1/{flag=1;print;next} /^>/{flag=0} flag' GCF_041296385.1_UR_Acoe_1.0_genomic.fna >> aphelocoma_coerulescens.fna

# Cyanocitta cristata
awk '/^>CM100531.1/{flag=1;print;next} /^>/{flag=0} flag' GCA_046129655.1_bCyaCrs1.hap1_genomic.fna > cyanocitta_cristata.fna
awk '/^>CM100532.1/{flag=1;print;next} /^>/{flag=0} flag' GCA_046129655.1_bCyaCrs1.hap1_genomic.fna >> cyanocitta_cristata.fna
```

Bacterial taxa

Enterococcus faecalis	
Corynebacterium ciconiae	
Escherichia coli	
Bacteroides caccae	
Bacteroides uniformis	
Oceanisphaera sp012518835	
Akkermansia muciniphila	
Phocaeicola vulgatus	
Mediterranea pullorum	
Sellimonas intestinalis	
Scybalousia sp900543675	
Phocaeicola vulgatus	
Mammaliicoccus sciuri	
Bradyrhizobium sp001908235	
Clostridium_H haemolyticum	
Lactobacillus johnsonii	
Sutterella sp905215795	

### Simulation of longitudinal progression of mice

Prepare the working directory
```sh
mkdir holosimulator_mice
cd holosimulator_mice
```

#### Fetch and subset host genome

Download the genomes and check contig headers.
```sh
# Corvus hawaiiensis
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/947/179/515/GCF_947179515.1_mApoSyl1.1/GCF_947179515.1_mApoSyl1.1_genomic.fna.gz
gunzip GCF_947179515.1_mApoSyl1.1_genomic.fna.gz
grep ">" GCF_947179515.1_mApoSyl1.1_genomic.fna | head
```

Subset genome to the two first chromosomes to avoid unnecessary data volumes. 
```sh
# Corvus hawaiiensis
awk '/^>NC_067472.1/{flag=1;print;next} /^>/{flag=0} flag' GCF_947179515.1_mApoSyl1.1_genomic.fna > apodemus_sylvaticus.fna
awk '/^>NC_067473.1/{flag=1;print;next} /^>/{flag=0} flag' GCF_947179515.1_mApoSyl1.1_genomic.fna >> apodemus_sylvaticus.fna
```

Bacterial taxa

Helicobacter_C japonicus
Yersinia bercovieri
Duncaniella dubosii
Choladocola sp009774145
UBA3282 sp910575245
Duncaniella sp002494015	
Enterococcus_B hirae	
Merdisoma sp910574255	
Muribaculum intestinale	
Eubacterium_R sp011958665	
CAG-873 sp910584905	
Paramuribaculum sp001689565	
Akkermansia muciniphila	
Sporofaciens sp000403455	
Ligilactobacillus murinus
Bacteroides sp910585845	
Ligilactobacillus apodemi	
Novisyntrophococcus sp900539175	
Brachyspira sp910578695	