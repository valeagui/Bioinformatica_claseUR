# Clase de filogenómica

Primero entramos al cluster de la universidad y corremos salloc. 

## Obtención de los datos: 
Descargamos directamente en el cluster los datos correspondientes (https://academic.oup.com/sysbio/article/65/6/1057/2281640) Estos hacen referencia a los transcriptomas de 23 especies de vertebrados. 

## Activación del ambiente: 
Usamos este código para la activación: ```conda activate biopt```.Posteriormente clonamos el repositorio (https://github.com/iirisarri/EGT3C_2019_PHYLO.git) con ```git clone```. 
Debido a que los datos se encuentran comprimidos con .tar, utilizando ```tar -xvf archivo.tar```, los descomprimimos para poder manejarlos.

## Grupos ortólogos: 
Utilizando el programa de OrthoFinder buscamos las regiones ortólogas entre los transcriptomas de la carpeta, generando un archivo de secuencia para cada uno de ellos: ```orthofinder -os -M msa -S blast -f vertebrate_proteomes```

Entrando a la carpeta generada de Orthogroups, se encuentra una lista de secuencias que pertenecen a un mismo orthogroup bajo el nombre de _Orthogroups.tvs_. Las secuencias se encuentran en la carpeta de _Orthogroup_Sequence_.

Posteriormente, cambiamos los nombres del patrón de género_especie_gen_xxx a género_especie con: ```sed -E -i 's/_GENE_[0-9]+//g' Orthogroups_copia.tsv``` y ```for f in *fa; do sed -E 's/_GENE.+//g' $f > out; mv out $f; done```. 
![image](https://github.com/valeagui/Bioinformatica_claseUR/assets/127573975/580b7097-6f97-493c-961b-392320e3c637)

## Revisar la calidad:
Con el fin de eliminar los errores generados por el ensamblaje o la anotación, utilizamos el programa de PREQUAL, el cual identifica fragmentos que no comparten homología para posteriormente ser removidas con: ```for f in *fa; do prequal $f ; done```, generando archivo _.filtered_

## Alineamiento múltiple:
Utilizando el programa de MAFFT generamos el alineamiento de cada orthogroup con: ```for f in *filtered; do mafft $f > $f.mafft; done```

![image](https://github.com/valeagui/Bioinformatica_claseUR/assets/127573975/01636b23-9b2e-42e2-9e0f-140e471629c2)
![image](https://github.com/valeagui/Bioinformatica_claseUR/assets/127573975/c298f5ab-dac2-4b38-816e-234142687741)

## Alignment trimming: 
En este caso borramos las regiones que tienen muchos gaps, lo realizamos eliminando las posiciones con más de 80% de gaps con: ```for f in *mafft; do bmge -i $f -t AA -g 0.8 -h 1 -w 1 -of $f.g08.fas; done```. 

![image](https://github.com/valeagui/Bioinformatica_claseUR/assets/127573975/18023960-ebd4-4d01-9bb4-8b3b836c5c64)
![image](https://github.com/valeagui/Bioinformatica_claseUR/assets/127573975/bbd21912-30a5-40d1-b625-60ba9ce912af)

## Alineamiento concatenado:
Copiamos el archivo de _FASconCAT-G_v1.05.1.pl_ a la carpeta. Con FASconCAT unimos todas las secuencias de los archivos por separado con ```perl ./FASconCAT-G_v1.05.1.pl -l -s```. Con esto generamos la matriz de diferencias y el archivo con las particiones. 
![image](https://github.com/valeagui/Bioinformatica_claseUR/assets/127573975/582c2875-b6bd-4a05-86bb-ca3e8fc47402)

## Árbol de Máxima verosimilitud:
Utilizamos el programa de IQTREE. En la primera parte generamos el alineamiento concatenado con un solo locus que evoluciona de una misma forma. Generamos un sbatch con 10 núcleos y 20 minutos de tiempo y con la línea de código:
![image](https://github.com/valeagui/Bioinformatica_claseUR/assets/127573975/97b54f13-2a17-4edf-8caf-4a74f78ec992)

Posteriormente, generamos otro sbatch en el cual el código le indica al programa donde inicia y termina cada gen de cada una de las particiones con: ```module load iqtree iqtree -s FcC_supermatrix.fas -spp FcC_supermatrix_partition.txt -m TEST -msub nuclear -merit AICc -bb 1000 -alrt 1000 -nt 1 -bnni -pre partitioned```. Aca se generaron los archivos de partitioned.treefile y  unpartitioned.treefile.

## Análisis de coalescencia:
Para realizar el árbol usamos el programa de ASTRAL. Primero generamos los árboles de genes: ```for f in *filtered.mafft.g08.fas; do iqtree -s $f -m TEST -msub nuclear -merit AICc -nt 1; done```. Luego colocamos todos los árboles en un mismo archivo ```cat *filtered.mafft.g08.fas.treefile > my_gene_trees.tre```

Corrimos el ASTRAL con: ```java -jar /home/bio.pt/data/Astral/astral.5.7.8.jar -i my_gene_trees.tre -o species_tree_ASTRAL.tre 2> out.log```. Generando el archivo de _species_tree_ASTRAL.tre_

Para ver los árboles usamos FigTree y lo enraizamos al punto medio, pasamos los archivos al computador con:

```scp -r -i bio.pt.pem -P 53841 bio.pt@loginpub-hpc.urosario.edu.co:/home/bio.pt/data/Vale_Aguilar/Filogenomica/phylogenomics/vertebrate_proteomes/OrthoFinder/Results_May12/Orthogroup_Sequences/concatenado/arboles /mnt/c/Users/valea/OneDrive/Documents/Bioinformatica```

### Árbol de particiones:
![image](https://github.com/valeagui/Bioinformatica_claseUR/assets/127573975/00888252-5508-4c69-b73c-058daff9d3e2)
[Texto del enlace](partitioned.treefile)

### Árbol sin particiones:
![image](https://github.com/valeagui/Bioinformatica_claseUR/assets/127573975/1c1235e2-0d90-4aed-9fa0-e09473035957)

En el caso del árbol con particiones y sin particiones, las relaciones entre las ramas es la misma, junto con una diferencia en el número de sustituciones similar. Los cambios más pronunciados se ven en el valor de soporte de cada rama (Bootstrap), estos suelen ser más altos en el caso del árbol generado con particiones. Esto puedo indicar que el realizar las particiones, delimitando cada uno de los genes, puede llegar a ser más efectivo en el momento de la construcción del árbol.

Dichas diferencias pueden ser causadas porque al usar las particiones se puede reflejar de una mejor forma las diferencias en el grado de las fuerzas o eventos evolutivos, lo que hace el análisis más sensible y exacto, en comparación al árbol generado sin las particiones que no toma en cuenta estas diferencias y trata toda la secuencia de la misma forma.

En ambos árboles se puede decir que se forman 2 clados principales. El clado superior se puede considerar monofilético porque agrupa peces óseos de la clase actinpterygii. El clado inferior o grupo inferior se puede considerar casi como polifilético debido a que, el grupo de los peces que pertenecen a la clase sarcopteryggi se dividen en dos, y en la “mitad” de estos se encuentran especies agrupadas en reptiles, mamíferos y anfibios.

### Árbol de ASTRA: 
![image](https://github.com/valeagui/Bioinformatica_claseUR/assets/127573975/04ef3c53-813d-4ee9-b20c-297133cbd1f9)

En este caso se busca el punto en el que las secuencias diferentes concuerdan o encuentran un ancestro en común, retrocediendo desde las especies hacia el ancestro. Por la forma en la que está construida el árbol se puede decir que tiene un mayor soporte como se ve reflejado en la mayoría de los Bootstrap a excepción de algunas ramas.

En comparación con los árboles de máxima verosimilitud cambian, se podría decir de forma significativa, las relaciones entre las especies. Esto puede tener varias explicaciones, una de ellas es que los árboles de coalescencia son más sensibles a los eventos de reticulación o hibridación y coalescencia, en la historia evolutiva de las especies y de alguna forma esto se tiene en cuenta en la formación del árbol.

En este caso se puede decir que se forman aproximadamente dos grupos o clados. El grupo de arriba podemos decir que es monofilético porque agrupa todas las especies, en este caso, anfibios, reptiles y mamíferos, juntos. El grupo de la mitad es un grupo parafilético porque no agrupa todos los miembros de la clase de sarcopteryggi. Por último, el grupo de abajo es polifilético porque dentro de las especies de la clase Actinopterygii incluye un individuo de la clase Chondrichthyes y Sarcopteryggi. 
