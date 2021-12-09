# dataset-edna
By [Diana LaScala-Gruenewald](https://github.com/dianalg)

## Introduction
**Rationale:**

DNA derived data are increasingly being used to document taxon 
occurrences. To ensure these data are useful to the broadest possible 
community, [GBIF](https://www.gbif.org/) published a guide entitled "[Publishing DNA-derived 
data through biodiversity data platforms](https://docs.gbif-uat.org/publishing-dna-derived-data/1.0/en/)." 
This guide is supported by the [DNA derived data extension](https://tools.gbif.org/dwca-validator/extension.do?id=http://rs.gbif.org/terms/1.0/DNADerivedData) 
for [Darwin Core](https://dwc.tdwg.org/), which incorporates [MIxS](https://gensc.org/mixs/) 
terms into the Darwin Core standard. 

This use case draws on both the guide and the extension to illustrate 
how to incorporate a DNA derived data extension file into a Darwin Core
archive. 

For further information on this use case and the DNA Derived data extension
in general, see the recording of the [OBIS Webinar on Genetic Data](https://obis.org/2021/10/13/gendatawebinar/).

**Project abstract:**

The example data employed in this use case are from marine 
filtered seawater samples collected at a nearshore station in 
Monterey Bay, California, USA. They were collected by CTD
rosette and filtered by a peristaltic pump system. Subsequently, 
they underwent metabarcoding for the 18S V9 region. The resulting
ASVs, their assigned taxonomy, and the metadata associated with their
collection are the input data for the conversion scripts 
presented here.

A selection of samples from this collection were included in the 
publication "[Environmental DNA reveals seasonal shifts and potential 
interactions in a marine community](https://www.nature.com/articles/s41467-019-14105-1)"
 which was published with open access in *Nature Communications* in 2020.

**Contacts:**
- Francisco Chavez - Principle Investigator ([chfr@mbari.org](chfr@mbari.org))
- Kathleen Pitz - Research Associate ([kpitz@mbari.org](kpitz@mbari.org))
- Diana LaScala-Gruenewald - Point of Contact ([dianalg@mbari.org](dianalg@mbari.org))

## Published data
- [GBIF](https://www.gbif.org/dataset/e0b59ee7-19ae-4eb0-9217-33317fb50d47)
- [OBIS](https://obis.org/dataset/62b97724-da17-4ca7-9b26-b2a22aeaab51)

## Repo structure
```
.
+-- README.md                   :Description of this repository
+-- LICENSE                     :Repository license
+-- .gitignore                  :Files and directories to be ignored by git
|
+-- raw
|   +-- asv_table.csv           :Source data containing ASV sequences and number of reads
|   +-- taxa_table.csv          :Source data containing taxon matches for each ASV
|   +-- metadata_table.csv      :Source data containing metadata about samples (e.g. collection information)
|
+-- src
|   +-- conversion_code.py      :Darwin Core mapping script
|   +-- conversion_code.ipynb   :Darwin Core mapping Jupyter Notebook
|
+-- processed
|   +-- occurrence.csv          :Occurrence file, generated by conversion_code
|   +-- dna_extension.csv       :DNA Derived Data Extension file, generated by conversion_code

```

