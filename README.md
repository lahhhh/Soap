# SOAP, a desktop platform for single-cell omics data analysis
## Version 1.0
- Soap is a Windows application with a graphic user interface, designed for data analysis of single-cell RNA-seq or ATAC-seq

- All operations are provided via GUI.

- If you want to add new functions, just add a new action and implement corresponding slots somewhere. 

- This project is built in MSVC2022 community edition.

## Installation
- 1. Download the release zip package.
- 2. unzip it wherever you want (at your convenience).
- 3. ok

## Vignette
- 1. double click the `soap.exe`, the main interface will appear.
- 2. find the `File` menu at top left (use your mouse).
- 3. choose which type you want to analyse (depend on your data type).
- 4. click the corresponding action (like `scRNA-seq` >> `Load 10X`).
- 5. choose the right file in the dialog (for example, if you choose 10x scRNA-seq, you should choose `barcodes.tsv.gz` for `Barcodes File`
, `features.tsv.gz` for `Feature File` and `matrix.mtx.gz` for `Matrix File`).
- 6. waiting for data loading (usually not long, data loading will fail if the file is broken).
- 7. analyse through the right-click menu.

## System Requirement
- 1. Windows operating system.
- 2. Intel CPU (because we use Intel Math Kernel Library), I don't know if it can run appropriately on AMD CPU (maybe need openblas).
- 3. CPU should support AVX2 archive (CPU produced after 2012).
- 4. enough memory (the more cells you analyze, the more memory you need).


If you meet bugs in using SOAP, please kindly report the issue.

Contact the author: liuanhang1997@qq.com

