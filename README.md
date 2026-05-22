![MaraudersGenoMap-logo](MaraudersGenoMap.png)

# Marauders GenoMap

Marauders GenoMap is a Python/PyQt6 graphical interface for downloading SRA data, performing genomic/transcriptomic assembly, searching for protein domains using HMMER, and exploring results in folders organized by SRA accession number.

The program was designed for local execution. The GUI calls scripts from the `scripts/` folder, and each step creates subfolders within the provided SRA code folder, for example: `SRRXXXXXX/01_QC_Reports`, `SRRXXXXXX/05_Assembly_Results`, and `SRRXXXXXX/02_HMMER_Results`.

## Docker and Singularity/Apptainer

The project includes files to generate a local image with the main dependencies:

- `Dockerfile`: Docker image based on `debian:trixie-slim`.

- `.dockerignore`: Prevents copying virtual environments, SRA data, FASTQ, BAM/SAM, FASTA, results, and caches into the image.

- `Makefile`: Shortcuts for build and execution.

The image includes MEGAHIT, SPAdes, and Trinity, as well as Salmon, SRA Toolkit, FastQC, Trimmomatic, BBMap, MultiQC, Prodigal, HMMER, SeqTK, Bowtie2, Samtools, PyQt6, Pandas, Matplotlib, and Biopython.

### Obtain Docker image from DockerHub

To download the latest version, run the following command in your terminal:

Bash

```
docker pull evomol/marauders-genomap
```

And run the image:

Bash

```
docker run --rm -it \
    --network host \
    --user $(id -u):$(id -g) \
    -e DISPLAY=$DISPLAY \
    -e XAUTHORITY=/tmp/.Xauthority \
    -e HOME=/tmp \
    -e QT_X11_NO_MITSHM=1 \
    -v /tmp/.X11-unix:/tmp/.X11-unix:rw \
    -v "${XAUTHORITY:-$HOME/.Xauthority}":/tmp/.Xauthority:ro \
    -v "$(pwd)":/data \
    -w /data \
    evomol/marauders-genomap:latest bash -c "python3 /app/marauders.py"
```

### Build Docker image locally

Bash

```
make docker-build
```

This target executes the build with a defined name/tag:

Bash

```
docker build --no-cache -t marauders-genomap .
```

Using `-t marauders-genomap:latest` is important because it provides a predictable name for the final image. Without this tag, Docker still creates the image, but it may remain unnamed, making subsequent execution difficult.

The `--no-cache` option prevents reusing cached layers when the base or its dependencies change.

To perform a faster build reusing the cache:

Bash

```
make docker-build-cache
```

This target is equivalent to:

Bash

```
docker build -t marauders-genomap .
```

After the build, it is normal for at least two images to appear when listing with `docker images` or graphical tools:

- `debian:trixie-slim`: Base image downloaded by Docker, small, used to build the final image;

- `marauders-genomap:latest`: Final Marauders GenoMap image; this is the image that should be executed.

The `debian:trixie-slim` image can remain on the system without issue. It occupies little space and accelerates future builds. If you want to remove it, use:

Bash

```
docker rmi debian:trixie-slim
```

### Running the GUI with Docker

Before opening the interface, allow local access to the X11 server:

Bash

```
xhost +local:docker
```

Then, enter the project folder where the `Makefile` is located and start the final `marauders-genomap:latest` image:

Bash

```
cd /path/to/your/project/folder
make docker-run
```

The `make docker-run` command must be executed within the Marauders GenoMap folder. The `make` command looks for the `Makefile` in the current directory, and this target mounts the current folder (`$PWD`) to `/app` inside the container.

If you want to start the image from any other directory, use the full Docker command pointing to the project folder:

Bash

```
docker run --rm -it \
  --network host \
  --user $(id -u):$(id -g) \
  -e DISPLAY=$DISPLAY \
  -e XAUTHORITY=/tmp/.Xauthority \
  -e HOME=/tmp \
  -e QT_X11_NO_MITSHM=1 \
  -v /tmp/.X11-unix:/tmp/.X11-unix:rw \
  -v "${XAUTHORITY:-$HOME/.Xauthority}:/tmp/.Xauthority:ro" \
  -v "/path/to/your/project/folder":/app \
  -w /app \
  marauders-genomap:latest
```

The `make docker-run` target executes the same concept, automatically using the current folder (you must be in the Marauders folder).

`docker-run` mounts the current folder to `/app`, so downloaded files and results remain in the project directory rather than in a temporary container layer.

When finished, you can revoke permission:

Bash

```
xhost -local:docker
```

### Using with Singularity/Apptainer

The project does not need to maintain a separate `Singularity.def` if the Docker image is published on Docker Hub or another compatible registry. Singularity/Apptainer can download a Docker image directly and convert it to `.sif` locally.

With Singularity:

Bash

```
singularity pull marauders-genomap.sif docker://evomol/marauders-genomap
```

Then run:

Bash

```
singularity run \
  --writable-tmpfs \
  --bind "$PWD":/app \
  --bind /tmp/.X11-unix:/tmp/.X11-unix \
  --bind "${XAUTHORITY:-$HOME/.Xauthority}":/tmp/.Xauthority:ro \
  --env DISPLAY=$DISPLAY,QT_X11_NO_MITSHM=1,XAUTHORITY=/tmp/.Xauthority,HOME=/tmp \
  marauders-genomap.sif
```

With Apptainer:

Bash

```
apptainer pull marauders-genomap.sif docker://evomol/marauders-genomap
```

Then run:

Bash

```
apptainer run \
  --writable-tmpfs \
  --bind "$PWD":/app \
  --bind /tmp/.X11-unix:/tmp/.X11-unix \
  --bind "${XAUTHORITY:-$HOME/.Xauthority}":/tmp/.Xauthority:ro \
  --env DISPLAY=$DISPLAY,QT_X11_NO_MITSHM=1,XAUTHORITY=/tmp/.Xauthority,HOME=/tmp \
  marauders-genomap.sif
```

Docker Hub hosts Docker images. It is not, in general, a native repository for `.sif` files, but it can serve the same Docker image to both Docker and Singularity/Apptainer users.

## How to run locally without a container

If the dependencies are already installed on the system:

Bash

```
make run
```

The `install` target creates a local command called `marauders`, pointing to the current folder:

Bash

```
make install
```

To remove:

Bash

```
make uninstall
```

## Program Tab Flow

The program has 5 main tabs. The SRA ID provided in the first tab is used by the subsequent steps to automatically locate inputs and results.

### 1. Download

Function: download raw SRA data and prepare compressed FASTQ files.

Script used:

- `scripts/SRAget.sh`

What the tab does:

- Receives the SRA accession;

- Allows selecting the sample type;

- For transcriptomics/metatranscriptomics, allows choosing paired-end or single-end layout;

- Creates a folder with the accession name;

- Executes `prefetch`;

- Executes `fasterq-dump` or `fastq-dump` as a fallback;

- Compresses FASTQ files with `pigz`.

Expected files after this step:

- `SRA_ID/SRA_ID_1.fastq.gz` and `SRA_ID/SRA_ID_2.fastq.gz` for paired-end;

- `SRA_ID/SRA_ID.fastq.gz` for single-end.

### 2. Assembly

Function: execute quality control, trimming, normalization, MultiQC, and assembly.

Script used:

- `scripts/run_assembly.sh`

What the tab does:

- Selects the assembler: MEGAHIT, SPAdes, or Trinity;

- Adjusts threads and RAM;

- Shows a recommended GB per thread;

- Allows removing raw files and/or the trimming folder at the end;

- Runs FastQC, Trimmomatic, BBNorm, MultiQC, and the chosen assembler.

Interface rules:

- MEGAHIT is recommended for metagenomes and fast execution.

- SPAdes is available for isolate genomics.

- Trinity is available for transcriptomics and metatranscriptomics.

Folders created within the SRA folder:

- `01_QC_Reports`

- `02_Trimmed_Reads`

- `03_Normalized_Reads`

- `04_MultiQC_Report`

- `05_Assembly_Results`

Main outputs:

- MEGAHIT: `05_Assembly_Results/MEGAHIT_*/final.contigs.fa`

- SPAdes: `05_Assembly_Results/SPADES_*/scaffolds.fasta` or `contigs.fasta`

- Trinity: `*.Trinity.fasta` files

### 3. Protein Search

Function: predict proteins from contigs and search for domains using an HMM profile.

Script used:

- `scripts/ProtSearch.sh`

What the tab does:

- Receives a `.hmm` file;

- Automatically detects the contigs file generated by the selected assembler;

- Chooses the Prodigal mode:
  
  - `single` for genomics;
  
  - `meta` for metagenomics/metatranscriptomics;

- Executes Prodigal to generate predicted proteins;

- Executes `hmmsearch` to generate the hits table.

Folders created within the SRA folder:

- `01_Predicted_Proteins`

- `02_HMMER_Results`

Main outputs:

- `01_Predicted_Proteins/predicted_proteins.faa`

- `02_HMMER_Results/*_domain_hits.tbl`

### 4. Get Results

Function: extract protein sequences associated with the hits detected by HMMER.

Script used:

- `scripts/get_Seq_results.sh`

What the tab does:

- Automatically locates the `.tbl` table in `02_HMMER_Results`;

- Locates `01_Predicted_Proteins/predicted_proteins.faa`;

- Extracts unique IDs from the first column of the HMMER table;

- Uses `seqtk subseq` to recover the corresponding sequences.

Main outputs:

- Temporary `protein_ids_to_extract.txt`;

- `lectin_hits.faa` with the extracted sequences.

### 5. Analysis & Cutting

Function: visualize hits, select targets, and extract/cut specific sequences.

Scripts used:

- This tab uses internal `marauders.py` functions;

- It does not call an external script directly.

What the tab does:

- Automatically loads a `.tbl` table from the `SRA_ID/02_HMMER_Results` folder;

- Also allows manual selection of a `.tbl` table;

- Filters hits with `full_E-value < 0.0001`;

- Shows the top 25 hits;

- Clicking a row transfers the `target_name` to the Severus Snap(e) extractor;

- Automatically locates `SRA_ID/01_Predicted_Proteins/predicted_proteins.faa` or accepts a manual FASTA;

- Saves sequences found by the provided term;

- When possible, cuts the region using coordinates present in the FASTA description;

- Generates a sequence length histogram with Matplotlib.

## Scripts in the `scripts/` folder

### `SRAget.sh`

Downloads and extracts SRA reads. Uses `prefetch`, `fasterq-dump`, fallback with `fastq-dump`, and compression with `pigz`.

### `run_assembly.sh`

Executes the assembly pipeline. Uses FastQC, Trimmomatic, BBNorm, MultiQC, and one of the assemblers: MEGAHIT, SPAdes, or Trinity.

### `ProtSearch.sh`

Predicts proteins with Prodigal and searches for domains with HMMER (`hmmsearch`).

### `get_Seq_results.sh`

Extracts protein IDs from the HMMER table and recovers corresponding sequences with SeqTK.

### `map_reads.sh`

Auxiliary script to map reads against a reference FASTA. Uses Bowtie2 and Samtools, creating results in `06_Mapping_Results`. Currently, it is not called directly by any of the 5 main GUI tabs.

### `setup.sh`

Historical local installation script via `apt`/`pip`. With Docker or Singularity, dependency installation is handled by the image, so this script tends to be unnecessary for container users.

## Output Structure

A complete flow for an SRA accession tends to create a structure like:

Plaintext

```
SRA_ID/
├── 01_QC_Reports/
├── 02_Trimmed_Reads/
├── 03_Normalized_Reads/
├── 04_MultiQC_Report/
├── 05_Assembly_Results/
├── 01_Predicted_Proteins/
├── 02_HMMER_Results/
└── lectin_hits.faa
```

Some folders may vary according to the chosen options, sample type, and step executed.

## Notes

- The program depends on an X11 graphical display when executed in Docker.

- Results are generated within the project folder when execution uses `make docker-run`, as the current directory is mounted to `/app`.

- Trinity significantly increases the image size but allows keeping the transcriptomics workflow in the same image as MEGAHIT and SPAdes.

- If Trinity complains that `salmon` is not installed, rebuild the updated image. `Salmon` is a runtime dependency of Trinity and must be present in the Docker/Singularity image.

- MultiQC warnings about missing `assets/js/packages/*.js` files may appear in some Debian-packaged versions. When the log ends with `MultiQC complete`, these warnings do not interrupt the pipeline; they only affect features embedded in the HTML report.

## Development Team

Marauder's GenoMap was developed by Djorkaeff Oliveira, Rodrigo Orvate, João Pedro Lemos e Silva Rodrigues, and João Paulo MS Lima.
[EvoMol-Lab](evomol-lab.imd.ufrn.br), BioMe, UFRN, Brazil.
