# Marauders GenoMap

Marauders GenoMap e uma interface grafica em Python/PyQt6 para baixar dados SRA, executar montagem genomica/transcriptomica, buscar dominios proteicos com HMMER e explorar os resultados gerados em pastas organizadas por accession SRA.

O programa foi pensado para execucao local. A GUI chama os scripts da pasta `scripts/` e cada etapa cria subpastas dentro da pasta do codigo SRA informado, por exemplo `SRRXXXXXX/01_QC_Reports`, `SRRXXXXXX/05_Assembly_Results` e `SRRXXXXXX/02_HMMER_Results`.

## Docker e Singularity/Apptainer

O projeto inclui arquivos para gerar uma imagem local com as dependencias principais:

- `Dockerfile`: imagem Docker baseada em `debian:trixie-slim`.
- `.dockerignore`: evita copiar ambientes virtuais, dados SRA, FASTQ, BAM/SAM, FASTA, resultados e caches para dentro da imagem.
- `Makefile`: atalhos para build e execucao.

A imagem inclui MEGAHIT, SPAdes e Trinity, alem de Salmon, SRA Toolkit, FastQC, Trimmomatic, BBMap, MultiQC, Prodigal, HMMER, SeqTK, Bowtie2, Samtools, PyQt6, Pandas, Matplotlib e Biopython.

### Obter imagem Docker no DockerHub

Para baixar a versão mais recente, execute o comando abaixo no seu terminal:

```bash
docker pull evomol/marauders-genomap
```

e rodar a imagem:

```bash
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

### Gerar imagem Docker local

```bash
make docker-build
```

Esse alvo executa o build com nome/tag definido:

```bash
docker build --no-cache -t marauders-genomap .
```

O uso de `-t marauders-genomap:latest` e importante porque da um nome previsivel para a imagem final. Sem essa tag, o Docker ainda cria a imagem, mas ela pode ficar sem nome claro, dificultando a execucao depois.

O `--no-cache` evita reaproveitar camadas antigas quando a base ou as dependencias mudam.

Para fazer um build mais rapido reaproveitando cache:

```bash
make docker-build-cache
```

Esse alvo equivale a:

```bash
docker build -t marauders-genomap .
```

Depois do build, e normal aparecerem pelo menos duas imagens ao listar com `docker images` ou ferramentas graficas:

- `debian:trixie-slim`: imagem base baixada pelo Docker, pequena, usada para construir a imagem final;
- `marauders-genomap:latest`: imagem final do Marauders GenoMap, esta e a imagem que deve ser executada.

A imagem `debian:trixie-slim` pode ficar no sistema sem problema. Ela ocupa pouco espaco e acelera builds futuros. Se quiser remover, use:

```bash
docker rmi debian:trixie-slim
```

### Rodar a GUI com Docker

Antes de abrir a interface, permita acesso local ao servidor X11:

```bash
xhost +local:docker
```

Depois, entre na pasta do projeto, onde esta o `Makefile`, e inicie a imagem final `marauders-genomap:latest`:

```bash
cd /caminho/para/sua/pasta/de/projetos
make docker-run
```

O `make docker-run` precisa ser executado dentro da pasta do Marauders GenoMap. O `make` procura o `Makefile` no diretorio atual, e esse alvo monta a pasta atual (`$PWD`) em `/app` dentro do container.

Se quiser iniciar a imagem de qualquer outro diretorio, use o comando Docker completo apontando para a pasta do projeto:

```bash
docker run --rm -it \
  --network host \
  --user $(id -u):$(id -g) \
  -e DISPLAY=$DISPLAY \
  -e XAUTHORITY=/tmp/.Xauthority \
  -e HOME=/tmp \
  -e QT_X11_NO_MITSHM=1 \
  -v /tmp/.X11-unix:/tmp/.X11-unix:rw \
  -v "${XAUTHORITY:-$HOME/.Xauthority}:/tmp/.Xauthority:ro" \
  -v "/caminho/para/sua/pasta/de/projetos":/app \
  -w /app \
  marauders-genomap:latest
```

O alvo `make docker-run` executa a mesma ideia, usando automaticamente a pasta atual (deve estar na pasta do Marauders).

O `docker-run` monta a pasta atual em `/app`, entao os arquivos baixados e resultados permanecem no diretorio do projeto, nao dentro de uma camada temporaria do container.

Ao terminar, voce pode revogar a permissao:

```bash
xhost -local:docker
```

### Usar com Singularity/Apptainer

O projeto nao precisa manter um `Singularity.def` separado se a imagem Docker for publicada no Docker Hub ou em outro registry compativel. Singularity/Apptainer conseguem baixar uma imagem Docker diretamente e converter para `.sif` localmente.

Com Singularity:

```bash
singularity pull marauders-genomap.sif docker://evomol/marauders-genomap
```

Depois rode:

```bash
singularity run \
  --writable-tmpfs \
  --bind "$PWD":/app \
  --bind /tmp/.X11-unix:/tmp/.X11-unix \
  --bind "${XAUTHORITY:-$HOME/.Xauthority}":/tmp/.Xauthority:ro \
  --env DISPLAY=$DISPLAY,QT_X11_NO_MITSHM=1,XAUTHORITY=/tmp/.Xauthority,HOME=/tmp \
  marauders-genomap.sif
```

Com Apptainer:

```bash
apptainer pull marauders-genomap.sif docker://evomol/marauders-genomap
```

Depois rode:

```bash
apptainer run \
  --writable-tmpfs \
  --bind "$PWD":/app \
  --bind /tmp/.X11-unix:/tmp/.X11-unix \
  --bind "${XAUTHORITY:-$HOME/.Xauthority}":/tmp/.Xauthority:ro \
  --env DISPLAY=$DISPLAY,QT_X11_NO_MITSHM=1,XAUTHORITY=/tmp/.Xauthority,HOME=/tmp \
  marauders-genomap.sif
```

O Docker Hub hospeda imagens Docker. Ele não é, em geral, um repositorio nativo de arquivos `.sif`, mas pode servir a mesma imagem Docker para usuarios Docker e tambem para usuarios Singularity/Apptainer.

## Como Executar Localmente Sem Container

Se as dependencias ja estiverem instaladas no sistema:

```bash
make run
```

O alvo `install` cria um comando local chamado `marauders`, apontando para a pasta atual:

```bash
make install
```

Para remover:

```bash
make uninstall
```

## Fluxo das Abas do Programa

O programa possui 5 abas principais. O ID SRA informado na primeira aba e usado pelas demais etapas para localizar entradas e resultados automaticamente.

### 1. Download

Funcao: baixar os dados brutos do SRA e preparar os FASTQ compactados.

Script usado:

- `scripts/SRAget.sh`

O que a aba faz:

- recebe o accession SRA;
- permite escolher o tipo de amostra;
- para transcriptomica/metatranscriptomica, permite escolher layout paired-end ou single-end;
- cria uma pasta com o nome do accession;
- executa `prefetch`;
- executa `fasterq-dump` ou `fastq-dump` como fallback;
- compacta os FASTQ com `pigz`.

Arquivos esperados apos a etapa:

- `SRA_ID/SRA_ID_1.fastq.gz` e `SRA_ID/SRA_ID_2.fastq.gz` para paired-end;
- `SRA_ID/SRA_ID.fastq.gz` para single-end.

### 2. Assembly

Funcao: executar controle de qualidade, trimming, normalizacao, MultiQC e montagem.

Script usado:

- `scripts/run_assembly.sh`

O que a aba faz:

- seleciona o montador: MEGAHIT, SPAdes ou Trinity;
- ajusta threads e RAM;
- mostra uma recomendacao de GB por thread;
- permite remover arquivos brutos e/ou a pasta de trimming ao final;
- roda FastQC, Trimmomatic, BBNorm, MultiQC e o montador escolhido.

Regras da interface:

- MEGAHIT e recomendado para metagenomas e execucao rapida;
- SPAdes fica disponivel para genomica de isolados;
- Trinity fica disponivel para transcriptomica e metatranscriptomica.

Pastas criadas dentro da pasta SRA:

- `01_QC_Reports`
- `02_Trimmed_Reads`
- `03_Normalized_Reads`
- `04_MultiQC_Report`
- `05_Assembly_Results`

Saidas principais:

- MEGAHIT: `05_Assembly_Results/MEGAHIT_*/final.contigs.fa`
- SPAdes: `05_Assembly_Results/SPADES_*/scaffolds.fasta` ou `contigs.fasta`
- Trinity: arquivos `*.Trinity.fasta`

### 3. Protein Search

Funcao: predizer proteinas a partir dos contigs e buscar dominios usando um perfil HMM.

Script usado:

- `scripts/ProtSearch.sh`

O que a aba faz:

- recebe um arquivo `.hmm`;
- detecta automaticamente o arquivo de contigs gerado pelo montador selecionado;
- escolhe o modo do Prodigal:
  - `single` para genomica;
  - `meta` para metagenomica/metatranscriptomica;
- executa Prodigal para gerar proteinas preditas;
- executa `hmmsearch` para gerar a tabela de hits.

Pastas criadas dentro da pasta SRA:

- `01_Predicted_Proteins`
- `02_HMMER_Results`

Saidas principais:

- `01_Predicted_Proteins/predicted_proteins.faa`
- `02_HMMER_Results/*_domain_hits.tbl`

### 4. Get Results

Funcao: extrair sequencias proteicas associadas aos hits detectados pelo HMMER.

Script usado:

- `scripts/get_Seq_results.sh`

O que a aba faz:

- localiza automaticamente a tabela `.tbl` em `02_HMMER_Results`;
- localiza `01_Predicted_Proteins/predicted_proteins.faa`;
- extrai IDs unicos da primeira coluna da tabela HMMER;
- usa `seqtk subseq` para recuperar as sequencias correspondentes.

Saidas principais:

- `protein_ids_to_extract.txt` temporario;
- `lectin_hits.faa` com as sequencias extraidas.

### 5. Analise & Cortes

Funcao: visualizar hits, selecionar alvos e extrair/cortar sequencias especificas.

Scripts usados:

- essa aba usa funcoes internas do `marauders.py`;
- nao chama um script externo diretamente.

O que a aba faz:

- carrega automaticamente uma tabela `.tbl` da pasta `SRA_ID/02_HMMER_Results`;
- tambem permite selecionar uma tabela `.tbl` manualmente;
- filtra hits com `full_E-value < 0.0001`;
- mostra os 25 primeiros hits;
- ao clicar em uma linha, transfere o `target_name` para o extrator Severus Snap(e);
- localiza automaticamente `SRA_ID/01_Predicted_Proteins/predicted_proteins.faa` ou aceita um FASTA manual;
- salva sequencias encontradas pelo termo informado;
- quando possivel, corta a regiao usando coordenadas presentes na descricao FASTA;
- gera histograma de tamanhos das sequencias com Matplotlib.

## Scripts da Pasta `scripts/`

### `SRAget.sh`

Baixa e extrai leituras SRA. Usa `prefetch`, `fasterq-dump`, fallback com `fastq-dump` e compactacao com `pigz`.

### `run_assembly.sh`

Executa o pipeline de montagem. Usa FastQC, Trimmomatic, BBNorm, MultiQC e um dos montadores: MEGAHIT, SPAdes ou Trinity.

### `ProtSearch.sh`

Prediz proteinas com Prodigal e busca dominios com HMMER (`hmmsearch`).

### `get_Seq_results.sh`

Extrai IDs de proteinas da tabela HMMER e recupera sequencias correspondentes com SeqTK.

### `map_reads.sh`

Script auxiliar para mapear reads contra um FASTA de referencia. Usa Bowtie2 e Samtools, criando resultados em `06_Mapping_Results`. Atualmente nao e chamado diretamente por nenhuma das 5 abas principais da GUI.

### `setup.sh`

Script historico de instalacao local via `apt`/`pip`. Com Docker ou Singularity, a instalacao das dependencias passa a ser feita pela imagem, entao esse script tende a ser desnecessario para usuarios de container.

## Estrutura de Saida

Um fluxo completo para um accession SRA tende a criar uma estrutura como:

```text
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

Algumas pastas podem variar conforme as opcoes escolhidas, o tipo de amostra e a etapa executada.

## Observacoes

- O programa depende de display grafico X11 quando executado em Docker.
- Os resultados sao gerados dentro da pasta do projeto quando a execucao usa o `make docker-run`, pois o diretorio atual e montado em `/app`.
- Trinity aumenta significativamente o tamanho da imagem, mas permite manter o fluxo de transcriptomica na mesma imagem que MEGAHIT e SPAdes.
- Se o Trinity reclamar que `salmon` nao esta instalado, reconstrua a imagem atualizada. O `salmon` e uma dependencia de runtime do Trinity e deve estar presente na imagem Docker/Singularity.
- Avisos do MultiQC sobre arquivos `assets/js/packages/*.js` ausentes podem aparecer em algumas versoes empacotadas pelo Debian. Quando o log termina com `MultiQC complete`, esses avisos nao interrompem o pipeline; eles afetam apenas recursos embutidos no relatorio HTML.
