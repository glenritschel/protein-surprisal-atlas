# Notebooks

The pipeline runs as three Colab notebooks, in order. Each clones this repository and installs its own
dependencies, so they can be run standalone.

| Notebook | Stage | Runtime | Persistence |
|---|---|---|---|
| [`01_pilot_gpu.ipynb`](01_pilot_gpu.ipynb) | 500-protein pilot: scoring → exact-vs-approx validation → controls → analysis → pilot associations | **GPU** | in-session |
| [`02_proteome_scoring_gpu.ipynb`](02_proteome_scoring_gpu.ipynb) | Score the full reviewed proteome (~20k), partitioned + resumable | **GPU** | Google Drive |
| [`03_proteome_annotations_cpu.ipynb`](03_proteome_annotations_cpu.ipynb) | Annotate the scored proteome + run proteome-wide associations | **CPU** | Google Drive |

Open in Colab:

- [![Open 01 in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/glenritschel/protein-surprisal-atlas/blob/main/notebooks/01_pilot_gpu.ipynb) &nbsp; `01_pilot_gpu.ipynb`
- [![Open 02 in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/glenritschel/protein-surprisal-atlas/blob/main/notebooks/02_proteome_scoring_gpu.ipynb) &nbsp; `02_proteome_scoring_gpu.ipynb`
- [![Open 03 in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/glenritschel/protein-surprisal-atlas/blob/main/notebooks/03_proteome_annotations_cpu.ipynb) &nbsp; `03_proteome_annotations_cpu.ipynb`

## How to run

**Order.** Run `01` to reproduce and validate the method on the pilot. Run `02` then `03` to build the
full proteome atlas. `03` depends on the scored output from `02`.

**Runtime type.** `01` and `02` need a **GPU** (Runtime → Change runtime type → GPU); `03` is CPU-only.

**Google Drive.** `02` and `03` mount Google Drive and read/write under `MyDrive/psa_proteome/`. This is
what makes the multi-hour runs **resumable**: the proteome scores, and the UniRef/MobiDB annotation
caches, persist there. If Colab disconnects, re-run — each notebook rebuilds and continues from the
cached state on Drive rather than starting over.

**DepMap requires a manual download.** figshare (which hosts the DepMap files) blocks Colab's shared IP
ranges, so the two DepMap CSVs cannot be downloaded from within Colab. Download them once in your own
browser and place them on Drive at `MyDrive/psa_proteome/data/external/`:

- `CRISPRGeneEffect.csv` — https://ndownloader.figshare.com/files/51064667
- `CRISPRInferredCommonEssentials.csv` — https://ndownloader.figshare.com/files/51064916

Notebook `03` verifies these are present before integrating; if they are missing it prints instructions
and can proceed with DepMap skipped.

## A note on committing notebooks

Commit these with outputs **cleared** (Kernel → Restart & Clear Output, or `nbstripout`). They
regenerate their own figures, tables, and logs when run; committing embedded outputs bloats the
repository and serves no purpose.
