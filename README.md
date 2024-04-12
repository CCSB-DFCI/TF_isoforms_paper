[![DOI](https://zenodo.org/badge/761455587.svg)](https://zenodo.org/doi/10.5281/zenodo.10961030)

This is the code for the preprint **Widespread variation in molecular interactions and regulatory properties among transcription factor isoforms** https://www.biorxiv.org/content/10.1101/2024.03.12.584681v3

The python version used to produce the figures is listed in `.python-version` (3.8.13). [pyenv](https://github.com/pyenv/pyenv) can be used to install that specific version. The dependencies are listed in `requirements.txt` and can be installed with `pip -r requirements.txt` after setting up a virtual environment.

The supplementary data tables are available in `supp/` or on the biorxiv page.

To run the notebooks, you will need to download the source data which are
available here: https://zenodo.org/records/10957151
```
wget https://zenodo.org/records/10957151/files/TF_isoforms_paper_data_dir.tar.gz
tar -xzf TF_isoforms_paper_data_dir.tar.gz
```
Any original analysis should use the data from the supplementary tables rather
than these source files.