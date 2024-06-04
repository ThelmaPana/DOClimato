# DOClimato

A ML-based climatology of DOC.

## Main idea

Boosted regression trees are used to relate DOC observations to various environmental climatologies. Inferred relationships were extrapolated to the entire globe to compute DOC climatologies and associated uncertainties.

Overall, 8 climatologies are generated:

-   annual climatologies in 4  layers: surface (0 - 10 m), epipelagic (10 - 200 m), mesopelagic (200 - 1000 m) and bathypelagic (\> 1000 m)

```{=html}
<!-- -->
```
-   seasonal climatologies in the surface layer

## Repo organisation

### Data

Where data lives.

### Figures

Generated figures.

### Output

Generated climatologies.

### Scripts

-   `00.prepare_env`: download, clean and format environmental data

-   `01.prepare_doc`: clean and format DOC data

-   `02.assemble`: assemble environmental and DOC data

-   `03a.doc_surf_ann_fit` to `06a.doc_bathy_ann_fit`: fit BRT models for annual climatologies in surface, epipelagic, mesopelagic and bathypelagic layers

-   `03b.doc_surf_ann_assess` to `06a.doc_bathy_ann_assess`: assess fitting of BRT models for annual climatologies in surface, epipelagic, mesopelagic and bathypelagic layers

-   `07a.doc_surf_seas_fit`: fit BRT models for seasonal surface climatologies

-   `07b.doc_surf_seas_assess`: assess fitting of BRT models for seasonal surface climatologies

-   `08.total_doc_content.R`: compute total DOC content by integrating projections across all oceans

-   `09.paper_figs.R`: generate all figures for the paper
