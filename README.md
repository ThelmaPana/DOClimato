# DOClimato

A ML-based climatology of DOC.

## Main idea

The goal is to relate DOC observations to a set of climatologies of environmental variables in order to predict DOC values in the global ocean.

## Scripts

-   `00.prepare_env`: download, clean and format environmental data

-   `01.prepare_doc`: clean and format DOC data

-   `02.assemble`: assemble environmental and DOC data

-   `03.fit_model`: train a XGBoost model to relate DOC to env data

-   `04.pred_quality`: assess the prediction quality of the model

-   `05.doc_projection`: generate global maps of DOC predictions
