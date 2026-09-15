# Drive notebook inventory

The project Drive archive contains historical exploratory notebooks. They are retained as provenance references, not as the authoritative publication workflow, because several use obsolete Drive-mounted paths, earlier split definitions, or superseded performance claims.

| Notebook | Drive ID | Status |
|---|---|---|
| 01_IEDB_data_process | `1uVwGNUFipHcr-AePiZxjxJlw7yd38u8w` | Historical curation notebook; source/provenance reference |
| 02_baseline_modelling_iedb_ifng | `1QydXFUGIu8WgkTOLpa62UDXdCcoXs8w9` | Superseded baseline |
| 02a_multi_model_pipeline_iedb_ifng_clean | `1lvG-mPfB1LfKxD2e0ByLINXXIGP3F6-s` | Superseded pipeline |
| 03_nn_baseline_iedb_ifng | `1mDke8woD8T96LLsn6y8jGAOLV3FeHNR9` | Exploratory neural baseline |
| 04_transformer_baseline_iedb_ifng | `1z_CH5riYtyLOJ_obI2hzMUEfcQ81M0yN` | Exploratory transformer baseline; not the final frozen ESM-2 analysis |
| 05_model_comparison_error_analysis_iedb_ifng_fixed | `13mDMPi372BeYjEs9VMg1iUlX4TcNNRCD` | Historical comparison/error analysis |
| 06_grouped_cv_gridsearch_iedb_ifng | `1Dc_NkRL2uHJmx924xXe8zAH5XMky-Njr` | Historical tuning notebook |
| 07_calibration_thresholding_iedb_ifng | `1fxmq3Fn11G0MyPor2oob0lPddBQPxS1N` | Historical calibration notebook |
| 08_regularised_logistic_comparison_clean | `1DlF4NN1f_e-CoGPj6Vb8RAqaCS0vyyKk` | Historical model comparison |
| 09_stacked_ensemble_lr_rf_xgb_clean | `1MH44oGxzZtC-nr1BUsbr6p66-cCh-DE3` | Historical stacking notebook |
| 10_calibrated_boosted_peptide_group_cv_platt_fixed | `1BiEPbVVt6mln9BtZEzY1TXhtEIpoXp5g` | Historical calibration notebook |

The authoritative notebooks are the six repository notebooks under `notebooks/`. They consume the fixed, validated artifacts under `data/` and produce the registered outputs under `results/`.
