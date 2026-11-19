# New Optimized layouts (NN≈8)

See **`docs/NewFourArrays_Summary.md`** for the full summary (methods, soft-distance tables, regenerate commands) and **`decoys/README.md`** for decoys.

| Array | Layout | Real NN μ | Decoy→real μ |
|-------|--------|-----------|--------------|
| Zipcodes flip | `Zipcodes/output_235x363_flip_g10_h5_d12-16_far10/` | 7.0 | 7.3 |
| HSDB colored | `HSDB/outputs/235x363_g10_h5_d12-16_far8_24+12_nn8_colored/` | 8.0 | 7.9 |
| Dual shared_pool minH6 | `Dual_barcodes/layout_235x363_dual_shared_pool_g10_minH6.csv` | **8.6** | **8.6** |
| Dual fixed_bar2 8×8 | `Dual_barcodes/layout_235x363_dual_fixed_bar2_8x8_g10_nn8.csv` | 7.0 | 7.0 |

All: 235×363, 36 nt, G&lt;10%, HP≤5. Soft Zip/HSDB: **12/14/16/16**. HSDB uses 16-color bar₂ libraries for NN≈8. Dual shared_pool packs halves to **min Hamming ≥ 6**.
