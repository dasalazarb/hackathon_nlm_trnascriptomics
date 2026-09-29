# Blood–MSG hackathon pipeline — flat inputs

This version matches the current `/filesystems/` layout shown in the screenshot.

## Expected structure

```text
/filesystems/
├── RawCountFile_filtered_blood.txt
├── RawCountFile_filtered_msg.txt
├── sampletable_blood.txt
├── sampletable_msg.txt
├── sampletable_SjD_grouping_blood.txt
├── sampletable_SjD_grouping_msg.txt
├── D_DE_neg-pos_blood.txt             # optional
├── D_DE_neg-pos_msg.txt               # optional
├── config_paths.R
├── scripts/
│   ├── 00_pair_and_validate.R
│   └── 01_hackathon_paired_blood_msg.R
└── output/
```

The `inputs/` folder is not used by these scripts.

## Configuration

Normally keep:

```r
PROJECT_ROOT <- "/filesystems"
```

If the root changes later, change only that line in `config_paths.R`.

## Run Stage 00

```bash
cd /filesystems
Rscript scripts/00_pair_and_validate.R
```

Then inspect:

```bash
cat output/00_pairing/pairing_summary.txt
```

## Run Stage 01

Only after Stage 00 passes:

```bash
Rscript scripts/01_hackathon_paired_blood_msg.R
```

Then inspect:

```bash
cat output/01_hackathon/hackathon_summary.txt
```

The two `sampletable_*` files retain the historical `neg/pos` condition. The two `sampletable_SjD_grouping_*` files supply `SjD/nonSjD/HV` for the hackathon.
