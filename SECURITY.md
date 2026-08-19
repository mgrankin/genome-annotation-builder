# Security and genomic privacy

Human genomic data is identifying and sensitive. Do not open an issue with a
real VCF, BAM/CRAM, FASTQ, `.genome` bundle, sample name, access URL, credential,
or clinical result.

The repository uses a default-deny `.gitignore`, but Git ignore rules are not a
security boundary. Before every publication or pull request, run:

```bash
python scripts/privacy_scan.py --tracked --deny-file .privacy-denylist
git status --ignored
```

Build releases into a new empty audit directory and scan the archive payloads
with the same deployment-specific deny-list:

```bash
python -m build --outdir /path/to/new/release-audit
python scripts/verify_distributions.py /path/to/new/release-audit \
  --deny-file .privacy-denylist
```

`.privacy-denylist` is optional, local, and ignored. Put one deployment-specific
identifier per line in it; the scan reports matches without echoing the private
value.

Use only deterministic synthetic fixtures in tests. Store local configuration,
inputs, resource snapshots, work directories, and outputs outside the tracked
tree or in the ignored locations shown in `config/config.example.yaml`.

The builder never uploads data. Network access is limited to explicitly
enabled, checksum-pinned public resources. Disable it with
`resources.allow_download: false` after the resource cache has been populated.

Report a software vulnerability privately through the repository hosting
platform's security advisory mechanism.
